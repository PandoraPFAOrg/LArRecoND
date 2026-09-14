/**
 *  @file   LArRecoMP/test/PandoraInterface.cc
 *
 *  @brief  Implementation of the lar reco mp application
 *
 *  $Log: $
 */

#include "TFile.h"
#include "TTree.h"

#include "TGeoBBox.h"
#include "TGeoManager.h"
#include "TGeoMatrix.h"
#include "TGeoShape.h"
#include "TGeoVolume.h"

#ifdef USE_EDEPSIM
#include "TG4PrimaryVertex.h"
#endif

#include "Api/PandoraApi.h"
#include "Geometry/LArTPC.h"
#include "Helpers/XmlHelper.h"
#include "Managers/GeometryManager.h"
#include "Managers/PluginManager.h"
#include "Xml/tinyxml.h"

#include "larpandoracontent/LArContent.h"
#include "larpandoracontent/LArControlFlow/MasterAlgorithm.h"
#include "larpandoracontent/LArControlFlow/MultiPandoraApi.h"
#include "larpandoracontent/LArObjects/LArCaloHit.h"
#include "larpandoracontent/LArObjects/LArMCParticle.h"
#include "larpandoracontent/LArPlugins/LArPseudoLayerPlugin.h"
#include "larpandoracontent/LArPlugins/LArRotationalTransformationPlugin.h"

#ifdef LIBTORCH_DL
#include "larpandoradlcontent/LArDLContent.h"
#endif

#include "LArNDContent.h"
#include "LArNDGeomSimple.h"
#include "LArRay.h"
#include "PandoraInterface.h"

#ifdef MONITORING
#include "TApplication.h"
#endif

#include <algorithm>
#include <cmath>
#include <getopt.h>
#include <iostream>
#include <memory>
#include <random>
#include <string>
#include <vector>

using namespace pandora;
using namespace lar_nd_reco;

int main(int argc, char *argv[])
{
    int errorNo(0);
    const Pandora *pPrimaryPandora(nullptr);

    try
    {
        Parameters parameters;

        if (!ParseCommandLine(argc, argv, parameters))
            return 1;

#ifdef MONITORING
        TApplication *pTApplication = new TApplication("LArReco", &argc, argv);
        pTApplication->SetReturnFromRun(kTRUE);
#endif

        pPrimaryPandora = new pandora::Pandora();

        if (!pPrimaryPandora)
            throw StatusCodeException(STATUS_CODE_FAILURE);

        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, LArContent::RegisterAlgorithms(*pPrimaryPandora));
#ifdef LIBTORCH_DL
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, LArDLContent::RegisterAlgorithms(*pPrimaryPandora));
#endif
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, LArContent::RegisterBasicPlugins(*pPrimaryPandora));

        if (parameters.m_use3D)
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, LArNDContent::RegisterAlgorithms(*pPrimaryPandora));

        MultiPandoraApi::AddPrimaryPandoraInstance(pPrimaryPandora);

        LArNDGeomSimple simpleGeom;
        CreateGeometry(parameters, pPrimaryPandora, simpleGeom);
        ProcessExternalParameters(parameters, pPrimaryPandora);
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::SetPseudoLayerPlugin(*pPrimaryPandora, new lar_content::LArPseudoLayerPlugin));
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=,
            PandoraApi::SetLArTransformationPlugin(*pPrimaryPandora, new lar_content::LArRotationalTransformationPlugin));
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::ReadSettings(*pPrimaryPandora, parameters.m_settingsFile));

        ProcessEvents(parameters, pPrimaryPandora, simpleGeom);
    }
    catch (const StatusCodeException &statusCodeException)
    {
        std::cerr << "Pandora StatusCodeException: " << statusCodeException.ToString() << statusCodeException.GetBackTrace() << std::endl;
        errorNo = 1;
    }
    catch (...)
    {
        std::cerr << "Unknown exception: " << std::endl;
        errorNo = 1;
    }

    MultiPandoraApi::DeletePandoraInstances(pPrimaryPandora);
    return errorNo;
}

//------------------------------------------------------------------------------------------------------------------------------------------

namespace lar_nd_reco
{

void CreateGeometry(const Parameters &parameters, const Pandora *const pPrimaryPandora, LArNDGeomSimple &geom)
{
    // Get the geometry info from the appropriate ROOT file
    TFile *fileSource = TFile::Open(parameters.m_geomFileName.c_str(), "READ");
    if (!fileSource)
    {
        std::cout << "Error in CreateGeometry(): can't open file " << parameters.m_geomFileName << std::endl;
        return;
    }

    TGeoManager *pSimGeom = dynamic_cast<TGeoManager *>(fileSource->Get(parameters.m_geomManagerName.c_str()));
    if (!pSimGeom)
    {
        std::cout << "Could not find the geometry manager named " << parameters.m_geomManagerName << std::endl;
        fileSource->Close();
        return;
    }

    // Go through the geometry and find the paths to the nodes we are interested in
    std::vector<std::vector<unsigned int>> nodePaths; // Store the daughter indices in the path to the node
    std::vector<unsigned int> currentPath;
    const std::string nameToFind = parameters.m_useModularGeometry ? parameters.m_sensitiveDetName : parameters.m_geometryVolName;
    RecursiveGeometrySearch(pSimGeom, nameToFind, nodePaths, currentPath);
    std::cout << "Found " << nodePaths.size() << " matches for volumes containing the name " << nameToFind << std::endl;

    // Navigate to each node and use them to build the pandora geometry
    for (unsigned int n = 0; n < nodePaths.size(); ++n)
    {
        const TGeoNode *pTopNode = pSimGeom->GetCurrentNode();
        // We have to multiply together matrices at each depth to convert local coordinates to the world volume
        std::unique_ptr<TGeoHMatrix> pVolMatrix = std::make_unique<TGeoHMatrix>(*pTopNode->GetMatrix());
        for (unsigned int d = 0; d < nodePaths.at(n).size(); ++d)
        {
            pSimGeom->CdDown(nodePaths.at(n).at(d));
            const TGeoNode *pNode = pSimGeom->GetCurrentNode();
            std::unique_ptr<TGeoHMatrix> pMatrix = std::make_unique<TGeoHMatrix>(*pNode->GetMatrix());
            pVolMatrix->Multiply(pMatrix.get());
        }
        const TGeoNode *pTargetNode = pSimGeom->GetCurrentNode();

        MakePandoraTPC(pPrimaryPandora, parameters, geom, pVolMatrix, pTargetNode, n);

        for (const unsigned int &daughter : nodePaths.at(n))
        {
            (void)daughter;
            pSimGeom->CdUp();
        }
    }
    std::cout << "Created " << nodePaths.size() << " TPCs" << std::endl;

    fileSource->Close();
}

//------------------------------------------------------------------------------------------------------------------------------------------

void RecursiveGeometrySearch(TGeoManager *pSimGeom, const std::string &targetName, std::vector<std::vector<unsigned int>> &nodePaths,
    std::vector<unsigned int> &currentPath)
{
    const std::string nodeName{pSimGeom->GetCurrentNode()->GetName()};
    if (nodeName.find(targetName) != std::string::npos)
    {
        nodePaths.emplace_back(currentPath);
    }
    else
    {
        for (unsigned int i = 0; i < pSimGeom->GetCurrentNode()->GetNdaughters(); ++i)
        {
            pSimGeom->CdDown(i);
            currentPath.emplace_back(i);
            RecursiveGeometrySearch(pSimGeom, targetName, nodePaths, currentPath);
            pSimGeom->CdUp();
            currentPath.pop_back();
        }
    }
    return;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void MakePandoraTPC(const pandora::Pandora *const pPrimaryPandora, const Parameters &parameters, LArNDGeomSimple &geom,
    const std::unique_ptr<TGeoHMatrix> &pVolMatrix, const TGeoNode *pTargetNode, const unsigned int tpcNumber)
{
    // Get the BBox dimensions from the placement volume, which is assumed to be a cube
    TGeoVolume *pCurrentVol = pTargetNode->GetVolume();
    TGeoShape *pCurrentShape = pCurrentVol->GetShape();
    // pCurrentShape->InspectShape();
    TGeoBBox *pBox = dynamic_cast<TGeoBBox *>(pCurrentShape);

    // Now can get origin/width data from the BBox
    const double dx = pBox->GetDX() * parameters.m_lengthScale; // Note these are the half widths
    const double dy = pBox->GetDY() * parameters.m_lengthScale;
    const double dz = pBox->GetDZ() * parameters.m_lengthScale;
    const double *pOrigin = pBox->GetOrigin();

    // std::cout << "Origin = (" << pOrigin[0] << ", " << pOrigin[1] << ", " << pOrigin[2] << ")" << std::endl;

    // Translate local origin to global coordinates
    double level1[3] = {0.0, 0.0, 0.0};
    pTargetNode->LocalToMasterVect(pOrigin, level1);

    // std::cout << "Level1 = (" << level1[0] << ", " << level1[1] << ", " << level1[2] << ")" << std::endl;

    // Can now create a geometry using the found parameters
    PandoraApi::Geometry::LArTPC::Parameters geoparameters;

    try
    {
        const double *pVolTrans = pVolMatrix->GetTranslation();
        const double centreX = (level1[0] + pVolTrans[0]) * parameters.m_lengthScale;
        const double centreY = (level1[1] + pVolTrans[1]) * parameters.m_lengthScale;
        const double centreZ = (level1[2] + pVolTrans[2]) * parameters.m_lengthScale;
        geoparameters.m_centerX = centreX;
        geoparameters.m_centerY = centreY;
        geoparameters.m_centerZ = centreZ;
        geoparameters.m_widthX = dx * 2.0;
        geoparameters.m_widthY = dy * 2.0;
        geoparameters.m_widthZ = dz * 2.0;

        // ATTN: parameters past here taken from uboone
        geoparameters.m_larTPCVolumeId = tpcNumber;
        geoparameters.m_wirePitchU = 0.300000011921;
        geoparameters.m_wirePitchV = 0.300000011921;
        geoparameters.m_wirePitchW = 0.300000011921;
        geoparameters.m_wireAngleU = 1.04719758034;
        geoparameters.m_wireAngleV = -1.04719758034;
        geoparameters.m_wireAngleW = 0.0;
        geoparameters.m_sigmaUVW = 1;
        geoparameters.m_isDriftInPositiveX = tpcNumber % 2;

        geom.AddTPC(centreX - dx, centreX + dx, centreY - dy, centreY + dy, centreZ - dz, centreZ + dz, tpcNumber);

        std::cout << "Creating TPC: " << centreX - dx << ", " << centreX + dx << ", " << centreY - dy << ", " << centreY + dy << ", "
                  << centreZ - dz << ", " << centreZ + dz << std::endl;
    }
    catch (const pandora::StatusCodeException &)
    {
        std::cout << "CreatePandoraLArTPCs - invalid tpc parameter provided" << std::endl;
    }

    try
    {
        PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::Geometry::LArTPC::Create(*pPrimaryPandora, geoparameters));
    }
    catch (const pandora::StatusCodeException &)
    {
        std::cout << "CreatePandoraLArTPCs - unable to create tpc, insufficient or "
                     "invalid information supplied"
                  << std::endl;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ProcessEvents(const Parameters &parameters, const Pandora *const pPrimaryPandora, const LArNDGeomSimple &geom)
{
    if (parameters.m_dataFormat == Parameters::LArNDFormat::EDepSim)
    {
#ifdef USE_EDEPSIM
        ProcessEDepSimEvents(parameters, pPrimaryPandora, geom);
#endif
    }
    else if (parameters.m_dataFormat == Parameters::LArNDFormat::SED)
    {
        ProcessSEDEvents(parameters, pPrimaryPandora, geom);
    }
    else
    {
        ProcessSPEvents(parameters, pPrimaryPandora, geom);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ProcessSPEvents(const Parameters &parameters, const Pandora *const pPrimaryPandora, const LArNDGeomSimple &geom)
{

    TFile *fileSource = TFile::Open(parameters.m_inputFileName.c_str(), "READ");
    if (!fileSource)
    {
        std::cout << "Error in ProcessSPEvents(): can't open file " << parameters.m_inputFileName << std::endl;
        return;
    }

    TTree *ndsptree = dynamic_cast<TTree *>(fileSource->Get(parameters.m_inputTreeName.c_str()));
    if (!ndsptree)
    {
        std::cout << "Could not find the event tree " << parameters.m_inputTreeName << std::endl;
        fileSource->Close();
        return;
    }

    std::unique_ptr<LArSP> larsp =
        parameters.m_dataFormat == Parameters::LArNDFormat::SPMC ? std::make_unique<LArSPMC>(ndsptree) : std::make_unique<LArSP>(ndsptree);

    // Factory for creating LArCaloHits
    lar_content::LArCaloHitFactory m_larCaloHitFactory;

    // Voxel width
    const float voxelWidth(parameters.m_voxelWidth);

    // Total number of entries in the TTree
    const int nEntries(ndsptree->GetEntries());

    // Starting event
    const int startEvt = parameters.m_nEventsToSkip > 0 ? parameters.m_nEventsToSkip : 0;
    // Number of events to process, up to nEntries
    const int nProcess = parameters.m_nEventsToProcess > 0 ? parameters.m_nEventsToProcess : nEntries;
    // End event, up to nEntries
    const int endEvt = (startEvt + nProcess) < nEntries ? startEvt + nProcess : nEntries;

    std::cout << "Start event is " << startEvt << " and end event is " << endEvt - 1 << std::endl;

    for (int iEvt = startEvt; iEvt < endEvt; iEvt++)
    {
        if (parameters.m_shouldDisplayEventNumber)
            std::cout << std::endl << "   PROCESSING EVENT: " << iEvt << std::endl << std::endl;

        ndsptree->GetEntry(iEvt);

        // Stop processing the event if we have too many space points: reco takes too long
        const int nSP = larsp->m_x->size();
        if (parameters.m_maxMergedVoxels > 0 && nSP > parameters.m_maxMergedVoxels)
        {
            std::cout << "SKIPPING EVENT: number of space points " << nSP << " > " << parameters.m_maxMergedVoxels << std::endl;
            continue;
        }
        // Also stop processing the event if it has too few hits (we can't make CaloHits for essentially empty events)
        if (nSP < parameters.m_minNSpacePoints)
        {
            std::cout << "SKIPPING EVENT: number of space points " << nSP << " < " << parameters.m_minNSpacePoints << std::endl;
            continue;
        }

        // Some truth information first
        if (parameters.m_dataFormat == Parameters::LArNDFormat::SPMC)
        {
            LArSPMC *larspmc = dynamic_cast<LArSPMC *>(larsp.get());
            CreateSPMCParticles(*larspmc, pPrimaryPandora, parameters);
        }

        int hitCounter(0);

        // Loop over the space points and make them into caloHits
        for (size_t isp = 0; isp < nSP; ++isp)
        {
            const float voxelX = (*larsp->m_x)[isp];
            const float voxelY = (*larsp->m_y)[isp];
            const float voxelZ = (*larsp->m_z)[isp];
            const float voxelE = (*larsp->m_charge)[isp];

	    const float voxelT0 = (*larsp->m_t0)[isp];

            // Skip this hit if its coordinates or energy are NaNs
            if (std::isnan(voxelX) || std::isnan(voxelY) || std::isnan(voxelZ) || std::isnan(voxelE))
            {
                std::cout << "Ignoring hit " << isp << " which contains NaNs: (" << voxelX << ", " << voxelY << ", " << voxelZ
                          << "), E = " << voxelE << std::endl;
                continue;
            }

            const pandora::CartesianVector voxelPos(voxelX, voxelY, voxelZ);
            const float MipE{0.00075};
            const float voxelMipEquivalentE = voxelE / MipE;
            const int tpcID(geom.GetTPCNumber(voxelPos));
            lar_content::LArCaloHitParameters caloHitParameters;
            caloHitParameters.m_positionVector = voxelPos;
            caloHitParameters.m_expectedDirection = pandora::CartesianVector(0.f, 0.f, 1.f);
            caloHitParameters.m_cellNormalVector = pandora::CartesianVector(0.f, 0.f, 1.f);
            caloHitParameters.m_cellGeometry = pandora::RECTANGULAR;
            caloHitParameters.m_cellSize0 = voxelWidth;
            caloHitParameters.m_cellSize1 = voxelWidth;
            caloHitParameters.m_cellThickness = voxelWidth;
            caloHitParameters.m_nCellRadiationLengths = 1.f;
            caloHitParameters.m_nCellInteractionLengths = 1.f;
            caloHitParameters.m_time = voxelT0;
            caloHitParameters.m_inputEnergy = voxelE;
            caloHitParameters.m_mipEquivalentEnergy = voxelMipEquivalentE;
            caloHitParameters.m_electromagneticEnergy = voxelE;
            caloHitParameters.m_hadronicEnergy = voxelE;
            caloHitParameters.m_isDigital = false;
            caloHitParameters.m_hitType = pandora::TPC_3D;
            caloHitParameters.m_hitRegion = pandora::SINGLE_REGION;
            caloHitParameters.m_layer = 0;
            caloHitParameters.m_isInOuterSamplingLayer = false;
            caloHitParameters.m_pParentAddress = (void *)(static_cast<uintptr_t>(++hitCounter));
            caloHitParameters.m_larTPCVolumeId = tpcID < 0 ? 0 : tpcID;
            caloHitParameters.m_daughterVolumeId = 0;

            // Only used for truth
            long trackID{0};
            float energyFrac{0.f};
            // Set calo hit to MCParticle relation using trackID
            if (parameters.m_dataFormat == Parameters::LArNDFormat::SPMC)
            {
                LArSPMC *larspmc = dynamic_cast<LArSPMC *>(larsp.get());
                const std::vector<float> mcContribs = (*larspmc->m_hit_packetFrac)[isp];
                const int biggestContribIndex = std::distance(mcContribs.begin(), std::max_element(mcContribs.begin(), mcContribs.end()));
                const std::vector<long> hitPartIDVect = (*larspmc->m_hit_particleID)[isp];
                trackID = (hitPartIDVect.size() > biggestContribIndex) ? hitPartIDVect[biggestContribIndex] : 0;

                // Due to the merging of hits, the contributions can sometimes add up to more than 1.
                // Normalise first
                const float sum = std::accumulate(mcContribs.begin(), mcContribs.end(), 0.f);
                energyFrac = (biggestContribIndex < mcContribs.size() && std::abs(sum) > 0.0) ? mcContribs[biggestContribIndex] / sum : 0.f;
                // Make sure the energy fraction is not larger than 1
                if (energyFrac > 1.f + std::numeric_limits<float>::epsilon())
                    energyFrac = 1.f;

                if (std::find(larspmc->m_mcp_id->begin(), larspmc->m_mcp_id->end(), trackID) == larspmc->m_mcp_id->end())
                    std::cout << "Problem? Could not find MC particle with file ID " << trackID << std::endl;
            }

            if (parameters.m_use3D)
                PANDORA_THROW_RESULT_IF(
                    pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::CaloHit::Create(*pPrimaryPandora, caloHitParameters, m_larCaloHitFactory));

            if (parameters.m_dataFormat == Parameters::LArNDFormat::SPMC)
                PandoraApi::SetCaloHitToMCParticleRelationship(*pPrimaryPandora, (void *)((intptr_t)hitCounter), (void *)((intptr_t)trackID), energyFrac);

            if (parameters.m_useLArTPC)
            {
                // Create LArCaloHits for U, V and W views assuming x is the common drift coordinate
                const float x0_cm(voxelPos.GetX());
                const float y0_cm(voxelPos.GetY());
                const float z0_cm(voxelPos.GetZ());

                // U view
                lar_content::LArCaloHitParameters caloHitPars_UView(caloHitParameters);
                caloHitPars_UView.m_hitType = pandora::TPC_VIEW_U;
                caloHitPars_UView.m_pParentAddress = (void *)(intptr_t(++hitCounter));
                const float upos_cm(pPrimaryPandora->GetPlugins()->GetLArTransformationPlugin()->YZtoU(y0_cm, z0_cm));
                caloHitPars_UView.m_positionVector = pandora::CartesianVector(x0_cm, 0.f, upos_cm);

                PANDORA_THROW_RESULT_IF(
                    pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::CaloHit::Create(*pPrimaryPandora, caloHitPars_UView, m_larCaloHitFactory));
                if (parameters.m_dataFormat == Parameters::LArNDFormat::SPMC)
                    PandoraApi::SetCaloHitToMCParticleRelationship(
                        *pPrimaryPandora, (void *)((intptr_t)hitCounter), (void *)((intptr_t)trackID), energyFrac);

                // V view
                lar_content::LArCaloHitParameters caloHitPars_VView(caloHitParameters);
                caloHitPars_VView.m_hitType = pandora::TPC_VIEW_V;
                caloHitPars_VView.m_pParentAddress = (void *)(intptr_t(++hitCounter));
                const float vpos_cm(pPrimaryPandora->GetPlugins()->GetLArTransformationPlugin()->YZtoV(y0_cm, z0_cm));
                caloHitPars_VView.m_positionVector = pandora::CartesianVector(x0_cm, 0.f, vpos_cm);
                PANDORA_THROW_RESULT_IF(
                    pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::CaloHit::Create(*pPrimaryPandora, caloHitPars_VView, m_larCaloHitFactory));
                if (parameters.m_dataFormat == Parameters::LArNDFormat::SPMC)
                    PandoraApi::SetCaloHitToMCParticleRelationship(
                        *pPrimaryPandora, (void *)((intptr_t)hitCounter), (void *)((intptr_t)trackID), energyFrac);
                // W view
                lar_content::LArCaloHitParameters caloHitPars_WView(caloHitParameters);
                caloHitPars_WView.m_hitType = pandora::TPC_VIEW_W;
                caloHitPars_WView.m_pParentAddress = (void *)(intptr_t(++hitCounter));
                const float wpos_cm(pPrimaryPandora->GetPlugins()->GetLArTransformationPlugin()->YZtoW(y0_cm, z0_cm));
                caloHitPars_WView.m_positionVector = pandora::CartesianVector(x0_cm, 0.f, wpos_cm);

                PANDORA_THROW_RESULT_IF(
                    pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::CaloHit::Create(*pPrimaryPandora, caloHitPars_WView, m_larCaloHitFactory));
                if (parameters.m_dataFormat == Parameters::LArNDFormat::SPMC)
                    PandoraApi::SetCaloHitToMCParticleRelationship(
                        *pPrimaryPandora, (void *)((intptr_t)hitCounter), (void *)((intptr_t)trackID), energyFrac);
            }

        } // end space point loop

        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::ProcessEvent(*pPrimaryPandora));
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::Reset(*pPrimaryPandora));
    } // end event loop

    fileSource->Close();
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CreateSPMCParticles(const LArSPMC &larspmc, const pandora::Pandora *const pPrimaryPandora, const Parameters &parameters)
{
    lar_content::LArMCParticleFactory mcParticleFactory;

    const int nNeutrinos(larspmc.m_nuPDG->size());
    std::cout << "Read in " << nNeutrinos << " true neutrinos" << std::endl;

    // Create MC neutrinos. Keep track of the vertex ID's of the neutrinos
    std::map<long, int> vertexIdToIndex;

    for (size_t i = 0; i < nNeutrinos; ++i)
    {
        // Unique vertex ID for the neutrino
        const long vertexID = (*larspmc.m_vertex_id)[i];
        // Store the vertex ID for this neutrino entry
        vertexIdToIndex[vertexID] = i;

        const int neutrinoPDG = (*larspmc.m_nuPDG)[i];
        const std::string reaction = GetNuanceReaction((*larspmc.m_ccnc)[i], (*larspmc.m_mode)[i]);
        const int nuanceCode = GetNuanceCode(reaction);
        const float nuVtxX = (*larspmc.m_nuvtxx)[i] * parameters.m_lengthScale;
        const float nuVtxY = (*larspmc.m_nuvtxy)[i] * parameters.m_lengthScale;
        const float nuVtxZ = (*larspmc.m_nuvtxz)[i] * parameters.m_lengthScale;

        const float nuE = (*larspmc.m_nue)[i] * parameters.m_energyScale;
        const float nuPx = (*larspmc.m_nupx)[i];
        const float nuPy = (*larspmc.m_nupy)[i];
        const float nuPz = (*larspmc.m_nupz)[i];

        lar_content::LArMCParticleParameters mcNeutrinoParameters;
        mcNeutrinoParameters.m_nuanceCode = nuanceCode;
        mcNeutrinoParameters.m_process = lar_content::MC_PROC_INCIDENT_NU;

        mcNeutrinoParameters.m_energy = nuE;
        mcNeutrinoParameters.m_visibleEnergy = nuE;
        mcNeutrinoParameters.m_isCC = ((*larspmc.m_ccnc)[i] == 1);
        mcNeutrinoParameters.m_momentum = pandora::CartesianVector(nuPx, nuPy, nuPz);
        mcNeutrinoParameters.m_vertex = pandora::CartesianVector(nuVtxX, nuVtxY, nuVtxZ);
        mcNeutrinoParameters.m_endDirection = mcNeutrinoParameters.m_momentum;
        mcNeutrinoParameters.m_endpoint = pandora::CartesianVector(nuVtxX, nuVtxY, nuVtxZ);

        mcNeutrinoParameters.m_particleId = neutrinoPDG;
        mcNeutrinoParameters.m_mcParticleType = pandora::MC_3D;
        mcNeutrinoParameters.m_pParentAddress = (void *)((intptr_t)vertexID);

        // INFO: Neutrino has no trajactory points.
        mcNeutrinoParameters.m_nTrajPoints = 0;
        mcNeutrinoParameters.m_trajPoints = {};

        PANDORA_THROW_RESULT_IF(
            pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::MCParticle::Create(*pPrimaryPandora, mcNeutrinoParameters, mcParticleFactory));
    }

    // Create a map for associating the unique (file-based) MC IDs with their corresponding (vertex_id, mcp_idLocal) pairs.
    // This is needed to find the unique ID for parent MC particles, where we only know their vertex_id & mcp_idLocal's.
    // The "mcp_mother" local ID doesn't have a corresponding unique "mcp_file_mother" stored in the input ROOT file.
    // We need to do this before creating the MC particles since the ancestry could be stored in any order
    std::map<std::pair<long, long>, long> mcIDMap;

    for (size_t i = 0; i < larspmc.m_mcp_id->size(); ++i)
    {
        // Store the corresponding unique (file-based) MC ID for each <vertex_id, mcp_idLocal> key pair
        const long mcpID = (*larspmc.m_mcp_id)[i];
        const long mcpVertexID = (*larspmc.m_mcp_vertex_id)[i];
        const long mcpIDLocal = (*larspmc.m_mcp_idLocal)[i];
        mcIDMap[std::make_pair(mcpVertexID, mcpIDLocal)] = mcpID;
    }

    // Create MC particles
    for (size_t i = 0; i < larspmc.m_mcp_id->size(); ++i)
    {
        // LArMCParticle parameters
        lar_content::LArMCParticleParameters mcParticleParameters;

        // Initial momentum and energy in GeV
        const float px = (*larspmc.m_mcp_px)[i] * parameters.m_energyScale;
        const float py = (*larspmc.m_mcp_py)[i] * parameters.m_energyScale;
        const float pz = (*larspmc.m_mcp_pz)[i] * parameters.m_energyScale;
        const float energy = (*larspmc.m_mcp_energy)[i] * parameters.m_energyScale;
        mcParticleParameters.m_energy = energy;
        mcParticleParameters.m_visibleEnergy = energy;
        mcParticleParameters.m_momentum = pandora::CartesianVector(px, py, pz);
        mcParticleParameters.m_endDirection = mcParticleParameters.m_momentum;

        // Particle codes
        mcParticleParameters.m_particleId = (*larspmc.m_mcp_pdg)[i];
        mcParticleParameters.m_mcParticleType = pandora::MC_3D;
        mcParticleParameters.m_isCC = false;

        // Neutrino info
        const long mcpVertexID = (*larspmc.m_mcp_vertex_id)[i];
        const int nuIndex = vertexIdToIndex[mcpVertexID];
        const std::string reaction = GetNuanceReaction((*larspmc.m_ccnc)[nuIndex], (*larspmc.m_mode)[nuIndex]);
        mcParticleParameters.m_nuanceCode = GetNuanceCode(reaction);

        // Unique file-based ID for this MC particle
        const long mcpID = (*larspmc.m_mcp_id)[i];
        mcParticleParameters.m_pParentAddress = (void *)((intptr_t)mcpID);

        // Start and end points in cm
        const float startx = (*larspmc.m_mcp_startx)[i] * parameters.m_lengthScale;
        const float starty = (*larspmc.m_mcp_starty)[i] * parameters.m_lengthScale;
        const float startz = (*larspmc.m_mcp_startz)[i] * parameters.m_lengthScale;
        mcParticleParameters.m_vertex = pandora::CartesianVector(startx, starty, startz);

        const float endx = (*larspmc.m_mcp_endx)[i] * parameters.m_lengthScale;
        const float endy = (*larspmc.m_mcp_endy)[i] * parameters.m_lengthScale;
        const float endz = (*larspmc.m_mcp_endz)[i] * parameters.m_lengthScale;
        mcParticleParameters.m_endpoint = pandora::CartesianVector(endx, endy, endz);

        // Process ID
        mcParticleParameters.m_process = lar_content::MC_PROC_UNKNOWN;

        // Empty trajectory info for now.
        mcParticleParameters.m_nTrajPoints = 0;
        mcParticleParameters.m_trajPoints = {};

        // Create MCParticle
        try
        {
            PANDORA_THROW_RESULT_IF(
                pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::MCParticle::Create(*pPrimaryPandora, mcParticleParameters, mcParticleFactory));
        }
        catch (const pandora::StatusCodeException &)
        {
            std::cout << "Unable to create MCParticle " << i << " : invalid info supplied, e.g. non-unique trackID or NaNs" << std::endl;
            continue;
        }

        // Set parent relationship. For the parent, use its <vertex_id, mcp_idLocal> pair to get its unique ID
        const long mcpMotherID = (*larspmc.m_mcp_mother)[i];
        const std::pair<long, long> parentPair = std::make_pair(mcpVertexID, mcpMotherID);
        const long mcpParentID = (mcIDMap.find(parentPair) != mcIDMap.end()) ? mcIDMap.at(parentPair) : mcpMotherID;

        if (mcpParentID == -1) // link to mc neutrino
        {
            PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=,
                PandoraApi::SetMCParentDaughterRelationship(*pPrimaryPandora, (void *)((intptr_t)mcpVertexID), (void *)((intptr_t)mcpID)));
        }
        else
        {
            PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=,
                PandoraApi::SetMCParentDaughterRelationship(*pPrimaryPandora, (void *)((intptr_t)mcpParentID), (void *)((intptr_t)mcpID)));
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

#ifdef USE_EDEPSIM
void ProcessEDepSimEvents(const Parameters &parameters, const Pandora *const pPrimaryPandora, const LArNDGeomSimple &geom)
{

    TFile *fileSource = TFile::Open(parameters.m_inputFileName.c_str(), "READ");
    if (!fileSource)
    {
        std::cout << "Error in ProcessEDepSimEvents(): can't open file " << parameters.m_inputFileName << std::endl;
        return;
    }

    TTree *pEDepSimTree = dynamic_cast<TTree *>(fileSource->Get(parameters.m_inputTreeName.c_str()));
    if (!pEDepSimTree)
    {
        std::cout << "Could not find the event tree " << parameters.m_inputTreeName << std::endl;
        fileSource->Close();
        return;
    }

    TG4Event *pEDepSimEvent(nullptr);
    pEDepSimTree->SetBranchAddress("Event", &pEDepSimEvent);

    // Factory for creating LArCaloHits
    lar_content::LArCaloHitFactory m_larCaloHitFactory;

    const LArGrid grid = parameters.m_useModularGeometry ? MakeVoxelisationGrid(geom, parameters) : MakeVoxelisationGrid(pPrimaryPandora, parameters);

    std::cout << "Total grid volume: bot = " << grid.m_bottom << "\n top = " << grid.m_top << std::endl;
    std::cout << "Making voxels with size " << grid.m_binWidths << std::endl;

    // Total number of entries in the TTree
    const int nEntries(pEDepSimTree->GetEntries());

    // Starting event
    const int startEvt = parameters.m_nEventsToSkip > 0 ? parameters.m_nEventsToSkip : 0;
    // Number of events to process, up to nEntries
    const int nProcess = parameters.m_nEventsToProcess > 0 ? parameters.m_nEventsToProcess : nEntries;
    // End event, up to nEntries
    const int endEvt = (startEvt + nProcess) < nEntries ? startEvt + nProcess : nEntries;

    std::cout << "Start event is " << startEvt << " and end event is " << endEvt - 1 << std::endl;

    for (int iEvt = startEvt; iEvt < endEvt; iEvt++)
    {
        if (parameters.m_shouldDisplayEventNumber)
            std::cout << std::endl << "   PROCESSING EVENT: " << iEvt << std::endl << std::endl;

        pEDepSimTree->GetEntry(iEvt);

        if (!pEDepSimEvent)
            return;

        // Create MCParticles from Geant4 trajectories
        const MCParticleEnergyMap MCEnergyMap = CreateEDepSimMCParticles(*pEDepSimEvent, pPrimaryPandora, parameters);

        int hitCounter{0};

        // Loop over (EDep) hits, which are stored in the hit segment detectors.
        // Only process hits from the detector we are interested in
        for (TG4HitSegmentDetectors::iterator detector = pEDepSimEvent->SegmentDetectors.begin();
             detector != pEDepSimEvent->SegmentDetectors.end(); ++detector)
        {
            if (detector->first.find(parameters.m_sensitiveDetName) == std::string::npos)
            {
                std::cout << "Skipping sensitive detector " << detector->first << "; expecting " << parameters.m_sensitiveDetName << std::endl;
                continue;
            }

            std::cout << "Show hits for " << detector->first << " (" << detector->second.size() << " hits)" << std::endl;
            std::cout << "                                 " << std::endl;

            LArVoxelList voxelList;

            // Loop over hit segments and create voxels from them
            for (TG4HitSegment &g4Hit : detector->second)
            {
                const TLorentzVector &hitStart = g4Hit.GetStart();
                const TLorentzVector &hitStop = g4Hit.GetStop();
                const pandora::CartesianVector start(hitStart.X(), hitStart.Y(), hitStart.Z());
                const pandora::CartesianVector end(hitStop.X(), hitStop.Y(), hitStop.Z());
                const float energy = g4Hit.GetEnergyDeposit();
                const int g4id = g4Hit.GetContributors()[0];

                const LArHitInfo hitInfo(start, end, energy, g4id, parameters.m_lengthScale, parameters.m_energyScale);
                const LArVoxelList currentVoxelList = MakeVoxels(hitInfo, grid, parameters, geom);

                for (const LArVoxel &voxel : currentVoxelList)
                    voxelList.emplace_back(voxel);
            }

            std::cout << "Produced " << voxelList.size() << " voxels from " << detector->second.size() << " hit segments." << std::endl;

            // Merge voxels with the same IDs
            const LArVoxelList mergedVoxels = MergeSameVoxels(voxelList);

            std::cout << "Produced " << mergedVoxels.size() << " merged voxels from " << voxelList.size() << " voxels." << std::endl;
            voxelList.clear();

            // Stop processing the event if we have too many voxels: reco takes too long
            if (parameters.m_maxMergedVoxels > 0 && mergedVoxels.size() > parameters.m_maxMergedVoxels)
            {
                std::cout << "SKIPPING EVENT: number of merged voxels " << mergedVoxels.size() << " > " << parameters.m_maxMergedVoxels << std::endl;
                break;
            }

            MakeCaloHitsFromVoxels(mergedVoxels, MCEnergyMap, pPrimaryPandora, parameters, hitCounter);
        } // end segment detector loop

        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::ProcessEvent(*pPrimaryPandora));
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::Reset(*pPrimaryPandora));
    }

    // Close input file
    fileSource->Close();
}

//------------------------------------------------------------------------------------------------------------------------------------------

MCParticleEnergyMap CreateEDepSimMCParticles(const TG4Event &event, const pandora::Pandora *const pPrimaryPandora, const Parameters &parameters)
{
    // Create map of trackID and energy
    MCParticleEnergyMap energyMap;

    if (!pPrimaryPandora)
    {
        std::cout << "Could not create MC particles, since pPrimaryPandora is null" << std::endl;
        return energyMap;
    }

    lar_content::LArMCParticleFactory mcParticleFactory;

    // Loop over the initial primary neutrinos, storing their IDs and vertex positions inside vectors
    // since we need these to work out the associated neutrino ancestors for all MC trajectories
    std::vector<pandora::CartesianVector> neutrinoVertices;
    std::vector<int> neutrinoIDVector, nuanceCodeVector;

    for (size_t i = 0; i < event.Primaries.size(); ++i)
    {
        const TG4PrimaryVertex &g4PrimaryVtx = event.Primaries[i];
        const TLorentzVector neutrinoVtx = g4PrimaryVtx.GetPosition() * parameters.m_lengthScale;
        std::cout << "Neutrino vertex = " << neutrinoVtx.X() << ", " << neutrinoVtx.Y() << ", " << neutrinoVtx.Z() << std::endl;

        const std::string reaction(g4PrimaryVtx.GetReaction());
        const int nuanceCode = GetNuanceCode(reaction);

        // Get the primary vertex particle information
        if (g4PrimaryVtx.Informational.size() > 0)
        {
            const TG4PrimaryVertex &g4Info = g4PrimaryVtx.Informational[0];

            // Get the first primary particle, which should be the neutrino.
            // Other primaries would be nuclei etc.
            if (g4Info.Particles.size() > 0)
            {
                const TG4PrimaryParticle &g4Primary = g4Info.Particles[0];

                // The primary neutrinoIDs are usually the same value in a full spill event, e.g. -2.
                // Introduce an artificial offset (the primary vertex number "i") to give unique IDs
                const int neutrinoID = g4Primary.GetTrackId() - i;
                const int neutrinoPDG = g4Primary.GetPDGCode();
                const TLorentzVector neutrinoP4 = g4Primary.GetMomentum() * parameters.m_energyScale;

                std::cout << "Neutrino ID = " << neutrinoID << ", PDG = " << neutrinoPDG << ", E = " << neutrinoP4.E()
                          << ", px = " << neutrinoP4.Px() << ", py = " << neutrinoP4.Py() << ", pz = " << neutrinoP4.Pz() << std::endl;

                lar_content::LArMCParticleParameters mcNeutrinoParameters;
                mcNeutrinoParameters.m_nuanceCode = nuanceCode;
                mcNeutrinoParameters.m_isCC = (g4PrimaryVtx.GetReaction().find("CC") != std::string::npos);
                mcNeutrinoParameters.m_process = lar_content::MC_PROC_INCIDENT_NU;
                mcNeutrinoParameters.m_energy = neutrinoP4.E();
                mcNeutrinoParameters.m_visibleEnergy = neutrinoP4.E();
                mcNeutrinoParameters.m_momentum = pandora::CartesianVector(neutrinoP4.Px(), neutrinoP4.Py(), neutrinoP4.Pz());
                mcNeutrinoParameters.m_endDirection = mcNeutrinoParameters.m_momentum;
                mcNeutrinoParameters.m_vertex = pandora::CartesianVector(neutrinoVtx.X(), neutrinoVtx.Y(), neutrinoVtx.Z());
                mcNeutrinoParameters.m_endpoint = pandora::CartesianVector(neutrinoVtx.X(), neutrinoVtx.Y(), neutrinoVtx.Z());
                mcNeutrinoParameters.m_particleId = neutrinoPDG;
                mcNeutrinoParameters.m_mcParticleType = pandora::MC_3D;
                mcNeutrinoParameters.m_pParentAddress = (void *)((intptr_t)neutrinoID);
                mcNeutrinoParameters.m_nTrajPoints = 0;
                mcNeutrinoParameters.m_trajPoints = {};

                PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=,
                    PandoraApi::MCParticle::Create(*pPrimaryPandora, mcNeutrinoParameters, mcParticleFactory));

                // Keep track of neutrino vertex, ID and reaction code
                neutrinoVertices.emplace_back(mcNeutrinoParameters.m_vertex.Get());
                neutrinoIDVector.emplace_back(neutrinoID);
                nuanceCodeVector.emplace_back(nuanceCode);
            }
        }
    }

    std::cout << "Creating MC Particles from the Geant4 event trajectories" << std::endl;

    // Loop over trajectories
    std::cout << "Number of trajectories = " << event.Trajectories.size() << std::endl;

    // Keep track of the primary Nuance codes for the trajectories using a map[trackID] container.
    // The trackIDs and their parentIDs will be in cascading historical order in the trajectory loop,
    // meaning that a given trajectory's parentID will have been previously stored in the map
    std::map<int, int> trajNuanceCodes;

    for (const TG4Trajectory &g4Traj : event.Trajectories)
    {
        // LArMCParticle parameters
        lar_content::LArMCParticleParameters mcParticleParameters;

        // Initial momentum and energy in GeV (Geant4 uses MeV)
        const TLorentzVector initMtm(g4Traj.GetInitialMomentum() * parameters.m_energyScale);
        const float energy(initMtm.E());
        mcParticleParameters.m_energy = energy;
        mcParticleParameters.m_visibleEnergy = energy;
        mcParticleParameters.m_momentum = pandora::CartesianVector(initMtm.X(), initMtm.Y(), initMtm.Z());

        // Particle codes
        mcParticleParameters.m_particleId = g4Traj.GetPDGCode();
        mcParticleParameters.m_mcParticleType = pandora::MC_3D;
        mcParticleParameters.m_isCC = false;

        // Set unique parent integer address using trackID
        const int trackID = g4Traj.GetTrackId();
        mcParticleParameters.m_pParentAddress = (void *)((intptr_t)trackID);

        // Start and end points in cm (Geant4 uses mm)
        const std::vector<TG4TrajectoryPoint> trajPoints = g4Traj.Points;
        const int nPoints(trajPoints.size());
        mcParticleParameters.m_nTrajPoints = nPoints;
        CartesianPointVector trajPointsVector({});

        if (nPoints > 1)
        {
            const TG4TrajectoryPoint start = trajPoints[0];
            const TLorentzVector vertex = start.GetPosition() * parameters.m_lengthScale;
            mcParticleParameters.m_vertex = pandora::CartesianVector(vertex.X(), vertex.Y(), vertex.Z());

            const TG4TrajectoryPoint end = trajPoints[nPoints - 1];
            const TLorentzVector endPos = end.GetPosition() * parameters.m_lengthScale;
            mcParticleParameters.m_endpoint = pandora::CartesianVector(endPos.X(), endPos.Y(), endPos.Z());
            mcParticleParameters.m_endDirection = pandora::CartesianVector(end.GetMomentum().X(), end.GetMomentum().Y(), end.GetMomentum().Z());

            // Process ID
            mcParticleParameters.m_process = start.GetProcess();

            for (const TG4TrajectoryPoint &trajPoint : trajPoints)
            {
                const TLorentzVector pos = trajPoint.GetPosition() * parameters.m_lengthScale;
                const pandora::CartesianVector cartPos(pos.X(), pos.Y(), pos.Z());
                trajPointsVector.emplace_back(cartPos);
            }
        }
        else
        {
            // Should not reach here, but set sensible values just in case
            mcParticleParameters.m_vertex = pandora::CartesianVector(0.f, 0.f, 0.f);
            mcParticleParameters.m_endpoint = pandora::CartesianVector(0.f, 0.f, 0.f);
            mcParticleParameters.m_endDirection = pandora::CartesianVector(0.f, 0.f, 0.f);
            mcParticleParameters.m_process = lar_content::MC_PROC_UNKNOWN;
        }

        mcParticleParameters.m_trajPoints = trajPointsVector;

        // Set parent relationship and nuance interaction code
        mcParticleParameters.m_nuanceCode = 0;
        const int trajParentID = g4Traj.GetParentId();

        int parentID{trajParentID};
        if (trajParentID < 0) // link to MC neutrino
        {
            // In full spill events, GetParentId() will always return -1 for those particles originating from
            // any neutrino interaction vertex. So we need to compare and match this track's vertex position
            // with the list of primary neutrino vertices to get the required unique neutrino parent ID entry
            // and Nuance reaction code. This also works for single neutrino events

            for (size_t iV = 0; iV < neutrinoVertices.size(); ++iV)
            {
                // Compare the distance squared between the trajectory and neutrino vertices
                const pandora::CartesianVector nuVtx = neutrinoVertices[iV];
                if (nuVtx.GetDistanceSquared(mcParticleParameters.m_vertex.Get()) < std::numeric_limits<float>::epsilon())
                {
                    // Vertex positions match
                    parentID = neutrinoIDVector[iV];
                    const int nuanceValue = nuanceCodeVector[iV];
                    mcParticleParameters.m_nuanceCode = nuanceValue;
                    trajNuanceCodes[trackID] = nuanceValue;
                    break;
                }
            }
        }
        else
        {
            // Retrieve the Nuance code using its parentID entry
            const int nuanceValue = trajNuanceCodes.find(parentID) != trajNuanceCodes.end() ? trajNuanceCodes[parentID] : 0;
            // Set the Nuance code for this trackID; any secondary particle will then retrieve this value
            trajNuanceCodes[trackID] = nuanceValue;
            mcParticleParameters.m_nuanceCode = nuanceValue;
        }

        // Create MCParticle
        PANDORA_THROW_RESULT_IF(
            pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::MCParticle::Create(*pPrimaryPandora, mcParticleParameters, mcParticleFactory));

        // Store the parentID, which will recursively find the primary neutrino if required
        PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=,
            PandoraApi::SetMCParentDaughterRelationship(*pPrimaryPandora, (void *)((intptr_t)parentID), (void *)((intptr_t)trackID)));

        // Store particle energy for given trackID
        energyMap[trackID] = energy;
    }

    return energyMap;
}

#endif

//------------------------------------------------------------------------------------------------------------------------------------------

void ProcessSEDEvents(const Parameters &parameters, const Pandora *const pPrimaryPandora, const LArNDGeomSimple &geom)
{
    std::cout << "About to process SED events" << std::endl;
    TFile *fileSource = TFile::Open(parameters.m_inputFileName.c_str(), "READ");
    if (!fileSource)
    {
        std::cout << "Error in ProcessSEDEvents(): can't open file " << parameters.m_inputFileName << std::endl;
        return;
    }

    TTree *ndsim = dynamic_cast<TTree *>(fileSource->Get(parameters.m_inputTreeName.c_str()));
    if (!ndsim)
    {
        std::cout << "Could not find the event tree " << parameters.m_inputTreeName << std::endl;
        fileSource->Close();
        return;
    }

    const LArSED larsed(ndsim);

    const LArGrid grid = parameters.m_useModularGeometry ? MakeVoxelisationGrid(geom, parameters) : MakeVoxelisationGrid(pPrimaryPandora, parameters);

    std::cout << "Total grid volume: bot = " << grid.m_bottom << "\n top = " << grid.m_top << std::endl;
    std::cout << "Making voxels with size " << grid.m_binWidths << std::endl;

    // Total number of entries in the TTree
    const int nEntries(ndsim->GetEntries());

    // Starting event
    const int startEvt = parameters.m_nEventsToSkip > 0 ? parameters.m_nEventsToSkip : 0;
    // Number of events to process, up to nEntries
    const int nProcess = parameters.m_nEventsToProcess > 0 ? parameters.m_nEventsToProcess : nEntries;
    // End event, up to nEntries
    const int endEvt = (startEvt + nProcess) < nEntries ? startEvt + nProcess : nEntries;

    std::cout << "Start event is " << startEvt << " and end event is " << endEvt - 1 << std::endl;

    for (int iEvt = startEvt; iEvt < endEvt; iEvt++)
    {
        if (parameters.m_shouldDisplayEventNumber)
            std::cout << std::endl << "   PROCESSING EVENT: " << iEvt << std::endl << std::endl;

        ndsim->GetEntry(iEvt);

        // Create MCParticles from Geant4 trajectories
        MCParticleEnergyMap MCEnergyMap;
        for (size_t imcp = 0; imcp < larsed.m_mcp_id->size(); ++imcp)
        {
            MCEnergyMap[(*larsed.m_mcp_id)[imcp]] = (*larsed.m_mcp_energy)[imcp];
        }
        CreateSEDMCParticles(larsed, pPrimaryPandora, parameters);

        LArVoxelList voxelList;

        // Loop over the energy deposits and create voxels
        for (size_t ised = 0; ised < larsed.m_sed_det->size(); ++ised)
        {
            if ((*larsed.m_sed_det)[ised] == parameters.m_sensitiveDetName) // usually volTPCActive
            {
                const float startx = (*larsed.m_sed_startx)[ised];
                const float starty = (*larsed.m_sed_starty)[ised];
                const float startz = (*larsed.m_sed_startz)[ised];
                const float endx = (*larsed.m_sed_endx)[ised];
                const float endy = (*larsed.m_sed_endy)[ised];
                const float endz = (*larsed.m_sed_endz)[ised];
                const float energy = (*larsed.m_sed_energy)[ised] * 1e-3; // sed_energy is in MeV, convert it to GeV
                const int g4id = std::abs((*larsed.m_sed_id)[ised]);

                const pandora::CartesianVector start(startx, starty, startz);
                const pandora::CartesianVector end(endx, endy, endz);

                const LArHitInfo hitInfo(start, end, energy, g4id, parameters.m_lengthScale, parameters.m_energyScale);
                const LArVoxelList currentVoxelList = MakeVoxels(hitInfo, grid, parameters, geom);

                for (const LArVoxel &voxel : currentVoxelList)
                    voxelList.emplace_back(voxel);
            }
        }

        std::cout << "Produced " << voxelList.size() << " voxels from " << larsed.m_sed_det->size() << " hit segments." << std::endl;

        // Merge voxels with the same IDs
        const LArVoxelList mergedVoxels = MergeSameVoxels(voxelList);

        std::cout << "Produced " << mergedVoxels.size() << " merged voxels from " << voxelList.size() << " voxels." << std::endl;
        voxelList.clear();

        // Stop processing the event if we have too many voxels: reco takes too long
        if (parameters.m_maxMergedVoxels > 0 && mergedVoxels.size() > parameters.m_maxMergedVoxels)
        {
            std::cout << "SKIPPING EVENT: number of merged voxels " << mergedVoxels.size() << " > " << parameters.m_maxMergedVoxels << std::endl;
            break;
        }

        int hitCounter{0};
        MakeCaloHitsFromVoxels(mergedVoxels, MCEnergyMap, pPrimaryPandora, parameters, hitCounter);

        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::ProcessEvent(*pPrimaryPandora));
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::Reset(*pPrimaryPandora));
    } // end event loop

    fileSource->Close();
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CreateSEDMCParticles(const LArSED &larsed, const pandora::Pandora *const pPrimaryPandora, const Parameters &parameters)
{

    lar_content::LArMCParticleFactory mcParticleFactory;

    const int nuidoffset(100000000);

    std::cout << "Read in " << larsed.m_nuPDG->size() << " true neutrinos" << std::endl;

    // Create MC neutrinos
    for (size_t i = 0; i < larsed.m_nuPDG->size(); ++i)
    {
        const int neutrinoID = nuidoffset + i;
        const int neutrinoPDG = (*larsed.m_nuPDG)[i];
        const std::string reaction = GetNuanceReaction((*larsed.m_ccnc)[i], (*larsed.m_mode)[i]);
        const int nuanceCode = GetNuanceCode(reaction);

        const float nuVtxX = (*larsed.m_nuvtxx)[i] * parameters.m_lengthScale;
        const float nuVtxY = (*larsed.m_nuvtxy)[i] * parameters.m_lengthScale;
        const float nuVtxZ = (*larsed.m_nuvtxz)[i] * parameters.m_lengthScale;

        const float nuE = (*larsed.m_enu)[i] * parameters.m_energyScale;
        const float nuPx = nuE * (*larsed.m_nu_dcosx)[i];
        const float nuPy = nuE * (*larsed.m_nu_dcosy)[i];
        const float nuPz = nuE * (*larsed.m_nu_dcosz)[i];

        lar_content::LArMCParticleParameters mcNeutrinoParameters;
        mcNeutrinoParameters.m_nuanceCode = nuanceCode;
        mcNeutrinoParameters.m_process = lar_content::MC_PROC_INCIDENT_NU;
        mcNeutrinoParameters.m_isCC = (*larsed.m_ccnc)[i] == 0;

        mcNeutrinoParameters.m_energy = nuE;
        mcNeutrinoParameters.m_visibleEnergy = nuE;
        mcNeutrinoParameters.m_momentum = pandora::CartesianVector(nuPx, nuPy, nuPz);
        mcNeutrinoParameters.m_vertex = pandora::CartesianVector(nuVtxX, nuVtxY, nuVtxZ);
        mcNeutrinoParameters.m_endpoint = pandora::CartesianVector(nuVtxX, nuVtxY, nuVtxZ);
        mcNeutrinoParameters.m_endDirection = mcNeutrinoParameters.m_momentum;

        mcNeutrinoParameters.m_particleId = neutrinoPDG;
        mcNeutrinoParameters.m_mcParticleType = pandora::MC_3D;
        mcNeutrinoParameters.m_pParentAddress = (void *)((intptr_t)neutrinoID);
        mcNeutrinoParameters.m_nTrajPoints = 0;
        mcNeutrinoParameters.m_trajPoints = {};

        PANDORA_THROW_RESULT_IF(
            pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::MCParticle::Create(*pPrimaryPandora, mcNeutrinoParameters, mcParticleFactory));
    }

    // Create MC particles
    for (size_t i = 0; i < larsed.m_mcp_id->size(); ++i)
    {
        // LArMCParticle parameters
        lar_content::LArMCParticleParameters mcParticleParameters;

        // Initial momentum and energy in GeV
        const float px = (*larsed.m_mcp_px)[i] * parameters.m_energyScale;
        const float py = (*larsed.m_mcp_py)[i] * parameters.m_energyScale;
        const float pz = (*larsed.m_mcp_pz)[i] * parameters.m_energyScale;
        const float energy = (*larsed.m_mcp_energy)[i] * parameters.m_energyScale;
        mcParticleParameters.m_energy = energy;
        mcParticleParameters.m_visibleEnergy = energy;
        mcParticleParameters.m_momentum = pandora::CartesianVector(px, py, pz);

        // Particle codes
        mcParticleParameters.m_particleId = (*larsed.m_mcp_pdg)[i];
        mcParticleParameters.m_mcParticleType = pandora::MC_3D;
        mcParticleParameters.m_isCC = false;

        // Neutrino info
        const int nuid = (*larsed.m_mcp_nuid)[i];
        const int neutrinoID = nuid + nuidoffset;
        const std::string reaction = GetNuanceReaction((*larsed.m_ccnc)[nuid], (*larsed.m_mode)[nuid]);
        mcParticleParameters.m_nuanceCode = GetNuanceCode(reaction);

        // Set unique parent integer address using trackID
        const int trackID = (*larsed.m_mcp_id)[i];
        mcParticleParameters.m_pParentAddress = (void *)((intptr_t)trackID);

        // Start and end points in cm
        const float startx = (*larsed.m_mcp_startx)[i] * parameters.m_lengthScale;
        const float starty = (*larsed.m_mcp_starty)[i] * parameters.m_lengthScale;
        const float startz = (*larsed.m_mcp_startz)[i] * parameters.m_lengthScale;
        mcParticleParameters.m_vertex = pandora::CartesianVector(startx, starty, startz);

        const float endx = (*larsed.m_mcp_endx)[i] * parameters.m_lengthScale;
        const float endy = (*larsed.m_mcp_endy)[i] * parameters.m_lengthScale;
        const float endz = (*larsed.m_mcp_endz)[i] * parameters.m_lengthScale;
        mcParticleParameters.m_endpoint = pandora::CartesianVector(endx, endy, endz);
        mcParticleParameters.m_endDirection = mcParticleParameters.m_momentum;

        // Process ID
        mcParticleParameters.m_process = lar_content::MC_PROC_UNKNOWN;

        mcParticleParameters.m_nTrajPoints = 0;
        mcParticleParameters.m_trajPoints = {};

        // Create MCParticle
        PANDORA_THROW_RESULT_IF(
            pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::MCParticle::Create(*pPrimaryPandora, mcParticleParameters, mcParticleFactory));

        // Set parent relationships
        const int parentID = (*larsed.m_mcp_mother)[i];

        if (parentID == 0) // link to mc neutrino
        {
            PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=,
                PandoraApi::SetMCParentDaughterRelationship(*pPrimaryPandora, (void *)((intptr_t)neutrinoID), (void *)((intptr_t)trackID)));
        }
        else
        {
            PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=,
                PandoraApi::SetMCParentDaughterRelationship(*pPrimaryPandora, (void *)((intptr_t)parentID), (void *)((intptr_t)trackID)));
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

int GetNuanceCode(const std::string &reaction)
{
    // The GENIE reaction string (also stored by edep-sim) is created using
    // https://github.com/GENIE-MC/Generator/blob/master/src/Framework/Interaction/Interaction.cxx#L249
    // String format is "nu:PDGId;tgt:PDGId;N:PDGId;proc:interactionType,scattering;", e.g.
    //                  "nu:14;tgt:1000180400;N:2112;proc:Weak[CC],QES;"

    // GENIE scattering codes:
    // https://github.com/GENIE-MC/Generator/blob/master/src/Framework/Interaction/ScatteringType.h

    // Nuance codes:
    // https://internal.dunescience.org/doxygen/MCNeutrino_8h_source.html

    // GENIE conversion code for RooTracker output files:
    // https://github.com/GENIE-MC/Generator/blob/master/src/contrib/t2k/neut_code_from_rootracker.C
    // Similar code is available here (Neut reaction code):
    // https://internal.dunescience.org/doxygen/namespacegenie_1_1utils_1_1ghep.html

    // For now, just set the basic reaction types, excluding any specific final states:
    // https://github.com/GENIE-MC/Generator/blob/master/src/contrib/t2k/neut_code_from_rootracker.C#L276
    int code(1000);

    const bool is_cc = (reaction.find("Weak[CC]") != std::string::npos); // weak charged-current
    const bool is_nc = (reaction.find("Weak[NC]") != std::string::npos); // weak neutral-current
    // const bool is_charm = (reaction.find("charm")    != std::string::npos); // charm production
    const bool is_qel = (reaction.find("QES") != std::string::npos);   // quasi-elastic scattering
    const bool is_dis = (reaction.find("DIS") != std::string::npos);   // deep inelastic scattering
    const bool is_res = (reaction.find("RES") != std::string::npos);   // resonance
    const bool is_cohpi = (reaction.find("COH") != std::string::npos); // coherent pi
    const bool is_ve = (reaction.find("NuEEL") != std::string::npos);  // nu e elastic
    const bool is_imd = (reaction.find("IMD") != std::string::npos);   // inverse mu decay
    const bool is_mec = (reaction.find("MEC") != std::string::npos);   // meson exchange current

    if (is_qel)
    {
        code = 0;
        if (is_cc)
            code = 1001;
        else if (is_nc)
            code = 1002;
    }
    else if (is_dis)
    {
        code = 2;
        if (is_cc)
            code = 1091;
        else if (is_nc)
            code = 1092;
    }
    else if (is_res)
        code = 1;
    else if (is_cohpi)
    {
        code = 3;
        if (is_qel)
            code = 4;
    }
    else if (is_ve)
        code = 1098;
    else if (is_imd)
        code = 1099;
    else if (is_mec)
        code = 10;

    // std::cout << "Reaction " << reaction << " has code = " << code << std::endl;

    return code;
}

std::string GetNuanceReaction(const int ccnc, const int mode)
{
    std::string reaction("");

    if (mode == 0)
    {
        reaction = "QES";
    }
    else if (mode == 2)
    {
        reaction = "DIS";
    }
    else if (mode == 1)
    {
        reaction = "RES";
    }
    else if (mode == 3)
    {
        reaction = "COH";
    }
    else if (mode == 5)
    {
        reaction = "NuEEL";
    }
    else if (mode == 6)
    {
        reaction = "IMD";
    }
    else if (mode == 10)
    {
        reaction = "MEC";
    }

    if (ccnc == 0)
    {
        reaction += "Weak[CC]";
    }
    else if (ccnc == 1)
    {
        reaction += "Weak[NC]";
    }

    return reaction;
}

//------------------------------------------------------------------------------------------------------------------------------------------

LArGrid MakeVoxelisationGrid(const pandora::Pandora *const pPrimaryPandora, const Parameters &parameters)
{
    // Detector volume for voxelising the hits
    const GeometryManager *geom = pPrimaryPandora->GetGeometry();
    const LArTPC &tpc = geom->GetLArTPC();

    const float botX = tpc.GetCenterX() - 0.5 * tpc.GetWidthX();
    const float botY = tpc.GetCenterY() - 0.5 * tpc.GetWidthY();
    const float botZ = tpc.GetCenterZ() - 0.5 * tpc.GetWidthZ();
    const float topX = botX + tpc.GetWidthX();
    const float topY = botY + tpc.GetWidthY();
    const float topZ = botZ + tpc.GetWidthZ();
    const float voxelWidth(parameters.m_voxelWidth);

    return LArGrid(pandora::CartesianVector(botX, botY, botZ), pandora::CartesianVector(topX, topY, topZ),
        pandora::CartesianVector(voxelWidth, voxelWidth, voxelWidth));
}

//------------------------------------------------------------------------------------------------------------------------------------------

LArGrid MakeVoxelisationGrid(const LArNDGeomSimple &geom, const Parameters &parameters)
{
    double minX{0.f};
    double maxX{0.f};
    double minY{0.f};
    double maxY{0.f};
    double minZ{0.f};
    double maxZ{0.f};
    const float voxelWidth(parameters.m_voxelWidth);
    geom.GetSurroundingBox(minX, maxX, minY, maxY, minZ, maxZ);

    return LArGrid(pandora::CartesianVector(minX, minY, minZ), pandora::CartesianVector(maxX, maxY, maxZ),
        pandora::CartesianVector(voxelWidth, voxelWidth, voxelWidth));
}

//------------------------------------------------------------------------------------------------------------------------------------------

LArVoxelList MakeVoxels(const LArHitInfo &hitInfo, const LArGrid &grid, const Parameters &parameters, const LArNDGeomSimple &geom)
{
    // Code based on
    // https://github.com/chenel/larcv2/tree/edepsim-formattruth/larcv/app/Supera/Voxelize.cxx
    // which is made available under the MIT license (which is fully compatible with Pandora's GPLv3 license)

    LArVoxelList currentVoxelList;

    // Start and end positions
    const pandora::CartesianVector start(hitInfo.m_start);
    const pandora::CartesianVector stop(hitInfo.m_stop);

    // Direction vector and hit segment length
    const pandora::CartesianVector dir = stop - start;
    const float hitLength(dir.GetMagnitude());

    // Check hit length is greater than epsilon limit
    if (hitLength < std::numeric_limits<float>::epsilon())
        return currentVoxelList;

    // Hit segment total energy in GeV (Geant4 uses MeV)
    const float g4HitEnergy(hitInfo.m_energy);

    // Check hit energy is greater than epsilon limit
    if (g4HitEnergy < std::numeric_limits<float>::epsilon())
        return currentVoxelList;

    // Get the trackID of the (main) contributing particle.
    // ATTN: this can very rarely be more than one track
    const int trackID = hitInfo.m_trackID;

    // Define ray trajectory, which checks dirMag (hitLength) >= epsilon limit
    const pandora::CartesianVector dirNorm = dir.GetUnitVector();
    LArRay ray(start, dirNorm);

    // We need to shuffle along the hit segment path and create voxels as we go.
    // There are 4 cases for the start and end points inside the voxelisation region.
    // Case 1: start & stop are both inside the voxelisation boundary
    // Case 2: start & stop are both outside, but path direction intersects boundary
    // Case 3: start is inside boundary, stop = intersection at region boundary
    // Case 4: end is inside boundary, start = intersection at region boundary

    double t0(0.0), t1(0.0);
    pandora::CartesianVector point1(0.f, 0.f, 0.f), point2(0.f, 0.f, 0.f);

    // Check if the start and end points are inside the voxelisation region
    const bool inStart = grid.Inside(start);
    const bool inStop = grid.Inside(stop);

    if (inStart && inStop)
    {
        // Case 1: Start and end points are inside boundary
        point1 = start;
        point2 = stop;
    }
    else if (!inStart && !inStop)
    {
        // Case 2: Start and end points are outside boundary
        if (grid.Intersect(ray, t0, t1))
        {
            point1 = ray.GetPoint(t0);
            point2 = ray.GetPoint(t1);
        }
        else
            return currentVoxelList;
    }
    else if (inStart && !inStop)
    {
        // Case 3: Start inside boundary
        point1 = start;
        if (grid.Intersect(ray, t0, t1))
            point2 = ray.GetPoint(t1);
        else
            return currentVoxelList;
    }
    else if (!inStart && inStop)
    {
        // Case 4: End inside boundary
        point2 = stop;
        if (grid.Intersect(ray, t0, t1))
            point1 = ray.GetPoint(t0);
        else
            return currentVoxelList;
    }

    // Now create voxels between point1 and point2.
    // Ray direction will be the same, but update starting point
    ray.UpdateOrigin(point1);

    bool shuffle(true);

    // Keep track of total voxel path length so far
    float totalPath(0.0);
    int loop(0);

    while (shuffle)
    {
        // Get point along path to define voxel bin (bottom corner)
        const pandora::CartesianVector voxelPoint = ray.GetPoint(parameters.m_voxelPathShift);

        // Grid 3d bin containing this point; 4th element is the total bin number
        const LongBin4Array gridBins = grid.GetBinIndices(voxelPoint);
        const long voxelID = gridBins[3];
        const long xBin = gridBins[0];
        const long yBin = gridBins[1];
        const long zBin = gridBins[2];

        // Voxel bottom and top corners
        const pandora::CartesianVector voxBot = grid.GetPoint(xBin, yBin, zBin);
        const pandora::CartesianVector voxTop = grid.GetPoint(xBin + 1, yBin + 1, zBin + 1);

        // Voxel box
        const LArBox vBox(voxBot, voxTop);

        // Get ray intersections with this box: t0 and t1 are set as the start
        // and end intersection pathlengths relative to the current ray point.
        // If we can't find t0 and t1, then stop shuffling along the path
        if (!vBox.Intersect(ray, t0, t1))
            shuffle = false;

        // Voxel extent = intersection path difference
        double dL(t1 - t0);
        // For the first path length, use the distance from the
        // starting ray point to the 2nd intersection t1
        if (loop == 0)
            dL = t1;

        // Stop processing if we are not moving along the path
        if (dL < parameters.m_voxelPathShift)
            shuffle = false;

        totalPath += dL;

        // Stop adding voxels if we have enough
        if (totalPath > hitLength)
        {
            shuffle = false;
            // Adjust final path according to hit segment total length
            dL = hitLength - totalPath + dL;
        }

        // Voxel energy (GeV) using path length fraction w.r.t hit length.
        // Here, hitLength is guaranteed to be greater than zero
        const float voxelEnergy(g4HitEnergy * dL / hitLength);

        if (parameters.m_useModularGeometry)
        {
            // If using modular geometry we need to assign the tpc number
            const int tpcID(geom.GetTPCNumber(voxelPoint));
            if (tpcID != -1)
            {
                const LArVoxel voxel(voxelID, voxelEnergy, voxBot, trackID, tpcID);
                currentVoxelList.emplace_back(voxel);
            }
            else
                std::cout << "Hit not in TPC: " << voxelPoint << std::endl;
        }
        else
        {
            const LArVoxel voxel(voxelID, voxelEnergy, voxBot, trackID);
            currentVoxelList.emplace_back(voxel);
        }

        // Update ray starting position using intersection path difference
        const pandora::CartesianVector newStart = ray.GetPoint(dL);
        ray.UpdateOrigin(newStart);
        loop++;
    }

    return currentVoxelList;
}

//------------------------------------------------------------------------------------------------------------------------------------------

LArVoxelList MergeSameVoxels(const LArVoxelList &voxelList)
{
    std::cout << "Merging voxels with the same IDs" << std::endl;
    LArVoxelList mergedVoxels;

    const int nVoxels = voxelList.size();
    std::vector<bool> processed(nVoxels, false);

    for (int i = 0; i < nVoxels; i++)
    {
        // Skip voxel if it was already used in a merge
        if (processed[i])
            continue;

        LArVoxel voxel1 = voxelList[i];
        float voxE1 = voxel1.m_energyInVoxel;
        std::map<int, float> trackIDToEnergy;
        trackIDToEnergy[voxel1.m_trackID] = voxE1;

        // Loop over other voxels (from i+1) and check if we have an ID match.
        // If so, add their energies and only store the combined voxel at the end
        for (int j = i + 1; j < nVoxels; j++)
        {
            // Skip voxel if it was already used in a merge
            if (processed[j])
                continue;

            const LArVoxel voxel2 = voxelList[j];
            const int trackid2 = voxel2.m_trackID;
            const float voxE2 = voxel2.m_energyInVoxel;
            if (voxel2.m_voxelID == voxel1.m_voxelID)
            {
                // IDs match. Add energy and set processed integer
                voxE1 += voxE2;
                processed[j] = true;
                // Amend the true particle contribution map
                if (trackIDToEnergy.count(trackid2) != 0)
                    trackIDToEnergy[trackid2] += voxE2;
                else
                    trackIDToEnergy[trackid2] = voxE2;
            }
        }

        // Add combined (or untouched) voxel to the merged list
        voxel1.SetEnergy(voxE1);
        // Update the track ID if necessary
        if (trackIDToEnergy.size() > 1)
        {
            float highestEnergy{0.f};
            int bestTrackID{-1};
            // std::cout << "Merged voxel had contributions from " << trackIDToEnergy.size() << " particles" << std::endl;
            for (auto const &pair : trackIDToEnergy)
            {
                // std::cout << " - " << pair.first << " with energy " << pair.second << std::endl;
                if (pair.second > highestEnergy)
                {
                    highestEnergy = pair.second;
                    bestTrackID = pair.first;
                }
            }
            // std::cout << " = chose track id " << bestTrackID << std::endl;
            voxel1.SetTrackID(bestTrackID);
        }

        mergedVoxels.emplace_back(voxel1);

        // We have processed the ith voxel
        processed[i] = true;
    }

    return mergedVoxels;
}

//------------------------------------------------------------------------------------------------------------------------------------------

LArVoxelProjectionList MergeSameProjections(const LArVoxelProjectionList &hits)
{
    LArVoxelProjectionList outputHits;
    std::vector<bool> areUsed(hits.size(), false);

    for (unsigned int vp1 = 0; vp1 < hits.size(); ++vp1)
    {
        if (areUsed.at(vp1))
            continue;

        LArVoxelProjection voxProj1{hits.at(vp1)};
        std::map<int, float> trackIDToEnergy;
        trackIDToEnergy[voxProj1.m_trackID] = voxProj1.m_energy;
        for (unsigned int vp2 = vp1 + 1; vp2 < hits.size(); ++vp2)
        {
            if (areUsed.at(vp2))
                continue;

            const LArVoxelProjection &voxProj2 = hits.at(vp2);
            if ((voxProj1.m_wire != voxProj2.m_wire) || (voxProj1.m_drift != voxProj2.m_drift))
                continue;

            // Add the energy, but keep track of the highest energy contributor
            voxProj1.m_energy += voxProj2.m_energy;
            if (trackIDToEnergy.count(voxProj2.m_trackID) != 0)
                trackIDToEnergy[voxProj2.m_trackID] += voxProj2.m_energy;
            else
                trackIDToEnergy[voxProj2.m_trackID] = voxProj2.m_energy;

            areUsed.at(vp2) = true;
        }
        // Add the hit to the output
        areUsed.at(vp1) = true;

        // Update the track ID if necessary
        if (trackIDToEnergy.size() > 1)
        {
            float highestEnergy{0.f};
            int bestTrackID{-1};
            for (auto const &pair : trackIDToEnergy)
            {
                if (pair.second > highestEnergy)
                {
                    highestEnergy = pair.second;
                    bestTrackID = pair.first;
                }
            }
            voxProj1.m_trackID = bestTrackID;
            voxProj1.m_parentVoxelID = hits.at(bestTrackID).m_parentVoxelID;
        }
        outputHits.emplace_back(voxProj1);
    }

    std::cout << outputHits.size() << " projected hits remain after merging" << std::endl;
    return outputHits;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void MakeCaloHitsFromVoxels(const LArVoxelList &voxels, const MCParticleEnergyMap &mcEnergyMap,
    const pandora::Pandora *const pPrimaryPandora, const Parameters &parameters, int &hitCounter)
{

    // Factory for creating LArCaloHits
    lar_content::LArCaloHitFactory m_larCaloHitFactory;
    const float voxelWidth(parameters.m_voxelWidth);
    const float MipE = 0.00075;
    lar_content::LArCaloHitParameters caloHitParameters = MakeDefaultCaloHitParams(voxelWidth);

    if (parameters.m_use3D)
    {
        for (unsigned int v = 0; v < voxels.size(); ++v)
        {
            const LArVoxel voxel = voxels.at(v);
            const pandora::CartesianVector voxelPos(voxel.m_voxelPosVect);
            const float voxelE = voxel.m_energyInVoxel;
            const float voxelMipEquivalentE = voxelE / MipE;

            if (voxelMipEquivalentE < parameters.m_minVoxelMipEquivE)
                continue;

            // Modify the important fields
            caloHitParameters.m_positionVector = voxelPos;
            caloHitParameters.m_inputEnergy = voxelE;
            caloHitParameters.m_mipEquivalentEnergy = voxelMipEquivalentE;
            caloHitParameters.m_electromagneticEnergy = voxelE;
            caloHitParameters.m_hadronicEnergy = voxelE;
            caloHitParameters.m_pParentAddress = (void *)(static_cast<uintptr_t>(++hitCounter));
            caloHitParameters.m_larTPCVolumeId = voxel.m_tpcID;

            PANDORA_THROW_RESULT_IF(
                pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::CaloHit::Create(*pPrimaryPandora, caloHitParameters, m_larCaloHitFactory));

            // Set calo hit voxel to MCParticle relation using trackID
            const int trackID = voxel.m_trackID;
            const float energyFrac = GetMCEnergyFraction(mcEnergyMap, voxelE, trackID);
            PandoraApi::SetCaloHitToMCParticleRelationship(*pPrimaryPandora, (void *)((intptr_t)hitCounter), (void *)((intptr_t)trackID), energyFrac);
        }
    }

    // Treat the 3 x 2D view case separately as we need to merge hits
    // on some of the views depending on the geometry
    if (parameters.m_useLArTPC)
    {
        LArVoxelProjectionList voxelProjectionsU;
        LArVoxelProjectionList voxelProjectionsV;
        LArVoxelProjectionList voxelProjectionsW;

        for (unsigned int v = 0; v < voxels.size(); ++v)
        {
            const LArVoxel &voxel = voxels.at(v);

            const pandora::CartesianVector voxelPos = voxel.m_voxelPosVect;
            const float uPos(pPrimaryPandora->GetPlugins()->GetLArTransformationPlugin()->YZtoU(voxelPos.GetY(), voxelPos.GetZ()));
            voxelProjectionsU.emplace_back(LArVoxelProjection(
                voxel.m_energyInVoxel, uPos, voxelPos.GetX(), pandora::TPC_VIEW_U, voxel.m_voxelID, voxel.m_trackID, voxel.m_tpcID));

            const float vPos(pPrimaryPandora->GetPlugins()->GetLArTransformationPlugin()->YZtoV(voxelPos.GetY(), voxelPos.GetZ()));
            voxelProjectionsV.emplace_back(LArVoxelProjection(
                voxel.m_energyInVoxel, vPos, voxelPos.GetX(), pandora::TPC_VIEW_V, voxel.m_voxelID, voxel.m_trackID, voxel.m_tpcID));

            const float wPos(pPrimaryPandora->GetPlugins()->GetLArTransformationPlugin()->YZtoW(voxelPos.GetY(), voxelPos.GetZ()));
            voxelProjectionsW.emplace_back(LArVoxelProjection(
                voxel.m_energyInVoxel, wPos, voxelPos.GetX(), pandora::TPC_VIEW_W, voxel.m_voxelID, voxel.m_trackID, voxel.m_tpcID));
        }

        std::vector<LArVoxelProjectionList> viewProjections;
        viewProjections.emplace_back(MergeSameProjections(voxelProjectionsU));
        viewProjections.emplace_back(MergeSameProjections(voxelProjectionsV));
        viewProjections.emplace_back(MergeSameProjections(voxelProjectionsW));

        voxelProjectionsU.clear();
        voxelProjectionsV.clear();
        voxelProjectionsW.clear();

        for (const LArVoxelProjectionList &view : viewProjections)
        {
            for (const LArVoxelProjection &hit : view)
            {
                const float voxelE = hit.m_energy;
                const float voxelMipEquivalentE = voxelE / MipE;

                if (voxelMipEquivalentE < parameters.m_minVoxelMipEquivE)
                    continue;

                // Modify the important fields
                caloHitParameters.m_positionVector = pandora::CartesianVector(hit.m_drift, 0.f, hit.m_wire);
                caloHitParameters.m_inputEnergy = voxelE;
                caloHitParameters.m_mipEquivalentEnergy = voxelMipEquivalentE;
                caloHitParameters.m_electromagneticEnergy = voxelE;
                caloHitParameters.m_hadronicEnergy = voxelE;
                caloHitParameters.m_pParentAddress = (void *)(static_cast<uintptr_t>(++hitCounter));
                caloHitParameters.m_hitType = hit.m_view;
                caloHitParameters.m_larTPCVolumeId = hit.m_tpcID;

                // Create LArCaloHits for U, V and W views
                PANDORA_THROW_RESULT_IF(
                    pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::CaloHit::Create(*pPrimaryPandora, caloHitParameters, m_larCaloHitFactory));

                // Set calo hit voxel to MCParticle relation using trackID
                const int trackID = hit.m_trackID;
                const float energyFrac = GetMCEnergyFraction(mcEnergyMap, voxelE, trackID);
                PandoraApi::SetCaloHitToMCParticleRelationship(*pPrimaryPandora, (void *)((intptr_t)hitCounter), (void *)((intptr_t)trackID), energyFrac);
            } // end voxel projection loop
        } // end view loop
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

float GetMCEnergyFraction(const MCParticleEnergyMap &mcEnergyMap, const float voxelE, const int trackID)
{
    // Find the energy fraction: voxelHitE/MCParticleE
    float energyFrac(0.f), MCEnergy(0.f);
    MCParticleEnergyMap::const_iterator mapIter = mcEnergyMap.find(trackID);
    if (mapIter != mcEnergyMap.end())
        MCEnergy = mapIter->second;

    if (MCEnergy > 0.0)
        energyFrac = voxelE / MCEnergy;

    return energyFrac;
}

//------------------------------------------------------------------------------------------------------------------------------------------

lar_content::LArCaloHitParameters MakeDefaultCaloHitParams(float voxelWidth)
{
    lar_content::LArCaloHitParameters caloHitParameters;
    caloHitParameters.m_positionVector = pandora::CartesianVector(0.f, 0.f, 1.f);
    caloHitParameters.m_expectedDirection = pandora::CartesianVector(0.f, 0.f, 1.f);
    caloHitParameters.m_cellNormalVector = pandora::CartesianVector(0.f, 0.f, 1.f);
    caloHitParameters.m_cellGeometry = pandora::RECTANGULAR;
    caloHitParameters.m_cellSize0 = voxelWidth;
    caloHitParameters.m_cellSize1 = voxelWidth;
    caloHitParameters.m_cellThickness = voxelWidth;
    caloHitParameters.m_nCellRadiationLengths = 1.f;
    caloHitParameters.m_nCellInteractionLengths = 1.f;
    caloHitParameters.m_time = -1.f;
    caloHitParameters.m_inputEnergy = 0.f;
    caloHitParameters.m_mipEquivalentEnergy = 0.f;
    caloHitParameters.m_electromagneticEnergy = 0.f;
    caloHitParameters.m_hadronicEnergy = 0.f;
    caloHitParameters.m_isDigital = false;
    caloHitParameters.m_hitType = pandora::TPC_3D;
    caloHitParameters.m_hitRegion = pandora::SINGLE_REGION;
    caloHitParameters.m_layer = 0;
    caloHitParameters.m_isInOuterSamplingLayer = false;
    caloHitParameters.m_pParentAddress = (void *)(static_cast<uintptr_t>(0));
    caloHitParameters.m_larTPCVolumeId = 0;
    caloHitParameters.m_daughterVolumeId = 0;
    return caloHitParameters;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool ParseCommandLine(int argc, char *argv[], Parameters &parameters)
{
    if (1 == argc)
        return PrintOptions();

    int cOpt(0);

    std::string recoOption("");
    std::string viewOption("both");
    std::string formatOption("SP");
    std::string inputTreeName("");
    std::string geomManagerName("");
    std::string geomVolName("");
    std::string sensDetName("");

    while ((cOpt = getopt(argc, argv, "r:i:e:k:f:g:t:v:d:n:s:j:w:m:b:c:MpNh")) != -1)
    {
        switch (cOpt)
        {
            case 'r':
                recoOption = optarg;
                break;
            case 'i':
                parameters.m_settingsFile = optarg;
                break;
            case 'e':
                parameters.m_inputFileName = optarg;
                break;
            case 'k':
                inputTreeName = optarg;
                break;
            case 'f':
                formatOption = optarg;
                break;
            case 'g':
                parameters.m_geomFileName = optarg;
                break;
            case 't':
                geomManagerName = optarg;
                break;
            case 'v':
                geomVolName = optarg;
                break;
            case 'd':
                sensDetName = optarg;
                break;
            case 'M':
                parameters.m_useModularGeometry = true;
                break;
            case 'n':
                parameters.m_nEventsToProcess = atoi(optarg);
                break;
            case 's':
                parameters.m_nEventsToSkip = atoi(optarg);
                break;
            case 'p':
                parameters.m_printOverallRecoStatus = true;
                break;
            case 'j':
                viewOption = optarg;
                break;
            case 'w':
                parameters.m_voxelWidth = atof(optarg);
                break;
            case 'm':
                parameters.m_maxMergedVoxels = atoi(optarg);
                break;
            case 'b':
                parameters.m_minNSpacePoints = atoi(optarg);
                break;
            case 'c':
                parameters.m_minVoxelMipEquivE = atof(optarg);
                break;
            case 'N':
                parameters.m_shouldDisplayEventNumber = true;
                break;
            case 'h':
            default:
                return PrintOptions();
        }
    }

    ProcessViewOption(viewOption, parameters);
    const bool gotFormat = ProcessFormatOption(formatOption, inputTreeName, geomManagerName, geomVolName, sensDetName, parameters);
    const bool gotRecoOpt = ProcessRecoOption(recoOption, parameters);
    const bool passed = gotFormat && gotRecoOpt;
    if (!passed)
    {
        return PrintOptions();
    }
    return passed;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool PrintOptions()
{
    std::cout << std::endl
              << "./bin/PandoraInterface " << std::endl
              << "    -r RecoOption          (required) [Full, AllHitsCR, AllHitsNu, CRRemHitsSliceCR, CRRemHitsSliceNu, AllHitsSliceCR, AllHitsSliceNu]"
              << std::endl
              << "    -i Settings            (required) [Run xml file for setting up the Pandora algorithms]" << std::endl
              << "    -e EventsFile          (required) [Events input data ROOT file]" << std::endl
              << "    -g GeometryFile        (required) [ROOT file containing the TGeoManager geometry]" << std::endl
              << "    -f DataFormat          (optional) [SP (SpacePoint default), SPMC (SpacePoint MC), EDepSim (rooTracker) or SED (LArSoft-like)]"
              << std::endl
              << "    -k EventsTreeName      (optional) [Name of the input events ROOT TTree (default = events)]" << std::endl
              << "    -t TGeoManagerName     (optional) [TGeoManager name (default = Default)]" << std::endl
              << "    -v geometryVolName     (optional) [ND LAr physical volume name (default = volArgonCubeCryostat_PV)]" << std::endl
              << "    -d sensitiveDetName    (optional) [ND LAr sensitive detector name (default = volTPCActive)]" << std::endl
              << "    -M                     (optional) [Use modular geometry that makes each TPC active volume separately (default = false)]"
              << std::endl
              << "    -j Projection          (optional) [Both (default), 3D or LArTPC (2D projections only)]" << std::endl
              << "    -n NEventsToProcess    (optional) [Number of events to process]" << std::endl
              << "    -s NEventsToSkip       (optional) [Number of events to skip in the event input file]" << std::endl
              << "    -p                     (optional) [Print status]" << std::endl
              << "    -N                     (optional) [Print event numbers]" << std::endl
              << "    -w width               (optional) [Voxel bin width (cm), default = 0.4 cm]" << std::endl
              << "    -m maxMergedVoxels     (optional) [Skip events that have N(space points) or N(merged voxels) > maxMergedVoxels (default = no events skipped)]"
              << std::endl
              << "    -b minNSpacePoints     (optional) [Skip events that have N(space points) < minNSpacePoints (default < 2)]" << std::endl
              << "    -c minMipEquivE        (optional) [Minimum MIP equivalent energy, default = 0.3]" << std::endl
              << std::endl;

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ProcessViewOption(const std::string &viewOption, Parameters &parameters)
{
    std::string chosenViewOption(viewOption);
    std::transform(chosenViewOption.begin(), chosenViewOption.end(), chosenViewOption.begin(), ::tolower);

    if (chosenViewOption == "both")
    {
        // Both LArTPC 2D projections and 3D hits (default)
        std::cout << "Using LArTPC projections _and_ 3D hits" << std::endl;
        parameters.m_useLArTPC = true;
        parameters.m_use3D = true;
    }
    else if (chosenViewOption == "3d")
    {
        // 3D hits only
        std::cout << "Using 3D hits" << std::endl;
        parameters.m_useLArTPC = false;
        parameters.m_use3D = true;
    }
    else if (chosenViewOption == "lartpc")
    {
        // LArTPC 2D projected hits only
        std::cout << "Using LArTPC projections" << std::endl;
        parameters.m_useLArTPC = true;
        parameters.m_use3D = false;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool ProcessRecoOption(const std::string &recoOption, Parameters &parameters)
{
    std::string chosenRecoOption(recoOption);
    std::transform(chosenRecoOption.begin(), chosenRecoOption.end(), chosenRecoOption.begin(), ::tolower);
    bool processed(true);

    if ("full" == chosenRecoOption)
    {
        parameters.m_shouldRunAllHitsCosmicReco = true;
        parameters.m_shouldRunStitching = true;
        parameters.m_shouldRunCosmicHitRemoval = true;
        parameters.m_shouldRunSlicing = true;
        parameters.m_shouldRunNeutrinoRecoOption = true;
        parameters.m_shouldRunCosmicRecoOption = true;
        parameters.m_shouldPerformSliceId = true;
    }
    else if ("allhitscr" == chosenRecoOption)
    {
        parameters.m_shouldRunAllHitsCosmicReco = true;
        parameters.m_shouldRunStitching = true;
        parameters.m_shouldRunCosmicHitRemoval = false;
        parameters.m_shouldRunSlicing = false;
        parameters.m_shouldRunNeutrinoRecoOption = false;
        parameters.m_shouldRunCosmicRecoOption = false;
        parameters.m_shouldPerformSliceId = false;
    }
    else if ("nostitchingcr" == chosenRecoOption)
    {
        parameters.m_shouldRunAllHitsCosmicReco = false;
        parameters.m_shouldRunStitching = false;
        parameters.m_shouldRunCosmicHitRemoval = false;
        parameters.m_shouldRunSlicing = false;
        parameters.m_shouldRunNeutrinoRecoOption = false;
        parameters.m_shouldRunCosmicRecoOption = true;
        parameters.m_shouldPerformSliceId = false;
    }
    else if ("allhitsnu" == chosenRecoOption)
    {
        parameters.m_shouldRunAllHitsCosmicReco = false;
        parameters.m_shouldRunStitching = false;
        parameters.m_shouldRunCosmicHitRemoval = false;
        parameters.m_shouldRunSlicing = false;
        parameters.m_shouldRunNeutrinoRecoOption = true;
        parameters.m_shouldRunCosmicRecoOption = false;
        parameters.m_shouldPerformSliceId = false;
    }
    else if ("crremhitsslicecr" == chosenRecoOption)
    {
        parameters.m_shouldRunAllHitsCosmicReco = true;
        parameters.m_shouldRunStitching = true;
        parameters.m_shouldRunCosmicHitRemoval = true;
        parameters.m_shouldRunSlicing = true;
        parameters.m_shouldRunNeutrinoRecoOption = false;
        parameters.m_shouldRunCosmicRecoOption = true;
        parameters.m_shouldPerformSliceId = false;
    }
    else if ("crremhitsslicenu" == chosenRecoOption)
    {
        parameters.m_shouldRunAllHitsCosmicReco = true;
        parameters.m_shouldRunStitching = true;
        parameters.m_shouldRunCosmicHitRemoval = true;
        parameters.m_shouldRunSlicing = true;
        parameters.m_shouldRunNeutrinoRecoOption = true;
        parameters.m_shouldRunCosmicRecoOption = false;
        parameters.m_shouldPerformSliceId = false;
    }
    else if ("allhitsslicecr" == chosenRecoOption)
    {
        parameters.m_shouldRunAllHitsCosmicReco = false;
        parameters.m_shouldRunStitching = false;
        parameters.m_shouldRunCosmicHitRemoval = false;
        parameters.m_shouldRunSlicing = true;
        parameters.m_shouldRunNeutrinoRecoOption = false;
        parameters.m_shouldRunCosmicRecoOption = true;
        parameters.m_shouldPerformSliceId = false;
    }
    else if ("allhitsslicenu" == chosenRecoOption)
    {
        parameters.m_shouldRunAllHitsCosmicReco = false;
        parameters.m_shouldRunStitching = false;
        parameters.m_shouldRunCosmicHitRemoval = false;
        parameters.m_shouldRunSlicing = true;
        parameters.m_shouldRunNeutrinoRecoOption = true;
        parameters.m_shouldRunCosmicRecoOption = false;
        parameters.m_shouldPerformSliceId = false;
    }
    else
    {
        std::cout << "Unrecognized reconstruction option: " << recoOption << std::endl;
        processed = false;
    }

    return processed;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool ProcessFormatOption(const std::string &formatOption, const std::string &inputTreeName, const std::string &geomManagerName,
    const std::string &geomVolName, const std::string &sensDetName, Parameters &parameters)
{
    std::string chosenFormatOption(formatOption);
    std::transform(chosenFormatOption.begin(), chosenFormatOption.end(), chosenFormatOption.begin(), ::tolower);
    bool processed(true);

    if (chosenFormatOption == "sp" || chosenFormatOption == "spmc")
    {
        // Space point ROOT format (data = default or MC)
        parameters.m_dataFormat = Parameters::LArNDFormat::SP;
        if (chosenFormatOption == "spmc")
            parameters.m_dataFormat = Parameters::LArNDFormat::SPMC;
        // Set the event input TTree name
        parameters.m_inputTreeName = inputTreeName.empty() ? "events" : inputTreeName;
        // Set the TGeoManager name
        parameters.m_geomManagerName = geomManagerName.empty() ? "Default" : geomManagerName;
        // Set geometry volume name
        parameters.m_geometryVolName = geomVolName.empty() ? "volArgonCubeCryostat_PV" : geomVolName;
        // Set the sensitive detector name
        parameters.m_sensitiveDetName = sensDetName.empty() ? "volTPCActive" : sensDetName;
        // All lengths are already in cm, so don't rescale
        parameters.m_lengthScale = 1.0f;
        // All energies are already in GeV, so don't rescale
        parameters.m_energyScale = 1.0f;
    }
    else if (chosenFormatOption == "edepsim")
    {
        // Assume EDepSim rooTracker format
        parameters.m_dataFormat = Parameters::LArNDFormat::EDepSim;
        // Set the event input TTree name
        parameters.m_inputTreeName = inputTreeName.empty() ? "EDepSimEvents" : inputTreeName;
        // Set the TGeoManager name
        parameters.m_geomManagerName = geomManagerName.empty() ? "EDepSimGeometry" : geomManagerName;
        // Set geometry volume name
        parameters.m_geometryVolName = geomVolName.empty() ? "volArgonCubeDetector_PV_0" : geomVolName;
        // Set the sensitive detector name
        parameters.m_sensitiveDetName = sensDetName.empty() ? "ArgonCube" : sensDetName;
        // All lengths are in mm, so we need to convert them to cm
        parameters.m_lengthScale = parameters.m_mm2cm;
        // All energies are in MeV, so we need to convert them to GeV
        parameters.m_energyScale = parameters.m_MeV2GeV;
    }
    else if (chosenFormatOption == "sed")
    {
        // LArSoft-type SimEnergyDeposit (SED) ROOT format
        parameters.m_dataFormat = Parameters::LArNDFormat::SED;
        // Set the event input TTree name
        parameters.m_inputTreeName = inputTreeName.empty() ? "simdump/ndsim" : inputTreeName;
        // Set the TGeoManager name
        parameters.m_geomManagerName = geomManagerName.empty() ? "Default" : geomManagerName;
        // Set geometry volume name
        parameters.m_geometryVolName = geomVolName.empty() ? "volArgonCubeCryostat_PV" : geomVolName;
        // Set the sensitive detector name if not set
        parameters.m_sensitiveDetName = sensDetName.empty() ? "volTPCActive" : sensDetName;
        // All lengths are already in cm, so don't rescale
        parameters.m_lengthScale = 1.0f;
        // All energies are already in GeV, so don't rescale
        parameters.m_energyScale = 1.0f;
    }
    else
    {
        std::cout << "Unrecognized data format option: " << formatOption << std::endl;
        processed = false;
    }
    return processed;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ProcessExternalParameters(const Parameters &parameters, const Pandora *const pPandora)
{
    auto *const pEventSteeringParameters = new lar_content::MasterAlgorithm::ExternalSteeringParameters;
    pEventSteeringParameters->m_shouldRunAllHitsCosmicReco = parameters.m_shouldRunAllHitsCosmicReco;
    pEventSteeringParameters->m_shouldRunStitching = parameters.m_shouldRunStitching;
    pEventSteeringParameters->m_shouldRunCosmicHitRemoval = parameters.m_shouldRunCosmicHitRemoval;
    pEventSteeringParameters->m_shouldRunSlicing = parameters.m_shouldRunSlicing;
    pEventSteeringParameters->m_shouldRunNeutrinoRecoOption = parameters.m_shouldRunNeutrinoRecoOption;
    pEventSteeringParameters->m_shouldRunCosmicRecoOption = parameters.m_shouldRunCosmicRecoOption;
    pEventSteeringParameters->m_shouldPerformSliceId = parameters.m_shouldPerformSliceId;
    pEventSteeringParameters->m_printOverallRecoStatus = parameters.m_printOverallRecoStatus;

    // LArMaster or LArMasterThreeD algorithms
    if (!parameters.m_use3D)
    {
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::SetExternalParameters(*pPandora, "LArMaster", pEventSteeringParameters));
    }
    else
    {
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::SetExternalParameters(*pPandora, "LArMasterThreeD", pEventSteeringParameters));
    }

#ifdef LIBTORCH_DL
    auto *const pEventSettingsParametersCopy = new lar_content::MasterAlgorithm::ExternalSteeringParameters(*pEventSteeringParameters);
    PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=,
        pandora::ExternallyConfiguredAlgorithm::SetExternalParameters(*pPandora, "LArDLMaster", pEventSettingsParametersCopy));
#endif
}

} // namespace lar_nd_reco
