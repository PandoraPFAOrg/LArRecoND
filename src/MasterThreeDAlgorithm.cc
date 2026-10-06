/**
 *  @file   src/MasterThreeDAlgorithm.cc
 *
 *  @brief  Implementation of the 3D master algorithm class.
 *
 *  $Log: $
 */

#include "Api/PandoraApi.h"

#include "Pandora/AlgorithmHeaders.h"

#include "LArNDContent.h"
#include "MasterThreeDAlgorithm.h"

#include "larpandoracontent/LArContent.h"
#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArFileHelper.h"
#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"
#include "larpandoracontent/LArHelpers/LArStitchingHelper.h"

#include "larpandoracontent/LArObjects/LArCaloHit.h"
#include "larpandoracontent/LArObjects/LArMCParticle.h"

#include "larpandoracontent/LArPlugins/LArPseudoLayerPlugin.h"
#include "larpandoracontent/LArPlugins/LArRotationalTransformationPlugin.h"

#include "larpandoracontent/LArUtility/PfoMopUpBaseAlgorithm.h"

#ifdef LIBTORCH_DL
#include "larpandoradlcontent/LArDLContent.h"
#endif

using namespace pandora;

namespace lar_content
{

MasterThreeDAlgorithm::MasterThreeDAlgorithm() :
    m_shouldRunRockMus_Xworkers(false),
    m_tagRockMuons(false)
{
}

StatusCode MasterThreeDAlgorithm::Run()
{

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->Reset());

    if (!m_workerInstancesInitialized)
    {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->InitializeWorkerInstances(m_workerToLArTPCMap));
    }

    if (m_passMCParticlesToWorkerInstances)
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->CopyMCParticles());

    PfoToFloatMap stitchedPfosToX0Map;
    VolumeIdToHitListMap volumeIdToHitListMap;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->GetVolumeIdToHitListMap(volumeIdToHitListMap));

    if (m_shouldRunAllHitsCosmicReco)
    {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->RunCosmicRayReconstruction(volumeIdToHitListMap));

        PfoToLArTPCMap pfoToLArTPCMap;
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->RecreateCosmicRayPfos(pfoToLArTPCMap));

        if (m_shouldRunStitching)
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->StitchCosmicRayPfos(pfoToLArTPCMap, stitchedPfosToX0Map));
    }

    if (m_shouldRunCosmicHitRemoval)
    {
        PfoList clearCosmicRayPfos, ambiguousPfos;
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->TagCosmicRayPfos(stitchedPfosToX0Map, clearCosmicRayPfos, ambiguousPfos));

        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->RunCosmicRayHitRemoval(ambiguousPfos));
    }

    SliceVector sliceVector;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->RunSlicing(volumeIdToHitListMap, sliceVector));

    if (m_shouldRunNeutrinoRecoOption || m_shouldRunCosmicRecoOption)
    {
        SliceHypotheses nuSliceHypotheses, crSliceHypotheses;
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->RunSliceReconstruction(sliceVector, nuSliceHypotheses, crSliceHypotheses));
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->SelectBestSliceHypotheses(nuSliceHypotheses, crSliceHypotheses));
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------
StatusCode MasterThreeDAlgorithm::TagCosmicRayPfos(const PfoToFloatMap &stitchedPfosToX0Map, PfoList &clearCosmicRayPfos, PfoList &ambiguousPfos) const
{

    PfoList ambiguousPfos_wRock;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, MasterAlgorithm::TagCosmicRayPfos(stitchedPfosToX0Map, clearCosmicRayPfos, ambiguousPfos_wRock));

    if (!m_tagRockMuons)
    {
        for (const Pfo *const pPfo : ambiguousPfos_wRock)
            ambiguousPfos.push_back(pPfo);

        return STATUS_CODE_SUCCESS;
    }

    for (RockMuonTaggingTool *const pRockMuonTaggingTool : m_rockMuonTaggingToolVector)
        pRockMuonTaggingTool->FindAmbiguousPfos(ambiguousPfos_wRock, ambiguousPfos, this);

    // filter ambiguousPfos form rocks

    const PfoList *pRecreatedCRPfos(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::GetCurrentPfoList(this->GetPandora(), pRecreatedCRPfos));

    for (const Pfo *const pPfo : *pRecreatedCRPfos)
    {
        const bool isClearRock(ambiguousPfos.end() == std::find(ambiguousPfos.begin(), ambiguousPfos.end(), pPfo));
        PandoraContentApi::ParticleFlowObject::Metadata metadata;
        metadata.m_propertiesToAdd["IsClearCosmic"] = (isClearRock ? 1.f : 0.f); // TODO maybe decouple labels? isClearRock, good for now
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ParticleFlowObject::AlterMetadata(*this, pPfo, metadata));

        if (isClearRock)
            clearCosmicRayPfos.push_back(pPfo);
    }

    if (m_visualizeOverallRecoStatus)
    {
        PANDORA_MONITORING_API(VisualizeParticleFlowObjects(this->GetPandora(), &clearCosmicRayPfos, "ClearCRPfos", RED));
        PANDORA_MONITORING_API(VisualizeParticleFlowObjects(this->GetPandora(), &ambiguousPfos, "AmbiguousPfos", BLUE));
        PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
    }

    return STATUS_CODE_SUCCESS;
}
//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode MasterThreeDAlgorithm::RunCosmicRayReconstruction(const VolumeIdToHitListMap &volumeIdToHitListMap) const
{
    for (const Pandora *const pCRWorker : m_crWorkerInstances)
    {
        const LArTPC &worker_larTPC(pCRWorker->GetGeometry()->GetLArTPC());
        const unsigned int worker_id = worker_larTPC.GetLArTPCVolumeId();

        const auto worker_it = m_workerToLArTPCMap.find(worker_id);
        if (worker_it == m_workerToLArTPCMap.end())
        {
            std::cout << "Problem: worker_id " << worker_id << "not present in the m_workerToLArTPCMap. Skipping this worker.\n";
            continue;
        }

        // loop over worker's TPCs
        for (const pandora::LArTPC *pLArTPC : worker_it->second)
        {
            const unsigned int larTPC_id = (*pLArTPC).GetLArTPCVolumeId();

            // get all TPC's hits
            VolumeIdToHitListMap::const_iterator iter(volumeIdToHitListMap.find(larTPC_id));

            if (volumeIdToHitListMap.end() == iter)
                continue;

            // copy hits into the worker
            std::cout << "Copying " << iter->second.m_allHitList.size() << " hits from LArTPC " << larTPC_id << " to worker " << worker_id << "\n";
            for (const CaloHit *const pCaloHit : iter->second.m_allHitList)
                PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->Copy(pCRWorker, pCaloHit));
        }

        if (m_printOverallRecoStatus)
            std::cout << "Running cosmic-ray reconstruction worker instance " << worker_id << " of " << m_crWorkerInstances.size() << std::endl;

        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::ProcessEvent(*pCRWorker));
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------
StatusCode MasterThreeDAlgorithm::RunCosmicRayHitRemoval(const PfoList &ambiguousPfos) const
{
    PfoList allPfosToDelete;
    LArPfoHelper::GetAllConnectedPfos(ambiguousPfos, allPfosToDelete);

    for (const Pfo *const pPfoToDelete : allPfosToDelete)
    {
        const ClusterList clusterList(pPfoToDelete->GetClusterList());
        const VertexList vertexList(pPfoToDelete->GetVertexList());

        // ATTN: If an ambiguous pfo has been stitched, reset the calo hit positions in preparation for subsequent algorithm chains
        if (LArStitchingHelper::HasPfoBeenStitched(pPfoToDelete))
        {
            CaloHitList caloHitList2D;
            LArPfoHelper::GetCaloHits(pPfoToDelete, TPC_VIEW_U, caloHitList2D);
            LArPfoHelper::GetCaloHits(pPfoToDelete, TPC_VIEW_V, caloHitList2D);
            LArPfoHelper::GetCaloHits(pPfoToDelete, TPC_VIEW_W, caloHitList2D);
            LArPfoHelper::GetCaloHits(pPfoToDelete, TPC_3D, caloHitList2D);
            LArPfoHelper::GetIsolatedCaloHits(pPfoToDelete, TPC_VIEW_U, caloHitList2D);
            LArPfoHelper::GetIsolatedCaloHits(pPfoToDelete, TPC_VIEW_V, caloHitList2D);
            LArPfoHelper::GetIsolatedCaloHits(pPfoToDelete, TPC_VIEW_W, caloHitList2D);
            LArPfoHelper::GetIsolatedCaloHits(pPfoToDelete, TPC_3D, caloHitList2D);

            for (const CaloHit *const pCaloHit : caloHitList2D)
            {
                PandoraContentApi::CaloHit::Metadata metadata;
                metadata.m_x0 = 0.f;
                PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::CaloHit::AlterMetadata(*this, pCaloHit, metadata));
            }
        }

        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Delete(*this, pPfoToDelete));
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Delete(*this, &clusterList));
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Delete(*this, &vertexList));
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode MasterThreeDAlgorithm::RunSlicing(const VolumeIdToHitListMap &volumeIdToHitListMap, SliceVector &sliceVector) const
{
    std::cout << "There are " << volumeIdToHitListMap.size() << " volumes" << std::endl;
    for (const VolumeIdToHitListMap::value_type &mapEntry : volumeIdToHitListMap)
    {
        std::cout << "- Volume has " << mapEntry.second.m_allHitList.size() << " hits" << std::endl;
        for (const CaloHit *const pCaloHit : (m_shouldRemoveOutOfTimeHits ? mapEntry.second.m_truncatedHitList : mapEntry.second.m_allHitList))
        {
            if (!PandoraContentApi::IsAvailable(*this, pCaloHit))
                continue;

            if (m_shouldRunSlicing)
            {
                PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->Copy(m_pSlicingWorkerInstance, pCaloHit));
            }
            else
            {
                if (sliceVector.empty())
                    sliceVector.push_back(CaloHitList());
                sliceVector.back().push_back(pCaloHit);
            }
        }
    }

    if (m_shouldRunSlicing)
    {
        if (m_printOverallRecoStatus)
            std::cout << "Running slicing worker instance" << std::endl;

        const PfoList *pSlicePfos(nullptr);
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::ProcessEvent(*m_pSlicingWorkerInstance));
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::GetCurrentPfoList(*m_pSlicingWorkerInstance, pSlicePfos));

        if (m_visualizeOverallRecoStatus)
        {
            PANDORA_MONITORING_API(VisualizeParticleFlowObjects(this->GetPandora(), pSlicePfos, "OnePfoPerSlice", BLUE));
            PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
        }

        for (const Pfo *const pSlicePfo : *pSlicePfos)
        {
            sliceVector.push_back(CaloHitList());
            LArPfoHelper::GetCaloHits(pSlicePfo, TPC_VIEW_U, sliceVector.back());
            LArPfoHelper::GetCaloHits(pSlicePfo, TPC_VIEW_V, sliceVector.back());
            LArPfoHelper::GetCaloHits(pSlicePfo, TPC_VIEW_W, sliceVector.back());
            LArPfoHelper::GetCaloHits(pSlicePfo, TPC_3D, sliceVector.back());
        }
    }

    if (m_printOverallRecoStatus)
        std::cout << "Identified " << sliceVector.size() << " slice(s)" << std::endl;

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode MasterThreeDAlgorithm::Recreate(const ParticleFlowObject *const pInputPfo, const ParticleFlowObject *const pNewParentPfo, PfoList &newPfoList) const
{
    ClusterList inputClusterList2D, inputClusterList3D, totalClusterList, newClusterList;

    // 2D clusters
    LArPfoHelper::GetTwoDClusterList(pInputPfo, inputClusterList2D);
    totalClusterList.insert(totalClusterList.end(), inputClusterList2D.begin(), inputClusterList2D.end());

    // 3D clusters
    LArPfoHelper::GetThreeDClusterList(pInputPfo, inputClusterList3D);
    totalClusterList.insert(totalClusterList.end(), inputClusterList3D.begin(), inputClusterList3D.end());

    for (const Cluster *const pInputCluster : totalClusterList)
    {
        CaloHitList inputCaloHitList, newCaloHitList, newIsolatedCaloHitList;
        pInputCluster->GetOrderedCaloHitList().FillCaloHitList(inputCaloHitList);

        for (const CaloHit *const pInputCaloHit : inputCaloHitList)
            newCaloHitList.push_back(static_cast<const CaloHit *>(pInputCaloHit->GetParentAddress()));

        for (const CaloHit *const pInputCaloHit : pInputCluster->GetIsolatedCaloHitList())
            newIsolatedCaloHitList.push_back(static_cast<const CaloHit *>(pInputCaloHit->GetParentAddress()));

        if (!newCaloHitList.empty())
            newClusterList.push_back(this->CreateCluster(pInputCluster, newCaloHitList, newIsolatedCaloHitList));
    }

    VertexList newVertexList;

    for (const Vertex *const pInputVertex : pInputPfo->GetVertexList())
        newVertexList.push_back(this->CreateVertex(pInputVertex));

    const ParticleFlowObject *const pNewPfo(this->CreatePfo(pInputPfo, newClusterList, newVertexList));
    newPfoList.push_back(pNewPfo);

    if (pNewParentPfo)
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SetPfoParentDaughterRelationship(*this, pNewParentPfo, pNewPfo))

    for (const ParticleFlowObject *const pInputDaughterPfo : pInputPfo->GetDaughterPfoList())
        this->Recreate(pInputDaughterPfo, pNewPfo, newPfoList);

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

const Pandora *MasterThreeDAlgorithm::CreateWorkerInstance(
    const LArTPC &larTPC, const DetectorGapList &gapList, const std::string &settingsFile, const std::string &name) const
{
    // The Pandora instance
    const Pandora *const pPandora(new Pandora(name));
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, LArContent::RegisterAlgorithms(*pPandora));
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, LArContent::RegisterBasicPlugins(*pPandora));
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, LArNDContent::RegisterAlgorithms(*pPandora));
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::SetPseudoLayerPlugin(*pPandora, new lar_content::LArPseudoLayerPlugin));
    PANDORA_THROW_RESULT_IF(
        STATUS_CODE_SUCCESS, !=, PandoraApi::SetLArTransformationPlugin(*pPandora, new lar_content::LArRotationalTransformationPlugin));
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->RegisterCustomContent(pPandora));
    MultiPandoraApi::AddDaughterPandoraInstance(&(this->GetPandora()), pPandora);

    // The LArTPC
    PandoraApi::Geometry::LArTPC::Parameters larTPCParameters;
    larTPCParameters.m_larTPCVolumeId = larTPC.GetLArTPCVolumeId();
    larTPCParameters.m_centerX = larTPC.GetCenterX();
    larTPCParameters.m_centerY = larTPC.GetCenterY();
    larTPCParameters.m_centerZ = larTPC.GetCenterZ();
    larTPCParameters.m_widthX = larTPC.GetWidthX();
    larTPCParameters.m_widthY = larTPC.GetWidthY();
    larTPCParameters.m_widthZ = larTPC.GetWidthZ();
    larTPCParameters.m_wirePitchU = larTPC.GetWirePitchU();
    larTPCParameters.m_wirePitchV = larTPC.GetWirePitchV();
    larTPCParameters.m_wirePitchW = larTPC.GetWirePitchW();
    larTPCParameters.m_wireAngleU = larTPC.GetWireAngleU();
    larTPCParameters.m_wireAngleV = larTPC.GetWireAngleV();
    larTPCParameters.m_wireAngleW = larTPC.GetWireAngleW();
    larTPCParameters.m_sigmaUVW = larTPC.GetSigmaUVW();
    larTPCParameters.m_isDriftInPositiveX = larTPC.IsDriftInPositiveX();
    this->AppendReadoutVolumeParameters(larTPC, 0, larTPCParameters.m_readoutVolumeParametersVector);
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::Geometry::LArTPC::Create(*pPandora, larTPCParameters));

    // The Gaps
    for (const DetectorGap *const pGap : gapList)
    {
        const BoxGap *const pBoxGap( dynamic_cast<const BoxGap *>(pGap) );

        if ((pBoxGap) &&
            (pBoxGap->GetSide1().GetX() < larTPC.GetWidthX()) &&
            (pBoxGap->GetSide2().GetY() <= larTPC.GetWidthY()) &&
            (pBoxGap->GetSide3().GetZ() < larTPC.GetWidthZ()))
        {
            PandoraApi::Geometry::BoxGap::Parameters boxGapParameters;
            boxGapParameters.m_side1 = pBoxGap->GetSide1();
            boxGapParameters.m_side2 = pBoxGap->GetSide2();
            boxGapParameters.m_side3 = pBoxGap->GetSide3();
            boxGapParameters.m_vertex = pBoxGap->GetVertex();

            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::Geometry::BoxGap::Create(*pPandora, boxGapParameters));
        }
    }

    // Configuration
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::ReadSettings(*pPandora, settingsFile));
    return pPandora;
}

//------------------------------------------------------------------------------------------------------------------------------------------

const Pandora *MasterThreeDAlgorithm::CreateWorkerInstance(const LArTPCMap &larTPCMap, const DetectorGapList &gapList,
    const std::string &settingsFile, const std::string &name, const unsigned int id = 0) const
{
    if (larTPCMap.empty())
    {
        std::cout << "MasterThreeDAlgorithm::CreateWorkerInstance - no LArTPC details provided" << std::endl;
        throw StatusCodeException(STATUS_CODE_NOT_INITIALIZED);
    }

    // The Pandora instance
    const Pandora *const pPandora(new Pandora(name));
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, LArContent::RegisterAlgorithms(*pPandora));
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, LArContent::RegisterBasicPlugins(*pPandora));
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, LArNDContent::RegisterAlgorithms(*pPandora));
    // Deep Learning algorithms (e.g. vertexing)
#ifdef LIBTORCH_DL
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, LArDLContent::RegisterAlgorithms(*pPandora));
#endif
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::SetPseudoLayerPlugin(*pPandora, new lar_content::LArPseudoLayerPlugin));
    PANDORA_THROW_RESULT_IF(
        STATUS_CODE_SUCCESS, !=, PandoraApi::SetLArTransformationPlugin(*pPandora, new lar_content::LArRotationalTransformationPlugin));
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->RegisterCustomContent(pPandora));
    MultiPandoraApi::AddDaughterPandoraInstance(&(this->GetPandora()), pPandora);

    // The Parent LArTPC
    const LArTPC *const pFirstLArTPC(larTPCMap.begin()->second);
    float parentMinX(pFirstLArTPC->GetCenterX() - 0.5f * pFirstLArTPC->GetWidthX());
    float parentMaxX(pFirstLArTPC->GetCenterX() + 0.5f * pFirstLArTPC->GetWidthX());
    float parentMinY(pFirstLArTPC->GetCenterY() - 0.5f * pFirstLArTPC->GetWidthY());
    float parentMaxY(pFirstLArTPC->GetCenterY() + 0.5f * pFirstLArTPC->GetWidthY());
    float parentMinZ(pFirstLArTPC->GetCenterZ() - 0.5f * pFirstLArTPC->GetWidthZ());
    float parentMaxZ(pFirstLArTPC->GetCenterZ() + 0.5f * pFirstLArTPC->GetWidthZ());

    std::vector<int> lof_tpcs = {};
    for (const LArTPCMap::value_type &mapEntry : larTPCMap)
    {
        const LArTPC *const pLArTPC(mapEntry.second);
        lof_tpcs.push_back(pLArTPC->GetLArTPCVolumeId());
        parentMinX = std::min(parentMinX, pLArTPC->GetCenterX() - 0.5f * pLArTPC->GetWidthX());
        parentMaxX = std::max(parentMaxX, pLArTPC->GetCenterX() + 0.5f * pLArTPC->GetWidthX());
        parentMinY = std::min(parentMinY, pLArTPC->GetCenterY() - 0.5f * pLArTPC->GetWidthY());
        parentMaxY = std::max(parentMaxY, pLArTPC->GetCenterY() + 0.5f * pLArTPC->GetWidthY());
        parentMinZ = std::min(parentMinZ, pLArTPC->GetCenterZ() - 0.5f * pLArTPC->GetWidthZ());
        parentMaxZ = std::max(parentMaxZ, pLArTPC->GetCenterZ() + 0.5f * pLArTPC->GetWidthZ());
    }

    std::cout << "Creating worker instance " << name << " with boundaries X = (" << parentMinX << ", " << parentMaxX << ") , Y = ("
              << parentMinY << ", " << parentMaxY << ") , Z = (" << parentMinZ << ", " << parentMaxZ << ")";
    std::cout << " worker's tpcs : ";
    for (const int tpc_id : lof_tpcs)
        std::cout << tpc_id << ",";
    std::cout << "\n";

    PandoraApi::Geometry::LArTPC::Parameters larTPCParameters;
    larTPCParameters.m_larTPCVolumeId = id;
    larTPCParameters.m_centerX = 0.5f * (parentMaxX + parentMinX);
    larTPCParameters.m_centerY = 0.5f * (parentMaxY + parentMinY);
    larTPCParameters.m_centerZ = 0.5f * (parentMaxZ + parentMinZ);
    larTPCParameters.m_widthX = parentMaxX - parentMinX;
    larTPCParameters.m_widthY = parentMaxY - parentMinY;
    larTPCParameters.m_widthZ = parentMaxZ - parentMinZ;
    larTPCParameters.m_wirePitchU = std::max(pFirstLArTPC->GetWirePitchU(), pFirstLArTPC->GetWirePitchV());
    larTPCParameters.m_wirePitchV = std::max(pFirstLArTPC->GetWirePitchU(), pFirstLArTPC->GetWirePitchV());
    larTPCParameters.m_wirePitchW = pFirstLArTPC->GetWirePitchW();
    larTPCParameters.m_wireAngleU = pFirstLArTPC->GetWireAngleU();
    larTPCParameters.m_wireAngleV = pFirstLArTPC->GetWireAngleV();
    larTPCParameters.m_wireAngleW = pFirstLArTPC->GetWireAngleW();
    larTPCParameters.m_sigmaUVW = pFirstLArTPC->GetSigmaUVW();
    larTPCParameters.m_isDriftInPositiveX = pFirstLArTPC->IsDriftInPositiveX();

    // Merge in readout volumes from every child TPC, keeping ids unique, and record the offset used for each original LArTPC so hits copied
    // into this worker can be remapped consistently in Copy()
    m_daughterVolumeIdOffsetMap.clear();
    unsigned int idOffset(0);
    for (const LArTPCMap::value_type &mapEntry : larTPCMap)
    {
        const LArTPC *const pLArTPC(mapEntry.second);
        // ATTN: This gets rebuilt on each call at startup (slicing), but it's always consistent, so it's fine to do this.
        m_daughterVolumeIdOffsetMap[pLArTPC->GetLArTPCVolumeId()] = idOffset;
        this->AppendReadoutVolumeParameters(*pLArTPC, idOffset, larTPCParameters.m_readoutVolumeParametersVector);
        idOffset += pLArTPC->GetReadoutVolumes().size();
    }

    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::Geometry::LArTPC::Create(*pPandora, larTPCParameters));

    // The Gaps
    for (const DetectorGap *const pGap : gapList)
    {
        const BoxGap *const pBoxGap( dynamic_cast<const BoxGap *>(pGap) );

        if (pBoxGap)
        {
            PandoraApi::Geometry::BoxGap::Parameters boxGapParameters;
            boxGapParameters.m_side1 = pBoxGap->GetSide1();
            boxGapParameters.m_side2 = pBoxGap->GetSide2();
            boxGapParameters.m_side3 = pBoxGap->GetSide3();
            boxGapParameters.m_vertex = pBoxGap->GetVertex();

            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::Geometry::BoxGap::Create(*pPandora, boxGapParameters));
        }
    }

    // Configuration
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::ReadSettings(*pPandora, settingsFile));
    return pPandora;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode MasterThreeDAlgorithm::InitializeWorkerInstances(WorkerToLArTPCMap &workerToLArTPCMap)
{
    // ATTN Used to be in the regular Initialize callback, but detector gap list cannot be extracted in client app before the first event
    if (m_workerInstancesInitialized)
        return STATUS_CODE_ALREADY_INITIALIZED;

    try
    {
        const LArTPCMap &larTPCMap(this->GetPandora().GetGeometry()->GetLArTPCMap());
        const DetectorGapList &gapList(this->GetPandora().GetGeometry()->GetDetectorGapList());

        if (m_shouldRunRockMus_Xworkers)
        {
            // Group TPCs (proxy for drift volumes here) that share the same XY coordinates
            // in a unique "columnar" worker instance.
            // In the next lines, a map with the following structure is created:
            // (x0, y0) <---> (0: drift_volume_z0, 1: drift_volume_z1, ...)
            // (x1, y1) <---> (0: drift_volume_z0, 1: drift_volume_z1, ...)
            // where (x, y) represents the shared coordinates of drift volumes in the same column.

            std::map<std::pair<float, float>, LArTPCMap> XYgrouped;
            const unsigned int FIRST_TPC_ID = 0;
            for (const LArTPCMap::value_type &mapEntry : larTPCMap)
            {
                const LArTPC &tpc(*(mapEntry.second));
                const auto key_xy = std::make_pair(tpc.GetCenterX(), tpc.GetCenterY());
                auto it_tpcMap = XYgrouped.find(key_xy);

                if (it_tpcMap != XYgrouped.end())
                {
                    // get the current group of tpcs correspondint to key_xy
                    LArTPCMap &tpcMap = it_tpcMap->second;
                    // get the tpc id of the most recent added tpc - the newest tpc
                    // id is the most recent + 1
                    const unsigned int current_max_id = tpcMap.empty() ? FIRST_TPC_ID : tpcMap.rbegin()->first + 1;
                    XYgrouped[key_xy].emplace(current_max_id, &tpc);
                }
                else // new column of tpcs
                {
                    XYgrouped[key_xy].emplace(FIRST_TPC_ID, &tpc);
                }
            }
            // now that we have grouped drift volumes along xy in a map
            // create a worker instance for each group
            unsigned int worker_id = 0;
            for (const auto &[xy, submap] : XYgrouped)
            {
                // const auto& [x, y] = xy;
                m_crWorkerInstances.push_back(
                    this->CreateWorkerInstance(submap, gapList, m_crSettingsFile, "CRWorkerInstance" + std::to_string(worker_id), worker_id));

                // loop over the group of TPCs along the same XY
                // and fill the map worker <---> vector<TPCs>
                for (const LArTPCMap::value_type &mapEntry : submap)
                    workerToLArTPCMap[worker_id].push_back(mapEntry.second);

                worker_id++;
            }
        }
        else
        { // dafault: crate 1 worker instance per drift volume
            for (const LArTPCMap::value_type &mapEntry : larTPCMap)
            {
                const unsigned int volumeId(mapEntry.second->GetLArTPCVolumeId());
                m_crWorkerInstances.push_back(
                    this->CreateWorkerInstance(*(mapEntry.second), gapList, m_crSettingsFile, "CRWorkerInstance" + std::to_string(volumeId)));

                workerToLArTPCMap[volumeId].push_back(mapEntry.second);
            }
        }

        if (m_shouldRunSlicing)
            m_pSlicingWorkerInstance = this->CreateWorkerInstance(larTPCMap, gapList, m_slicingSettingsFile, "SlicingWorker");

        if (m_shouldRunNeutrinoRecoOption)
            m_pSliceNuWorkerInstance = this->CreateWorkerInstance(larTPCMap, gapList, m_nuSettingsFile, "SliceNuWorker");

        if (m_shouldRunCosmicRecoOption)
            m_pSliceCRWorkerInstance = this->CreateWorkerInstance(larTPCMap, gapList, m_crSettingsFile, "SliceCRWorker");
    }
    catch (const StatusCodeException &statusCodeException)
    {
        std::cout << "MasterAlgorithm: Exception during initialization of worker instances " << statusCodeException.ToString() << std::endl;
        return statusCodeException.GetStatusCode();
    }

    m_workerInstancesInitialized = true;
    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode MasterThreeDAlgorithm::GetVolumeIdToHitListMap(VolumeIdToHitListMap &volumeIdToHitListMap) const
{
    const LArTPCMap &larTPCMap(this->GetPandora().GetGeometry()->GetLArTPCMap());
    const unsigned int nLArTPCs(larTPCMap.size());

    const CaloHitList *pCaloHitList(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_inputHitListName, pCaloHitList));

    for (const CaloHit *const pCaloHit : *pCaloHitList)
    {
        const LArCaloHit *const pLArCaloHit(dynamic_cast<const LArCaloHit *>(pCaloHit));

        if (!pLArCaloHit && (1 != nLArTPCs))
            return STATUS_CODE_INVALID_PARAMETER;

        const unsigned int volumeId(pLArCaloHit ? pLArCaloHit->GetLArTPCVolumeId() : 0);
        const LArTPC *const pLArTPC(larTPCMap.at(volumeId));

        LArTPCHitList &larTPCHitList(volumeIdToHitListMap[volumeId]);
        larTPCHitList.m_allHitList.push_back(pCaloHit);

        if (((pCaloHit->GetPositionVector().GetX() >= (pLArTPC->GetCenterX() - 0.5f * pLArTPC->GetWidthX())) &&
                (pCaloHit->GetPositionVector().GetX() <= (pLArTPC->GetCenterX() + 0.5f * pLArTPC->GetWidthX()))))
        {
            larTPCHitList.m_truncatedHitList.push_back(pCaloHit);
        }
        else
            std::cout << "Hit of type " << pCaloHit->GetHitType() << " outside TPC " << volumeId << "? "
                      << pCaloHit->GetPositionVector().GetX() << ", " << pLArTPC->GetCenterX() - 0.5f * pLArTPC->GetWidthX() << ", "
                      << pLArTPC->GetCenterX() + 0.5f * pLArTPC->GetWidthX() << std::endl;
    }

    return STATUS_CODE_SUCCESS;
}

StatusCode MasterThreeDAlgorithm::ReadSettings(const pandora::TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "ShouldRunRockMus_Xworkers", m_shouldRunRockMus_Xworkers));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "TagRockMuons", m_tagRockMuons));

    if (m_shouldRunCosmicHitRemoval)
    {
        AlgorithmToolVector algorithmToolVector;
        PANDORA_RETURN_RESULT_IF(
            STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithmToolList(*this, xmlHandle, "RockMuonTaggingTools", algorithmToolVector));

        for (AlgorithmTool *const pAlgorithmTool : algorithmToolVector)
        {
            RockMuonTaggingTool *const pRockMuonTaggingTool(dynamic_cast<RockMuonTaggingTool *>(pAlgorithmTool));
            if (!pRockMuonTaggingTool)
                return STATUS_CODE_INVALID_PARAMETER;
            m_rockMuonTaggingToolVector.push_back(pRockMuonTaggingTool);
        }
    }

    return MasterAlgorithm::ReadSettings(xmlHandle);
}

} // namespace lar_content
