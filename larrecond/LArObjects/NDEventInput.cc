/**
 *  @file   larrecond/LArObjects/NDEventInput.cc
 *
 *  @brief  Implementation of the ND event input
 *
 *  $Log: $
 */

#include <array>
#include <utility>

#include "Api/PandoraApi.h"

#include "Managers/PluginManager.h"
#include "Pandora/PandoraSettings.h"

#include "larpandoracontent/LArControlFlow/MultiPandoraApi.h"
#include "larpandoracontent/LArPlugins/LArRotationalTransformationPlugin.h"

#include "larrecond/LArObjects/NDEventInput.h"

#include "TFile.h"
#include "TTree.h"

namespace
{
// Mip-equivalent energy of a voxel (GeV)
constexpr float MipE{0.00075f};
} // namespace

namespace lar_nd_reco
{

NDEventInput::NDEventInput() :
    m_pPandora{nullptr},
    m_pMainPandora{nullptr},
    m_parameters(),
    m_tpcGeom(),
    m_pFile{nullptr},
    m_pTree{nullptr},
    m_nEntries{0},
    m_LArCaloHitFactory()
{
}

NDEventInput::NDEventInput(const pandora::Pandora *pPandora, const NDParameters &parameters, const NDSimpleTPCGeom &tpcGeom) :
    m_pPandora{pPandora},
    m_pMainPandora{nullptr},
    m_parameters(parameters),
    m_tpcGeom(tpcGeom),
    m_pFile{nullptr},
    m_pTree{nullptr},
    m_nEntries{0},
    m_LArCaloHitFactory()
{
    this->Initialize();
}

void NDEventInput::Initialize()
{
    if (this->ShouldDisplayInfo())
    {
        std::cout << "Constructing NDEventInput for " << m_pPandora->GetName() << ", " << m_parameters.m_inputFileName << ", "
                  << m_parameters.m_inputTreeName << std::endl;
    }

    // Get the main Pandora instance pointer
    m_pMainPandora = MultiPandoraApi::GetPrimaryPandoraInstance(m_pPandora);
    if (this->ShouldDisplayInfo())
        std::cout << "Main Pandora instance is " << m_pMainPandora->GetName() << std::endl;

    m_pFile = TFile::Open(m_parameters.m_inputFileName.c_str(), "READ");
    if (m_pFile)
    {
        m_pTree = dynamic_cast<TTree *>(m_pFile->Get(m_parameters.m_inputTreeName.c_str()));
        if (!m_pTree)
        {
            std::cout << "Error in NDEventInput: can't find tree " << m_parameters.m_inputTreeName << std::endl;
            m_pFile->Close();
            delete m_pFile;
            m_pFile = nullptr;
            throw pandora::StatusCodeException(pandora::STATUS_CODE_NOT_FOUND);
        }
        else
            m_nEntries = m_pTree->GetEntries();
    }
    else
    {
        std::cout << "Error in NDEventInput: can't open file " << m_parameters.m_inputFileName << std::endl;
        throw pandora::StatusCodeException(pandora::STATUS_CODE_NOT_FOUND);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

NDEventInput::~NDEventInput()
{
    if (m_pFile)
    {
        m_pFile->Close();
        delete m_pFile;
    }
}

bool NDEventInput::ShouldDisplayInfo() const
{
    return m_pPandora && m_pPandora->GetSettings() && m_pPandora->GetSettings()->ShouldDisplayAlgorithmInfo();
}

//------------------------------------------------------------------------------------------------------------------------------------------

lar_content::LArCaloHitParameters NDEventInput::MakeCaloHitParameters(const pandora::CartesianVector &voxelPos, const float voxelE, const int tpcID) const
{
    const float voxelWidth(m_parameters.m_voxelWidth);

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
    caloHitParameters.m_time = 0.f;
    caloHitParameters.m_inputEnergy = voxelE;
    caloHitParameters.m_mipEquivalentEnergy = voxelE / MipE;
    caloHitParameters.m_electromagneticEnergy = voxelE;
    caloHitParameters.m_hadronicEnergy = voxelE;
    caloHitParameters.m_isDigital = false;
    caloHitParameters.m_hitType = pandora::TPC_3D;
    caloHitParameters.m_hitRegion = pandora::SINGLE_REGION;
    caloHitParameters.m_layer = 0;
    caloHitParameters.m_isInOuterSamplingLayer = false;
    caloHitParameters.m_larTPCVolumeId = tpcID < 0 ? 0 : tpcID;
    caloHitParameters.m_daughterVolumeId = 0;

    return caloHitParameters;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void NDEventInput::CreateVoxelCaloHits(const lar_content::LArCaloHitParameters &caloHitParameters, const pandora::CartesianVector &voxelPos,
    int &hitCounter, const MCLink *pMCLink) const
{
    // Create a CaloHit for the main Pandora instance and a copy for this instance, with an optional MC link
    const auto createCaloHit = [&](lar_content::LArCaloHitParameters hitParameters, const int hitID)
    {
        void *const pHitAddress = (void *)(static_cast<intptr_t>(hitID));
        hitParameters.m_pParentAddress = pHitAddress;

        PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::CaloHit::Create(*m_pMainPandora, hitParameters, m_LArCaloHitFactory));
        PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::CaloHit::Create(*m_pPandora, hitParameters, m_LArCaloHitFactory));

        if (pMCLink)
        {
            void *const pTrackAddress = (void *)(static_cast<intptr_t>(pMCLink->m_trackID));
            PandoraApi::SetCaloHitToMCParticleRelationship(*m_pMainPandora, pHitAddress, pTrackAddress, pMCLink->m_energyFrac);
            PandoraApi::SetCaloHitToMCParticleRelationship(*m_pPandora, pHitAddress, pTrackAddress, pMCLink->m_energyFrac);
        }
    };

    // ATTN: A hit ID is always reserved for the 3D hit, even if the 3D hit is not made
    const int hitID3D = ++hitCounter;

    if (m_parameters.m_use3D)
        createCaloHit(caloHitParameters, hitID3D);

    if (m_parameters.m_useLArTPC)
    {
        // Create LArCaloHits for U, V and W views assuming x is the common drift coordinate
        const float x0_cm(voxelPos.GetX());
        const float y0_cm(voxelPos.GetY());
        const float z0_cm(voxelPos.GetZ());

        const auto *const pTransformationPlugin(m_pPandora->GetPlugins()->GetLArTransformationPlugin());

        const std::array<std::pair<pandora::HitType, float>, 3> views{{{pandora::TPC_VIEW_U, pTransformationPlugin->YZtoU(y0_cm, z0_cm)},
            {pandora::TPC_VIEW_V, pTransformationPlugin->YZtoV(y0_cm, z0_cm)}, {pandora::TPC_VIEW_W, pTransformationPlugin->YZtoW(y0_cm, z0_cm)}}};

        for (const auto &view : views)
        {
            lar_content::LArCaloHitParameters viewParameters(caloHitParameters);
            viewParameters.m_hitType = view.first;
            viewParameters.m_positionVector = pandora::CartesianVector(x0_cm, 0.f, view.second);
            createCaloHit(viewParameters, ++hitCounter);
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

} // namespace lar_nd_reco
