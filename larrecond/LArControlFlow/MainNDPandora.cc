/**
 *  @file   larrecond/LArControlFlow/MainNDPandora.cc
 *
 *  @brief  Implementation of the main Pandora ND instance
 *
 *  $Log: $
 */

#include "Api/PandoraApi.h"
#include "Pandora/PandoraSettings.h"
#include "Pandora/StatusCodes.h"
#include "PandoraMonitoringApi.h"

#include "larpandoracontent/LArContent.h"
#include "larpandoracontent/LArControlFlow/MasterAlgorithm.h"
#include "larpandoracontent/LArControlFlow/MultiPandoraApi.h"
#include "larpandoracontent/LArHelpers/LArFileHelper.h"
#include "larpandoracontent/LArPlugins/LArPseudoLayerPlugin.h"
#include "larpandoracontent/LArPlugins/LArRotationalTransformationPlugin.h"

#ifdef LIBTORCH_DL
#include "larpandoradlcontent/LArDLContent.h"
#endif

#include "larrecond/LArControlFlow/MainNDPandora.h"
#include "larrecond/LArHelpers/NDGeometryHelper.h"
#include "larrecond/LArNDContent.h"
#include "larrecond/LArObjects/LArNDEventInput.h"
#include "larrecond/LArObjects/TMSEventInput.h"

using namespace pandora;

namespace
{

// ATTN: The wire parameters below are placeholders for the ND (voxel-based) geometry,
// taken from the MicroBooNE LArTPC definition.
constexpr float WirePitch{0.300000011921f};  ///< Wire pitch for the U, V and W views
constexpr float WireAngleU{1.04719758034f};  ///< Wire angle for the U view (rad), +60 degrees
constexpr float WireAngleV{-1.04719758034f}; ///< Wire angle for the V view (rad), -60 degrees
constexpr float WireAngleW{0.f};             ///< Wire angle for the W view (rad)
constexpr float SigmaUVW{1.f};               ///< Wire position uncertainty

/**
 *  @brief  Make the Pandora LArTPC parameters from a simple ND TPC box
 *
 *  @param  tpc The simple TPC description
 *
 *  @return The Pandora LArTPC parameters
 */
PandoraApi::Geometry::LArTPC::Parameters CreateLArTPCParameters(const lar_nd_reco::NDSimpleTPC &tpc)
{
    PandoraApi::Geometry::LArTPC::Parameters geoParameters;

    geoParameters.m_centerX = 0.5 * (tpc.m_x_min + tpc.m_x_max);
    geoParameters.m_centerY = 0.5 * (tpc.m_y_min + tpc.m_y_max);
    geoParameters.m_centerZ = 0.5 * (tpc.m_z_min + tpc.m_z_max);
    geoParameters.m_widthX = (tpc.m_x_max - tpc.m_x_min);
    geoParameters.m_widthY = (tpc.m_y_max - tpc.m_y_min);
    geoParameters.m_widthZ = (tpc.m_z_max - tpc.m_z_min);

    geoParameters.m_larTPCVolumeId = tpc.m_TPC_ID;
    geoParameters.m_wirePitchU = WirePitch;
    geoParameters.m_wirePitchV = WirePitch;
    geoParameters.m_wirePitchW = WirePitch;
    geoParameters.m_wireAngleU = WireAngleU;
    geoParameters.m_wireAngleV = WireAngleV;
    geoParameters.m_wireAngleW = WireAngleW;
    geoParameters.m_sigmaUVW = SigmaUVW;

    // ATTN: Assumes the drift direction alternates with the TPC ID (even IDs drift in -x, odd IDs in +x)
    geoParameters.m_isDriftInPositiveX = (tpc.m_TPC_ID % 2) != 0;

    return geoParameters;
}

} // namespace

namespace lar_nd_reco
{

//------------------------------------------------------------------------------------------------------------------------------------------

MainNDPandora::MainNDPandora(const std::string &mainName, const NDParameters &mainParameters) :
    m_mainParameters(mainParameters),
    m_pMainPandora{nullptr},
    m_NDParametersMap{},
    m_NDSimpleTPCGeomMap{},
    m_instanceOrder{}
{
    m_pMainPandora = this->CreatePandoraInstance(mainName, mainParameters);
    MultiPandoraApi::AddPrimaryPandoraInstance(m_pMainPandora);
    MultiPandoraApi::SetVolumeId(m_pMainPandora, NDParameters::VolType::Main);
}

//------------------------------------------------------------------------------------------------------------------------------------------

MainNDPandora::~MainNDPandora()
{
    // Delete all Pandora instances, including the main one
    if (this->ShouldDisplayInfo())
        std::cout << "DeletePandoraInstances" << std::endl;
    MultiPandoraApi::DeletePandoraInstances(m_pMainPandora);

    for (const auto &pair : m_NDEventInputMap)
    {
        const NDEventInput *pEventInput = pair.second;
        delete pEventInput;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void MainNDPandora::AddPandoraInstance(const std::string &name, const NDParameters &parameters)
{
    const pandora::Pandora *pPandora = this->CreatePandoraInstance(name, parameters);

    // Make this a daughter of the main instance
    MultiPandoraApi::AddDaughterPandoraInstance(m_pMainPandora, pPandora);

    // Also add this as a primary instance for the MultiPandoraApi,
    // which is used to add daughter instances inside MasterThreeDAlgorithm etc.
    MultiPandoraApi::AddPrimaryPandoraInstance(pPandora);
    // Also set the volume Id (VolType enum)
    MultiPandoraApi::SetVolumeId(pPandora, parameters.m_volType);

    // Add parameters to the internal map
    m_NDParametersMap.insert(std::make_pair(name, parameters));

    // Keep track of the instance insertion order
    m_instanceOrder.emplace_back(name);
}

//------------------------------------------------------------------------------------------------------------------------------------------

const pandora::Pandora *MainNDPandora::CreatePandoraInstance(const std::string &name, const NDParameters &parameters) const
{
    std::cout << "Creating Pandora instance " << name << std::endl;

    const pandora::Pandora *pPandora = new pandora::Pandora(name);

    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, LArContent::RegisterAlgorithms(*pPandora));
#ifdef LIBTORCH_DL
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, LArDLContent::RegisterAlgorithms(*pPandora));
#endif
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, LArContent::RegisterBasicPlugins(*pPandora));

    if (parameters.m_use3D)
    {
        std::cout << "Registered 3D algorithms" << std::endl;
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, LArNDContent::RegisterAlgorithms(*pPandora));
    }

    return pPandora;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void MainNDPandora::CreatePandoraTPCs()
{
    std::cout << "CreatePandoraTPCs" << std::endl;

    // Keep track of the total number of TPCs being made, since we need
    // to use this to set unique TPC IDs for all Pandora instances
    int nTotalTPCs{0};

    // Loop over the Pandora instances
    const PandoraInstanceList pandoraList = MultiPandoraApi::GetDaughterPandoraInstanceList(m_pMainPandora);
    for (const auto *pPandora : pandoraList)
    {
        const std::string pandoraName = pPandora->GetName();
        const NDParameters &parameters = this->GetPandoraParameters(pandoraName);

        std::cout << "Creating TPCs for " << pandoraName << std::endl;

        const NDGeometryHelper geometry(parameters.m_geomFileName, parameters.m_geomManagerName);

        // TPC name from the parameters
        const std::string tpcName = parameters.m_tpcName;

        // Simple object storing the TPC boxes
        const NDSimpleTPCGeom tpcGeom = geometry.GetTPCGeometry(tpcName, nTotalTPCs);
        // Store this in the TPC geometry map (for processing the events)
        m_NDSimpleTPCGeomMap.insert(std::make_pair(pandoraName, tpcGeom));

        nTotalTPCs += tpcGeom.m_TPCs.size();
        std::cout << "nTotalTPCs = " << nTotalTPCs << std::endl;

        // Create the TPCs for the given & main Pandora instances
        std::map<unsigned int, NDSimpleTPC>::const_iterator tpcIter;
        for (tpcIter = tpcGeom.m_TPCs.begin(); tpcIter != tpcGeom.m_TPCs.end(); ++tpcIter)
        {
            const NDSimpleTPC &tpc = tpcIter->second;
            const PandoraApi::Geometry::LArTPC::Parameters geoParameters(CreateLArTPCParameters(tpc));

            try
            {
                // Create TPC for both pPandora and the main Pandora instance
                PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::Geometry::LArTPC::Create(*pPandora, geoParameters));
                PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::Geometry::LArTPC::Create(*m_pMainPandora, geoParameters));
                std::cout << "Created TPC " << tpc.m_TPC_ID << " for " << pPandora->GetName() << " & " << m_pMainPandora->GetName() << std::endl;
            }
            catch (const pandora::StatusCodeException &)
            {
                std::cout << "CreatePandoraTPCs - unable to create TPC, insufficient or invalid information supplied" << std::endl;
            }
        } // TPCs

        // TPC transformation plugins for this Pandora instance
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::SetPseudoLayerPlugin(*pPandora, new lar_content::LArPseudoLayerPlugin));
        PANDORA_THROW_RESULT_IF(
            STATUS_CODE_SUCCESS, !=, PandoraApi::SetLArTransformationPlugin(*pPandora, new lar_content::LArRotationalTransformationPlugin));
    } // Pandora instances

    // Main instance TPC transformation plugins
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::SetPseudoLayerPlugin(*m_pMainPandora, new lar_content::LArPseudoLayerPlugin));
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=,
        PandoraApi::SetLArTransformationPlugin(*m_pMainPandora, new lar_content::LArRotationalTransformationPlugin));
}

//------------------------------------------------------------------------------------------------------------------------------------------

void MainNDPandora::CreatePandoraDetectorGaps()
{
    // To be implemented
}

//------------------------------------------------------------------------------------------------------------------------------------------

void MainNDPandora::ConfigurePandoraInstances()
{
    std::cout << "ConfigurePandoraInstances" << std::endl;

    this->ProcessExternalParameters(m_pMainPandora, m_mainParameters);

    std::string settingsFilePath(lar_content::LArFileHelper::FindFileInPath(m_mainParameters.m_settingsFile, "FW_SEARCH_PATH"));
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::ReadSettings(*m_pMainPandora, settingsFilePath));

    const PandoraInstanceList pandoraList = MultiPandoraApi::GetDaughterPandoraInstanceList(m_pMainPandora);
    for (const auto *pPandora : pandoraList)
    {
        const NDParameters &parameters = this->GetPandoraParameters(pPandora->GetName());
        this->ProcessExternalParameters(pPandora, parameters);
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::ReadSettings(*pPandora, parameters.m_settingsFile));
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void MainNDPandora::ProcessExternalParameters(const pandora::Pandora *const pPandora, const NDParameters &parameters) const
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

    // LArMasterThreeD or LArMaster algorithms
    if (parameters.m_use3D)
    {
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::SetExternalParameters(*pPandora, "LArMasterThreeD", pEventSteeringParameters));
    }
    else
    {
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::SetExternalParameters(*pPandora, "LArMaster", pEventSteeringParameters));
    }

#ifdef LIBTORCH_DL
    auto *const pEventSettingsParametersCopy = new lar_content::MasterAlgorithm::ExternalSteeringParameters(*pEventSteeringParameters);
    PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=,
        pandora::ExternallyConfiguredAlgorithm::SetExternalParameters(*pPandora, "LArDLMaster", pEventSettingsParametersCopy));
#endif
}

//------------------------------------------------------------------------------------------------------------------------------------------

void MainNDPandora::ConfigureEventInputs()
{
    // Create NDEventInputs for each Pandora daughter instance
    const PandoraInstanceList pandoraList = MultiPandoraApi::GetDaughterPandoraInstanceList(m_pMainPandora);
    for (const auto *pPandora : pandoraList)
    {
        const std::string pandoraName = pPandora->GetName();
        const NDParameters &parameters = this->GetPandoraParameters(pandoraName);

        // The TPC geometry
        const NDSimpleTPCGeom &tpcGeom = this->GetNDSimpleTPCGeom(pandoraName);

        // Create the NDEventInput pointer
        const NDEventInput *pEventInput{nullptr};
        const int volType = parameters.m_volType;

        if (volType == NDParameters::VolType::LArND)
            pEventInput = new LArNDEventInput(pPandora, parameters, tpcGeom);
        else if (volType == NDParameters::VolType::TMS)
            pEventInput = new TMSEventInput(pPandora, parameters, tpcGeom);
        else
            throw pandora::StatusCodeException(pandora::STATUS_CODE_INVALID_PARAMETER);

        // Insert into the event input map
        m_NDEventInputMap.insert(std::make_pair(pandoraName, pEventInput));
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void MainNDPandora::ProcessEvents()
{
    if (this->ShouldDisplayInfo())
        std::cout << "MainNDPandora::ProcessEvents()" << std::endl;

    if (m_NDEventInputMap.empty())
    {
        std::cout << "MainNDPandora::ProcessEvents() - no NDEventInputs found, nothing to process" << std::endl;
        return;
    }

    // Find the total number of events
    const NDEventInput *pFirstEventInput = m_NDEventInputMap.begin()->second;
    const int nEntries = (pFirstEventInput != nullptr) ? pFirstEventInput->GetNEntries() : 0;

    // Check all ND event inputs have the same number of events
    NDEventInputMap::const_iterator iter;
    for (iter = m_NDEventInputMap.begin(); iter != m_NDEventInputMap.end(); ++iter)
    {
        const NDEventInput *pEventInput = iter->second;
        const int checkNEntries = (pEventInput != nullptr) ? pEventInput->GetNEntries() : 0;
        if (checkNEntries != nEntries)
            std::cout << "Warning. NDEventInput " << iter->first << " has " << checkNEntries << " entries != " << nEntries << std::endl;
    }

    // Starting event
    const int startEvt = m_mainParameters.m_nEventsToSkip > 0 ? m_mainParameters.m_nEventsToSkip : 0;
    // Number of events to process, up to nEntries
    const int nProcess = m_mainParameters.m_nEventsToProcess > 0 ? m_mainParameters.m_nEventsToProcess : nEntries;
    // End event, up to nEntries
    const int endEvt = (startEvt + nProcess) < nEntries ? startEvt + nProcess : nEntries;

    std::cout << "Start event is " << startEvt << " and end event is " << endEvt - 1 << std::endl;
    int nTotalHits{0}, nTotalTracks{0};

    // Main event loop, which assumes the same event numbers are used for all Pandora inputs
    for (int iEvt = startEvt; iEvt < endEvt; iEvt++)
    {
        std::cout << std::endl << "   PROCESSING EVENT: " << iEvt << std::endl << std::endl;

        // Make MCParticles, CaloHits & Tracks for each Pandora event input.
        // Keep track of how many CaloHits & Tracks we are making in order to set unique
        // parent addresses for the main Pandora instance (copied to the specific instance).
        // MCParticles should use true IDs to set their parent addresses
        for (iter = m_NDEventInputMap.begin(); iter != m_NDEventInputMap.end(); ++iter)
        {
            const NDEventInput *pEventInput = iter->second;
            if (!pEventInput)
                continue;

            if (this->ShouldDisplayInfo())
                std::cout << "NDEventInput name = " << iter->first << std::endl;
            pEventInput->CreateMCParticles(iEvt);
            nTotalHits = pEventInput->CreateCaloHits(iEvt, nTotalHits);
            nTotalTracks = pEventInput->CreateTracks(iEvt, nTotalTracks);

            if (this->ShouldDisplayInfo())
                std::cout << "nTotalHits = " << nTotalHits << ", nTotalTracks = " << nTotalTracks << std::endl;

            // Run the algorithms for this Pandora event instance
            const pandora::Pandora *pPandora = pEventInput->GetPandoraInstance();
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::ProcessEvent(*pPandora));
        }

        // Further process the event using the main instance, which should include copying
        // the PFOs from the daughter instances to the main one using CreateMainNDListsAlgorithm
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::ProcessEvent(*m_pMainPandora));

        // Reset for the next event
        for (iter = m_NDEventInputMap.begin(); iter != m_NDEventInputMap.end(); ++iter)
        {
            const NDEventInput *pEventInput = iter->second;
            if (!pEventInput)
                continue;

            const pandora::Pandora *pPandora = pEventInput->GetPandoraInstance();
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::Reset(*pPandora));
        }
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::Reset(*m_pMainPandora));
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool MainNDPandora::ShouldDisplayInfo() const
{
    return m_pMainPandora && m_pMainPandora->GetSettings() && m_pMainPandora->GetSettings()->ShouldDisplayAlgorithmInfo();
}

//------------------------------------------------------------------------------------------------------------------------------------------

const NDParameters &MainNDPandora::GetPandoraParameters(const std::string &name) const
{
    NDParametersMap::const_iterator iter = m_NDParametersMap.find(name);

    if (m_NDParametersMap.end() == iter)
        throw pandora::StatusCodeException(pandora::STATUS_CODE_NOT_FOUND);

    return iter->second;
}

//------------------------------------------------------------------------------------------------------------------------------------------

const NDSimpleTPCGeom &MainNDPandora::GetNDSimpleTPCGeom(const std::string &name) const
{
    NDSimpleTPCGeomMap::const_iterator iter = m_NDSimpleTPCGeomMap.find(name);

    if (m_NDSimpleTPCGeomMap.end() == iter)
        throw pandora::StatusCodeException(pandora::STATUS_CODE_NOT_FOUND);

    return iter->second;
}

//------------------------------------------------------------------------------------------------------------------------------------------

const NDEventInput *MainNDPandora::GetNDEventInput(const std::string &name) const
{
    NDEventInputMap::const_iterator iter = m_NDEventInputMap.find(name);

    if (m_NDEventInputMap.end() == iter)
        throw pandora::StatusCodeException(pandora::STATUS_CODE_NOT_FOUND);

    return iter->second;
}

} // namespace lar_nd_reco
