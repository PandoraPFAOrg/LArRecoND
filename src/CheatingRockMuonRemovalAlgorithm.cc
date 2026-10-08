/**
 *  @file   src/CheatingRockMuonRemovalAlgorithm.cc
 *
 *  @brief  Implementation of the cheating rock muon removal algorithm class.
 *
 *  $Log: $
 */

#include "Objects/CaloHit.h"
#include "Objects/MCParticle.h"
#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"
#include "larpandoracontent/LArHelpers/LArVertexHelper.h"

#include "CheatingRockMuonRemovalAlgorithm.h"

using namespace pandora;

namespace lar_content
{

bool CheatingRockMuonRemovalAlgorithm::IsRockMuon(const Pandora &pandora, const MCParticle *const pMCParticle)
{
    // Get the parent MCParticle.
    // This means things like delta rays, showers etc are all traced back to the original muon, and classified as a rock muon.
    const auto parentMCParticle = LArMCParticleHelper::GetParentMCParticle(pMCParticle);

    const bool isNuMu = (std::abs(parentMCParticle->GetParticleId()) == 14);
    const bool vertexInsideFV = LArVertexHelper::IsInFiducialVolume(pandora, parentMCParticle->GetVertex(), "dune_nd");

    return isNuMu && !vertexInsideFV;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CheatingRockMuonRemovalAlgorithm::Run()
{
    // Load the input lists
    const CaloHitList *pCaloHitList(nullptr);

    if (PandoraContentApi::GetList(*this, m_inputCaloHitListName, pCaloHitList) != STATUS_CODE_SUCCESS)
    {
        std::cout << "CheatingRockMuonRemovalAlgorithm: No calo hit list found with name: " << m_inputCaloHitListName << ", skipping!" << std::endl;
        return STATUS_CODE_FAILURE;
    }

    // Create the output lists
    CaloHitList pRockMuonCaloHitList;
    CaloHitList pNeutrinoCaloHitList;

    // Process the hits
    for (const CaloHit *const pCaloHit : *pCaloHitList)
    {
        try
        {
            const MCParticle *const pMCParticle = MCParticleHelper::GetMainMCParticle(pCaloHit);

            if (this->IsRockMuon(this->GetPandora(), pMCParticle))
                pRockMuonCaloHitList.push_back(pCaloHit);
            else
                pNeutrinoCaloHitList.push_back(pCaloHit);
        }
        catch (StatusCodeException &)
        {
            // INFO: There may be small number of hits with no associated MC... we don't want to break on them.
            continue;
        }
    }

    if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
        std::cout << "CheatingRockMuonRemovalAlgorithm: Rock muon hits: " << pRockMuonCaloHitList.size() << ", Neutrino hits: " << pNeutrinoCaloHitList.size() << std::endl;

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList(*this, pRockMuonCaloHitList, m_rockMuonCaloHitListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList(*this, pNeutrinoCaloHitList, m_neutrinoCaloHitListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<CaloHit>(*this, m_rockMuonCaloHitListName));

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CheatingRockMuonRemovalAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "InputCaloHitListName", m_inputCaloHitListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "RockMuonCaloHitListName", m_rockMuonCaloHitListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "NeutrinoCaloHitListName", m_neutrinoCaloHitListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "MCParticleListName", m_inputMCParticleListName));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
