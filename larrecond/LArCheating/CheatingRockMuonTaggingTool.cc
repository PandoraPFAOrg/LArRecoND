/**
 *  @file   larrecond/LArCheating/CheatingRockMuonTaggingTool.cc
 *
 *  @brief  Implementation of the cheating rock muon tagging tool class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArCheating/CheatingSliceIdBaseTool.h"
#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

#include "larrecond/LArCheating/CheatingRockMuonRemovalAlgorithm.h"
#include "larrecond/LArCheating/CheatingRockMuonTaggingTool.h"

using namespace pandora;

namespace lar_content
{

CheatingRockMuonTaggingTool::CheatingRockMuonTaggingTool() :
    m_maxRockMuonFraction(0.25f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CheatingRockMuonTaggingTool::FindAmbiguousPfos(const PfoList &parentRockMuonPfos, PfoList &ambiguousPfos, const MasterAlgorithm *const /*pAlgorithm*/)
{
    if (this->GetPandora().GetSettings()->ShouldDisplayAlgorithmInfo())
        std::cout << "----> Running Algorithm Tool: " << this->GetInstanceName() << ", " << this->GetType() << std::endl;

    PfoList ambiguousParentPfos;

    for (const Pfo *const pParentRockMuonPfo : parentRockMuonPfos)
    {
        PfoList downstreamPfos;
        LArPfoHelper::GetAllDownstreamPfos(pParentRockMuonPfo, downstreamPfos);

        float thisRockMuWeight(0.f), thisTotalWeight(0.f);
        const auto isRockMu = [this](const MCParticle *const pMCParticle)
        { return CheatingRockMuonRemovalAlgorithm::IsRockMuon(this->GetPandora(), pMCParticle); };
        CheatingSliceIdBaseTool::GetTargetParticleWeight(&downstreamPfos, thisRockMuWeight, thisTotalWeight, isRockMu);

        if ((thisTotalWeight > 0.f) && ((thisRockMuWeight / thisTotalWeight) < m_maxRockMuonFraction))
            ambiguousParentPfos.push_back(pParentRockMuonPfo);
    }

    LArPfoHelper::GetAllConnectedPfos(ambiguousParentPfos, ambiguousPfos);
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CheatingRockMuonTaggingTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MaxRockMuonFraction", m_maxRockMuonFraction));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
