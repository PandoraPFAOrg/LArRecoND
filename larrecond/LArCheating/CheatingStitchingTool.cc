/**
 *  @file   larrecond/LArCheating/CheatingStitchingTool.cc
 *
 *  @brief  Implementation of the cheating stitching tool class.
 *
 *  $Log: $
 */

#include <larpandoracontent/LArHelpers/LArMCParticleHelper.h>
#include <larpandoracontent/LArHelpers/LArPfoHelper.h>
#include <larpandoracontent/LArThreeDReco/LArPfoStitching/StitchingPfoOperations.h>

#include "larrecond/LArCheating/CheatingStitchingTool.h"

using namespace pandora;

namespace lar_content
{
CheatingStitchingTool::CheatingStitchingTool()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CheatingStitchingTool::RunStitching(const Algorithm *const pAlgorithm, const StitchingPfoOperations *const pStitchingOperations,
    const PfoList *const /*pMultiPfoList*/, PfoToLArTPCMap &pfoToLArTPCMap, PfoToFloatMap & /*stitchedPfosToX0Map*/)
{
    if (PandoraContentApi::GetSettings(*pAlgorithm)->ShouldDisplayAlgorithmInfo())
        std::cout << "----> Running Algorithm Tool: " << this->GetInstanceName() << ", " << this->GetType() << std::endl;

    if (this->GetPandora().GetGeometry()->GetLArTPCMap().size() < 2)
        return;

    if (pfoToLArTPCMap.empty())
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

    const PfoList *pPfoList{nullptr};
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*pAlgorithm, pPfoList));

    std::map<const MCParticle *, std::map<unsigned int, std::vector<const ParticleFlowObject *>>> mcParticleToLArTPCToPfosMap;

    // Find all PFOs associated to each MCParticle
    for (const ParticleFlowObject *const pPfo : *pPfoList)
    {
        CaloHitList allHits;
        LArPfoHelper::GetAllCaloHits2D(pPfo, allHits);
        std::map<const MCParticle *, unsigned int> mcParticleToHitCountMap;

        for (const CaloHit *const pCaloHit : allHits)
        {
            float maxWeight{std::numeric_limits<float>::lowest()};
            const MCParticle *pMainMC{nullptr};
            for (const auto &[pMC, weight] : pCaloHit->GetMCParticleWeightMap())
            {
                if (weight > maxWeight)
                {
                    pMainMC = pMC;
                    maxWeight = weight;
                }
            }

            // Consider only the MCParticle with the most weight associated to this hit
            if (pMainMC)
                ++mcParticleToHitCountMap[pMainMC];
        }

        // Find the MCParticle with the most hits associated to this PFO
        const MCParticle *pMCParticle{nullptr};
        size_t maxHits(0);

        for (const auto &mcParticleHitCountPair : mcParticleToHitCountMap)
        {
            if (mcParticleHitCountPair.second > maxHits)
            {
                maxHits = mcParticleHitCountPair.second;
                pMCParticle = mcParticleHitCountPair.first;
            }
        }

        // Add this PFO to the map for this MCParticle
        if (nullptr == pMCParticle)
            continue;

        const unsigned int tpcVolumeId = pfoToLArTPCMap.at(pPfo)->GetLArTPCVolumeId();
        mcParticleToLArTPCToPfosMap[pMCParticle][tpcVolumeId].push_back(pPfo);
    }

    // Now, join together any PFOs from different TPCs associated to the same MCParticle
    for (const auto &mcParticlePfosPair : mcParticleToLArTPCToPfosMap)
    {
        const std::map<unsigned int, std::vector<const ParticleFlowObject *>> &tpcToPfosMap = mcParticlePfosPair.second;

        // Nothing to stitch if only one TPC involved
        if (tpcToPfosMap.size() < 2)
            continue;

        // Sort the PFOs in each TPC by number of hits.
        // Since this is JUST cheated stitching, we will take the largest PFO in each TPC as the "main" one to stitch to.
        std::map<unsigned int, const ParticleFlowObject *> tpcToMainPfoMap;
        for (const auto &tpcPfosPair : tpcToPfosMap)
        {
            const unsigned int tpcVolumeId = tpcPfosPair.first;
            const std::vector<const ParticleFlowObject *> &pfosInTpc = tpcPfosPair.second;
            const ParticleFlowObject *pLargestPfo{nullptr};

            size_t maxHits(0);
            for (const ParticleFlowObject *const pPfo : pfosInTpc)
            {
                const size_t nHits = LArPfoHelper::GetNumberOfThreeDHits(pPfo) + LArPfoHelper::GetNumberOfTwoDHits(pPfo);
                if (nHits > maxHits)
                {
                    maxHits = nHits;
                    pLargestPfo = pPfo;
                }
            }

            tpcToMainPfoMap[tpcVolumeId] = pLargestPfo;
        }

        // Stitch all other PFOs in this TPC to the main one
        auto tpcToMainPfoIter = tpcToMainPfoMap.begin();
        const ParticleFlowObject *const pReferencePfo = tpcToMainPfoIter->second;
        ++tpcToMainPfoIter;

        for (; tpcToMainPfoIter != tpcToMainPfoMap.end(); ++tpcToMainPfoIter)
        {
            const ParticleFlowObject *const pCurrentPfo = tpcToMainPfoIter->second;

            // Stitch the current PFO to the reference PFO
            pStitchingOperations->StitchPfos(pReferencePfo, pCurrentPfo, pfoToLArTPCMap);
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CheatingStitchingTool::ReadSettings(const TiXmlHandle /*xmlHandle*/)
{

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
