/**
 *  @file   src/PfoThreeDHitAssignmentAlgorithm.cc
 *
 *  @brief
 *
 *  $Log: $
 */

#include "Api/PandoraApi.h"
#include "Objects/CartesianVector.h"
#include "Objects/ParticleFlowObject.h"
#include "Pandora/AlgorithmHeaders.h"

#include "Pandora/PandoraInternal.h"
#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

#include "PfoThreeDHitAssignmentAlgorithm.h"

#include <limits>

using namespace pandora;

namespace lar_content
{

PfoThreeDHitAssignmentAlgorithm::PfoThreeDHitAssignmentAlgorithm() :
    m_inputCaloHitList3DName{""}
{
}

pandora::StatusCode PfoThreeDHitAssignmentAlgorithm::Run()
{
    if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
        std::cout << "----> Running Algorithm: " << this->GetInstanceName() << ", " << this->GetType() << std::endl;

    const CaloHitList *pCaloHits3D{nullptr};
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_inputCaloHitList3DName, pCaloHits3D));

    CaloHitList availableHits;
    std::map<const CaloHit *, float> availableHitUPos, availableHitVPos, availableHitWPos;
    for (const CaloHit *pCaloHit : (*pCaloHits3D))
    {
        if (!PandoraContentApi::IsAvailable(*this, pCaloHit))
            continue;

        availableHits.emplace_back(pCaloHit);
        const CartesianVector pos3D = pCaloHit->GetPositionVector();
        availableHitUPos[pCaloHit] = PandoraContentApi::GetPlugins(*this)->GetLArTransformationPlugin()->YZtoU(pos3D.GetY(), pos3D.GetZ());
        availableHitVPos[pCaloHit] = PandoraContentApi::GetPlugins(*this)->GetLArTransformationPlugin()->YZtoV(pos3D.GetY(), pos3D.GetZ());
        availableHitWPos[pCaloHit] = PandoraContentApi::GetPlugins(*this)->GetLArTransformationPlugin()->YZtoW(pos3D.GetY(), pos3D.GetZ());
    }

    // Maps to keep track of which 3D hit matches to a 2D hit in a given pfo
    std::map<const ParticleFlowObject *, std::string> pfoToClusterListName;
    std::map<const CaloHit *, const ParticleFlowObject *> availableHitToPfoU, availableHitToPfoV, availableHitToPfoW;

    using HitPositionMap = std::map<std::pair<float, float>, const CaloHit *>;
    std::map<const ParticleFlowObject *, HitPositionMap> pfoToUHitMap, pfoToVHitMap, pfoToWHitMap;

    for (unsigned int i = 0; i < m_inputPfoListNames.size(); ++i)
    {
        const std::string pfoListName(m_inputPfoListNames.at(i));
        const PfoList *pPfoList{nullptr};
        PANDORA_THROW_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetList(*this, pfoListName, pPfoList));
        if (!pPfoList || pPfoList->empty())
        {
            if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
                std::cout << "PfoThreeDHitAssignment: couldn't find pfo list " << pfoListName << std::endl;
            continue;
        }

        for (const ParticleFlowObject *pPfo : (*pPfoList))
        {
            pfoToClusterListName[pPfo] = m_outputClusterListNames.at(i);

            CaloHitList theseCaloHitsU, theseCaloHitsV, theseCaloHitsW;
            LArPfoHelper::GetCaloHits(pPfo, TPC_VIEW_U, theseCaloHitsU);
            LArPfoHelper::GetCaloHits(pPfo, TPC_VIEW_V, theseCaloHitsV);
            LArPfoHelper::GetCaloHits(pPfo, TPC_VIEW_W, theseCaloHitsW);

            for (const CaloHit *const pHit : theseCaloHitsU)
            {
                const CartesianVector pos = pHit->GetPositionVector();
                pfoToUHitMap[pPfo][std::make_pair(pos.GetX(), pos.GetZ())] = pHit;
            }

            for (const CaloHit *const pHit : theseCaloHitsV)
            {
                const CartesianVector pos = pHit->GetPositionVector();
                pfoToVHitMap[pPfo][std::make_pair(pos.GetX(), pos.GetZ())] = pHit;
            }

            for (const CaloHit *const pHit : theseCaloHitsW)
            {
                const CartesianVector pos = pHit->GetPositionVector();
                pfoToWHitMap[pPfo][std::make_pair(pos.GetX(), pos.GetZ())] = pHit;
            }
        }
    }

    for (const CaloHit *const pHit3D : availableHits)
    {
        const auto pos3D = pHit3D->GetPositionVector();
        const auto xPos = pos3D.GetX();
        const auto uPos = availableHitUPos.at(pHit3D);
        const auto vPos = availableHitVPos.at(pHit3D);
        const auto wPos = availableHitWPos.at(pHit3D);

        for (const auto &pfoMapPair : pfoToUHitMap)
        {
            const ParticleFlowObject *const pPfo = pfoMapPair.first;
            const HitPositionMap &uHitMap = pfoMapPair.second;

            if (uHitMap.count(std::make_pair(xPos, uPos)) != 0)
            {
                availableHitToPfoU[pHit3D] = pPfo;
                break;
            }
        }

        for (const auto &pfoMapPair : pfoToVHitMap)
        {
            const ParticleFlowObject *const pPfo = pfoMapPair.first;
            const HitPositionMap &vHitMap = pfoMapPair.second;

            if (vHitMap.count(std::make_pair(xPos, vPos)) != 0)
            {
                availableHitToPfoV[pHit3D] = pPfo;
                break;
            }
        }

        for (const auto &pfoMapPair : pfoToWHitMap)
        {
            const ParticleFlowObject *const pPfo = pfoMapPair.first;
            const HitPositionMap &wHitMap = pfoMapPair.second;

            if (wHitMap.count(std::make_pair(xPos, wPos)) != 0)
            {
                availableHitToPfoW[pHit3D] = pPfo;
                break;
            }
        }
    }

    CaloHitList threeDHitsMatchedToOnePfo;
    CaloHitList threeDHitsMatchedToMultiPfos;

    std::map<const CaloHit *, PfoSet> hits3DToPfosSets;
    for (const CaloHit *const pCaloHit3D : availableHits)
    {
        PfoSet matchedPfos;

        auto uIter = availableHitToPfoU.find(pCaloHit3D);
        if (uIter != availableHitToPfoU.end())
            matchedPfos.insert(uIter->second);

        auto vIter = availableHitToPfoV.find(pCaloHit3D);
        if (vIter != availableHitToPfoV.end())
            matchedPfos.insert(vIter->second);

        auto wIter = availableHitToPfoW.find(pCaloHit3D);
        if (wIter != availableHitToPfoW.end())
            matchedPfos.insert(wIter->second);

        if (!hits3DToPfosSets.count(pCaloHit3D))
            hits3DToPfosSets[pCaloHit3D] = matchedPfos;

        const unsigned int nPfos(matchedPfos.size());
        if (0 == nPfos)
            continue;
        else if (1 == nPfos)
            threeDHitsMatchedToOnePfo.emplace_back(pCaloHit3D);
        else
            threeDHitsMatchedToMultiPfos.emplace_back(pCaloHit3D);
    }

    // Swap sets to vectors for easier access
    std::map<const CaloHit *, PfoVector> hits3DToPfos;
    for (const auto &hitToPfos : hits3DToPfosSets)
        hits3DToPfos[hitToPfos.first] = PfoVector(hitToPfos.second.begin(), hitToPfos.second.end());

    std::map<const ParticleFlowObject *, CaloHitList> pfoToHits;

    // Deal with the unambiguous 3D hits first
    for (const CaloHit *const pCaloHit3D : threeDHitsMatchedToOnePfo)
    {
        const ParticleFlowObject *const pPfo = hits3DToPfos[pCaloHit3D][0];
        if (!pfoToHits.count(pPfo))
            pfoToHits[pPfo] = CaloHitList();
        pfoToHits[pPfo].emplace_back(pCaloHit3D);
    }

    // If a 3D hit has 2D hits in more than one pfo we have a decision to make
    for (const CaloHit *const pCaloHit3D : threeDHitsMatchedToMultiPfos)
    {
        if (!hits3DToPfos.count(pCaloHit3D))
            continue;

        const PfoVector &pfoList = hits3DToPfos[pCaloHit3D];

        int bestIndex{-999};
        int mostHits{0};
        int biggestPfoIndex{-999};
        float minDist = std::numeric_limits<float>::max();

        for (unsigned int iPfo = 0; iPfo < pfoList.size(); ++iPfo)
        {
            const ParticleFlowObject *const pPfo = pfoList[iPfo];

            const int n2DHits = LArPfoHelper::GetNumberOfTwoDHits(pPfo);
            if (mostHits < n2DHits)
            {
                mostHits = n2DHits;
                biggestPfoIndex = iPfo;
            }

            if (!pfoToHits.count(pPfo))
                continue;

            const CartesianVector hit3DPos = pCaloHit3D->GetPositionVector();
            for (const CaloHit *const pPfoHit : pfoToHits.at(pPfo))
            {
                const float dist = (hit3DPos - pPfoHit->GetPositionVector()).GetMagnitude();
                if (minDist > dist)
                {
                    minDist = dist;
                    bestIndex = iPfo;
                }
            }
        }

        // If no best index was found, use the pfo with the most 2D hits
        if (bestIndex < 0)
            bestIndex = biggestPfoIndex;

        if (!pfoToHits.count(pfoList[bestIndex]))
            pfoToHits[pfoList[bestIndex]] = CaloHitList();
        pfoToHits[pfoList[bestIndex]].emplace_back(pCaloHit3D);
    }

    // Assign the hits
    for (auto const &pfoPair : pfoToHits)
        AddHitsToPfo(pfoPair.first, pfoPair.second, pfoToClusterListName.at(pfoPair.first));

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PfoThreeDHitAssignmentAlgorithm::AddHitsToPfo(const ParticleFlowObject *pPfo, const CaloHitList &hits, const std::string listName) const
{
    ClusterList clusters3D;
    LArPfoHelper::GetThreeDClusterList(pPfo, clusters3D);
    const unsigned int nClusters(clusters3D.size());

    if (0 == nClusters)
    {
        const ClusterList *pClusterList(nullptr);
        std::string clusterListName;
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::CreateTemporaryListAndSetCurrent(*this, pClusterList, clusterListName));

        PandoraContentApi::Cluster::Parameters clusterParams;
        clusterParams.m_caloHitList.insert(clusterParams.m_caloHitList.end(), hits.begin(), hits.end());

        const Cluster *pCluster3D(nullptr);
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Cluster::Create(*this, clusterParams, pCluster3D));
        if (!pCluster3D || !pClusterList || pClusterList->empty())
            throw StatusCodeException(STATUS_CODE_FAILURE);

        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList<Cluster>(*this, listName));
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::AddToPfo(*this, pPfo, pCluster3D));
    }
    else if (1 == nClusters)
    {
        const Cluster *const pCluster3D = *(clusters3D.begin());
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::AddToCluster(*this, pCluster3D, &hits));
    }
    else
        throw StatusCodeException(STATUS_CODE_FAILURE);
}

//------------------------------------------------------------------------------------------------------------------------------------------

pandora::StatusCode PfoThreeDHitAssignmentAlgorithm::ReadSettings(const pandora::TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "InputCaloHitList3DName", m_inputCaloHitList3DName));
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadVectorOfValues(xmlHandle, "InputPfoListNames", m_inputPfoListNames));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadVectorOfValues(xmlHandle, "OutputClusterListNames", m_outputClusterListNames));

    if (m_inputPfoListNames.size() != m_outputClusterListNames.size())
    {
        std::cout << "LArPfoThreeDHitAssignment: number of input pfo list names must be the same as the number of output cluster list names"
                  << std::endl;
        return STATUS_CODE_FAILURE;
    }

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
