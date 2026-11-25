/**
 *  @file   larrecond/LArThreeDReco/LArTrackReco/MergeClearTracksThreeDAlgorithm.cc
 *
 *  @brief  Implementation of the 3D track merging algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArPointingClusterHelper.h"

#include "larpandoracontent/LArObjects/LArThreeDSlidingFitResult.h"

#include "larrecond/LArThreeDReco/LArTrackReco/MergeClearTracksThreeDAlgorithm.h"

using namespace pandora;

namespace lar_content
{

MergeClearTracksThreeDAlgorithm::MergeClearTracksThreeDAlgorithm() :
    m_slidingFitWindow(10),
    m_maxGapLengthCut(25.f),
    m_maxGapTransverseCut(1.f),
    m_minCosThetaCut(0.96f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode MergeClearTracksThreeDAlgorithm::Run()
{
    bool madeMerges{true};
    while (madeMerges)
    {
        const ClusterList *pClusterList3D{nullptr};
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_inputClusterListName, pClusterList3D));

        madeMerges = this->FindMerges(pClusterList3D);
    }
    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool MergeClearTracksThreeDAlgorithm::FindMerges(const ClusterList *const pClusterList) const
{
    // Use a cluster vector so that we can sort track-like clusters by the number of hits
    ClusterVector sortedClusters;
    for (const Cluster *const pCluster : *pClusterList)
    {
        if (MU_MINUS != pCluster->GetParticleId())
            continue;

        sortedClusters.emplace_back(pCluster);
    }
    if (sortedClusters.size() < 2)
        return false;

    std::sort(sortedClusters.begin(), sortedClusters.end(), LArClusterHelper::SortByNHits);

    ClusterMergeMap mergeCandidates;
    ClusterFloatMap mergeClosestDistance;

    for (ClusterVector::iterator iter1 = sortedClusters.begin(); iter1 != sortedClusters.end(); ++iter1)
    {
        for (ClusterVector::iterator iter2 = iter1 + 1; iter2 != sortedClusters.end(); ++iter2)
        {
            this->CanMergeClusters(*iter1, *iter2, mergeCandidates, mergeClosestDistance);
        }
    }

    bool madeMerges(!mergeCandidates.empty());
    if (madeMerges)
    {
        ClusterSet usedClusters;
        for (auto const &pair : mergeCandidates)
        {
            if (!usedClusters.count(pair.second) && !usedClusters.count(pair.first))
            {
                this->MergeClusters(pair.first, pair.second);
                usedClusters.insert(pair.first);
                usedClusters.insert(pair.second);
            }
        }
        return true;
    }
    else
        return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void MergeClearTracksThreeDAlgorithm::CanMergeClusters(const Cluster *const pLargeCluster, const Cluster *const pSmallCluster,
    ClusterMergeMap &mergeCandidates, ClusterFloatMap &mergeDistance) const
{
    try
    {
        const LArPointingCluster largePointingCluster(pLargeCluster, m_slidingFitWindow, LArGeometryHelper::GetWireZPitch(this->GetPandora()));
        const LArPointingCluster smallPointingCluster(pSmallCluster, m_slidingFitWindow, LArGeometryHelper::GetWireZPitch(this->GetPandora()));

        LArPointingCluster::Vertex largeClusterVertex, smallClusterVertex;
        LArPointingClusterHelper::GetClosestVertices(largePointingCluster, smallPointingCluster, largeClusterVertex, smallClusterVertex);

        float longitudinalGap(std::numeric_limits<float>::max());
        float transverseGap(std::numeric_limits<float>::max());
        LArPointingClusterHelper::GetImpactParameters(largeClusterVertex, smallClusterVertex, longitudinalGap, transverseGap);
        longitudinalGap = std::fabs(longitudinalGap);
        transverseGap = std::fabs(transverseGap);

        const CartesianVector &largeClusterDirection(largeClusterVertex.GetDirection());
        const CartesianVector &smallClusterDirection(smallClusterVertex.GetDirection());
        const float cosAngle = largeClusterDirection.GetCosOpeningAngle(smallClusterDirection * -1.0);

        if (longitudinalGap <= m_maxGapLengthCut && transverseGap <= m_maxGapTransverseCut && cosAngle >= m_minCosThetaCut)
        {
            if (!mergeCandidates.count(pLargeCluster))
            {
                mergeCandidates.insert(std::make_pair(pLargeCluster, pSmallCluster));
                mergeDistance.insert(std::make_pair(pLargeCluster, longitudinalGap));
            }
            else
            {
                if (longitudinalGap < mergeDistance[pLargeCluster])
                {
                    mergeCandidates[pLargeCluster] = pSmallCluster;
                    mergeDistance[pLargeCluster] = longitudinalGap;
                }
            }
        }
    }
    catch (const StatusCodeException &)
    {
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void MergeClearTracksThreeDAlgorithm::MergeClusters(const Cluster *const pLargeCluster, const Cluster *const pSmallCluster) const
{
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=,
        PandoraContentApi::MergeAndDeleteClusters(*this, pLargeCluster, pSmallCluster, m_inputClusterListName, m_inputClusterListName));
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode MergeClearTracksThreeDAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "SlidingFitWindow", m_slidingFitWindow));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MaxGapLengthCut", m_maxGapLengthCut));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MaxGapTransverseCut", m_maxGapTransverseCut));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinCosThetaCut", m_minCosThetaCut));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "InputClusterListName", m_inputClusterListName));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
