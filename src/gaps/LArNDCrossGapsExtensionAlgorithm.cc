/**
 *  @file   larpandoracontent/LArTwoDReco/LArClusterAssociation/LArNDCrossGapsExtensionAlgorithm.cc
 *
 *  @brief  Implementation of the cross gaps extension algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArPointingClusterHelper.h"

#include "LArNDCrossGapsExtensionAlgorithm.h"

using namespace pandora;

namespace lar_content
{

LArNDCrossGapsExtensionAlgorithm::LArNDCrossGapsExtensionAlgorithm() :
    m_minClusterLength(5.f),
    m_minGapFraction(0.5f),
    m_maxGapTolerance(2.f),
    m_maxTransverseDisplacement(2.5f),
    m_maxRelativeAngle(10.f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArNDCrossGapsExtensionAlgorithm::GetListOfCleanClusters(const ClusterList *const pClusterList, ClusterVector &clusterVector) const
{
    // ATTN May want to opt-out completely if no gap information available
    // if (PandoraContentApi::GetGeometry(*this)->GetDetectorGapList().empty())
    //     return;

    for (const Cluster *const pCluster : *pClusterList)
    {
        if (LArClusterHelper::GetLengthSquared(pCluster) < m_minClusterLength * m_minClusterLength)
            continue;

        clusterVector.push_back(pCluster);
    }

    std::sort(clusterVector.begin(), clusterVector.end(), LArClusterHelper::SortByNHits);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArNDCrossGapsExtensionAlgorithm::FillClusterAssociationMatrix(const ClusterVector &clusterVector, ClusterAssociationMatrix &clusterAssociationMatrix) const
{
    // Build lists of pointing clusters in proximity to gaps
    LArPointingClusterList innerPointingClusterList, outerPointingClusterList;
    this->BuildPointingClusterList(clusterVector, innerPointingClusterList, outerPointingClusterList);

    // Form associations between pairs of pointing clusters
    for (const LArPointingCluster &pointingClusterInner : innerPointingClusterList)
    {
        const LArPointingCluster::Vertex &pointingVertexInner(pointingClusterInner.GetInnerVertex());

        for (const LArPointingCluster &pointingClusterOuter : outerPointingClusterList)
        {
            const LArPointingCluster::Vertex &pointingVertexOuter(pointingClusterOuter.GetOuterVertex());

            if( pointingClusterOuter.GetCluster() == pointingClusterInner.GetCluster() )
                continue;

            if (!this->IsAssociated(pointingVertexInner, pointingVertexOuter))
                continue;

            if (!this->IsAcrossGap(pointingVertexOuter.GetPosition(), pointingVertexInner.GetPosition(), 
                        LArClusterHelper::GetClusterHitType(pointingClusterInner.GetCluster())))
                continue;

            const Cluster *const pClusterInner(pointingClusterInner.GetCluster());
            const Cluster *const pClusterOuter(pointingClusterOuter.GetCluster());

            const float lengthSquaredInner(LArClusterHelper::GetLengthSquared(pClusterInner));
            const float lengthSquaredOuter(LArClusterHelper::GetLengthSquared(pClusterOuter));

            (void)clusterAssociationMatrix[pClusterInner].insert(ClusterAssociationMap::value_type(pClusterOuter,
                ClusterAssociation(ClusterAssociation::INNER, ClusterAssociation::OUTER, ClusterAssociation::STRONG, lengthSquaredOuter)));
            (void)clusterAssociationMatrix[pClusterOuter].insert(ClusterAssociationMap::value_type(pClusterInner,
                ClusterAssociation(ClusterAssociation::OUTER, ClusterAssociation::INNER, ClusterAssociation::STRONG, lengthSquaredInner)));
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArNDCrossGapsExtensionAlgorithm::BuildPointingClusterList(const ClusterVector &clusterVector,
    LArPointingClusterList &innerPointingClusterList, LArPointingClusterList &outerPointingClusterList) const
{
    // Convert each input cluster into a pointing cluster
    LArPointingClusterList pointingClusterList;

    for (const Cluster *const pCluster : clusterVector)
    {
        try
        {
            pointingClusterList.push_back(LArPointingCluster(pCluster));
        }
        catch (StatusCodeException &)
        {
        }
    }

    // Identify clusters adjacent to detector gaps
    this->BuildPointingClusterList(true, pointingClusterList, innerPointingClusterList);
    this->BuildPointingClusterList(false, pointingClusterList, outerPointingClusterList);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArNDCrossGapsExtensionAlgorithm::BuildPointingClusterList(
    const bool useInner, const LArPointingClusterList &inputPointingClusterList, LArPointingClusterList &outputPointingClusterList) const
{
    for (const LArPointingCluster &pointingCluster : inputPointingClusterList)
    {
        const LArPointingCluster::Vertex &pointingVertex(useInner ? pointingCluster.GetInnerVertex() : pointingCluster.GetOuterVertex());
        const HitType hitType(LArClusterHelper::GetClusterHitType(pointingCluster.GetCluster()));

        if (LArGeometryHelper::IsInGap(this->GetPandora(), pointingVertex.GetPosition(), hitType, m_maxGapTolerance))
            outputPointingClusterList.push_back(pointingCluster);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArNDCrossGapsExtensionAlgorithm::IsAssociated(const LArPointingCluster::Vertex &pointingVertex1, const LArPointingCluster::Vertex &pointingVertex2) const
{
    const float ratio{LArGeometryHelper::GetWirePitchRatio(this->GetPandora(), LArClusterHelper::GetClusterHitType(pointingVertex1.GetCluster()))};
    const float maxTransverseDisplacementAdjusted{ratio * m_maxTransverseDisplacement};
    const float maxLongitudinalDisplacement((pointingVertex2.GetPosition() - pointingVertex1.GetPosition()).GetMagnitude());

    const bool isAssociated1(LArPointingClusterHelper::IsEmission(pointingVertex1.GetPosition(), pointingVertex2, -1.f,
        maxLongitudinalDisplacement + 1.f, maxTransverseDisplacementAdjusted, m_maxRelativeAngle));
    const bool isAssociated2(LArPointingClusterHelper::IsEmission(pointingVertex2.GetPosition(), pointingVertex1, -1.f,
        maxLongitudinalDisplacement + 1.f, maxTransverseDisplacementAdjusted, m_maxRelativeAngle));

    return (isAssociated1 && isAssociated2);
}

//------------------------------------------------------------------------------------------------------------------------------------------
float LArNDCrossGapsExtensionAlgorithm::CalculateGapPathLength(const pandora::Pandora &pandora, const pandora::CartesianVector& startPoint,
            const pandora::CartesianVector& endPoint, const pandora::HitType hitType)
{
    const CartesianVector delta(endPoint - startPoint);
    const float length(delta.GetMagnitude());

    if( length < std::numeric_limits<float>::epsilon())
        return 0.f;

    std::vector<std::pair<float, float>> intervals;

    for( const DetectorGap *const pDetectorGap : pandora.GetGeometry()->GetDetectorGapList() )
    {
        const BoxGap *const pBoxGap = dynamic_cast<const BoxGap *>(pDetectorGap);

        if( (!pBoxGap) || (hitType != TPC_3D) )
            throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

        float tMin{0.f}, tMax{1.f};

        for( const CartesianVector &side : {pBoxGap->GetSide1(), pBoxGap->GetSide2(), pBoxGap->GetSide3()} )
        {
            const CartesianVector& axis{side.GetUnitVector()};
            const float sideLength{side.GetMagnitude()};
            const float p0{(startPoint - pBoxGap->GetVertex()).GetDotProduct(axis)};
            const float d{ delta.GetDotProduct(axis) };

            if (std::fabs(d) < std::numeric_limits<float>::epsilon())
            {
                // Unlikely to cross the current gap if the projection into
                // a face of the gap is almost 0 and start position is before or
                // after the gap vertex
                if ((p0 < 0.f) || (p0 > sideLength))
                {
                    tMax = tMin;
                    break;
                }

                continue;
            }

            // The portion of the distance between start and
            // end point within the gap is: 0 <= p0 + t*d <= sideLength
            // Solve for t.
            float t1{-p0/d}, t2{(sideLength-p0)/d};

            if( t1 > t2)
                std::swap(t1, t2);

            tMin = std::max(tMin, t1);
            tMax = std::min(tMax, t2);

            if( tMax <= tMin )
                break;
        }

        if( tMax > tMin )
            intervals.emplace_back(tMin, tMax);
    }

    std::sort(intervals.begin(), intervals.end());

    float gapFraction{0.f}, currentStart{0.f}, currentEnd{0.f};

    for( const auto [start, end] : intervals )
    {
        if(start > currentEnd)
        {
            gapFraction += (currentEnd - currentStart);
            currentStart = start;
            currentEnd = end;
        }
        else
        {
            currentEnd = std::max(currentEnd, end);
        }
    }

    gapFraction += currentEnd - currentStart;

    return gapFraction * length;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArNDCrossGapsExtensionAlgorithm::IsAcrossGap(const CartesianVector& pOuter, const CartesianVector& pInner, const HitType hitType) const
{
    const float length{ (pInner - pOuter).GetMagnitude() };

    if( length < std::numeric_limits<float>::epsilon() )
        return false;

    const float gapLength{ CalculateGapPathLength(this->GetPandora(), pOuter, pInner, hitType) };

    if( (gapLength / length) < m_minGapFraction )
        return false;

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArNDCrossGapsExtensionAlgorithm::FillClusterMergeMap(const ClusterAssociationMatrix &inputAssociationMatrix, ClusterMergeMap &clusterMergeMap) const
{
    // Decide which associations will become merges
    // To make the merge A <-> B, both A -> B and B -> A must be strong associations
    // with the largest figures of merit of all the A -> X and B -> Y associations

    // First step: remove double-counting from the map of associations
    // i.e. if the map has A <-> B, B <-> C, A <-> C, then remove A <-> C
    ClusterAssociationMatrix clusterAssociationMatrix;

    ClusterVector sortedInputClusters;
    for (const auto &mapEntry : inputAssociationMatrix)
        sortedInputClusters.push_back(mapEntry.first);
    std::sort(sortedInputClusters.begin(), sortedInputClusters.end(), LArClusterHelper::SortByNHits);

    for (const Cluster *const pCluster1 : sortedInputClusters)
    {
        const ClusterAssociationMap &associationMap1(inputAssociationMatrix.at(pCluster1));

        for (const Cluster *const pCluster2 : sortedInputClusters)
        {
            if (pCluster1 == pCluster2)
                continue;

            const ClusterAssociationMap &associationMap2(inputAssociationMatrix.at(pCluster2));

            ClusterAssociationMap::const_iterator iter12 = associationMap1.find(pCluster2);
            if (associationMap1.end() == iter12)
                continue;

            ClusterAssociationMap::const_iterator iter21 = associationMap2.find(pCluster1);
            if (associationMap2.end() == iter21)
                continue;

            const ClusterAssociation &association12(iter12->second);
            const ClusterAssociation &association21(iter21->second);

            bool isAssociated(true);

            ClusterVector sortedAssociationClusters;
            for (const auto &mapEntry : associationMap1)
                sortedAssociationClusters.push_back(mapEntry.first);
            std::sort(sortedAssociationClusters.begin(), sortedAssociationClusters.end(), LArClusterHelper::SortByNHits);

            for (const Cluster *const pCluster3 : sortedAssociationClusters)
            {
                const ClusterAssociation &association13(associationMap1.at(pCluster3));

                ClusterAssociationMap::const_iterator iter23 = associationMap2.find(pCluster3);
                if (associationMap2.end() == iter23)
                    continue;

                const ClusterAssociation &association23(iter23->second);

                if (association12.GetParent() == association13.GetParent() && association23.GetParent() == association21.GetParent() &&
                    association13.GetDaughter() != association23.GetDaughter())
                {
                    isAssociated = false;
                    break;
                }
            }

            if (isAssociated)
            {
                (void)clusterAssociationMatrix[pCluster1].insert(ClusterAssociationMap::value_type(pCluster2, association12));
                (void)clusterAssociationMatrix[pCluster2].insert(ClusterAssociationMap::value_type(pCluster1, association21));
            }
        }
    }

    // Second step: find the best associations A -> X and B -> Y
    ClusterAssociationMatrix intermediateAssociationMatrix;

    ClusterVector sortedClusters;
    for (const auto &mapEntry : clusterAssociationMatrix)
        sortedClusters.push_back(mapEntry.first);
    std::sort(sortedClusters.begin(), sortedClusters.end(), LArClusterHelper::SortByNHits);

    for (const Cluster *const pParentCluster : sortedClusters)
    {
        const ClusterAssociationMap &clusterAssociationMap(clusterAssociationMatrix.at(pParentCluster));

        const Cluster *pBestClusterInner(nullptr);
        ClusterAssociation bestAssociationInner(ClusterAssociation::UNDEFINED, ClusterAssociation::UNDEFINED, ClusterAssociation::NONE, 0.f);

        const Cluster *pBestClusterOuter(nullptr);
        ClusterAssociation bestAssociationOuter(ClusterAssociation::UNDEFINED, ClusterAssociation::UNDEFINED, ClusterAssociation::NONE, 0.f);

        ClusterVector sortedAssociationClusters;
        for (const auto &mapEntry : clusterAssociationMap)
            sortedAssociationClusters.push_back(mapEntry.first);
        std::sort(sortedAssociationClusters.begin(), sortedAssociationClusters.end(), LArClusterHelper::SortByNHits);

        for (const Cluster *const pDaughterCluster : sortedAssociationClusters)
        {
            const ClusterAssociation &clusterAssociation(clusterAssociationMap.at(pDaughterCluster));

            // Inner associations
            if (clusterAssociation.GetParent() == ClusterAssociation::INNER)
            {
                if (clusterAssociation.GetFigureOfMerit() > bestAssociationInner.GetFigureOfMerit())
                {
                    bestAssociationInner = clusterAssociation;
                    pBestClusterInner = pDaughterCluster;
                }
            }

            // Outer associations
            if (clusterAssociation.GetParent() == ClusterAssociation::OUTER)
            {
                if (clusterAssociation.GetFigureOfMerit() > bestAssociationOuter.GetFigureOfMerit())
                {
                    bestAssociationOuter = clusterAssociation;
                    pBestClusterOuter = pDaughterCluster;
                }
            }
        }

        if (pBestClusterInner)
            (void)intermediateAssociationMatrix[pParentCluster].insert(ClusterAssociationMap::value_type(pBestClusterInner, bestAssociationInner));

        if (pBestClusterOuter)
            (void)intermediateAssociationMatrix[pParentCluster].insert(ClusterAssociationMap::value_type(pBestClusterOuter, bestAssociationOuter));
    }

    // Third step: make the merge if A -> X and B -> Y is in fact A -> B and B -> A
    ClusterVector intermediateSortedClusters;
    for (const auto &mapEntry : intermediateAssociationMatrix)
        intermediateSortedClusters.push_back(mapEntry.first);
    std::sort(intermediateSortedClusters.begin(), intermediateSortedClusters.end(), LArClusterHelper::SortByNHits);

    for (const Cluster *const pParentCluster : intermediateSortedClusters)
    {
        const ClusterAssociationMap &parentAssociationMap(intermediateAssociationMatrix.at(pParentCluster));

        ClusterVector sortedAssociationClusters;
        for (const auto &mapEntry : parentAssociationMap)
            sortedAssociationClusters.push_back(mapEntry.first);
        std::sort(sortedAssociationClusters.begin(), sortedAssociationClusters.end(), LArClusterHelper::SortByNHits);

        for (const Cluster *const pDaughterCluster : sortedAssociationClusters)
        {
            const ClusterAssociation &parentToDaughterAssociation(parentAssociationMap.at(pDaughterCluster));

            ClusterAssociationMatrix::const_iterator iter5 = intermediateAssociationMatrix.find(pDaughterCluster);

            if (intermediateAssociationMatrix.end() == iter5)
                continue;

            const ClusterAssociationMap &daughterAssociationMap(iter5->second);

            ClusterAssociationMap::const_iterator iter6 = daughterAssociationMap.find(pParentCluster);

            if (daughterAssociationMap.end() == iter6)
                continue;

            const ClusterAssociation &daughterToParentAssociation(iter6->second);

            if (parentToDaughterAssociation.GetParent() == daughterToParentAssociation.GetDaughter() &&
                parentToDaughterAssociation.GetDaughter() == daughterToParentAssociation.GetParent())
            {
                ClusterList &parentList(clusterMergeMap[pParentCluster]);

                if (parentList.end() == std::find(parentList.begin(), parentList.end(), pDaughterCluster))
                    parentList.push_back(pDaughterCluster);
            }
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode LArNDCrossGapsExtensionAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinClusterLength", m_minClusterLength));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinGapFraction", m_minGapFraction));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MaxGapTolerance", m_maxGapTolerance));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "MaxTransverseDisplacement", m_maxTransverseDisplacement));

    float maxCosRelativeAngle(std::cos(m_maxRelativeAngle * M_PI / 180.0));
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MaxCosRelativeAngle", maxCosRelativeAngle));
    m_maxRelativeAngle = (180.0 / M_PI) * std::acos(maxCosRelativeAngle);

    return ClusterExtensionAlgorithm::ReadSettings(xmlHandle);
}

} // namespace lar_content
