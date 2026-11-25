/**
 *  @file   larrecond/LArVertex/CandidateVertexCreationThreeDAlgorithm.cc
 *
 *  @brief  Implementation of the candidate vertex creation algorithm class.
 *
 *  $Log: $
 */

#include <utility>

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"

#include "larrecond/LArVertex/CandidateVertexCreationThreeDAlgorithm.h"

using namespace pandora;

namespace lar_content
{

CandidateVertexCreationThreeDAlgorithm::CandidateVertexCreationThreeDAlgorithm() :
    m_replaceCurrentVertexList(true),
    m_slidingFitWindow(20),
    m_minClusterCaloHits(5),
    m_minClusterLengthSquared(3.f * 3.f),
    m_chiSquaredCut(2.f),
    m_enableEndpointCandidates(true),
    m_maxEndpointXDiscrepancy(4.f),
    m_enableCrossingCandidates(false),
    m_nMaxCrossingCandidates(500),
    m_maxCrossingDiscrepancy(0.5f),
    m_extrapolationNSteps(200),
    m_extrapolationStepSize(0.1f),
    m_maxCrossingSeparationSquared(2.f * 2.f),
    m_minNearbyCrossingDistanceSquared(0.5f * 0.5f),
    m_reducedCandidates(false),
    m_selectionCutFactorMax(2.f),
    m_nClustersPassingMaxCutsPar(26.f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CandidateVertexCreationThreeDAlgorithm::Run()
{
    try
    {
        ClusterVector clusterVector3D;
        this->SelectClusters(clusterVector3D);

        const VertexList *pVertexList(nullptr);
        std::string temporaryListName;
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::CreateTemporaryListAndSetCurrent(*this, pVertexList, temporaryListName));

        if (m_enableEndpointCandidates)
        {
            this->CreateEndpointVertices(clusterVector3D);
        }

        if (m_enableCrossingCandidates)
            this->CreateCrossingCandidates(clusterVector3D);

        if (!m_inputVertexListName.empty())
            this->AddInputVertices();

        if (!pVertexList->empty())
        {
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList<Vertex>(*this, m_outputVertexListName));

            if (m_replaceCurrentVertexList)
                PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<Vertex>(*this, m_outputVertexListName));
        }
    }
    catch (StatusCodeException &statusCodeException)
    {
        this->TidyUp();
        throw statusCodeException;
    }

    this->TidyUp();

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CandidateVertexCreationThreeDAlgorithm::SelectClusters(ClusterVector &clusterVector3D)
{
    if (!clusterVector3D.empty())
        throw StatusCodeException(STATUS_CODE_FAILURE);

    const ClusterList *pClusterList(nullptr);
    PANDORA_THROW_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetList(*this, m_inputClusterListName, pClusterList));

    if (!pClusterList || pClusterList->empty())
    {
        if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
            std::cout << "CandidateVertexCreationThreeDAlgorithm: unable to find cluster list " << m_inputClusterListName << std::endl;

        throw StatusCodeException(STATUS_CODE_FAILURE);
    }

    if (TPC_3D != LArClusterHelper::GetClusterHitType(*(pClusterList->begin())))
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

    ClusterVector sortedClusters(pClusterList->begin(), pClusterList->end());
    std::sort(sortedClusters.begin(), sortedClusters.end(), LArClusterHelper::SortByNHits);

    unsigned int nClustersPassingMaxCuts(0);
    if (m_reducedCandidates)
    {
        for (const Cluster *const pCluster : sortedClusters)
        {
            float selectionCutFactor(1.f);

            if (pCluster->GetParticleId() == E_MINUS)
                selectionCutFactor = m_selectionCutFactorMax;

            if (pCluster->GetNCaloHits() < m_minClusterCaloHits * selectionCutFactor)
                continue;

            if (LArClusterHelper::GetLengthSquared(pCluster) < m_minClusterLengthSquared * selectionCutFactor * selectionCutFactor)
                continue;

            nClustersPassingMaxCuts++;
        }
    }

    for (const Cluster *const pCluster : sortedClusters)
    {
        float selectionCutFactor(1.f);

        if (pCluster->GetParticleId() == E_MINUS && m_reducedCandidates)
        {
            selectionCutFactor = (m_selectionCutFactorMax + 1.f) * 0.5f +
                (m_selectionCutFactorMax - 1.f) * 0.5f * std::tanh(static_cast<float>(nClustersPassingMaxCuts) - m_nClustersPassingMaxCutsPar);
        }

        if (pCluster->GetNCaloHits() < m_minClusterCaloHits * selectionCutFactor)
            continue;

        if (LArClusterHelper::GetLengthSquared(pCluster) < m_minClusterLengthSquared * selectionCutFactor * selectionCutFactor)
            continue;

        try
        {
            this->AddToSlidingFitCache(pCluster);
            clusterVector3D.push_back(pCluster);
        }
        catch (StatusCodeException &statusCodeException)
        {
            if (STATUS_CODE_FAILURE == statusCodeException.GetStatusCode())
                throw statusCodeException;
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CandidateVertexCreationThreeDAlgorithm::CreateEndpointVertices(const ClusterVector &clusterVector3D) const
{
    for (const Cluster *const pCluster : clusterVector3D)
    {
        const ThreeDSlidingFitResult &fitResult(this->GetCachedSlidingFitResult(pCluster));

        const CartesianVector minLayerPosition(fitResult.GetGlobalMinLayerPosition());
        const CartesianVector maxLayerPosition(fitResult.GetGlobalMaxLayerPosition());

        PandoraContentApi::Vertex::Parameters parameters;
        parameters.m_position = minLayerPosition;
        parameters.m_vertexLabel = VERTEX_INTERACTION;
        parameters.m_vertexType = VERTEX_3D;

        const Vertex *pMinVertex(nullptr);
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Vertex::Create(*this, parameters, pMinVertex));

        parameters.m_position = maxLayerPosition;
        const Vertex *pMaxVertex(nullptr);
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Vertex::Create(*this, parameters, pMaxVertex));
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CandidateVertexCreationThreeDAlgorithm::CreateCrossingCandidates(const ClusterVector &clusterVector3D) const
{
    CartesianPointVector crossings3D;
    this->FindCrossingPoints(clusterVector3D, crossings3D);

    unsigned int nCrossingCandidates(0);
    this->CreateCrossingVertices(crossings3D, nCrossingCandidates);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CandidateVertexCreationThreeDAlgorithm::FindCrossingPoints(const ClusterVector &clusterVector, CartesianPointVector &crossingPoints) const
{
    ClusterToSpacepointsMap clusterToSpacepointsMap;

    for (const Cluster *const pCluster : clusterVector)
    {
        ClusterToSpacepointsMap::iterator mapIter(clusterToSpacepointsMap.emplace(pCluster, CartesianPointVector()).first);
        this->GetSpacepoints(pCluster, mapIter->second);
    }

    for (const Cluster *const pCluster1 : clusterVector)
    {
        for (const Cluster *const pCluster2 : clusterVector)
        {
            if (pCluster1 == pCluster2)
                continue;

            this->FindCrossingPoints(clusterToSpacepointsMap.at(pCluster1), clusterToSpacepointsMap.at(pCluster2), crossingPoints);
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CandidateVertexCreationThreeDAlgorithm::GetSpacepoints(const Cluster *const pCluster, CartesianPointVector &spacepoints) const
{
    LArClusterHelper::GetCoordinateVector(pCluster, spacepoints);
    /*
    const ThreeDSlidingFitResult &fitResult(this->GetCachedSlidingFitResult(pCluster));
    const float minLayerRL(fitResult.GetL(fitResult.GetMinLayer()));
    const float maxLayerRL(fitResult.GetL(fitResult.GetMaxLayer()));

    for (unsigned int iStep = 0; iStep < m_extrapolationNSteps; ++iStep)
    {
        const float deltaRL(static_cast<float>(iStep) * m_extrapolationStepSize);

        CartesianVector positionPositive(0.f, 0.f, 0.f), positionNegative(0.f, 0.f, 0.f);
        fitResult.GetExtrapolatedPosition(maxLayerRL + deltaRL, positionPositive);
        fitResult.GetExtrapolatedPosition(minLayerRL - deltaRL, positionNegative);

        spacepoints.push_back(positionPositive);
        spacepoints.push_back(positionNegative);
    }
*/
    std::sort(spacepoints.begin(), spacepoints.end(), LArClusterHelper::SortCoordinatesByPosition);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CandidateVertexCreationThreeDAlgorithm::FindCrossingPoints(
    const CartesianPointVector &spacepoints1, const CartesianPointVector &spacepoints2, CartesianPointVector &crossingPoints) const
{
    bool bestCrossingFound(false);
    float bestSeparationSquared(m_maxCrossingSeparationSquared);
    CartesianVector bestPosition1(0.f, 0.f, 0.f), bestPosition2(0.f, 0.f, 0.f);

    for (const CartesianVector &position1 : spacepoints1)
    {
        for (const CartesianVector &position2 : spacepoints2)
        {
            const float separationSquared((position1 - position2).GetMagnitudeSquared());

            if (separationSquared < bestSeparationSquared)
            {
                bestCrossingFound = true;
                bestSeparationSquared = separationSquared;
                bestPosition1 = position1;
                bestPosition2 = position2;
            }
        }
    }

    if (bestCrossingFound)
    {
        bool alreadyPopulated(false);

        for (const CartesianVector &existingPosition : crossingPoints)
        {
            if (((existingPosition - bestPosition1).GetMagnitudeSquared() < m_minNearbyCrossingDistanceSquared) ||
                ((existingPosition - bestPosition2).GetMagnitudeSquared() < m_minNearbyCrossingDistanceSquared))
            {
                alreadyPopulated = true;
                break;
            }
        }

        if (!alreadyPopulated)
        {
            crossingPoints.push_back(bestPosition1);
            crossingPoints.push_back(bestPosition2);
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CandidateVertexCreationThreeDAlgorithm::CreateCrossingVertices(const CartesianPointVector &crossingPoints3D, unsigned int &nCrossingCandidates) const
{

    for (const CartesianVector &position3D : crossingPoints3D)
    {
        if (nCrossingCandidates > m_nMaxCrossingCandidates)
            return;

        PandoraContentApi::Vertex::Parameters parameters;
        parameters.m_position = position3D;
        parameters.m_vertexLabel = VERTEX_INTERACTION;
        parameters.m_vertexType = VERTEX_3D;

        const Vertex *pVertex(nullptr);
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Vertex::Create(*this, parameters, pVertex));
        ++nCrossingCandidates;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CandidateVertexCreationThreeDAlgorithm::AddInputVertices() const
{
    const VertexList *pInputVertexList{nullptr};
    try
    { // ATTN - No guarantee the list has been initialised, but silent failure here is ok
        PandoraContentApi::GetList(*this, m_inputVertexListName, pInputVertexList);
        if (!pInputVertexList)
            return;

        for (const Vertex *pInputVertex : *pInputVertexList)
        {
            PandoraContentApi::Vertex::Parameters parameters;
            parameters.m_position = pInputVertex->GetPosition();
            parameters.m_vertexLabel = VERTEX_INTERACTION;
            parameters.m_vertexType = VERTEX_3D;

            const Vertex *pVertex(nullptr);
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Vertex::Create(*this, parameters, pVertex));
        }
    }
    catch (const StatusCodeException &)
    {
        return;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CandidateVertexCreationThreeDAlgorithm::AddToSlidingFitCache(const Cluster *const pCluster)
{
    const float slidingFitPitch(LArGeometryHelper::GetWireZPitch(this->GetPandora()));
    const ThreeDSlidingFitResult slidingFitResult(pCluster, m_slidingFitWindow, slidingFitPitch);

    if (!m_slidingFitResultMap.insert(ThreeDSlidingFitResultMap::value_type(pCluster, slidingFitResult)).second)
        throw StatusCodeException(STATUS_CODE_FAILURE);
}

//------------------------------------------------------------------------------------------------------------------------------------------

const ThreeDSlidingFitResult &CandidateVertexCreationThreeDAlgorithm::GetCachedSlidingFitResult(const Cluster *const pCluster) const
{
    ThreeDSlidingFitResultMap::const_iterator iter = m_slidingFitResultMap.find(pCluster);

    if (m_slidingFitResultMap.end() == iter)
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

    return iter->second;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CandidateVertexCreationThreeDAlgorithm::TidyUp()
{
    m_slidingFitResultMap.clear();
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CandidateVertexCreationThreeDAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "InputClusterListName", m_inputClusterListName));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "InputVertexListName", m_inputVertexListName));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "OutputVertexListName", m_outputVertexListName));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "ReplaceCurrentVertexList", m_replaceCurrentVertexList));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "SlidingFitWindow", m_slidingFitWindow));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinClusterCaloHits", m_minClusterCaloHits));

    float minClusterLength = std::sqrt(m_minClusterLengthSquared);
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinClusterLength", minClusterLength));
    m_minClusterLengthSquared = minClusterLength * minClusterLength;

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "ChiSquaredCut", m_chiSquaredCut));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "EnableEndpointCandidates", m_enableEndpointCandidates));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "MaxEndpointXDiscrepancy", m_maxEndpointXDiscrepancy));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "EnableCrossingCandidates", m_enableCrossingCandidates));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "NMaxCrossingCandidates", m_nMaxCrossingCandidates));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MaxCrossingDiscrepancy", m_maxCrossingDiscrepancy));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "ExtrapolationNSteps", m_extrapolationNSteps));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "ExtrapolationStepSize", m_extrapolationStepSize));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "ReducedCandidates", m_reducedCandidates));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "SelectionCutFactorMax", m_selectionCutFactorMax));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "NClustersPassingMaxCutsPar", m_nClustersPassingMaxCutsPar));

    float maxCrossingSeparation = std::sqrt(m_maxCrossingSeparationSquared);
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MaxCrossingSeparation", maxCrossingSeparation));
    m_maxCrossingSeparationSquared = maxCrossingSeparation * maxCrossingSeparation;

    float minNearbyCrossingDistance = std::sqrt(m_minNearbyCrossingDistanceSquared);
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "MinNearbyCrossingDistance", minNearbyCrossingDistance));
    m_minNearbyCrossingDistanceSquared = minNearbyCrossingDistance * minNearbyCrossingDistance;

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
