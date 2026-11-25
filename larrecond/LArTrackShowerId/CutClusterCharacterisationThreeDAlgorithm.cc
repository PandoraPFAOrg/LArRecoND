/**
 *  @file   larrecond/LArTrackShowerId/CutClusterCharacterisationThreeDAlgorithm.cc
 *
 *  @brief  Implementation of the cut based cluster characterisation algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"

#include "larpandoracontent/LArObjects/LArThreeDSlidingFitResult.h"
#include "larpandoracontent/LArTrackShowerId/CutClusterCharacterisationAlgorithm.h"

#include "larrecond/LArTrackShowerId/CutClusterCharacterisationThreeDAlgorithm.h"

using namespace pandora;

namespace lar_content
{

CutClusterCharacterisationThreeDAlgorithm::CutClusterCharacterisationThreeDAlgorithm() :
    m_slidingFitWindow(10),
    m_minCaloHitsCut(6),
    m_maxShowerLengthCut(80.f),
    m_pathLengthRatioCut(1.005f),
    m_rTWidthRatioCut(0.05f),
    m_vertexDistanceRatioCut(0.5f),
    m_showerWidthRatioCut(0.35f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool CutClusterCharacterisationThreeDAlgorithm::IsClearTrack(const Cluster *const pCluster) const
{
    if (pCluster->GetNCaloHits() < m_minCaloHitsCut)
        return false;

    float straightLineLength(-1.f), integratedPathLength(-1.f);
    float rTMinFit1(+std::numeric_limits<float>::max()), rTMaxFit1(-std::numeric_limits<float>::max());
    float rTMinFit2(+std::numeric_limits<float>::max()), rTMaxFit2(-std::numeric_limits<float>::max());

    try
    {
        // We can actually do a sliding cone fit and then extract the sliding fit result. We'll use the cones later

        const ThreeDSlidingFitResult slidingFitResult(pCluster, m_slidingFitWindow, LArGeometryHelper::GetWireZPitch(this->GetPandora()));

        const CartesianVector globalMinLayerPosition(slidingFitResult.GetGlobalMinLayerPosition());
        straightLineLength = (slidingFitResult.GetGlobalMaxLayerPosition() - globalMinLayerPosition).GetMagnitude();

        integratedPathLength = 0.f;
        CartesianVector previousFitPosition(globalMinLayerPosition);

        for (int layer = slidingFitResult.GetMinLayer(); layer < slidingFitResult.GetMaxLayer(); ++layer)
        {
            CartesianVector thisFitPosition(0.f, 0.f, 0.f);
            slidingFitResult.GetGlobalFitPosition(layer, thisFitPosition);
            integratedPathLength += (thisFitPosition - previousFitPosition).GetMagnitude();
            previousFitPosition = thisFitPosition;
        }

        // Find the width from the sliding fits along the track in the secondary and tertiary axes
        const TwoDSlidingFitResult &firstFitResult2D(slidingFitResult.GetFirstFitResult());
        for (const auto &mapEntry : firstFitResult2D.GetLayerFitResultMap())
        {
            rTMinFit1 = std::min(rTMinFit1, static_cast<float>(mapEntry.second.GetFitT()));
            rTMaxFit1 = std::max(rTMaxFit1, static_cast<float>(mapEntry.second.GetFitT()));
        }
        const TwoDSlidingFitResult &secondFitResult2D(slidingFitResult.GetSecondFitResult());
        for (const auto &mapEntry : secondFitResult2D.GetLayerFitResultMap())
        {
            rTMinFit2 = std::min(rTMinFit2, static_cast<float>(mapEntry.second.GetFitT()));
            rTMaxFit2 = std::max(rTMaxFit2, static_cast<float>(mapEntry.second.GetFitT()));
        }
    }
    catch (const StatusCodeException &)
    {
        return false;
    }

    const float deltaTFit1{rTMaxFit1 - rTMinFit1};
    const float deltaTFit2{rTMaxFit2 - rTMinFit2};
    const float deltaT{deltaTFit1 > deltaTFit2 ? deltaTFit1 : deltaTFit2};

    //    std::cout << pCluster << " " << pCluster->GetNCaloHits() << ": " << straightLineLength << " (" << m_maxShowerLengthCut << ") :: "
    //              << integratedPathLength / straightLineLength << " (" << m_pathLengthRatioCut << ") :: "
    //              << deltaT / straightLineLength << " (" << m_rTWidthRatioCut << ")" << std::endl;

    if (straightLineLength < std::numeric_limits<float>::epsilon())
        return false;

    if (straightLineLength > m_maxShowerLengthCut)
        return true;

    if ((integratedPathLength < std::numeric_limits<float>::epsilon()) || (integratedPathLength / straightLineLength > m_pathLengthRatioCut))
        return false;

    if (deltaT / straightLineLength > m_rTWidthRatioCut)
        return false;

    const float vertexDistance(CutClusterCharacterisationAlgorithm::GetVertexDistance(this, pCluster));

    if ((vertexDistance > std::numeric_limits<float>::epsilon()) && ((vertexDistance / straightLineLength) > m_vertexDistanceRatioCut))
        return false;

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CutClusterCharacterisationThreeDAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "SlidingFitWindow", m_slidingFitWindow));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinCaloHitsCut", m_minCaloHitsCut));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MaxShowerLengthCut", m_maxShowerLengthCut));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "PathLengthRatioCut", m_pathLengthRatioCut));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "RTWidthRatioCut", m_rTWidthRatioCut));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "VertexDistanceRatioCut", m_vertexDistanceRatioCut));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "ShowerWidthRatioCut", m_showerWidthRatioCut));

    return ClusterCharacterisationBaseAlgorithm::ReadSettings(xmlHandle);
}

} // namespace lar_content
