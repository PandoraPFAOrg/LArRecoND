/**
 *  @file   src/SimpleClusterCreationThreeDAlgorithm.cc
 *
 *  @brief  Implementation of the 3D cluster creation algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "SimpleClusterCreationThreeDAlgorithm.h"

using namespace pandora;

namespace lar_content
{

SimpleClusterCreationThreeDAlgorithm::SimpleClusterCreationThreeDAlgorithm()
{
    m_clusteringWindowSquared = 0.25f;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void SimpleClusterCreationThreeDAlgorithm::BuildAssociationMap(const CaloHitVector &caloHitVector, HitAssociationMap &hitAssociationMap) const
{
    this->BuildAssociationMapKDTree<3>(caloHitVector, hitAssociationMap);
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode SimpleClusterCreationThreeDAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, SimpleClusterCreationAlgorithm::ReadSettings(xmlHandle));

    if (!m_inputCaloHitListName.empty() || !m_outputClusterListName.empty())
    {
        std::cout << "SimpleClusterCreationThreeDAlgorithm: 2D calo hit list XML settings are set for 3D clustering!" << std::endl;
        return STATUS_CODE_INVALID_PARAMETER;
    }

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "InputCaloHitListName3D", m_inputCaloHitListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "OutputClusterListName3D", m_outputClusterListName));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
