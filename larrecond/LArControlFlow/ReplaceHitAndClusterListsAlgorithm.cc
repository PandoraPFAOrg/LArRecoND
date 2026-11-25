/**
 *  @file   larrecond/LArControlFlow/ReplaceHitAndClusterListsAlgorithm.cc
 *
 *  @brief  Implementation of the cluster creation algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larrecond/LArControlFlow/ReplaceHitAndClusterListsAlgorithm.h"

using namespace pandora;

namespace lar_content
{

ReplaceHitAndClusterListsAlgorithm::ReplaceHitAndClusterListsAlgorithm()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ReplaceHitAndClusterListsAlgorithm::Run()
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<Cluster>(*this, m_inputClusterListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<CaloHit>(*this, m_inputCaloHitListName));

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ReplaceHitAndClusterListsAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "InputClusterListName", m_inputClusterListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "InputCaloHitListName", m_inputCaloHitListName));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
