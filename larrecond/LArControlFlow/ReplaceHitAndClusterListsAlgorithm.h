/**
 *  @file   larrecond/LArControlFlow/ReplaceHitAndClusterListsAlgorithm.h
 *
 *  @brief  Header file for the class to replace current hit and cluster lists.
 *
 *  $Log: $
 */
#ifndef LAR_REPLACE_HIT_AND_CLUSTER_LISTS_ALGORITHM_H
#define LAR_REPLACE_HIT_AND_CLUSTER_LISTS_ALGORITHM_H 1

#include <unordered_map>

#include "Pandora/Algorithm.h"

namespace lar_content
{

/**
 *  @brief  ReplaceHitAndClusterListsAlgorithm class
 */
class ReplaceHitAndClusterListsAlgorithm : public pandora::Algorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    ReplaceHitAndClusterListsAlgorithm();

private:
    pandora::StatusCode Run();

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    std::string m_inputClusterListName; ///< Name of the cluster list
    std::string m_inputCaloHitListName; ///< Name of the calo hit list
};

} // namespace lar_content

#endif // #ifndef LAR_REPLACE_HIT_AND_CLUSTER_LISTS_ALGORITHM_H
