/**
 *  @file   include/SimpleClusterCreationThreeDAlgorithm.h
 *
 *  @brief  Header file for the 3D cluster creation algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_SIMPLE_CLUSTER_CREATION_THREE_D_ALGORITHM_H
#define LAR_SIMPLE_CLUSTER_CREATION_THREE_D_ALGORITHM_H 1

#include "larpandoracontent/LArTwoDReco/LArClusterCreation/SimpleClusterCreationAlgorithm.h"

namespace lar_content
{

/**
 *  @brief  SimpleClusterCreationThreeDAlgorithm class: cluster creation for a named list of 3D hits, searching for nearby hits in (x,y,z)
 */
class SimpleClusterCreationThreeDAlgorithm : public SimpleClusterCreationAlgorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    SimpleClusterCreationThreeDAlgorithm();

private:
    void BuildAssociationMap(const pandora::CaloHitVector &caloHitVector, HitAssociationMap &hitAssociationMap) const;

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);
};

} // namespace lar_content

#endif // #ifndef LAR_SIMPLE_CLUSTER_CREATION_THREE_D_ALGORITHM_H
