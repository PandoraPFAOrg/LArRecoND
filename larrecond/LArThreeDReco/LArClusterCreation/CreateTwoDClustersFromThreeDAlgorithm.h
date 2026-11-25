/**
 *  @file   larrecond/LArThreeDReco/LArClusterCreation/CreateTwoDClustersFromThreeDAlgorithm.h
 *
 *  @brief  Header file for the 3D to 2D cluster creation alg
 *
 *  $Log: $
 */
#ifndef LAR_CREATE_TWO_D_CLUSTERS_FROM_THREE_D_ALGORITHM_H
#define LAR_CREATE_TWO_D_CLUSTERS_FROM_THREE_D_ALGORITHM_H 1

#include <unordered_map>

#include "Api/PandoraContentApi.h"
#include "Pandora/Algorithm.h"

namespace lar_content
{

/**
 *  @brief  CreateTwoDClustersFromThreeDAlgorithm class
 */
class CreateTwoDClustersFromThreeDAlgorithm : public pandora::Algorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    CreateTwoDClustersFromThreeDAlgorithm();

private:
    pandora::StatusCode Run();

    typedef std::unordered_map<const pandora::CaloHit *, pandora::CaloHitList> HitAssociationMap;
    typedef std::unordered_map<const pandora::CaloHit *, bool> HitUsedMap;
    typedef std::map<std::pair<int, int>, const pandora::CaloHit *> SpatialHitMap;

    /**
     *  @brief Get the hits in each view matching to the clustered 3D hits
     *
     *  @param pCaloHit3D a pointer to the 3D hit
     *  @param hitMap a spatial map of 2D hits
     *  @param associatedHits reference to empty hit list to be filled with 2D hits
     *  @param usedHits2D list of already used 2D hits
     *  @param hitType the type of hits in the hit list
     */
    void GetAssociatedTwoDHit(const pandora::CaloHit *const pCaloHit3D, const SpatialHitMap &hitMap, pandora::CaloHitList &associatedHits,
        HitUsedMap &usedHits2D, const pandora::HitType &hitType) const;

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    std::vector<std::string> m_inputCaloHitListNames2D; ///< Names of the input 2D hit lists
    std::string m_inputClusterListName3D;               ///< Name of the input 3D cluster list
    std::vector<std::string> m_outputClusterListNames;  ///< Names of the output 2D cluster lists
};

} // namespace lar_content

#endif // #ifndef LAR_CREATE_TWO_D_CLUSTERS_FROM_THREE_D_ALGORITHM_H
