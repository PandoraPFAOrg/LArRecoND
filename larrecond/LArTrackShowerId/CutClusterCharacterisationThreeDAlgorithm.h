/**
 *  @file   larrecond/LArTrackShowerId/CutClusterCharacterisationThreeDAlgorithm.h
 *
 *  @brief  Header file for the cut based cluster characterisation algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_CUT_CLUSTER_CHARACTERISATION_THREE_D_ALGORITHM_H
#define LAR_CUT_CLUSTER_CHARACTERISATION_THREE_D_ALGORITHM_H 1

#include "larpandoracontent/LArTrackShowerId/ClusterCharacterisationBaseAlgorithm.h"

namespace lar_content
{

/**
 *  @brief  CutClusterCharacterisationThreeDAlgorithm class
 */
class CutClusterCharacterisationThreeDAlgorithm : public ClusterCharacterisationBaseAlgorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    CutClusterCharacterisationThreeDAlgorithm();

private:
    /**
     *  @brief  Whether cluster is identified as a clear track
     *
     *  @param  pCluster address of the relevant cluster
     *
     *  @return whether the cluster is a clear track
     */
    virtual bool IsClearTrack(const pandora::Cluster *const pCluster) const;

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    unsigned int m_slidingFitWindow; ///< The layer window for the sliding linear fits
    unsigned int m_minCaloHitsCut;   ///< The minimum number of calo hits to qualify as a track
    float m_maxShowerLengthCut;      ///< The maximum cluster length to qualify as a shower
    float m_pathLengthRatioCut;      ///< The maximum ratio of path length to straight line length to qualify as a track
    float m_rTWidthRatioCut;         ///< The maximum ratio of transverse fit position width to straight line length to qualify as a track
    float m_vertexDistanceRatioCut;  ///< The maximum ratio of vertex separation to straight line length to qualify as a track
    float m_showerWidthRatioCut;     ///< The maximum ratio of shower fit width to straight line length to qualify as a track
};

} // namespace lar_content

#endif // #ifndef LAR_CUT_CLUSTER_CHARACTERISATION_THREE_D_ALGORITHM_H
