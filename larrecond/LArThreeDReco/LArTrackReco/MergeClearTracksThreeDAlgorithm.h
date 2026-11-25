/**
 *  @file   larrecond/LArThreeDReco/LArTrackReco/MergeClearTracksThreeDAlgorithm.h
 *
 *  @brief  Header file for the clear 3D track merging algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_MERGE_CLEAR_TRACKS_THREE_D_ALGORITHM_H
#define LAR_MERGE_CLEAR_TRACKS_THREE_D_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include "larpandoracontent/LArObjects/LArThreeDSlidingFitResult.h"

namespace lar_content
{

/**
 *  @brief  MergeClearTracksThreeDAlgorithm class
 */
class MergeClearTracksThreeDAlgorithm : public pandora::Algorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    MergeClearTracksThreeDAlgorithm();

    pandora::StatusCode Run();

private:
    typedef std::map<const pandora::Cluster *const, float> ClusterFloatMap;
    typedef std::map<const pandora::Cluster *const, const pandora::Cluster *> ClusterMergeMap;

    /**
     *  @brief  Look for possible merges between hit-sorted track-like clusters
     *  @param  pClusterList the list of clusters
     *
     *  @return whether we have found possible cluster merges
     */
    bool FindMerges(const pandora::ClusterList *const pClusterList) const;

    /**
     *  @brief  Try to merge small clusters with large ones depending on their relative distances
     *  @param  pLargeCluster the large cluster
     *  @param  pSmallCluster the small cluster (reduced number of hits)
     *  @param  mergeCandidates large-small cluster pair that can be merged
     *  @param  mergeDistances large-small cluster separation within distance and angle cuts
     */
    void CanMergeClusters(const pandora::Cluster *const pLargeCluster, const pandora::Cluster *const pSmallCluster,
        ClusterMergeMap &mergeCandidates, ClusterFloatMap &mergeDistances) const;

    /**
     *  @brief  Merge the small cluster with the large one
     *  @param  pLargeCluster the large cluster (with the greater number of hits)
     *  @param  pSmallCluster the smaller cluster
     */
    void MergeClusters(const pandora::Cluster *const pLargeCluster, const pandora::Cluster *const pSmallCluster) const;

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    float m_slidingFitWindow;           ///< Number of layers in the sliding fit window
    float m_maxGapLengthCut;            ///< The maximum distance between to tracks to allow merging
    float m_maxGapTransverseCut;        ///< The maximum transverse distance between tracks to allow merging
    float m_minCosThetaCut;             ///< The minimum angle cosine between track segment directions
    std::string m_inputClusterListName; ///< The name of the input 3D cluster list
};

} // namespace lar_content

#endif // #ifndef LAR_MERGE_CLEAR_TRACKS_THREE_D_ALGORITHM_H
