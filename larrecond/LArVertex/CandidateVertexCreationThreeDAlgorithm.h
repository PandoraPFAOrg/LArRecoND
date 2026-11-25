/**
 *  @file   larrecond/LArVertex/CandidateVertexCreationThreeDAlgorithm.h
 *
 *  @brief  Header file for the candidate vertex creation algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_CANDIDATE_VERTEX_CREATION_THREE_D_ALGORITHM_H
#define LAR_CANDIDATE_VERTEX_CREATION_THREE_D_ALGORITHM_H 1

#include <unordered_map>

#include "Pandora/Algorithm.h"

#include "larpandoracontent/LArObjects/LArThreeDSlidingFitResult.h"

namespace lar_content
{

/**
 *  @brief  CandidateVertexCreationThreeDAlgorithm::Algorithm class
 */
class CandidateVertexCreationThreeDAlgorithm : public pandora::Algorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    CandidateVertexCreationThreeDAlgorithm();

private:
    pandora::StatusCode Run();

    /**
     *  @brief  Select a subset of input clusters (contained in the input list names) for processing in this algorithm
     *
     *  @param  clusterVector3D to receive the selected 3D clusters
     */
    void SelectClusters(pandora::ClusterVector &clusterVector3D);

    /**
     *  @brief  Create candidate vertex positions by looking at the extreme positions of the clusters
     *
     *  @param  clusterVector3D the 3D clusters
     */
    void CreateEndpointVertices(const pandora::ClusterVector &clusterVector3D) const;

    /**
     *  @brief  Find crossing points between 3D clusters
     *
     *  @param  clusterVector3D the 3D clusters
     */
    void CreateCrossingCandidates(const pandora::ClusterVector &clusterVector3D) const;

    /**
     *  @brief  Identify where clusters plausibly cross in 3D
     *
     *  @param  clusterVector3D the input 3D clusters
     *  @param  crossingPoints to receive the 3D crossing points
     */
    void FindCrossingPoints(const pandora::ClusterVector &clusterVector3D, pandora::CartesianPointVector &crossingPoints) const;

    /**
     *  @brief  Get a list of spacepoints representing cluster 3D hit positions
     *
     *  @param  pCluster address of the cluster
     *  @param  spacePoints to receive the list of spacepoints
     */
    void GetSpacepoints(const pandora::Cluster *const pCluster, pandora::CartesianPointVector &spacePoints) const;

    /**
     *  @brief  Identify where clusters plausibly cross in 3D
     *
     *  @param  spacepoints1 space points for cluster 1
     *  @param  spacepoints2 space points for cluster 2
     *  @param  crossingPoints to receive the list of plausible 3D crossing points
     */
    void FindCrossingPoints(const pandora::CartesianPointVector &spacepoints1, const pandora::CartesianPointVector &spacepoints2,
        pandora::CartesianPointVector &crossingPoints) const;

    /**
     *  @brief  Attempt to create candidate vertex positions, using 3D crossing points
     *
     *  @param  crossingPoints3D the crossing points in 3D
     *  @param  nCrossingCandidates to count the number of crossing candidates created
     */
    void CreateCrossingVertices(const pandora::CartesianPointVector &crossingPoints3D, unsigned int &nCrossingCandidates) const;

    /**
     *  @brief  Add candidate vertices from any input vertices
     */
    void AddInputVertices() const;

    /**
     *  @brief  Creates a 2D sliding fit of a cluster and stores it for later use
     *
     *  @param  pCluster address of the relevant cluster
     */
    void AddToSlidingFitCache(const pandora::Cluster *const pCluster);

    /**
     *  @brief  Get a sliding fit result from the algorithm cache
     *
     *  @param  pCluster address of the relevant cluster
     *
     *  @return 3d sliding fit result object
     */
    const ThreeDSlidingFitResult &GetCachedSlidingFitResult(const pandora::Cluster *const pCluster) const;

    /**
     *  @brief  Clear relevant algorithm member variables between events
     */
    void TidyUp();

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    typedef std::unordered_map<const pandora::Cluster *, pandora::CartesianPointVector> ClusterToSpacepointsMap;

    std::string m_inputClusterListName; ///< The list of cluster list name
    std::string m_inputVertexListName;  ///< The list name for existing candidate vertices
    std::string m_outputVertexListName; ///< The name under which to save the output vertex list
    bool m_replaceCurrentVertexList;    ///< Whether to replace the current vertex list with the output list

    unsigned int m_slidingFitWindow;                 ///< The layer window for the sliding linear fits
    ThreeDSlidingFitResultMap m_slidingFitResultMap; ///< The sliding fit result map

    unsigned int m_minClusterCaloHits; ///< The min number of hits in base cluster selection method
    float m_minClusterLengthSquared;   ///< The min length (squared) in base cluster selection method
    float m_chiSquaredCut;             ///< The chi squared cut (accept only 3D vertex positions with values below cut)

    bool m_enableEndpointCandidates; ///< Whether to create endpoint-based candidates
    float m_maxEndpointXDiscrepancy; ///< The max cluster endpoint discrepancy

    bool m_enableCrossingCandidates;          ///< Whether to create crossing vertex candidates
    unsigned int m_nMaxCrossingCandidates;    ///< The max number of crossing candidates to create
    float m_maxCrossingDiscrepancy;           ///< The max cluster endpoint discrepancy
    unsigned int m_extrapolationNSteps;       ///< Number of extrapolation steps, at each end of cluster, of specified size
    float m_extrapolationStepSize;            ///< The extrapolation step size in cm
    float m_maxCrossingSeparationSquared;     ///< The separation (squared) between spacepoints below which a crossing can be identified
    float m_minNearbyCrossingDistanceSquared; ///< The minimum allowed distance between identified crossing positions

    bool m_reducedCandidates;           ///< Whether to reduce the number of candidates
    float m_selectionCutFactorMax;      ///< Maximum factor to multiply the base cluster selection cuts
    float m_nClustersPassingMaxCutsPar; ///< Parameter for number of clusters passing the max base cluster selection cuts
};

} // namespace lar_content

#endif // #ifndef LAR_CANDIDATE_VERTEX_CREATION_THREE_D_ALGORITHM_H
