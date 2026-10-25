/**
 *  @file   larpandoracontent/LArVertex/VertexSelectionBaseAlgorithm.h
 *
 *  @brief  Header file for the vertex selection base algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_LArNDVERTEX_SELECTION_BASE_ALGORITHM_H
#define LAR_LArNDVERTEX_SELECTION_BASE_ALGORITHM_H 1

#include "larpandoracontent/LArVertex/VertexSelectionBaseAlgorithm.h"

namespace lar_content
{

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  VertexSelectionBaseAlgorithm class
 */
class LArNDVertexSelectionBaseAlgorithm : public VertexSelectionBaseAlgorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    LArNDVertexSelectionBaseAlgorithm();

protected:
    /**
     *  @brief  Whether the vertex lies on a hit in the specified view
     *
     *  @param  pVertex the address of the vertex
     *  @param  hitType the relevant hit type
     *  @param  kdTree the relevant kd tree
     *
     *  @return boolean
     */
    bool IsVertexOnHit(const pandora::Vertex *const pVertex, const pandora::HitType hitType, HitKDTree2D &kdTree) const;

    /**
     *  @brief  Filter the input list of vertices to obtain a reduced number of vertex candidates
     *
     *  @param  pInputVertexList the address of the input vertex list
     *  @param  kdTreeU the kd tree for u hits
     *  @param  kdTreeV the kd tree for v hits
     *  @param  kdTreeW the kd tree for w hits
     *  @param  filteredVertices to receive the filtered vertex list
     */
    void FilterVertexList(const pandora::VertexList *const pInputVertexList, HitKDTree2D &kdTreeU, HitKDTree2D &kdTreeV,
        HitKDTree2D &kdTreeW, pandora::VertexVector &filteredVertices) const override;

    /**
     *  @brief  Calculate the energy of a vertex candidate by summing values from all three planes
     *
     *  @param  pVertex the address of the vertex
     *  @param  kdTreeMap the map of 2D hit kd trees
     *
     *  @return the summed vertex energy
     */
    float GetVertexEnergy(const pandora::Vertex *const pVertex, const KDTreeMap &kdTreeMap) const;

    /**
     *  @brief  Sort vertices by increasing z position
     *
     *  @param  pLhs address of the lhs vertex
     *  @param  pRhs address of the rhs vertex
     *
     *  @return whether lhs should precedes rhs
     */
    static bool SortByVertexZPosition(const pandora::Vertex *const pLhs, const pandora::Vertex *const pRhs);

private:
    float m_maxOnHitDisplacement;               ///< Max hit-vertex displacement for declaring vertex to lie on a hit in each view

    bool m_isEmptyViewAcceptable;               ///< Whether views entirely empty of hits are classed as 'acceptable' for candidate filtration
    unsigned int m_minVertexAcceptableViews;    ///< The minimum number of views in which a candidate must sit on/near a hit or in a gap (or view can be empty)

};

} // namespace lar_content

#endif // #ifndef LAR_LArNDVERTEX_SELECTION_BASE_ALGORITHM_H
