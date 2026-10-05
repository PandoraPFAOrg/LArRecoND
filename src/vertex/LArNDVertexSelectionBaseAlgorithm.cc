/**
 *  @file   larpandoracontent/LArVertex/VertexSelectionBaseAlgorithm.cc
 *
 *  @brief  Implementation of the vertex selection base algorithm class.
 *
 *  $Log: $
 */
#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"

#include "larpandoracontent/LArUtility/KDTreeLinkerAlgoT.h"

#include "LArNDVertexSelectionBaseAlgorithm.h"

using namespace pandora;

namespace lar_content
{

LArNDVertexSelectionBaseAlgorithm::LArNDVertexSelectionBaseAlgorithm() :
    m_maxOnHitDisplacement(1.f),
    m_isEmptyViewAcceptable(true),
    m_minVertexAcceptableViews(3)
{
}

void LArNDVertexSelectionBaseAlgorithm::FilterVertexList(const VertexList *const pInputVertexList, HitKDTree2D &kdTreeU, HitKDTree2D &kdTreeV,
    HitKDTree2D &kdTreeW, VertexVector &filteredVertices) const
{
    for (const Vertex *const pVertex : *pInputVertexList)
    {
        unsigned int nAcceptableViews(0);

        if ((m_isEmptyViewAcceptable && kdTreeU.empty()) || this->IsVertexOnHit(pVertex, TPC_VIEW_U, kdTreeU))
            ++nAcceptableViews;

        if ((m_isEmptyViewAcceptable && kdTreeV.empty()) || this->IsVertexOnHit(pVertex, TPC_VIEW_V, kdTreeV))
            ++nAcceptableViews;

        if ((m_isEmptyViewAcceptable && kdTreeW.empty()) || this->IsVertexOnHit(pVertex, TPC_VIEW_W, kdTreeW))
            ++nAcceptableViews;

        if (nAcceptableViews >= m_minVertexAcceptableViews)
            filteredVertices.push_back(pVertex);
    }

    std::sort(filteredVertices.begin(), filteredVertices.end(), SortByVertexZPosition);
}

//------------------------------------------------------------------------------------------------------------------------------------------

float LArNDVertexSelectionBaseAlgorithm::GetVertexEnergy(const pandora::Vertex *const pVertex, const KDTreeMap &kdTreeMap) const
{
    float totalEnergy(0.f);

    if (this->IsVertexOnHit(pVertex, TPC_VIEW_U, kdTreeMap.at(TPC_VIEW_U)))
        totalEnergy += this->VertexHitEnergy(pVertex, TPC_VIEW_U, kdTreeMap.at(TPC_VIEW_U));

    if (this->IsVertexOnHit(pVertex, TPC_VIEW_V, kdTreeMap.at(TPC_VIEW_V)))
        totalEnergy += this->VertexHitEnergy(pVertex, TPC_VIEW_V, kdTreeMap.at(TPC_VIEW_V));

    if (this->IsVertexOnHit(pVertex, TPC_VIEW_W, kdTreeMap.at(TPC_VIEW_W)))
        totalEnergy += this->VertexHitEnergy(pVertex, TPC_VIEW_W, kdTreeMap.at(TPC_VIEW_W));

    return totalEnergy;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArNDVertexSelectionBaseAlgorithm::IsVertexOnHit(const Vertex *const pVertex, const HitType hitType, HitKDTree2D &kdTree) const
{
    const CartesianVector vertexPosition2D(LArGeometryHelper::ProjectPosition(this->GetPandora(), pVertex->GetPosition(), hitType));
    KDTreeBox searchRegionHits = build_2d_kd_search_region(vertexPosition2D, m_maxOnHitDisplacement, m_maxOnHitDisplacement);

    HitKDNode2DList found;
    kdTree.search(searchRegionHits, found);

    return (!found.empty());
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArNDVertexSelectionBaseAlgorithm::SortByVertexZPosition(const pandora::Vertex *const pLhs, const pandora::Vertex *const pRhs)
{
    const CartesianVector deltaPosition(pRhs->GetPosition() - pLhs->GetPosition());

    if (std::fabs(deltaPosition.GetZ()) > std::numeric_limits<float>::epsilon())
        return (deltaPosition.GetZ() > std::numeric_limits<float>::epsilon());

    if (std::fabs(deltaPosition.GetX()) > std::numeric_limits<float>::epsilon())
        return (deltaPosition.GetX() > std::numeric_limits<float>::epsilon());

    // ATTN No way to distinguish between vertices if still have a tie in y coordinate
    return (deltaPosition.GetY() > std::numeric_limits<float>::epsilon());
}


} // namespace lar_content
