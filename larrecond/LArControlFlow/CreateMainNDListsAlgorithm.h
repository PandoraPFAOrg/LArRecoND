/**
 *  @file   larrecond/LArControlFlow/CreateMainNDListsAlgorithm.h
 *
 *  @brief  Header file for creating various lists for the MainND instance
 *
 *  $Log: $
 */
#ifndef LAR_CREATE_MAIN_ND_LISTS_ALGORITHM_H
#define LAR_CREATE_MAIN_ND_LISTS_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

namespace lar_content
{

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  CreateMainNDListsAlgorithm class
 */
class CreateMainNDListsAlgorithm : public pandora::Algorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    CreateMainNDListsAlgorithm();

private:
    pandora::StatusCode Reset();
    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    /**
     *  @brief Create the various lists
     */
    pandora::StatusCode CreateLists();

    /**
     *  @brief Create the main PFO list
     */
    pandora::StatusCode CreateMainPFOList();

    /**
     *  @brief  Recreate a specified list of pfos in the main pandora instance
     *
     *  @param  inputPfoList the input pfo list
     *  @param  newPfoList to receive the list of new pfos
     */
    pandora::StatusCode Recreate(const pandora::PfoList &inputPfoList, pandora::PfoList &newPfoList) const;

    /**
     *  @brief  Recreate a specified pfo in the main pandora instance
     *
     *  @param  pInputPfo the input pfo
     *  @param  pNewParentPfo the new parent of the new output pfo (nullptr if none)
     *  @param  newPfoList to receive the list of new pfos
     */
    pandora::StatusCode Recreate(const pandora::ParticleFlowObject *const pInputPfo, const pandora::ParticleFlowObject *const pNewParentPfo,
        pandora::PfoList &newPfoList) const;

    /**
     *  @brief  Create a new calo hit in the main pandora instance, based upon the provided input calo hit
     *
     *  @param  pInputCaloHit the address of the input calo hit
     *  @param  pParentCaloHit the address of the parent calo hit
     *
     *  @return the address of the new calo hit
     */
    const pandora::CaloHit *CreateCaloHit(const pandora::CaloHit *const pInputCaloHit) const;

    /**
     *  @brief  Create a new cluster in the current pandora instance, based upon the provided input cluster
     *
     *  @param  pInputCluster the address of the input cluster
     *  @param  newCaloHitList the list of calo hits for the new cluster
     *  @param  newIsolatedCaloHitList the list of isolated calo hits for the new cluster
     *
     *  @return the address of the new cluster
     */
    const pandora::Cluster *CreateCluster(const pandora::Cluster *const pInputCluster, const pandora::CaloHitList &newCaloHitList,
        const pandora::CaloHitList &newIsolatedCaloHitList) const;

    /**
     *  @brief  Create a new vertex in the main pandora instance, based upon the provided input vertex
     *
     *  @param  pInputVertex the address of the input vertex
     *
     *  @return the address of the new vertex
     */
    const pandora::Vertex *CreateVertex(const pandora::Vertex *const pInputVertex) const;

    /**
     *  @brief  Create a new pfo in the main pandora instance, based upon the provided input pfo
     *
     *  @param  pInputPfo the address of the input pfo
     *  @param  newClusterList the list of clusters for the new pfo
     *  @param  newVertexList the list of vertices for the new pfo
     *
     *  @return the address of the new pfo
     */
    const pandora::ParticleFlowObject *CreatePfo(const pandora::ParticleFlowObject *const pInputPfo,
        const pandora::ClusterList &newClusterList, const pandora::VertexList &newVertexList) const;

    std::string m_inputCaloHitListName;     ///< The input calo hit list name
    std::string m_outputCaloHitListNameU;   ///< The output calo hit list name for TPC_VIEW_U hits
    std::string m_outputCaloHitListNameV;   ///< The output calo hit list name for TPC_VIEW_V hits
    std::string m_outputCaloHitListNameW;   ///< The output calo hit list name for TPC_VIEW_W hits
    std::string m_outputCaloHitListName3D;  ///< The output calo hit list name for TPC_3D hits
    std::string m_recreatedPfoListName;     ///< The output recreated pfo list name
    std::string m_recreatedClusterListName; ///< The output recreated cluster list name
    std::string m_recreatedVertexListName;  ///< The output recreated vertex list name
};

} // namespace lar_content

#endif // #ifndef LAR_PRE_PROCESSING_THREE_D_ALGORITHM_H
