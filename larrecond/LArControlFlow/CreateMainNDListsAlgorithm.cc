/**
 *  @file   larrecond/LArControlFlow/CreateMainNDListsAlgorithm.cc
 *
 *  @brief  Implementation of creating various lists for the MainND instance
 *
 *  $Log: $
 */

#include "Api/PandoraApi.h"
#include "Api/PandoraContentApi.h"
#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArControlFlow/MultiPandoraApi.h"

#include "larrecond/LArControlFlow/CreateMainNDListsAlgorithm.h"
using namespace pandora;

namespace lar_content
{

CreateMainNDListsAlgorithm::CreateMainNDListsAlgorithm() :
    m_inputCaloHitListName("Input"),
    m_outputCaloHitListNameU("CaloHitListU"),
    m_outputCaloHitListNameV("CaloHitListV"),
    m_outputCaloHitListNameW("CaloHitListW"),
    m_outputCaloHitListName3D("CaloHitList3D"),
    m_recreatedPfoListName("MainPFOs"),
    m_recreatedClusterListName("MainClusters"),
    m_recreatedVertexListName("MainVertices")
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CreateMainNDListsAlgorithm::Reset()
{
    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CreateMainNDListsAlgorithm::Run()
{
    if (this->GetPandora().GetSettings()->ShouldDisplayAlgorithmInfo())
        std::cout << "Running CreateMainNDListsAlgorithm for Pandora instance " << this->GetPandora().GetName() << std::endl;

    return this->CreateLists();
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CreateMainNDListsAlgorithm::CreateLists()
{
    const CaloHitList *pCaloHitList(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_inputCaloHitListName, pCaloHitList));

    if (pCaloHitList->empty())
        return STATUS_CODE_SUCCESS;

    CaloHitList caloHitListU, caloHitListV, caloHitListW, caloHitList3D;

    for (const CaloHit *const pCaloHit : *pCaloHitList)
    {
        if (TPC_VIEW_U == pCaloHit->GetHitType())
        {
            caloHitListU.push_back(pCaloHit);
        }
        else if (TPC_VIEW_V == pCaloHit->GetHitType())
        {
            caloHitListV.push_back(pCaloHit);
        }
        else if (TPC_VIEW_W == pCaloHit->GetHitType())
        {
            caloHitListW.push_back(pCaloHit);
        }
        else if (TPC_3D == pCaloHit->GetHitType())
        {
            caloHitList3D.push_back(pCaloHit);
        }
    }

    // Save these hit lists
    if (!caloHitListU.empty() && !m_outputCaloHitListNameU.empty())
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList(*this, caloHitListU, m_outputCaloHitListNameU));

    if (!caloHitListV.empty() && !m_outputCaloHitListNameV.empty())
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList(*this, caloHitListV, m_outputCaloHitListNameV));

    if (!caloHitListW.empty() && !m_outputCaloHitListNameW.empty())
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList(*this, caloHitListW, m_outputCaloHitListNameW));

    if (!caloHitList3D.empty() && !m_outputCaloHitListName3D.empty())
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList(*this, caloHitList3D, m_outputCaloHitListName3D));

    // Save PFOs
    return this->CreateMainPFOList();
}

StatusCode CreateMainNDListsAlgorithm::CreateMainPFOList()
{
    // Loop over the Pandora instances, retrieve their PFOs, and store them in a combined, temporary list.
    // This list will then be deep-copied to the main Pandora instance
    PfoList inputPfoList;
    const PandoraInstanceList pandoraList = MultiPandoraApi::GetDaughterPandoraInstanceList(&this->GetPandora());
    for (const auto *pPandora : pandoraList)
    {
        const std::string pandoraName = pPandora->GetName();
        const PfoList *pPfoList{nullptr};
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::GetCurrentPfoList(*pPandora, pPfoList));

        if (!pPfoList)
            continue;

        for (const ParticleFlowObject *pPfo : *pPfoList)
            inputPfoList.emplace_back(pPfo);
    }

    // Save the main PFO list to the main instance by creating new PFO deep copies.
    // Based on LArContent/LArControlFlow/MasterAlgorithm::Recreate(oldPFOs, newPFOs)
    if (this->GetPandora().GetSettings()->ShouldDisplayAlgorithmInfo())
        std::cout << "Main Pandora PFO list size = " << inputPfoList.size() << std::endl;

    PfoList copiedPfoList;
    return this->Recreate(inputPfoList, copiedPfoList);
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CreateMainNDListsAlgorithm::Recreate(const PfoList &inputPfoList, PfoList &newPfoList) const
{
    if (inputPfoList.empty())
        return STATUS_CODE_SUCCESS;

    std::string clusterListName;
    const ClusterList *pClusterList(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::CreateTemporaryListAndSetCurrent(*this, pClusterList, clusterListName));

    std::string vertexListName;
    const VertexList *pVertexList(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::CreateTemporaryListAndSetCurrent(*this, pVertexList, vertexListName));

    std::string pfoListName;
    const PfoList *pPfoList(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::CreateTemporaryListAndSetCurrent(*this, pPfoList, pfoListName));

    for (const Pfo *const pPfo : inputPfoList)
    {
        if (pPfo->GetParentPfoList().empty())
            this->Recreate(pPfo, nullptr, newPfoList);
    }

    if (!pClusterList->empty())
    {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList<Cluster>(*this, m_recreatedClusterListName));
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<Cluster>(*this, m_recreatedClusterListName));
    }

    if (!pVertexList->empty())
    {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList<Vertex>(*this, m_recreatedVertexListName));
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<Vertex>(*this, m_recreatedVertexListName));
    }

    if (!pPfoList->empty())
    {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList<ParticleFlowObject>(*this, m_recreatedPfoListName));
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<ParticleFlowObject>(*this, m_recreatedPfoListName));
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CreateMainNDListsAlgorithm::Recreate(
    const ParticleFlowObject *const pInputPfo, const ParticleFlowObject *const pNewParentPfo, PfoList &newPfoList) const
{
    // Get the 2D & 3D clusters for this given input PFO
    ClusterList inputClusterList2D, inputClusterList3D, newClusterList;
    LArPfoHelper::GetTwoDClusterList(pInputPfo, inputClusterList2D);
    LArPfoHelper::GetThreeDClusterList(pInputPfo, inputClusterList3D);

    // 2D clusters
    for (const Cluster *const pInputCluster : inputClusterList2D)
    {
        CaloHitList inputCaloHitList, newCaloHitList, newIsolatedCaloHitList;
        pInputCluster->GetOrderedCaloHitList().FillCaloHitList(inputCaloHitList);

        // Create hit deep copies
        for (const CaloHit *const pInputCaloHit : inputCaloHitList)
            newCaloHitList.push_back(this->CreateCaloHit(pInputCaloHit));

        for (const CaloHit *const pInputCaloHit : pInputCluster->GetIsolatedCaloHitList())
            newIsolatedCaloHitList.push_back(this->CreateCaloHit(pInputCaloHit));

        // Create cluster deep copies
        if (!newCaloHitList.empty())
            newClusterList.push_back(this->CreateCluster(pInputCluster, newCaloHitList, newIsolatedCaloHitList));
    }

    // 3D clusters
    for (const Cluster *const pInputCluster : inputClusterList3D)
    {
        CaloHitList inputCaloHitList, newCaloHitList, newIsolatedCaloHitList;
        pInputCluster->GetOrderedCaloHitList().FillCaloHitList(inputCaloHitList);

        // Create hit deep copies
        for (const CaloHit *const pInputCaloHit : inputCaloHitList)
            newCaloHitList.push_back(this->CreateCaloHit(pInputCaloHit));

        for (const CaloHit *const pInputCaloHit : pInputCluster->GetIsolatedCaloHitList())
            newIsolatedCaloHitList.push_back(this->CreateCaloHit(pInputCaloHit));

        // Create cluster deep copies
        if (!newCaloHitList.empty())
            newClusterList.push_back(this->CreateCluster(pInputCluster, newCaloHitList, newIsolatedCaloHitList));
    }

    // Create vertex deep copies
    VertexList newVertexList;
    for (const Vertex *const pInputVertex : pInputPfo->GetVertexList())
        newVertexList.push_back(this->CreateVertex(pInputVertex));

    // Create PFO deep copy
    const ParticleFlowObject *const pNewPfo(this->CreatePfo(pInputPfo, newClusterList, newVertexList));
    newPfoList.push_back(pNewPfo);

    if (pNewParentPfo)
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SetPfoParentDaughterRelationship(*this, pNewParentPfo, pNewPfo))

    // Recursive deep copies for daughter PFOs
    for (const ParticleFlowObject *const pInputDaughterPfo : pInputPfo->GetDaughterPfoList())
        this->Recreate(pInputDaughterPfo, pNewPfo, newPfoList);

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

const CaloHit *CreateMainNDListsAlgorithm::CreateCaloHit(const CaloHit *const pInputCaloHit) const
{
    PandoraContentApi::CaloHit::Parameters parameters;
    parameters.m_positionVector = pInputCaloHit->GetPositionVector();
    parameters.m_expectedDirection = pInputCaloHit->GetExpectedDirection();
    parameters.m_cellNormalVector = pInputCaloHit->GetCellNormalVector();
    parameters.m_cellGeometry = pInputCaloHit->GetCellGeometry();
    parameters.m_cellSize0 = pInputCaloHit->GetCellSize0();
    parameters.m_cellSize1 = pInputCaloHit->GetCellSize1();
    parameters.m_cellThickness = pInputCaloHit->GetCellThickness();
    parameters.m_nCellRadiationLengths = pInputCaloHit->GetNCellRadiationLengths();
    parameters.m_nCellInteractionLengths = pInputCaloHit->GetNCellInteractionLengths();
    parameters.m_time = pInputCaloHit->GetTime();
    parameters.m_inputEnergy = pInputCaloHit->GetInputEnergy();
    parameters.m_mipEquivalentEnergy = pInputCaloHit->GetMipEquivalentEnergy();
    parameters.m_electromagneticEnergy = pInputCaloHit->GetElectromagneticEnergy();
    parameters.m_hadronicEnergy = pInputCaloHit->GetHadronicEnergy();
    parameters.m_isDigital = pInputCaloHit->IsDigital();
    parameters.m_hitType = pInputCaloHit->GetHitType();
    parameters.m_hitRegion = pInputCaloHit->GetHitRegion();
    parameters.m_layer = pInputCaloHit->GetLayer();
    parameters.m_isInOuterSamplingLayer = pInputCaloHit->IsInOuterSamplingLayer();
    parameters.m_pParentAddress = pInputCaloHit->GetParentAddress();

    const CaloHit *pNewCaloHit(nullptr);
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::CaloHit::Create(*this, parameters, pNewCaloHit));

    PandoraContentApi::CaloHit::Metadata metadata;
    metadata.m_isIsolated = pInputCaloHit->IsIsolated();
    metadata.m_isPossibleMip = pInputCaloHit->IsPossibleMip();
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::CaloHit::AlterMetadata(*this, pNewCaloHit, metadata));

    return pNewCaloHit;
}

//------------------------------------------------------------------------------------------------------------------------------------------

const Cluster *CreateMainNDListsAlgorithm::CreateCluster(
    const Cluster *const pInputCluster, const CaloHitList &newCaloHitList, const CaloHitList &newIsolatedCaloHitList) const
{
    PandoraContentApi::Cluster::Parameters parameters;
    parameters.m_caloHitList = newCaloHitList;
    parameters.m_isolatedCaloHitList = newIsolatedCaloHitList;

    const Cluster *pNewCluster(nullptr);
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Cluster::Create(*this, parameters, pNewCluster));

    PandoraContentApi::Cluster::Metadata metadata;
    metadata.m_particleId = pInputCluster->GetParticleId();
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Cluster::AlterMetadata(*this, pNewCluster, metadata));

    return pNewCluster;
}

//------------------------------------------------------------------------------------------------------------------------------------------

const Vertex *CreateMainNDListsAlgorithm::CreateVertex(const Vertex *const pInputVertex) const
{
    PandoraContentApi::Vertex::Parameters parameters;
    parameters.m_position = pInputVertex->GetPosition();
    parameters.m_vertexLabel = pInputVertex->GetVertexLabel();
    parameters.m_vertexType = pInputVertex->GetVertexType();

    const Vertex *pNewVertex(nullptr);
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Vertex::Create(*this, parameters, pNewVertex));

    return pNewVertex;
}

//------------------------------------------------------------------------------------------------------------------------------------------

const ParticleFlowObject *CreateMainNDListsAlgorithm::CreatePfo(
    const ParticleFlowObject *const pInputPfo, const ClusterList &newClusterList, const VertexList &newVertexList) const
{
    PandoraContentApi::ParticleFlowObject::Parameters parameters;
    parameters.m_particleId = pInputPfo->GetParticleId();
    parameters.m_charge = pInputPfo->GetCharge();
    parameters.m_mass = pInputPfo->GetMass();
    parameters.m_energy = pInputPfo->GetEnergy();
    parameters.m_momentum = pInputPfo->GetMomentum();
    parameters.m_clusterList = newClusterList;
    parameters.m_trackList.clear();
    parameters.m_vertexList = newVertexList;
    parameters.m_propertiesToAdd = pInputPfo->GetPropertiesMap();

    const ParticleFlowObject *pNewPfo(nullptr);
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ParticleFlowObject::Create(*this, parameters, pNewPfo));

    return pNewPfo;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CreateMainNDListsAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "InputCaloHitListName", m_inputCaloHitListName));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "OutputCaloHitListNameU", m_outputCaloHitListNameU));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "OutputCaloHitListNameV", m_outputCaloHitListNameV));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "OutputCaloHitListNameW", m_outputCaloHitListNameW));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "OutputCaloHitListName3D", m_outputCaloHitListName3D));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "RecreatedPfoListName", m_recreatedPfoListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "RecreatedClusterListName", m_recreatedClusterListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "RecreatedVertexListName", m_recreatedVertexListName));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
