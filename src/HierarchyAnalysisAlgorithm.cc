/**
 *  @file   src/HierarchyAnalysisAlgorithm.cc
 *
 *  @brief  Implementation of the hierarchy analysis output algorithm
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "HierarchyAnalysisAlgorithm.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

#include "LArT0Helper.h"

#include "TFile.h"
#include "TTree.h"

using namespace pandora;

namespace lar_content
{

HierarchyAnalysisAlgorithm::HierarchyAnalysisAlgorithm() :
    m_count{-1},
    m_event{-1},
    m_run{0},
    m_subRun{0},
    m_unixTime{0},
    m_unixTimeUsec{0},
    m_startTime{0},
    m_endTime{0},
    m_triggers{0},
    m_nhits{0},
    m_mcIDs{nullptr},
    m_mcLocalIDs{nullptr},
    m_eventFileName{""},
    m_eventTreeName{"events"},
    m_eventLeafName{"event"},
    m_runLeafName{"run"},
    m_subRunLeafName{"subrun"},
    m_unixTimeLeafName{"unix_ts"},
    m_unixTimeUsecLeafName{"unix_ts_usec"},
    m_startTimeLeafName{"event_start_t"},
    m_endTimeLeafName{"event_end_t"},
    m_triggersLeafName{"triggers"},
    m_nhitsLeafName{"nhits"},
    m_mcIdLeafName{"mcp_id"},
    m_mcLocalIdLeafName{"mcp_idLocal"},
    m_eventsToSkip{0},
    m_minHitsToSkip{2},
    m_eventFile{nullptr},
    m_eventTree{nullptr},
    m_caloHitListName{"CaloHitList2D"},
    m_pfoListName{"RecreatedPfos"},
    m_minTrackScore{0.5f},
    m_analysisFileName{"LArRecoND.root"},
    m_analysisTreeName{"LArRecoND"},
    m_foldToPrimaries{false},
    m_foldToLeadingShowers{false},
    m_foldDynamic{true},
    m_minPurity{0.5f},
    m_minCompleteness{0.1f},
    m_minRecoHits{15},
    m_minRecoHitsPerView{5},
    m_minRecoGoodViews{2},
    m_removeRecoNeutrons{true},
    m_selectRecoHits{true},
    m_storeClusterRecoHits{true},
    m_gotMCEventInput{false},
    m_mcIdMap{}
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

HierarchyAnalysisAlgorithm::~HierarchyAnalysisAlgorithm()
{
    // Save the analysis output ROOT file. Always recreate this
    PANDORA_MONITORING_API(SaveTree(this->GetPandora(), m_analysisTreeName.c_str(), m_analysisFileName.c_str(), "RECREATE"));

    // Cleanup ROOT file used for the event numbers
    if (m_eventFile && m_eventFile->IsOpen())
    {
        delete m_eventTree;
        m_eventTree = nullptr;
    }
    delete m_eventFile;
    m_eventFile = nullptr;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode HierarchyAnalysisAlgorithm::Run()
{
    // Increment the algorithm run count
    ++m_count;

    // Set the event run number and trigger timing info, as well as the unique-local MCParticle Id map
    this->SetEventRunMCIdInfo();

    // Need to use 2D calo hit list for now since LArHierarchyHelper::MCHierarchy::IsReconstructable()
    // checks for minimum number of hits in the U, V & W views only, which will fail for 3D
    const CaloHitList *pCaloHitList(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_caloHitListName, pCaloHitList));
    const MCParticleList *pMCParticleList(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pMCParticleList));
    const PfoList *pPfoList(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_pfoListName, pPfoList));

    LArHierarchyHelper::FoldingParameters foldParameters;
    if (m_foldToPrimaries)
        foldParameters.m_foldToTier = true;
    else if (m_foldToLeadingShowers)
        foldParameters.m_foldToLeadingShowers = true;
    else if (m_foldDynamic)
        foldParameters.m_foldDynamic = true;

    const LArHierarchyHelper::MCHierarchy::ReconstructabilityCriteria recoCriteria(
        m_minRecoHits, m_minRecoHitsPerView, m_minRecoGoodViews, m_removeRecoNeutrons);

    LArHierarchyHelper::MCHierarchy mcHierarchy(recoCriteria);
    LArHierarchyHelper::FillMCHierarchy(*pMCParticleList, *pCaloHitList, foldParameters, mcHierarchy);
    LArHierarchyHelper::RecoHierarchy recoHierarchy;
    LArHierarchyHelper::FillRecoHierarchy(*pPfoList, foldParameters, recoHierarchy);

    const LArHierarchyHelper::QualityCuts quality(m_minPurity, m_minCompleteness, m_selectRecoHits);
    LArHierarchyHelper::MatchInfo matchInfo(mcHierarchy, recoHierarchy, quality);
    LArHierarchyHelper::MatchHierarchies(matchInfo);
    matchInfo.Print(mcHierarchy);

    // Analysis PFO & matched reco-MC output
    this->EventAnalysisOutput(matchInfo);

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void HierarchyAnalysisAlgorithm::SetEventRunMCIdInfo()
{
    // Set the event and run numbers as well as the trigger timing.
    // Also fill the map linking the local & unique MC particle Ids.
    // MCParticles use unique MC Ids, but CAFs need the local ones as well
    m_mcIdMap.clear();

    if (m_eventTree)
    {
        // Sets m_event, m_run, m_subRun, m_unixTime, m_unixTimeUsec, m_startTime, m_endTime & m_triggers
        const int iEntry = m_count + m_eventsToSkip;
        m_eventTree->GetEntry(iEntry);

        // Check if we should actually be pointing to a higher event number due to skipping some events in Pandora
        if (m_nhits < m_minHitsToSkip)
        {
            // Skip ahead as many events as we need to
            int thisHits = m_nhits;
            while (thisHits < m_minHitsToSkip)
            {
                m_count += 1;
                const int newEntry = m_count + m_eventsToSkip;
                m_eventTree->GetEntry(newEntry);
                thisHits = m_nhits;
            }
        }

        // Fill the Id map
        if (m_gotMCEventInput)
        {
            for (size_t i = 0; i < m_mcIDs->size(); i++)
                m_mcIdMap[(*m_mcIDs)[i]] = (*m_mcLocalIDs)[i];
        }
    }
    else
        // Use the algorithm run count number
        m_event = m_count;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void HierarchyAnalysisAlgorithm::EventAnalysisOutput(const LArHierarchyHelper::MatchInfo &matchInfo) const
{
    // For storing various reconstructed PFO quantities in the given event
    int sliceId{-1};
    // Slice & cluster IDs, and number of hits
    IntVector sliceIdVect, clusterIdVect, n3DHitsVect, nUHitsVect, nVHitsVect, nWHitsVect;
    // Cluster isShower, isRecoPrimary & reco PDG hypothesis, as well as the track score
    IntVector isShowerVect, isRecoPrimaryVect, recoPDGVect;
    FloatVector trackScoreVect;
    // Reco neutrino vertex
    FloatVector nuVtxXVect, nuVtxYVect, nuVtxZVect;
    // Cluster start, end, direction, PCA axis lengths and total hit energy
    FloatVector startXVect, startYVect, startZVect, endXVect, endYVect, endZVect;
    FloatVector t0Vect;
    FloatVector dirXVect, dirYVect, dirZVect, centroidXVect, centroidYVect, centroidZVect;
    FloatVector primaryLVect, secondaryLVect, tertiaryLVect, energyVect;
    // Best matched MC info
    IntVector matchVect, mcPDGVect, nSharedHitsVect, isPrimaryVect;
    FloatVector completenessVect, purityVect;
    // MC matched energy, momentum, vertex and end position
    FloatVector mcEVect, mcPxVect, mcPyVect, mcPzVect;
    FloatVector mcVtxXVect, mcVtxYVect, mcVtxZVect, mcEndXVect, mcEndYVect, mcEndZVect;
    // MC neutrino parent info
    IntVector mcNuPDGVect, mcNuCodeVect;
    FloatVector mcNuVtxXVect, mcNuVtxYVect, mcNuVtxZVect;
    FloatVector mcNuEVect, mcNuPxVect, mcNuPyVect, mcNuPzVect;
    // Long integers for the MC IDs: vertex, unique and local trajectories
    std::vector<long> mcNuIdVect, mcIdVect, mcLocalIdVect, mcParentIdVect;
    //MC pfo parent info
    IntVector mcParentPDGVect;

    // Hit info for each reconstructed PFO. Since we can't store vectors of vectors, the size of
    // these vectors = n3DHits*nPFOs, whereas all of the above vectors have size = nPFOs.
    // The hit vector entries follow the order PFO1[n3DHits1], PFO2[n3DHits2], PFO3[n3DHits3] etc.
    // The sliceId & clusterId vectors keep track of where a given hit comes from
    IntVector recoHitIdVect, recoHitSliceIdVect, recoHitClusterIdVect;
    FloatVector recoHitXVect, recoHitYVect, recoHitZVect, recoHitEVect, recoHitT0Vect;

    // Get the list of root MCParticles for the MC truth matching
    MCParticleList rootMCParticles;
    matchInfo.GetRootMCParticles(rootMCParticles);

    // Get reconstructed root PFOs (neutrinos)
    PfoList rootPfos;
    const LArHierarchyHelper::RecoHierarchy &recoHierarchy{matchInfo.GetRecoHierarchy()};
    recoHierarchy.GetRootPfos(rootPfos);

    // Loop over the root PFOs
    for (const ParticleFlowObject *const pRoot : rootPfos)
    {
        // Slice id = root PFO number
        ++sliceId;

        // Get (first) root vertex
        const VertexList &rootVertices{pRoot->GetVertexList()};
        const Vertex *pRootVertex = (rootVertices.size() > 0) ? (*rootVertices.begin()) : nullptr;
        const float max{std::numeric_limits<float>::max()};
        const CartesianVector rootRecoVtx = (pRootVertex != nullptr) ? pRootVertex->GetPosition() : CartesianVector(max, max, max);

        // Get reco nodes for each root PFO
        LArHierarchyHelper::RecoHierarchy::NodeVector recoNodes;
        recoHierarchy.GetFlattenedNodes(pRoot, recoNodes);

        // Cluster id for the given slice
        int clusterId{-1};

        // Loop over the reco nodes
        for (const LArHierarchyHelper::RecoHierarchy::Node *pRecoNode : recoNodes)
        {
            // Get the list of PFOs for each node
            const PfoList recoParticles = pRecoNode->GetRecoParticles();

            // Is this a primary node (hierarchy tier = 1)?
            const int isRecoPrimary = (pRecoNode->GetHierarchyTier() == 1) ? 1 : 0;

            // Get individual PFOs
            for (const ParticleFlowObject *pPfo : recoParticles)
            {
                // Get the 3D cluster
                const Cluster *pCluster3D = this->GetCluster(pPfo, TPC_3D);
                if (!pCluster3D)
                    continue;

                // Make sure the cluster has some hits
                const int n3DHits(pCluster3D->GetNCaloHits());
                if (n3DHits <= 1)
                    continue;

                // Increment clusterId
                clusterId++;

                // Find first and last cluster hit points
                CartesianVector first(max, max, max), last(max, max, max);
                LArClusterHelper::GetExtremalCoordinates(pCluster3D, first, last);
                // Get the (first) vertex if it exists, otherwise use the first hit position
                const VertexList &vertices{pPfo->GetVertexList()};
                const CartesianVector vertex = (vertices.size() > 0) ? (*vertices.begin())->GetPosition() : first;

                // Principal component analysis of the cluster
                CartesianPointVector pointVector;
                LArClusterHelper::GetCoordinateVector(pCluster3D, pointVector);
                // Use cluster local vertex for relative axis directions
                const LArShowerPCA pca = LArPfoHelper::GetPrincipalComponents(pointVector, vertex);

                // Centroid, primary axis direction and lengths
                const CartesianVector centroid{pca.GetCentroid()};
                const CartesianVector primaryAxis{pca.GetPrimaryAxis()};
                const float primaryLength{pca.GetPrimaryLength()};
                const float secondaryLength{pca.GetSecondaryLength()};
                const float tertiaryLength{pca.GetTertiaryLength()};

                // Cluster start is assumed to be the vertex. Set the end point as either the
                // first or last hit point that has the largest squared distance from the vertex
                const float firstDistSq = first.GetDistanceSquared(vertex);
                const float lastDistSq = last.GetDistanceSquared(vertex);
                const CartesianVector endPoint = (lastDistSq > firstDistSq) ? last : first;

                // Use the primary axis to set the direction.
                // Reverse this if (endPoint - vertex) dot primaryAxis < 0
                const CartesianVector displacement{endPoint - vertex};
                const CartesianVector direction = (displacement.GetDotProduct(primaryAxis) < 0.0) ? primaryAxis * (-1.0) : primaryAxis;

                // Cluster deposited energy (all hits assume EM energy = hadronic energy)
                const float clusterEnergy{pCluster3D->GetElectromagneticEnergy()};

                // Get number of hits in the U, V and W views (if they exist)
                const Cluster *pClusterU = this->GetCluster(pPfo, TPC_VIEW_U);
                const int nUHits = (pClusterU != nullptr) ? pClusterU->GetNCaloHits() : 0;
                const Cluster *pClusterV = this->GetCluster(pPfo, TPC_VIEW_V);
                const int nVHits = (pClusterV != nullptr) ? pClusterV->GetNCaloHits() : 0;
                const Cluster *pClusterW = this->GetCluster(pPfo, TPC_VIEW_W);
                const int nWHits = (pClusterW != nullptr) ? pClusterW->GetNCaloHits() : 0;

                // Store quantities in the vectors
                sliceIdVect.emplace_back(sliceId);
                // Neutrino reco vertex
                nuVtxXVect.emplace_back(rootRecoVtx.GetX());
                nuVtxYVect.emplace_back(rootRecoVtx.GetY());
                nuVtxZVect.emplace_back(rootRecoVtx.GetZ());

                // Cluster Id
                clusterIdVect.emplace_back(clusterId);

                // Number of hits in the cluster (by views)
                n3DHitsVect.emplace_back(n3DHits);
                nUHitsVect.emplace_back(nUHits);
                nVHitsVect.emplace_back(nVHits);
                nWHitsVect.emplace_back(nWHits);

                // Save the track score, getting the appropriate metadata
                // see e.g. https://github.com/PandoraPFA/larpandora/blob/develop/larpandora/LArPandoraInterface/LArPandoraOutput.cxx#L325 for similar
                const auto &properties = pPfo->GetPropertiesMap();
                float trackScore{-1.f};
                const auto iterTrackScore(properties.find("TrackScore"));
                if (iterTrackScore != properties.end())
                {
                    trackScore = iterTrackScore->second;
                }
                trackScoreVect.emplace_back(trackScore);

                // Define isShower based on track score
                const int isShower = (trackScore >= m_minTrackScore) ? 0 : 1;
                isShowerVect.emplace_back(isShower);

                // Set reco PDG hypothesis, e.g track = muon, shower = electron.
                // Since all PFOs are tracks for now, this will always be muon
                const int recoPDG = (isShower == 0) ? MU_MINUS : E_MINUS;
                recoPDGVect.emplace_back(recoPDG);

                // Is this a reconstructed primary PFO?
                isRecoPrimaryVect.emplace_back(isRecoPrimary);

                // Cluster vertex, end and direction (from PCA)
                startXVect.emplace_back(vertex.GetX());
                startYVect.emplace_back(vertex.GetY());
                startZVect.emplace_back(vertex.GetZ());
                endXVect.emplace_back(endPoint.GetX());
                endYVect.emplace_back(endPoint.GetY());
                endZVect.emplace_back(endPoint.GetZ());
                dirXVect.emplace_back(direction.GetX());
                dirYVect.emplace_back(direction.GetY());
                dirZVect.emplace_back(direction.GetZ());
                // Cluster centroid and axis lengths (from PCA)
                centroidXVect.emplace_back(centroid.GetX());
                centroidYVect.emplace_back(centroid.GetY());
                centroidZVect.emplace_back(centroid.GetZ());
                primaryLVect.emplace_back(primaryLength);
                secondaryLVect.emplace_back(secondaryLength);
                tertiaryLVect.emplace_back(tertiaryLength);

                // Cluster energy (sum over all hits)
                energyVect.emplace_back(clusterEnergy);

		pandora::CaloHitList calo3DHitList;
		LArClusterHelper::GetAllHits(pCluster3D, calo3DHitList);
		float clusterT0 = LArT0Helper::GetMedianTime(calo3DHitList);
		t0Vect.emplace_back(clusterT0);

                if (m_storeClusterRecoHits)
                {
                    // Store 3D reco hit information for this cluster/PFO. Vector sizes = nPFOs*n3DHits not nPFOs.
                    // The hit vector entries follow the order PFO1[n3DHits1], PFO2[n3DHits2], PFO3[n3DHits3] etc.
                    // Sort hits using their positions
                    calo3DHitList.sort(LArClusterHelper::SortHitsByPosition);

                    for (const auto *pCalo3DHit : calo3DHitList)
                    {
                        const int hitId = reinterpret_cast<intptr_t>(pCalo3DHit->GetParentAddress());
                        const CartesianVector hitPos = pCalo3DHit->GetPositionVector();
                        const float hitE = pCalo3DHit->GetInputEnergy();
                        recoHitIdVect.emplace_back(hitId);
                        recoHitSliceIdVect.emplace_back(sliceId);
                        recoHitClusterIdVect.emplace_back(clusterId);
                        recoHitXVect.emplace_back(hitPos.GetX());
                        recoHitYVect.emplace_back(hitPos.GetY());
                        recoHitZVect.emplace_back(hitPos.GetZ());
                        recoHitEVect.emplace_back(hitE);
			recoHitT0Vect.emplace_back(pCalo3DHit->GetTime());
                    }
                }

                // Find best-matched MC particle for this reconstructed PFO cluster
                const HierarchyAnalysisAlgorithm::RecoMCMatch bestMatch = GetRecoMCMatch(pPfo, pRecoNode, matchInfo, rootMCParticles);

                // Best matched MC particle
                const MCParticle *pLeadingMC = bestMatch.m_pLeadingMC;
                const int gotMatch = (pLeadingMC != nullptr) ? 1 : 0;
                const int mcPDG = (pLeadingMC != nullptr) ? pLeadingMC->GetParticleId() : 0;
                // Unique and local MC Ids
                const long mcId = (pLeadingMC != nullptr) ? reinterpret_cast<intptr_t>(pLeadingMC->GetUid()) : 0;
                const long mcLocalId = (m_mcIdMap.find(mcId) != m_mcIdMap.end()) ? m_mcIdMap.at(mcId) : mcId;
                const int isPrimary = (pLeadingMC != nullptr && LArMCParticleHelper::IsPrimary(pLeadingMC)) ? 1 : 0;
                const float mcEnergy = (pLeadingMC != nullptr) ? pLeadingMC->GetEnergy() : 0.f;
                const CartesianVector mcMomentum = (pLeadingMC != nullptr) ? pLeadingMC->GetMomentum() : CartesianVector(0.f, 0.f, 0.f);
                const CartesianVector mcVertex = (pLeadingMC != nullptr) ? pLeadingMC->GetVertex() : CartesianVector(max, max, max);
                const CartesianVector mcEndPoint = (pLeadingMC != nullptr) ? pLeadingMC->GetEndpoint() : CartesianVector(max, max, max);

                // Retrieve MC parent info
                int mcParentPDG{-999}, mcParentId{-999};
                if (pLeadingMC)
                {
                    const MCParticleList &parentList{pLeadingMC->GetParentList()};
                    if (!parentList.empty())
                    {
                        const MCParticle *pParent{parentList.front()};
                        mcParentPDG = (pParent != nullptr) ? pParent->GetParticleId() : -999;
                        mcParentId = (pParent != nullptr) ? reinterpret_cast<intptr_t>(pParent->GetUid()) : -999;
                    }
                }

                mcParentPDGVect.emplace_back(mcParentPDG);
                mcParentIdVect.emplace_back(mcParentId);

                // MC neutrino parent info, including Nuance interaction code
                const MCParticle *pNuRoot = bestMatch.m_pNuRoot;
                const int mcNuPDG = (pNuRoot != nullptr) ? pNuRoot->GetParticleId() : 0;
                // Neutrino Id = unique vertex Id
                const long mcNuId = (pNuRoot != nullptr) ? reinterpret_cast<intptr_t>(pNuRoot->GetUid()) : 0;
                const int mcNuCode = (dynamic_cast<const LArMCParticle *>(pNuRoot) != nullptr) ? LArMCParticleHelper::GetNuanceCode(pNuRoot) : 0;
                const CartesianVector mcNuVertex = (pNuRoot != nullptr) ? pNuRoot->GetVertex() : CartesianVector(max, max, max);
                const float mcNuEnergy = (pNuRoot != nullptr) ? pNuRoot->GetEnergy() : 0.f;
                const CartesianVector mcNuMomentum = (pNuRoot != nullptr) ? pNuRoot->GetMomentum() : CartesianVector(0.f, 0.f, 0.f);

                matchVect.emplace_back(gotMatch);
                mcPDGVect.emplace_back(mcPDG);
                mcIdVect.emplace_back(mcId);
                mcLocalIdVect.emplace_back(mcLocalId);
                isPrimaryVect.emplace_back(isPrimary);
                nSharedHitsVect.emplace_back(bestMatch.m_nSharedHits);
                completenessVect.emplace_back(bestMatch.m_completeness);
                purityVect.emplace_back(bestMatch.m_purity);
                mcEVect.emplace_back(mcEnergy);
                mcPxVect.emplace_back(mcMomentum.GetX());
                mcPyVect.emplace_back(mcMomentum.GetY());
                mcPzVect.emplace_back(mcMomentum.GetZ());
                mcVtxXVect.emplace_back(mcVertex.GetX());
                mcVtxYVect.emplace_back(mcVertex.GetY());
                mcVtxZVect.emplace_back(mcVertex.GetZ());
                mcEndXVect.emplace_back(mcEndPoint.GetX());
                mcEndYVect.emplace_back(mcEndPoint.GetY());
                mcEndZVect.emplace_back(mcEndPoint.GetZ());
                mcNuPDGVect.emplace_back(mcNuPDG);
                mcNuIdVect.emplace_back(mcNuId);
                mcNuCodeVect.emplace_back(mcNuCode);
                mcNuVtxXVect.emplace_back(mcNuVertex.GetX());
                mcNuVtxYVect.emplace_back(mcNuVertex.GetY());
                mcNuVtxZVect.emplace_back(mcNuVertex.GetZ());
                mcNuEVect.emplace_back(mcNuEnergy);
                mcNuPxVect.emplace_back(mcNuMomentum.GetX());
                mcNuPyVect.emplace_back(mcNuMomentum.GetY());
                mcNuPzVect.emplace_back(mcNuMomentum.GetZ());

            } // Reco PFOs
        } // Reco nodes
    } // Root PFOs

    // Fill ROOT ntuple
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_analysisTreeName.c_str(), "event", m_event));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_analysisTreeName.c_str(), "run", m_run));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_analysisTreeName.c_str(), "subRun", m_subRun));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_analysisTreeName.c_str(), "unixTime", m_unixTime));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_analysisTreeName.c_str(), "unixTimeUsec", m_unixTimeUsec));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_analysisTreeName.c_str(), "startTime", m_startTime));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_analysisTreeName.c_str(), "endTime", m_endTime));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_analysisTreeName.c_str(), "triggers", m_triggers));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_analysisTreeName.c_str(), "sliceId", &sliceIdVect));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_analysisTreeName.c_str(), "nuVtxX", &nuVtxXVect));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_analysisTreeName.c_str(), "nuVtxY", &nuVtxYVect));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_analysisTreeName.c_str(), "nuVtxZ", &nuVtxZVect));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_analysisTreeName.c_str(), "clusterId", &clusterIdVect));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_analysisTreeName.c_str(), "n3DHits", &n3DHitsVect));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_analysisTreeName.c_str(), "nUHits", &nUHitsVect));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_analysisTreeName.c_str(), "nVHits", &nVHitsVect));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_analysisTreeName.c_str(), "nWHits", &nWHitsVect));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_analysisTreeName.c_str(), "isShower", &isShowerVect));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_analysisTreeName.c_str(), "trackScore", &trackScoreVect));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_analysisTreeName.c_str(), "recoPDG", &recoPDGVect));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_analysisTreeName.c_str(), "isRecoPrimary", &isRecoPrimaryVect));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_analysisTreeName.c_str(), "startX", &startXVect));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_analysisTreeName.c_str(), "startY", &startYVect));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_analysisTreeName.c_str(), "startZ", &startZVect));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_analysisTreeName.c_str(), "endX", &endXVect));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_analysisTreeName.c_str(), "endY", &endYVect));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_analysisTreeName.c_str(), "endZ", &endZVect));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_analysisTreeName.c_str(), "t0", &t0Vect));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_analysisTreeName.c_str(), "dirX", &dirXVect));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_analysisTreeName.c_str(), "dirY", &dirYVect));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_analysisTreeName.c_str(), "dirZ", &dirZVect));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_analysisTreeName.c_str(), "centroidX", &centroidXVect));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_analysisTreeName.c_str(), "centroidY", &centroidYVect));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_analysisTreeName.c_str(), "centroidZ", &centroidZVect));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_analysisTreeName.c_str(), "length1", &primaryLVect));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_analysisTreeName.c_str(), "length2", &secondaryLVect));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_analysisTreeName.c_str(), "length3", &tertiaryLVect));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_analysisTreeName.c_str(), "energy", &energyVect));
    if (m_storeClusterRecoHits)
    {
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_analysisTreeName.c_str(), "recoHitId", &recoHitIdVect));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_analysisTreeName.c_str(), "recoHitSliceId", &recoHitSliceIdVect));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_analysisTreeName.c_str(), "recoHitClusterId", &recoHitClusterIdVect));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_analysisTreeName.c_str(), "recoHitX", &recoHitXVect));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_analysisTreeName.c_str(), "recoHitY", &recoHitYVect));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_analysisTreeName.c_str(), "recoHitZ", &recoHitZVect));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_analysisTreeName.c_str(), "recoHitE", &recoHitEVect));
	PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_analysisTreeName.c_str(), "recoHitT0", &recoHitT0Vect));
    }
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_analysisTreeName.c_str(), "gotMatch", &matchVect));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_analysisTreeName.c_str(), "mcPDG", &mcPDGVect));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_analysisTreeName.c_str(), "mcId", &mcIdVect));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_analysisTreeName.c_str(), "mcLocalId", &mcLocalIdVect));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_analysisTreeName.c_str(), "isPrimary", &isPrimaryVect));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_analysisTreeName.c_str(), "nSharedHits", &nSharedHitsVect));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_analysisTreeName.c_str(), "completeness", &completenessVect));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_analysisTreeName.c_str(), "purity", &purityVect));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_analysisTreeName.c_str(), "mcEnergy", &mcEVect));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_analysisTreeName.c_str(), "mcPx", &mcPxVect));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_analysisTreeName.c_str(), "mcPy", &mcPyVect));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_analysisTreeName.c_str(), "mcPz", &mcPzVect));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_analysisTreeName.c_str(), "mcVtxX", &mcVtxXVect));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_analysisTreeName.c_str(), "mcVtxY", &mcVtxYVect));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_analysisTreeName.c_str(), "mcVtxZ", &mcVtxZVect));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_analysisTreeName.c_str(), "mcEndX", &mcEndXVect));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_analysisTreeName.c_str(), "mcEndY", &mcEndYVect));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_analysisTreeName.c_str(), "mcEndZ", &mcEndZVect));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_analysisTreeName.c_str(), "mcNuPDG", &mcNuPDGVect));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_analysisTreeName.c_str(), "mcNuId", &mcNuIdVect));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_analysisTreeName.c_str(), "mcNuCode", &mcNuCodeVect));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_analysisTreeName.c_str(), "mcNuVtxX", &mcNuVtxXVect));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_analysisTreeName.c_str(), "mcNuVtxY", &mcNuVtxYVect));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_analysisTreeName.c_str(), "mcNuVtxZ", &mcNuVtxZVect));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_analysisTreeName.c_str(), "mcNuE", &mcNuEVect));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_analysisTreeName.c_str(), "mcNuPx", &mcNuPxVect));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_analysisTreeName.c_str(), "mcNuPy", &mcNuPyVect));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_analysisTreeName.c_str(), "mcNuPz", &mcNuPzVect));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_analysisTreeName.c_str(), "mcParentPDG", &mcParentPDGVect));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_analysisTreeName.c_str(), "mcParentId", &mcParentIdVect));

    PANDORA_MONITORING_API(FillTree(this->GetPandora(), m_analysisTreeName.c_str()));
}

//------------------------------------------------------------------------------------------------------------------------------------------

const Cluster *HierarchyAnalysisAlgorithm::GetCluster(const ParticleFlowObject *pPfo, const HitType hitType) const
{
    // Predicate for finding the specific cluster type from a given PFO's list of clusters (3D or 2D views)
    const auto predicate = [&hitType](const Cluster *pCluster) { return LArClusterHelper::GetClusterHitType(pCluster) == hitType; };

    const ClusterList &clusters{pPfo->GetClusterList()};
    const Cluster *pCluster{nullptr};
    ClusterList::const_iterator iter = std::find_if(clusters.begin(), clusters.end(), predicate);
    if (iter != clusters.end())
        pCluster = *iter;
    return pCluster;
}

//------------------------------------------------------------------------------------------------------------------------------------------

const HierarchyAnalysisAlgorithm::RecoMCMatch HierarchyAnalysisAlgorithm::GetRecoMCMatch(const pandora::ParticleFlowObject *pPfo,
    const LArHierarchyHelper::RecoHierarchy::Node *pRecoNode, const LArHierarchyHelper::MatchInfo &matchInfo, MCParticleList &rootMCParticles) const
{
    int nSharedHits{0};
    float completeness{0.f}, purity{0.f};

    const MCParticle *pRootNu{nullptr}, *pLeadingMC{nullptr};

    // Loop over the root (neutrino) MC particles
    for (const MCParticle *const pMCRoot : rootMCParticles)
    {
        // Loop over the possible matches
        const LArHierarchyHelper::MCMatchesVector &matches{matchInfo.GetMatches(pMCRoot)};

        for (const LArHierarchyHelper::MCMatches &match : matches)
        {
            // MC node
            const LArHierarchyHelper::MCHierarchy::Node *pMCNode{match.GetMC()};

            // Reco matches
            const LArHierarchyHelper::RecoHierarchy::NodeVector &nodeVector{match.GetRecoMatches()};

            // See if the current recoNode is in the reco matches vector
            if (std::find(nodeVector.begin(), nodeVector.end(), pRecoNode) != nodeVector.end())
            {
                // The MC matching uses nodes which can have more than 1 folded particle/PFO.
                // The PFO passed to this function will be part of the matched reco node
                // once the code reaches here, but this reco node can have more than 1 PFO via folding.
                // For completeness & purity, only use the required PFO hits and not all hits from
                // the reco node (which will include other PFOs if they are present).
                // Copies and adapts code from GetSharedHits(), GetCompleteness() & GetPurity() from
                // LArHierarchyHelper::MCMatches to avoid duplicate shared hits, completeness & purities
                // for the case when reco nodes have more than 1 PFO. Otherwise each PFO from
                // the same reco node will always return the same match quantities

                // Get hits for this best matched MC node
                CaloHitList mcHits = pMCNode->GetCaloHits();

                // Get all the reco hits for the passed PFO only, not the whole reco node
                CaloHitList pfoHits;
                LArPfoHelper::GetAllCaloHits2D(pPfo, pfoHits);

                // Get the selected hits of the reco node, such that each reco hit has a corresponding MC hit
                CaloHitList selectedRecoHits = match.GetSelectedRecoHits(pRecoNode);

                // Sort the hit lists
                mcHits.sort();
                pfoHits.sort();
                selectedRecoHits.sort();

                // Get the PFO reco hits that are also part of the selected reco list (such that each reco hit has a corresponding MC hit)
                CaloHitVector selectedPfoHits;
                std::set_intersection(
                    pfoHits.begin(), pfoHits.end(), selectedRecoHits.begin(), selectedRecoHits.end(), std::back_inserter(selectedPfoHits));

                // Find the overlap of these MC & selected PFO reco hit lists
                CaloHitVector intersection;
                std::set_intersection(mcHits.begin(), mcHits.end(), selectedPfoHits.begin(), selectedPfoHits.end(), std::back_inserter(intersection));

                const int currentSharedHits = intersection.size();

                if (currentSharedHits > nSharedHits)
                {
                    nSharedHits = currentSharedHits;
                    completeness = (mcHits.size() > 0) ? nSharedHits / static_cast<float>(mcHits.size()) : 0.0;
                    purity = (selectedPfoHits.size() > 0) ? nSharedHits / static_cast<float>(selectedPfoHits.size()) : 0.0;
                    pRootNu = pMCRoot;
                    pLeadingMC = pMCNode->GetLeadingMCParticle();
                }
            } // Find recoNode
        } // Match loop
    } // Root MC particles

    const HierarchyAnalysisAlgorithm::RecoMCMatch info(pRootNu, pLeadingMC, nSharedHits, completeness, purity);
    return info;
}

//------------------------------------------------------------------------------------------------------------------------------------------

HierarchyAnalysisAlgorithm::RecoMCMatch::RecoMCMatch(const pandora::MCParticle *pNuRoot, const pandora::MCParticle *pLeadingMC,
    const int nSharedHits, const float completeness, const float purity) :
    m_pNuRoot(pNuRoot),
    m_pLeadingMC(pLeadingMC),
    m_nSharedHits(nSharedHits),
    m_completeness(completeness),
    m_purity(purity)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode HierarchyAnalysisAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "EventFileName", m_eventFileName));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "EventTreeName", m_eventTreeName));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "EventLeafName", m_eventLeafName));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "RunLeafName", m_runLeafName));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "SubRunLeafName", m_subRunLeafName));
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "UnixTimeLeafName", m_unixTimeLeafName));
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "UnixTimeUsecLeafName", m_unixTimeUsecLeafName));
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "StartTimeLeafName", m_startTimeLeafName));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "EndTimeLeafName", m_endTimeLeafName));
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "TriggersLeafName", m_triggersLeafName));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "NHitsLeafName", m_nhitsLeafName));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MCIdLeafName", m_mcIdLeafName));
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MCLocalIdLeafName", m_mcLocalIdLeafName));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "EventsToSkip", m_eventsToSkip));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinHitsToSkip", m_minHitsToSkip));

    // Setup the event ROOT file
    if (m_eventFileName.size() > 0)
    {
        m_eventFile = TFile::Open(m_eventFileName.c_str(), "READ");
        if (m_eventFile && m_eventFile->IsOpen())
        {
            m_eventTree = dynamic_cast<TTree *>(m_eventFile->Get(m_eventTreeName.c_str()));
            if (m_eventTree)
            {
                // Only enable the event and run number leaves as well as the trigger timing.
                // Also enable the vertex_id leaf
                m_eventTree->SetBranchStatus("*", 0);
                m_eventTree->SetBranchStatus(m_eventLeafName.c_str(), 1);
                m_eventTree->SetBranchStatus(m_runLeafName.c_str(), 1);
                m_eventTree->SetBranchStatus(m_subRunLeafName.c_str(), 1);
                m_eventTree->SetBranchStatus(m_unixTimeLeafName.c_str(), 1);
                m_eventTree->SetBranchStatus(m_unixTimeUsecLeafName.c_str(), 1);
                m_eventTree->SetBranchStatus(m_startTimeLeafName.c_str(), 1);
                m_eventTree->SetBranchStatus(m_endTimeLeafName.c_str(), 1);
                m_eventTree->SetBranchStatus(m_triggersLeafName.c_str(), 1);
                m_eventTree->SetBranchStatus(m_nhitsLeafName.c_str(), 1);
                m_eventTree->SetBranchAddress(m_eventLeafName.c_str(), &m_event);
                m_eventTree->SetBranchAddress(m_runLeafName.c_str(), &m_run);
                m_eventTree->SetBranchAddress(m_subRunLeafName.c_str(), &m_subRun);
                m_eventTree->SetBranchAddress(m_unixTimeLeafName.c_str(), &m_unixTime);
                m_eventTree->SetBranchAddress(m_unixTimeUsecLeafName.c_str(), &m_unixTimeUsec);
                m_eventTree->SetBranchAddress(m_startTimeLeafName.c_str(), &m_startTime);
                m_eventTree->SetBranchAddress(m_endTimeLeafName.c_str(), &m_endTime);
                m_eventTree->SetBranchAddress(m_triggersLeafName.c_str(), &m_triggers);
                m_eventTree->SetBranchAddress(m_nhitsLeafName.c_str(), &m_nhits);

                // Check if we have MC branches
                if (m_eventTree->GetBranch(m_mcIdLeafName.c_str()) && m_eventTree->GetBranch(m_mcLocalIdLeafName.c_str()))
                {
                    m_gotMCEventInput = true;
                    m_eventTree->SetBranchStatus(m_mcIdLeafName.c_str(), 1);
                    m_eventTree->SetBranchStatus(m_mcLocalIdLeafName.c_str(), 1);
                    m_eventTree->SetBranchAddress(m_mcIdLeafName.c_str(), &m_mcIDs);
                    m_eventTree->SetBranchAddress(m_mcLocalIdLeafName.c_str(), &m_mcLocalIDs);
                }
            }
        }
    }

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "CaloHitListName", m_caloHitListName));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "PfoListName", m_pfoListName));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinTrackScore", m_minTrackScore));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "AnalysisFileName", m_analysisFileName));
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "AnalysisTreeName", m_analysisTreeName));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "FoldToPrimaries", m_foldToPrimaries));
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "FoldToLeadingShowers", m_foldToLeadingShowers));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "FoldDynamic", m_foldDynamic));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinPurity", m_minPurity));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinCompleteness", m_minCompleteness));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinRecoHits", m_minRecoHits));
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinRecoHitsPerView", m_minRecoHitsPerView));
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinRecoGoodViews", m_minRecoGoodViews));
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "RemoveRecoNeutrons", m_removeRecoNeutrons));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "SelectRecoHits", m_selectRecoHits));
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "StoreClusterRecoHits", m_storeClusterRecoHits));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
