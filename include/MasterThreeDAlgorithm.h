/**
 *  @file   include/MasterThreeDAlgorithm.h
 *
 *  @brief  Header file for the master algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_MASTER_THREE_D_ALGORITHM_H
#define LAR_MASTER_THREE_D_ALGORITHM_H 1

#include "Pandora/AlgorithmTool.h"
#include "Pandora/ExternallyConfiguredAlgorithm.h"

#include "larpandoracontent/LArControlFlow/MasterAlgorithm.h"
#include "larpandoracontent/LArControlFlow/MultiPandoraApi.h"
#include "larpandoracontent/LArObjects/LArCaloHit.h"

#include "RockMuonTaggingTool.h"

#include <unordered_map>

namespace lar_content
{

typedef std::unordered_map<unsigned int, std::vector<const pandora::LArTPC *>> WorkerToLArTPCMap;
//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  MasterThreeDAlgorithm class
 */
class MasterThreeDAlgorithm : public MasterAlgorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    MasterThreeDAlgorithm();

protected:
    pandora::StatusCode Run() override;

    /**
     *  @brief  Copy mc particles in the named input list to the given pandora instance.
     *
     *  @param  pPandora the pandora instance to which to copy the mc particles
     */
    pandora::StatusCode CopyMCParticles(const pandora::Pandora *pPandora) const;


    /**
     * @brief  Run cosmic-ray reconstruction, then recreate the resulting pfos in the current pandora instance.
     *
     * @param  volumeIdToHitListMap the LArTPC Volume ID to Hit List Map
     * @param  pfoToLArTPCMap To receive the mapping from pfos
     */
    pandora::StatusCode RunCosmicRayReconstructionThenRecreate(const VolumeIdToHitListMap &volumeIdToHitListMap, PfoToLArTPCMap &pfoToLArTPCMap) const;

    /**
     *  @brief  Tag clear, unambiguous cosmic-ray pfos
     *
     *  @param  stitchedPfosToX0Map a map of cosmic-ray pfos that have been stitched between lar tpcs to the X0 shift
     *  @param  clearCosmicRayPfos to receive the list of clear cosmic-ray pfos
     *  @param  ambiguousPfos to receive the list of ambiguous cosmic-ray pfos for further analysis
     */
    pandora::StatusCode TagCosmicRayPfos(
        const PfoToFloatMap &stitchedPfosToX0Map, pandora::PfoList &clearCosmicRayPfos, pandora::PfoList &ambiguousPfos) const override;

    /**
     *  @brief  Run cosmic-ray hit removal, freeing hits in ambiguous pfos for further processing
     *
     *  @param  ambiguousPfos the list of ambiguous cosmic-ray pfos
     */
    pandora::StatusCode RunCosmicRayHitRemoval(const pandora::PfoList &ambiguousPfos) const;

    /**
     *  @brief  Run the event slicing procedures, dividing available hits up into distinct 3D regions
     *
     *  @param  volumeIdToHitListMap the volume id to hit list map
     *  @param  sliceVector to receive the populated slice vector
     *
     *  @return whether slicing could be run
     */
    pandora::StatusCode RunSlicing(const VolumeIdToHitListMap &volumeIdToHitListMap, SliceVector &sliceVector) const;

    /**
     *  @brief  Recreate a specified pfo in the current pandora instance
     *
     *  @param  pInputPfo the input pfo
     *  @param  pNewParentPfo the new parent of the new output pfo (nullptr if none)
     *  @param  newPfoList to receive the list of new pfos
     */
    pandora::StatusCode Recreate(const pandora::ParticleFlowObject *const pInputPfo, const pandora::ParticleFlowObject *const pNewParentPfo,
        pandora::PfoList &newPfoList) const override;

    /**
     *  @brief  Create a pandora worker instance to handle a single LArTPC
     *
     *  @param  larTPC the lar tpc
     *  @param  gapList the gap list
     *  @param  settingsFile the pandora settings file
     *  @param  name the pandora instance name
     *
     *  @return the address of the pandora instance
     */
    const pandora::Pandora *CreateWorkerInstance(const pandora::LArTPC &larTPC, const pandora::DetectorGapList &gapList,
        const std::string &settingsFile, const std::string &name) const;

    /**
     *  @brief  Create a pandora worker instance to handle a number of LArTPCs
     *
     *  @param  larTPCMap the lar tpc map
     *  @param  gapList the gap list
     *  @param  settingsFile the pandora settings file
     *  @param  name the pandora instance name
     *  @param  id for the created worker instance (i.e. larTPCParameters.m_larTPCVolumeId)   
     *
     *  @return the address of the pandora instance
     */
    const pandora::Pandora *CreateWorkerInstance(const pandora::LArTPCMap &larTPCMap, const pandora::DetectorGapList &gapList,
        const std::string &settingsFile, const std::string &name, const unsigned int id) const;

    /**
     *  @brief  Initialize pandora worker instances
     *
     *  @param workerToLArTPCMap to map each worker instance to the list of TPCs it acts on
     */
    pandora::StatusCode InitializeWorkerInstances(WorkerToLArTPCMap &workerToLArTPCMap);

    /**
     *  @brief  Get the mapping from lar tpc volume id to lists of all hits, and truncated hits
     *
     *  @param  volumeIdToHitListMap to receive the populated volume id to hit list map
     *
     *  @return status code
     */
    pandora::StatusCode GetVolumeIdToHitListMap(VolumeIdToHitListMap &volumeIdToHitListMap) const;

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle) override;

    typedef std::vector<RockMuonTaggingTool *> RockMuonTaggingToolVector;

    bool m_shouldRunRockMus_Xworkers;                      ///< Whether to run rock muons reconstruction using a columnar X worker
    bool m_tagRockMuons;                                   ///< bool to activate tagging of rock muons
    RockMuonTaggingToolVector m_rockMuonTaggingToolVector; ///< The cosmic-ray tagging tool vector
    WorkerToLArTPCMap m_workerToLArTPCMap;                 ///< mapping between worker instances and LArTPCs
};

} // namespace lar_content

#endif // #ifndef LAR_MASTER_THREE_D_ALGORITHM_H
