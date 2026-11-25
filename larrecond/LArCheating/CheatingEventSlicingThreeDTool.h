/**
 *  @file   larrecond/LArCheating/CheatingEventSlicingThreeDTool.h
 *
 *  @brief  Header file for the 3D event slicing tool class.
 *
 *  $Log: $
 */
#ifndef LAR_CHEATING_EVENT_SLICING_THREE_D_TOOL_H
#define LAR_CHEATING_EVENT_SLICING_THREE_D_TOOL_H 1

#include <unordered_map>

#include "Pandora/Algorithm.h"
#include "Pandora/AlgorithmTool.h"

#include "larpandoracontent/LArObjects/LArPointingCluster.h"
#include "larpandoracontent/LArObjects/LArThreeDSlidingConeFitResult.h"

#include "larrecond/LArControlFlow/EventSlicingThreeDBaseTool.h"
#include "larrecond/LArControlFlow/SlicingThreeDAlgorithm.h"
#include "larrecond/LArObjects/LArSlice3D.h"

namespace lar_content
{

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  CheatingEventSlicingThreeDTool class
 */
class CheatingEventSlicingThreeDTool : public EventSlicingThreeDBaseTool
{
public:
    /**
     *  @brief  Default constructor
     */
    CheatingEventSlicingThreeDTool();

    /**
     *  @brief  Run the 3D slicing tool
     *
     *  @param  pAlgorithm the address of the calling algorithm
     *  @param  caloHitListNames the hit type to calo hit list name map
     *  @param  clusterListNames the hit type to cluster list name map
     *  @param  sliceList to receive the populated slice list
     */
    void RunSlicing(const pandora::Algorithm *const pAlgorithm, const HitTypeToNameMap &caloHitListNames,
        const HitTypeToNameMap &clusterListNames, Slice3DList &sliceList);

private:
    typedef std::unordered_map<const pandora::MCParticle *, Slice3D> MCParticleToSlice3DMap;

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    /**
     *  @brief  Initialize the map from parent mc particles to slice objects
     *
     *  @param  pAlgorithm address of the calling algorithm
     *  @param  caloHitListNames the hit type to calo hit list name map
     *  @param  mcParticleToSliceMap to receive the parent mc particle to slice map
     */
    void InitializeMCParticleToSlice3DMap(const pandora::Algorithm *const pAlgorithm, const HitTypeToNameMap &caloHitListNames,
        MCParticleToSlice3DMap &mcParticleToSliceMap) const;

    /**
     *  @brief  Fill slices using hits from a specified view
     *
     *  @param  pAlgorithm address of the calling algorithm
     *  @param  hitType the hit type (i.e. view)
     *  @param  caloHitListNames the hit type to calo hit list name map
     *  @param  mcParticleToSliceMap to receive the parent mc particle to slice map
     */
    void FillSlices(const pandora::Algorithm *const pAlgorithm, const pandora::HitType hitType, const HitTypeToNameMap &caloHitListNames,
        MCParticleToSlice3DMap &mcParticleToSliceMap) const;
};

} // namespace lar_content

#endif // #ifndef LAR_CHEATING_EVENT_SLICING_THREE_D_TOOL_H
