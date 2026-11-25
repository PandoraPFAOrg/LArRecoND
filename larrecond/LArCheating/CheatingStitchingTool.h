/**
 *  @file   larrecond/LArCheating/CheatingStitchingTool.h
 *
 *  @brief  Header file for the cheating stitching tool class.
 *
 *  $Log: $
 */
#ifndef LAR_CHEATING_STITCHING_TOOL_H
#define LAR_CHEATING_STITCHING_TOOL_H 1

#include "larpandoracontent/LArThreeDReco/LArPfoStitching/StitchingBaseTool.h"
#include "larpandoracontent/LArThreeDReco/LArPfoStitching/StitchingPfoOperations.h"

namespace lar_content
{

/**
 *  @brief  CheatingStitchingTool class
 */
class CheatingStitchingTool : public StitchingBaseTool
{
public:
    /**
     *  @brief  Default constructor
     */
    CheatingStitchingTool();

protected:
    /**
     *  @brief  Run the cheating stitching tool, stitching together pfos from different LArTPCs based on MCParticle associations.
     *
     *  @param  pAlgorithm address of the calling algorithm
     *  @param  pStitchingOperations address of the calling algorithm's stitching operations implementation
     *  @param  pMultiPfoList the input pfo list to be stitched (unused)
     *  @param  pfoToLArTPCMap to receive the map from pfos to LArTPCs
     *  @param  stitchedPfosToX0Map to receive the map from stitched pfos to X0 positions (unused)
     */
    void RunStitching(const pandora::Algorithm *const pAlgorithm, const StitchingPfoOperations *const pStitchingOperations,
        const pandora::PfoList *const pMultiPfoList, PfoToLArTPCMap &pfoToLArTPCMap, PfoToFloatMap &stitchedPfosToX0Map);

private:
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);
};

} // namespace lar_content

#endif // #ifndef LAR_CHEATING_STITCHING_TOOL_H
