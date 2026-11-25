/**
 *  @file   larrecond/LArControlFlow/EventSlicingThreeDBaseTool.h
 *
 *  @brief  Header file for the event slicing tool base class.
 *
 *  $Log: $
 */
#ifndef LAR_EVENT_SLICING_THREE_D_BASE_TOOL_H
#define LAR_EVENT_SLICING_THREE_D_BASE_TOOL_H 1

#include "Pandora/Algorithm.h"
#include "Pandora/AlgorithmTool.h"

#include "larrecond/LArObjects/LArSlice3D.h"

namespace lar_content
{

typedef std::map<pandora::HitType, std::string> HitTypeToNameMap;
/**
 *  @brief  EventSlicingThreeDBaseTool class
 */
class EventSlicingThreeDBaseTool : public pandora::AlgorithmTool
{
public:
    /**
     *  @brief  Run the slicing tool
     *
     *  @param  pAlgorithm address of the calling algorithm
     *  @param  caloHitListNames the hit type to calo hit list name map
     *  @param  clusterListNames the hit type to cluster list name map
     *  @param  sliceList to receive the populated slice list
     */
    virtual void RunSlicing(const pandora::Algorithm *const pAlgorithm, const HitTypeToNameMap &caloHitListNames,
        const HitTypeToNameMap &clusterListNames, Slice3DList &sliceList) = 0;
};

} // namespace lar_content

#endif // #ifndef LAR_EVENT_SLICING_THREE_D_BASE_TOOL_H
