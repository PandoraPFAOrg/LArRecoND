/**
 *  @file   larrecond/LArControlFlow/SlicingThreeDAlgorithm.h
 *
 *  @brief  Header file for the 3D slicing algorithm.
 *
 *  $Log: $
 */
#ifndef LAR_SLICING_THREE_D_ALGORITHM_H
#define LAR_SLICING_THREE_D_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include "larpandoracontent/LArControlFlow/SlicingAlgorithm.h"

#include "larrecond/LArControlFlow/EventSlicingThreeDBaseTool.h"

namespace lar_content
{

typedef std::map<pandora::HitType, std::string> HitTypeToNameMap;

//------------------------------------------------------------------------------------------------------------------------------------------

class EventSlicingThreeDTool;

/**
 *  @brief  SlicingThreeDAlgorithm class
 */
class SlicingThreeDAlgorithm : public SlicingAlgorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    SlicingThreeDAlgorithm();

private:
    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    EventSlicingThreeDBaseTool *m_pEventSlicingTool; ///< The address of the event slicing tool
    std::string m_slicingListDeletionAlgorithm;      ///< The name of the slicing list deletion algorithm

    HitTypeToNameMap m_caloHitListNames; ///< The hit type to calo hit list name map
    HitTypeToNameMap m_clusterListNames; ///< The hit type to cluster list name map

    std::string m_sliceClusterListName; ///< The name of the output slice cluster list
    std::string m_slicePfoListName;     ///< The name of the output slice pfo list
};

//------------------------------------------------------------------------------------------------------------------------------------------

} // namespace lar_content

#endif // #ifndef LAR_SLICING_THREE_D_ALGORITHM_H
