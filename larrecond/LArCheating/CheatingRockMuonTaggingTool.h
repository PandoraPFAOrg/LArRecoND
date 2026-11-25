/**
 *  @file   larrecond/LArCheating/CheatingRockMuonTaggingTool.h
 *
 *  @brief  Header file for the cheating rock muon tagging tool class.
 *
 *  $Log: $
 */
#ifndef LAR_CHEATING_ROCK_MUON_TAGGING_TOOL_H
#define LAR_CHEATING_ROCK_MUON_TAGGING_TOOL_H 1

#include "larpandoracontent/LArControlFlow/CosmicRayTaggingBaseTool.h"
#include "larpandoracontent/LArControlFlow/MasterAlgorithm.h"

namespace lar_content
{

/**
 *  @brief  CheatingRockMuonTaggingTool class
 */
class CheatingRockMuonTaggingTool : public CosmicRayTaggingBaseTool
{
public:
    /**
     *  @brief  Default constructor
     */
    CheatingRockMuonTaggingTool();

    /**
     *  @brief  Find and remove out of fiducial volume rock muon pfos from the input pfo list.
     *
     *  @param  parentRockMuonPfos the input list of rock muon pfos
     *  @param  ambiguousPfos to receive the list of ambiguous rock muon pfos
     *  @param  pAlgorithm address of the calling algorithm
     */
    void FindAmbiguousPfos(const pandora::PfoList &parentRockMuonPfos, pandora::PfoList &ambiguousPfos, const MasterAlgorithm *const pAlgorithm);

private:
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    float m_maxRockMuonFraction; ///< The maximum rock muon fraction for a pfo to be declared an ambiguous rock muon
};

} // namespace lar_content

#endif // #ifndef LAR_CHEATING_ROCK_MUON_TAGGING_TOOL_H
