/** 
 * @file   include/RockMuonTaggingTool.h
 *
 *  @brief  Header file for the rock muon tagging tool class.
 *
 *  $Log: $
 */

#ifndef LAR_ROCK_MUON_TAGGING_TOOL_H
#define LAR_ROCK_MUON_TAGGING_TOOL_H 1

#include "larpandoracontent/LArControlFlow/CosmicRayTaggingTool.h"
#include "larpandoracontent/LArControlFlow/MasterAlgorithm.h"

#include <unordered_map>

namespace lar_content
{

/**
 *  @brief  RockMuonTaggingAlgorithm class
 */

class RockMuonTaggingTool : public CosmicRayTaggingTool
{
public:
    /**
     *  @brief  Default constructor
     */
    RockMuonTaggingTool();

    void FindAmbiguousPfos(const pandora::PfoList &parentCosmicRayPfos, pandora::PfoList &ambiguousPfos, const MasterAlgorithm *const pAlgorithm) override;

    /**
     *  @brief Tag clear rock muons and exclude it from ambiguousPfos 
     *
     *  @param  ambiguousPfos input pfos list
     */
    pandora::StatusCode TagRockMuonPfos(pandora::PfoList &ambiguousPfos) const;

private:
    /**
     *  @brief  Find which detector face a 3D point is close to, within the configured margins
     *
     *  @param  x point x coordinate
     *  @param  y point y coordinate
     *  @param  z point z coordinate
     *
     *  @return 0 = left (Xa), 1 = right (Xc), 2 = bottom (Yb), 3 = top (Yt), 4 = upstream (Zu), 5 = downstream (Zd),
     *          or -1 if the point is not close to any face. If the point is close to more than one face (edge or corner),
     *          the first match in the order X, Y, Z is returned.
     */
    int GetNearFace(const float x, const float y, const float z) const;

    /**
     *  @brief  Check if each candidate is throughgoing (i.e emerging and exiting from any of the detector boundaies)
     *
     *  @param  candidate input
     */
    bool CheckIfThroughgoing(const CRCandidate &candidate) const;

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle) override;

    float m_marginX; ///< the minimum distance from the dector X-face to define a fiducial volume for tagging
    float m_marginY; ///< the minimum distance from the dector Y-face to define a fiducial volume for tagging
    float m_marginZ; ///< the minimum distance from the dector Z-face to define a fiducial volume for tagging
};

} // namespace lar_content
#endif // #ifndef LAR_ROCK_MUON_TAGGING_TOOL_H
