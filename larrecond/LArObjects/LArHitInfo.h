/**
 *  @file   larrecond/LArObjects/LArHitInfo.h
 *
 *  @brief  Header file for storing hit information
 *
 *  $Log: $
 */
#ifndef PANDORA_LAR_HIT_INFO_H
#define PANDORA_LAR_HIT_INFO_H 1

#include "Pandora/PandoraInputTypes.h"

namespace lar_nd_reco
{

class LArHitInfo
{
public:
    /**
     *  @brief  Constructor
     *
     *  @param  start The starting point of the hit step
     *  @param  stop The end point of the hit step
     *  @param  energy The total energy of the hit step
     *  @param  trackID The ID of the main contributing particle
     *  @param  lengthScale Scaling factor to use cm length dimensions
     *  @param  energyScale Scaling factor to use GeV energies
     */
    LArHitInfo(const pandora::CartesianVector &start, const pandora::CartesianVector &stop, const float energy, const int trackID,
        const float lengthScale, const float energyScale);

    pandora::CartesianVector m_start; ///< Starting point of the hit step
    pandora::CartesianVector m_stop;  ///< End point of the hit step
    float m_energy;                   ///< The total energy of the step
    int m_trackID;                    ///< The ID of the (main) contributing particle
};

inline LArHitInfo::LArHitInfo(const pandora::CartesianVector &start, const pandora::CartesianVector &stop, const float energy,
    const int trackID, const float lengthScale, const float energyScale) :
    m_start(start * lengthScale),
    m_stop(stop * lengthScale),
    m_energy(energy * energyScale),
    m_trackID(trackID)
{
}
} // namespace lar_nd_reco

#endif
