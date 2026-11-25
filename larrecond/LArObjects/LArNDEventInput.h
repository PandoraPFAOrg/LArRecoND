/**
 *  @file   larrecond/LArObjects/LArNDEventInput.h
 *
 *  @brief  Header file for the LArND event processor
 *
 *  $Log: $
 */
#ifndef PANDORA_LARND_EVENT_INPUT_H
#define PANDORA_LARND_EVENT_INPUT_H 1

#include <memory>

#include "larrecond/LArObjects/LArSP.h"
#include "larrecond/LArObjects/LArSPMC.h"
#include "larrecond/LArObjects/NDEventInput.h"

namespace lar_nd_reco
{

class LArNDEventInput : public NDEventInput
{
public:
    /**
     *  @brief  Constructor
     *
     *  @param  pPandora the address of the Pandora instance
     *  @param  parameters the application parameters
     *  @param  tpcGeom the TPC geometry
     */
    LArNDEventInput(const pandora::Pandora *pPandora, const NDParameters &parameters, const NDSimpleTPCGeom &tpcGeom);

    /**
     *  @brief  Destructor
     */
    virtual ~LArNDEventInput();

    /**
     *  @brief  Create the MCParticles
     *
     *  @param  eventIndex the event number
     */
    virtual void CreateMCParticles(int eventIndex) const override;

    /**
     *  @brief  Create the CaloHits
     *
     *  @param  eventIndex the event number
     *  @param  startHitIndex the starting hit index
     *
     *  @return the last hit index
     */
    virtual int CreateCaloHits(int eventIndex, int startHitIndex = 0) const override;

    /**
     *  @brief  Create the Tracks
     *
     *  @param  eventIndex the event number
     *  @param  startTrackIndex the starting track index
     *
     *  @return the last track index
     */
    virtual int CreateTracks(int eventIndex, int startTrackIndex = 0) const override;

private:
    std::unique_ptr<LArSP> m_larsp;
};

} // namespace lar_nd_reco

#endif
