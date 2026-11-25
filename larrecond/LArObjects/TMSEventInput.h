/**
 *  @file   larrecond/LArObjects/TMSEventInput.h
 *
 *  @brief  Header file for the TMS event processor
 *
 *  $Log: $
 */
#ifndef PANDORA_TMS_EVENT_INPUT_H
#define PANDORA_TMS_EVENT_INPUT_H 1

#include <memory>

#include "larrecond/LArObjects/NDEventInput.h"
#include "larrecond/LArObjects/TMS.h"

namespace lar_nd_reco
{

class TMSEventInput : public NDEventInput
{
public:
    /**
     *  @brief  Constructor
     *
     *  @param  pPandora the address of the Pandora instance
     *  @param  parameters the application parameters
     *  @param  tpcGeom the TPC geometry
     */
    TMSEventInput(const pandora::Pandora *pPandora, const NDParameters &parameters, const NDSimpleTPCGeom &tpcGeom);

    /**
     *  @brief  Destructor
     */
    virtual ~TMSEventInput();

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
    std::unique_ptr<TMS> m_tms;
};

} // namespace lar_nd_reco

#endif
