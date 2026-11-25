/**
 *  @file   larrecond/LArObjects/NDEventInput.h
 *
 *  @brief  Header file for the ND event input
 *
 *  $Log: $
 */
#ifndef PANDORA_ND_EVENT_INPUT_H
#define PANDORA_ND_EVENT_INPUT_H 1

#include "Pandora/Pandora.h"

#include "larpandoracontent/LArObjects/LArCaloHit.h"

#include "larrecond/LArObjects/NDParameters.h"
#include "larrecond/LArObjects/NDSimpleTPCGeom.h"

class TFile;
class TTree;

namespace lar_nd_reco
{

class NDEventInput
{
public:
    /**
     *  @brief  Default constructor
     */
    NDEventInput();

    /**
     *  @brief  Constructor
     *
     *  @param  pPandora the address of the Pandora instance
     *  @param  parameters the application parameters
     *  @param  tpcGeom the TPC geometry
     */
    NDEventInput(const pandora::Pandora *pPandora, const NDParameters &parameters, const NDSimpleTPCGeom &tpcGeom);

    /**
     *  @brief  Destructor
     */
    virtual ~NDEventInput();

    /**
     *  @brief  Create the MCParticles
     *
     *  @param  eventIndex the event number
     */
    virtual void CreateMCParticles(int eventIndex) const = 0;

    /**
     *  @brief  Create the CaloHits
     *
     *  @param  eventIndex the event number
     *  @param  startHitIndex the starting hit index
     *
     *  @return the last hit index
     */
    virtual int CreateCaloHits(int eventIndex, int startHitIndex = 0) const = 0;

    /**
     *  @brief  Create the Tracks
     *
     *  @param  eventIndex the event number
     *  @param  startTrackIndex the starting track index
     *
     *  @return the last track index
     */
    virtual int CreateTracks(int eventIndex, int startTrackIndex = 0) const = 0;

    /**
     *  @brief  Get the Pandora instance
     *
     *  @return the address of the Pandora instance
     */
    const pandora::Pandora *GetPandoraInstance() const
    {
        return m_pPandora;
    }

    /**
     *  @brief  Get the total number of event entries
     *
     *  @return the total number of events
     */
    int GetNEntries() const
    {
        return m_nEntries;
    }

protected:
    /**
     *  @brief  Initialization
     */
    void Initialize();

    /**
     *  @brief  Link between a CaloHit and the MCParticle that contributed most to it
     */
    struct MCLink
    {
        long m_trackID;     ///< The (file-based) ID of the MCParticle
        float m_energyFrac; ///< The fraction of the hit energy from the MCParticle
    };

    /**
     *  @brief  Make the default 3D CaloHit parameters for a voxel (the parent address is not set)
     *
     *  @param  voxelPos The voxel position (cm)
     *  @param  voxelE The voxel energy (GeV)
     *  @param  tpcID The TPC number of the voxel (negative values are mapped to TPC 0)
     *
     *  @return The CaloHit parameters
     */
    lar_content::LArCaloHitParameters MakeCaloHitParameters(const pandora::CartesianVector &voxelPos, const float voxelE, const int tpcID) const;

    /**
     *  @brief  Create the CaloHits for a voxel: a 3D hit (if enabled) and U, V and W view hits (if enabled).
     *          Each hit is created in both the main and this Pandora instance, with a unique parent address.
     *
     *  @param  caloHitParameters The default parameters from MakeCaloHitParameters
     *  @param  voxelPos The voxel position (cm)
     *  @param  hitCounter The running hit counter, which is incremented for each hit (or reserved 3D hit) made
     *  @param  pMCLink The MC link to set for each hit, or nullptr if there is no truth information
     */
    void CreateVoxelCaloHits(const lar_content::LArCaloHitParameters &caloHitParameters, const pandora::CartesianVector &voxelPos,
        int &hitCounter, const MCLink *pMCLink) const;

    /**
     *  @brief  Whether verbose output is enabled for this Pandora instance
     *
     *  @return True if the Pandora settings request algorithm info to be displayed
     */
    bool ShouldDisplayInfo() const;

    const pandora::Pandora *m_pPandora;
    const pandora::Pandora *m_pMainPandora;
    NDParameters m_parameters;
    NDSimpleTPCGeom m_tpcGeom;
    TFile *m_pFile;
    TTree *m_pTree;
    int m_nEntries;
    lar_content::LArCaloHitFactory m_LArCaloHitFactory;
};

} // namespace lar_nd_reco

#endif
