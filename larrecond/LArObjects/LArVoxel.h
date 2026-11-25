/**
 *  @file   larrecond/LArObjects/LArVoxel.h
 *
 *  @brief  Header file for LArVoxel and LArVoxelProjection
 *
 *  $Log: $
 */
#ifndef PANDORA_LAR_VOXEL_H
#define PANDORA_LAR_VOXEL_H 1

#include <vector>

#include "Pandora/PandoraInputTypes.h"

namespace lar_nd_reco
{

class LArVoxel
{
public:
    /**
     *  @brief  Constructor
     *
     *  @param  voxelID Total bin number for the voxel (long integer, since it can be > 2^31)
     *  @param  energyInVoxel The total deposited energy in the voxel (GeV)
     *  @param  voxelPosVect Voxel position, set as the first corner of the voxel bin
     *  @param  trackID The Geant4 ID of the (first) contributing track to this voxel
     */
    LArVoxel(const long voxelID, const float energyInVoxel, const pandora::CartesianVector &voxelPosVect, const int trackID);

    /**
     *  @brief  Constructor
     *
     *  @param  voxelID Total bin number for the voxel (long integer, since it can be > 2^31)
     *  @param  energyInVoxel The total deposited energy in the voxel (GeV)
     *  @param  voxelPosVect Voxel position, set as the first corner of the voxel bin
     *  @param  trackID The Geant4 ID of the (first) contributing track to this voxel
     *  @param  tpcID ID of the tpc containing this voxel
     */
    LArVoxel(const long voxelID, const float energyInVoxel, const pandora::CartesianVector &voxelPosVect, const int trackID, const unsigned int tpcID);

    /**
     *  @brief  Set voxel energy (GeV)
     *
     *  @param  E voxel energy
     */
    void SetEnergy(const float E);

    /**
     *  @brief  Set voxel track id
     *
     *  @param  trackid is the track id
     */
    void SetTrackID(const int trackid);

    long m_voxelID;                          ///< The long integer ID of the voxel (can be larger than 2^31)
    float m_energyInVoxel;                   ///< The energy in the voxel (GeV)
    pandora::CartesianVector m_voxelPosVect; ///< Position vector (x,y,z) of the first voxel corner
    int m_trackID;                           ///< The Geant4 ID of the (first) contributing track to this voxel
    int m_tpcID;
};

typedef std::vector<LArVoxel> LArVoxelList;

//------------------------------------------------------------------------------------------------------------------------------------------

inline LArVoxel::LArVoxel(const long voxelID, const float energyInVoxel, const pandora::CartesianVector &voxelPosVect, const int trackID) :
    m_voxelID(voxelID),
    m_energyInVoxel(energyInVoxel),
    m_voxelPosVect(voxelPosVect),
    m_trackID(trackID),
    m_tpcID(0)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline LArVoxel::LArVoxel(const long voxelID, const float energyInVoxel, const pandora::CartesianVector &voxelPosVect, const int trackID,
    const unsigned int tpcID) :
    m_voxelID(voxelID),
    m_energyInVoxel(energyInVoxel),
    m_voxelPosVect(voxelPosVect),
    m_trackID(trackID),
    m_tpcID(tpcID)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void LArVoxel::SetEnergy(const float E)
{
    m_energyInVoxel = E;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void LArVoxel::SetTrackID(const int trackid)
{
    m_trackID = trackid;
}

//------------------------------------------------------------------------------------------------------------------------------------------

class LArVoxelProjection
{
public:
    /**
     *  @brief  Constructor
     *
     *  @param  energy the energy deposited in the voxel before projection
     *  @param  w the coordinate in the wire direction after projection
     *  @param  x the drift coordinate
     *  @param  view the readout view - typically U, V, or W
     *  @param  parentid the id of the parent voxel
     *  @param  trackid id of the true particle that produced the energy in the voxel before projection
     */
    LArVoxelProjection(const float energy, const float w, const float x, const pandora::HitType &view, const int parentid, const int trackid);

    /**
     *  @brief  Constructor
     *
     *  @param  energy the energy deposited in the voxel before projection
     *  @param  w the coordinate in the wire direction after projection
     *  @param  x the drift coordinate
     *  @param  view the readout view - typically U, V, or W
     *  @param  parentid the id of the parent voxel
     *  @param  trackid id of the true particle that produced the energy in the voxel before projection
     *  @param  tpcid id of the TPC containing the voxel
     */
    LArVoxelProjection(const float energy, const float w, const float x, const pandora::HitType &view, const int parentid,
        const int trackid, const unsigned int tpcid);

    float m_energy;          ///< energy deposited in the voxel
    float m_wire;            ///< projected wire coordinate of the voxel
    float m_drift;           ///< drift coordinate (x coordinate of the voxel)
    pandora::HitType m_view; ///< hit type to label the U, V and W views
    int m_parentVoxelID;     ///< id of the parent voxel
    int m_trackID;           ///< true particle responsible for the majority of the energy
    int m_tpcID;             ///< id of the TPC containing the voxel
};

typedef std::vector<LArVoxelProjection> LArVoxelProjectionList;

//------------------------------------------------------------------------------------------------------------------------------------------

inline LArVoxelProjection::LArVoxelProjection(
    const float energy, const float w, const float x, const pandora::HitType &view, const int parentid, const int trackid) :
    m_energy(energy),
    m_wire(w),
    m_drift(x),
    m_view(view),
    m_parentVoxelID(parentid),
    m_trackID(trackid),
    m_tpcID(0)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline LArVoxelProjection::LArVoxelProjection(const float energy, const float w, const float x, const pandora::HitType &view,
    const int parentid, const int trackid, const unsigned int tpcid) :
    m_energy(energy),
    m_wire(w),
    m_drift(x),
    m_view(view),
    m_parentVoxelID(parentid),
    m_trackID(trackid),
    m_tpcID(tpcid)
{
}

} // namespace lar_nd_reco

#endif
