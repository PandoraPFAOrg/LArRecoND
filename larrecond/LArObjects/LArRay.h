/**
 *  @file   larrecond/LArObjects/LArRay.h
 *
 *  @brief  Header file for the voxel tracing rays LArRay
 *
 *  $Log: $
 */
#ifndef PANDORA_LAR_RAY_H
#define PANDORA_LAR_RAY_H 1

#include <array>

#include "Pandora/PandoraInputTypes.h"

namespace lar_nd_reco
{

typedef std::array<int, 3> SignComponentArray;

class LArRay
{
public:
    /**
     *  @brief  Constructor
     *
     *  @param  origin The starting point of the ray
     *  @param  dir The direction vector of the ray
     */
    LArRay(const pandora::CartesianVector &origin, const pandora::CartesianVector &dir);

    /**
     *  @brief  Update ray starting position
     *
     *  @param  newOrigin The new starting position (CartesianVector)
     */
    void UpdateOrigin(const pandora::CartesianVector &newOrigin);

    /**
     *  @brief  Get position along ray from its start point
     *
     *  @param  length Distance (cm) along the ray from its starting point
     *
     *  @return The position as a CartesianVector
     */
    pandora::CartesianVector GetPoint(const double length) const;

    pandora::CartesianVector m_origin; ///< Starting point of the ray
    pandora::CartesianVector m_dir;    ///< Ray direction
    pandora::CartesianVector m_invDir; ///< Reciprocal direction vector
    SignComponentArray m_sign;         ///< Sign of reciprocal direction vector components
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline LArRay::LArRay(const pandora::CartesianVector &origin, const pandora::CartesianVector &dir) :
    m_origin(origin),
    m_dir(dir),
    m_invDir({0.f, 0.f, 0.f}),
    m_sign({0, 0, 0})
{
    const float dx(m_dir.GetX());
    const float dy(m_dir.GetY());
    const float dz(m_dir.GetZ());

    // maxVal is needed to set the inverse direction components for parallel lines
    const float maxVal(std::numeric_limits<float>::max());
    const float invdx(std::fabs(dx) > 0.0 ? 1.0 / dx : maxVal);
    const float invdy(std::fabs(dy) > 0.0 ? 1.0 / dy : maxVal);
    const float invdz(std::fabs(dz) > 0.0 ? 1.0 / dz : maxVal);

    m_invDir.SetValues(invdx, invdy, invdz);

    m_sign[0] = (m_invDir.GetX() < 0);
    m_sign[1] = (m_invDir.GetY() < 0);
    m_sign[2] = (m_invDir.GetZ() < 0);
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void LArRay::UpdateOrigin(const pandora::CartesianVector &newOrigin)
{
    m_origin = newOrigin;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::CartesianVector LArRay::GetPoint(const double length) const
{
    return (m_origin + m_dir * length);
}

//------------------------------------------------------------------------------------------------------------------------------------------

} // namespace lar_nd_reco

#endif
