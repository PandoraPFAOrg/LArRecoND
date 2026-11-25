/**
 *  @file   larrecond/LArObjects/LArBox.h
 *
 *  @brief  Header file for LArBox, an element in the voxelisation grid
 *
 *  $Log: $
 */
#ifndef PANDORA_LAR_BOX_H
#define PANDORA_LAR_BOX_H 1

#include "Pandora/PandoraInputTypes.h"

#include "larrecond/LArObjects/LArRay.h"

namespace lar_nd_reco
{

class LArBox
{
public:
    /**
     *  @brief  Constructor
     *
     *  @param  bottom The bottom box corner
     *  @param  top The top box corner
     */
    LArBox(const pandora::CartesianVector &bottom, const pandora::CartesianVector &top);

    /**
     *  @brief  Find box intersection lengths for given ray
     *
     *  @param  ray The ray (starting point and direction)
     *  @param  t0 The first intersection length along ray from its starting point (double precision)
     *  @param  t1 The second intersection length along ray from its starting point (double precision)
     *
     *  @return Success/failure of finding both intersection lengths
     */
    bool Intersect(const LArRay &ray, double &t0, double &t1) const;

    /**
     *  @brief  Check if the given point is inside the box
     *
     *  @param  point The point (CartesianVector)
     *
     *  @return Success/failure
     */
    bool Inside(const pandora::CartesianVector &point) const;

    pandora::CartesianVector m_bottom; ///< The bottom corner of the box
    pandora::CartesianVector m_top;    ///< The top corner of the box
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline LArBox::LArBox(const pandora::CartesianVector &bottom, const pandora::CartesianVector &top) :
    m_bottom(bottom),
    m_top(top)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArBox::Intersect(const LArRay &ray, double &t0, double &t1) const
{
    // Brian Smits ray-box intersection algorithm with improvements from Amy Williams et al. Code based on
    // https://github.com/chenel/larcv2/tree/edepsim-formattruth/larcv/app/Supera/Voxelize.cxx (MIT license)
    // which uses code from i) "An Efficient and Robust Ray–Box Intersection Algorithm", Amy Williams et al. (2004):
    // https://doi.org/10.1145/1198555.1198748 and ii) "Efficiency Issues for Ray Tracing", Brian Smits (1998):
    // https://doi.org/10.1080/10867651.1998.10487488. The Smits-Williams GPLv3 licensed code is available from:
    // https://www.scratchapixel.com/lessons/3d-basic-rendering/minimal-ray-tracer-rendering-simple-shapes/ray-box-intersection
    const double x(ray.m_origin.GetX());
    const double y(ray.m_origin.GetY());

    const double invDirX(ray.m_invDir.GetX());
    const double invDirY(ray.m_invDir.GetY());

    double tMin(0.0), tMax(0.0), tyMin(0.0), tyMax(0.0), tzMin(0.0), tzMax(0.0);

    if (ray.m_sign[0] == 0)
    {
        tMin = (m_bottom.GetX() - x) * invDirX;
        tMax = (m_top.GetX() - x) * invDirX;
    }
    else
    {
        tMin = (m_top.GetX() - x) * invDirX;
        tMax = (m_bottom.GetX() - x) * invDirX;
    }

    if (ray.m_sign[1] == 0)
    {
        tyMin = (m_bottom.GetY() - y) * invDirY;
        tyMax = (m_top.GetY() - y) * invDirY;
    }
    else
    {
        tyMin = (m_top.GetY() - y) * invDirY;
        tyMax = (m_bottom.GetY() - y) * invDirY;
    }

    if ((tMin > tyMax) || (tyMin > tMax))
        return false;

    if (tyMin > tMin)
        tMin = tyMin;

    if (tyMax < tMax)
        tMax = tyMax;

    const double z(ray.m_origin.GetZ());
    const double invDirZ(ray.m_invDir.GetZ());
    if (ray.m_sign[2] == 0)
    {
        tzMin = (m_bottom.GetZ() - z) * invDirZ;
        tzMax = (m_top.GetZ() - z) * invDirZ;
    }
    else
    {
        tzMin = (m_top.GetZ() - z) * invDirZ;
        tzMax = (m_bottom.GetZ() - z) * invDirZ;
    }

    if ((tMin > tzMax) || (tzMin > tMax))
        return false;

    if (tzMin > tMin)
        tMin = tzMin;

    if (tzMax < tMax)
        tMax = tzMax;

    // First and second intersection path lengths
    t0 = tMin;
    t1 = tMax;
    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArBox::Inside(const pandora::CartesianVector &point) const
{
    const float x = point.GetX();
    const float y = point.GetY();
    const float z = point.GetZ();

    if (x > m_bottom.GetX() && x < m_top.GetX() && y > m_bottom.GetY() && y < m_top.GetY() && z > m_bottom.GetZ() && z < m_top.GetZ())
    {
        return true;
    }

    return false;
}

} // namespace lar_nd_reco

#endif
