/**
 *  @file   larrecond/LArObjects/LArGrid.h
 *
 *  @brief  Header file for the voxelisation grid
 *
 *  $Log: $
 */
#ifndef PANDORA_LAR_GRID_H
#define PANDORA_LAR_GRID_H 1

#include <array>

#include "Pandora/PandoraInputTypes.h"

#include "larrecond/LArObjects/LArBox.h"

namespace lar_nd_reco
{

typedef std::array<long, 3> LongBin3Array;
typedef std::array<long, 4> LongBin4Array;

class LArGrid : public LArBox
{
public:
    /**
     *  @brief  Constructor
     *
     *  @param  bottom The bottom starting corner of the grid (CartesianVector)
     *  @param  top The top end corner of the grid (CartesianVector)
     *  @param  binWidths The regular bin widths stored as a CartesianVector(dx, dy, dz)
     */
    LArGrid(const pandora::CartesianVector &bottom, const pandora::CartesianVector &top, const pandora::CartesianVector &binWidths);

    /**
     *  @brief  Get the (x,y,z,total) bin indices for the given point. Need long integers, since total bin can be > 2^31
     *
     *  @param  point The CartesianVector position
     *
     *  @return The array of (x,y,z,total) bin long integer indices
     */
    LongBin4Array GetBinIndices(const pandora::CartesianVector &point) const;

    /**
     *  @brief  Get the position for the given x, y and z bins
     *
     *  @param  xBin The x bin index (long integer)
     *  @param  yBin The y bin index (long integer)
     *  @param  zBin The z bin index (long integer)
     *
     *  @return The CartesianVector position
     */
    pandora::CartesianVector GetPoint(const long xBin, const long yBin, const long zBin) const;

    /**
     *  @brief  Get the position for the given array of bin indices (x,y,z,total)
     *
     *  @param  bins The (x,y,z,total) bins as an array of long integers (since total can be > 2^31)
     *
     *  @return The CartesianVector position
     */
    pandora::CartesianVector GetPoint(const LongBin4Array &bins) const;

    pandora::CartesianVector m_bottom;    ///< The bottom corner of the box
    pandora::CartesianVector m_top;       ///< The top corner of the box
    pandora::CartesianVector m_binWidths; ///< The bin widths (dx, dy, dz)
    LongBin3Array m_nBins;                ///< The (x,y,z) bin indices (long integers)
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline LArGrid::LArGrid(const pandora::CartesianVector &bottom, const pandora::CartesianVector &top, const pandora::CartesianVector &binWidths) :
    LArBox(bottom, top),
    m_bottom(bottom),
    m_top(top),
    m_binWidths(binWidths),
    m_nBins({0, 0, 0})
{
    if (binWidths.GetX() <= 0 || binWidths.GetY() <= 0 || binWidths.GetZ() <= 0)
        throw pandora::StatusCodeException(pandora::STATUS_CODE_OUT_OF_RANGE);

    const long NxBins = static_cast<long>((top.GetX() - bottom.GetX()) / binWidths.GetX());
    const long NyBins = static_cast<long>((top.GetY() - bottom.GetY()) / binWidths.GetY());
    const long NzBins = static_cast<long>((top.GetZ() - bottom.GetZ()) / binWidths.GetZ());
    m_nBins = {NxBins, NyBins, NzBins};
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline long GetBinIndex(const double x, const double botX, const double binWidth, const long maxBin)
{
    long index(0);
    if (binWidth > 0.0)
        index = static_cast<long>(((x - botX) / binWidth));

    if (index < 0)
        index = 0;
    else if (index >= maxBin)
        index = maxBin - 1;

    return index;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline LongBin4Array LArGrid::GetBinIndices(const pandora::CartesianVector &point) const
{
    // Bin widths should always be non-zero. Need to use long integers to avoid
    // integer overflow for the total bin, which can be larger than 2^31
    const long xBin = GetBinIndex(point.GetX(), m_bottom.GetX(), m_binWidths.GetX(), m_nBins[0]);
    const long yBin = GetBinIndex(point.GetY(), m_bottom.GetY(), m_binWidths.GetY(), m_nBins[1]);
    const long zBin = GetBinIndex(point.GetZ(), m_bottom.GetZ(), m_binWidths.GetZ(), m_nBins[2]);

    const long totBin = (zBin * m_nBins[1] + yBin) * m_nBins[0] + xBin;

    LongBin4Array binIndices = {xBin, yBin, zBin, totBin};
    return binIndices;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::CartesianVector LArGrid::GetPoint(const long xBin, const long yBin, const long zBin) const
{
    const float x = m_bottom.GetX() + xBin * m_binWidths.GetX();
    const float y = m_bottom.GetY() + yBin * m_binWidths.GetY();
    const float z = m_bottom.GetZ() + zBin * m_binWidths.GetZ();

    return pandora::CartesianVector(x, y, z);
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::CartesianVector LArGrid::GetPoint(const LongBin4Array &bins) const
{
    return GetPoint(bins[0], bins[1], bins[2]);
}

} // namespace lar_nd_reco

#endif
