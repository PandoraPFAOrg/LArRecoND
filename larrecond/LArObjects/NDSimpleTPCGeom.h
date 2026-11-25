/**
 *  @file   larrecond/LArObjects/NDSimpleTPCGeom.h
 *
 *  @brief  Header file for storing TPC-like geometry for the ND
 *
 *  $Log: $
 */
#ifndef PANDORA_ND_SIMPLE_TPC_GEOM_H
#define PANDORA_ND_SIMPLE_TPC_GEOM_H 1

#include "Pandora/PandoraInputTypes.h"

namespace lar_nd_reco
{

class NDSimpleTPC
{
public:
    /**
     *  @brief  default constructor
     */
    NDSimpleTPC();

    /**
     *  @brief  constructor with coordinate limits and id
     *
     *  @param  x_min minimum x value of the TPC box
     *  @param  x_max maximum x value of the TPC box
     *  @param  y_min minimum y value of the TPC box
     *  @param  y_max maximum y value of the TPC box
     *  @param  z_min minimum z value of the TPC box
     *  @param  z_max maximum z value of the TPC box
     *  @param  tpcID unique id of the tpc
     */
    NDSimpleTPC(const double x_min, const double x_max, const double y_min, const double y_max, const double z_min, const double z_max, const int tpcID);

    /**
     *  @brief  Find out if a 3D point is inside a given TPC
     *
     *  @param  pos the 3D position to query
     *  @param  ignoreY flag to ignore the y coordinate (for 2D points)
     *  @return whether the point is inside or not
     */
    bool IsInTPC(const pandora::CartesianVector &pos, bool ignoreY = false) const;

    double m_x_min; ///< minimum x value of the TPC cubioid
    double m_x_max; ///< maximum x value of the TPC cubioid
    double m_y_min; ///< minimum y value of the TPC cubioid
    double m_y_max; ///< maximum y value of the TPC cubioid
    double m_z_min; ///< minimum z value of the TPC cubioid
    double m_z_max; ///< maximum z value of the TPC cubioid
    int m_TPC_ID;   ///< unique id of the TPC
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline NDSimpleTPC::NDSimpleTPC() :
    m_x_min{0.},
    m_x_max{0.},
    m_y_min{0.},
    m_y_max{0.},
    m_z_min{0.},
    m_z_max{0.},
    m_TPC_ID{-1}
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline NDSimpleTPC::NDSimpleTPC(const double x_min, const double x_max, const double y_min, const double y_max, const double z_min,
    const double z_max, const int tpcID) :
    m_x_min{x_min},
    m_x_max{x_max},
    m_y_min{y_min},
    m_y_max{y_max},
    m_z_min{z_min},
    m_z_max{z_max},
    m_TPC_ID{tpcID}
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline bool NDSimpleTPC::IsInTPC(const pandora::CartesianVector &pos, bool ignoreY) const
{
    const double epsilon{1.0e-3};
    const double m_x_min_eps{m_x_min - epsilon};
    const double m_x_max_eps{m_x_max + epsilon};
    const double m_y_min_eps{m_y_min - epsilon};
    const double m_y_max_eps{m_y_max + epsilon};
    const double m_z_min_eps{m_z_min - epsilon};
    const double m_z_max_eps{m_z_max + epsilon};
    const double x{pos.GetX()};
    const double y{pos.GetY()};
    const double z{pos.GetZ()};

    if (x <= m_x_min_eps || x >= m_x_max_eps)
        return false;
    if (!ignoreY && (y <= m_y_min_eps || y >= m_y_max_eps))
        return false;
    if (z <= m_z_min_eps || z >= m_z_max_eps)
        return false;

    return true;
}

class NDSimpleTPCGeom
{
public:
    /**
     *  @brief  default constructor
     */
    NDSimpleTPCGeom();

    /**
     *  @brief  Get the TPC number from a 3D position
     *
     *  @param  position 3d position to query
     *  @param  ignoreY boolean flag to ignore the y coordinate (for 2D points)
     *  @return tpc number as unsigned int
     */
    int GetTPCNumber(const pandora::CartesianVector &position, bool ignoreY = false) const;

    /**
     *  @brief  Get the module number from a 3D position
     *
     *  @param  position 3d position to query
     *  @param  ignoreY boolean flag to ignore the y coordinate (for 2D points)
     *  @return module number as unsigned int
     */
    int GetModuleNumber(const pandora::CartesianVector &position, bool ignoreY = false) const;

    /**
     *  @brief  Add a TPC to the geometry
     *
     *  @param  min_x minimum x position of the TPC
     *  @param  max_x maximum x position of the TPC
     *  @param  min_y minimum y position of the TPC
     *  @param  max_y maximum y position of the TPC
     *  @param  min_z minimum z position of the TPC
     *  @param  max_z maximum z position of the TPC
     *  @param  tpcID tpc id as int
     */
    void AddTPC(const double min_x, const double max_x, const double min_y, const double max_y, const double min_z, const double max_z, const int tpcID);

    /**
     *  @brief  Get the box surrounding all TPCs
     *
     *  @param  min_x minimum x position of the TPCs
     *  @param  max_x maximum x position of the TPCs
     *  @param  min_y minimum y position of the TPCs
     *  @param  max_y maximum y position of the TPCs
     *  @param  min_z minimum z position of the TPCs
     *  @param  max_z maximum z position of the TPCs
     */
    void GetSurroundingBox(double &min_x, double &max_x, double &min_y, double &max_y, double &min_z, double &max_z) const;

    std::map<unsigned int, NDSimpleTPC> m_TPCs; ///< map of the TPCs keyed on their unique id

private:
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline NDSimpleTPCGeom::NDSimpleTPCGeom()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline int NDSimpleTPCGeom::GetTPCNumber(const pandora::CartesianVector &position, bool ignoreY) const
{
    for (auto const &tpc : m_TPCs)
    {
        if (tpc.second.IsInTPC(position, ignoreY))
            return tpc.first;
    }
    return -1;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline int NDSimpleTPCGeom::GetModuleNumber(const pandora::CartesianVector &position, bool ignoreY) const
{
    for (auto const &tpc : m_TPCs)
    {
        if (tpc.second.IsInTPC(position, ignoreY))
            return tpc.first / 2;
    }
    return -1;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void NDSimpleTPCGeom::AddTPC(
    const double min_x, const double max_x, const double min_y, const double max_y, const double min_z, const double max_z, const int tpcID)
{
    if (m_TPCs.count(tpcID))
        std::cout << "NDSimpleTPCGeom: trying to add another TPC with tpc id " << tpcID << "! Doing nothing. " << std::endl;
    else
        m_TPCs[tpcID] = NDSimpleTPC(min_x, max_x, min_y, max_y, min_z, max_z, tpcID);
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void NDSimpleTPCGeom::GetSurroundingBox(double &min_x, double &max_x, double &min_y, double &max_y, double &min_z, double &max_z) const
{
    min_x = std::numeric_limits<double>::max();
    min_y = std::numeric_limits<double>::max();
    min_z = std::numeric_limits<double>::max();
    max_x = -std::numeric_limits<double>::max();
    max_y = -std::numeric_limits<double>::max();
    max_z = -std::numeric_limits<double>::max();

    for (auto const &tpc : m_TPCs)
    {
        min_x = tpc.second.m_x_min < min_x ? tpc.second.m_x_min : min_x;
        min_y = tpc.second.m_y_min < min_y ? tpc.second.m_y_min : min_y;
        min_z = tpc.second.m_z_min < min_z ? tpc.second.m_z_min : min_z;
        max_x = tpc.second.m_x_max > max_x ? tpc.second.m_x_max : max_x;
        max_y = tpc.second.m_y_max > max_y ? tpc.second.m_y_max : max_y;
        max_z = tpc.second.m_z_max > max_z ? tpc.second.m_z_max : max_z;
    }
}

} // namespace lar_nd_reco

#endif
