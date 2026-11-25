/**
 *  @file   larrecond/LArObjects/LArNDGeomSimple.h
 *
 *  @brief  Header file for storing geometry information for ND LAr
 *
 *  $Log: $
 */
#ifndef PANDORA_LAR_ND_GEOM_SIMPLE_H
#define PANDORA_LAR_ND_GEOM_SIMPLE_H 1

#include "Pandora/PandoraInputTypes.h"

namespace lar_nd_reco
{

class LArNDTPCSimple
{
public:
    /**
     *  @brief  default constructor
     */
    LArNDTPCSimple();

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
    LArNDTPCSimple(const double x_min, const double x_max, const double y_min, const double y_max, const double z_min, const double z_max,
        const int tpcID);

    bool IsInTPC(const pandora::CartesianVector &pos) const;

    double m_x_min; ///< minimum x value of the TPC cubioid
    double m_x_max; ///< maximum x value of the TPC cubioid
    double m_y_min; ///< minimum y value of the TPC cubioid
    double m_y_max; ///< maximum y value of the TPC cubioid
    double m_z_min; ///< minimum z value of the TPC cubioid
    double m_z_max; ///< maximum z value of the TPC cubioid
    int m_TPC_ID;   ///< unique id of the TPC
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline LArNDTPCSimple::LArNDTPCSimple() :
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

inline LArNDTPCSimple::LArNDTPCSimple(const double x_min, const double x_max, const double y_min, const double y_max, const double z_min,
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

inline bool LArNDTPCSimple::IsInTPC(const pandora::CartesianVector &pos) const
{
    const double epsilon{1.0};
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
    if (y <= m_y_min_eps || y >= m_y_max_eps)
        return false;
    if (z <= m_z_min_eps || z >= m_z_max_eps)
        return false;

    return true;
}

class LArNDGeomSimple
{
public:
    /**
     *  @brief  default constructor
     */
    LArNDGeomSimple();

    /**
     *  @brief  Get the TPC number from a 3D position
     *
     *  @param  position 3d position to query
     *
     *  @return tpc number as unsigned int
     */
    int GetTPCNumber(const pandora::CartesianVector &position) const;

    /**
     *  @brief  Get the module number from a 3D position
     *
     *  @param  position 3d position to query
     *
     *  @return module number as unsigned int
     */
    int GetModuleNumber(const pandora::CartesianVector &position) const;

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

    std::map<unsigned int, LArNDTPCSimple> m_TPCs; ///< map of the TPCs keyed on their unique id

private:
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline LArNDGeomSimple::LArNDGeomSimple()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline int LArNDGeomSimple::GetTPCNumber(const pandora::CartesianVector &position) const
{
    for (auto const &tpc : m_TPCs)
    {
        if (tpc.second.IsInTPC(position))
            return tpc.first;
    }
    return -1;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline int LArNDGeomSimple::GetModuleNumber(const pandora::CartesianVector &position) const
{
    for (auto const &tpc : m_TPCs)
    {
        if (tpc.second.IsInTPC(position))
            return tpc.first / 2;
    }
    return -1;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void LArNDGeomSimple::AddTPC(
    const double min_x, const double max_x, const double min_y, const double max_y, const double min_z, const double max_z, const int tpcID)
{
    if (m_TPCs.count(tpcID))
        std::cout << "LArNDGeomSimple: trying to add another TPC with tpc id " << tpcID << "! Doing nothing. " << std::endl;
    else
        m_TPCs[tpcID] = LArNDTPCSimple(min_x, max_x, min_y, max_y, min_z, max_z, tpcID);
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void LArNDGeomSimple::GetSurroundingBox(double &min_x, double &max_x, double &min_y, double &max_y, double &min_z, double &max_z) const
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
