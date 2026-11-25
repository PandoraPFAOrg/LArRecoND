/**
 *  @file   larrecond/LArHelpers/NDGeometryHelper.h
 *
 *  @brief  Header file for the ND geometry helper
 *
 *  $Log: $
 */
#ifndef PANDORA_ND_GEOMETRY_HELPER_H
#define PANDORA_ND_GEOMETRY_HELPER_H 1

#include "larrecond/LArObjects/NDSimpleTPCGeom.h"

#include <memory>
#include <string>
#include <vector>

class TFile;
class TGeoHMatrix;
class TGeoManager;
class TGeoNode;

namespace lar_nd_reco
{

class NDGeometryHelper
{
public:
    /**
     *  @brief  Constructor
     *
     *  @param  geomInputFile The ROOT file containing the geometry
     *  @param  geomManagerName The TGeoManager name
     *  @param  lengthScale The length scale to have the dimensions in cm (default = 1)
     */
    NDGeometryHelper(const std::string &geomInputFile, const std::string &geomManagerName, const float lengthScale = 1.0f);

    /**
     *  @brief  Destructor
     */
    virtual ~NDGeometryHelper();

    /**
     *  @brief  Get the TPC geometry
     *
     *  @param  tpcName The common names of the TPC volumes
     *  @param  tpcIDOffset The integer offset to make all the TPC IDs unique
     *  @return the NDSimpleTPCGeom object containing the TPC geometry
     */
    NDSimpleTPCGeom GetTPCGeometry(const std::string &tpcName, int tpcIDOffset) const;

private:
    /**
     *  @brief  Recursively search for volumes with the target name
     *
     *  @param  targetName the volume name that we want to find
     *  @param  nodePaths daughter indices to recreate the path to the target nodes
     *  @param  currentPath path to the current position in the geometry
     */
    void RecursiveGeometrySearch(
        const std::string &targetName, std::vector<std::vector<unsigned int>> &nodePaths, std::vector<unsigned int> &currentPath) const;

    /**
     *  @brief  Create and register a tpc in pandora
     *
     *  @param  geom Simple representation of the geometry for assigning TPC numbers
     *  @param  pVolMatrix matrix required to convert TPC coordinates to world
     *  @param  pTargetNode pointer to the TPC geometry node
     *  @param  tpcNumber the number for the TPC volume
     */
    void AddTPC(NDSimpleTPCGeom &geom, const std::unique_ptr<TGeoHMatrix> &pVolMatrix, const TGeoNode *pTargetNode, const unsigned int tpcNumber) const;

    TFile *m_pGeomFile;          ///< The geometry ROOT file
    TGeoManager *m_pGeomManager; ///< The TGeoManager name
    float m_lengthScale;         ///< The length scale for cm dimensions
};

} // namespace lar_nd_reco

#endif
