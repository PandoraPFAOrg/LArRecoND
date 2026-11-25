/**
 *  @file   larrecond/LArHelpers/NDGeometryHelper.cc
 *
 *  @brief  Implementation of the ND TPC Geometry helper
 *
 *  $Log: $
 */

#include "Api/PandoraApi.h"

#include "TFile.h"
#include "TGeoBBox.h"
#include "TGeoManager.h"
#include "TGeoMatrix.h"
#include "TGeoShape.h"
#include "TGeoVolume.h"

#include "larrecond/LArHelpers/NDGeometryHelper.h"

namespace lar_nd_reco
{

NDGeometryHelper::NDGeometryHelper(const std::string &geomFileName, const std::string &geomManagerName, const float lengthScale) :
    m_pGeomFile{nullptr},
    m_pGeomManager{nullptr},
    m_lengthScale{lengthScale}
{
    m_pGeomFile = TFile::Open(geomFileName.c_str(), "READ");
    if (m_pGeomFile)
    {
        m_pGeomManager = dynamic_cast<TGeoManager *>(m_pGeomFile->Get(geomManagerName.c_str()));
        if (!m_pGeomManager)
            std::cout << "Error in NDGeometryHelper: can't find TGeoManager " << geomManagerName << std::endl;
    }
    else
        std::cout << "Error in NDGeometryHelper: can't open file " << geomFileName << std::endl;
}

//------------------------------------------------------------------------------------------------------------------------------------------

NDGeometryHelper::~NDGeometryHelper()
{
    if (m_pGeomFile)
        m_pGeomFile->Close();
}

//------------------------------------------------------------------------------------------------------------------------------------------

NDSimpleTPCGeom NDGeometryHelper::GetTPCGeometry(const std::string &tpcName, int tpcIDOffset) const
{
    NDSimpleTPCGeom geom;

    if (!m_pGeomManager)
        return geom;

    // Go through the geometry and find the paths to the nodes we are interested in
    std::vector<std::vector<unsigned int>> nodePaths; // Store the daughter indices in the path to the node
    std::vector<unsigned int> currentPath;
    this->RecursiveGeometrySearch(tpcName, nodePaths, currentPath);
    std::cout << "Found " << nodePaths.size() << " matches for volumes containing the name " << tpcName << std::endl;

    // Navigate to each node and use them to build the pandora geometry
    for (unsigned int n = 0; n < nodePaths.size(); ++n)
    {
        const TGeoNode *pTopNode = m_pGeomManager->GetCurrentNode();
        // We have to multiply together matrices at each depth to convert local coordinates to the world volume
        std::unique_ptr<TGeoHMatrix> pVolMatrix = std::make_unique<TGeoHMatrix>(*pTopNode->GetMatrix());
        for (unsigned int d = 0; d < nodePaths.at(n).size(); ++d)
        {
            m_pGeomManager->CdDown(nodePaths.at(n).at(d));
            const TGeoNode *pNode = m_pGeomManager->GetCurrentNode();
            std::unique_ptr<TGeoHMatrix> pMatrix = std::make_unique<TGeoHMatrix>(*pNode->GetMatrix());
            pVolMatrix->Multiply(pMatrix.get());
        }
        const TGeoNode *pTargetNode = m_pGeomManager->GetCurrentNode();

        const int tpcNumber = n + tpcIDOffset;
        this->AddTPC(geom, pVolMatrix, pTargetNode, tpcNumber);

        for (const unsigned int &daughter : nodePaths.at(n))
        {
            (void)daughter;
            m_pGeomManager->CdUp();
        }
    }

    std::cout << "Created " << nodePaths.size() << " TPCs" << std::endl;
    return geom;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void NDGeometryHelper::RecursiveGeometrySearch(
    const std::string &targetName, std::vector<std::vector<unsigned int>> &nodePaths, std::vector<unsigned int> &currentPath) const
{
    const std::string nodeName{m_pGeomManager->GetCurrentNode()->GetName()};
    if (nodeName.find(targetName) != std::string::npos)
    {
        nodePaths.emplace_back(currentPath);
    }
    else
    {
        for (unsigned int i = 0; i < m_pGeomManager->GetCurrentNode()->GetNdaughters(); ++i)
        {
            m_pGeomManager->CdDown(i);
            currentPath.emplace_back(i);
            this->RecursiveGeometrySearch(targetName, nodePaths, currentPath);
            m_pGeomManager->CdUp();
            currentPath.pop_back();
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void NDGeometryHelper::AddTPC(
    NDSimpleTPCGeom &geom, const std::unique_ptr<TGeoHMatrix> &pVolMatrix, const TGeoNode *pTargetNode, const unsigned int tpcNumber) const
{
    // Get the BBox dimensions from the placement volume, which is assumed to be a cube
    TGeoVolume *pCurrentVol = pTargetNode->GetVolume();
    TGeoShape *pCurrentShape = (pCurrentVol != nullptr) ? pCurrentVol->GetShape() : nullptr;
    // pCurrentShape->InspectShape();
    TGeoBBox *pBox = dynamic_cast<TGeoBBox *>(pCurrentShape);

    // Now can get origin/width data from the BBox
    const double dx = (pBox != nullptr) ? pBox->GetDX() * m_lengthScale : 0.0; // Note these are the half widths
    const double dy = (pBox != nullptr) ? pBox->GetDY() * m_lengthScale : 0.0;
    const double dz = (pBox != nullptr) ? pBox->GetDZ() * m_lengthScale : 0.0;
    const double *pOrigin = (pBox != nullptr) ? pBox->GetOrigin() : nullptr;

    // Translate local origin to global coordinates
    double level1[3] = {0.0, 0.0, 0.0};
    if (pTargetNode && pOrigin)
        pTargetNode->LocalToMasterVect(pOrigin, level1);

    const double *pVolTrans = pVolMatrix->GetTranslation();
    const double centreX = (pVolTrans != nullptr) ? (level1[0] + pVolTrans[0]) * m_lengthScale : 0.0;
    const double centreY = (pVolTrans != nullptr) ? (level1[1] + pVolTrans[1]) * m_lengthScale : 0.0;
    const double centreZ = (pVolTrans != nullptr) ? (level1[2] + pVolTrans[2]) * m_lengthScale : 0.0;

    geom.AddTPC(centreX - dx, centreX + dx, centreY - dy, centreY + dy, centreZ - dz, centreZ + dz, tpcNumber);

    std::cout << "Creating TPC: " << centreX - dx << ", " << centreX + dx << ", " << centreY - dy << ", " << centreY + dy << ", "
              << centreZ - dz << ", " << centreZ + dz << std::endl;
}

} // namespace lar_nd_reco
