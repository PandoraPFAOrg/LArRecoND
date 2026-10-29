/**
 *  @file   src/RockMuonTaggingTool.cc
 *
 *  @brief  Implementation of the rock muon tagging tool class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArControlFlow/CosmicRayTaggingTool.h"
#include <memory>

#include "RockMuonTaggingTool.h"

using namespace pandora;

namespace lar_content
{

RockMuonTaggingTool::RockMuonTaggingTool() :
    m_marginX(5.f), // [cm]
    m_marginY(5.f),
    m_marginZ(5.f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------
bool RockMuonTaggingTool::IsOutsideBox(const float x, const float y, const float z) const
{
    const float BoxXmin = m_face_Xa + m_marginX;
    const float BoxXmax = m_face_Xc - m_marginX;

    const float BoxYmin = m_face_Yb + m_marginY;
    const float BoxYmax = m_face_Yt - m_marginY;

    const float BoxZmin = m_face_Zu + m_marginZ;
    const float BoxZmax = m_face_Zd - m_marginZ;

    const bool IsOutRangeX = (x < BoxXmin || x > BoxXmax);
    const bool IsOutRangeY = (y < BoxYmin || y > BoxYmax);
    const bool IsOutRangeZ = (z < BoxZmin || z > BoxZmax);

    return (IsOutRangeZ || IsOutRangeY || IsOutRangeX);
}

//------------------------------------------------------------------------------------------------------------------------------------------
bool RockMuonTaggingTool::CheckIfThroughgoing(const CRCandidate &candidate) const
{
    return ((candidate.m_endPoint1.GetX() != std::numeric_limits<float>::max()) && // candidates whose sliding fit fails, have default values set as inf
        IsOutsideBox(candidate.m_endPoint1.GetX(), candidate.m_endPoint1.GetY(), candidate.m_endPoint1.GetZ()) &&
        IsOutsideBox(candidate.m_endPoint2.GetX(), candidate.m_endPoint2.GetY(), candidate.m_endPoint2.GetZ()));
}

//----------------------------------------------------- -------------------------------------------------------------------------------------
void RockMuonTaggingTool::FindAmbiguousPfos(const PfoList &parentCosmicRayPfos, PfoList &ambiguousPfos, const MasterAlgorithm *const /*pAlgorithm*/)
{
    if (this->GetPandora().GetSettings()->ShouldDisplayAlgorithmInfo())
        std::cout << "----> Running Algorithm Tool: " << this->GetInstanceName() << ", " << this->GetType() << std::endl;

    const LArTPCMap &larTPCMap(this->GetPandora().GetGeometry()->GetLArTPCMap());
    const LArTPC *const pFirstLArTPC(larTPCMap.begin()->second);

    float parentMinX(pFirstLArTPC->GetCenterX() - 0.5f * pFirstLArTPC->GetWidthX());
    float parentMaxX(pFirstLArTPC->GetCenterX() + 0.5f * pFirstLArTPC->GetWidthX());
    float parentMinY(pFirstLArTPC->GetCenterY() - 0.5f * pFirstLArTPC->GetWidthY());
    float parentMaxY(pFirstLArTPC->GetCenterY() + 0.5f * pFirstLArTPC->GetWidthY());
    float parentMinZ(pFirstLArTPC->GetCenterZ() - 0.5f * pFirstLArTPC->GetWidthZ());
    float parentMaxZ(pFirstLArTPC->GetCenterZ() + 0.5f * pFirstLArTPC->GetWidthZ());

    for (const LArTPCMap::value_type &mapEntry : larTPCMap)
    {
        const LArTPC *const pLArTPC(mapEntry.second);
        parentMinX = std::min(parentMinX, pLArTPC->GetCenterX() - 0.5f * pLArTPC->GetWidthX());
        parentMaxX = std::max(parentMaxX, pLArTPC->GetCenterX() + 0.5f * pLArTPC->GetWidthX());
        parentMinY = std::min(parentMinY, pLArTPC->GetCenterY() - 0.5f * pLArTPC->GetWidthY());
        parentMaxY = std::max(parentMaxY, pLArTPC->GetCenterY() + 0.5f * pLArTPC->GetWidthY());
        parentMinZ = std::min(parentMinZ, pLArTPC->GetCenterZ() - 0.5f * pLArTPC->GetWidthZ());
        parentMaxZ = std::max(parentMaxZ, pLArTPC->GetCenterZ() + 0.5f * pLArTPC->GetWidthZ());
    }

    m_face_Xa = parentMinX;
    m_face_Xc = parentMaxX;
    m_face_Yb = parentMinY;
    m_face_Yt = parentMaxY;
    m_face_Zu = parentMinZ;
    m_face_Zd = parentMaxZ;

    if (this->GetPandora().GetSettings()->ShouldDisplayAlgorithmInfo())
    {
        std::cout << "RockMuonTaggingTool: Detector boundaries to tag Clear Rock Muons : " <<
                "X [" << m_face_Xa + m_marginX << ","  << m_face_Xc - m_marginX << "] " <<
                "Y [" << m_face_Yb + m_marginY << ", " << m_face_Yt - m_marginY << "] " <<
                "Z [" << m_face_Zu + m_marginZ << ", " << m_face_Zd - m_marginZ << "]" << std::endl;
    }

    int nof_tagged_rockmus = 0;

    // Filter out left over rock muons in the ambigious Pfos list....
    for (const ParticleFlowObject *const pPfo : parentCosmicRayPfos)
    {
        const CRCandidate candidate = CRCandidate(this->GetPandora(), pPfo, 0);
        // Note: 0 is a placeholder for pfoSliceId (unused in CheckIfThroughGoing).
        // CRCandidate is used (instead of pfo) because it runs sliding and gives track start/stop.

        const bool is_throughgoing = this->CheckIfThroughgoing(candidate);

        if (!is_throughgoing)
        {
            // Clean this clear rock mu from ambiguousPfos
            ambiguousPfos.push_back(pPfo);
        }
        else
        {
            nof_tagged_rockmus++;
        }
    }

    if (this->GetPandora().GetSettings()->ShouldDisplayAlgorithmInfo())
        std::cout << "RockMuonTaggingTool: Rock Muon tagged: " << nof_tagged_rockmus << std::endl;
}

//------------------------------------------------------------------------------------------------------------------------------------------
StatusCode RockMuonTaggingTool::ReadSettings(const TiXmlHandle xmlHandle)
{

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MarginZ", m_marginZ));

    return CosmicRayTaggingTool::ReadSettings(xmlHandle);
}

//------------------------------------------------------------------------------------------------------------------------------------------
} // namespace lar_content
