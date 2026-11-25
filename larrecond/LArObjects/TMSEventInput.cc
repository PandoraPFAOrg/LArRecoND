/**
 *  @file   larrecond/LArObjects/TMSEventInput.cc
 *
 *  @brief  Implementation of the TMS event input
 *
 *  $Log: $
 */

#include "Api/PandoraApi.h"
#include "Managers/PluginManager.h"

#include "larpandoracontent/LArObjects/LArCaloHit.h"
#include "larpandoracontent/LArObjects/LArMCParticle.h"
#include "larpandoracontent/LArPlugins/LArPseudoLayerPlugin.h"
#include "larpandoracontent/LArPlugins/LArRotationalTransformationPlugin.h"

#include "larrecond/LArObjects/TMSEventInput.h"

namespace lar_nd_reco
{

TMSEventInput::TMSEventInput(const pandora::Pandora *pPandora, const NDParameters &parameters, const NDSimpleTPCGeom &tpcGeom) :
    NDEventInput(pPandora, parameters, tpcGeom),
    m_tms{}
{
    if (this->ShouldDisplayInfo())
        std::cout << "TMSEventInput nEntries = " << m_nEntries << std::endl;
    m_tms = std::make_unique<TMS>(m_pTree);
}

//------------------------------------------------------------------------------------------------------------------------------------------

TMSEventInput::~TMSEventInput()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TMSEventInput::CreateMCParticles(int) const
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

int TMSEventInput::CreateCaloHits(int eventIndex, int startHitIndex) const
{
    if (this->ShouldDisplayInfo())
        std::cout << "TMSEventInput eventIndex = " << eventIndex << ", start = " << startHitIndex << std::endl;

    m_tms->GetEntry(eventIndex);

    const int nHits = m_tms->m_NRecoHits;
    if (this->ShouldDisplayInfo())
        std::cout << "TMS nHits = " << nHits << ", startHitIndex = " << startHitIndex << std::endl;

    int hitCounter{startHitIndex};

    // Loop over the hits and make them into CaloHits
    for (int iHit = 0; iHit < nHits; ++iHit)
    {
        //const float voxelX = m_tms->m_RecoHitX[iHit] * m_parameters.m_lengthScale;
        //const float voxelY = m_tms->m_RecoHitY[iHit] * m_parameters.m_lengthScale;
        //const float voxelZ = m_tms->m_RecoHitZ[iHit] * m_parameters.m_lengthScale;
        //const float voxelE = m_tms->m_RecoHitE[iHit] * m_parameters.m_energyScale;

        const float voxelX = m_tms->m_TrueHitX[iHit] * m_parameters.m_lengthScale;
        const float voxelY = m_tms->m_TrueHitY[iHit] * m_parameters.m_lengthScale;
        const float voxelZ = m_tms->m_TrueHitZ[iHit] * m_parameters.m_lengthScale;
        const float voxelE = m_tms->m_TrueHitE[iHit] * m_parameters.m_energyScale;

        //std::cout<<"CaloHit1: "<<voxelX1<<", "<<voxelY1<<", "<<voxelZ1<<", "<<voxelE1<<std::endl;

        // Skip this hit if its coordinates or energy are NaNs
        if (std::isnan(voxelX) || std::isnan(voxelY) || std::isnan(voxelZ) || std::isnan(voxelE))
        {
            std::cout << "Ignoring hit " << iHit << " which contains NaNs: (" << voxelX << ", " << voxelY << ", " << voxelZ
                      << "), E = " << voxelE << std::endl;
            continue;
        }

        const pandora::CartesianVector voxelPos(voxelX, voxelY, voxelZ);
        const bool ignoreY{true};
        const int tpcID(m_tpcGeom.GetTPCNumber(voxelPos, ignoreY));

        if (tpcID < 0)
            std::cout << "TMSEventInput::CreateCaloHits: hit " << iHit << " is outside of TPC volume: (" << voxelX << ", " << voxelY << ", "
                      << voxelZ << ")" << std::endl;

        const lar_content::LArCaloHitParameters caloHitParameters(this->MakeCaloHitParameters(voxelPos, voxelE, tpcID));

        // TODO: Add the MC link once the TMS truth information is read
        this->CreateVoxelCaloHits(caloHitParameters, voxelPos, hitCounter, nullptr);
    } // hit loop

    if (this->ShouldDisplayInfo())
        std::cout << "startHitIndex = " << startHitIndex << ", hitCounter = " << hitCounter << std::endl;
    return hitCounter;
}

//------------------------------------------------------------------------------------------------------------------------------------------

int TMSEventInput::CreateTracks(int, int) const
{
    return 0;
}

} // namespace lar_nd_reco
