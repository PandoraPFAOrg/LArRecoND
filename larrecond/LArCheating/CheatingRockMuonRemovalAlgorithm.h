/**
 *  @file   larrecond/LArCheating/CheatingRockMuonRemovalAlgorithm.h
 *
 *  @brief  Header file for the cheating rock muon removal algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_CHEATING_ROCK_MUON_REMOVAL_ALGORITHM_H
#define LAR_CHEATING_ROCK_MUON_REMOVAL_ALGORITHM_H 1

#include "Objects/MCParticle.h"
#include "Pandora/Algorithm.h"
#include "Pandora/Pandora.h"

namespace lar_content
{

/**
 *  @brief  CheatingRockMuonRemovalAlgorithm class
 */
class CheatingRockMuonRemovalAlgorithm : public pandora::Algorithm
{
public:
    /**
     *  @brief  Is the specified MC particle a rock muon?
     *
     *  @param  pandora the pandora instance
     *  @param  pMCParticle the MC particle pointer
     *
     *  @return whether the MC particle is a rock muon
     */
    static bool IsRockMuon(const pandora::Pandora &pandora, const pandora::MCParticle *const pMCParticle);

private:
    pandora::StatusCode Run();

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    std::string m_inputCaloHitListName;    ///< Name of the calo hit list to load
    std::string m_rockMuonCaloHitListName; ///< Name of the rock muon calo hit list to create
    std::string m_neutrinoCaloHitListName; ///< Name of the neutrino calo hit list to create
    std::string m_inputMCParticleListName; ///< Name of the MC particle list to use
};

} // namespace lar_content

#endif // #ifndef LAR_CHEATING_ROCK_MUON_REMOVAL_ALGORITHM_H
