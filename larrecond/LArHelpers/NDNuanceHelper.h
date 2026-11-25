/**
 *  @file   larrecond/LArHelpers/NDNuanceHelper.h
 *
 *  @brief  Header file for the Nuance helper functions
 *
 *  $Log: $
 */
#ifndef PANDORA_ND_NUANCE_HELPER_H
#define PANDORA_ND_NUANCE_HELPER_H 1

#include <string>

namespace lar_nd_reco
{

/**
 *  @brief  Convert the GENIE neutrino reaction string to a Nuance-like integer code
 *
 *  @param  reaction The neutrino reaction string
 *
 *  @return The reaction integer code
 */
int GetNuanceCode(const std::string &reaction);

/**
 *  @brief  Get Nuance reaction string
 *
 *  @param  ccnc Integer specifying CCNC reaction
 *  @param  mode Integer specifying general physics mode
 *
 *  @return Name of the Nuance reaction
 */
std::string GetNuanceReaction(const int ccnc, const int mode);

} // namespace lar_nd_reco

#endif
