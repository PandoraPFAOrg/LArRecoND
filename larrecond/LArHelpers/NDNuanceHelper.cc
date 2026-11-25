/**
 *  @file   larrecond/LArHelpers/NDNuanceHelper.cc
 *
 *  @brief  Implementation of the Nuance helper functions
 *
 *  $Log: $
 */

#include <map>
#include <string>

#include "larrecond/LArHelpers/NDNuanceHelper.h"

namespace lar_nd_reco
{

int GetNuanceCode(const std::string &reaction)
{
    // The GENIE reaction string (also stored by edep-sim) is created using
    // https://github.com/GENIE-MC/Generator/blob/master/src/Framework/Interaction/Interaction.cxx#L249
    // String format is "nu:PDGId;tgt:PDGId;N:PDGId;proc:interactionType,scattering;", e.g.
    //                  "nu:14;tgt:1000180400;N:2112;proc:Weak[CC],QES;"

    // GENIE scattering codes:
    // https://github.com/GENIE-MC/Generator/blob/master/src/Framework/Interaction/ScatteringType.h

    // Nuance codes:
    // https://internal.dunescience.org/doxygen/MCNeutrino_8h_source.html

    // GENIE conversion code for RooTracker output files:
    // https://github.com/GENIE-MC/Generator/blob/master/src/contrib/t2k/neut_code_from_rootracker.C
    // Similar code is available here (Neut reaction code):
    // https://internal.dunescience.org/doxygen/namespacegenie_1_1utils_1_1ghep.html

    // For now, just set the basic reaction types, excluding any specific final states:
    // https://github.com/GENIE-MC/Generator/blob/master/src/contrib/t2k/neut_code_from_rootracker.C#L276
    int code(1000);

    const bool is_cc = (reaction.find("Weak[CC]") != std::string::npos); // weak charged-current
    const bool is_nc = (reaction.find("Weak[NC]") != std::string::npos); // weak neutral-current
    // const bool is_charm = (reaction.find("charm")    != std::string::npos); // charm production
    const bool is_qel = (reaction.find("QES") != std::string::npos);   // quasi-elastic scattering
    const bool is_dis = (reaction.find("DIS") != std::string::npos);   // deep inelastic scattering
    const bool is_res = (reaction.find("RES") != std::string::npos);   // resonance
    const bool is_cohpi = (reaction.find("COH") != std::string::npos); // coherent pi
    const bool is_ve = (reaction.find("NuEEL") != std::string::npos);  // nu e elastic
    const bool is_imd = (reaction.find("IMD") != std::string::npos);   // inverse mu decay
    const bool is_mec = (reaction.find("MEC") != std::string::npos);   // meson exchange current

    if (is_qel)
    {
        code = 0;
        if (is_cc)
            code = 1001;
        else if (is_nc)
            code = 1002;
    }
    else if (is_dis)
    {
        code = 2;
        if (is_cc)
            code = 1091;
        else if (is_nc)
            code = 1092;
    }
    else if (is_res)
        code = 1;
    else if (is_cohpi)
    {
        code = 3;
        if (is_qel)
            code = 4;
    }
    else if (is_ve)
        code = 1098;
    else if (is_imd)
        code = 1099;
    else if (is_mec)
        code = 10;

    return code;
}

//------------------------------------------------------------------------------------------------------------------------------------------

std::string GetNuanceReaction(const int ccnc, const int mode)
{
    // GENIE interaction mode -> reaction process string
    static const std::map<int, std::string> modeToProcess = {{0, "QES"}, {1, "RES"}, {2, "DIS"}, {3, "COH"}, {5, "NuEEL"}, {6, "IMD"}, {10, "MEC"}};

    std::string reaction;

    const auto modeIter = modeToProcess.find(mode);
    if (modeIter != modeToProcess.end())
        reaction = modeIter->second;

    if (ccnc == 0)
        reaction += "Weak[CC]";
    else if (ccnc == 1)
        reaction += "Weak[NC]";

    return reaction;
}

} // namespace lar_nd_reco
