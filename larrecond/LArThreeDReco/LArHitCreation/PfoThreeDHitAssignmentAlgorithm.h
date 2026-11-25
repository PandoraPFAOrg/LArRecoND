/**
 *  @file   larrecond/LArThreeDReco/LArHitCreation/PfoThreeDHitAssignmentAlgorithm.h
 *
 *  @brief  Header file for the cluster creation algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_PFO_THREE_D_HIT_ASSIGNMENT_ALGORITHM_H
#define LAR_PFO_THREE_D_HIT_ASSIGNMENT_ALGORITHM_H 1

#include "Objects/CaloHit.h"
#include "Objects/ParticleFlowObject.h"
#include "Pandora/Algorithm.h"

namespace lar_content
{

/**
 *  @brief  PfoThreeDHitAssignmentAlgorithm class
 */
class PfoThreeDHitAssignmentAlgorithm : public pandora::Algorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    PfoThreeDHitAssignmentAlgorithm();

private:
    pandora::StatusCode Run();

    /**
     *  @brief  Assign hits to the pfos
     *
     *  @param  pPfo pointer to the pfo
     *  @param  hits reference to the hit list to add to the pfo
     *  @param  listName name of the ouput 3D cluster list
     */
    void AddHitsToPfo(const pandora::ParticleFlowObject *pPfo, const pandora::CaloHitList &hits, const std::string listName) const;

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    std::string m_inputCaloHitList3DName;              ///< Name of the input 3D calo hit list
    std::vector<std::string> m_inputPfoListNames;      ///< Name of the input pfo list(s)
    std::vector<std::string> m_outputClusterListNames; ///< Name of the output cluster list(s)
};

} // namespace lar_content

#endif // #ifndef LAR_PFO_THREE_D_HIT_ASSIGNMENT_ALGORITHM_H
