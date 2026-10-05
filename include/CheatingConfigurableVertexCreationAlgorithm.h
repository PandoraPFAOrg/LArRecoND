/**
 *  @file   CheatingConfigurableVertexCreationAlgorithm.h
 *
 *  @brief  Header file for the configurable cheating vertex creation algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_CHEATING_CONFIGURABLE_VERTEX_CREATION_ALGORITHM_H
#define LAR_CHEATING_CONFIGURABLE_VERTEX_CREATION_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

namespace lar_content
{

/**
 *  @brief  CheatingConfigurableVertexCreationAlgorithm::Algorithm class
 */
class CheatingConfigurableVertexCreationAlgorithm : public pandora::Algorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    CheatingConfigurableVertexCreationAlgorithm();

private:
    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    std::string m_outputVertexListName; ///< The name under which to save the output vertex list
    bool m_replaceCurrentVertexList;    ///< Whether to replace the current vertex list with the output list
    float m_vertexXCorrection;          ///< The vertex x correction, added to reported mc neutrino endpoint x value, in cm
    bool m_useNeutrinoInteractionVertex;///< Whether to use the neutrino interaction vertex. If set to false, use the start of the highest momentum MC particle. 
};

} // namespace lar_content

#endif // #ifndef LAR_CHEATING_CONFIGURABLE_VERTEX_CREATION_ALGORITHM_H
