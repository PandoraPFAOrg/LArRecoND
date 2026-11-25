/**
 *  @file   larrecond/LArObjects/NDParameters.h
 *
 *  @brief  Header file for ND parameters
 *
 *  $Log: $
 */
#ifndef PANDORA_ND_PARAMETERS_H
#define PANDORA_ND_PARAMETERS_H 1

#include <string>

namespace lar_nd_reco
{

class NDParameters
{
public:
    // Near Detector volume type
    enum VolType
    {
        Main = 0,
        LArND = 1,
        TMS = 2
    };

    // Data format enumeration
    enum DataFormat
    {
        SP = 0,
        SPMC = 1,
        TMSMC = 2
    };

    /**
     *  @brief  Default constructor
     */
    NDParameters();

    /**
     *  @brief  Set various parameters depending on what reco option we have
     *
     *  @param  recoOption the reco option string
     *  @return if the option has been recognised
     *
     */
    inline bool SetRecoOption(const std::string &recoOption);

    /**
     *  @brief  Set the projection view option
     *
     *  @param  viewOption the view option string
     *
     */
    inline void SetViewOption(const std::string &viewOption);

    /**
     *  @brief  Get the volume enum integer given the name
     *
     *  @param  volName the volume type name
     *  @return the volume enumeration integer
     */
    inline NDParameters::VolType GetVolEnum(const std::string &volName);

    /**
     *  @brief  Get the data format enum integer given the name
     *
     *  @param  formatName the data format name
     *  @return the data format enumeration integer
     */
    inline NDParameters::DataFormat GetDataEnum(const std::string &formatName);

    std::string m_settingsFile;    ///< The path to the Pandora settings file
    std::string m_inputFileName;   ///< The path to the input file containing the events
    std::string m_inputTreeName;   ///< The name of the event TTree
    std::string m_geomFileName;    ///< The ROOT file name containing the TGeoManager info
    std::string m_geomManagerName; ///< The name of the TGeoManager
    std::string m_tpcName;         ///< The name of the TPC volume nodes

    VolType m_volType;       ///< The Near Detector volume Id
    DataFormat m_dataFormat; ///< The input data format
    bool m_use3D;            ///< Create 3D LArCaloHits
    bool m_useLArTPC;        ///< Create LArTPC-like LArCaloHits with U,V,W views
    float m_voxelWidth;      ///< Voxel box width (cm), default = 0.4 cm
    float m_lengthScale;     ///< The scaling factor to set all lengths to cm
    float m_energyScale;     ///< The scaling factor to set all energies to GeV

    int m_nEventsToProcess; ///< The number of events to process (default all events)
    int m_nEventsToSkip;    ///< The number of events to skip
    int m_maxNHits;         ///< The max number of hits to process in a given event (default all)
    int m_minNHits;         ///< The minimum number of hits for processing an event (default = 2)

    bool m_shouldRunAllHitsCosmicReco;  ///< Whether to run all hits cosmic-ray reconstruction
    bool m_shouldRunStitching;          ///< Whether to stitch cosmic-ray muons crossing between volumes
    bool m_shouldRunCosmicHitRemoval;   ///< Whether to remove hits from tagged cosmic-rays
    bool m_shouldRunSlicing;            ///< Whether to slice events into separate regions for processing
    bool m_shouldRunNeutrinoRecoOption; ///< Whether to run neutrino reconstruction for each slice
    bool m_shouldRunCosmicRecoOption;   ///< Whether to run cosmic-ray reconstruction for each slice
    bool m_shouldPerformSliceId;        ///< Whether to identify slices and select most appropriate pfos
    bool m_printOverallRecoStatus;      ///< Whether to print current operation status messages

    // Internal maps associating strings to the enumerations
    typedef std::map<std::string, NDParameters::VolType> VolTypeMap;
    typedef std::map<std::string, NDParameters::DataFormat> DataFormatMap;
    VolTypeMap m_VolTypeMap;
    DataFormatMap m_DataFormatMap;
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline NDParameters::NDParameters() :
    m_settingsFile(""),
    m_inputFileName(""),
    m_inputTreeName(""),
    m_geomFileName(""),
    m_geomManagerName(""),
    m_tpcName(""),
    m_volType(NDParameters::VolType::Main),
    m_dataFormat(NDParameters::DataFormat::SP),
    m_use3D(true),
    m_useLArTPC(true),
    m_voxelWidth(0.4f),
    m_lengthScale(1.f),
    m_energyScale(1.f),
    m_nEventsToProcess(-1),
    m_nEventsToSkip(0),
    m_maxNHits(-1),
    m_minNHits(2),
    m_shouldRunAllHitsCosmicReco(true),
    m_shouldRunStitching(true),
    m_shouldRunCosmicHitRemoval(true),
    m_shouldRunSlicing(true),
    m_shouldRunNeutrinoRecoOption(true),
    m_shouldRunCosmicRecoOption(true),
    m_shouldPerformSliceId(true),
    m_printOverallRecoStatus(false),
    m_VolTypeMap{},
    m_DataFormatMap{}
{
    m_VolTypeMap.insert(std::make_pair("Main", VolType::Main));
    m_VolTypeMap.insert(std::make_pair("LArND", VolType::LArND));
    m_VolTypeMap.insert(std::make_pair("TMS", VolType::TMS));

    m_DataFormatMap.insert(std::make_pair("SP", DataFormat::SP));
    m_DataFormatMap.insert(std::make_pair("SPMC", DataFormat::SPMC));
    m_DataFormatMap.insert(std::make_pair("TMSMC", DataFormat::TMSMC));
}

inline bool NDParameters::SetRecoOption(const std::string &recoOption)
{
    std::cout << "NDParameters::SetRecoOption is " << recoOption << std::endl;
    std::string chosenRecoOption(recoOption);
    std::transform(chosenRecoOption.begin(), chosenRecoOption.end(), chosenRecoOption.begin(), ::tolower);
    bool processed(true);

    if ("full" == chosenRecoOption)
    {
        m_shouldRunAllHitsCosmicReco = true;
        m_shouldRunStitching = true;
        m_shouldRunCosmicHitRemoval = true;
        m_shouldRunSlicing = true;
        m_shouldRunNeutrinoRecoOption = true;
        m_shouldRunCosmicRecoOption = true;
        m_shouldPerformSliceId = true;
    }
    else if ("allhitscr" == chosenRecoOption)
    {
        m_shouldRunAllHitsCosmicReco = true;
        m_shouldRunStitching = true;
        m_shouldRunCosmicHitRemoval = false;
        m_shouldRunSlicing = false;
        m_shouldRunNeutrinoRecoOption = false;
        m_shouldRunCosmicRecoOption = false;
        m_shouldPerformSliceId = false;
    }
    else if ("nostitchingcr" == chosenRecoOption)
    {
        m_shouldRunAllHitsCosmicReco = false;
        m_shouldRunStitching = false;
        m_shouldRunCosmicHitRemoval = false;
        m_shouldRunSlicing = false;
        m_shouldRunNeutrinoRecoOption = false;
        m_shouldRunCosmicRecoOption = true;
        m_shouldPerformSliceId = false;
    }
    else if ("allhitsnu" == chosenRecoOption)
    {
        m_shouldRunAllHitsCosmicReco = false;
        m_shouldRunStitching = false;
        m_shouldRunCosmicHitRemoval = false;
        m_shouldRunSlicing = false;
        m_shouldRunNeutrinoRecoOption = true;
        m_shouldRunCosmicRecoOption = false;
        m_shouldPerformSliceId = false;
    }
    else if ("crremhitsslicecr" == chosenRecoOption)
    {
        m_shouldRunAllHitsCosmicReco = true;
        m_shouldRunStitching = true;
        m_shouldRunCosmicHitRemoval = true;
        m_shouldRunSlicing = true;
        m_shouldRunNeutrinoRecoOption = false;
        m_shouldRunCosmicRecoOption = true;
        m_shouldPerformSliceId = false;
    }
    else if ("crremhitsslicenu" == chosenRecoOption)
    {
        m_shouldRunAllHitsCosmicReco = true;
        m_shouldRunStitching = true;
        m_shouldRunCosmicHitRemoval = true;
        m_shouldRunSlicing = true;
        m_shouldRunNeutrinoRecoOption = true;
        m_shouldRunCosmicRecoOption = false;
        m_shouldPerformSliceId = false;
    }
    else if ("allhitsslicecr" == chosenRecoOption)
    {
        m_shouldRunAllHitsCosmicReco = false;
        m_shouldRunStitching = false;
        m_shouldRunCosmicHitRemoval = false;
        m_shouldRunSlicing = true;
        m_shouldRunNeutrinoRecoOption = false;
        m_shouldRunCosmicRecoOption = true;
        m_shouldPerformSliceId = false;
    }
    else if ("allhitsslicenu" == chosenRecoOption)
    {
        m_shouldRunAllHitsCosmicReco = false;
        m_shouldRunStitching = false;
        m_shouldRunCosmicHitRemoval = false;
        m_shouldRunSlicing = true;
        m_shouldRunNeutrinoRecoOption = true;
        m_shouldRunCosmicRecoOption = false;
        m_shouldPerformSliceId = false;
    }
    else
    {
        std::cout << "Unrecognized reconstruction option: " << recoOption << std::endl;
        std::cout << "Choose: Full, AllHitsSliceNu, AllHitsNu, AllHitsCR, CRRemHitsSliceCR, CRRemHitsSliceNu, or AllHitsSliceCR" << std::endl;
        processed = false;
    }

    return processed;
}

inline void NDParameters::SetViewOption(const std::string &viewOption)
{
    std::cout << "NDParameters::SetViewOption is " << viewOption << std::endl;
    std::string chosenViewOption(viewOption);
    std::transform(chosenViewOption.begin(), chosenViewOption.end(), chosenViewOption.begin(), ::tolower);

    if (chosenViewOption == "3d")
    {
        // 3D hits only
        std::cout << "Using 3D hits" << std::endl;
        m_useLArTPC = false;
        m_use3D = true;
    }
    else if (chosenViewOption == "lartpc")
    {
        // LArTPC 2D projected hits only
        std::cout << "Using LArTPC projections" << std::endl;
        m_useLArTPC = true;
        m_use3D = false;
    }
    else
    {
        // Both LArTPC 2D projections and 3D hits (default)
        std::cout << "Using LArTPC projections _and_ 3D hits" << std::endl;
        m_useLArTPC = true;
        m_use3D = true;
    }
}

inline NDParameters::VolType NDParameters::GetVolEnum(const std::string &volName)
{
    NDParameters::VolTypeMap::const_iterator iter = m_VolTypeMap.find(volName);
    if (m_VolTypeMap.end() == iter)
        throw pandora::StatusCodeException(pandora::STATUS_CODE_NOT_FOUND);

    return iter->second;
}

inline NDParameters::DataFormat NDParameters::GetDataEnum(const std::string &volName)
{
    NDParameters::DataFormatMap::const_iterator iter = m_DataFormatMap.find(volName);
    if (m_DataFormatMap.end() == iter)
        throw pandora::StatusCodeException(pandora::STATUS_CODE_NOT_FOUND);

    return iter->second;
}

} // namespace lar_nd_reco

#endif
