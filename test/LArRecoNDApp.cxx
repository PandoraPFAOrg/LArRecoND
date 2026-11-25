/**
 *  @file   test/LArRecoNDApp.cxx
 *
 *  @brief  Implementation of the LArRecoND application
 *
 *  $Log: $
 */

#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <vector>

#include "Helpers/XmlHelper.h"
#include "Xml/tinyxml.h"
#include <Pandora/StatusCodes.h>

#include "larpandoracontent/LArHelpers/LArFileHelper.h"

#include "larrecond/LArControlFlow/MainNDPandora.h"
#include "larrecond/LArObjects/NDParameters.h"

#include "LArRecoNDApp.h"

#ifdef MONITORING
#include "TApplication.h"
#endif

using namespace pandora;
using namespace lar_nd_reco;

typedef std::map<std::string, std::string> InstanceMap;

struct LArRecoNDConfig
{
    // Default config file name
    std::string configFileName = "LArND_TMS.xml";
    // Input file names.
    InstanceMap inputFiles;
    // Geometry file names.
    InstanceMap geometryFiles;
    // Whether the help message was requested.
    bool showHelp = false;
    // TODO: Expand this to match all the expected CLI options for "Main" in the XML.
};

std::tuple<std::string, std::string> ParseInstanceInput(const std::string &input)
{
    const size_t sepPos = input.find_first_of(":");

    if (sepPos == std::string::npos)
    {
        return std::make_tuple("", "");
    }

    std::string instanceName = input.substr(0, sepPos);
    std::string filePath = input.substr(sepPos + 1);
    return std::make_tuple(instanceName, filePath);
}

LArRecoNDConfig ParseCommandLine(int argc, char *argv[])
{
    LArRecoNDConfig options;

    for (int i = 1; i < argc; ++i)
    {
        std::string arg = argv[i];

        if ((arg == "-c" || arg == "--config") && i + 1 < argc)
        {
            options.configFileName = argv[++i];
        }
        else if ((arg == "-i" || arg == "--input") && i + 1 < argc)
        {
            std::string inputArg = argv[++i];
            auto [instanceName, filePath] = ParseInstanceInput(inputArg);

            if (!instanceName.empty() && !filePath.empty())
                options.inputFiles[instanceName] = filePath;
            else
                std::cerr << "Warning: Malformed input argument '" << inputArg << "'. Expected format is InstanceName:FilePath" << std::endl;
        }
        else if ((arg == "-g" || arg == "--geometry") && i + 1 < argc)
        {
            std::string geomArg = argv[++i];
            auto [instanceName, filePath] = ParseInstanceInput(geomArg);

            if (!instanceName.empty() && !filePath.empty())
                options.geometryFiles[instanceName] = filePath;
            else
                std::cerr << "Warning: Malformed geometry argument '" << geomArg << "'. Expected format is InstanceName:FilePath" << std::endl;
        }
        else if (arg == "-h" || arg == "--help")
        {
            std::cout << "Usage: " << argv[0] << " [options]\n"
                      << "Options:\n"
                      << "  -c, --config <file>       Specify the XML config file (default: LArND_TMS.xml)\n"
                      << "  -i, --input <InstanceName:FilePath>  Specify input file for a Pandora instance\n"
                      << "  -g, --geometry <InstanceName:FilePath>  Specify geometry file for a Pandora instance\n"
                      << "  -h, --help                Show this help message\n";
            options.showHelp = true;
            return options;
        }
        else
        {
            std::cerr << "Warning: Unknown or incomplete command line argument: " << arg << std::endl;
        }
    }

    return options;
}

std::string GetValue(InstanceMap &cliOptions, const TiXmlHandle &xmlHandle, const std::string &xmlKey, const std::string &instanceName)
{
    // Check if the CLI option is provided for this instance
    auto cliIt = cliOptions.find(instanceName);
    if (cliIt != cliOptions.end())
        return cliIt->second;

    // If not, check the XML configuration
    TiXmlHandle instanceHandle = xmlHandle.FirstChildElement(instanceName.c_str());
    if (instanceHandle.Element())
    {
        std::string value;
        if (XmlHelper::ReadValue(instanceHandle, xmlKey, value) == STATUS_CODE_SUCCESS)
            return value;
    }

    // If neither is provided, raise an error.
    throw StatusCodeException(STATUS_CODE_NOT_FOUND);
}

int main(int argc, char *argv[])
{
    int errorNo(0);

    LArRecoNDConfig cliOptions = ParseCommandLine(argc, argv);
    if (cliOptions.showHelp)
        return 0;

    const std::string configFileName = cliOptions.configFileName;
    const std::string configFilePath(lar_content::LArFileHelper::FindFileInPath(configFileName, "FW_SEARCH_PATH"));
    std::cout << "XML ConfigFileName = " << configFilePath << std::endl;

#ifdef MONITORING
    std::unique_ptr<TApplication> pTApplication;
#endif

    try
    {
#ifdef MONITORING
        pTApplication = std::make_unique<TApplication>("LArRecoNDApp", &argc, argv);
        pTApplication->SetReturnFromRun(kTRUE);
#endif

        TiXmlDocument xmlDocument;
        if (!xmlDocument.LoadFile(configFilePath.c_str()))
        {
            std::cerr << "Error: Failed to load XML config file: " << configFilePath << std::endl;
            throw StatusCodeException(STATUS_CODE_FAILURE);
        }

        // INFO: The XML document will be the full file...we want inside the <pandora> element.
        const TiXmlHandle xmlDocumentHandle(&xmlDocument);
        const TiXmlHandle xmlHandle(TiXmlHandle(xmlDocumentHandle.FirstChildElement().Element()));

        // Get the main Pandora instance handle
        TiXmlHandle mainHandle(xmlHandle.FirstChildElement("Main"));
        if (!mainHandle.Element())
        {
            std::cerr << "Error: Could not find 'Main' element in XML config." << std::endl;
            throw StatusCodeException(STATUS_CODE_FAILURE);
        }

        NDParameters mainParameters;
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(mainHandle, "Settings", mainParameters.m_settingsFile));
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(mainHandle, "EventsToProcess", mainParameters.m_nEventsToProcess));
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(mainHandle, "EventsToSkip", mainParameters.m_nEventsToSkip));
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(mainHandle, "MaxNHits", mainParameters.m_maxNHits));
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(mainHandle, "MinNHits", mainParameters.m_minNHits));

        std::string viewOption, recoOption;
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(mainHandle, "ViewOption", viewOption));
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(mainHandle, "RecoOption", recoOption));

        mainParameters.SetViewOption(viewOption);
        if (!mainParameters.SetRecoOption(recoOption))
            throw StatusCodeException(STATUS_CODE_FAILURE);

        MainNDPandora mainND("Main", mainParameters);

        std::vector<std::string> instances;
        if (mainHandle.FirstChildElement("Instances").Element())
            XmlHelper::ReadVectorOfValues(mainHandle, "Instances", instances);
        else
            std::cout << "No additional Pandora instances specified in the XML config file." << std::endl;

        // Add the other Pandora instances
        for (const std::string &instanceName : instances)
        {
            std::cout << "Setting up Pandora instance : " << instanceName << std::endl;

            // Check an instance with this name exists in the XML file
            if (!xmlHandle.FirstChildElement(instanceName.c_str()).Element())
            {
                std::cerr << "Error: Could not find '" << instanceName << "' element in XML config." << std::endl;
                throw StatusCodeException(STATUS_CODE_FAILURE);
            }

            const auto instanceHandle = xmlHandle.FirstChildElement(instanceName.c_str());

            NDParameters NDPars(mainParameters);

            std::string volTypeStr;
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(instanceHandle, "VolType", volTypeStr));
            NDPars.m_volType = NDPars.GetVolEnum(volTypeStr);

            std::string dataFormatStr;
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(instanceHandle, "DataFormat", dataFormatStr));
            NDPars.m_dataFormat = NDPars.GetDataEnum(dataFormatStr);

            std::string settingsFile;
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(instanceHandle, "Settings", settingsFile));
            NDPars.m_settingsFile = lar_content::LArFileHelper::FindFileInPath(settingsFile, "FW_SEARCH_PATH");

            // Process any options that may have CLI overrides, otherwise use the XML values
            NDPars.m_inputFileName = GetValue(cliOptions.inputFiles, instanceHandle, "InputFile", instanceName);
            NDPars.m_geomFileName = GetValue(cliOptions.geometryFiles, instanceHandle, "GeomFile", instanceName);

            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(instanceHandle, "InputTree", NDPars.m_inputTreeName));
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(instanceHandle, "GeomManager", NDPars.m_geomManagerName));
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(instanceHandle, "TPCName", NDPars.m_tpcName));
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(instanceHandle, "LengthScale", NDPars.m_lengthScale));
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(instanceHandle, "EnergyScale", NDPars.m_energyScale));

            // This instance may have different reco options compared to the main one.
            std::string instanceRecoOption;
            if (XmlHelper::ReadValue(instanceHandle, "RecoOption", instanceRecoOption) == STATUS_CODE_SUCCESS)
            {
                if (!NDPars.SetRecoOption(instanceRecoOption))
                    throw StatusCodeException(STATUS_CODE_FAILURE);
            }

            // Instance may use different projection views as well
            std::string instanceViewOption;
            if (XmlHelper::ReadValue(instanceHandle, "ViewOption", instanceViewOption) == STATUS_CODE_SUCCESS)
                NDPars.SetViewOption(instanceViewOption);

            // Set the event info using the main Pandora instance.
            // This assumes the input files for each Pandora instance have
            // the same event number ordering
            NDPars.m_nEventsToProcess = mainParameters.m_nEventsToProcess;
            NDPars.m_nEventsToSkip = mainParameters.m_nEventsToSkip;
            NDPars.m_maxNHits = mainParameters.m_maxNHits;
            NDPars.m_minNHits = mainParameters.m_minNHits;

            std::cout << "Settings file = " << NDPars.m_settingsFile << std::endl;
            std::cout << "VolType = " << NDPars.m_volType << ", DataFormat = " << NDPars.m_dataFormat << std::endl;
            std::cout << "Input file = " << NDPars.m_inputFileName << ", tree = " << NDPars.m_inputTreeName << std::endl;
            std::cout << "Geometry file = " << NDPars.m_geomFileName << ", manager = " << NDPars.m_geomManagerName << std::endl;
            std::cout << "TPC volume name = " << NDPars.m_tpcName << std::endl;
            std::cout << "Length scale (cm) = " << NDPars.m_lengthScale << ", Energy scale (GeV) = " << NDPars.m_energyScale << std::endl;

            mainND.AddPandoraInstance(instanceName, NDPars);
        }

        // Create the TPC geometry using the tpcNames in the NDParameters geometry ROOT file.
        // This creates TPCs for both the main & added Pandora instances
        mainND.CreatePandoraTPCs();

        // Create detector gaps
        mainND.CreatePandoraDetectorGaps();

        // Setup external & algorithm parameters for all Pandora instances
        mainND.ConfigurePandoraInstances();

        // Setup the event inputs
        mainND.ConfigureEventInputs();

        // Process the events
        mainND.ProcessEvents();

        std::cout << "Done" << std::endl;
    }
    catch (const StatusCodeException &statusCodeException)
    {
        std::cerr << "Pandora StatusCodeException: " << statusCodeException.ToString() << statusCodeException.GetBackTrace() << std::endl;
        errorNo = 1;
    }
    catch (...)
    {
        std::cerr << "Misconfigurated XML file/parameters or unknown exception" << std::endl;
        errorNo = 1;
    }

    return errorNo;
}
