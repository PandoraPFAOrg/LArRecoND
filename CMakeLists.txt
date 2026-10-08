cmake_minimum_required(VERSION 3.20 FATAL_ERROR)
project(LArRecoND VERSION 01.05.00 LANGUAGES CXX)

# Automatically add the path of the linked libraries to the RPATH
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# Don't replace the build-tree RPATH with the install-tree RPATH
set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)

# Enforce out-of-source builds, a good practice from the original file
if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
    message(FATAL_ERROR "LArRecoND requires an out-of-source build.")
endif()

if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    set(CMAKE_INSTALL_PREFIX "${CMAKE_BINARY_DIR}/install" CACHE PATH "Default install path for ${PROJECT_NAME}" FORCE)
endif()

# Build Options
option(PANDORA_LIBTORCH "Build with LibTorch-dependent libraries" OFF)
option(PANDORA_MONITORING "Build with PandoraMonitoring support" ON)
option(LArRecoND_BUILD_DOCS "Build documentation for ${PROJECT_NAME}" OFF)

# Dependencies
if (NOT TARGET PandoraPFA::PandoraSDK)
    find_package(PandoraSDK 05.00.00 REQUIRED)
endif()
if (NOT TARGET PandoraPFA::LArContent)
    find_package(LArContent 05.00.00 REQUIRED)
endif()

if(PANDORA_LIBTORCH)
    find_package(LArDLContent 05.00.00 REQUIRED)
endif()

if(PANDORA_MONITORING)
    if (NOT TARGET PandoraPFA::PandoraMonitoring)
        find_package(PandoraMonitoring 05.00.00 REQUIRED)
    endif()
    find_package(ROOT 6.18.04 REQUIRED COMPONENTS Eve Geom RGL EG)

    # By marking ROOT includes as system headers, we suppress any warnings from them.
    include_directories(SYSTEM ${ROOT_INCLUDE_DIRS})
endif()

if(USE_EDEPSIM)
    find_package(EDepSim)
endif()

# Common compile options
set(LAR_RECO_COMPILE_OPTIONS
    -Wall
    -Wextra
    -Werror
    -pedantic
    -Wno-long-long
    -Wno-sign-compare
    -Wshadow
    -fno-strict-aliasing
)

# Build products
file(GLOB_RECURSE LAR_RECO_ND_SRCS CONFIGURE_DEPENDS "${PROJECT_SOURCE_DIR}/src/*.cc")

add_library(${PROJECT_NAME} SHARED ${LAR_RECO_ND_SRCS})

set_target_properties(${PROJECT_NAME} PROPERTIES
    CXX_STANDARD 17
    CXX_STANDARD_REQUIRED ON
    VERSION ${PROJECT_VERSION}
    SOVERSION ${PROJECT_VERSION_MAJOR}
)

target_include_directories(${PROJECT_NAME} PUBLIC ${PROJECT_SOURCE_DIR}/include)

target_compile_options(${PROJECT_NAME} PRIVATE ${LAR_RECO_COMPILE_OPTIONS})

target_link_libraries(${PROJECT_NAME} PRIVATE
    PandoraPFA::PandoraSDK
    PandoraPFA::LArContent
)

if(PANDORA_LIBTORCH)
    target_link_libraries(${PROJECT_NAME} PRIVATE PandoraPFA::LArDLContent)
    target_compile_definitions(${PROJECT_NAME} PRIVATE LIBTORCH_DL=1)
endif()

if(PANDORA_MONITORING)
    target_link_libraries(${PROJECT_NAME} PRIVATE PandoraPFA::PandoraMonitoring)
    target_compile_definitions(${PROJECT_NAME} PRIVATE MONITORING)
endif()

if(USE_EDEPSIM)
    target_compile_definitions(${PROJECT_NAME} PRIVATE USE_EDEPSIM)
endif()

# Executables
set(LAR_RECO_EXECUTABLES
    PandoraInterface
    PandoraOuterface
)

foreach(executable_name IN LISTS LAR_RECO_EXECUTABLES)
    add_executable(${executable_name} test/${executable_name}.cxx)

    set_target_properties(${executable_name} PROPERTIES
        CXX_STANDARD 17
        CXX_STANDARD_REQUIRED ON
        BUILD_RPATH "$<TARGET_FILE_DIR:${PROJECT_NAME}>"
    )

    target_compile_options(${executable_name} PRIVATE ${LAR_RECO_COMPILE_OPTIONS})

    target_link_libraries(${executable_name} PRIVATE
        ${PROJECT_NAME}
        PandoraPFA::PandoraSDK
        PandoraPFA::LArContent
    )

    if(PANDORA_LIBTORCH)
        target_link_libraries(${executable_name} PRIVATE PandoraPFA::LArDLContent)
        target_compile_definitions(${executable_name} PRIVATE LIBTORCH_DL=1)
    endif()

    if(PANDORA_MONITORING)
        target_link_libraries(${executable_name} PRIVATE
            PandoraPFA::PandoraMonitoring
            ROOT::Eve
            ROOT::Geom
            ROOT::RGL
            ROOT::EG
        )
        target_compile_definitions(${executable_name} PRIVATE MONITORING)
    endif()

    if(USE_EDEPSIM)
        target_compile_definitions(${executable_name} PRIVATE USE_EDEPSIM)
    endif()
endforeach()

if(PANDORA_MONITORING AND USE_EDEPSIM)
    target_link_libraries(PandoraInterface PRIVATE EDepSim::edepsim_io)
endif()

# Optional documents
if(LArRecoND_BUILD_DOCS)
    add_subdirectory(doc)
endif()

# Installation
install(DIRECTORY include/ DESTINATION include COMPONENT Development FILES_MATCHING PATTERN "*.h")
install(TARGETS ${PROJECT_NAME} DESTINATION lib COMPONENT Runtime)
foreach(executable_name IN LISTS LAR_RECO_EXECUTABLES)
    install(TARGETS ${executable_name} DESTINATION bin PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)
endforeach()
