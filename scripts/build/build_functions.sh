#!/bin/bash

# Shared error handler
RunWithError() {
    "$@"
    local status=$?
    if [ $status -ne 0 ]; then
        echo "Error running command: $@"
        exit $status
    fi
}

# Shared PANDORA_PROJECT_DIR inference logic
GetPandoraProjectDir() {
    if [ -z "$PANDORA_PROJECT_DIR" ]; then
        if [ -d "$PWD/LArRecoND" ]; then
            export PANDORA_PROJECT_DIR=$PWD
        elif [ -d "$PWD/../LArRecoND" ]; then
            export PANDORA_PROJECT_DIR=$PWD/..
        elif [ -d "$PWD/../../LArRecoND" ]; then
            export PANDORA_PROJECT_DIR=$PWD/../..
        else
            echo "PANDORA_PROJECT_DIR is not set and cannot be inferred from the current working directory."
            echo "Please set the PANDORA_PROJECT_DIR environment variable to the Pandora project directory."
            exit 1
        fi
    fi
}

# Define any common build flags here.
COMMON_CMAKE_FLAGS="-DCMAKE_BUILD_TYPE=${PANDORA_BUILD_TYPE:-RelWithDebInfo}"
COMMON_CMAKE_FLAGS="${COMMON_CMAKE_FLAGS} -DCMAKE_INSTALL_LIBDIR=lib"
COMMON_CMAKE_FLAGS="${COMMON_CMAKE_FLAGS} -DCMAKE_INSTALL_PREFIX=install"
COMMON_CMAKE_FLAGS="${COMMON_CMAKE_FLAGS} -DCMAKE_INSTALL_RPATH_USE_LINK_PATH=TRUE"
MAKE_COMMAND="cmake --build . --target install -j4"

BuildPandoraPFA() {

    if [ -d "$PANDORA_PROJECT_DIR/PandoraPFA" ]; then
        echo "PandoraPFA directory already exists. Skipping!"
        return
    fi

    RunWithError cd $PANDORA_PROJECT_DIR
    RunWithError git clone https://github.com/PandoraPFA/PandoraPFA.git
    RunWithError cd PandoraPFA
    RunWithError git checkout $PANDORA_PFA_VERSION
}

BuildEigen3() {
    RunWithError cd $PANDORA_PROJECT_DIR

    # Download Eigen tarball, rather than git clone, as the repo history is massive.
    if [ ! -d "eigen" ]; then
        echo "Downloading Eigen ${EIGEN_VERSION}..."
        EIGEN_TARBALL="eigen-${EIGEN_VERSION}.tar.gz"
        EIGEN_URL="https://gitlab.com/libeigen/eigen/-/archive/${EIGEN_VERSION}/${EIGEN_TARBALL}"

        RunWithError curl -L -o $EIGEN_TARBALL $EIGEN_URL
        RunWithError tar -xzf $EIGEN_TARBALL
        RunWithError rm $EIGEN_TARBALL
        RunWithError mv eigen-${EIGEN_VERSION} eigen
    fi

    RunWithError cd eigen
    RunWithError mkdir build
    RunWithError cd build
    EIGEN_BUILD_FLAGS="-DBUILD_TESTING=OFF"
    EIGEN_BUILD_FLAGS="${EIGEN_BUILD_FLAGS} -DEIGEN_BUILD_TESTING=OFF"
    EIGEN_BUILD_FLAGS="${EIGEN_BUILD_FLAGS} -DEIGEN_BUILD_DOC=OFF"
    EIGEN_BUILD_FLAGS="${EIGEN_BUILD_FLAGS} -DEIGEN_BUILD_BLAS=OFF"
    EIGEN_BUILD_FLAGS="${EIGEN_BUILD_FLAGS} -DEIGEN_BUILD_LAPACK=OFF"
    EIGEN_BUILD_FLAGS="${EIGEN_BUILD_FLAGS} -DEIGEN_BUILD_DEMOS=OFF"
    EIGEN_BUILD_FLAGS="${EIGEN_BUILD_FLAGS} -DCMAKE_INSTALL_PREFIX=install"
    EIGEN_BUILD_FLAGS="${EIGEN_BUILD_FLAGS} -DCMAKE_INSTALL_LIBDIR=lib"
    RunWithError cmake $EIGEN_BUILD_FLAGS ..
    RunWithError $MAKE_COMMAND
}

BuildPandoraSDK() {

    if [ -d "$PANDORA_PROJECT_DIR/PandoraSDK" ]; then
        echo "PandoraSDK directory already exists. Skipping!"
        return
    fi

    RunWithError cd $PANDORA_PROJECT_DIR
    RunWithError git clone https://github.com/PandoraPFA/PandoraSDK.git
    RunWithError cd PandoraSDK
    RunWithError git checkout $PANDORA_SDK_VERSION
    RunWithError mkdir build
    RunWithError cd build
    RunWithError cmake $COMMON_CMAKE_FLAGS -DCMAKE_MODULE_PATH=$PANDORA_PROJECT_DIR/PandoraPFA/cmakemodules ..
    RunWithError $MAKE_COMMAND
}

BuildPandoraMonitoring() {

    if [ -d "$PANDORA_PROJECT_DIR/PandoraMonitoring" ]; then
        echo "PandoraMonitoring directory already exists. Skipping!"
        return
    fi

    RunWithError cd $PANDORA_PROJECT_DIR
    RunWithError git clone https://github.com/PandoraPFA/PandoraMonitoring.git
    RunWithError cd PandoraMonitoring
    RunWithError git checkout $PANDORA_MONITORING_VERSION
    RunWithError mkdir build
    RunWithError cd build
    RunWithError cmake $COMMON_CMAKE_FLAGS -DCMAKE_MODULE_PATH="$PANDORA_PROJECT_DIR/PandoraPFA/cmakemodules;$ROOTSYS/etc/cmake" \
        -DPandoraSDK_DIR=$PANDORA_PROJECT_DIR/PandoraSDK/build/install/lib/cmake/PandoraSDK ..
    RunWithError $MAKE_COMMAND
}

BuildLArContent() {

    if [ -d "$PANDORA_PROJECT_DIR/LArContent" ]; then
        echo "LArContent directory already exists. Skipping!"
        return
    fi

    RunWithError cd $PANDORA_PROJECT_DIR
    RunWithError git clone https://github.com/PandoraPFA/LArContent.git
    RunWithError cd LArContent
    RunWithError git checkout $PANDORA_LAR_CONTENT_VERSION
    RunWithError mkdir build
    RunWithError cd build
    RunWithError cmake $COMMON_CMAKE_FLAGS \
        -DPANDORA_MONITORING=ON \
        -DPandoraSDK_DIR=$PANDORA_PROJECT_DIR/PandoraSDK/build/install/lib/cmake/PandoraSDK \
        -DPandoraMonitoring_DIR=$PANDORA_PROJECT_DIR/PandoraMonitoring/build/install/lib/cmake/PandoraMonitoring \
        "$@" ..
    RunWithError $MAKE_COMMAND
}

BuildLArRecoND() {
    # INFO: Assumed to be cloned already...you have this script!
    RunWithError cd $PANDORA_PROJECT_DIR/LArRecoND
    RunWithError mkdir build
    RunWithError cd build
    RunWithError cmake $COMMON_CMAKE_FLAGS \
        -DPANDORA_MONITORING=ON \
        -DPandoraSDK_DIR=$PANDORA_PROJECT_DIR/PandoraSDK/build/install/lib/cmake/PandoraSDK \
        -DPandoraMonitoring_DIR=$PANDORA_PROJECT_DIR/PandoraMonitoring/build/install/lib/cmake/PandoraMonitoring \
        -DLArContent_DIR=$PANDORA_PROJECT_DIR/LArContent/build/install/lib/cmake/LArContent \
        "$@" ..
    RunWithError $MAKE_COMMAND
}

GetMachineLearningData() {

    if [ -d "$PANDORA_PROJECT_DIR/LArMachineLearningData" ]; then
        echo "LArMachineLearningData directory already exists. Skipping!"
        return
    fi

    RunWithError cd $PANDORA_PROJECT_DIR
    RunWithError git clone https://github.com/PandoraPFA/LArMachineLearningData.git
    RunWithError cd LArMachineLearningData
    RunWithError git checkout $PANDORA_LAR_MLDATA_VERSION

    # INFO: Get the training files.
    #       Since the above lines will fail if the files already exist, this should be safe.
    ./download.sh dunend
    ./download.sh dune lbl
    ./download.sh uboone
}
