################################################################################################
#
# CMakeLists.txt for WAC package.
#
# Creates  libraries with dictionary and various executables
#
################################################################################################
cmake_minimum_required(VERSION 3.0 FATAL_ERROR)

################################################################################################
# Project WAC
################################################################################################
project(WAC VERSION 1.0 DESCRIPTION "Wayne State Code" LANGUAGES CXX)

set(WAC_ROOT "${WAC_ROOT}")

#option(USE_ECCENTRICITY "Use eccentricity lib" OFF)

################################################################################################
# You need to tell CMake where to find the ROOT installation. This can be done in a number of ways:
#   - ROOT built with classic configure/make use the provided $ROOTSYS/etc/cmake/FindROOT.cmake
#   - ROOT built with CMake. Add in CMAKE_PREFIX_PATH the installation prefix for ROOT
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})

#---Locate the ROOT package and defines a number of variables (e.g. ROOT_INCLUDE_DIRS)
find_package(ROOT REQUIRED COMPONENTS MathCore RIO Hist Tree Net )
find_library(EG EG )
find_library(PYTHIA8_LIB pythia8 PATHS $ENV{PYTHIA8}/lib)
find_library(EGPYTHIA8 EGPythia8)


#---Define useful ROOT functions and macros (e.g. ROOT_GENERATE_DICTIONARY)
include(${ROOT_USE_FILE})
include_directories(${CMAKE_SOURCE_DIR} ${ROOT_INCLUDE_DIRS} "Base" "Cluster" "CollisionGeometry" "Ampt" "Epos" "Pythia" "Eccentricity" "Blastwave" 
"$ENV{ROOTSYS}/include" 
"$ENV{PYTHIA8}/include" 
"$ENV{PYTHIA8}/include/Pythia8")
add_definitions(${ROOT_CXX_FLAGS})


add_library(EG SHARED IMPORTED)
add_library(PYTHIA8 SHARED IMPORTED)
add_library(EGPYTHIA8 SHARED IMPORTED)
add_library(EGLib SHARED IMPORTED)
add_library(PythiaLib SHARED IMPORTED)
add_library(EGPythiaLib SHARED IMPORTED)
set_target_properties(EGLib PROPERTIES IMPORTED_LOCATION ${EG})
set_target_properties(PythiaLib PROPERTIES IMPORTED_LOCATION ${PYTHIA8_LIB})
set_target_properties(EGPythiaLib PROPERTIES IMPORTED_LOCATION ${EGPYTHIA8})


add_subdirectory(Base)
#add_subdirectory(Cluster)
#add_subdirectory(CollisionGeometry)
add_subdirectory(Pythia)
#add_subdirectory(Ampt)
add_subdirectory(Blastwave)
#add_subdirectory(Eccentricity)
#add_subdirectory(Epos)
#add_subdirectory(LambdaAnalysis)
#add_subdirectory(Music)
#add_subdirectory(Urqmd)
#add_subdirectory(Plots)

