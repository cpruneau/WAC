################################################################################################
# Project WAC/Base
################################################################################################

#include(${ROOT_USE_FILE})
#include_directories(${CMAKE_SOURCE_DIR} ${ROOT_INCLUDE_DIRS})
#add_definitions(${ROOT_CXX_FLAGS})

ROOT_GENERATE_DICTIONARY(G__Base Timer.hpp CanvasCollection.hpp CanvasConfiguration.hpp TaskConfiguration.hpp EnhancedGraph.hpp Factory.hpp Collection.hpp FunctionCollection.hpp GraphConfiguration.hpp HistogramCollection.hpp  Histograms.hpp RandomGenerators.hpp RapidityGenerator.hpp Task.hpp TaskCollection.hpp Property.hpp MessageLogger.hpp Plotter.hpp StatStudyHistograms.hpp StatStudyPlotter.hpp StatStudyModel.hpp  Event.hpp EventFilter.hpp  EventHistos.hpp EventLoop.hpp GeneralizedGaussianConfiguration.hpp GeneralizedGaussianProfile.hpp GeneratorConfiguration.hpp  TwoPartCorrelationAnalyzer.hpp Particle.hpp ParticleFilter.hpp ParticleFilterAliceV0.hpp ParticleHistos.hpp ParticleAnalyzer.hpp ParticleAnalyzerConfiguration.hpp ParticlePlotter.hpp  ParticlePairCombinedHistos.hpp ParticlePairDerivedHistos.hpp ParticlePairFilter.hpp ParticlePairHistos.hpp ParticlePairAnalyzerConfiguration.hpp TrackAndPairConfiguration.hpp NuDynTask.hpp NuDynHistos.hpp NuDynDerivedHistos.hpp NuDynConfiguration.hpp NuDynPlotter.hpp  CollisionGeometry.hpp CollisionGeometryConfiguration.hpp  CollisionGeometryAnalyzer.hpp CollisionGeometryGenerator.hpp CollisionGeometryHistograms.hpp GeometryMoments.hpp  Nucleon.hpp Nucleus.hpp NucleusGenerator.hpp PTCorrelator.hpp PTHistos.hpp RadialBoostConfiguration.hpp RadialBoostHistos.hpp RadialBoostTask.hpp TotalXsection.hpp SubSampleStatCalculator.hpp GlobalAnalyzer.hpp GlobalHistos.hpp GlobalAnalyzerConfiguration.hpp ParticleFilterAliceV0.hpp
LINKDEF BaseLinkDef.h)

################################################################################################
# Create a shared library with geneated dictionary
################################################################################################
add_compile_options(-Wall -Wextra -pedantic)
add_library(Base SHARED Timer.cpp  CanvasCollection.cpp CanvasConfiguration.cpp TaskConfiguration.cpp EnhancedGraph.cpp  Factory.cpp Collection.cpp FunctionCollection.cpp  GraphConfiguration.cpp HistogramCollection.cpp Histograms.cpp RandomGenerators.cpp RapidityGenerator.cpp Task.cpp TaskCollection.cpp Property.cpp MessageLogger.cpp Plotter.cpp  StatStudyHistograms.cpp StatStudyPlotter.cpp StatStudyModel.cpp  Event.cpp EventFilter.cpp  EventHistos.cpp EventLoop.cpp GeneralizedGaussianConfiguration.cpp GeneralizedGaussianProfile.cpp GeneratorConfiguration.cpp TwoPartCorrelationAnalyzer.cpp Particle.cpp  ParticleFilter.cpp ParticleFilterAliceV0.cpp  ParticleHistos.cpp ParticleAnalyzer.cpp  ParticleAnalyzerConfiguration.cpp  ParticlePlotter.cpp ParticlePairCombinedHistos.cpp ParticlePairDerivedHistos.cpp ParticlePairFilter.cpp ParticlePairHistos.cpp ParticlePairAnalyzerConfiguration.cpp TrackAndPairConfiguration.cpp  NuDynTask.cpp NuDynHistos.cpp NuDynDerivedHistos.cpp NuDynConfiguration.cpp  NuDynPlotter.cpp CollisionGeometry.cpp CollisionGeometryConfiguration.cpp CollisionGeometryAnalyzer.cpp CollisionGeometryGenerator.cpp CollisionGeometryHistograms.cpp GeometryMoments.cpp  Nucleon.cpp Nucleus.cpp NucleusGenerator.cpp PTCorrelator.cpp PTHistos.cpp RadialBoostConfiguration.cpp RadialBoostHistos.cpp RadialBoostTask.cpp TotalXsection.cpp SubSampleStatCalculator.cpp GlobalAnalyzer.cpp GlobalHistos.cpp  GlobalAnalyzerConfiguration.cpp ParticleFilterAliceV0.cpp G__Base.cxx)

target_link_libraries(Base ${ROOT_LIBRARIES} ${EXTRA_LIBS} )
target_include_directories(Base  PUBLIC Base ${EXTRA_INCLUDES} ) 

# optimization for big histogram access within the pair inner loop
# WARNING: the sumw2 structure is not used -> estimate uncertainties using sub-samples
add_definitions(-DOPTIMIZEADDBINCONTENT)

add_executable(RunCollisionGeometrySimulation RunCollisionGeometrySimulation.cpp)
target_include_directories(RunCollisionGeometrySimulation PUBLIC Base ${EXTRA_INCLUDES})
target_link_libraries(RunCollisionGeometrySimulation Base  EGLib   ${ROOT_LIBRARIES} ${EXTRA_LIBS} ${CMAKE_DL_LIBS})


install(FILES  "${CMAKE_CURRENT_BINARY_DIR}/libBase.rootmap" "${CMAKE_CURRENT_BINARY_DIR}/libBase_rdict.pcm" DESTINATION "$ENV{WAC_LIB}")
install(TARGETS Base  LIBRARY DESTINATION "$ENV{WAC_LIB}")

install(TARGETS RunCollisionGeometrySimulation DESTINATION "$ENV{WAC_BIN}") 

