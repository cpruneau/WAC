################################################################################################
# Project WAC/Base
################################################################################################

#include(${ROOT_USE_FILE})
#include_directories(${CMAKE_SOURCE_DIR} ${ROOT_INCLUDE_DIRS})
#add_definitions(${ROOT_CXX_FLAGS})

ROOT_GENERATE_DICTIONARY(G__Base CanvasCollection.hpp CanvasConfiguration.hpp TaskConfiguration.hpp EnhancedGraph.hpp Factory.hpp FunctionCollection.hpp GraphConfiguration.hpp HistogramCollection.hpp  Histograms.hpp RandomGenerators.hpp RapidityGenerator.hpp Task.hpp TaskCollection.hpp Property.hpp MessageLogger.hpp Plotter.hpp StatStudyHistograms.hpp StatStudyPlotter.hpp StatStudyModel.hpp AnalysisConfiguration.hpp Event.hpp EventFilter.hpp  EventHistos.hpp EventLoop.hpp GeneralizedGaussianConfiguration.hpp GeneralizedGaussianProfile.hpp GeneratorConfiguration.hpp  TwoPartCorrelationAnalyzer.hpp Particle.hpp ParticleFilter.hpp ParticleHistos.hpp ParticlePairCombinedHistos.hpp ParticlePairDerivedHistos.hpp ParticlePairFilter.hpp ParticlePairHistos.hpp TrackAndPairConfiguration.hpp NuDynTask.hpp NuDynHistos.hpp NuDynDerivedHistos.hpp PTCorrelator.hpp PTHistos.hpp LINKDEF BaseLinkDef.h)

################################################################################################
# Create a shared library with geneated dictionary
################################################################################################
add_compile_options(-Wall -Wextra -pedantic)
add_library(Base SHARED CanvasCollection.cpp CanvasConfiguration.cpp TaskConfiguration.cpp EnhancedGraph.cpp  Factory.cpp FunctionCollection.cpp  GraphConfiguration.cpp HistogramCollection.cpp Histograms.cpp RandomGenerators.cpp RapidityGenerator.cpp Task.cpp TaskCollection.cpp Property.cpp MessageLogger.cpp Plotter.cpp  StatStudyHistograms.cpp StatStudyPlotter.cpp StatStudyModel.cpp AnalysisConfiguration.cpp Event.cpp EventFilter.cpp  EventHistos.cpp EventLoop.cpp GeneralizedGaussianConfiguration.cpp GeneralizedGaussianProfile.cpp GeneratorConfiguration.cpp TwoPartCorrelationAnalyzer.cpp Particle.cpp  ParticleFilter.cpp ParticleHistos.cpp ParticlePairCombinedHistos.cpp ParticlePairDerivedHistos.cpp ParticlePairFilter.cpp ParticlePairHistos.cpp TrackAndPairConfiguration.cpp  NuDynTask.cpp NuDynHistos.cpp NuDynDerivedHistos.cpp PTCorrelator.cpp PTHistos.cpp G__Base.cxx)
target_link_libraries(Base ${ROOT_LIBRARIES} ${EXTRA_LIBS} )
target_include_directories(Base  PUBLIC Base ${EXTRA_INCLUDES} ) 

install(FILES  "${CMAKE_CURRENT_BINARY_DIR}/libBase.rootmap" "${CMAKE_CURRENT_BINARY_DIR}/libBase_rdict.pcm" DESTINATION "$ENV{WAC_LIB}")
install(TARGETS Base  LIBRARY DESTINATION "$ENV{WAC_LIB}")

