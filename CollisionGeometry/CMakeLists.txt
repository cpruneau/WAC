################################################################################################
# Project WAC/CollisionGeometry
################################################################################################

#include(${ROOT_USE_FILE})
#include_directories(${CMAKE_SOURCE_DIR} ${ROOT_INCLUDE_DIRS})


ROOT_GENERATE_DICTIONARY(G__CollisionGeometry CollisionGeometry.hpp CollisionGeometryConfiguration.hpp  CollisionGeometryAnalyzer.hpp CollisionGeometryGenerator.hpp CollisionGeometryHistograms.hpp Nucleon.hpp Nucleus.hpp NucleusGenerator.hpp  LINKDEF CollisionGeometryLinkDef.h)

################################################################################################
# Create a shared library with geneated dictionary
################################################################################################
add_compile_options(-Wall -Wextra -pedantic)
add_library(CollisionGeometry SHARED CollisionGeometry.cpp CollisionGeometryConfiguration.cpp CollisionGeometryAnalyzer.cpp CollisionGeometryGenerator.cpp CollisionGeometryHistograms.cpp Nucleon.cpp Nucleus.cpp NucleusGenerator.cpp  G__CollisionGeometry.cxx)
target_link_libraries(CollisionGeometry Base  EGLib PythiaLib EGPythiaLib  ${ROOT_LIBRARIES} ${EXTRA_LIBS} )
target_include_directories(CollisionGeometry  PUBLIC  Base CollisionGeometry   ${EXTRA_INCLUDES})

add_executable(RunCollisionGeometrySimulation RunCollisionGeometrySimulation.cpp)
target_include_directories(RunCollisionGeometrySimulation PUBLIC Base CollisionGeometry ${EXTRA_INCLUDES})
target_link_libraries(RunCollisionGeometrySimulation Base CollisionGeometry  EGLib CollisionGeometryLib   ${ROOT_LIBRARIES} ${EXTRA_LIBS} ${CMAKE_DL_LIBS})

install(FILES  "${CMAKE_CURRENT_BINARY_DIR}/libCollisionGeometry.rootmap" "${CMAKE_CURRENT_BINARY_DIR}/libCollisionGeometry_rdict.pcm" DESTINATION /Users/claudeapruneau/opt/WAC/lib)
install(TARGETS CollisionGeometry  LIBRARY DESTINATION /Users/claudeapruneau/opt/WAC/lib) 
install(TARGETS RunCollisionGeometrySimulation DESTINATION /Users/claudeapruneau/opt/WAC/bin) 

