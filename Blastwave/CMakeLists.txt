################################################################################################
# Project WAC/Blastwave
################################################################################################


ROOT_GENERATE_DICTIONARY(G__Blastwave  BWModelConfiguration.hpp BWModelEventGenerator.hpp  LINKDEF BWLinkDef.h)

################################################################################################
# Create a shared library with geneated dictionary
################################################################################################
add_library(Blastwave SHARED  BWModelConfiguration.cpp BWModelEventGenerator.cpp G__Blastwave.cxx)
target_link_libraries(Blastwave Base  ${ROOT_LIBRARIES} ${EXTRA_LIBS} )
target_include_directories(Blastwave  PUBLIC Base Blastwave ${EXTRA_INCLUDES} ) 

# RunBWModelSimulation app
add_executable(RunBWModelSimulation RunBWModelSimulation.cpp)
target_include_directories(RunBWModelSimulation PUBLIC Base Blastwave ${EXTRA_INCLUDES})
target_link_libraries(RunBWModelSimulation Base Blastwave ${CMAKE_DL_LIBS})

install(FILES  "${CMAKE_CURRENT_BINARY_DIR}/libBlastwave.rootmap" "${CMAKE_CURRENT_BINARY_DIR}/libBlastwave_rdict.pcm" DESTINATION "$ENV{WAC_LIB}")
install(TARGETS Blastwave  LIBRARY DESTINATION "$ENV{WAC_LIB}") 
install(TARGETS RunBWModelSimulation DESTINATION "$ENV{WAC_BIN}") 

