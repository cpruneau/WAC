################################################################################################
# Project WAC/HadronGas
################################################################################################

ROOT_GENERATE_DICTIONARY(G__HadronGas Hadron.hpp HadronGasConfiguration.hpp  HadronGas.hpp HadronGasHistograms.hpp HadronGasParticleHistograms.hpp  ChemicalPotential.hpp LINKDEF HadronGasLinkDef.h)

################################################################################################
# Create a shared library with geneated dictionary
################################################################################################
add_compile_options(-Wall -Wextra -pedantic )
add_library(HadronGas SHARED Hadron.cpp HadronGasConfiguration.cpp  HadronGas.cpp HadronGasHistograms.cpp HadronGasParticleHistograms.cpp Arsenal.cpp ChemicalPotential.cpp  Table2D.cpp G__HadronGas.cxx)
target_link_libraries(HadronGas  Base  ${ROOT_LIBRARIES} MathMore ${EXTRA_LIBS} )
target_include_directories(HadronGas  PUBLIC Base HadronGas ${EXTRA_INCLUDES})

#add_executable(RunHadronGas RunHadronGas.cpp)
#target_include_directories(RunHadronGas PUBLIC HadronGas  Base ${EXTRA_INCLUDES})
#target_link_libraries(RunHadronGas HadronGas  Base ${ROOT_LIBRARIES} MathMore ${EXTRA_LIBS} ${CMAKE_DL_LIBS})

################################################################################################
# Install
################################################################################################

install(FILES  "${CMAKE_CURRENT_BINARY_DIR}/libHadronGas.rootmap" "${CMAKE_CURRENT_BINARY_DIR}/libHadronGas_rdict.pcm" DESTINATION "$ENV{WAC_LIB}")
install(TARGETS HadronGas  LIBRARY DESTINATION "$ENV{WAC_LIB}")
#install(TARGETS RunHadronGas DESTINATION "$ENV{WAC_BIN}") 

