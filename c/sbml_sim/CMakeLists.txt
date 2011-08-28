SET(LIBRARY_OUTPUT_PATH ${EXECUTABLE_OUTPUT_PATH})

 INCLUDE_DIRECTORIES( BEFORE
    ${COPASIAPI_SOURCE_DIR}/libSBML
    ${COPASIAPI_SOURCE_DIR}/libSBML/include
    ${COPASIAPI_SOURCE_DIR}/libSBML/src
    ${COPASIAPI_SOURCE_DIR}/c/sbml_sim
    ${COPASIAPI_SOURCE_DIR}/muparser
    ${COPASIAPI_SOURCE_DIR}/qwt
    ${COPASIAPI_SOURCE_DIR}/qwtplot3d/include
    ${COPASIAPI_SOURCE_DIR}/optim
    ${COPASIAPI_SOURCE_DIR}/lapack/INCLUDE
    ${COPASIAPI_SOURCE_DIR}/cvode260/include
    ${TINKERCELL_BINARY_DIR}/cvode260/include
    ${COPASIAPI_SOURCE_DIR}/cvode260/include/cvode
    ${TINKERCELL_BINARY_DIR}/cvode260/include/cvode
    ${COPASIAPI_SOURCE_DIR}/cvode260/include/sundials
    ${TINKERCELL_BINARY_DIR}/cvode260/include/sundials
    ${COPASIAPI_SOURCE_DIR}/cvode260/include/nvector
    ${TINKERCELL_BINARY_DIR}/cvode260/include/nvector
 )
 
IF (WIN32)
   ADD_DEFINITIONS(-D_WIN32)
ENDIF (WIN32)

ADD_LIBRARY( sbml_sim 
	SHARED 
	sbml_sim.cpp )

TARGET_LINK_LIBRARIES( sbml_sim 
   sbml_sim 
   sbml 
   ode
   ssa 
   muparser
)
 
#IF ( WIN32 )
#  INSTALL(TARGETS sbml_sim DESTINATION bin/../)
#ELSE ( WIN32 )
#  INSTALL(TARGETS sbml_sim LIBRARY DESTINATION bin/../)
#ENDIF( WIN32 )

