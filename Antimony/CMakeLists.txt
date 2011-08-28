SET(LIBRARY_OUTPUT_PATH ${EXECUTABLE_OUTPUT_PATH})

INCLUDE_DIRECTORIES( BEFORE
	${COPASIAPI_SOURCE_DIR}/libSBML/include/
	${COPASIAPI_SOURCE_DIR}/libSBML/include/sbml/
	${COPASIAPI_SOURCE_DIR}/libSBML/src/
	${COPASIAPI_SOURCE_DIR}/Antimony/src/
)

SET( ANTIMONY_SRC
#    src/antimony2sbml.cpp
    src/antimony_api.cpp
    src/Antimony.cpp
    src/antimony.tab.cpp
    src/dnastrand.cpp
    src/event.cpp
    src/formula.cpp
    src/module.cpp
    src/reactantlist.cpp
    src/reaction.cpp
    src/registry.cpp
#    src/sbml2antimony.cpp
    src/sbmlx.cpp
    src/stringx.cpp
#    src/testantimony.cpp
    src/typex.cpp
    src/userfunction.cpp
    src/variable.cpp
)
  
#ADD_DEFINITIONS(-DNSBML)
ADD_DEFINITIONS(-DNCELLML)
ADD_DEFINITIONS(-DLIB_EXPORTS=1)
ADD_DEFINITIONS(-D__cplusplus)

ADD_LIBRARY ( antimony
   SHARED
#  STATIC
  ${ANTIMONY_SRC}
)


TARGET_LINK_LIBRARIES( antimony
  antimony
  sbml
)

IF ( WIN32 )
  INSTALL(TARGETS antimony DESTINATION bin/../)
ELSE ( WIN32 )
  INSTALL(TARGETS antimony LIBRARY DESTINATION bin/../)
ENDIF( WIN32 )

