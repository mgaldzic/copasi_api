SET(LIBRARY_OUTPUT_PATH ${EXECUTABLE_OUTPUT_PATH})

#----------- include folders ----------------------

INCLUDE_DIRECTORIES( 
	BEFORE
	${COPASIAPI_SOURCE_DIR}/expat
	${COPASIAPI_SOURCE_DIR}/expat/lib
	${COPASIAPI_SOURCE_DIR}/libSBML/include
	${COPASIAPI_SOURCE_DIR}/libSBML/src/
	${COPASIAPI_SOURCE_DIR}/libSBML/src/xml 
	${COPASIAPI_SOURCE_DIR}/libSBML/src/util 
	${COPASIAPI_SOURCE_DIR}/libSBML/src/validator 
	${COPASIAPI_SOURCE_DIR}/libSBML/src/validator/constraints/
	${COPASIAPI_SOURCE_DIR}/libSBML/src/units 
	${COPASIAPI_SOURCE_DIR}/libSBML/src/sbml
	${COPASIAPI_SOURCE_DIR}/libSBML/src/sbml/layout 
	${COPASIAPI_SOURCE_DIR}/libSBML/src/math 
	${COPASIAPI_SOURCE_DIR}/libSBML/src/common 
	${COPASIAPI_SOURCE_DIR}/libSBML/src/compress 
	${COPASIAPI_SOURCE_DIR}/libSBML/src/annotation
)

#------ DEFINITIONS --------------

#ADD_DEFINITIONS(-DLIBSBML_STATIC)
ADD_DEFINITIONS(-DLIBSBML_EXPORTS)
ADD_DEFINITIONS(-DLIBLAX_EXPORTS)

ADD_DEFINITIONS(-DPACKAGE_NAME=\"libSBML\")
ADD_DEFINITIONS(-DNDEBUG)
#ADD_DEFINITIONS(-DUSE_LIBXML)
ADD_DEFINITIONS(-DUSE_EXPAT)
ADD_DEFINITIONS(-DUSE_LAYOUT)

IF (WIN32)
	ADD_DEFINITIONS(-DWIN32)
	ADD_DEFINITIONS(-D_WIN32)
	ADD_DEFINITIONS(-D_WINDOWS)
	ADD_DEFINITIONS(-D_USRDLL)
	ADD_DEFINITIONS(-D_EXPORTS)
	#ADD_DEFINITIONS(-DIN_LIBXML)
	IF (MINGW)
		ADD_DEFINITIONS(-DCYGWIN)
	ENDIF(MINGW)
ENDIF (WIN32)

#------ SBML source files --------------

FILE( GLOB SBML_SRC_FILES
	${COPASIAPI_SOURCE_DIR}/libSBML/src/xml/XML*.cpp 
	${COPASIAPI_SOURCE_DIR}/libSBML/src/xml/XML*.c
	${COPASIAPI_SOURCE_DIR}/libSBML/src/xml/Expat*.cpp 
	${COPASIAPI_SOURCE_DIR}/libSBML/src/xml/Expat*.c
	${COPASIAPI_SOURCE_DIR}/libSBML/src/util/*.cpp 
	${COPASIAPI_SOURCE_DIR}/libSBML/src/util/*.c
	${COPASIAPI_SOURCE_DIR}/libSBML/src/validator/*.cpp 
	${COPASIAPI_SOURCE_DIR}/libSBML/src/validator/*.c
	${COPASIAPI_SOURCE_DIR}/libSBML/src/validator/constraints/*.cpp
	${COPASIAPI_SOURCE_DIR}/libSBML/src/validator/constraints/*.c
	${COPASIAPI_SOURCE_DIR}/libSBML/src/units/*.cpp 
	${COPASIAPI_SOURCE_DIR}/libSBML/src/units/*.c
	${COPASIAPI_SOURCE_DIR}/libSBML/src/math/*.cpp
	${COPASIAPI_SOURCE_DIR}/libSBML/src/math/*.c
	${COPASIAPI_SOURCE_DIR}/libSBML/src/annotation/*.cpp 
	${COPASIAPI_SOURCE_DIR}/libSBML/src/annotation/*.c
	${COPASIAPI_SOURCE_DIR}/libSBML/src/sbml/*.cpp
	${COPASIAPI_SOURCE_DIR}/libSBML/src/sbml/*.c
	${COPASIAPI_SOURCE_DIR}/libSBML/src/sbml/layout/*.cpp
	${COPASIAPI_SOURCE_DIR}/libSBML/src/sbml/layout/*.c
	${COPASIAPI_SOURCE_DIR}/libSBML/src/common/*.cpp
	${COPASIAPI_SOURCE_DIR}/libSBML/src/common/*.c
	${COPASIAPI_SOURCE_DIR}/libSBML/src/compress/CompressCommon.cpp
	${COPASIAPI_SOURCE_DIR}/libSBML/src/compress/InputDecompressor.cpp
	${COPASIAPI_SOURCE_DIR}/libSBML/src/compress/OutputCompressor.cpp
)

ADD_DEFINITIONS(-DHAVE_EXPAT_CONFIG_H)
FILE( GLOB  EXPAT_SRC_FILES   ${COPASIAPI_SOURCE_DIR}/expat/lib/*.c)

#------ make the library --------------

ADD_LIBRARY( sbml
  SHARED
  #STATIC
  ${SBML_SRC_FILES}
  ${EXPAT_SRC_FILES}
)

TARGET_LINK_LIBRARIES( sbml
  sbml
)

#------ install the library in the main TinkerCell folder --------------

IF ( WIN32 )
   INSTALL(TARGETS sbml DESTINATION bin/../)
ELSE ( WIN32 )
   INSTALL(TARGETS sbml LIBRARY DESTINATION bin/../)
ENDIF( WIN32 )

