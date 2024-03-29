INCLUDE_DIRECTORIES( BEFORE
    ${COPASIAPI_SOURCE_DIR}/curl/include
    ${COPASIAPI_SOURCE_DIR}/expat
	${COPASIAPI_SOURCE_DIR}/expat/lib
	${COPASIAPI_SOURCE_DIR}/raptor/src
	${COPASIAPI_SOURCE_DIR}/raptor/librdfa
	${COPASIAPI_SOURCE_DIR}/raptor/utils
)

#ADD_DEFINITIONS(-D__cplusplus)
ADD_DEFINITIONS(-DHAVE_CONFIG_H)
ADD_DEFINITIONS(-DLIBRDFA_IN_RAPTOR)
ADD_DEFINITIONS(-DRAPTOR_INTERNAL=1)
ADD_DEFINITIONS(-DRAPTOR_V2_EXPERIMENTAL=1)

IF(WIN32)
	ADD_DEFINITIONS(-DWIN32)
	IF (MINGW)
		ADD_DEFINITIONS(-D__MINGW32__)
	ELSE(MINGW)
		ADD_DEFINITIONS(-Dstrcasecmp=_stricmp)
		ADD_DEFINITIONS(-Dstrncasecmp=_strnicmp)
		ADD_DEFINITIONS(-Dstrtok_r=strtok_s)
		ADD_DEFINITIONS(-Dgettimeofday=timeGetTime)
	ENDIF(MINGW)
ENDIF(WIN32)

FILE( GLOB RAPTOR_SRC src/*.c src/*.y librdfa/*.c)# utils/*.c)

ADD_DEFINITIONS(-DHAVE_EXPAT_CONFIG_H)
FILE( GLOB  EXPAT_SRC_FILES   ${COPASIAPI_SOURCE_DIR}/expat/lib/*.c)
#SET_SOURCE_FILES_PROPERTIES(${RAPTOR_SRC}  PROPERTIES LANGUAGE CXX)

ADD_LIBRARY( raptor
  SHARED
#  STATIC
  ${RAPTOR_SRC} 
  ${EXPAT_SRC_FILES}
)

#SET_TARGET_PROPERTIES(raptor PROPERTIES LINKER_LANGUAGE CXX)

TARGET_LINK_LIBRARIES( raptor
   raptor
   libcurl
)

#----------------------------------
#  INSTALL libraptor in . (CPack)
#----------------------------------

IF ( WIN32 )
   INSTALL(TARGETS raptor DESTINATION bin/../)
ELSE ( WIN32 )
   INSTALL(TARGETS raptor LIBRARY DESTINATION bin/../)
ENDIF( WIN32 )

