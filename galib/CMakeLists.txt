#---------------------------------------------------------------------------
# GAlib from lancet.mit.edu/ga/
#---------------------------------------------------------------------------

INCLUDE_DIRECTORIES( 
	BEFORE
	${COPASIAPI_SOURCE_DIR}/galib/
	${COPASIAPI_SOURCE_DIR}/galib/ga
)

FILE(GLOB GALIB_FILES ga/*.cpp)

ADD_LIBRARY( galib
  STATIC
  ${GALIB_FILES}
)
