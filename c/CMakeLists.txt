SET(LIBRARY_OUTPUT_PATH ${C_LIBRARY_OUTPUT_PATH})
ADD_DEFINITIONS(-DTC_EXPORTS)

OPTION( USE_COPASI "Use COPASI for simulations" ON )

IF ( UNIX )
    ADD_CUSTOM_COMMAND ( 
        TARGET tinkercellapi POST_BUILD
        COMMAND source doSwig.sh
    )
ENDIF ( UNIX )

#***************************************
#    ODE and SSA libraries
#***************************************

ADD_DEFINITIONS(-DNO_FPRINTF_OUTPUT)
ADD_LIBRARY( ode
  STATIC
  eigen.c
  cvodesim.c
  ${COPASIAPI_SOURCE_DIR}/optim/mat.c
  ${COPASIAPI_SOURCE_DIR}/optim/neldermead.c
)

TARGET_LINK_LIBRARIES( ode ode lapack blas f2c sundials_cvode sundials_nvecserial )

SET_TARGET_PROPERTIES(ode PROPERTIES LINKER_LANGUAGE C)

ADD_LIBRARY( ssa
	STATIC
	mtrand.c
	ssa.c
	cells_ssa.c
)

SET_TARGET_PROPERTIES(ssa PROPERTIES LINKER_LANGUAGE C)
IF ( WIN32 )
	  INSTALL(TARGETS ssa ode DESTINATION lib )
ELSE ( WIN32 )
	  INSTALL(TARGETS ssa ode ARCHIVE DESTINATION lib )
ENDIF( WIN32 )

#*******************************************************************************

INCLUDE_DIRECTORIES( BEFORE
    ${COPASIAPI_SOURCE_DIR}/c
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

ADD_SUBDIRECTORY( sbml_sim )

IF (USE_COPASI)	 #use COPASI for the numerical analysis (more features)
	 
	 #all the simulation work is handled by C++ plugins
	INSTALL( FILES 
			${COPASIAPI_SOURCE_DIR}/API/TC_structs.h
			${COPASIAPI_SOURCE_DIR}/copasi/copasi_api.h
			${COPASIAPI_SOURCE_DIR}/lapack/INCLUDE/f2c.h
			${COPASIAPI_SOURCE_DIR}/lapack/INCLUDE/clapack.h
			${COPASIAPI_SOURCE_DIR}/lapack/INCLUDE/blaswrap.h
	  DESTINATION c
	)

ELSE(USE_COPASI)   #use custom programs for the numerical analysis (less features, but faster)
	INCLUDE( CustomSimulationPrograms.cmake )
ENDIF (USE_COPASI)

INCLUDE_DIRECTORIES( BEFORE
	${COPASIAPI_SOURCE_DIR}/lpsolve
)

ADD_LIBRARY( lpsolve
  SHARED
  lpsolve.c
)

TARGET_LINK_LIBRARIES( lpsolve lpsolve lp_solve tinkercellapi )

IF( WIN32 )
  INSTALL(TARGETS lpsolve DESTINATION plugins/c/ )
ELSE( WIN32 )
  INSTALL(TARGETS lpsolve LIBRARY DESTINATION plugins/c/ )
ENDIF( WIN32 )



