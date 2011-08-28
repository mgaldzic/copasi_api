SET(LIBRARY_OUTPUT_PATH ${C_LIBRARY_OUTPUT_PATH})

ADD_LIBRARY( runsteadystate 
  SHARED
  runsteadystate.c
)
TARGET_LINK_LIBRARIES( runsteadystate runsteadystate tinkercellapi )

IF ( WIN32 )
  INSTALL(TARGETS runsteadystate DESTINATION plugins/c)
ELSE ( WIN32 )
  INSTALL(TARGETS runsteadystate LIBRARY DESTINATION plugins/c)
ENDIF( WIN32 )

ADD_LIBRARY( runcvode 
  SHARED
  runcvode.c
)

TARGET_LINK_LIBRARIES( runcvode runcvode tinkercellapi )

ADD_LIBRARY( runssa 
  SHARED
  runssa.c
)

TARGET_LINK_LIBRARIES( runssa runssa tinkercellapi )

ADD_LIBRARY( runloops 
  SHARED
  runLoops.c
)

TARGET_LINK_LIBRARIES( runloops runloops tinkercellapi )

ADD_LIBRARY( loadFullBindingKinetics 
  SHARED
  loadFullBindingKinetics.c
)

TARGET_LINK_LIBRARIES( loadFullBindingKinetics loadFullBindingKinetics tinkercellapi )

ADD_LIBRARY( runbistablega 
  SHARED
  runbistablega.c
)

TARGET_LINK_LIBRARIES( runbistablega runbistablega tinkercellapi )

ADD_LIBRARY( runvaluesattime 
  SHARED
  runvaluesattime.c
)

TARGET_LINK_LIBRARIES( runvaluesattime runvaluesattime tinkercellapi )

ADD_LIBRARY( tabasco_like 
  SHARED
  tabasco_like.c
)

TARGET_LINK_LIBRARIES( tabasco_like tabasco_like tinkercellapi )

ADD_LIBRARY( parametercorrelation 
  SHARED
  parametercorrelation.c
)

TARGET_LINK_LIBRARIES( parametercorrelation parametercorrelation tinkercellapi )

##
## CPack
##


IF ( WIN32 )
  INSTALL(TARGETS runssa DESTINATION plugins/c)
  INSTALL(TARGETS runcvode DESTINATION plugins/c)
  INSTALL(TARGETS runloops DESTINATION plugins/c)
  INSTALL(TARGETS loadFullBindingKinetics DESTINATION plugins/c)
#  INSTALL(TARGETS runbistablega DESTINATION plugins/c)
#  INSTALL(TARGETS runvaluesattime DESTINATION plugins/c)
#  INSTALL(TARGETS parametercorrelation DESTINATION plugins/c)
 # INSTALL(TARGETS tabasco_like DESTINATION plugins/c)
  INSTALL(TARGETS lpsolve DESTINATION plugins/c)
ELSE ( WIN32 )
  INSTALL(TARGETS runssa LIBRARY DESTINATION plugins/c)
  INSTALL(TARGETS runcvode LIBRARY DESTINATION plugins/c)
  INSTALL(TARGETS runloops LIBRARY DESTINATION plugins/c)
  INSTALL(TARGETS loadFullBindingKinetics LIBRARY DESTINATION plugins/c)
#  INSTALL(TARGETS runbistablega LIBRARY DESTINATION plugins/c)
#  INSTALL(TARGETS runvaluesattime LIBRARY DESTINATION plugins/c)
#  INSTALL(TARGETS parametercorrelation LIBRARY DESTINATION plugins/c)
 # INSTALL(TARGETS tabasco_like LIBRARY DESTINATION plugins/c)
  INSTALL(TARGETS lpsolve LIBRARY DESTINATION plugins/c)
ENDIF( WIN32 )

FILE( GLOB C_SRC "*.c" "*.h" ${COPASIAPI_SOURCE_DIR}/lapack/INCLUDE/*.h )

INSTALL( FILES ${C_SRC}
  DESTINATION c
  COMPONENT C_SOURCE
)

FILE( GLOB CVODE_HDRS ${COPASIAPI_SOURCE_DIR}/cvode260/include/cvode/*.h )

INSTALL( FILES ${CVODE_HDRS}
  DESTINATION c/cvode
  COMPONENT CVODE
)

FILE( GLOB 
SUNDIALS_HDRS ${COPASIAPI_SOURCE_DIR}/cvode260/include/sundials/*.h 
SUNDIALS_HDRS ${TINKERCELL_BINARY_DIR}/cvode260/include/sundials/*.h 
)

INSTALL( FILES ${SUNDIALS_HDRS}
  DESTINATION c/sundials
  COMPONENT CVODE
)

FILE( GLOB NVECTOR_HDRS ${COPASIAPI_SOURCE_DIR}/cvode260/include/nvector/*.h )

INSTALL( FILES ${NVECTOR_HDRS}
  DESTINATION c/nvector
  COMPONENT CVODE
)

#ADD_SUBDIRECTORY( sbml_sim )

