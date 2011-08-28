PROJECT( COPASIAPI )
CMAKE_MINIMUM_REQUIRED( VERSION 2.6 )
#SET(LIBRARY_OUTPUT_PATH ${EXECUTABLE_OUTPUT_PATH})

# The version number.
#set (COPASI_VERSION_MAJOR 1)
#set (COPASI_VERSION_MINOR 0)

#---------------------
#    OS STUFF
#---------------------
#IF (WIN32)
#   SET(BUILD_OS WIN32)
#ELSE(WIN32)    
#    IF (UNIX)
#        IF (APPLE)
#            SET(BUILD_OS MAC)
#        ELSEIF(APPLE)
#            SET(BUILD_OS LINUX)
#        ENDIF(APPLE)
#    ENDIF(UNIX)
#ENDIF(WIN32)


#-----------------------------------
#    INCLUDE DIRECTORIES
#-----------------------------------

INCLUDE_DIRECTORIES( 
    BEFORE
    ${COPASIAPI_SOURCE_DIR}/lapack/INCLUDE
    ${COPASIAPI_SOURCE_DIR}/muparser/
    ${COPASIAPI_SOURCE_DIR}/c/
    ${QT_INCLUDES}
    ${COPASIAPI_SOURCE_DIR}
)

#--------------------------------------------------------------------------
#  Special compile flags for Qt
#--------------------------------------------------------------------------
#SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DQT_NO_DEBUG -DQT_XML_LIB -DQT_GUI_LIB -DQT_CORE_LIB -DQT_THREAD_SUPPORT")

#FIND_PACKAGE( Qt4 )

FIND_PACKAGE(Qt4 REQUIRED)
ADD_DEFINITIONS(${QT_DEFINITIONS})

# Build only if the version of qt is newer than Qt4.5.
#IF( QT4_FOUND AND QT_VERSION_MINOR GREATER 5)
#  IF( QT_USE_FILE )
#          SET( QT_USE_QTCORE TRUE )
#	  SET( QT_USE_QTGUI FALSE )
#	  SET( QT_USE_QTXML FALSE )
#	  SET( QT_USE_QTOPENGL FALSE )
#          SET( QT_USE_QT3SUPPORT FALSE )
#     #SET( QT_USE_QTSVG TRUE )
     INCLUDE( ${QT_USE_FILE} )
#   ELSE( QT_USE_FILE )
#     SET( QT_LIBRARIES ${QT_QT_LIBRARY} )
#  ENDIF( QT_USE_FILE )
#ENDIF( QT4_FOUND AND QT_VERSION_MINOR GREATER 5)

#  All of these are third-party libraries
#OPTION( BUILD_LAPACK "check this if your system does not have lapack library" ON )
ADD_SUBDIRECTORY( copasi )  #main package
ADD_SUBDIRECTORY( muparser )  #math parser
ADD_SUBDIRECTORY( curl )      #needed by raptor
ADD_SUBDIRECTORY( raptor )    #RDF parser
#IF (BUILD_LAPACK)
	ADD_SUBDIRECTORY( lapack )     #linear algebra
#ENDIF(BUILD_LAPACK)
ADD_SUBDIRECTORY( libSBML )    #systems biology markup language parser
ADD_SUBDIRECTORY( Antimony )   #human-readable test-based reaction networks
ADD_SUBDIRECTORY( libstructural )   #network structural analysis based on stoichiometry matrix
ADD_SUBDIRECTORY( galib )   #genetic algorithms library



#INSTALL( FILES
#         ${COPASIAPI_SOURCE_DIR}/API/TC_structs.h
#         ${COPASIAPI_SOURCE_DIR}/copasi/copasi_api.h
#         ${COPASIAPI_SOURCE_DIR}/lapack/INCLUDE/f2c.h
#         ${COPASIAPI_SOURCE_DIR}/lapack/INCLUDE/clapack.h
#         ${COPASIAPI_SOURCE_DIR}/lapack/INCLUDE/blaswrap.h
#         DESTINATION c
#       )


#--------------------------
#  MAKE libcopasi
#--------------------------



