# Begin CVS Header
#   $Source: /fs/turing/cvs/copasi_dev/copasi/elementaryFluxModes/elementaryFluxModes.pro,v $
#   $Revision: 1.16 $
#   $Name: Build-33 $
#   $Author: heilmand $
#   $Date: 2010/08/02 15:12:41 $
# End CVS Header

# Copyright (C) 2010 by Pedro Mendes, Virginia Tech Intellectual 
# Properties, Inc., University of Heidelberg, and The University 
# of Manchester. 
# All rights reserved. 

# Copyright (C) 2008 by Pedro Mendes, Virginia Tech Intellectual
# Properties, Inc., EML Research, gGmbH, University of Heidelberg,
# and The University of Manchester.
# All rights reserved.

######################################################################
# Created by Stefan Hoops <shoops@vt.edu> Sat May 17 20:42:30 2003
######################################################################

LIB = elementaryFluxModes

# Input
HEADERS += CBitPatternMethod.h
SOURCES += CBitPatternMethod.cpp

HEADERS += CBitPatternTreeMethod.h
SOURCES += CBitPatternTreeMethod.cpp

HEADERS += CBitPatternTree.h
SOURCES += CBitPatternTree.cpp

HEADERS += CBitPatternTreeNode.h
SOURCES += CBitPatternTreeNode.cpp

HEADERS += CEFMAlgorithm.h
SOURCES += CEFMAlgorithm.cpp

HEADERS += CEFMMethod.h
SOURCES += CEFMMethod.cpp

HEADERS += CEFMProblem.h
SOURCES += CEFMProblem.cpp

HEADERS += CEFMTask.h
SOURCES += CEFMTask.cpp

HEADERS += CFluxMode.h
SOURCES += CFluxMode.cpp

HEADERS += CFluxScore.h
SOURCES += CFluxScore.cpp

contains(DEFINES, COPASI_SSA) {
   HEADERS += CSSAMethod.h
   SOURCES += CSSAMethod.cpp
}

HEADERS += CStepMatrix.h
SOURCES += CStepMatrix.cpp

HEADERS += CStepMatrixColumn.h
SOURCES += CStepMatrixColumn.cpp

HEADERS += CTableauLine.h
SOURCES += CTableauLine.cpp

HEADERS += CTableauMatrix.h
SOURCES += CTableauMatrix.cpp

HEADERS += CZeroSet.h
SOURCES += CZeroSet.cpp

include(../lib.pri)
include(../common.pri)
