# Begin CVS Header 
#   $Source: /fs/turing/cvs/copasi_dev/copasi/barChart/barChart.pro,v $ 
#   $Revision: 1.5 $ 
#   $Name: Build-33 $ 
#   $Author: shoops $ 
#   $Date: 2010/07/16 18:54:05 $ 
# End CVS Header 

# Copyright (C) 2010 by Pedro Mendes, Virginia Tech Intellectual 
# Properties, Inc., University of Heidelberg, and The University 
# of Manchester. 
# All rights reserved. 

# Copyright (C) 2008 by Pedro Mendes, Virginia Tech Intellectual 
# Properties, Inc., EML Research, gGmbH, University of Heidelberg, 
# and The University of Manchester. 
# All rights reserved. 

# Copyright (C) 2001 - 2007 by Pedro Mendes, Virginia Tech Intellectual 
# Properties, Inc. and EML Research, gGmbH. 
# All rights reserved. 

LIB = barChart

# Input
HEADERS +=  qwt3dBars.h \
  qwt3dScale.h \
  qwt3dPlot.h \
  qwt3dColor.h \
  qwt3dBaseWidget.h 
  

SOURCES += qwt3dBars.cpp \
  qwt3dScale.cpp \
  qwt3dPlot.cpp \
  qwt3dColor.cpp\
  qwt3dBaseWidget.cpp 

  
include(../lib.pri)

CONFIG += qt 

include(../common.pri)

#The following line was inserted by qt3to4
QT += xml opengl 
