# Begin CVS Header
#   $Source: /fs/turing/cvs/copasi_dev/copasi/plotUI/plotUI.pro,v $
#   $Revision: 1.10 $
#   $Name: Build-33 $
#   $Author: shoops $
#   $Date: 2010/09/07 16:33:27 $
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
# $Revision: 1.10 $ $Author: shoops $ $Date: 2010/09/07 16:33:27 $
######################################################################

LIB = plotUI

HEADERS += CopasiPlot.h \
           COutputHandlerPlot.h \
           curve2dwidget.h \
           HistoWidget.h \
           PlotSubwidget.h \
           plotwindow.h \
           scrollbar.h \
           scrollzoomer.h

SOURCES += CopasiPlot.cpp \
           COutputHandlerPlot.cpp \
           curve2dwidget.cpp \
           HistoWidget.cpp \
           PlotSubwidget.cpp \
           plotwindow.cpp \
           scrollbar.cpp \
           scrollzoomer.cpp

FORMS += PlotSubwidget.ui
FORMS +=  HistoWidget.ui
FORMS +=  curve2dwidget.ui

DISTFILES += QwtLicense

include(../lib.pri)

CONFIG += qt

include(../common.pri)


#The following line was inserted by qt3to4
QT +=  opengl 

#The following line is for implementing QtSvg module
QT += svg qt3support
