# Begin CVS Header
#   $Source: /fs/turing/cvs/copasi_dev/copasi/optimization/optimization.pro,v $
#   $Revision: 1.35 $
#   $Name: Build-33 $
#   $Author: shoops $
#   $Date: 2010/07/16 19:01:58 $
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

######################################################################
# $Revision: 1.35 $ $Author: shoops $ $Date: 2010/07/16 19:01:58 $
######################################################################

LIB = optimization

# Input
HEADERS += COptItem.h \
           COptMethod.h \
           COptMethodEP.h \
           COptMethodGA.h \
           COptMethodGASR.h \
           COptMethodHookeJeeves.h \
           COptMethodLevenbergMarquardt.h \
           COptMethodNelderMead.h \
           COptMethodPraxis.h \
           COptMethodPS.h \
           COptMethodSA.h \
           COptMethodSRES.h \
           COptMethodSteepestDescent.h \
           COptMethodStatistics.h \
           COptMethodTruncatedNewton.h \
           COptProblem.h \
           COptTask.h \
           CPraxis.h \
           CRandomSearch.h \
           CTruncatedNewton.h \
           FminBrent.h

SOURCES += COptItem.cpp \
           COptMethod.cpp \
           COptMethodEP.cpp \
           COptMethodGA.cpp \
           COptMethodGASR.cpp \
           COptMethodHookeJeeves.cpp \
           COptMethodLevenbergMarquardt.cpp \
           COptMethodNelderMead.cpp \
           COptMethodPraxis.cpp \
           COptMethodPS.cpp \
           COptMethodSA.cpp \
           COptMethodSRES.cpp \
           COptMethodSteepestDescent.cpp \
           COptMethodStatistics.cpp \
           COptMethodTruncatedNewton.cpp \
           COptProblem.cpp \
           COptTask.cpp \
           CPraxis.cpp \
           CRandomSearch.cpp \
           CTruncatedNewton.cpp \
           FminBrent.cpp

include(../lib.pri)
include(../common.pri)
