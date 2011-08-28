# Begin CVS Header
#   $Source: /fs/turing/cvs/copasi_dev/copasi/compareExpressions/compareExpressions.pro,v $
#   $Revision: 1.12 $
#   $Name: Build-33 $
#   $Author: shoops $
#   $Date: 2010/07/16 18:57:32 $
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
# $Revision: 1.12 $ $Author: shoops $ $Date: 2010/07/16 18:57:32 $
######################################################################

LIB = compareExpressions

#QMAKE_EXTRA_UNIX_TARGETS += check
#check.commands = cd unittests; $(QMAKE) $$(0);make

# Input
HEADERS += CNormalItem.h \
           CNormalItemPower.h\
           CNormalProduct.h \
           CNormalSum.h \
           CNormalLcm.h \
           CNormalFraction.h \
           CNormalBase.h \
           CNormalChoiceLogical.h \
           CNormalChoice.h \
           CNormalLogical.h \
           CNormalLogicalItem.h \
           CNormalFunction.h \
           CNormalCall.h \
           CNormalGeneralPower.h \
           CNormalTranslation.h\
           ConvertToCEvaluationNode.h \
           CEvaluationNodeNormalizer.h \
           compare_utilities.h

SOURCES += CNormalItem.cpp \
           CNormalItemPower.cpp \
           CNormalProduct.cpp \
           CNormalSum.cpp \
           CNormalLcm.cpp \
           CNormalFraction.cpp \
           CNormalChoice.cpp \
           CNormalChoiceLogical.cpp \
           CNormalLogical.cpp \
           CNormalLogicalItem.cpp \
      	   CNormalFunction.cpp \
      	   CNormalCall.cpp \
           CNormalGeneralPower.cpp \
           CNormalTranslation.cpp \
           CNormalBase.cpp \
           ConvertToCEvaluationNode.cpp \
           CEvaluationNodeNormalizer.cpp \
           compare_utilities.cpp


include(../lib.pri)
include(../common.pri)
