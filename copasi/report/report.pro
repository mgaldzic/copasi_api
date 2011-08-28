# Begin CVS Header 
#   $Source: /fs/turing/cvs/copasi_dev/copasi/report/report.pro,v $ 
#   $Revision: 1.23 $ 
#   $Name: Build-33 $ 
#   $Author: shoops $ 
#   $Date: 2010/09/22 13:21:10 $ 
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
# $Revision: 1.23 $ $Author: shoops $ $Date: 2010/09/22 13:21:10 $
######################################################################

LIB = report

# Input
HEADERS += ArtisticLicense.h \
           CArrayElementReference.h \
           CCopasiArray.h \
           CCopasiContainer.h \
           CCopasiObject.h \
           CCopasiObjectName.h \
           CCopasiObjectReference.h \
           CCopasiRootContainer.h \
           CCopasiStaticString.h \
           CCopasiTimer.h \
           CKeyFactory.h \
           COutputAssistant.h \
           CRenameHandler.h \
           CReportDefinition.h \
           CReportDefinitionVector.h \
           CReport.h

SOURCES += CArrayElementReference.cpp \
           CCopasiArray.cpp \
           CCopasiContainer.cpp \
           CCopasiObject.cpp \
           CCopasiObjectName.cpp \
           CCopasiRootContainer.cpp \
           CCopasiStaticString.cpp \
           CCopasiTimer.cpp \
           CKeyFactory.cpp \
           COutputAssistant.cpp \
           CReportDefinition.cpp \
           CReportDefinitionVector.cpp \
           CReport.cpp

include(../lib.pri)
include(../common.pri)
