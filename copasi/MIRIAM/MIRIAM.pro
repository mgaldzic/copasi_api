# Begin CVS Header
#   $Source: /fs/turing/cvs/copasi_dev/copasi/MIRIAM/MIRIAM.pro,v $
#   $Revision: 1.18 $
#   $Name: Build-33 $
#   $Author: shoops $
#   $Date: 2010/07/16 19:00:07 $
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

LIB = MIRIAM

CONFIG -= qt

DEPENDPATH += ..
INCLUDEPATH += ..

!win32:{
  QMAKE_CFLAGS_WARN_ON += -Wno-unused-parameter 
  QMAKE_CXXFLAGS_WARN_ON += -Wno-unused-parameter 
}

# Input
HEADERS += CBiologicalDescription.h
HEADERS += CCreator.h
HEADERS += CConstants.h
HEADERS += CMIRIAMResource.h
HEADERS += CModelMIRIAMInfo.h
HEADERS += CModified.h
HEADERS += CReference.h
HEADERS += CRaptorInit.h
HEADERS += CRDFGraph.h
HEADERS += CRDFGraphConverter.h
HEADERS += CRDFLiteral.h
HEADERS += CRDFNode.h
HEADERS += CRDFObject.h
HEADERS += CRDFParser.h
HEADERS += CRDFPredicate.h
HEADERS += CRDFSubject.h
HEADERS += CRDFTriplet.h
HEADERS += CRDFUtilities.h
HEADERS += CRDFWriter.h
HEADERS += WebServicesIssues/stdsoap2.h
HEADERS += WebServicesIssues/soapStub.h
HEADERS += WebServicesIssues/soapH.h
HEADERS += WebServicesIssues/soapMiriamWebServicesSoapBindingProxy.h

SOURCES += CBiologicalDescription.cpp
SOURCES += CCreator.cpp
SOURCES += CConstants.cpp
SOURCES += CMIRIAMResource.cpp
SOURCES += CModelMIRIAMInfo.cpp
SOURCES += CModified.cpp
SOURCES += CReference.cpp
SOURCES += CRaptorInit.cpp
SOURCES += CRDFGraph.cpp
SOURCES += CRDFGraphConverter.cpp
SOURCES += CRDFLiteral.cpp
SOURCES += CRDFNode.cpp
SOURCES += CRDFObject.cpp
SOURCES += CRDFParser.cpp
SOURCES += CRDFPredicate.cpp
SOURCES += CRDFSubject.cpp
SOURCES += CRDFTriplet.cpp
SOURCES += CRDFUtilities.cpp
SOURCES += CRDFWriter.cpp
SOURCES += WebServicesIssues/soapC.cpp
SOURCES += WebServicesIssues/soapMiriamWebServicesSoapBindingProxy.cpp
SOURCES += WebServicesIssues/stdsoap2.cpp

DISTFILES += WebServicesIssues/MiriamWebServicesSoapBinding.nsmap

include(../lib.pri)
include(../common.pri)
