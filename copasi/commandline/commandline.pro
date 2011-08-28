# Begin CVS Header
#   $Source: /fs/turing/cvs/copasi_dev/copasi/commandline/commandline.pro,v $
#   $Revision: 1.16 $
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

######################################################################
# $Revision: 1.16 $ $Author: shoops $ $Date: 2010/07/16 18:57:32 $
######################################################################

LIB = commandline

#Input
HEADERS += CConfigurationFile.h \
           COptionParser.h \
           COptions.h

SOURCES += CConfigurationFile.cpp \
           COptionParser.cpp \
           COptions.cpp

include(../lib.pri)
include(../common.pri)

contains(BUILD_PARSER, yes) {
  clo.target = COptionParser.cpp
  clo.depends = COptionParser.xml
  win32:{
    clo.commands = C:\cygwin\bin\bash ../../admin/clo++.sh $$clo.depends
    QMAKE_EXTRA_WIN_TARGETS += clo
  } else {
    clo.commands = ../../admin/clo++.sh $$clo.depends
    QMAKE_EXTRA_UNIX_TARGETS += clo
  }
}


DISTFILES += \
             COptionParser.xml
