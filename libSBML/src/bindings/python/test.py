#!/usr/bin/env python
##
## @file    test.py
## @brief   AutoRunner for Python test scripts
## @author  Akiya Jouraku
##
## $Id: test.py 10866 2010-01-29 19:52:27Z mhucka $
## $HeadURL: https://sbml.svn.sourceforge.net/svnroot/sbml/trunk/libsbml/src/bindings/python/test.py $
##
##<!---------------------------------------------------------------------------
## This file is part of libSBML.  Please visit http://sbml.org for more
## information about SBML, and the latest version of libSBML.
##
## Copyright 2005-2010 California Institute of Technology.
## Copyright 2002-2005 California Institute of Technology and
##                     Japan Science and Technology Corporation.
##
## This library is free software; you can redistribute it and/or modify it
## under the terms of the GNU Lesser General Public License as published by
## the Free Software Foundation.  A copy of the license agreement is provided
## in the file named "LICENSE.txt" included with this software distribution
## and also available online as http://sbml.org/software/libsbml/license.html
##----------------------------------------------------------------------- -->*/

import os
import sys
import re
import glob
import unittest

test_basedir = 'test'
test_subdirs = ('sbml','xml','math','annotation') 
test_files   = "/Test*.py"

def suite():
  suite = unittest.TestSuite()
  for subdir in test_subdirs :
    for file in glob.glob( test_basedir + '/' + subdir + '/' + test_files ) :
      module_name = re.compile(r"\.py$").sub('',os.path.basename(file)) 
      module = __import__(module_name)
      class_name = getattr(module, module_name)
      suite.addTest(unittest.makeSuite(class_name))
  return suite

if __name__ == "__main__":
  if unittest.TextTestRunner(verbosity=1).run(suite()).wasSuccessful() :
    sys.exit(0)
  else:
    sys.exit(1)
