#
# @file    TestRule.py
# @brief   Rule unit tests
#
# @author  Akiya Jouraku (Python conversion)
# @author  Ben Bornstein 
#
# $Id$
# $HeadURL$
#
# ====== WARNING ===== WARNING ===== WARNING ===== WARNING ===== WARNING ======
#
# DO NOT EDIT THIS FILE.
#
# This file was generated automatically by converting the file located at
# src/sbml/test/TestRule.c
# using the conversion program dev/utilities/translateTests/translateTests.pl.
# Any changes made here will be lost the next time the file is regenerated.
#
# -----------------------------------------------------------------------------
# This file is part of libSBML.  Please visit http://sbml.org for more
# information about SBML, and the latest version of libSBML.
#
# Copyright 2005-2010 California Institute of Technology.
# Copyright 2002-2005 California Institute of Technology and
#                     Japan Science and Technology Corporation.
# 
# This library is free software; you can redistribute it and/or modify it
# under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation.  A copy of the license agreement is provided
# in the file named "LICENSE.txt" included with this software distribution
# and also available online as http://sbml.org/software/libsbml/license.html
# -----------------------------------------------------------------------------

import sys
import unittest
import libsbml


class TestRule(unittest.TestCase):

  global R
  R = None

  def setUp(self):
    self.R = libsbml.AlgebraicRule(2,4)
    if (self.R == None):
      pass    
    pass  

  def tearDown(self):
    _dummyList = [ self.R ]; _dummyList[:] = []; del _dummyList
    pass  

  def test_Rule_init(self):
    self.assert_( self.R.getTypeCode() == libsbml.SBML_ALGEBRAIC_RULE )
    self.assert_( self.R.getMetaId() == "" )
    self.assert_( self.R.getNotes() == None )
    self.assert_( self.R.getAnnotation() == None )
    self.assert_( self.R.getFormula() == "" )
    self.assert_( self.R.getMath() == None )
    pass  

  def test_Rule_setFormula(self):
    formula =  "k1*X0";
    self.R.setFormula(formula)
    self.assert_(( formula == self.R.getFormula() ))
    self.assert_( self.R.isSetFormula() == True )
    if (self.R.getFormula() == formula):
      pass    
    self.R.setFormula(self.R.getFormula())
    self.assert_(( formula == self.R.getFormula() ))
    self.R.setFormula( "")
    self.assert_( self.R.isSetFormula() == False )
    if (self.R.getFormula() != None):
      pass    
    pass  

  def test_Rule_setMath(self):
    math = libsbml.parseFormula("1 + 1")
    self.R.setMath(math)
    self.assert_( self.R.getMath() != math )
    self.assertEqual( True, self.R.isSetMath() )
    self.R.setMath(self.R.getMath())
    self.assert_( self.R.getMath() != math )
    self.R.setMath(None)
    self.assertEqual( False, self.R.isSetMath() )
    if (self.R.getMath() != None):
      pass    
    pass  

def suite():
  suite = unittest.TestSuite()
  suite.addTest(unittest.makeSuite(TestRule))

  return suite

if __name__ == "__main__":
  if unittest.TextTestRunner(verbosity=1).run(suite()).wasSuccessful() :
    sys.exit(0)
  else:
    sys.exit(1)
