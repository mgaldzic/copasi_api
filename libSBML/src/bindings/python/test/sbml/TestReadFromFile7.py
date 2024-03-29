#
# @file    TestReadFromFile7.py
# @brief   Reads test-data/l2v3-all.xml into memory and tests it.
#
# @author  Akiya Jouraku (Python conversion)
# @author  Sarah Keating 
#
# $Id$
# $HeadURL$
#
# ====== WARNING ===== WARNING ===== WARNING ===== WARNING ===== WARNING ======
#
# DO NOT EDIT THIS FILE.
#
# This file was generated automatically by converting the file located at
# src/sbml/test/TestReadFromFile7.cpp
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


class TestReadFromFile7(unittest.TestCase):


  def test_read_l2v3_all(self):
    reader = libsbml.SBMLReader()
    filename = "../../sbml/test/test-data/"
    filename += "l2v3-all.xml"
    d = reader.readSBML(filename)
    if (d == None):
      pass    
    self.assert_( d.getLevel() == 2 )
    self.assert_( d.getVersion() == 3 )
    m = d.getModel()
    self.assert_( m != None )
    self.assert_( m.getId() ==  "l2v3_all" )
    self.assert_( m.getNumCompartments() == 1 )
    c = m.getCompartment(0)
    self.assert_( c != None )
    self.assert_( c.getId() ==  "a" )
    self.assert_( c.getCompartmentType() ==  "hh" )
    self.assert_( c.getSBOTerm() == 236 )
    self.assert_( c.getSBOTermID() ==  "SBO:0000236" )
    self.assert_( c.getSize() == 2.3 )
    self.assert_( m.getNumCompartmentTypes() == 1 )
    ct = m.getCompartmentType(0)
    self.assert_( ct != None )
    self.assert_( ct.getId() ==  "hh" )
    self.assert_( ct.getSBOTerm() == 236 )
    self.assert_( ct.getSBOTermID() ==  "SBO:0000236" )
    self.assert_( m.getNumSpeciesTypes() == 1 )
    st = m.getSpeciesType(0)
    self.assert_( st != None )
    self.assert_( st.getId() ==  "gg" )
    self.assert_( st.getName() ==  "dd" )
    self.assert_( st.getSBOTerm() == 236 )
    self.assert_( st.getSBOTermID() ==  "SBO:0000236" )
    lost = m.getListOfSpeciesTypes()
    st1 = lost.get(0)
    self.assert_( st1 == st )
    st1 = lost.get("gg")
    self.assert_( st1 == st )
    self.assert_( m.getNumConstraints() == 1 )
    con = m.getConstraint(0)
    self.assert_( con != None )
    ast = con.getMath()
    self.assert_((  "lt(x, 3)" == libsbml.formulaToString(ast) ))
    self.assert_( m.getNumEvents() == 1 )
    e = m.getEvent(0)
    self.assert_( e != None )
    self.assert_( e.getId() ==  "e1" )
    self.assert_( e.getSBOTerm() == 231 )
    self.assert_( e.getSBOTermID() ==  "SBO:0000231" )
    self.assertEqual( True, e.isSetDelay() )
    delay = e.getDelay()
    self.assert_( delay != None )
    self.assert_( delay.getSBOTerm() == 64 )
    self.assert_( delay.getSBOTermID() ==  "SBO:0000064" )
    ast = delay.getMath()
    self.assert_((  "p + 3" == libsbml.formulaToString(ast) ))
    self.assertEqual( True, e.isSetTrigger() )
    trigger = e.getTrigger()
    self.assert_( trigger != None )
    self.assert_( trigger.getSBOTerm() == 64 )
    self.assert_( trigger.getSBOTermID() ==  "SBO:0000064" )
    ast = trigger.getMath()
    self.assert_((  "lt(x, 3)" == libsbml.formulaToString(ast) ))
    loe = m.getListOfEvents()
    e1 = loe.get(0)
    self.assert_( e1 == e )
    e1 = loe.get("e1")
    self.assert_( e1 == e )
    self.assert_( e.getNumEventAssignments() == 1 )
    ea = e.getEventAssignment(0)
    self.assert_( ea != None )
    self.assert_( ea.getVariable() ==  "a" )
    self.assert_( ea.getSBOTerm() == 64 )
    self.assert_( ea.getSBOTermID() ==  "SBO:0000064" )
    ast = ea.getMath()
    self.assert_((  "x * p3" == libsbml.formulaToString(ast) ))
    loea = e.getListOfEventAssignments()
    ea1 = loea.get(0)
    self.assert_( ea1 == ea )
    ea1 = loea.get("a")
    self.assert_( ea1 == ea )
    self.assert_( m.getNumFunctionDefinitions() == 1 )
    fd = m.getFunctionDefinition(0)
    self.assert_( fd != None )
    self.assert_( fd.getId() ==  "fd" )
    self.assert_( fd.getSBOTerm() == 64 )
    self.assert_( fd.getSBOTermID() ==  "SBO:0000064" )
    ast = fd.getMath()
    self.assert_((  "lambda(x, pow(x, 3))" == libsbml.formulaToString(ast) ))
    lofd = m.getListOfFunctionDefinitions()
    fd1 = lofd.get(0)
    self.assert_( fd1 == fd )
    fd1 = lofd.get("fd")
    self.assert_( fd1 == fd )
    self.assert_( m.getNumInitialAssignments() == 1 )
    ia = m.getInitialAssignment(0)
    self.assert_( ia != None )
    self.assert_( ia.getSymbol() ==  "p1" )
    ast = ia.getMath()
    self.assert_((  "x * p3" == libsbml.formulaToString(ast) ))
    self.assert_( m.getNumRules() == 3 )
    alg = m.getRule(0)
    self.assert_( alg != None )
    self.assert_( alg.getSBOTerm() == 64 )
    self.assert_( alg.getSBOTermID() ==  "SBO:0000064" )
    ast = alg.getMath()
    self.assert_((  "pow(x, 3)" == libsbml.formulaToString(ast) ))
    ar = m.getRule(1)
    self.assert_( ar != None )
    self.assert_( ar.getVariable() ==  "p2" )
    self.assert_( ar.getSBOTerm() == 64 )
    self.assert_( ar.getSBOTermID() ==  "SBO:0000064" )
    ast = ar.getMath()
    self.assert_((  "x * p3" == libsbml.formulaToString(ast) ))
    rr = m.getRule(2)
    self.assert_( rr != None )
    self.assert_( rr.getVariable() ==  "p3" )
    self.assert_( rr.getSBOTerm() == 64 )
    self.assert_( rr.getSBOTermID() ==  "SBO:0000064" )
    ast = rr.getMath()
    self.assert_((  "p1 / p" == libsbml.formulaToString(ast) ))
    self.assert_( m.getNumSpecies() == 1 )
    s = m.getSpecies(0)
    self.assert_( s != None )
    self.assert_( s.getId() ==  "s" )
    self.assert_( s.getSpeciesType() ==  "gg" )
    self.assert_( s.getCompartment() ==  "a" )
    self.assert_( s.getSBOTerm() == 236 )
    self.assert_( s.getSBOTermID() ==  "SBO:0000236" )
    self.assertEqual( True, s.isSetInitialAmount() )
    self.assertEqual( False, s.isSetInitialConcentration() )
    self.assert_( s.getInitialAmount() == 0 )
    self.assert_( m.getNumReactions() == 1 )
    r = m.getReaction(0)
    self.assert_( r != None )
    self.assert_( r.getId() ==  "r" )
    self.assertEqual( False, r.getReversible() )
    self.assertEqual( True, r.getFast() )
    self.assertEqual( True, r.isSetKineticLaw() )
    kl = r.getKineticLaw()
    self.assert_( kl != None )
    self.assertEqual( True, kl.isSetMath() )
    ast = kl.getMath()
    self.assert_((  "s * k / p" == libsbml.formulaToString(ast) ))
    self.assert_( kl.getNumParameters() == 2 )
    p = kl.getParameter(0)
    self.assert_( p != None )
    self.assert_( p.getId() ==  "k" )
    self.assert_( p.getUnits() ==  "litre" )
    self.assert_( p.getValue() == 9 )
    ud = p.getDerivedUnitDefinition()
    self.assert_( ud.getNumUnits() == 1 )
    self.assert_( ud.getUnit(0).getKind() == libsbml.UNIT_KIND_LITRE )
    self.assert_( ud.getUnit(0).getExponent() == 1 )
    lop = kl.getListOfParameters()
    p1 = lop.get(0)
    self.assert_( p1 == p )
    p1 = lop.get("k")
    self.assert_( p1 == p )
    p = kl.getParameter(1)
    self.assert_( p != None )
    self.assert_( p.getId() ==  "k1" )
    self.assert_( p.getUnits() ==  "ud1" )
    self.assert_( p.getValue() == 9 )
    ud = p.getDerivedUnitDefinition()
    self.assert_( ud.getNumUnits() == 1 )
    self.assert_( ud.getUnit(0).getKind() == libsbml.UNIT_KIND_MOLE )
    self.assert_( ud.getUnit(0).getExponent() == 1 )
    self.assert_( r.getNumReactants() == 1 )
    self.assert_( r.getNumProducts() == 0 )
    self.assert_( r.getNumModifiers() == 0 )
    sr = r.getReactant(0)
    self.assert_( sr != None )
    self.assert_( sr.getSpecies() ==  "s" )
    self.assert_( sr.getSBOTerm() == 11 )
    self.assert_( sr.getSBOTermID() ==  "SBO:0000011" )
    stoich = sr.getStoichiometryMath()
    self.assert_( stoich != None )
    self.assert_( stoich.getSBOTerm() == 64 )
    self.assert_( stoich.getSBOTermID() ==  "SBO:0000064" )
    ast = stoich.getMath()
    self.assert_((  "s * p" == libsbml.formulaToString(ast) ))
    self.assert_( m.getNumUnitDefinitions() == 1 )
    ud = m.getUnitDefinition(0)
    self.assert_( ud != None )
    self.assert_( ud.getId() ==  "ud1" )
    loud = m.getListOfUnitDefinitions()
    ud1 = loud.get(0)
    self.assert_( ud1 == ud )
    ud1 = loud.get("ud1")
    self.assert_( ud1 == ud )
    self.assert_( ud.getNumUnits() == 1 )
    u = ud.getUnit(0)
    self.assert_( u != None )
    self.assert_( u.getKind() == libsbml.UNIT_KIND_MOLE )
    lou = ud.getListOfUnits()
    u1 = lou.get(0)
    self.assert_( u1 == u )
    d = None
    pass  

def suite():
  suite = unittest.TestSuite()
  suite.addTest(unittest.makeSuite(TestReadFromFile7))

  return suite

if __name__ == "__main__":
  if unittest.TextTestRunner(verbosity=1).run(suite()).wasSuccessful() :
    sys.exit(0)
  else:
    sys.exit(1)
