# @file    TestReadFromFile9.rb
# @brief   Reads tests/l3v1-new.xml into memory and tests it.
#
# @author  Akiya Jouraku (Ruby conversion)
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
# src/sbml/test/TestReadFromFile9.cpp
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
require 'test/unit'
require 'libSBML'

class TestReadFromFile9 < Test::Unit::TestCase

  def isnan(x)
    return (x != x)
  end
  @@SBML_INT_MAX = 2147483647
  def test_read_l3v1_new
    reader = LibSBML::SBMLReader.new()
    filename = "../../sbml/test/test-data/"
    filename += "l3v1-new.xml"
    d = reader.readSBML(filename)
    if (d == nil)
    end
    assert( d.getLevel() == 3 )
    assert( d.getVersion() == 1 )
    m = d.getModel()
    assert( m != nil )
    assert_equal true, m.isSetSubstanceUnits()
    assert_equal true, m.isSetTimeUnits()
    assert_equal true, m.isSetVolumeUnits()
    assert_equal true, m.isSetLengthUnits()
    assert_equal true, m.isSetAreaUnits()
    assert_equal true, m.isSetExtentUnits()
    assert_equal true, m.isSetConversionFactor()
    assert( m.getSubstanceUnits() ==  "mole" )
    assert( m.getTimeUnits() ==  "second" )
    assert( m.getVolumeUnits() ==  "litre" )
    assert( m.getLengthUnits() ==  "metre" )
    assert( m.getAreaUnits() ==  "metre" )
    assert( m.getExtentUnits() ==  "mole" )
    assert( m.getConversionFactor() ==  "p" )
    assert( m.getNumUnitDefinitions() == 2 )
    ud = m.getUnitDefinition(0)
    assert( ud.getNumUnits() == 1 )
    u = ud.getUnit(0)
    assert_equal false, u.isSetExponent()
    assert_equal false, u.isSetScale()
    assert_equal false, u.isSetMultiplier()
    assert_equal true, isnan(u.getExponentAsDouble())
    assert( u.getScale() == @@SBML_INT_MAX )
    assert_equal true, isnan(u.getMultiplier())
    ud = m.getUnitDefinition(1)
    assert( ud.getNumUnits() == 3 )
    u = ud.getUnit(0)
    assert_equal true, u.isSetExponent()
    assert_equal true, u.isSetScale()
    assert_equal true, u.isSetMultiplier()
    assert( u.getExponent() == -1 )
    assert( u.getExponentAsDouble() == -1 )
    assert( u.getScale() == 2 )
    assert( u.getMultiplier() == 1.3 )
    u = ud.getUnit(1)
    assert_equal true, u.isSetExponent()
    assert_equal true, u.isSetScale()
    assert_equal true, u.isSetMultiplier()
    assert( u.getExponentAsDouble() == 1.5 )
    assert( u.getScale() == 10 )
    assert( u.getMultiplier() == 0.5 )
    u = ud.getUnit(2)
    assert_equal true, u.isSetExponent()
    assert_equal true, u.isSetScale()
    assert_equal true, u.isSetMultiplier()
    assert( u.getExponent() == 1 )
    assert( u.getScale() == 0 )
    assert( u.getMultiplier() == 1 )
    assert( m.getNumCompartments() == 3 )
    c = m.getCompartment(0)
    assert_equal true, c.isSetSize()
    assert_equal true, c.isSetSpatialDimensions()
    assert_equal true, c.isSetConstant()
    assert( c.getId() ==  "comp" )
    assert( c.getSize() == 1e-14 )
    assert( c.getSpatialDimensions() == 3 )
    assert( c.getSpatialDimensionsAsDouble() == 3 )
    assert( c.getUnits() ==  "litre" )
    assert( c.getConstant() == true )
    c = m.getCompartment(1)
    assert_equal false, c.isSetSize()
    assert_equal false, c.isSetSpatialDimensions()
    assert_equal true, c.isSetConstant()
    assert( c.getId() ==  "comp1" )
    assert_equal true, isnan(c.getSize())
    assert_equal true, isnan(c.getSpatialDimensionsAsDouble())
    assert( c.getConstant() == false )
    c = m.getCompartment(2)
    assert_equal false, c.isSetSize()
    assert_equal true, c.isSetSpatialDimensions()
    assert_equal false, c.isSetConstant()
    assert( c.getId() ==  "comp2" )
    assert( c.getSpatialDimensionsAsDouble() == 4.6 )
    assert( m.getNumSpecies() == 2 )
    s = m.getSpecies(0)
    assert( s.getId() ==  "ES" )
    assert( s.getCompartment() ==  "comp" )
    assert_equal true, s.isSetConversionFactor()
    assert( s.getConversionFactor() ==  "p" )
    assert_equal true, s.isSetBoundaryCondition()
    assert( s.getBoundaryCondition() == false )
    assert_equal true, s.isSetHasOnlySubstanceUnits()
    assert( s.getHasOnlySubstanceUnits() == false )
    assert_equal true, s.isSetSubstanceUnits()
    assert( s.getSubstanceUnits() ==  "mole" )
    assert_equal true, s.isSetConstant()
    assert( s.getConstant() == false )
    assert_equal true, s.isSetInitialAmount()
    assert( s.getInitialAmount() == 0 )
    assert_equal false, s.isSetInitialConcentration()
    assert_equal true, isnan(s.getInitialConcentration())
    s = m.getSpecies(1)
    assert( s.getId() ==  "P" )
    assert( s.getCompartment() ==  "comp" )
    assert_equal false, s.isSetConversionFactor()
    assert( s.getConversionFactor() ==  "" )
    assert_equal false, s.isSetBoundaryCondition()
    assert_equal false, s.isSetHasOnlySubstanceUnits()
    assert_equal false, s.isSetSubstanceUnits()
    assert( s.getSubstanceUnits() ==  "" )
    assert_equal false, s.isSetConstant()
    assert_equal false, s.isSetInitialAmount()
    assert_equal true, isnan(s.getInitialAmount())
    assert_equal false, s.isSetInitialConcentration()
    assert_equal true, isnan(s.getInitialConcentration())
    assert( m.getNumParameters() == 3 )
    p = m.getParameter(0)
    assert( p.getId() ==  "Keq" )
    assert_equal true, p.isSetValue()
    assert( p.getValue() == 2.5 )
    assert_equal true, p.isSetUnits()
    assert( p.getUnits() ==  "dimensionless" )
    assert_equal true, p.isSetConstant()
    assert( p.getConstant() == true )
    p = m.getParameter(1)
    assert( p.getId() ==  "Keq1" )
    assert_equal false, p.isSetValue()
    assert_equal true, isnan(p.getValue())
    assert_equal false, p.isSetUnits()
    assert( p.getUnits() ==  "" )
    assert_equal true, p.isSetConstant()
    assert( p.getConstant() == false )
    p = m.getParameter(2)
    assert( p.getId() ==  "Keq2" )
    assert_equal false, p.isSetValue()
    assert_equal true, isnan(p.getValue())
    assert_equal false, p.isSetUnits()
    assert( p.getUnits() ==  "" )
    assert_equal false, p.isSetConstant()
    assert( m.getNumReactions() == 3 )
    r = m.getReaction(0)
    assert_equal true, r.isSetFast()
    assert( r.getFast() == false )
    assert_equal true, r.isSetReversible()
    assert( r.getReversible() == false )
    assert_equal true, r.isSetCompartment()
    assert( r.getCompartment() ==  "comp" )
    sr = r.getReactant(0)
    assert_equal true, sr.isSetConstant()
    assert( sr.getConstant() == true )
    assert_equal true, sr.isSetStoichiometry()
    assert( sr.getStoichiometry() == 1 )
    sr = r.getProduct(0)
    assert_equal true, sr.isSetConstant()
    assert( sr.getConstant() == false )
    assert_equal false, sr.isSetStoichiometry()
    assert_equal true, isnan(sr.getStoichiometry())
    kl = r.getKineticLaw()
    assert( kl.getNumLocalParameters() == 2 )
    assert( kl.getNumParameters() == 2 )
    p = kl.getParameter(0)
    assert_equal true, p.isSetUnits()
    assert( p.getUnits() ==  "per_second" )
    assert_equal true, p.isSetValue()
    assert( p.getValue() == 0.1 )
    lp = kl.getLocalParameter(0)
    assert_equal true, lp.isSetUnits()
    assert( lp.getUnits() ==  "per_second" )
    assert_equal true, lp.isSetValue()
    assert( lp.getValue() == 0.1 )
    lp = kl.getLocalParameter(1)
    assert_equal false, lp.isSetUnits()
    assert( lp.getUnits() ==  "" )
    assert_equal false, lp.isSetValue()
    assert_equal true, isnan(lp.getValue())
    r = m.getReaction(1)
    assert_equal true, r.isSetFast()
    assert( r.getFast() == true )
    assert_equal true, r.isSetReversible()
    assert( r.getReversible() == true )
    assert_equal false, r.isSetCompartment()
    assert( r.getCompartment() ==  "" )
    sr = r.getReactant(0)
    assert_equal false, sr.isSetConstant()
    r = m.getReaction(2)
    assert_equal false, r.isSetFast()
    assert_equal false, r.isSetReversible()
    assert_equal false, r.isSetCompartment()
    assert( r.getCompartment() ==  "" )
    d = nil
  end

end
