# @file    TestSpeciesReference.rb
# @brief   SpeciesReference unit tests
#
# @author  Akiya Jouraku (Ruby conversion)
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
# src/sbml/test/TestSpeciesReference.c
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

class TestSpeciesReference < Test::Unit::TestCase

  def setup
    @@sr = LibSBML::SpeciesReference.new(2,4)
    if (@@sr == nil)
    end
  end

  def teardown
    @@sr = nil
  end

  def test_SpeciesReference_create
    assert( @@sr.getTypeCode() == LibSBML::SBML_SPECIES_REFERENCE )
    assert( @@sr.getMetaId() == "" )
    assert( @@sr.getNotes() == nil )
    assert( @@sr.getAnnotation() == nil )
    assert( @@sr.getSpecies() == "" )
    assert( @@sr.getStoichiometry() == 1 )
    assert( @@sr.getStoichiometryMath() == nil )
    assert( @@sr.getDenominator() == 1 )
    assert_equal false, @@sr.isSetSpecies()
    assert_equal false, @@sr.isSetStoichiometryMath()
  end

  def test_SpeciesReference_createWithNS
    xmlns = LibSBML::XMLNamespaces.new()
    xmlns.add( "http://www.sbml.org", "testsbml")
    sbmlns = LibSBML::SBMLNamespaces.new(2,1)
    sbmlns.addNamespaces(xmlns)
    object = LibSBML::SpeciesReference.new(sbmlns)
    assert( object.getTypeCode() == LibSBML::SBML_SPECIES_REFERENCE )
    assert( object.getMetaId() == "" )
    assert( object.getNotes() == nil )
    assert( object.getAnnotation() == nil )
    assert( object.getLevel() == 2 )
    assert( object.getVersion() == 1 )
    assert( object.getNamespaces() != nil )
    assert( object.getNamespaces().getLength() == 2 )
    object = nil
  end

  def test_SpeciesReference_free_NULL
  end

  def test_SpeciesReference_setId
    species =  "X0";
    @@sr.setId(species)
    assert (( species == @@sr.getId() ))
    assert_equal true, @@sr.isSetId()
    if (@@sr.getId() == species)
    end
    @@sr.setId(@@sr.getId())
    assert (( species == @@sr.getId() ))
    @@sr.setId("")
    assert_equal false, @@sr.isSetId()
    if (@@sr.getId() != nil)
    end
  end

  def test_SpeciesReference_setSpecies
    species =  "X0";
    @@sr.setSpecies(species)
    assert (( species == @@sr.getSpecies() ))
    assert_equal true, @@sr.isSetSpecies()
    if (@@sr.getSpecies() == species)
    end
    @@sr.setSpecies(@@sr.getSpecies())
    assert (( species == @@sr.getSpecies() ))
    @@sr.setSpecies("")
    assert_equal false, @@sr.isSetSpecies()
    if (@@sr.getSpecies() != nil)
    end
  end

  def test_SpeciesReference_setStoichiometryMath
    math = LibSBML::parseFormula("k3 / k2")
    stoich = LibSBML::StoichiometryMath.new(2,4)
    stoich.setMath(math)
    @@sr.setStoichiometryMath(stoich)
    math1 = @@sr.getStoichiometryMath()
    assert( math1 != nil )
    formula = LibSBML::formulaToString(math1.getMath())
    assert( formula != nil )
    assert ((  "k3 / k2" == formula ))
    assert_equal true, @@sr.isSetStoichiometryMath()
  end

end
