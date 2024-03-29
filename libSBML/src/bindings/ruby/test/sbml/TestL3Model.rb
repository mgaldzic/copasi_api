# @file    TestL3Model.rb
# @brief   L3 Model unit tests
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
# src/sbml/test/TestL3Model.c
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

class TestL3Model < Test::Unit::TestCase

  def setup
    @@m = LibSBML::Model.new(3,1)
    if (@@m == nil)
    end
  end

  def teardown
    @@m = nil
  end

  def test_L3_Model_NS
    assert( @@m.getNamespaces() != nil )
    assert( @@m.getNamespaces().getLength() == 1 )
    assert ((     "http://www.sbml.org/sbml/level3/version1/core" == @@m.getNamespaces().getURI(0) ))
  end

  def test_L3_Model_areaUnits
    units =  "mole";
    assert_equal false, @@m.isSetAreaUnits()
    @@m.setAreaUnits(units)
    assert (( units == @@m.getAreaUnits() ))
    assert_equal true, @@m.isSetAreaUnits()
    if (@@m.getAreaUnits() == units)
    end
    @@m.unsetAreaUnits()
    assert_equal false, @@m.isSetAreaUnits()
    if (@@m.getAreaUnits() != nil)
    end
  end

  def test_L3_Model_conversionFactor
    units =  "mole";
    assert_equal false, @@m.isSetConversionFactor()
    @@m.setConversionFactor(units)
    assert (( units == @@m.getConversionFactor() ))
    assert_equal true, @@m.isSetConversionFactor()
    if (@@m.getConversionFactor() == units)
    end
    @@m.unsetConversionFactor()
    assert_equal false, @@m.isSetConversionFactor()
    if (@@m.getConversionFactor() != nil)
    end
  end

  def test_L3_Model_create
    assert( @@m.getTypeCode() == LibSBML::SBML_MODEL )
    assert( @@m.getMetaId() == "" )
    assert( @@m.getNotes() == nil )
    assert( @@m.getAnnotation() == nil )
    assert( @@m.getId() == "" )
    assert( @@m.getName() == "" )
    assert( @@m.getSubstanceUnits() == "" )
    assert( @@m.getTimeUnits() == "" )
    assert( @@m.getVolumeUnits() == "" )
    assert( @@m.getAreaUnits() == "" )
    assert( @@m.getLengthUnits() == "" )
    assert( @@m.getConversionFactor() == "" )
    assert_equal false, @@m.isSetId()
    assert_equal false, @@m.isSetName()
    assert_equal false, @@m.isSetSubstanceUnits()
    assert_equal false, @@m.isSetTimeUnits()
    assert_equal false, @@m.isSetVolumeUnits()
    assert_equal false, @@m.isSetAreaUnits()
    assert_equal false, @@m.isSetLengthUnits()
    assert_equal false, @@m.isSetConversionFactor()
  end

  def test_L3_Model_createWithNS
    xmlns = LibSBML::XMLNamespaces.new()
    xmlns.add( "http://www.sbml.org", "testsbml")
    sbmlns = LibSBML::SBMLNamespaces.new(3,1)
    sbmlns.addNamespaces(xmlns)
    m = LibSBML::Model.new(sbmlns)
    assert( m.getTypeCode() == LibSBML::SBML_MODEL )
    assert( m.getMetaId() == "" )
    assert( m.getNotes() == nil )
    assert( m.getAnnotation() == nil )
    assert( m.getLevel() == 3 )
    assert( m.getVersion() == 1 )
    assert( m.getNamespaces() != nil )
    assert( m.getNamespaces().getLength() == 2 )
    assert( m.getId() == "" )
    assert( m.getName() == "" )
    assert( m.getSubstanceUnits() == "" )
    assert( m.getTimeUnits() == "" )
    assert( m.getVolumeUnits() == "" )
    assert( m.getAreaUnits() == "" )
    assert( m.getLengthUnits() == "" )
    assert( m.getConversionFactor() == "" )
    assert_equal false, m.isSetId()
    assert_equal false, m.isSetName()
    assert_equal false, m.isSetSubstanceUnits()
    assert_equal false, m.isSetTimeUnits()
    assert_equal false, m.isSetVolumeUnits()
    assert_equal false, m.isSetAreaUnits()
    assert_equal false, m.isSetLengthUnits()
    assert_equal false, m.isSetConversionFactor()
    m = nil
  end

  def test_L3_Model_extentUnits
    units =  "mole";
    assert_equal false, @@m.isSetExtentUnits()
    @@m.setExtentUnits(units)
    assert (( units == @@m.getExtentUnits() ))
    assert_equal true, @@m.isSetExtentUnits()
    if (@@m.getExtentUnits() == units)
    end
    @@m.unsetExtentUnits()
    assert_equal false, @@m.isSetExtentUnits()
    if (@@m.getExtentUnits() != nil)
    end
  end

  def test_L3_Model_free_NULL
  end

  def test_L3_Model_id
    id =  "mitochondria";
    assert_equal false, @@m.isSetId()
    @@m.setId(id)
    assert (( id == @@m.getId() ))
    assert_equal true, @@m.isSetId()
    if (@@m.getId() == id)
    end
    @@m.unsetId()
    assert_equal false, @@m.isSetId()
    if (@@m.getId() != nil)
    end
  end

  def test_L3_Model_lengthUnits
    units =  "mole";
    assert_equal false, @@m.isSetLengthUnits()
    @@m.setLengthUnits(units)
    assert (( units == @@m.getLengthUnits() ))
    assert_equal true, @@m.isSetLengthUnits()
    if (@@m.getLengthUnits() == units)
    end
    @@m.unsetLengthUnits()
    assert_equal false, @@m.isSetLengthUnits()
    if (@@m.getLengthUnits() != nil)
    end
  end

  def test_L3_Model_name
    name =  "My_Favorite_Factory";
    assert_equal false, @@m.isSetName()
    @@m.setName(name)
    assert (( name == @@m.getName() ))
    assert_equal true, @@m.isSetName()
    if (@@m.getName() == name)
    end
    @@m.unsetName()
    assert_equal false, @@m.isSetName()
    if (@@m.getName() != nil)
    end
  end

  def test_L3_Model_substanceUnits
    units =  "mole";
    assert_equal false, @@m.isSetSubstanceUnits()
    @@m.setSubstanceUnits(units)
    assert (( units == @@m.getSubstanceUnits() ))
    assert_equal true, @@m.isSetSubstanceUnits()
    if (@@m.getSubstanceUnits() == units)
    end
    @@m.unsetSubstanceUnits()
    assert_equal false, @@m.isSetSubstanceUnits()
    if (@@m.getSubstanceUnits() != nil)
    end
  end

  def test_L3_Model_timeUnits
    units =  "mole";
    assert_equal false, @@m.isSetTimeUnits()
    @@m.setTimeUnits(units)
    assert (( units == @@m.getTimeUnits() ))
    assert_equal true, @@m.isSetTimeUnits()
    if (@@m.getTimeUnits() == units)
    end
    @@m.unsetTimeUnits()
    assert_equal false, @@m.isSetTimeUnits()
    if (@@m.getTimeUnits() != nil)
    end
  end

  def test_L3_Model_volumeUnits
    units =  "mole";
    assert_equal false, @@m.isSetVolumeUnits()
    @@m.setVolumeUnits(units)
    assert (( units == @@m.getVolumeUnits() ))
    assert_equal true, @@m.isSetVolumeUnits()
    if (@@m.getVolumeUnits() == units)
    end
    @@m.unsetVolumeUnits()
    assert_equal false, @@m.isSetVolumeUnits()
    if (@@m.getVolumeUnits() != nil)
    end
  end

end
