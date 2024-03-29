# @file    TestUnit_newSetters.rb
# @brief   Unit unit tests for new set function API
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
# src/sbml/test/TestUnit_newSetters.c
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

class TestUnit_newSetters < Test::Unit::TestCase

  def setup
    @@u = LibSBML::Unit.new(1,2)
    if (@@u == nil)
    end
  end

  def teardown
    @@u = nil
  end

  def test_Unit_removeScale
    i = @@u.setScale(2)
    assert( i == LibSBML::LIBSBML_OPERATION_SUCCESS )
    assert( @@u.getScale() == 2 )
    i = LibSBML::Unit.removeScale(@@u)
    assert( i == LibSBML::LIBSBML_OPERATION_SUCCESS )
    assert( @@u.getScale() == 0 )
    assert( @@u.getMultiplier() == 100 )
  end

  def test_Unit_setExponent1
    i = @@u.setExponent(2)
    assert( i == LibSBML::LIBSBML_OPERATION_SUCCESS )
    assert( @@u.getExponent() == 2 )
  end

  def test_Unit_setExponent2
    i = @@u.setExponent(2.0)
    assert( i == LibSBML::LIBSBML_OPERATION_SUCCESS )
    assert( @@u.getExponent() == 2 )
  end

  def test_Unit_setExponent3
    i = @@u.setExponent(2.2)
    assert( i == LibSBML::LIBSBML_INVALID_ATTRIBUTE_VALUE )
    assert( @@u.getExponent() == 1 )
  end

  def test_Unit_setKind1
    i = @@u.setKind(LibSBML::UnitKind_forName("cell"))
    assert( i == LibSBML::LIBSBML_INVALID_ATTRIBUTE_VALUE )
    assert_equal false, @@u.isSetKind()
  end

  def test_Unit_setKind2
    i = @@u.setKind(LibSBML::UnitKind_forName("litre"))
    assert( i == LibSBML::LIBSBML_OPERATION_SUCCESS )
    assert_equal true, @@u.isSetKind()
  end

  def test_Unit_setMultiplier1
    i = @@u.setMultiplier(2)
    assert( i == LibSBML::LIBSBML_UNEXPECTED_ATTRIBUTE )
    assert( @@u.getMultiplier() == 2 )
  end

  def test_Unit_setMultiplier2
    c = LibSBML::Unit.new(2,2)
    i = c.setMultiplier(4)
    assert( i == LibSBML::LIBSBML_OPERATION_SUCCESS )
    assert( c.getMultiplier() == 4 )
    c = nil
  end

  def test_Unit_setOffset1
    i = @@u.setOffset(2.0)
    assert( i == LibSBML::LIBSBML_UNEXPECTED_ATTRIBUTE )
    assert( @@u.getOffset() == 0 )
  end

  def test_Unit_setOffset2
    u1 = LibSBML::Unit.new(2,1)
    i = u1.setOffset(2.0)
    assert( i == LibSBML::LIBSBML_OPERATION_SUCCESS )
    assert( u1.getOffset() == 2 )
  end

  def test_Unit_setScale1
    i = @@u.setScale(2)
    assert( i == LibSBML::LIBSBML_OPERATION_SUCCESS )
    assert( @@u.getScale() == 2 )
  end

end
