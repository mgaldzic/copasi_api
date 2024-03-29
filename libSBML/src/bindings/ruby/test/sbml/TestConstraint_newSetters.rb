# @file    TestConstraint_newSetters.rb
# @brief   Constraint unit tests for new set function API
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
# src/sbml/test/TestConstraint_newSetters.c
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

class TestConstraint_newSetters < Test::Unit::TestCase

  def setup
    @@c = LibSBML::Constraint.new(2,4)
    if (@@c == nil)
    end
  end

  def teardown
    @@c = nil
  end

  def test_Constraint_setMath1
    math = LibSBML::parseFormula("2 * k")
    i = @@c.setMath(math)
    assert( i == LibSBML::LIBSBML_OPERATION_SUCCESS )
    assert( @@c.getMath() != math )
    assert_equal true, @@c.isSetMath()
    i = @@c.setMath(nil)
    assert( i == LibSBML::LIBSBML_OPERATION_SUCCESS )
    assert( @@c.getMath() == nil )
    assert_equal false, @@c.isSetMath()
    math = nil
  end

  def test_Constraint_setMath2
    math = LibSBML::ASTNode.new(LibSBML::AST_TIMES)
    i = @@c.setMath(math)
    assert( i == LibSBML::LIBSBML_INVALID_OBJECT )
    assert_equal false, @@c.isSetMath()
    math = nil
  end

  def test_Constraint_setMessage1
    node = LibSBML::XMLNode.new()
    i = @@c.setMessage(node)
    assert( i == LibSBML::LIBSBML_INVALID_OBJECT )
    assert( @@c.isSetMessage() == false )
    i = @@c.unsetMessage()
    assert( i == LibSBML::LIBSBML_OPERATION_SUCCESS )
    assert_equal false, @@c.isSetMessage()
    if (@@c.getMessage() != nil)
    end
    node = nil
  end

  def test_Constraint_setMessage2
    text = LibSBML::XMLNode.convertStringToXMLNode(" Some text ",nil)
    triple = LibSBML::XMLTriple.new("p", "http://www.w3.org/1999/xhtml", "")
    att = LibSBML::XMLAttributes.new()
    xmlns = LibSBML::XMLNamespaces.new()
    xmlns.add( "http://www.w3.org/1999/xhtml", "")
    p = LibSBML::XMLNode.new(triple,att,xmlns)
    p.addChild(text)
    triple1 = LibSBML::XMLTriple.new("message", "", "")
    att1 = LibSBML::XMLAttributes.new()
    node = LibSBML::XMLNode.new(triple1,att1)
    node.addChild(p)
    i = @@c.setMessage(node)
    assert( i == LibSBML::LIBSBML_OPERATION_SUCCESS )
    assert( @@c.isSetMessage() == true )
    i = @@c.unsetMessage()
    assert( i == LibSBML::LIBSBML_OPERATION_SUCCESS )
    assert_equal false, @@c.isSetMessage()
    if (@@c.getMessage() != nil)
    end
    node = nil
  end

end
