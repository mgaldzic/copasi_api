# @file    TestMathReadFromFile2.rb
# @brief   Tests for reading MathML from files into ASTNodes.
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
# src/math/test/TestReadFromFile2.cpp
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

class TestMathReadFromFile2 < Test::Unit::TestCase

  def test_read_MathML_2
    reader = LibSBML::SBMLReader.new()
    filename = "../../math/test/test-data/"
    filename += "mathML_2.xml"
    d = reader.readSBML(filename)
    if (d == nil)
    end
    m = d.getModel()
    assert( m != nil )
    assert( m.getNumFunctionDefinitions() == 2 )
    assert( m.getNumInitialAssignments() == 1 )
    assert( m.getNumRules() == 2 )
    fd = m.getFunctionDefinition(0)
    fd_math = fd.getMath()
    assert( fd_math.getType() == LibSBML::AST_LAMBDA )
    assert( fd_math.getNumChildren() == 1 )
    assert ((  "lambda()" == LibSBML::formulaToString(fd_math) ))
    child = fd_math.getChild(0)
    assert( child.getType() == LibSBML::AST_UNKNOWN )
    assert( child.getNumChildren() == 0 )
    assert ((  "" == LibSBML::formulaToString(child) ))
    fd = m.getFunctionDefinition(1)
    fd1_math = fd.getMath()
    assert( fd1_math.getType() == LibSBML::AST_LAMBDA )
    assert( fd1_math.getNumChildren() == 2 )
    assert ((                            "lambda(x, piecewise(p, leq(x, 4)))" == LibSBML::formulaToString(fd1_math) ))
    child1 = fd1_math.getRightChild()
    assert( child1.getType() == LibSBML::AST_FUNCTION_PIECEWISE )
    assert( child1.getNumChildren() == 2 )
    assert ((                                      "piecewise(p, leq(x, 4))" == LibSBML::formulaToString(child1) ))
    c1 = child1.getChild(0)
    assert( c1.getType() == LibSBML::AST_NAME )
    assert( c1.getNumChildren() == 0 )
    assert ((  "p" == LibSBML::formulaToString(c1) ))
    c2 = child1.getChild(1)
    assert( c2.getType() == LibSBML::AST_RELATIONAL_LEQ )
    assert( c2.getNumChildren() == 2 )
    assert ((  "leq(x, 4)" == LibSBML::formulaToString(c2) ))
    ia = m.getInitialAssignment(0)
    ia_math = ia.getMath()
    assert( ia_math.getType() == LibSBML::AST_FUNCTION_PIECEWISE )
    assert( ia_math.getNumChildren() == 4 )
    assert ((                      "piecewise(-x, lt(x, 0), 0, eq(x, 0))" == LibSBML::formulaToString(ia_math) ))
    child1 = ia_math.getChild(0)
    child2 = ia_math.getChild(1)
    child3 = ia_math.getChild(2)
    child4 = ia_math.getChild(3)
    assert( child1.getType() == LibSBML::AST_MINUS )
    assert( child1.getNumChildren() == 1 )
    assert ((  "-x" == LibSBML::formulaToString(child1) ))
    assert( child2.getType() == LibSBML::AST_RELATIONAL_LT )
    assert( child2.getNumChildren() == 2 )
    assert ((  "lt(x, 0)" == LibSBML::formulaToString(child2) ))
    assert( child3.getType() == LibSBML::AST_REAL )
    assert( child3.getNumChildren() == 0 )
    assert ((  "0" == LibSBML::formulaToString(child3) ))
    assert( child4.getType() == LibSBML::AST_RELATIONAL_EQ )
    assert( child4.getNumChildren() == 2 )
    assert ((  "eq(x, 0)" == LibSBML::formulaToString(child4) ))
    r = m.getRule(0)
    r_math = r.getMath()
    assert( r_math.getType() == LibSBML::AST_CONSTANT_TRUE )
    assert( r_math.getNumChildren() == 0 )
    assert ((  "true" == LibSBML::formulaToString(r_math) ))
    r = m.getRule(1)
    r1_math = r.getMath()
    assert( r1_math.getType() == LibSBML::AST_FUNCTION_LOG )
    assert( r1_math.getNumChildren() == 2 )
    assert ((  "log(3, x)" == LibSBML::formulaToString(r1_math) ))
    child1 = r1_math.getChild(0)
    child2 = r1_math.getChild(1)
    assert( child1.getType() == LibSBML::AST_REAL )
    assert( child1.getNumChildren() == 0 )
    assert ((  "3" == LibSBML::formulaToString(child1) ))
    assert( child2.getType() == LibSBML::AST_NAME )
    assert( child2.getNumChildren() == 0 )
    assert ((  "x" == LibSBML::formulaToString(child2) ))
    d = nil
  end

end
