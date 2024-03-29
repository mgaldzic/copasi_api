# @file    TestWriteMathML.rb
# @brief   Write MathML unit tests
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
# src/math/test/TestWriteMathML.cpp
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

class TestWriteMathML < Test::Unit::TestCase

  def MATHML_FOOTER
    return "</math>"
  end

  def MATHML_HEADER
    return "<math xmlns=\"http://www.w3.org/1998/Math/MathML\">\n"
  end

  def MATHML_HEADER_UNITS
    return "<math xmlns=\"http://www.w3.org/1998/Math/MathML\""
  end

  def MATHML_HEADER_UNITS2
    return " xmlns:sbml=\"http://www.sbml.org/sbml/level3/version1/core\">\n"
  end

  def XML_HEADER
    return "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"
  end

  def wrapMathML(s)
    r = XML_HEADER()
    r += MATHML_HEADER()
    r += s
    r += MATHML_FOOTER()
    return r
  end

  def wrapMathMLUnits(s)
    r = XML_HEADER()
    r += MATHML_HEADER_UNITS()
    r += MATHML_HEADER_UNITS2()
    r += s
    r += MATHML_FOOTER()
    return r
  end

  def util_NaN
    z = 0.0
    return 0.0/z
  end

  def util_PosInf
    z = 0.0
    return 1.0/z
  end

  def util_NegInf
    z = 0.0
    return -1.0/z
  end

  def equals(*x)
    case x.size
    when 2
      e, s = x
      return e == s
    when 1
      e, = x
      return e == @@oss.str()
    end
  end

  def setup
    @@n = nil
    @@s = nil
  end

  def teardown
    @@n = nil
    @@s = nil
  end

  def test_MathMLFormatter_ci
    expected = wrapMathML("  <ci> foo </ci>\n")
    @@n = LibSBML::parseFormula("foo")
    @@s = LibSBML::writeMathMLToString(@@n)
    assert_equal true, equals(expected,@@s)
  end

  def test_MathMLFormatter_cn_e_notation_1
    expected = wrapMathML("  <cn type=\"e-notation\"> 0 <sep/> 3 </cn>\n"  
    )
    @@n = LibSBML::parseFormula("0e3")
    @@s = LibSBML::writeMathMLToString(@@n)
    assert_equal true, equals(expected,@@s)
  end

  def test_MathMLFormatter_cn_e_notation_2
    expected = wrapMathML("  <cn type=\"e-notation\"> 2 <sep/> 3 </cn>\n"  
    )
    @@n = LibSBML::parseFormula("2e3")
    @@s = LibSBML::writeMathMLToString(@@n)
    assert_equal true, equals(expected,@@s)
  end

  def test_MathMLFormatter_cn_e_notation_3
    expected = wrapMathML("  <cn type=\"e-notation\"> 1234567.8 <sep/> 3 </cn>\n"  
    )
    @@n = LibSBML::parseFormula("1234567.8e3")
    @@s = LibSBML::writeMathMLToString(@@n)
    assert_equal true, equals(expected,@@s)
  end

  def test_MathMLFormatter_cn_e_notation_4
    expected = wrapMathML("  <cn type=\"e-notation\"> 6.0221367 <sep/> 23 </cn>\n"  
    )
    @@n = LibSBML::parseFormula("6.0221367e+23")
    @@s = LibSBML::writeMathMLToString(@@n)
    assert_equal true, equals(expected,@@s)
  end

  def test_MathMLFormatter_cn_e_notation_5
    expected = wrapMathML("  <cn type=\"e-notation\"> 4 <sep/> -6 </cn>\n"  
    )
    @@n = LibSBML::parseFormula(".000004")
    @@s = LibSBML::writeMathMLToString(@@n)
    assert_equal true, equals(expected,@@s)
  end

  def test_MathMLFormatter_cn_e_notation_6
    expected = wrapMathML("  <cn type=\"e-notation\"> 4 <sep/> -12 </cn>\n"  
    )
    @@n = LibSBML::parseFormula(".000004e-6")
    @@s = LibSBML::writeMathMLToString(@@n)
    assert_equal true, equals(expected,@@s)
  end

  def test_MathMLFormatter_cn_e_notation_7
    expected = wrapMathML("  <cn type=\"e-notation\"> -1 <sep/> -6 </cn>\n"  
    )
    @@n = LibSBML::parseFormula("-1e-6")
    @@s = LibSBML::writeMathMLToString(@@n)
    assert_equal true, equals(expected,@@s)
  end

  def test_MathMLFormatter_cn_integer
    expected = wrapMathML("  <cn type=\"integer\"> 5 </cn>\n")
    @@n = LibSBML::parseFormula("5")
    @@s = LibSBML::writeMathMLToString(@@n)
    assert_equal true, equals(expected,@@s)
  end

  def test_MathMLFormatter_cn_rational
    expected = wrapMathML("  <cn type=\"rational\"> 1 <sep/> 3 </cn>\n"  
    )
    @@n = LibSBML::ASTNode.new()
    @@n.setValue(1,3)
    @@s = LibSBML::writeMathMLToString(@@n)
    assert_equal true, equals(expected,@@s)
  end

  def test_MathMLFormatter_cn_real_1
    expected = wrapMathML("  <cn> 1.2 </cn>\n")
    @@n = LibSBML::parseFormula("1.2")
    @@s = LibSBML::writeMathMLToString(@@n)
    assert_equal true, equals(expected,@@s)
  end

  def test_MathMLFormatter_cn_real_2
    expected = wrapMathML("  <cn> 1234567.8 </cn>\n")
    @@n = LibSBML::parseFormula("1234567.8")
    @@s = LibSBML::writeMathMLToString(@@n)
    assert_equal true, equals(expected,@@s)
  end

  def test_MathMLFormatter_cn_real_3
    expected = wrapMathML("  <cn> -3.14 </cn>\n")
    @@n = LibSBML::parseFormula("-3.14")
    @@s = LibSBML::writeMathMLToString(@@n)
    assert_equal true, equals(expected,@@s)
  end

  def test_MathMLFormatter_cn_real_locale
    expected = wrapMathML("  <cn> 2.72 </cn>\n")
    @@n = LibSBML::parseFormula("2.72")
    @@s = LibSBML::writeMathMLToString(@@n)
    assert_equal true, equals(expected,@@s)
  end

  def test_MathMLFormatter_cn_units
    expected = wrapMathMLUnits("  <cn sbml:units=\"mole\"> 1.2 </cn>\n")
    @@n = LibSBML::parseFormula("1.2")
    @@n.setUnits("mole")
    @@s = LibSBML::writeMathMLToString(@@n)
    assert_equal true, equals(expected,@@s)
  end

  def test_MathMLFormatter_constant_exponentiale
    expected = wrapMathML("  <exponentiale/>\n")
    @@n = LibSBML::ASTNode.new(LibSBML::AST_CONSTANT_E)
    @@s = LibSBML::writeMathMLToString(@@n)
    assert_equal true, equals(expected,@@s)
  end

  def test_MathMLFormatter_constant_false
    expected = wrapMathML("  <false/>\n")
    @@n = LibSBML::ASTNode.new(LibSBML::AST_CONSTANT_FALSE)
    @@s = LibSBML::writeMathMLToString(@@n)
    assert_equal true, equals(expected,@@s)
  end

  def test_MathMLFormatter_constant_infinity
    expected = wrapMathML("  <infinity/>\n")
    @@n = LibSBML::ASTNode.new()
    @@n.setValue( util_PosInf() )
    @@s = LibSBML::writeMathMLToString(@@n)
    assert_equal true, equals(expected,@@s)
  end

  def test_MathMLFormatter_constant_infinity_neg
    expected = wrapMathML("  <apply> <minus/> <infinity/> </apply>\n"  
    )
    @@n = LibSBML::ASTNode.new()
    @@n.setValue(- util_PosInf())
    @@s = LibSBML::writeMathMLToString(@@n)
    assert_equal true, equals(expected,@@s)
  end

  def test_MathMLFormatter_constant_notanumber
    expected = wrapMathML("  <notanumber/>\n")
    @@n = LibSBML::ASTNode.new(LibSBML::AST_REAL)
    @@n.setValue( util_NaN() )
    @@s = LibSBML::writeMathMLToString(@@n)
    assert_equal true, equals(expected,@@s)
  end

  def test_MathMLFormatter_constant_true
    expected = wrapMathML("  <true/>\n")
    @@n = LibSBML::ASTNode.new(LibSBML::AST_CONSTANT_TRUE)
    @@s = LibSBML::writeMathMLToString(@@n)
    assert_equal true, equals(expected,@@s)
  end

  def test_MathMLFormatter_csymbol_avogadro
    expected = wrapMathML("  <csymbol encoding=\"text\" " + "definitionURL=\"http://www.sbml.org/sbml/symbols/avogadro\"> NA </csymbol>\n")
    @@n = LibSBML::ASTNode.new(LibSBML::AST_NAME_AVOGADRO)
    @@n.setName("NA")
    @@s = LibSBML::writeMathMLToString(@@n)
    assert_equal true, equals(expected,@@s)
  end

  def test_MathMLFormatter_csymbol_delay
    expected = wrapMathML("  <apply>\n" + 
    "    <csymbol encoding=\"text\" definitionURL=\"http://www.sbml.org/sbml/" + 
    "symbols/delay\"> my_delay </csymbol>\n" + 
    "    <ci> x </ci>\n" + 
    "    <cn> 0.1 </cn>\n" + 
    "  </apply>\n")
    @@n = LibSBML::parseFormula("delay(x, 0.1)")
    @@n.setName("my_delay")
    @@s = LibSBML::writeMathMLToString(@@n)
    assert_equal true, equals(expected,@@s)
  end

  def test_MathMLFormatter_csymbol_time
    expected = wrapMathML("  <csymbol encoding=\"text\" " + "definitionURL=\"http://www.sbml.org/sbml/symbols/time\"> t </csymbol>\n")
    @@n = LibSBML::ASTNode.new(LibSBML::AST_NAME_TIME)
    @@n.setName("t")
    @@s = LibSBML::writeMathMLToString(@@n)
    assert_equal true, equals(expected,@@s)
  end

  def test_MathMLFormatter_function_1
    expected = wrapMathML("  <apply>\n" + 
    "    <ci> foo </ci>\n" + 
    "    <cn type=\"integer\"> 1 </cn>\n" + 
    "    <cn type=\"integer\"> 2 </cn>\n" + 
    "    <cn type=\"integer\"> 3 </cn>\n" + 
    "  </apply>\n")
    @@n = LibSBML::parseFormula("foo(1, 2, 3)")
    @@s = LibSBML::writeMathMLToString(@@n)
    assert_equal true, equals(expected,@@s)
  end

  def test_MathMLFormatter_function_2
    expected = wrapMathML("  <apply>\n" + 
    "    <ci> foo </ci>\n" + 
    "    <cn type=\"integer\"> 1 </cn>\n" + 
    "    <cn type=\"integer\"> 2 </cn>\n" + 
    "    <apply>\n" + 
    "      <ci> bar </ci>\n" + 
    "      <ci> z </ci>\n" + 
    "    </apply>\n" + 
    "  </apply>\n")
    @@n = LibSBML::parseFormula("foo(1, 2, bar(z))")
    @@s = LibSBML::writeMathMLToString(@@n)
    assert_equal true, equals(expected,@@s)
  end

  def test_MathMLFormatter_lambda
    expected = wrapMathML("  <lambda>\n" + 
    "    <bvar>\n" + 
    "      <ci> x </ci>\n" + 
    "    </bvar>\n" + 
    "    <bvar>\n" + 
    "      <ci> y </ci>\n" + 
    "    </bvar>\n" + 
    "    <apply>\n" + 
    "      <root/>\n" + 
    "      <degree>\n" + 
    "        <cn type=\"integer\"> 2 </cn>\n" + 
    "      </degree>\n" + 
    "      <apply>\n" + 
    "        <plus/>\n" + 
    "        <apply>\n" + 
    "          <power/>\n" + 
    "          <ci> x </ci>\n" + 
    "          <cn type=\"integer\"> 2 </cn>\n" + 
    "        </apply>\n" + 
    "        <apply>\n" + 
    "          <power/>\n" + 
    "          <ci> y </ci>\n" + 
    "          <cn type=\"integer\"> 2 </cn>\n" + 
    "        </apply>\n" + 
    "      </apply>\n" + 
    "    </apply>\n" + 
    "  </lambda>\n")
    @@n = LibSBML::parseFormula("lambda(x, y, root(2, x^2 + y^2))")
    @@s = LibSBML::writeMathMLToString(@@n)
    assert_equal true, equals(expected,@@s)
  end

  def test_MathMLFormatter_lambda_no_bvars
    expected = wrapMathML("  <lambda>\n" + 
    "    <apply>\n" + 
    "      <plus/>\n" + 
    "      <cn type=\"integer\"> 2 </cn>\n" + 
    "      <cn type=\"integer\"> 2 </cn>\n" + 
    "    </apply>\n" + 
    "  </lambda>\n")
    @@n = LibSBML::parseFormula("lambda(2 + 2)")
    @@s = LibSBML::writeMathMLToString(@@n)
    assert_equal true, equals(expected,@@s)
  end

  def test_MathMLFormatter_log
    expected = wrapMathML("  <apply>\n" + 
    "    <log/>\n" + 
    "    <logbase>\n" + 
    "      <cn type=\"integer\"> 2 </cn>\n" + 
    "    </logbase>\n" + 
    "    <ci> N </ci>\n" + 
    "  </apply>\n")
    @@n = LibSBML::parseFormula("log(2, N)")
    @@s = LibSBML::writeMathMLToString(@@n)
    assert_equal true, equals(expected,@@s)
  end

  def test_MathMLFormatter_minus
    expected = wrapMathML("  <apply>\n" + 
    "    <minus/>\n" + 
    "    <cn type=\"integer\"> 1 </cn>\n" + 
    "    <cn type=\"integer\"> 2 </cn>\n" + 
    "  </apply>\n")
    @@n = LibSBML::parseFormula("1 - 2")
    @@s = LibSBML::writeMathMLToString(@@n)
    assert_equal true, equals(expected,@@s)
  end

  def test_MathMLFormatter_minus_unary_1
    expected = wrapMathML("  <cn type=\"integer\"> -2 </cn>\n"  
    )
    @@n = LibSBML::parseFormula("-2")
    @@s = LibSBML::writeMathMLToString(@@n)
    assert_equal true, equals(expected,@@s)
  end

  def test_MathMLFormatter_minus_unary_2
    expected = wrapMathML("  <apply>\n" + 
    "    <minus/>\n" + 
    "    <ci> a </ci>\n" + 
    "  </apply>\n")
    @@n = LibSBML::parseFormula("-a")
    @@s = LibSBML::writeMathMLToString(@@n)
    assert_equal true, equals(expected,@@s)
  end

  def test_MathMLFormatter_piecewise
    expected = wrapMathML("  <piecewise>\n" + 
    "    <piece>\n" + 
    "      <apply>\n" + 
    "        <minus/>\n" + 
    "        <ci> x </ci>\n" + 
    "      </apply>\n" + 
    "      <apply>\n" + 
    "        <lt/>\n" + 
    "        <ci> x </ci>\n" + 
    "        <cn type=\"integer\"> 0 </cn>\n" + 
    "      </apply>\n" + 
    "    </piece>\n" + 
    "    <piece>\n" + 
    "      <cn type=\"integer\"> 0 </cn>\n" + 
    "      <apply>\n" + 
    "        <eq/>\n" + 
    "        <ci> x </ci>\n"  + 
    "        <cn type=\"integer\"> 0 </cn>\n" + 
    "      </apply>\n" + 
    "    </piece>\n" + 
    "    <piece>\n" + 
    "      <ci> x </ci>\n" + 
    "      <apply>\n" + 
    "        <gt/>\n" + 
    "        <ci> x </ci>\n"  + 
    "        <cn type=\"integer\"> 0 </cn>\n" + 
    "      </apply>\n" + 
    "    </piece>\n" + 
    "  </piecewise>\n")
    f =  "piecewise(-x, lt(x, 0), 0, eq(x, 0), x, gt(x, 0))";
    @@n = LibSBML::parseFormula(f)
    @@s = LibSBML::writeMathMLToString(@@n)
    assert_equal true, equals(expected,@@s)
  end

  def test_MathMLFormatter_piecewise_otherwise
    expected = wrapMathML("  <piecewise>\n" + 
    "    <piece>\n" + 
    "      <cn type=\"integer\"> 0 </cn>\n" + 
    "      <apply>\n" + 
    "        <lt/>\n" + 
    "        <ci> x </ci>\n" + 
    "        <cn type=\"integer\"> 0 </cn>\n" + 
    "      </apply>\n" + 
    "    </piece>\n" + 
    "    <otherwise>\n" + 
    "      <ci> x </ci>\n"  + 
    "    </otherwise>\n" + 
    "  </piecewise>\n")
    @@n = LibSBML::parseFormula("piecewise(0, lt(x, 0), x)")
    @@s = LibSBML::writeMathMLToString(@@n)
    assert_equal true, equals(expected,@@s)
  end

  def test_MathMLFormatter_plus_binary
    expected = wrapMathML("  <apply>\n" + 
    "    <plus/>\n" + 
    "    <cn type=\"integer\"> 1 </cn>\n" + 
    "    <cn type=\"integer\"> 2 </cn>\n" + 
    "  </apply>\n")
    @@n = LibSBML::parseFormula("1 + 2")
    @@s = LibSBML::writeMathMLToString(@@n)
    assert_equal true, equals(expected,@@s)
  end

  def test_MathMLFormatter_plus_nary_1
    expected = wrapMathML("  <apply>\n" + 
    "    <plus/>\n" + 
    "    <cn type=\"integer\"> 1 </cn>\n" + 
    "    <cn type=\"integer\"> 2 </cn>\n" + 
    "    <cn type=\"integer\"> 3 </cn>\n" + 
    "  </apply>\n")
    @@n = LibSBML::parseFormula("1 + 2 + 3")
    @@s = LibSBML::writeMathMLToString(@@n)
    assert_equal true, equals(expected,@@s)
  end

  def test_MathMLFormatter_plus_nary_2
    expected = wrapMathML("  <apply>\n" + 
    "    <plus/>\n" + 
    "    <cn type=\"integer\"> 1 </cn>\n" + 
    "    <cn type=\"integer\"> 2 </cn>\n" + 
    "    <cn type=\"integer\"> 3 </cn>\n" + 
    "  </apply>\n")
    @@n = LibSBML::parseFormula("(1 + 2) + 3")
    @@s = LibSBML::writeMathMLToString(@@n)
    assert_equal true, equals(expected,@@s)
  end

  def test_MathMLFormatter_plus_nary_3
    expected = wrapMathML("  <apply>\n" + 
    "    <plus/>\n" + 
    "    <cn type=\"integer\"> 1 </cn>\n" + 
    "    <cn type=\"integer\"> 2 </cn>\n" + 
    "    <cn type=\"integer\"> 3 </cn>\n" + 
    "  </apply>\n")
    @@n = LibSBML::parseFormula("1 + (2 + 3)")
    @@s = LibSBML::writeMathMLToString(@@n)
    assert_equal true, equals(expected,@@s)
  end

  def test_MathMLFormatter_plus_nary_4
    expected = wrapMathML("  <apply>\n" + 
    "    <plus/>\n" + 
    "    <cn type=\"integer\"> 1 </cn>\n" + 
    "    <cn type=\"integer\"> 2 </cn>\n" + 
    "    <apply>\n" + 
    "      <times/>\n" + 
    "      <ci> x </ci>\n" + 
    "      <ci> y </ci>\n" + 
    "      <ci> z </ci>\n" + 
    "    </apply>\n" + 
    "    <cn type=\"integer\"> 3 </cn>\n" + 
    "  </apply>\n")
    @@n = LibSBML::parseFormula("1 + 2 + x * y * z + 3")
    @@s = LibSBML::writeMathMLToString(@@n)
    assert_equal true, equals(expected,@@s)
  end

  def test_MathMLFormatter_root
    expected = wrapMathML("  <apply>\n" + 
    "    <root/>\n" + 
    "    <degree>\n" + 
    "      <cn type=\"integer\"> 3 </cn>\n" + 
    "    </degree>\n" + 
    "    <ci> x </ci>\n" + 
    "  </apply>\n")
    @@n = LibSBML::parseFormula("root(3, x)")
    @@s = LibSBML::writeMathMLToString(@@n)
    assert_equal true, equals(expected,@@s)
  end

  def test_MathMLFormatter_sin
    expected = wrapMathML("  <apply>\n" + 
    "    <sin/>\n" + 
    "    <ci> x </ci>\n" + 
    "  </apply>\n")
    @@n = LibSBML::parseFormula("sin(x)")
    @@s = LibSBML::writeMathMLToString(@@n)
    assert_equal true, equals(expected,@@s)
  end

end
