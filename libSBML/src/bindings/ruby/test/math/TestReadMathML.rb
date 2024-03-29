# @file    TestReadMathML.rb
# @brief   Read MathML unit tests
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
# src/math/test/TestReadMathML.cpp
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

class TestReadMathML < Test::Unit::TestCase

  def MATHML_FOOTER
    return "</math>"
  end

  def MATHML_HEADER
    return "<math xmlns='http://www.w3.org/1998/Math/MathML'>\n"
  end

  def MATHML_HEADER_UNITS
    return "<math xmlns='http://www.w3.org/1998/Math/MathML'\n"
  end

  def MATHML_HEADER_UNITS2
    return " xmlns:sbml='http://www.sbml.org/sbml/level3/version1/core'>\n"
  end

  def XML_HEADER
    return "<?xml version='1.0' encoding='UTF-8'?>\n"
  end

  def isnan(x)
    return (x != x)
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

  def wrapXML(s)
    r = XML_HEADER()
    r += s
    return r
  end

  def util_isInf(*x)
    e, = x 
    return ( e == util_PosInf() || e == util_NegInf() )
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
    @@f = nil
  end

  def teardown
    @@n = nil
    @@f = nil
  end

  def test_element_abs
    s = wrapMathML("<apply><abs/><ci> x </ci></apply>")
    @@n = LibSBML::readMathMLFromString(s)
    assert( @@n != nil )
    @@f = LibSBML::formulaToString(@@n)
    assert ((  "abs(x)" == @@f ))
  end

  def test_element_and
    s = wrapMathML("<apply> <and/> <ci>a</ci> <ci>b</ci> <ci>c</ci> </apply>"  
    )
    @@n = LibSBML::readMathMLFromString(s)
    assert( @@n != nil )
    @@f = LibSBML::formulaToString(@@n)
    assert ((  "and(a, b, c)" == @@f ))
  end

  def test_element_arccos
    s = wrapMathML("<apply><arccos/><ci> x </ci></apply>")
    @@n = LibSBML::readMathMLFromString(s)
    assert( @@n != nil )
    @@f = LibSBML::formulaToString(@@n)
    assert ((  "acos(x)" == @@f ))
  end

  def test_element_arccosh
    s = wrapMathML("<apply><arccosh/><ci> x </ci></apply>")
    @@n = LibSBML::readMathMLFromString(s)
    assert( @@n != nil )
    @@f = LibSBML::formulaToString(@@n)
    assert ((  "arccosh(x)" == @@f ))
  end

  def test_element_arccot
    s = wrapMathML("<apply><arccot/><ci> x </ci></apply>")
    @@n = LibSBML::readMathMLFromString(s)
    assert( @@n != nil )
    @@f = LibSBML::formulaToString(@@n)
    assert ((  "arccot(x)" == @@f ))
  end

  def test_element_arccoth
    s = wrapMathML("<apply><arccoth/><ci> x </ci></apply>")
    @@n = LibSBML::readMathMLFromString(s)
    assert( @@n != nil )
    @@f = LibSBML::formulaToString(@@n)
    assert ((  "arccoth(x)" == @@f ))
  end

  def test_element_arccsc
    s = wrapMathML("<apply><arccsc/><ci> x </ci></apply>")
    @@n = LibSBML::readMathMLFromString(s)
    assert( @@n != nil )
    @@f = LibSBML::formulaToString(@@n)
    assert ((  "arccsc(x)" == @@f ))
  end

  def test_element_arccsch
    s = wrapMathML("<apply><arccsch/><ci> x </ci></apply>")
    @@n = LibSBML::readMathMLFromString(s)
    assert( @@n != nil )
    @@f = LibSBML::formulaToString(@@n)
    assert ((  "arccsch(x)" == @@f ))
  end

  def test_element_arcsec
    s = wrapMathML("<apply><arcsec/><ci> x </ci></apply>")
    @@n = LibSBML::readMathMLFromString(s)
    assert( @@n != nil )
    @@f = LibSBML::formulaToString(@@n)
    assert ((  "arcsec(x)" == @@f ))
  end

  def test_element_arcsech
    s = wrapMathML("<apply><arcsech/><ci> x </ci></apply>")
    @@n = LibSBML::readMathMLFromString(s)
    assert( @@n != nil )
    @@f = LibSBML::formulaToString(@@n)
    assert ((  "arcsech(x)" == @@f ))
  end

  def test_element_arcsin
    s = wrapMathML("<apply><arcsin/><ci> x </ci></apply>")
    @@n = LibSBML::readMathMLFromString(s)
    assert( @@n != nil )
    @@f = LibSBML::formulaToString(@@n)
    assert ((  "asin(x)" == @@f ))
  end

  def test_element_arcsinh
    s = wrapMathML("<apply><arcsinh/><ci> x </ci></apply>")
    @@n = LibSBML::readMathMLFromString(s)
    assert( @@n != nil )
    @@f = LibSBML::formulaToString(@@n)
    assert ((  "arcsinh(x)" == @@f ))
  end

  def test_element_arctan
    s = wrapMathML("<apply><arctan/><ci> x </ci></apply>")
    @@n = LibSBML::readMathMLFromString(s)
    assert( @@n != nil )
    @@f = LibSBML::formulaToString(@@n)
    assert ((  "atan(x)" == @@f ))
  end

  def test_element_arctanh
    s = wrapMathML("<apply><arctanh/><ci> x </ci></apply>")
    @@n = LibSBML::readMathMLFromString(s)
    assert( @@n != nil )
    @@f = LibSBML::formulaToString(@@n)
    assert ((  "arctanh(x)" == @@f ))
  end

  def test_element_bug_apply_ci_1
    s = wrapMathML("<apply>" + 
    "  <ci> Y </ci>" + 
    "  <cn> 1 </cn>" + 
    "</apply>")
    @@n = LibSBML::readMathMLFromString(s)
    assert( @@n != nil )
    assert( @@n.getType() == LibSBML::AST_FUNCTION )
    assert ((  "Y" == @@n.getName() ))
    assert( @@n.getNumChildren() == 1 )
    c = @@n.getLeftChild()
    assert( c != nil )
    assert( c.getType() == LibSBML::AST_REAL )
    assert( c.getReal() == 1 )
    assert( c.getNumChildren() == 0 )
  end

  def test_element_bug_apply_ci_2
    s = wrapMathML("<apply>" + 
    "  <ci> Y </ci>" + 
    "  <csymbol encoding='text' " + 
    "   definitionURL='http://www.sbml.org/sbml/symbols/time'> t </csymbol>" + 
    "</apply>")
    @@n = LibSBML::readMathMLFromString(s)
    assert( @@n != nil )
    assert( @@n.getType() == LibSBML::AST_FUNCTION )
    assert ((  "Y" == @@n.getName() ))
    assert( @@n.getNumChildren() == 1 )
    c = @@n.getLeftChild()
    assert( c != nil )
    assert( c.getType() == LibSBML::AST_NAME_TIME )
    assert ((  "t" == c.getName() ))
    assert( c.getNumChildren() == 0 )
  end

  def test_element_bug_cn_e_notation_1
    s = wrapMathML("<cn type='e-notation'> 2 <sep/> -8 </cn>")
    @@n = LibSBML::readMathMLFromString(s)
    assert( @@n != nil )
    assert( @@n.getType() == LibSBML::AST_REAL_E )
    assert( @@n.getMantissa() == 2.0 )
    assert( @@n.getExponent() == -8.0 )
    assert( @@n.getNumChildren() == 0 )
  end

  def test_element_bug_cn_e_notation_2
    s = wrapMathML("<cn type='e-notation'> -3 <sep/> 4 </cn>")
    @@n = LibSBML::readMathMLFromString(s)
    assert( @@n != nil )
    assert( @@n.getType() == LibSBML::AST_REAL_E )
    assert( @@n.getMantissa() == -3.0 )
    assert( @@n.getExponent() == 4.0 )
    assert( @@n.getNumChildren() == 0 )
  end

  def test_element_bug_cn_e_notation_3
    s = wrapMathML("<cn type='e-notation'> -6 <sep/> -1 </cn>")
    @@n = LibSBML::readMathMLFromString(s)
    assert( @@n != nil )
    assert( @@n.getType() == LibSBML::AST_REAL_E )
    assert( @@n.getMantissa() == -6.0 )
    assert( @@n.getExponent() == -1.0 )
    assert( @@n.getNumChildren() == 0 )
  end

  def test_element_bug_cn_integer_negative
    s = wrapMathML("<cn type='integer'> -7 </cn>")
    @@n = LibSBML::readMathMLFromString(s)
    assert( @@n != nil )
    assert( @@n.getType() == LibSBML::AST_INTEGER )
    assert( @@n.getInteger() == -7 )
    assert( @@n.getNumChildren() == 0 )
  end

  def test_element_bug_csymbol_1
    s = wrapMathML("<apply>" + 
    "  <gt/>" + 
    "  <csymbol encoding='text' " + 
    "    definitionURL='http://www.sbml.org/sbml/symbols/time'>time</csymbol>" + 
    "  <cn>5000</cn>" + 
    "</apply>")
    @@n = LibSBML::readMathMLFromString(s)
    assert( @@n != nil )
    assert( @@n.getType() == LibSBML::AST_RELATIONAL_GT )
    assert( @@n.getNumChildren() == 2 )
    c = @@n.getLeftChild()
    assert( c != nil )
    assert( c.getType() == LibSBML::AST_NAME_TIME )
    assert ((  "time" == c.getName() ))
    assert( c.getNumChildren() == 0 )
    c = @@n.getRightChild()
    assert( c != nil )
    assert( c.getType() == LibSBML::AST_REAL )
    assert( c.getReal() == 5000 )
    assert( c.getNumChildren() == 0 )
  end

  def test_element_bug_csymbol_delay_1
    s = wrapMathML("<apply>" + 
    "  <csymbol encoding='text' definitionURL='http://www.sbml.org/sbml/" + 
    "symbols/delay'> my_delay </csymbol>" + 
    "  <ci> x </ci>" + 
    "  <cn> 0.1 </cn>" + 
    "</apply>\n")
    @@n = LibSBML::readMathMLFromString(s)
    assert( @@n != nil )
    assert( @@n.getType() == LibSBML::AST_FUNCTION_DELAY )
    assert ((  "my_delay" == @@n.getName() ))
    assert( @@n.getNumChildren() == 2 )
    c = @@n.getLeftChild()
    assert( c != nil )
    assert( c.getType() == LibSBML::AST_NAME )
    assert ((  "x" == c.getName() ))
    assert( c.getNumChildren() == 0 )
    c = @@n.getRightChild()
    assert( c != nil )
    assert( c.getType() == LibSBML::AST_REAL )
    assert( c.getReal() == 0.1 )
    assert( c.getNumChildren() == 0 )
  end

  def test_element_bug_math_xmlns
    s = wrapXML("<foo:math xmlns:foo='http://www.w3.org/1998/Math/MathML'>" + 
    "  <foo:apply>" + 
    "    <foo:plus/> <foo:cn>1</foo:cn> <foo:cn>2</foo:cn>" + 
    "  </foo:apply>" + 
    "</foo:math>")
    @@n = LibSBML::readMathMLFromString(s)
    assert( @@n != nil )
    @@f = LibSBML::formulaToString(@@n)
    assert ((  "1 + 2" == @@f ))
  end

  def test_element_ceiling
    s = wrapMathML("<apply><ceiling/><cn> 1.6 </cn></apply>")
    @@n = LibSBML::readMathMLFromString(s)
    assert( @@n != nil )
    @@f = LibSBML::formulaToString(@@n)
    assert ((  "ceil(1.6)" == @@f ))
  end

  def test_element_ci
    s = wrapMathML("<ci> x </ci>")
    @@n = LibSBML::readMathMLFromString(s)
    assert( @@n != nil )
    assert( @@n.getType() == LibSBML::AST_NAME )
    assert ((  "x" == @@n.getName() ))
    assert( @@n.getNumChildren() == 0 )
  end

  def test_element_ci_definitionURL
    s = wrapMathML("<ci definitionURL=\"foobar\"> x </ci>")
    @@n = LibSBML::readMathMLFromString(s)
    assert( @@n != nil )
    assert( @@n.getType() == LibSBML::AST_NAME )
    assert ((  "x" == @@n.getName() ))
    assert( @@n.getNumChildren() == 0 )
    assert( @@n.getDefinitionURL().getValue(0) ==  "foobar" )
  end

  def test_element_ci_surrounding_spaces_bug
    s = wrapMathML("  <ci> s </ci>  ")
    @@n = LibSBML::readMathMLFromString(s)
    assert( @@n != nil )
    assert( @@n.getType() == LibSBML::AST_NAME )
    assert ((  "s" == @@n.getName() ))
    assert( @@n.getNumChildren() == 0 )
  end

  def test_element_cn_default
    s = wrapMathML("<cn> 12345.7 </cn>")
    @@n = LibSBML::readMathMLFromString(s)
    assert( @@n != nil )
    assert( @@n.getType() == LibSBML::AST_REAL )
    assert( @@n.getReal() == 12345.7 )
    assert( @@n.getNumChildren() == 0 )
  end

  def test_element_cn_e_notation
    s = wrapMathML("<cn type='e-notation'> 12.3 <sep/> 5 </cn>")
    @@n = LibSBML::readMathMLFromString(s)
    assert( @@n != nil )
    assert( @@n.getType() == LibSBML::AST_REAL_E )
    assert( @@n.getMantissa() == 12.3 )
    assert( @@n.getExponent() == 5 )
    assert( @@n.getNumChildren() == 0 )
  end

  def test_element_cn_integer
    s = wrapMathML("<cn type='integer'> 12345 </cn>")
    @@n = LibSBML::readMathMLFromString(s)
    assert( @@n != nil )
    assert( @@n.getType() == LibSBML::AST_INTEGER )
    assert( @@n.getInteger() == 12345 )
    assert( @@n.getNumChildren() == 0 )
  end

  def test_element_cn_rational
    s = wrapMathML("<cn type='rational'> 12342 <sep/> 2342342 </cn>"  
    )
    @@n = LibSBML::readMathMLFromString(s)
    assert( @@n != nil )
    assert( @@n.getType() == LibSBML::AST_RATIONAL )
    assert( @@n.getNumerator() == 12342 )
    assert( @@n.getDenominator() == 2342342 )
    assert( @@n.getNumChildren() == 0 )
  end

  def test_element_cn_real
    s = wrapMathML("<cn type='real'> 12345.7 </cn>")
    @@n = LibSBML::readMathMLFromString(s)
    assert( @@n != nil )
    assert( @@n.getType() == LibSBML::AST_REAL )
    assert( @@n.getReal() == 12345.7 )
    assert( @@n.getNumChildren() == 0 )
  end

  def test_element_cn_units
    s = wrapMathMLUnits("<cn sbml:units=\"mole\"> 12345.7 </cn>")
    @@n = LibSBML::readMathMLFromString(s)
    assert( @@n != nil )
    assert( @@n.getType() == LibSBML::AST_REAL )
    assert( @@n.getReal() == 12345.7 )
    assert( @@n.getUnits() ==  "mole"    )
    assert( @@n.getNumChildren() == 0 )
  end

  def test_element_constants_exponentiale
    s = wrapMathML("<exponentiale/>")
    @@n = LibSBML::readMathMLFromString(s)
    assert( @@n != nil )
    assert( @@n.getType() == LibSBML::AST_CONSTANT_E )
    assert( @@n.getNumChildren() == 0 )
  end

  def test_element_constants_false
    s = wrapMathML("<false/>")
    @@n = LibSBML::readMathMLFromString(s)
    assert( @@n != nil )
    assert( @@n.getType() == LibSBML::AST_CONSTANT_FALSE )
    assert( @@n.getNumChildren() == 0 )
  end

  def test_element_constants_infinity
    s = wrapMathML("<infinity/>")
    @@n = LibSBML::readMathMLFromString(s)
    assert( @@n != nil )
    assert( @@n.getType() == LibSBML::AST_REAL )
    assert( util_isInf(@@n.getReal()) == true )
    assert( @@n.getNumChildren() == 0 )
  end

  def test_element_constants_notanumber
    s = wrapMathML("<notanumber/>")
    @@n = LibSBML::readMathMLFromString(s)
    assert( @@n != nil )
    assert( @@n.getType() == LibSBML::AST_REAL )
    assert_equal true, isnan(@@n.getReal())
    assert( @@n.getNumChildren() == 0 )
  end

  def test_element_constants_pi
    s = wrapMathML("<pi/>")
    @@n = LibSBML::readMathMLFromString(s)
    assert( @@n != nil )
    assert( @@n.getType() == LibSBML::AST_CONSTANT_PI )
    assert( @@n.getNumChildren() == 0 )
  end

  def test_element_constants_true
    s = wrapMathML("<true/>")
    @@n = LibSBML::readMathMLFromString(s)
    assert( @@n != nil )
    assert( @@n.getType() == LibSBML::AST_CONSTANT_TRUE )
    assert( @@n.getNumChildren() == 0 )
  end

  def test_element_cos
    s = wrapMathML("<apply><cos/><ci> x </ci></apply>")
    @@n = LibSBML::readMathMLFromString(s)
    assert( @@n != nil )
    @@f = LibSBML::formulaToString(@@n)
    assert ((  "cos(x)" == @@f ))
  end

  def test_element_cosh
    s = wrapMathML("<apply><cosh/><ci> x </ci></apply>")
    @@n = LibSBML::readMathMLFromString(s)
    assert( @@n != nil )
    @@f = LibSBML::formulaToString(@@n)
    assert ((  "cosh(x)" == @@f ))
  end

  def test_element_cot
    s = wrapMathML("<apply><cot/><ci> x </ci></apply>")
    @@n = LibSBML::readMathMLFromString(s)
    assert( @@n != nil )
    @@f = LibSBML::formulaToString(@@n)
    assert ((  "cot(x)" == @@f ))
  end

  def test_element_coth
    s = wrapMathML("<apply><coth/><ci> x </ci></apply>")
    @@n = LibSBML::readMathMLFromString(s)
    assert( @@n != nil )
    @@f = LibSBML::formulaToString(@@n)
    assert ((  "coth(x)" == @@f ))
  end

  def test_element_csc
    s = wrapMathML("<apply><csc/><ci> x </ci></apply>")
    @@n = LibSBML::readMathMLFromString(s)
    assert( @@n != nil )
    @@f = LibSBML::formulaToString(@@n)
    assert ((  "csc(x)" == @@f ))
  end

  def test_element_csch
    s = wrapMathML("<apply><csch/><ci> x </ci></apply>")
    @@n = LibSBML::readMathMLFromString(s)
    assert( @@n != nil )
    @@f = LibSBML::formulaToString(@@n)
    assert ((  "csch(x)" == @@f ))
  end

  def test_element_csymbol_avogadro
    s = wrapMathML("<csymbol encoding='text' " + "definitionURL='http://www.sbml.org/sbml/symbols/avogadro'> NA </csymbol>")
    @@n = LibSBML::readMathMLFromString(s)
    assert( @@n != nil )
    assert( @@n.getType() == LibSBML::AST_NAME_AVOGADRO )
    assert ((  "NA" == @@n.getName() ))
    assert( @@n.getNumChildren() == 0 )
  end

  def test_element_csymbol_delay_1
    s = wrapMathML("<csymbol encoding='text' " + "definitionURL='http://www.sbml.org/sbml/symbols/delay'> delay </csymbol>")
    @@n = LibSBML::readMathMLFromString(s)
    assert( @@n != nil )
    assert( @@n.getType() == LibSBML::AST_FUNCTION_DELAY )
    assert ((  "delay" == @@n.getName() ))
    assert( @@n.getNumChildren() == 0 )
  end

  def test_element_csymbol_delay_2
    s = wrapMathML("<apply>" + 
    "  <csymbol encoding='text' definitionURL='http://www.sbml.org/sbml/" + 
    "symbols/delay'> my_delay </csymbol>" + 
    "  <ci> x </ci>" + 
    "  <cn> 0.1 </cn>" + 
    "</apply>\n")
    @@n = LibSBML::readMathMLFromString(s)
    assert( @@n != nil )
    @@f = LibSBML::formulaToString(@@n)
    assert ((  "my_delay(x, 0.1)" == @@f ))
  end

  def test_element_csymbol_delay_3
    s = wrapMathML("<apply>" + 
    "  <power/>" + 
    "  <apply>" + 
    "    <csymbol encoding='text' definitionURL='http://www.sbml.org/sbml/" + 
    "symbols/delay'> delay </csymbol>" + 
    "    <ci> P </ci>" + 
    "    <ci> delta_t </ci>" + 
    "  </apply>\n" + 
    "  <ci> q </ci>" + 
    "</apply>\n")
    @@n = LibSBML::readMathMLFromString(s)
    assert( @@n != nil )
    @@f = LibSBML::formulaToString(@@n)
    assert ((  "pow(delay(P, delta_t), q)" == @@f ))
  end

  def test_element_csymbol_time
    s = wrapMathML("<csymbol encoding='text' " + "definitionURL='http://www.sbml.org/sbml/symbols/time'> t </csymbol>")
    @@n = LibSBML::readMathMLFromString(s)
    assert( @@n != nil )
    assert( @@n.getType() == LibSBML::AST_NAME_TIME )
    assert ((  "t" == @@n.getName() ))
    assert( @@n.getNumChildren() == 0 )
  end

  def test_element_eq
    s = wrapMathML("<apply> <eq/> <ci>a</ci> <ci>b</ci> <ci>c</ci> </apply>"  
    )
    @@n = LibSBML::readMathMLFromString(s)
    assert( @@n != nil )
    @@f = LibSBML::formulaToString(@@n)
    assert ((  "eq(a, b, c)" == @@f ))
  end

  def test_element_exp
    s = wrapMathML("<apply><exp/><ci> x </ci></apply>")
    @@n = LibSBML::readMathMLFromString(s)
    assert( @@n != nil )
    @@f = LibSBML::formulaToString(@@n)
    assert ((  "exp(x)" == @@f ))
  end

  def test_element_factorial
    s = wrapMathML("<apply><factorial/><cn> 5 </cn></apply>")
    @@n = LibSBML::readMathMLFromString(s)
    assert( @@n != nil )
    @@f = LibSBML::formulaToString(@@n)
    assert ((  "factorial(5)" == @@f ))
  end

  def test_element_floor
    s = wrapMathML("<apply><floor/><cn> 1.2 </cn></apply>")
    @@n = LibSBML::readMathMLFromString(s)
    assert( @@n != nil )
    @@f = LibSBML::formulaToString(@@n)
    assert ((  "floor(1.2)" == @@f ))
  end

  def test_element_function_call_1
    s = wrapMathML("<apply> <ci> foo </ci> <ci> x </ci> </apply>")
    @@n = LibSBML::readMathMLFromString(s)
    assert( @@n != nil )
    @@f = LibSBML::formulaToString(@@n)
    assert ((  "foo(x)" == @@f ))
  end

  def test_element_function_call_2
    s = wrapMathML("<apply> <plus/> <cn> 1 </cn>" + 
    "                <apply> <ci> f </ci> <ci> x </ci> </apply>" + 
    "</apply>")
    @@n = LibSBML::readMathMLFromString(s)
    assert( @@n != nil )
    @@f = LibSBML::formulaToString(@@n)
    assert ((  "1 + f(x)" == @@f ))
  end

  def test_element_geq
    s = wrapMathML("<apply> <geq/> <cn>1</cn> <ci>x</ci> <cn>0</cn> </apply>"  
    )
    @@n = LibSBML::readMathMLFromString(s)
    assert( @@n != nil )
    @@f = LibSBML::formulaToString(@@n)
    assert ((  "geq(1, x, 0)" == @@f ))
  end

  def test_element_gt
    s = wrapMathML("<apply> <gt/> <infinity/>" + 
    "              <apply> <minus/> <infinity/> <cn>1</cn> </apply>" + 
    "</apply>")
    @@n = LibSBML::readMathMLFromString(s)
    assert( @@n != nil )
    @@f = LibSBML::formulaToString(@@n)
    assert ((  "gt(INF, INF - 1)" == @@f ))
  end

  def test_element_invalid_mathml
    invalid = wrapMathML("<lambda definitionURL=\"http://biomodels.net/SBO/#SBO:0000065\">" + 
    "<bvar>" + 
    "<ci>c</ci>" + 
    "</bvar>" + 
    "<apply>" + 
    "  <ci>c</ci>" + 
    "</apply>" + 
    "</lambda>\n")
    @@n = LibSBML::readMathMLFromString(nil)
    assert( @@n == nil )
    @@n = LibSBML::readMathMLFromString(invalid)
    assert( @@n == nil )
  end

  def test_element_lambda
    s = wrapMathML("<lambda>" + 
    "  <bvar> <ci>x</ci> </bvar>" + 
    "  <apply> <sin/>" + 
    "          <apply> <plus/> <ci>x</ci> <cn>1</cn> </apply>" + 
    "  </apply>" + 
    "</lambda>")
    @@n = LibSBML::readMathMLFromString(s)
    assert( @@n != nil )
    @@f = LibSBML::formulaToString(@@n)
    assert ((  "lambda(x, sin(x + 1))" == @@f ))
  end

  def test_element_leq
    s = wrapMathML("<apply> <leq/> <cn>0</cn> <ci>x</ci> <cn>1</cn> </apply>"  
    )
    @@n = LibSBML::readMathMLFromString(s)
    assert( @@n != nil )
    @@f = LibSBML::formulaToString(@@n)
    assert ((  "leq(0, x, 1)" == @@f ))
  end

  def test_element_ln
    s = wrapMathML("<apply><ln/><ci> a </ci></apply>")
    @@n = LibSBML::readMathMLFromString(s)
    assert( @@n != nil )
    @@f = LibSBML::formulaToString(@@n)
    assert ((  "log(a)" == @@f ))
  end

  def test_element_log_1
    s = wrapMathML("<apply> <log/> <logbase> <cn type='integer'> 3 </cn> </logbase>" + 
    "               <ci> x </ci>" + 
    "</apply>")
    @@n = LibSBML::readMathMLFromString(s)
    assert( @@n != nil )
    @@f = LibSBML::formulaToString(@@n)
    assert ((  "log(3, x)" == @@f ))
  end

  def test_element_log_2
    s = wrapMathML("<apply> <log/> <ci> x </ci> </apply>")
    @@n = LibSBML::readMathMLFromString(s)
    assert( @@n != nil )
    @@f = LibSBML::formulaToString(@@n)
    assert ((  "log10(x)" == @@f ))
  end

  def test_element_lt
    s = wrapMathML("<apply> <lt/> <apply> <minus/> <infinity/> <infinity/> </apply>" + 
    "              <cn>1</cn>" + 
    "</apply>")
    @@n = LibSBML::readMathMLFromString(s)
    assert( @@n != nil )
    @@f = LibSBML::formulaToString(@@n)
    assert ((  "lt(INF - INF, 1)" == @@f ))
  end

  def test_element_math
    s = wrapXML("<math xmlns='http://www.w3.org/1998/Math/MathML'/>"  
    )
    @@n = LibSBML::readMathMLFromString(s)
    assert( @@n != nil )
    assert( @@n.getType() == LibSBML::AST_UNKNOWN )
  end

  def test_element_neq
    s = wrapMathML("<apply> <neq/> <notanumber/> <notanumber/> </apply>"  
    )
    @@n = LibSBML::readMathMLFromString(s)
    assert( @@n != nil )
    @@f = LibSBML::formulaToString(@@n)
    assert ((  "neq(NaN, NaN)" == @@f ))
  end

  def test_element_not
    s = wrapMathML("<apply> <not/> <ci> TooShabby </ci> </apply>")
    @@n = LibSBML::readMathMLFromString(s)
    assert( @@n != nil )
    @@f = LibSBML::formulaToString(@@n)
    assert ((  "not(TooShabby)" == @@f ))
  end

  def test_element_operator_plus
    s = wrapMathML("<apply> <plus/> <cn> 1 </cn> <cn> 2 </cn> <cn> 3 </cn> </apply>"  
    )
    @@n = LibSBML::readMathMLFromString(s)
    assert( @@n != nil )
    @@f = LibSBML::formulaToString(@@n)
    assert ((  "1 + 2 + 3" == @@f ))
  end

  def test_element_operator_times
    s = wrapMathML("<apply> <times/> <ci> x </ci> <ci> y </ci> <ci> z </ci> </apply>"  
    )
    @@n = LibSBML::readMathMLFromString(s)
    assert( @@n != nil )
    @@f = LibSBML::formulaToString(@@n)
    assert ((  "x * y * z" == @@f ))
  end

  def test_element_or
    s = wrapMathML("<apply> <or/> <ci>a</ci> <ci>b</ci> <ci>c</ci> <ci>d</ci> </apply>"  
    )
    @@n = LibSBML::readMathMLFromString(s)
    assert( @@n != nil )
    @@f = LibSBML::formulaToString(@@n)
    assert ((  "or(a, b, c, d)" == @@f ))
  end

  def test_element_piecewise
    s = wrapMathML("<piecewise>" + 
    "  <piece>" + 
    "    <apply> <minus/> <ci>x</ci> </apply>" + 
    "    <apply> <lt/> <ci>x</ci> <cn>0</cn> </apply>" + 
    "  </piece>" + 
    "  <piece>" + 
    "    <cn>0</cn>" + 
    "    <apply> <eq/> <ci>x</ci> <cn>0</cn> </apply>" + 
    "  </piece>" + 
    "  <piece>" + 
    "    <ci>x</ci>" + 
    "    <apply> <gt/> <ci>x</ci> <cn>0</cn> </apply>" + 
    "  </piece>" + 
    "</piecewise>")
    @@n = LibSBML::readMathMLFromString(s)
    assert( @@n != nil )
    @@f = LibSBML::formulaToString(@@n)
    assert ((  "piecewise(-x, lt(x, 0), 0, eq(x, 0), x, gt(x, 0))" == @@f ))
  end

  def test_element_piecewise_otherwise
    s = wrapMathML("<piecewise>" + 
    "  <piece>" + 
    "    <cn>0</cn>" + 
    "    <apply> <lt/> <ci>x</ci> <cn>0</cn> </apply>" + 
    "  </piece>" + 
    "  <otherwise>" + 
    "    <ci>x</ci>" + 
    "  </otherwise>" + 
    "</piecewise>")
    @@n = LibSBML::readMathMLFromString(s)
    assert( @@n != nil )
    @@f = LibSBML::formulaToString(@@n)
    assert ((  "piecewise(0, lt(x, 0), x)" == @@f ))
  end

  def test_element_power
    s = wrapMathML("<apply><power/> <ci>x</ci> <cn>3</cn> </apply>")
    @@n = LibSBML::readMathMLFromString(s)
    assert( @@n != nil )
    @@f = LibSBML::formulaToString(@@n)
    assert ((  "pow(x, 3)" == @@f ))
  end

  def test_element_root_1
    s = wrapMathML("<apply> <root/> <degree> <cn type='integer'> 3 </cn> </degree>" + 
    "               <ci> a </ci>" + 
    "</apply>")
    @@n = LibSBML::readMathMLFromString(s)
    assert( @@n != nil )
    @@f = LibSBML::formulaToString(@@n)
    assert ((  "root(3, a)" == @@f ))
  end

  def test_element_root_2
    s = wrapMathML("<apply> <root/> <ci> a </ci> </apply>")
    @@n = LibSBML::readMathMLFromString(s)
    assert( @@n != nil )
    @@f = LibSBML::formulaToString(@@n)
    assert ((  "sqrt(a)" == @@f ))
  end

  def test_element_sec
    s = wrapMathML("<apply><sec/><ci> x </ci></apply>")
    @@n = LibSBML::readMathMLFromString(s)
    assert( @@n != nil )
    @@f = LibSBML::formulaToString(@@n)
    assert ((  "sec(x)" == @@f ))
  end

  def test_element_sech
    s = wrapMathML("<apply><sech/><ci> x </ci></apply>")
    @@n = LibSBML::readMathMLFromString(s)
    assert( @@n != nil )
    @@f = LibSBML::formulaToString(@@n)
    assert ((  "sech(x)" == @@f ))
  end

  def test_element_sin
    s = wrapMathML("<apply><sin/><ci> x </ci></apply>")
    @@n = LibSBML::readMathMLFromString(s)
    assert( @@n != nil )
    @@f = LibSBML::formulaToString(@@n)
    assert ((  "sin(x)" == @@f ))
  end

  def test_element_sinh
    s = wrapMathML("<apply><sinh/><ci> x </ci></apply>")
    @@n = LibSBML::readMathMLFromString(s)
    assert( @@n != nil )
    @@f = LibSBML::formulaToString(@@n)
    assert ((  "sinh(x)" == @@f ))
  end

  def test_element_tan
    s = wrapMathML("<apply><tan/><ci> x </ci></apply>")
    @@n = LibSBML::readMathMLFromString(s)
    assert( @@n != nil )
    @@f = LibSBML::formulaToString(@@n)
    assert ((  "tan(x)" == @@f ))
  end

  def test_element_tanh
    s = wrapMathML("<apply><tanh/><ci> x </ci></apply>")
    @@n = LibSBML::readMathMLFromString(s)
    assert( @@n != nil )
    @@f = LibSBML::formulaToString(@@n)
    assert ((  "tanh(x)" == @@f ))
  end

  def test_element_xor
    s = wrapMathML("<apply> <xor/> <ci>a</ci> <ci>b</ci> <ci>b</ci> <ci>a</ci> </apply>"  
    )
    @@n = LibSBML::readMathMLFromString(s)
    assert( @@n != nil )
    @@f = LibSBML::formulaToString(@@n)
    assert ((  "xor(a, b, b, a)" == @@f ))
  end

end
