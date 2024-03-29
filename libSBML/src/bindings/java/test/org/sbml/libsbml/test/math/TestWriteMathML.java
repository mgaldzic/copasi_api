/*
 * @file    TestWriteMathML.java
 * @brief   Write MathML unit tests
 *
 * @author  Akiya Jouraku (Java conversion)
 * @author  Ben Bornstein 
 *
 * $Id$
 * $HeadURL$
 *
 * ====== WARNING ===== WARNING ===== WARNING ===== WARNING ===== WARNING ======
 *
 * DO NOT EDIT THIS FILE.
 *
 * This file was generated automatically by converting the file located at
 * src/math/test/TestWriteMathML.cpp
 * using the conversion program dev/utilities/translateTests/translateTests.pl.
 * Any changes made here will be lost the next time the file is regenerated.
 *
 * -----------------------------------------------------------------------------
 * This file is part of libSBML.  Please visit http://sbml.org for more
 * information about SBML, and the latest version of libSBML.
 *
 * Copyright 2005-2010 California Institute of Technology.
 * Copyright 2002-2005 California Institute of Technology and
 *                     Japan Science and Technology Corporation.
 * 
 * This library is free software; you can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation.  A copy of the license agreement is provided
 * in the file named "LICENSE.txt" included with this software distribution
 * and also available online as http://sbml.org/software/libsbml/license.html
 * -----------------------------------------------------------------------------
 */

package org.sbml.libsbml.test.math;

import org.sbml.libsbml.*;

import java.io.File;
import java.lang.AssertionError;

public class TestWriteMathML {

  static void assertTrue(boolean condition) throws AssertionError
  {
    if (condition == true)
    {
      return;
    }
    throw new AssertionError();
  }

  static void assertEquals(Object a, Object b) throws AssertionError
  {
    if ( (a == null) && (b == null) )
    {
      return;
    }
    else if ( (a == null) || (b == null) )
    {
      throw new AssertionError();
    }
    else if (a.equals(b))
    {
      return;
    }

    throw new AssertionError();
  }

  static void assertNotEquals(Object a, Object b) throws AssertionError
  {
    if ( (a == null) && (b == null) )
    {
      throw new AssertionError();
    }
    else if ( (a == null) || (b == null) )
    {
      return;
    }
    else if (a.equals(b))
    {
      throw new AssertionError();
    }
  }

  static void assertEquals(boolean a, boolean b) throws AssertionError
  {
    if ( a == b )
    {
      return;
    }
    throw new AssertionError();
  }

  static void assertNotEquals(boolean a, boolean b) throws AssertionError
  {
    if ( a != b )
    {
      return;
    }
    throw new AssertionError();
  }

  static void assertEquals(int a, int b) throws AssertionError
  {
    if ( a == b )
    {
      return;
    }
    throw new AssertionError();
  }

  static void assertNotEquals(int a, int b) throws AssertionError
  {
    if ( a != b )
    {
      return;
    }
    throw new AssertionError();
  }
  private String S;
  private ASTNode N;

  public String MATHML_FOOTER()
  {
    return "</math>";
  }

  public String MATHML_HEADER()
  {
    return "<math xmlns=\"http://www.w3.org/1998/Math/MathML\">\n";
  }

  public String MATHML_HEADER_UNITS()
  {
    return "<math xmlns=\"http://www.w3.org/1998/Math/MathML\"";
  }

  public String MATHML_HEADER_UNITS2()
  {
    return " xmlns:sbml=\"http://www.sbml.org/sbml/level3/version1/core\">\n";
  }

  public String XML_HEADER()
  {
    return "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n";
  }

  public String wrapMathML(String s)
  {
    String r = XML_HEADER();
    r += MATHML_HEADER();
    r += s;
    r += MATHML_FOOTER();
    return r;
  }

  public String wrapMathMLUnits(String s)
  {
    String r = XML_HEADER();
    r += MATHML_HEADER_UNITS();
    r += MATHML_HEADER_UNITS2();
    r += s;
    r += MATHML_FOOTER();
    return r;
  }


  public double util_NaN()
  {
    double z = 0.0;
    return 0.0/z;
  }

  public double util_PosInf()
  {
    double z = 0.0;
    return 1.0/z;
  }

  public double util_NegInf()
  {
    double z = 0.0;
    return -1.0/z;
  }

  public boolean equals(String s1, String s2)
  {
    return s1.equals(s2);
  }

  protected void setUp() throws Exception
  {
    N = null;
    S = null;
  }

  protected void tearDown() throws Exception
  {
    S = null;
  }

  public void test_MathMLFormatter_ci()
  {
    String expected = wrapMathML("  <ci> foo </ci>\n");
    N = libsbml.parseFormula("foo");
    S = libsbml.writeMathMLToString(N);
    assertEquals( true, equals(expected,S) );
  }

  public void test_MathMLFormatter_cn_e_notation_1()
  {
    String expected = wrapMathML("  <cn type=\"e-notation\"> 0 <sep/> 3 </cn>\n"  
    );
    N = libsbml.parseFormula("0e3");
    S = libsbml.writeMathMLToString(N);
    assertEquals( true, equals(expected,S) );
  }

  public void test_MathMLFormatter_cn_e_notation_2()
  {
    String expected = wrapMathML("  <cn type=\"e-notation\"> 2 <sep/> 3 </cn>\n"  
    );
    N = libsbml.parseFormula("2e3");
    S = libsbml.writeMathMLToString(N);
    assertEquals( true, equals(expected,S) );
  }

  public void test_MathMLFormatter_cn_e_notation_3()
  {
    String expected = wrapMathML("  <cn type=\"e-notation\"> 1234567.8 <sep/> 3 </cn>\n"  
    );
    N = libsbml.parseFormula("1234567.8e3");
    S = libsbml.writeMathMLToString(N);
    assertEquals( true, equals(expected,S) );
  }

  public void test_MathMLFormatter_cn_e_notation_4()
  {
    String expected = wrapMathML("  <cn type=\"e-notation\"> 6.0221367 <sep/> 23 </cn>\n"  
    );
    N = libsbml.parseFormula("6.0221367e+23");
    S = libsbml.writeMathMLToString(N);
    assertEquals( true, equals(expected,S) );
  }

  public void test_MathMLFormatter_cn_e_notation_5()
  {
    String expected = wrapMathML("  <cn type=\"e-notation\"> 4 <sep/> -6 </cn>\n"  
    );
    N = libsbml.parseFormula(".000004");
    S = libsbml.writeMathMLToString(N);
    assertEquals( true, equals(expected,S) );
  }

  public void test_MathMLFormatter_cn_e_notation_6()
  {
    String expected = wrapMathML("  <cn type=\"e-notation\"> 4 <sep/> -12 </cn>\n"  
    );
    N = libsbml.parseFormula(".000004e-6");
    S = libsbml.writeMathMLToString(N);
    assertEquals( true, equals(expected,S) );
  }

  public void test_MathMLFormatter_cn_e_notation_7()
  {
    String expected = wrapMathML("  <cn type=\"e-notation\"> -1 <sep/> -6 </cn>\n"  
    );
    N = libsbml.parseFormula("-1e-6");
    S = libsbml.writeMathMLToString(N);
    assertEquals( true, equals(expected,S) );
  }

  public void test_MathMLFormatter_cn_integer()
  {
    String expected = wrapMathML("  <cn type=\"integer\"> 5 </cn>\n");
    N = libsbml.parseFormula("5");
    S = libsbml.writeMathMLToString(N);
    assertEquals( true, equals(expected,S) );
  }

  public void test_MathMLFormatter_cn_rational()
  {
    String expected = wrapMathML("  <cn type=\"rational\"> 1 <sep/> 3 </cn>\n"  
    );
    N = new ASTNode();
    N.setValue(((int) 1),3);
    S = libsbml.writeMathMLToString(N);
    assertEquals( true, equals(expected,S) );
  }

  public void test_MathMLFormatter_cn_real_1()
  {
    String expected = wrapMathML("  <cn> 1.2 </cn>\n");
    N = libsbml.parseFormula("1.2");
    S = libsbml.writeMathMLToString(N);
    assertEquals( true, equals(expected,S) );
  }

  public void test_MathMLFormatter_cn_real_2()
  {
    String expected = wrapMathML("  <cn> 1234567.8 </cn>\n");
    N = libsbml.parseFormula("1234567.8");
    S = libsbml.writeMathMLToString(N);
    assertEquals( true, equals(expected,S) );
  }

  public void test_MathMLFormatter_cn_real_3()
  {
    String expected = wrapMathML("  <cn> -3.14 </cn>\n");
    N = libsbml.parseFormula("-3.14");
    S = libsbml.writeMathMLToString(N);
    assertEquals( true, equals(expected,S) );
  }

  public void test_MathMLFormatter_cn_real_locale()
  {
    String expected = wrapMathML("  <cn> 2.72 </cn>\n");
    N = libsbml.parseFormula("2.72");
    S = libsbml.writeMathMLToString(N);
    assertEquals( true, equals(expected,S) );
  }

  public void test_MathMLFormatter_cn_units()
  {
    String expected = wrapMathMLUnits("  <cn sbml:units=\"mole\"> 1.2 </cn>\n");
    N = libsbml.parseFormula("1.2");
    N.setUnits("mole");
    S = libsbml.writeMathMLToString(N);
    assertEquals( true, equals(expected,S) );
  }

  public void test_MathMLFormatter_constant_exponentiale()
  {
    String expected = wrapMathML("  <exponentiale/>\n");
    N = new ASTNode(libsbml.AST_CONSTANT_E);
    S = libsbml.writeMathMLToString(N);
    assertEquals( true, equals(expected,S) );
  }

  public void test_MathMLFormatter_constant_false()
  {
    String expected = wrapMathML("  <false/>\n");
    N = new ASTNode(libsbml.AST_CONSTANT_FALSE);
    S = libsbml.writeMathMLToString(N);
    assertEquals( true, equals(expected,S) );
  }

  public void test_MathMLFormatter_constant_infinity()
  {
    String expected = wrapMathML("  <infinity/>\n");
    N = new ASTNode();
    N.setValue( util_PosInf() );
    S = libsbml.writeMathMLToString(N);
    assertEquals( true, equals(expected,S) );
  }

  public void test_MathMLFormatter_constant_infinity_neg()
  {
    String expected = wrapMathML("  <apply> <minus/> <infinity/> </apply>\n"  
    );
    N = new ASTNode();
    N.setValue(- util_PosInf());
    S = libsbml.writeMathMLToString(N);
    assertEquals( true, equals(expected,S) );
  }

  public void test_MathMLFormatter_constant_notanumber()
  {
    String expected = wrapMathML("  <notanumber/>\n");
    N = new ASTNode(libsbml.AST_REAL);
    N.setValue( util_NaN() );
    S = libsbml.writeMathMLToString(N);
    assertEquals( true, equals(expected,S) );
  }

  public void test_MathMLFormatter_constant_true()
  {
    String expected = wrapMathML("  <true/>\n");
    N = new ASTNode(libsbml.AST_CONSTANT_TRUE);
    S = libsbml.writeMathMLToString(N);
    assertEquals( true, equals(expected,S) );
  }

  public void test_MathMLFormatter_csymbol_avogadro()
  {
    String expected = wrapMathML("  <csymbol encoding=\"text\" " + "definitionURL=\"http://www.sbml.org/sbml/symbols/avogadro\"> NA </csymbol>\n");
    N = new ASTNode(libsbml.AST_NAME_AVOGADRO);
    N.setName("NA");
    S = libsbml.writeMathMLToString(N);
    assertEquals( true, equals(expected,S) );
  }

  public void test_MathMLFormatter_csymbol_delay()
  {
    String expected = wrapMathML("  <apply>\n" + 
    "    <csymbol encoding=\"text\" definitionURL=\"http://www.sbml.org/sbml/" + 
    "symbols/delay\"> my_delay </csymbol>\n" + 
    "    <ci> x </ci>\n" + 
    "    <cn> 0.1 </cn>\n" + 
    "  </apply>\n");
    N = libsbml.parseFormula("delay(x, 0.1)");
    N.setName("my_delay");
    S = libsbml.writeMathMLToString(N);
    assertEquals( true, equals(expected,S) );
  }

  public void test_MathMLFormatter_csymbol_time()
  {
    String expected = wrapMathML("  <csymbol encoding=\"text\" " + "definitionURL=\"http://www.sbml.org/sbml/symbols/time\"> t </csymbol>\n");
    N = new ASTNode(libsbml.AST_NAME_TIME);
    N.setName("t");
    S = libsbml.writeMathMLToString(N);
    assertEquals( true, equals(expected,S) );
  }

  public void test_MathMLFormatter_function_1()
  {
    String expected = wrapMathML("  <apply>\n" + 
    "    <ci> foo </ci>\n" + 
    "    <cn type=\"integer\"> 1 </cn>\n" + 
    "    <cn type=\"integer\"> 2 </cn>\n" + 
    "    <cn type=\"integer\"> 3 </cn>\n" + 
    "  </apply>\n");
    N = libsbml.parseFormula("foo(1, 2, 3)");
    S = libsbml.writeMathMLToString(N);
    assertEquals( true, equals(expected,S) );
  }

  public void test_MathMLFormatter_function_2()
  {
    String expected = wrapMathML("  <apply>\n" + 
    "    <ci> foo </ci>\n" + 
    "    <cn type=\"integer\"> 1 </cn>\n" + 
    "    <cn type=\"integer\"> 2 </cn>\n" + 
    "    <apply>\n" + 
    "      <ci> bar </ci>\n" + 
    "      <ci> z </ci>\n" + 
    "    </apply>\n" + 
    "  </apply>\n");
    N = libsbml.parseFormula("foo(1, 2, bar(z))");
    S = libsbml.writeMathMLToString(N);
    assertEquals( true, equals(expected,S) );
  }

  public void test_MathMLFormatter_lambda()
  {
    String expected = wrapMathML("  <lambda>\n" + 
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
    "  </lambda>\n");
    N = libsbml.parseFormula("lambda(x, y, root(2, x^2 + y^2))");
    S = libsbml.writeMathMLToString(N);
    assertEquals( true, equals(expected,S) );
  }

  public void test_MathMLFormatter_lambda_no_bvars()
  {
    String expected = wrapMathML("  <lambda>\n" + 
    "    <apply>\n" + 
    "      <plus/>\n" + 
    "      <cn type=\"integer\"> 2 </cn>\n" + 
    "      <cn type=\"integer\"> 2 </cn>\n" + 
    "    </apply>\n" + 
    "  </lambda>\n");
    N = libsbml.parseFormula("lambda(2 + 2)");
    S = libsbml.writeMathMLToString(N);
    assertEquals( true, equals(expected,S) );
  }

  public void test_MathMLFormatter_log()
  {
    String expected = wrapMathML("  <apply>\n" + 
    "    <log/>\n" + 
    "    <logbase>\n" + 
    "      <cn type=\"integer\"> 2 </cn>\n" + 
    "    </logbase>\n" + 
    "    <ci> N </ci>\n" + 
    "  </apply>\n");
    N = libsbml.parseFormula("log(2, N)");
    S = libsbml.writeMathMLToString(N);
    assertEquals( true, equals(expected,S) );
  }

  public void test_MathMLFormatter_minus()
  {
    String expected = wrapMathML("  <apply>\n" + 
    "    <minus/>\n" + 
    "    <cn type=\"integer\"> 1 </cn>\n" + 
    "    <cn type=\"integer\"> 2 </cn>\n" + 
    "  </apply>\n");
    N = libsbml.parseFormula("1 - 2");
    S = libsbml.writeMathMLToString(N);
    assertEquals( true, equals(expected,S) );
  }

  public void test_MathMLFormatter_minus_unary_1()
  {
    String expected = wrapMathML("  <cn type=\"integer\"> -2 </cn>\n"  
    );
    N = libsbml.parseFormula("-2");
    S = libsbml.writeMathMLToString(N);
    assertEquals( true, equals(expected,S) );
  }

  public void test_MathMLFormatter_minus_unary_2()
  {
    String expected = wrapMathML("  <apply>\n" + 
    "    <minus/>\n" + 
    "    <ci> a </ci>\n" + 
    "  </apply>\n");
    N = libsbml.parseFormula("-a");
    S = libsbml.writeMathMLToString(N);
    assertEquals( true, equals(expected,S) );
  }

  public void test_MathMLFormatter_piecewise()
  {
    String expected = wrapMathML("  <piecewise>\n" + 
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
    "  </piecewise>\n");
    String f =  "piecewise(-x, lt(x, 0), 0, eq(x, 0), x, gt(x, 0))";
    N = libsbml.parseFormula(f);
    S = libsbml.writeMathMLToString(N);
    assertEquals( true, equals(expected,S) );
  }

  public void test_MathMLFormatter_piecewise_otherwise()
  {
    String expected = wrapMathML("  <piecewise>\n" + 
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
    "  </piecewise>\n");
    N = libsbml.parseFormula("piecewise(0, lt(x, 0), x)");
    S = libsbml.writeMathMLToString(N);
    assertEquals( true, equals(expected,S) );
  }

  public void test_MathMLFormatter_plus_binary()
  {
    String expected = wrapMathML("  <apply>\n" + 
    "    <plus/>\n" + 
    "    <cn type=\"integer\"> 1 </cn>\n" + 
    "    <cn type=\"integer\"> 2 </cn>\n" + 
    "  </apply>\n");
    N = libsbml.parseFormula("1 + 2");
    S = libsbml.writeMathMLToString(N);
    assertEquals( true, equals(expected,S) );
  }

  public void test_MathMLFormatter_plus_nary_1()
  {
    String expected = wrapMathML("  <apply>\n" + 
    "    <plus/>\n" + 
    "    <cn type=\"integer\"> 1 </cn>\n" + 
    "    <cn type=\"integer\"> 2 </cn>\n" + 
    "    <cn type=\"integer\"> 3 </cn>\n" + 
    "  </apply>\n");
    N = libsbml.parseFormula("1 + 2 + 3");
    S = libsbml.writeMathMLToString(N);
    assertEquals( true, equals(expected,S) );
  }

  public void test_MathMLFormatter_plus_nary_2()
  {
    String expected = wrapMathML("  <apply>\n" + 
    "    <plus/>\n" + 
    "    <cn type=\"integer\"> 1 </cn>\n" + 
    "    <cn type=\"integer\"> 2 </cn>\n" + 
    "    <cn type=\"integer\"> 3 </cn>\n" + 
    "  </apply>\n");
    N = libsbml.parseFormula("(1 + 2) + 3");
    S = libsbml.writeMathMLToString(N);
    assertEquals( true, equals(expected,S) );
  }

  public void test_MathMLFormatter_plus_nary_3()
  {
    String expected = wrapMathML("  <apply>\n" + 
    "    <plus/>\n" + 
    "    <cn type=\"integer\"> 1 </cn>\n" + 
    "    <cn type=\"integer\"> 2 </cn>\n" + 
    "    <cn type=\"integer\"> 3 </cn>\n" + 
    "  </apply>\n");
    N = libsbml.parseFormula("1 + (2 + 3)");
    S = libsbml.writeMathMLToString(N);
    assertEquals( true, equals(expected,S) );
  }

  public void test_MathMLFormatter_plus_nary_4()
  {
    String expected = wrapMathML("  <apply>\n" + 
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
    "  </apply>\n");
    N = libsbml.parseFormula("1 + 2 + x * y * z + 3");
    S = libsbml.writeMathMLToString(N);
    assertEquals( true, equals(expected,S) );
  }

  public void test_MathMLFormatter_root()
  {
    String expected = wrapMathML("  <apply>\n" + 
    "    <root/>\n" + 
    "    <degree>\n" + 
    "      <cn type=\"integer\"> 3 </cn>\n" + 
    "    </degree>\n" + 
    "    <ci> x </ci>\n" + 
    "  </apply>\n");
    N = libsbml.parseFormula("root(3, x)");
    S = libsbml.writeMathMLToString(N);
    assertEquals( true, equals(expected,S) );
  }

  public void test_MathMLFormatter_sin()
  {
    String expected = wrapMathML("  <apply>\n" + 
    "    <sin/>\n" + 
    "    <ci> x </ci>\n" + 
    "  </apply>\n");
    N = libsbml.parseFormula("sin(x)");
    S = libsbml.writeMathMLToString(N);
    assertEquals( true, equals(expected,S) );
  }

  /**
   * Loads the SWIG-generated libSBML Java module when this class is
   * loaded, or reports a sensible diagnostic message about why it failed.
   */
  static
  {
    String varname;
    String shlibname;

    if (System.getProperty("mrj.version") != null)
    {
      varname = "DYLD_LIBRARY_PATH";    // We're on a Mac.
      shlibname = "libsbmlj.jnilib and/or libsbml.dylib";
    }
    else
    {
      varname = "LD_LIBRARY_PATH";      // We're not on a Mac.
      shlibname = "libsbmlj.so and/or libsbml.so";
    }

    try
    {
      System.loadLibrary("sbmlj");
      // For extra safety, check that the jar file is in the classpath.
      Class.forName("org.sbml.libsbml.libsbml");
    }
    catch (SecurityException e)
    {
      e.printStackTrace();
      System.err.println("Could not load the libSBML library files due to a"+
                         " security exception.\n");
      System.exit(1);
    }
    catch (UnsatisfiedLinkError e)
    {
      e.printStackTrace();
      System.err.println("Error: could not link with the libSBML library files."+
                         " It is likely\nyour " + varname +
                         " environment variable does not include the directories\n"+
                         "containing the " + shlibname + " library files.\n");
      System.exit(1);
    }
    catch (ClassNotFoundException e)
    {
      e.printStackTrace();
      System.err.println("Error: unable to load the file libsbmlj.jar."+
                         " It is likely\nyour -classpath option and CLASSPATH" +
                         " environment variable\n"+
                         "do not include the path to libsbmlj.jar.\n");
      System.exit(1);
    }
  }
}
