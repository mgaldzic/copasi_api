/*
 * @file    TestFunctionDefinition_newSetters.java
 * @brief   FunctionDefinition unit tests for new set function API
 *
 * @author  Akiya Jouraku (Java conversion)
 * @author  Sarah Keating 
 *
 * $Id$
 * $HeadURL$
 *
 * ====== WARNING ===== WARNING ===== WARNING ===== WARNING ===== WARNING ======
 *
 * DO NOT EDIT THIS FILE.
 *
 * This file was generated automatically by converting the file located at
 * src/sbml/test/TestFunctionDefinition_newSetters.c
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

package org.sbml.libsbml.test.sbml;

import org.sbml.libsbml.*;

import java.io.File;
import java.lang.AssertionError;

public class TestFunctionDefinition_newSetters {

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
  private FunctionDefinition E;

  protected void setUp() throws Exception
  {
    E = new  FunctionDefinition(2,4);
    if (E == null);
    {
    }
  }

  protected void tearDown() throws Exception
  {
    E = null;
  }

  public void test_FunctionDefinition_setId1()
  {
    String id =  "1e1";
    int i = E.setId(id);
    assertTrue( i == libsbml.LIBSBML_INVALID_ATTRIBUTE_VALUE );
    assertEquals( false, E.isSetId() );
  }

  public void test_FunctionDefinition_setId2()
  {
    String id =  "e1";
    int i = E.setId(id);
    assertTrue( i == libsbml.LIBSBML_OPERATION_SUCCESS );
    assertTrue(E.getId().equals(id));
    assertEquals( true, E.isSetId() );
    i = E.setId("");
    assertTrue( i == libsbml.LIBSBML_OPERATION_SUCCESS );
    assertEquals( false, E.isSetId() );
  }

  public void test_FunctionDefinition_setMath1()
  {
    ASTNode math = libsbml.parseFormula("2 * k");
    int i = E.setMath(math);
    assertTrue( i == libsbml.LIBSBML_OPERATION_SUCCESS );
    assertTrue( !E.getMath().equals(math) );
    assertEquals( true, E.isSetMath() );
    i = E.setMath(null);
    assertTrue( i == libsbml.LIBSBML_OPERATION_SUCCESS );
    assertTrue( E.getMath() == null );
    assertEquals( false, E.isSetMath() );
    math = null;
  }

  public void test_FunctionDefinition_setMath2()
  {
    ASTNode math = new  ASTNode(libsbml.AST_TIMES);
    int i = E.setMath(math);
    assertTrue( i == libsbml.LIBSBML_INVALID_OBJECT );
    assertEquals( false, E.isSetMath() );
    math = null;
  }

  public void test_FunctionDefinition_setName1()
  {
    String name =  "3Set_k2";
    int i = E.setName(name);
    assertTrue( i == libsbml.LIBSBML_OPERATION_SUCCESS );
    assertEquals( true, E.isSetName() );
  }

  public void test_FunctionDefinition_setName2()
  {
    String name =  "Set k2";
    int i = E.setName(name);
    assertTrue( i == libsbml.LIBSBML_OPERATION_SUCCESS );
    assertTrue(E.getName().equals(name));
    assertEquals( true, E.isSetName() );
    i = E.unsetName();
    assertTrue( i == libsbml.LIBSBML_OPERATION_SUCCESS );
    assertEquals( false, E.isSetName() );
  }

  public void test_FunctionDefinition_setName3()
  {
    int i = E.setName("");
    assertTrue( i == libsbml.LIBSBML_OPERATION_SUCCESS );
    assertEquals( false, E.isSetName() );
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
