/*
 * @file    TestUnitDefinition_newSetters.java
 * @brief   SBML UnitDefinition unit tests for new API
 *
 * @author  Akiya Jouraku (Java conversion)
 * @author  sarah Keating 
 *
 * $Id$
 * $HeadURL$
 *
 * ====== WARNING ===== WARNING ===== WARNING ===== WARNING ===== WARNING ======
 *
 * DO NOT EDIT THIS FILE.
 *
 * This file was generated automatically by converting the file located at
 * src/sbml/test/TestUnitDefinition_newSetters.c
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

public class TestUnitDefinition_newSetters {

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
  private UnitDefinition UD;

  protected void setUp() throws Exception
  {
    UD = new  UnitDefinition(2,4);
    if (UD == null);
    {
    }
  }

  protected void tearDown() throws Exception
  {
    UD = null;
  }

  public void test_UnitDefinition_addUnit1()
  {
    UnitDefinition m = new  UnitDefinition(2,2);
    Unit p = new  Unit(2,2);
    int i = m.addUnit(p);
    assertTrue( i == libsbml.LIBSBML_INVALID_OBJECT );
    p.setKind(libsbml.UNIT_KIND_MOLE);
    i = m.addUnit(p);
    assertTrue( i == libsbml.LIBSBML_OPERATION_SUCCESS );
    assertTrue( m.getNumUnits() == 1 );
    p = null;
    m = null;
  }

  public void test_UnitDefinition_addUnit2()
  {
    UnitDefinition m = new  UnitDefinition(2,2);
    Unit p = new  Unit(2,1);
    p.setKind(libsbml.UNIT_KIND_MOLE);
    int i = m.addUnit(p);
    assertTrue( i == libsbml.LIBSBML_VERSION_MISMATCH );
    assertTrue( m.getNumUnits() == 0 );
    p = null;
    m = null;
  }

  public void test_UnitDefinition_addUnit3()
  {
    UnitDefinition m = new  UnitDefinition(2,2);
    Unit p = new  Unit(1,2);
    p.setKind(libsbml.UNIT_KIND_MOLE);
    int i = m.addUnit(p);
    assertTrue( i == libsbml.LIBSBML_LEVEL_MISMATCH );
    assertTrue( m.getNumUnits() == 0 );
    p = null;
    m = null;
  }

  public void test_UnitDefinition_addUnit4()
  {
    UnitDefinition m = new  UnitDefinition(2,2);
    Unit p = null;
    int i = m.addUnit(p);
    assertTrue( i == libsbml.LIBSBML_OPERATION_FAILED );
    assertTrue( m.getNumUnits() == 0 );
    m = null;
  }

  public void test_UnitDefinition_createUnit()
  {
    UnitDefinition m = new  UnitDefinition(2,2);
    Unit p = m.createUnit();
    assertTrue( m.getNumUnits() == 1 );
    assertTrue( (p).getLevel() == 2 );
    assertTrue( (p).getVersion() == 2 );
    m = null;
  }

  public void test_UnitDefinition_setId1()
  {
    int i = UD.setId( "mmls");
    assertTrue( i == libsbml.LIBSBML_OPERATION_SUCCESS );
    assertTrue(UD.getId().equals( "mmls"));
    assertEquals( true, UD.isSetId() );
    i = UD.setId("");
    assertTrue( i == libsbml.LIBSBML_OPERATION_SUCCESS );
    assertEquals( false, UD.isSetId() );
    i = UD.setId( "123");
    assertTrue( i == libsbml.LIBSBML_INVALID_ATTRIBUTE_VALUE );
    assertEquals( false, UD.isSetId() );
  }

  public void test_UnitDefinition_setName1()
  {
    int i = UD.setName( "mmls");
    assertTrue( i == libsbml.LIBSBML_OPERATION_SUCCESS );
    assertTrue(UD.getName().equals( "mmls"));
    assertEquals( true, UD.isSetName() );
    i = UD.setName("");
    assertTrue( i == libsbml.LIBSBML_OPERATION_SUCCESS );
    assertEquals( false, UD.isSetName() );
    i = UD.setName( "123");
    assertTrue( i == libsbml.LIBSBML_OPERATION_SUCCESS );
    assertEquals( true, UD.isSetName() );
  }

  public void test_UnitDefinition_setName2()
  {
    int i = UD.setName( "mmls");
    assertTrue( i == libsbml.LIBSBML_OPERATION_SUCCESS );
    assertTrue(UD.getName().equals( "mmls"));
    assertEquals( true, UD.isSetName() );
    i = UD.unsetName();
    assertTrue( i == libsbml.LIBSBML_OPERATION_SUCCESS );
    assertEquals( false, UD.isSetName() );
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
