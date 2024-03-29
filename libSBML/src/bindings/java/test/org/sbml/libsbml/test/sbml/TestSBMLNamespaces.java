/*
 * @file    TestSBMLNamespaces.java
 * @brief   SBMLNamespaces unit tests
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
 * src/sbml/test/TestSBMLNamespaces.cpp
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

public class TestSBMLNamespaces {

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

  public void test_SBMLNamespaces_L1V1()
  {
    SBMLNamespaces sbml = new SBMLNamespaces(1,1);
    assertTrue( sbml.getLevel() == 1 );
    assertTrue( sbml.getVersion() == 1 );
    XMLNamespaces ns = sbml.getNamespaces();
    assertTrue( ns.getLength() == 1 );
    assertTrue( ns.getURI(0).equals( "http://www.sbml.org/sbml/level1") );
    assertTrue( ns.getPrefix(0).equals( "") );
    sbml = null;
  }

  public void test_SBMLNamespaces_L1V2()
  {
    SBMLNamespaces sbml = new SBMLNamespaces(1,2);
    assertTrue( sbml.getLevel() == 1 );
    assertTrue( sbml.getVersion() == 2 );
    XMLNamespaces ns = sbml.getNamespaces();
    assertTrue( ns.getLength() == 1 );
    assertTrue( ns.getURI(0).equals( "http://www.sbml.org/sbml/level1") );
    assertTrue( ns.getPrefix(0).equals( "") );
    sbml = null;
  }

  public void test_SBMLNamespaces_L2V1()
  {
    SBMLNamespaces sbml = new SBMLNamespaces(2,1);
    assertTrue( sbml.getLevel() == 2 );
    assertTrue( sbml.getVersion() == 1 );
    XMLNamespaces ns = sbml.getNamespaces();
    assertTrue( ns.getLength() == 1 );
    assertTrue( ns.getURI(0).equals( "http://www.sbml.org/sbml/level2") );
    assertTrue( ns.getPrefix(0).equals( "") );
    sbml = null;
  }

  public void test_SBMLNamespaces_L2V2()
  {
    SBMLNamespaces sbml = new SBMLNamespaces(2,2);
    assertTrue( sbml.getLevel() == 2 );
    assertTrue( sbml.getVersion() == 2 );
    XMLNamespaces ns = sbml.getNamespaces();
    assertTrue( ns.getLength() == 1 );
    assertTrue( ns.getURI(0).equals( "http://www.sbml.org/sbml/level2/version2") );
    assertTrue( ns.getPrefix(0).equals( "") );
    sbml = null;
  }

  public void test_SBMLNamespaces_L2V3()
  {
    SBMLNamespaces sbml = new SBMLNamespaces(2,3);
    assertTrue( sbml.getLevel() == 2 );
    assertTrue( sbml.getVersion() == 3 );
    XMLNamespaces ns = sbml.getNamespaces();
    assertTrue( ns.getLength() == 1 );
    assertTrue( ns.getURI(0).equals( "http://www.sbml.org/sbml/level2/version3") );
    assertTrue( ns.getPrefix(0).equals( "") );
    sbml = null;
  }

  public void test_SBMLNamespaces_L2V4()
  {
    SBMLNamespaces sbml = new SBMLNamespaces(2,4);
    assertTrue( sbml.getLevel() == 2 );
    assertTrue( sbml.getVersion() == 4 );
    XMLNamespaces ns = sbml.getNamespaces();
    assertTrue( ns.getLength() == 1 );
    assertTrue( ns.getURI(0).equals( "http://www.sbml.org/sbml/level2/version4") );
    assertTrue( ns.getPrefix(0).equals( "") );
    sbml = null;
  }

  public void test_SBMLNamespaces_getURI()
  {
    assertTrue( SBMLNamespaces.getSBMLNamespaceURI(1,1).equals(                             "http://www.sbml.org/sbml/level1") );
    assertTrue( SBMLNamespaces.getSBMLNamespaceURI(1,2).equals(                             "http://www.sbml.org/sbml/level1") );
    assertTrue( SBMLNamespaces.getSBMLNamespaceURI(2,1).equals(                             "http://www.sbml.org/sbml/level2") );
    assertTrue( SBMLNamespaces.getSBMLNamespaceURI(2,2).equals(                             "http://www.sbml.org/sbml/level2/version2") );
    assertTrue( SBMLNamespaces.getSBMLNamespaceURI(2,3).equals(                             "http://www.sbml.org/sbml/level2/version3") );
    assertTrue( SBMLNamespaces.getSBMLNamespaceURI(2,4).equals(                             "http://www.sbml.org/sbml/level2/version4") );
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
