/*
 * @file    TestListOf.java
 * @brief   ListOf unit tests
 *
 * @author  Akiya Jouraku (Java conversion)
 * @author  Ben Bornstein 
 *
 * $Id: TestListOf.java 11546 2010-07-23 02:32:42Z mhucka $
 * $HeadURL: https://sbml.svn.sourceforge.net/svnroot/sbml/trunk/libsbml/src/bindings/java/test/org/sbml/libsbml/test/sbml/TestListOf.java $
 *
 * ====== WARNING ===== WARNING ===== WARNING ===== WARNING ===== WARNING ======
 *
 * DO NOT EDIT THIS FILE.
 *
 * This file was generated automatically by converting the file located at
 * src/sbml/test/TestListOf.c
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

public class TestListOf {

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

  public void test_ListOf_clear()
  {
    ListOf lo = new  ListOf();
    SBase sp = new  Species(2,4);
    lo.append(sp);
    lo.append(sp);
    lo.append(sp);
    lo.append(sp);
    lo.append(sp);
    assertTrue( lo.size() == 5 );
    lo.clear(true);
    assertTrue( lo.size() == 0 );
    lo.append(sp);
    lo.append(sp);
    lo.append(sp);
    lo.append(sp);
    lo.appendAndOwn(sp);
    assertTrue( lo.size() == 5 );
    SBase elem;
    elem = lo.get(0);
    elem = null;
    elem = lo.get(1);
    elem = null;
    elem = lo.get(2);
    elem = null;
    elem = lo.get(3);
    elem = null;
    elem = lo.get(4);
    elem = null;
    lo.clear(false);
    assertTrue( lo.size() == 0 );
    lo = null;
  }

  public void test_ListOf_create()
  {
    ListOf lo = new  ListOf();
    assertTrue( lo.getTypeCode() == libsbml.SBML_LIST_OF );
    assertTrue( lo.getNotes() == null );
    assertTrue( lo.getAnnotation() == null );
    assertTrue( lo.getMetaId().equals("") == true );
    assertTrue( lo.size() == 0 );
    lo = null;
  }

  public void test_ListOf_free_NULL()
  {
  }

  public void test_ListOf_get()
  {
    ListOf lo = new  ListOf();
    assertTrue( lo.size() == 0 );
    SBase sp = new  Species(2,4);
    lo.append(sp);
    assertTrue( lo.size() == 1 );
    SBase elem = lo.get(1);
    assertTrue( !sp.equals(elem) );
    sp = null;
    lo = null;
  }

  public void test_ListOf_remove()
  {
    ListOf lo = new  ListOf();
    SBase sp = new  Species(2,4);
    assertTrue( lo.size() == 0 );
    lo.append(sp);
    lo.append(sp);
    lo.append(sp);
    lo.append(sp);
    lo.append(sp);
    assertTrue( lo.size() == 5 );
    SBase elem;
    elem = lo.remove(0);
    elem = null;
    elem = lo.remove(0);
    elem = null;
    elem = lo.remove(0);
    elem = null;
    elem = lo.remove(0);
    elem = null;
    elem = lo.remove(0);
    elem = null;
    assertTrue( lo.size() == 0 );
    lo.append(sp);
    lo.append(sp);
    lo.append(sp);
    lo.append(sp);
    lo.appendAndOwn(sp);
    assertTrue( lo.size() == 5 );
    lo = null;
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
