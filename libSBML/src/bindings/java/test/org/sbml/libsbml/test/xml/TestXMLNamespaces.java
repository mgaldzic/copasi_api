/*
 * @file    TestXMLNamespaces.java
 * @brief   XMLNamespaces unit tests
 *
 * @author  Akiya Jouraku (Java conversion)
 * @author  Michael Hucka <mhucka@caltech.edu> 
 *
 * $Id$
 * $HeadURL$
 *
 * ====== WARNING ===== WARNING ===== WARNING ===== WARNING ===== WARNING ======
 *
 * DO NOT EDIT THIS FILE.
 *
 * This file was generated automatically by converting the file located at
 * src/xml/test/TestXMLNamespaces.c
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

package org.sbml.libsbml.test.xml;

import org.sbml.libsbml.*;

import java.io.File;
import java.lang.AssertionError;

public class TestXMLNamespaces {

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
  private XMLNamespaces NS;

  protected void setUp() throws Exception
  {
    NS = new  XMLNamespaces();
    if (NS == null);
    {
    }
  }

  protected void tearDown() throws Exception
  {
    NS = null;
  }

  public void test_XMLNamespaces_add()
  {
    assertTrue( NS.getLength() == 0 );
    assertTrue( NS.isEmpty() == true );
    NS.add( "http://test1.org/", "test1");
    assertTrue( NS.getLength() == 1 );
    assertTrue( NS.isEmpty() == false );
    NS.add( "http://test2.org/", "test2");
    assertTrue( NS.getLength() == 2 );
    assertTrue( NS.isEmpty() == false );
    NS.add( "http://test1.org/", "test1a");
    assertTrue( NS.getLength() == 3 );
    assertTrue( NS.isEmpty() == false );
    NS.add( "http://test1.org/", "test1a");
    assertTrue( NS.getLength() == 3 );
    assertTrue( NS.isEmpty() == false );
    assertTrue( ! (NS.getIndex( "http://test1.org/") == -1) );
  }

  public void test_XMLNamespaces_add1()
  {
    assertTrue( NS.getLength() == 0 );
    assertTrue( NS.isEmpty() == true );
    int i = NS.add( "http://test1.org/", "test1");
    assertTrue( i == libsbml.LIBSBML_OPERATION_SUCCESS );
    assertTrue( NS.getLength() == 1 );
    assertTrue( NS.isEmpty() == false );
  }

  public void test_XMLNamespaces_baseline()
  {
    assertTrue( NS.getLength() == 0 );
    assertTrue( NS.isEmpty() == true );
  }

  public void test_XMLNamespaces_clear()
  {
    NS.add( "http://test1.org/", "test1");
    NS.add( "http://test2.org/", "test2");
    NS.add( "http://test3.org/", "test3");
    NS.add( "http://test4.org/", "test4");
    NS.add( "http://test5.org/", "test5");
    assertTrue( NS.getLength() == 5 );
    int i = NS.clear();
    assertTrue( i == libsbml.LIBSBML_OPERATION_SUCCESS );
    assertTrue( NS.getLength() == 0 );
  }

  public void test_XMLNamespaces_get()
  {
    NS.add( "http://test1.org/", "test1");
    NS.add( "http://test2.org/", "test2");
    NS.add( "http://test3.org/", "test3");
    NS.add( "http://test4.org/", "test4");
    NS.add( "http://test5.org/", "test5");
    NS.add( "http://test6.org/", "test6");
    NS.add( "http://test7.org/", "test7");
    NS.add( "http://test8.org/", "test8");
    NS.add( "http://test9.org/", "test9");
    assertTrue( NS.getLength() == 9 );
    assertTrue( NS.getIndex( "http://test1.org/") == 0 );
    assertTrue( !NS.getPrefix(1).equals( "test2") == false );
    assertTrue( !NS.getPrefix( "http://test1.org/").equals(		      "test1") == false );
    assertTrue( !NS.getURI(1).equals( "http://test2.org/") == false );
    assertTrue( !NS.getURI( "test2").equals(		      "http://test2.org/") == false );
    assertTrue( NS.getIndex( "http://test1.org/") == 0 );
    assertTrue( NS.getIndex( "http://test2.org/") == 1 );
    assertTrue( NS.getIndex( "http://test5.org/") == 4 );
    assertTrue( NS.getIndex( "http://test9.org/") == 8 );
    assertTrue( NS.getIndex( "http://testX.org/") == -1 );
    assertTrue( NS.hasURI( "http://test1.org/") != false );
    assertTrue( NS.hasURI( "http://test2.org/") != false );
    assertTrue( NS.hasURI( "http://test5.org/") != false );
    assertTrue( NS.hasURI( "http://test9.org/") != false );
    assertTrue( NS.hasURI( "http://testX.org/") == false );
    assertTrue( NS.getIndexByPrefix( "test1") == 0 );
    assertTrue( NS.getIndexByPrefix( "test5") == 4 );
    assertTrue( NS.getIndexByPrefix( "test9") == 8 );
    assertTrue( NS.getIndexByPrefix( "testX") == -1 );
    assertTrue( NS.hasPrefix( "test1") != false );
    assertTrue( NS.hasPrefix( "test5") != false );
    assertTrue( NS.hasPrefix( "test9") != false );
    assertTrue( NS.hasPrefix( "testX") == false );
    assertTrue( NS.hasNS( "http://test1.org/", "test1") != false );
    assertTrue( NS.hasNS( "http://test5.org/", "test5") != false );
    assertTrue( NS.hasNS( "http://test9.org/", "test9") != false );
    assertTrue( NS.hasNS( "http://testX.org/", "testX") == false );
  }

  public void test_XMLNamespaces_remove()
  {
    NS.add( "http://test1.org/", "test1");
    NS.add( "http://test2.org/", "test2");
    NS.add( "http://test3.org/", "test3");
    NS.add( "http://test4.org/", "test4");
    NS.add( "http://test5.org/", "test5");
    assertTrue( NS.getLength() == 5 );
    NS.remove(4);
    assertTrue( NS.getLength() == 4 );
    NS.remove(3);
    assertTrue( NS.getLength() == 3 );
    NS.remove(2);
    assertTrue( NS.getLength() == 2 );
    NS.remove(1);
    assertTrue( NS.getLength() == 1 );
    NS.remove(0);
    assertTrue( NS.getLength() == 0 );
    NS.add( "http://test1.org/", "test1");
    NS.add( "http://test2.org/", "test2");
    NS.add( "http://test3.org/", "test3");
    NS.add( "http://test4.org/", "test4");
    NS.add( "http://test5.org/", "test5");
    assertTrue( NS.getLength() == 5 );
    NS.remove(0);
    assertTrue( NS.getLength() == 4 );
    NS.remove(0);
    assertTrue( NS.getLength() == 3 );
    NS.remove(0);
    assertTrue( NS.getLength() == 2 );
    NS.remove(0);
    assertTrue( NS.getLength() == 1 );
    NS.remove(0);
    assertTrue( NS.getLength() == 0 );
  }

  public void test_XMLNamespaces_remove1()
  {
    NS.add( "http://test1.org/", "test1");
    NS.add( "http://test2.org/", "test2");
    assertTrue( NS.getLength() == 2 );
    int i = NS.remove(4);
    assertTrue( i == libsbml.LIBSBML_INDEX_EXCEEDS_SIZE );
    assertTrue( NS.getLength() == 2 );
    i = NS.remove( "test4");
    assertTrue( i == libsbml.LIBSBML_INDEX_EXCEEDS_SIZE );
    assertTrue( NS.getLength() == 2 );
    i = NS.remove(1);
    assertTrue( i == libsbml.LIBSBML_OPERATION_SUCCESS );
    assertTrue( NS.getLength() == 1 );
    i = NS.remove( "test1");
    assertTrue( i == libsbml.LIBSBML_OPERATION_SUCCESS );
    assertTrue( NS.getLength() == 0 );
  }

  public void test_XMLNamespaces_remove_by_prefix()
  {
    NS.add( "http://test1.org/", "test1");
    NS.add( "http://test2.org/", "test2");
    NS.add( "http://test3.org/", "test3");
    NS.add( "http://test4.org/", "test4");
    NS.add( "http://test5.org/", "test5");
    assertTrue( NS.getLength() == 5 );
    NS.remove( "test1");
    assertTrue( NS.getLength() == 4 );
    NS.remove( "test2");
    assertTrue( NS.getLength() == 3 );
    NS.remove( "test3");
    assertTrue( NS.getLength() == 2 );
    NS.remove( "test4");
    assertTrue( NS.getLength() == 1 );
    NS.remove( "test5");
    assertTrue( NS.getLength() == 0 );
    NS.add( "http://test1.org/", "test1");
    NS.add( "http://test2.org/", "test2");
    NS.add( "http://test3.org/", "test3");
    NS.add( "http://test4.org/", "test4");
    NS.add( "http://test5.org/", "test5");
    assertTrue( NS.getLength() == 5 );
    NS.remove( "test5");
    assertTrue( NS.getLength() == 4 );
    NS.remove( "test4");
    assertTrue( NS.getLength() == 3 );
    NS.remove( "test3");
    assertTrue( NS.getLength() == 2 );
    NS.remove( "test2");
    assertTrue( NS.getLength() == 1 );
    NS.remove( "test1");
    assertTrue( NS.getLength() == 0 );
    NS.add( "http://test1.org/", "test1");
    NS.add( "http://test2.org/", "test2");
    NS.add( "http://test3.org/", "test3");
    NS.add( "http://test4.org/", "test4");
    NS.add( "http://test5.org/", "test5");
    assertTrue( NS.getLength() == 5 );
    NS.remove( "test3");
    assertTrue( NS.getLength() == 4 );
    NS.remove( "test1");
    assertTrue( NS.getLength() == 3 );
    NS.remove( "test4");
    assertTrue( NS.getLength() == 2 );
    NS.remove( "test5");
    assertTrue( NS.getLength() == 1 );
    NS.remove( "test2");
    assertTrue( NS.getLength() == 0 );
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
