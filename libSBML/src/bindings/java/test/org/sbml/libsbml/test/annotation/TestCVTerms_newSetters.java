/*
 * @file    TestCVTerms_newSetters.java
 * @brief   CVTerms unit tests
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
 * src/annotation/test/TestCVTerms_newSetters.c
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

package org.sbml.libsbml.test.annotation;

import org.sbml.libsbml.*;

import java.io.File;
import java.lang.AssertionError;

public class TestCVTerms_newSetters {

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

  public void test_CVTerm_addResource()
  {
    CVTerm term = new  CVTerm(libsbml.MODEL_QUALIFIER);
    String resource =  "GO6666";
    XMLAttributes xa;
    assertTrue( term != null );
    assertTrue( term.getQualifierType() == libsbml.MODEL_QUALIFIER );
    int i = term.addResource( "");
    assertTrue( i == libsbml.LIBSBML_OPERATION_FAILED );
    xa = term.getResources();
    assertTrue( xa.getLength() == 0 );
    i = term.addResource(resource);
    assertTrue( i == libsbml.LIBSBML_OPERATION_SUCCESS );
    xa = term.getResources();
    assertTrue( xa.getLength() == 1 );
    assertTrue(xa.getName(0).equals( "rdf:resource"));
    assertTrue(xa.getValue(0).equals( "GO6666"));
    term = null;
  }

  public void test_CVTerm_removeResource()
  {
    CVTerm term = new  CVTerm(libsbml.MODEL_QUALIFIER);
    String resource =  "GO6666";
    XMLAttributes xa;
    assertTrue( term != null );
    assertTrue( term.getQualifierType() == libsbml.MODEL_QUALIFIER );
    term.addResource(resource);
    xa = term.getResources();
    assertTrue( xa.getLength() == 1 );
    int i = term.removeResource( "CCC");
    assertTrue( i == libsbml.LIBSBML_INVALID_ATTRIBUTE_VALUE );
    xa = term.getResources();
    assertTrue( xa.getLength() == 1 );
    i = term.removeResource(resource);
    assertTrue( i == libsbml.LIBSBML_OPERATION_SUCCESS );
    xa = term.getResources();
    assertTrue( xa.getLength() == 0 );
    term = null;
  }

  public void test_CVTerm_setBiolQualifierType()
  {
    CVTerm term = new  CVTerm(libsbml.BIOLOGICAL_QUALIFIER);
    assertTrue( term != null );
    assertTrue( term.getQualifierType() == libsbml.BIOLOGICAL_QUALIFIER );
    assertTrue( term.getModelQualifierType() == libsbml.BQM_UNKNOWN );
    assertTrue( term.getBiologicalQualifierType() == libsbml.BQB_UNKNOWN );
    int i = term.setBiologicalQualifierType(libsbml.BQB_IS);
    assertTrue( i == libsbml.LIBSBML_OPERATION_SUCCESS );
    assertTrue( term.getQualifierType() == libsbml.BIOLOGICAL_QUALIFIER );
    assertTrue( term.getBiologicalQualifierType() == libsbml.BQB_IS );
    assertTrue( term.getModelQualifierType() == libsbml.BQM_UNKNOWN );
    i = term.setQualifierType(libsbml.MODEL_QUALIFIER);
    assertTrue( i == libsbml.LIBSBML_OPERATION_SUCCESS );
    assertTrue( term.getQualifierType() == libsbml.MODEL_QUALIFIER );
    assertTrue( term.getModelQualifierType() == libsbml.BQM_UNKNOWN );
    assertTrue( term.getBiologicalQualifierType() == libsbml.BQB_UNKNOWN );
    i = term.setBiologicalQualifierType(libsbml.BQB_IS);
    assertTrue( i == libsbml.LIBSBML_INVALID_ATTRIBUTE_VALUE );
    assertTrue( term.getQualifierType() == libsbml.MODEL_QUALIFIER );
    assertTrue( term.getModelQualifierType() == libsbml.BQM_UNKNOWN );
    assertTrue( term.getBiologicalQualifierType() == libsbml.BQB_UNKNOWN );
    term = null;
  }

  public void test_CVTerm_setModelQualifierType()
  {
    CVTerm term = new  CVTerm(libsbml.MODEL_QUALIFIER);
    assertTrue( term != null );
    assertTrue( term.getQualifierType() == libsbml.MODEL_QUALIFIER );
    assertTrue( term.getModelQualifierType() == libsbml.BQM_UNKNOWN );
    assertTrue( term.getBiologicalQualifierType() == libsbml.BQB_UNKNOWN );
    int i = term.setModelQualifierType(libsbml.BQM_IS);
    assertTrue( i == libsbml.LIBSBML_OPERATION_SUCCESS );
    assertTrue( term.getQualifierType() == libsbml.MODEL_QUALIFIER );
    assertTrue( term.getModelQualifierType() == libsbml.BQM_IS );
    assertTrue( term.getBiologicalQualifierType() == libsbml.BQB_UNKNOWN );
    i = term.setQualifierType(libsbml.BIOLOGICAL_QUALIFIER);
    assertTrue( i == libsbml.LIBSBML_OPERATION_SUCCESS );
    assertTrue( term.getQualifierType() == libsbml.BIOLOGICAL_QUALIFIER );
    assertTrue( term.getModelQualifierType() == libsbml.BQM_UNKNOWN );
    assertTrue( term.getBiologicalQualifierType() == libsbml.BQB_UNKNOWN );
    i = term.setModelQualifierType(libsbml.BQM_IS);
    assertTrue( i == libsbml.LIBSBML_INVALID_ATTRIBUTE_VALUE );
    assertTrue( term.getQualifierType() == libsbml.BIOLOGICAL_QUALIFIER );
    assertTrue( term.getBiologicalQualifierType() == libsbml.BQB_UNKNOWN );
    assertTrue( term.getModelQualifierType() == libsbml.BQM_UNKNOWN );
    term = null;
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
