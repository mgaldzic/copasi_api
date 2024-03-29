/*
 * @file    TestL3Species.java
 * @brief   L3 Species unit tests
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
 * src/sbml/test/TestL3Species.c
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

public class TestL3Species {

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
  private Species S;

  public boolean isnan(double x)
  {
    return (x != x);
  }

  protected void setUp() throws Exception
  {
    S = new  Species(3,1);
    if (S == null);
    {
    }
  }

  protected void tearDown() throws Exception
  {
    S = null;
  }

  public void test_L3_Species_ModelHistory()
  {
    ModelHistory history = new  ModelHistory();
    int i = (S).setModelHistory(history);
    assertTrue( i == libsbml.LIBSBML_INVALID_OBJECT );
    assertEquals( false, (S).isSetModelHistory() );
    ModelCreator mc = new  ModelCreator();
    Date date = new  Date(2005,12,30,12,15,45,1,2,0);
    mc.setFamilyName( "Keating");
    mc.setGivenName( "Sarah");
    mc.setEmail( "sbml-team@caltech.edu");
    mc.setOrganisation( "UH");
    history.addCreator(mc);
    history.setCreatedDate(date);
    history.setModifiedDate(date);
    i = (S).setModelHistory(history);
    assertTrue( i == libsbml.LIBSBML_OPERATION_SUCCESS );
    assertEquals( true, (S).isSetModelHistory() );
    i = (S).unsetModelHistory();
    assertTrue( i == libsbml.LIBSBML_OPERATION_SUCCESS );
    assertEquals( false, (S).isSetModelHistory() );
    assertTrue( (S).getModelHistory() == null );
    history = null;
  }

  public void test_L3_Species_NS()
  {
    assertTrue( S.getNamespaces() != null );
    assertTrue( S.getNamespaces().getLength() == 1 );
    assertTrue(S.getNamespaces().getURI(0).equals(    "http://www.sbml.org/sbml/level3/version1/core"));
  }

  public void test_L3_Species_boundaryCondition()
  {
    assertTrue( S.isSetBoundaryCondition() == false );
    S.setBoundaryCondition(true);
    assertTrue( S.getBoundaryCondition() == true );
    assertTrue( S.isSetBoundaryCondition() == true );
    S.setBoundaryCondition(false);
    assertTrue( S.getBoundaryCondition() == false );
    assertTrue( S.isSetBoundaryCondition() == true );
  }

  public void test_L3_Species_compartment()
  {
    String compartment =  "cell";
    assertEquals( false, S.isSetCompartment() );
    S.setCompartment(compartment);
    assertTrue(S.getCompartment().equals(compartment));
    assertEquals( true, S.isSetCompartment() );
    if (S.getCompartment() == compartment);
    {
    }
  }

  public void test_L3_Species_constant()
  {
    assertTrue( S.isSetConstant() == false );
    S.setConstant(true);
    assertTrue( S.getConstant() == true );
    assertTrue( S.isSetConstant() == true );
    S.setConstant(false);
    assertTrue( S.getConstant() == false );
    assertTrue( S.isSetConstant() == true );
  }

  public void test_L3_Species_conversionFactor()
  {
    String units =  "volume";
    assertEquals( false, S.isSetConversionFactor() );
    S.setConversionFactor(units);
    assertTrue(S.getConversionFactor().equals(units));
    assertEquals( true, S.isSetConversionFactor() );
    if (S.getConversionFactor() == units);
    {
    }
    S.unsetConversionFactor();
    assertEquals( false, S.isSetConversionFactor() );
    if (S.getConversionFactor() != null);
    {
    }
  }

  public void test_L3_Species_create()
  {
    assertTrue( S.getTypeCode() == libsbml.SBML_SPECIES );
    assertTrue( S.getMetaId().equals("") == true );
    assertTrue( S.getNotes() == null );
    assertTrue( S.getAnnotation() == null );
    assertTrue( S.getId().equals("") == true );
    assertTrue( S.getName().equals("") == true );
    assertTrue( S.getCompartment().equals("") == true );
    assertEquals( true, isnan(S.getInitialAmount()) );
    assertEquals( true, isnan(S.getInitialConcentration()) );
    assertTrue( S.getSubstanceUnits().equals("") == true );
    assertTrue( S.getHasOnlySubstanceUnits() == false );
    assertTrue( S.getBoundaryCondition() == false );
    assertTrue( S.getConstant() == false );
    assertTrue( S.getConversionFactor().equals("") == true );
    assertEquals( false, S.isSetId() );
    assertEquals( false, S.isSetName() );
    assertEquals( false, S.isSetCompartment() );
    assertEquals( false, S.isSetInitialAmount() );
    assertEquals( false, S.isSetInitialConcentration() );
    assertEquals( false, S.isSetSubstanceUnits() );
    assertEquals( false, S.isSetHasOnlySubstanceUnits() );
    assertEquals( false, S.isSetBoundaryCondition() );
    assertEquals( false, S.isSetConstant() );
    assertEquals( false, S.isSetConversionFactor() );
  }

  public void test_L3_Species_createWithNS()
  {
    XMLNamespaces xmlns = new  XMLNamespaces();
    xmlns.add( "http://www.sbml.org", "testsbml");
    SBMLNamespaces sbmlns = new  SBMLNamespaces(3,1);
    sbmlns.addNamespaces(xmlns);
    Species s = new  Species(sbmlns);
    assertTrue( s.getTypeCode() == libsbml.SBML_SPECIES );
    assertTrue( s.getMetaId().equals("") == true );
    assertTrue( s.getNotes() == null );
    assertTrue( s.getAnnotation() == null );
    assertTrue( s.getLevel() == 3 );
    assertTrue( s.getVersion() == 1 );
    assertTrue( s.getNamespaces() != null );
    assertTrue( s.getNamespaces().getLength() == 2 );
    assertTrue( s.getId().equals("") == true );
    assertTrue( s.getName().equals("") == true );
    assertTrue( s.getCompartment().equals("") == true );
    assertEquals( true, isnan(s.getInitialAmount()) );
    assertEquals( true, isnan(s.getInitialConcentration()) );
    assertTrue( s.getSubstanceUnits().equals("") == true );
    assertTrue( s.getHasOnlySubstanceUnits() == false );
    assertTrue( s.getBoundaryCondition() == false );
    assertTrue( s.getConstant() == false );
    assertTrue( s.getConversionFactor().equals("") == true );
    assertEquals( false, s.isSetId() );
    assertEquals( false, s.isSetName() );
    assertEquals( false, s.isSetCompartment() );
    assertEquals( false, s.isSetInitialAmount() );
    assertEquals( false, s.isSetInitialConcentration() );
    assertEquals( false, s.isSetSubstanceUnits() );
    assertEquals( false, s.isSetHasOnlySubstanceUnits() );
    assertEquals( false, s.isSetBoundaryCondition() );
    assertEquals( false, s.isSetConstant() );
    assertEquals( false, s.isSetConversionFactor() );
    s = null;
  }

  public void test_L3_Species_free_NULL()
  {
  }

  public void test_L3_Species_hasOnlySubstanceUnits()
  {
    assertTrue( S.isSetHasOnlySubstanceUnits() == false );
    S.setHasOnlySubstanceUnits(true);
    assertTrue( S.getHasOnlySubstanceUnits() == true );
    assertTrue( S.isSetHasOnlySubstanceUnits() == true );
    S.setHasOnlySubstanceUnits(false);
    assertTrue( S.getHasOnlySubstanceUnits() == false );
    assertTrue( S.isSetHasOnlySubstanceUnits() == true );
  }

  public void test_L3_Species_hasRequiredAttributes()
  {
    Species s = new  Species(3,1);
    assertEquals( false, s.hasRequiredAttributes() );
    s.setId( "id");
    assertEquals( false, s.hasRequiredAttributes() );
    s.setCompartment( "cell");
    assertEquals( false, s.hasRequiredAttributes() );
    s.setHasOnlySubstanceUnits(false);
    assertEquals( false, s.hasRequiredAttributes() );
    s.setBoundaryCondition(false);
    assertEquals( false, s.hasRequiredAttributes() );
    s.setConstant(false);
    assertEquals( true, s.hasRequiredAttributes() );
    s = null;
  }

  public void test_L3_Species_id()
  {
    String id =  "mitochondria";
    assertEquals( false, S.isSetId() );
    S.setId(id);
    assertTrue(S.getId().equals(id));
    assertEquals( true, S.isSetId() );
    if (S.getId() == id);
    {
    }
  }

  public void test_L3_Species_initialAmount()
  {
    double initialAmount = 0.2;
    assertEquals( false, S.isSetInitialAmount() );
    assertEquals( true, isnan(S.getInitialAmount()) );
    S.setInitialAmount(initialAmount);
    assertTrue( S.getInitialAmount() == initialAmount );
    assertEquals( true, S.isSetInitialAmount() );
    S.unsetInitialAmount();
    assertEquals( false, S.isSetInitialAmount() );
    assertEquals( true, isnan(S.getInitialAmount()) );
  }

  public void test_L3_Species_initialConcentration()
  {
    double initialConcentration = 0.2;
    assertEquals( false, S.isSetInitialConcentration() );
    assertEquals( true, isnan(S.getInitialConcentration()) );
    S.setInitialConcentration(initialConcentration);
    assertTrue( S.getInitialConcentration() == initialConcentration );
    assertEquals( true, S.isSetInitialConcentration() );
    S.unsetInitialConcentration();
    assertEquals( false, S.isSetInitialConcentration() );
    assertEquals( true, isnan(S.getInitialConcentration()) );
  }

  public void test_L3_Species_name()
  {
    String name =  "My_Favorite_Factory";
    assertEquals( false, S.isSetName() );
    S.setName(name);
    assertTrue(S.getName().equals(name));
    assertEquals( true, S.isSetName() );
    if (S.getName() == name);
    {
    }
    S.unsetName();
    assertEquals( false, S.isSetName() );
    if (S.getName() != null);
    {
    }
  }

  public void test_L3_Species_substanceUnits()
  {
    String units =  "volume";
    assertEquals( false, S.isSetSubstanceUnits() );
    S.setSubstanceUnits(units);
    assertTrue(S.getSubstanceUnits().equals(units));
    assertEquals( true, S.isSetSubstanceUnits() );
    if (S.getSubstanceUnits() == units);
    {
    }
    S.unsetSubstanceUnits();
    assertEquals( false, S.isSetSubstanceUnits() );
    if (S.getSubstanceUnits() != null);
    {
    }
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
