///  @file    TestModel_newSetters.cs
///  @brief   Model unit tests for new set function API
///  @author  Frank Bergmann (Csharp conversion)
///  @author  Akiya Jouraku (Csharp conversion)
///  @author  Sarah Keating 
/// 
///  $Id: TestModel_newSetters.cs 11545 2010-07-23 02:19:10Z mhucka $
///  $HeadURL: https://sbml.svn.sourceforge.net/svnroot/sbml/trunk/libsbml/src/bindings/csharp/test/sbml/TestModel_newSetters.cs $
/// 
///  ====== WARNING ===== WARNING ===== WARNING ===== WARNING ===== WARNING ======
/// 
///  DO NOT EDIT THIS FILE.
/// 
///  This file was generated automatically by converting the file located at
///  src/sbml/test/TestModel_newSetters.c
///  using the conversion program dev/utilities/translateTests/translateTests.pl.
///  Any changes made here will be lost the next time the file is regenerated.
/// 
///  -----------------------------------------------------------------------------
///  This file is part of libSBML.  Please visit http://sbml.org for more
///  information about SBML, and the latest version of libSBML.
/// 
///  Copyright 2005-2010 California Institute of Technology.
///  Copyright 2002-2005 California Institute of Technology and
///                      Japan Science and Technology Corporation.
///  
///  This library is free software; you can redistribute it and/or modify it
///  under the terms of the GNU Lesser General Public License as published by
///  the Free Software Foundation.  A copy of the license agreement is provided
///  in the file named "LICENSE.txt" included with this software distribution
///  and also available online as http://sbml.org/software/libsbml/license.html
///  -----------------------------------------------------------------------------


namespace LibSBMLCSTest {

  using libsbml;

  using System;

  using System.IO;

  public class TestModel_newSetters {
    public class AssertionError : System.Exception 
    {
      public AssertionError() : base()
      {
        
      }
    }


    static void assertTrue(bool condition)
    {
      if (condition == true)
      {
        return;
      }
      throw new AssertionError();
    }

    static void assertEquals(object a, object b)
    {
      if ( (a == null) && (b == null) )
      {
        return;
      }
      else if ( (a == null) || (b == null) )
      {
        throw new AssertionError();
      }
      else if (a.Equals(b))
      {
        return;
      }
  
      throw new AssertionError();
    }

    static void assertNotEquals(object a, object b)
    {
      if ( (a == null) && (b == null) )
      {
        throw new AssertionError();
      }
      else if ( (a == null) || (b == null) )
      {
        return;
      }
      else if (a.Equals(b))
      {
        throw new AssertionError();
      }
    }

    static void assertEquals(bool a, bool b)
    {
      if ( a == b )
      {
        return;
      }
      throw new AssertionError();
    }

    static void assertNotEquals(bool a, bool b)
    {
      if ( a != b )
      {
        return;
      }
      throw new AssertionError();
    }

    static void assertEquals(int a, int b)
    {
      if ( a == b )
      {
        return;
      }
      throw new AssertionError();
    }

    static void assertNotEquals(int a, int b)
    {
      if ( a != b )
      {
        return;
      }
      throw new AssertionError();
    }

    private Model M;

    public void setUp()
    {
      M = new  Model(2,4);
      if (M == null);
      {
      }
    }

    public void tearDown()
    {
      M = null;
    }

    public void test_Model_addCompartment1()
    {
      Model m = new  Model(2,2);
      Compartment c = new  Compartment(2,2);
      int i = m.addCompartment(c);
      assertTrue( i == libsbml.LIBSBML_INVALID_OBJECT );
      c.setId( "c");
      i = m.addCompartment(c);
      assertTrue( i == libsbml.LIBSBML_OPERATION_SUCCESS );
      assertTrue( m.getNumCompartments() == 1 );
      c = null;
      m = null;
    }

    public void test_Model_addCompartment2()
    {
      Model m = new  Model(2,2);
      Compartment c = new  Compartment(2,1);
      c.setId( "c");
      int i = m.addCompartment(c);
      assertTrue( i == libsbml.LIBSBML_VERSION_MISMATCH );
      assertTrue( m.getNumCompartments() == 0 );
      c = null;
      m = null;
    }

    public void test_Model_addCompartment3()
    {
      Model m = new  Model(2,2);
      Compartment c = new  Compartment(1,2);
      c.setId( "c");
      int i = m.addCompartment(c);
      assertTrue( i == libsbml.LIBSBML_LEVEL_MISMATCH );
      assertTrue( m.getNumCompartments() == 0 );
      c = null;
      m = null;
    }

    public void test_Model_addCompartment4()
    {
      Model m = new  Model(2,2);
      Compartment c = null;
      int i = m.addCompartment(c);
      assertTrue( i == libsbml.LIBSBML_OPERATION_FAILED );
      assertTrue( m.getNumCompartments() == 0 );
      m = null;
    }

    public void test_Model_addCompartment5()
    {
      Model m = new  Model(2,2);
      Compartment c = new  Compartment(2,2);
      c.setId( "c");
      Compartment c1 = new  Compartment(2,2);
      c1.setId( "c");
      int i = m.addCompartment(c);
      assertTrue( i == libsbml.LIBSBML_OPERATION_SUCCESS );
      assertTrue( m.getNumCompartments() == 1 );
      i = m.addCompartment(c1);
      assertTrue( i == libsbml.LIBSBML_DUPLICATE_OBJECT_ID );
      assertTrue( m.getNumCompartments() == 1 );
      c = null;
      c1 = null;
      m = null;
    }

    public void test_Model_addCompartmentType1()
    {
      Model m = new  Model(2,2);
      CompartmentType ct = new  CompartmentType(2,2);
      int i = m.addCompartmentType(ct);
      assertTrue( i == libsbml.LIBSBML_INVALID_OBJECT );
      ct.setId( "ct");
      i = m.addCompartmentType(ct);
      assertTrue( i == libsbml.LIBSBML_OPERATION_SUCCESS );
      assertTrue( m.getNumCompartmentTypes() == 1 );
      ct = null;
      m = null;
    }

    public void test_Model_addCompartmentType2()
    {
      Model m = new  Model(2,2);
      CompartmentType ct = new  CompartmentType(2,3);
      ct.setId( "ct");
      int i = m.addCompartmentType(ct);
      assertTrue( i == libsbml.LIBSBML_VERSION_MISMATCH );
      assertTrue( m.getNumCompartmentTypes() == 0 );
      ct = null;
      m = null;
    }

    public void test_Model_addCompartmentType3()
    {
      Model m = new  Model(2,2);
      CompartmentType ct = null;
      int i = m.addCompartmentType(ct);
      assertTrue( i == libsbml.LIBSBML_OPERATION_FAILED );
      assertTrue( m.getNumCompartmentTypes() == 0 );
      m = null;
    }

    public void test_Model_addCompartmentType4()
    {
      Model m = new  Model(2,2);
      CompartmentType ct = new  CompartmentType(2,2);
      ct.setId( "ct");
      CompartmentType ct1 = new  CompartmentType(2,2);
      ct1.setId( "ct");
      int i = m.addCompartmentType(ct);
      assertTrue( i == libsbml.LIBSBML_OPERATION_SUCCESS );
      assertTrue( m.getNumCompartmentTypes() == 1 );
      i = m.addCompartmentType(ct1);
      assertTrue( i == libsbml.LIBSBML_DUPLICATE_OBJECT_ID );
      assertTrue( m.getNumCompartmentTypes() == 1 );
      ct = null;
      ct1 = null;
      m = null;
    }

    public void test_Model_addConstraint1()
    {
      Model m = new  Model(2,2);
      Constraint c = new  Constraint(2,2);
      int i = m.addConstraint(c);
      assertTrue( i == libsbml.LIBSBML_INVALID_OBJECT );
      c.setMath(libsbml.parseFormula("a+b"));
      i = m.addConstraint(c);
      assertTrue( i == libsbml.LIBSBML_OPERATION_SUCCESS );
      assertTrue( m.getNumConstraints() == 1 );
      c = null;
      m = null;
    }

    public void test_Model_addConstraint2()
    {
      Model m = new  Model(2,2);
      Constraint c = new  Constraint(2,3);
      c.setMath(libsbml.parseFormula("a+b"));
      int i = m.addConstraint(c);
      assertTrue( i == libsbml.LIBSBML_VERSION_MISMATCH );
      assertTrue( m.getNumConstraints() == 0 );
      c = null;
      m = null;
    }

    public void test_Model_addConstraint3()
    {
      Model m = new  Model(2,2);
      Constraint c = null;
      int i = m.addConstraint(c);
      assertTrue( i == libsbml.LIBSBML_OPERATION_FAILED );
      assertTrue( m.getNumConstraints() == 0 );
      m = null;
    }

    public void test_Model_addEvent1()
    {
      Model m = new  Model(2,2);
      Event e = new  Event(2,2);
      Trigger t = new  Trigger(2,2);
      int i = m.addEvent(e);
      assertTrue( i == libsbml.LIBSBML_INVALID_OBJECT );
      e.setTrigger(t);
      i = m.addEvent(e);
      assertTrue( i == libsbml.LIBSBML_INVALID_OBJECT );
      e.createEventAssignment();
      i = m.addEvent(e);
      assertTrue( i == libsbml.LIBSBML_OPERATION_SUCCESS );
      assertTrue( m.getNumEvents() == 1 );
      e = null;
      m = null;
    }

    public void test_Model_addEvent2()
    {
      Model m = new  Model(2,2);
      Event e = new  Event(2,1);
      Trigger t = new  Trigger(2,1);
      e.setTrigger(t);
      e.createEventAssignment();
      int i = m.addEvent(e);
      assertTrue( i == libsbml.LIBSBML_VERSION_MISMATCH );
      assertTrue( m.getNumEvents() == 0 );
      e = null;
      m = null;
    }

    public void test_Model_addEvent3()
    {
      Model m = new  Model(2,2);
      Event e = null;
      int i = m.addEvent(e);
      assertTrue( i == libsbml.LIBSBML_OPERATION_FAILED );
      assertTrue( m.getNumEvents() == 0 );
      m = null;
    }

    public void test_Model_addEvent4()
    {
      Model m = new  Model(2,2);
      Event e = new  Event(2,2);
      Trigger t = new  Trigger(2,2);
      e.setId( "e");
      e.setTrigger(t);
      e.createEventAssignment();
      Event e1 = new  Event(2,2);
      e1.setId( "e");
      e1.setTrigger(t);
      e1.createEventAssignment();
      int i = m.addEvent(e);
      assertTrue( i == libsbml.LIBSBML_OPERATION_SUCCESS );
      assertTrue( m.getNumEvents() == 1 );
      i = m.addEvent(e1);
      assertTrue( i == libsbml.LIBSBML_DUPLICATE_OBJECT_ID );
      assertTrue( m.getNumEvents() == 1 );
      e = null;
      e1 = null;
      m = null;
    }

    public void test_Model_addFunctionDefinition1()
    {
      Model m = new  Model(2,2);
      FunctionDefinition fd = new  FunctionDefinition(2,2);
      int i = m.addFunctionDefinition(fd);
      assertTrue( i == libsbml.LIBSBML_INVALID_OBJECT );
      fd.setId( "fd");
      i = m.addFunctionDefinition(fd);
      assertTrue( i == libsbml.LIBSBML_INVALID_OBJECT );
      fd.setMath(libsbml.parseFormula("fd"));
      i = m.addFunctionDefinition(fd);
      assertTrue( i == libsbml.LIBSBML_OPERATION_SUCCESS );
      assertTrue( m.getNumFunctionDefinitions() == 1 );
      fd = null;
      m = null;
    }

    public void test_Model_addFunctionDefinition2()
    {
      Model m = new  Model(2,2);
      FunctionDefinition fd = new  FunctionDefinition(2,1);
      fd.setId( "fd");
      fd.setMath(libsbml.parseFormula("fd"));
      int i = m.addFunctionDefinition(fd);
      assertTrue( i == libsbml.LIBSBML_VERSION_MISMATCH );
      assertTrue( m.getNumFunctionDefinitions() == 0 );
      fd = null;
      m = null;
    }

    public void test_Model_addFunctionDefinition3()
    {
      Model m = new  Model(2,2);
      FunctionDefinition fd = null;
      int i = m.addFunctionDefinition(fd);
      assertTrue( i == libsbml.LIBSBML_OPERATION_FAILED );
      assertTrue( m.getNumFunctionDefinitions() == 0 );
      m = null;
    }

    public void test_Model_addFunctionDefinition4()
    {
      Model m = new  Model(2,2);
      FunctionDefinition fd = new  FunctionDefinition(2,2);
      fd.setId( "fd");
      fd.setMath(libsbml.parseFormula("fd"));
      FunctionDefinition fd1 = new  FunctionDefinition(2,2);
      fd1.setId( "fd");
      fd1.setMath(libsbml.parseFormula("fd"));
      int i = m.addFunctionDefinition(fd);
      assertTrue( i == libsbml.LIBSBML_OPERATION_SUCCESS );
      assertTrue( m.getNumFunctionDefinitions() == 1 );
      i = m.addFunctionDefinition(fd1);
      assertTrue( i == libsbml.LIBSBML_DUPLICATE_OBJECT_ID );
      assertTrue( m.getNumFunctionDefinitions() == 1 );
      fd = null;
      fd1 = null;
      m = null;
    }

    public void test_Model_addInitialAssignment1()
    {
      Model m = new  Model(2,2);
      InitialAssignment ia = new  InitialAssignment(2,2);
      int i = m.addInitialAssignment(ia);
      assertTrue( i == libsbml.LIBSBML_INVALID_OBJECT );
      ia.setSymbol( "i");
      i = m.addInitialAssignment(ia);
      assertTrue( i == libsbml.LIBSBML_INVALID_OBJECT );
      ia.setMath(libsbml.parseFormula("gg"));
      i = m.addInitialAssignment(ia);
      assertTrue( i == libsbml.LIBSBML_OPERATION_SUCCESS );
      assertTrue( m.getNumInitialAssignments() == 1 );
      ia = null;
      m = null;
    }

    public void test_Model_addInitialAssignment2()
    {
      Model m = new  Model(2,2);
      InitialAssignment ia = new  InitialAssignment(2,3);
      ia.setSymbol( "i");
      ia.setMath(libsbml.parseFormula("gg"));
      int i = m.addInitialAssignment(ia);
      assertTrue( i == libsbml.LIBSBML_VERSION_MISMATCH );
      assertTrue( m.getNumInitialAssignments() == 0 );
      ia = null;
      m = null;
    }

    public void test_Model_addInitialAssignment3()
    {
      Model m = new  Model(2,2);
      InitialAssignment ia = null;
      int i = m.addInitialAssignment(ia);
      assertTrue( i == libsbml.LIBSBML_OPERATION_FAILED );
      assertTrue( m.getNumInitialAssignments() == 0 );
      m = null;
    }

    public void test_Model_addInitialAssignment4()
    {
      Model m = new  Model(2,2);
      InitialAssignment ia = new  InitialAssignment(2,2);
      ia.setSymbol( "ia");
      ia.setMath(libsbml.parseFormula("a+b"));
      InitialAssignment ia1 = new  InitialAssignment(2,2);
      ia1.setSymbol( "ia");
      ia1.setMath(libsbml.parseFormula("a+b"));
      int i = m.addInitialAssignment(ia);
      assertTrue( i == libsbml.LIBSBML_OPERATION_SUCCESS );
      assertTrue( m.getNumInitialAssignments() == 1 );
      i = m.addInitialAssignment(ia1);
      assertTrue( i == libsbml.LIBSBML_DUPLICATE_OBJECT_ID );
      assertTrue( m.getNumInitialAssignments() == 1 );
      ia = null;
      ia1 = null;
      m = null;
    }

    public void test_Model_addParameter1()
    {
      Model m = new  Model(2,2);
      Parameter p = new  Parameter(2,2);
      int i = m.addParameter(p);
      assertTrue( i == libsbml.LIBSBML_INVALID_OBJECT );
      p.setId( "p");
      i = m.addParameter(p);
      assertTrue( i == libsbml.LIBSBML_OPERATION_SUCCESS );
      assertTrue( m.getNumParameters() == 1 );
      p = null;
      m = null;
    }

    public void test_Model_addParameter2()
    {
      Model m = new  Model(2,2);
      Parameter p = new  Parameter(2,1);
      p.setId( "p");
      int i = m.addParameter(p);
      assertTrue( i == libsbml.LIBSBML_VERSION_MISMATCH );
      assertTrue( m.getNumParameters() == 0 );
      p = null;
      m = null;
    }

    public void test_Model_addParameter3()
    {
      Model m = new  Model(2,2);
      Parameter p = new  Parameter(1,2);
      p.setId( "p");
      int i = m.addParameter(p);
      assertTrue( i == libsbml.LIBSBML_LEVEL_MISMATCH );
      assertTrue( m.getNumParameters() == 0 );
      p = null;
      m = null;
    }

    public void test_Model_addParameter4()
    {
      Model m = new  Model(2,2);
      Parameter p = null;
      int i = m.addParameter(p);
      assertTrue( i == libsbml.LIBSBML_OPERATION_FAILED );
      assertTrue( m.getNumParameters() == 0 );
      m = null;
    }

    public void test_Model_addParameter5()
    {
      Model m = new  Model(2,2);
      Parameter p = new  Parameter(2,2);
      p.setId( "p");
      Parameter p1 = new  Parameter(2,2);
      p1.setId( "p");
      int i = m.addParameter(p);
      assertTrue( i == libsbml.LIBSBML_OPERATION_SUCCESS );
      assertTrue( m.getNumParameters() == 1 );
      i = m.addParameter(p1);
      assertTrue( i == libsbml.LIBSBML_DUPLICATE_OBJECT_ID );
      assertTrue( m.getNumParameters() == 1 );
      p = null;
      p1 = null;
      m = null;
    }

    public void test_Model_addReaction1()
    {
      Model m = new  Model(2,2);
      Reaction r = new  Reaction(2,2);
      int i = m.addReaction(r);
      assertTrue( i == libsbml.LIBSBML_INVALID_OBJECT );
      r.setId( "r");
      i = m.addReaction(r);
      assertTrue( i == libsbml.LIBSBML_OPERATION_SUCCESS );
      assertTrue( m.getNumReactions() == 1 );
      r = null;
      m = null;
    }

    public void test_Model_addReaction2()
    {
      Model m = new  Model(2,2);
      Reaction r = new  Reaction(2,1);
      r.setId( "r");
      int i = m.addReaction(r);
      assertTrue( i == libsbml.LIBSBML_VERSION_MISMATCH );
      assertTrue( m.getNumReactions() == 0 );
      r = null;
      m = null;
    }

    public void test_Model_addReaction3()
    {
      Model m = new  Model(2,2);
      Reaction r = new  Reaction(1,2);
      r.setId( "r");
      int i = m.addReaction(r);
      assertTrue( i == libsbml.LIBSBML_LEVEL_MISMATCH );
      assertTrue( m.getNumReactions() == 0 );
      r = null;
      m = null;
    }

    public void test_Model_addReaction4()
    {
      Model m = new  Model(2,2);
      Reaction r = null;
      int i = m.addReaction(r);
      assertTrue( i == libsbml.LIBSBML_OPERATION_FAILED );
      assertTrue( m.getNumReactions() == 0 );
      m = null;
    }

    public void test_Model_addReaction5()
    {
      Model m = new  Model(2,2);
      Reaction r = new  Reaction(2,2);
      r.setId( "r");
      Reaction r1 = new  Reaction(2,2);
      r1.setId( "r");
      int i = m.addReaction(r);
      assertTrue( i == libsbml.LIBSBML_OPERATION_SUCCESS );
      assertTrue( m.getNumReactions() == 1 );
      i = m.addReaction(r1);
      assertTrue( i == libsbml.LIBSBML_DUPLICATE_OBJECT_ID );
      assertTrue( m.getNumReactions() == 1 );
      r = null;
      r1 = null;
      m = null;
    }

    public void test_Model_addRule1()
    {
      Model m = new  Model(2,2);
      Rule r = new  AssignmentRule(2,2);
      int i = m.addRule(r);
      assertTrue( i == libsbml.LIBSBML_INVALID_OBJECT );
      r.setVariable( "f");
      i = m.addRule(r);
      assertTrue( i == libsbml.LIBSBML_INVALID_OBJECT );
      r.setMath(libsbml.parseFormula("a-n"));
      i = m.addRule(r);
      assertTrue( i == libsbml.LIBSBML_OPERATION_SUCCESS );
      assertTrue( m.getNumRules() == 1 );
      r = null;
      m = null;
    }

    public void test_Model_addRule2()
    {
      Model m = new  Model(2,2);
      Rule r = new  AssignmentRule(2,1);
      r.setVariable( "f");
      r.setMath(libsbml.parseFormula("a-n"));
      int i = m.addRule(r);
      assertTrue( i == libsbml.LIBSBML_VERSION_MISMATCH );
      assertTrue( m.getNumRules() == 0 );
      r = null;
      m = null;
    }

    public void test_Model_addRule3()
    {
      Model m = new  Model(2,2);
      Rule r = new  AssignmentRule(1,2);
      r.setVariable( "f");
      r.setMath(libsbml.parseFormula("a-n"));
      int i = m.addRule(r);
      assertTrue( i == libsbml.LIBSBML_LEVEL_MISMATCH );
      assertTrue( m.getNumRules() == 0 );
      r = null;
      m = null;
    }

    public void test_Model_addRule4()
    {
      Model m = new  Model(2,2);
      Rule r = null;
      int i = m.addRule(r);
      assertTrue( i == libsbml.LIBSBML_OPERATION_FAILED );
      assertTrue( m.getNumRules() == 0 );
      m = null;
    }

    public void test_Model_addRule5()
    {
      Model m = new  Model(2,2);
      Rule ar = new  AssignmentRule(2,2);
      ar.setVariable( "ar");
      ar.setMath(libsbml.parseFormula("a-j"));
      Rule ar1 = new  AssignmentRule(2,2);
      ar1.setVariable( "ar");
      ar1.setMath(libsbml.parseFormula("a-j"));
      int i = m.addRule(ar);
      assertTrue( i == libsbml.LIBSBML_OPERATION_SUCCESS );
      assertTrue( m.getNumRules() == 1 );
      i = m.addRule(ar1);
      assertTrue( i == libsbml.LIBSBML_DUPLICATE_OBJECT_ID );
      assertTrue( m.getNumRules() == 1 );
      ar = null;
      ar1 = null;
      m = null;
    }

    public void test_Model_addSpecies1()
    {
      Model m = new  Model(2,2);
      Species s = new  Species(2,2);
      int i = m.addSpecies(s);
      assertTrue( i == libsbml.LIBSBML_INVALID_OBJECT );
      s.setId( "s");
      i = m.addSpecies(s);
      assertTrue( i == libsbml.LIBSBML_INVALID_OBJECT );
      s.setCompartment( "c");
      i = m.addSpecies(s);
      assertTrue( i == libsbml.LIBSBML_OPERATION_SUCCESS );
      assertTrue( m.getNumSpecies() == 1 );
      s = null;
      m = null;
    }

    public void test_Model_addSpecies2()
    {
      Model m = new  Model(2,2);
      Species s = new  Species(2,1);
      s.setId( "s");
      s.setCompartment( "c");
      int i = m.addSpecies(s);
      assertTrue( i == libsbml.LIBSBML_VERSION_MISMATCH );
      assertTrue( m.getNumSpecies() == 0 );
      s = null;
      m = null;
    }

    public void test_Model_addSpecies3()
    {
      Model m = new  Model(2,2);
      Species s = new  Species(1,2);
      s.setId( "s");
      s.setCompartment( "c");
      s.setInitialAmount(2);
      int i = m.addSpecies(s);
      assertTrue( i == libsbml.LIBSBML_LEVEL_MISMATCH );
      assertTrue( m.getNumSpecies() == 0 );
      s = null;
      m = null;
    }

    public void test_Model_addSpecies4()
    {
      Model m = new  Model(2,2);
      Species s = null;
      int i = m.addSpecies(s);
      assertTrue( i == libsbml.LIBSBML_OPERATION_FAILED );
      assertTrue( m.getNumSpecies() == 0 );
      m = null;
    }

    public void test_Model_addSpecies5()
    {
      Model m = new  Model(2,2);
      Species s = new  Species(2,2);
      s.setId( "s");
      s.setCompartment( "c");
      Species s1 = new  Species(2,2);
      s1.setId( "s");
      s1.setCompartment( "c");
      int i = m.addSpecies(s);
      assertTrue( i == libsbml.LIBSBML_OPERATION_SUCCESS );
      assertTrue( m.getNumSpecies() == 1 );
      i = m.addSpecies(s1);
      assertTrue( i == libsbml.LIBSBML_DUPLICATE_OBJECT_ID );
      assertTrue( m.getNumSpecies() == 1 );
      s = null;
      s1 = null;
      m = null;
    }

    public void test_Model_addSpeciesType1()
    {
      Model m = new  Model(2,2);
      SpeciesType st = new  SpeciesType(2,2);
      int i = m.addSpeciesType(st);
      assertTrue( i == libsbml.LIBSBML_INVALID_OBJECT );
      st.setId( "st");
      i = m.addSpeciesType(st);
      assertTrue( i == libsbml.LIBSBML_OPERATION_SUCCESS );
      assertTrue( m.getNumSpeciesTypes() == 1 );
      st = null;
      m = null;
    }

    public void test_Model_addSpeciesType2()
    {
      Model m = new  Model(2,2);
      SpeciesType st = new  SpeciesType(2,3);
      st.setId( "st");
      int i = m.addSpeciesType(st);
      assertTrue( i == libsbml.LIBSBML_VERSION_MISMATCH );
      assertTrue( m.getNumSpeciesTypes() == 0 );
      st = null;
      m = null;
    }

    public void test_Model_addSpeciesType3()
    {
      Model m = new  Model(2,2);
      SpeciesType st = null;
      int i = m.addSpeciesType(st);
      assertTrue( i == libsbml.LIBSBML_OPERATION_FAILED );
      assertTrue( m.getNumSpeciesTypes() == 0 );
      m = null;
    }

    public void test_Model_addSpeciesType4()
    {
      Model m = new  Model(2,2);
      SpeciesType st = new  SpeciesType(2,2);
      st.setId( "st");
      SpeciesType st1 = new  SpeciesType(2,2);
      st1.setId( "st");
      int i = m.addSpeciesType(st);
      assertTrue( i == libsbml.LIBSBML_OPERATION_SUCCESS );
      assertTrue( m.getNumSpeciesTypes() == 1 );
      i = m.addSpeciesType(st1);
      assertTrue( i == libsbml.LIBSBML_DUPLICATE_OBJECT_ID );
      assertTrue( m.getNumSpeciesTypes() == 1 );
      st = null;
      st1 = null;
      m = null;
    }

    public void test_Model_addUnitDefinition1()
    {
      Model m = new  Model(2,2);
      UnitDefinition ud = new  UnitDefinition(2,2);
      int i = m.addUnitDefinition(ud);
      assertTrue( i == libsbml.LIBSBML_INVALID_OBJECT );
      ud.createUnit();
      i = m.addUnitDefinition(ud);
      assertTrue( i == libsbml.LIBSBML_INVALID_OBJECT );
      ud.setId( "ud");
      i = m.addUnitDefinition(ud);
      assertTrue( i == libsbml.LIBSBML_OPERATION_SUCCESS );
      assertTrue( m.getNumUnitDefinitions() == 1 );
      ud = null;
      m = null;
    }

    public void test_Model_addUnitDefinition2()
    {
      Model m = new  Model(2,2);
      UnitDefinition ud = new  UnitDefinition(2,1);
      ud.createUnit();
      ud.setId( "ud");
      int i = m.addUnitDefinition(ud);
      assertTrue( i == libsbml.LIBSBML_VERSION_MISMATCH );
      assertTrue( m.getNumUnitDefinitions() == 0 );
      ud = null;
      m = null;
    }

    public void test_Model_addUnitDefinition3()
    {
      Model m = new  Model(2,2);
      UnitDefinition ud = new  UnitDefinition(1,2);
      ud.createUnit();
      ud.setId( "ud");
      int i = m.addUnitDefinition(ud);
      assertTrue( i == libsbml.LIBSBML_LEVEL_MISMATCH );
      assertTrue( m.getNumUnitDefinitions() == 0 );
      ud = null;
      m = null;
    }

    public void test_Model_addUnitDefinition4()
    {
      Model m = new  Model(2,2);
      UnitDefinition ud = null;
      int i = m.addUnitDefinition(ud);
      assertTrue( i == libsbml.LIBSBML_OPERATION_FAILED );
      assertTrue( m.getNumUnitDefinitions() == 0 );
      m = null;
    }

    public void test_Model_addUnitDefinition5()
    {
      Model m = new  Model(2,2);
      UnitDefinition ud = new  UnitDefinition(2,2);
      ud.setId( "ud");
      ud.createUnit();
      UnitDefinition ud1 = new  UnitDefinition(2,2);
      ud1.setId( "ud");
      ud1.createUnit();
      int i = m.addUnitDefinition(ud);
      assertTrue( i == libsbml.LIBSBML_OPERATION_SUCCESS );
      assertTrue( m.getNumUnitDefinitions() == 1 );
      i = m.addUnitDefinition(ud1);
      assertTrue( i == libsbml.LIBSBML_DUPLICATE_OBJECT_ID );
      assertTrue( m.getNumUnitDefinitions() == 1 );
      ud = null;
      ud1 = null;
      m = null;
    }

    public void test_Model_createCompartment()
    {
      Model m = new  Model(2,2);
      Compartment p = m.createCompartment();
      assertTrue( m.getNumCompartments() == 1 );
      assertTrue( (p).getLevel() == 2 );
      assertTrue( (p).getVersion() == 2 );
      m = null;
    }

    public void test_Model_createCompartmentType()
    {
      Model m = new  Model(2,2);
      CompartmentType p = m.createCompartmentType();
      assertTrue( m.getNumCompartmentTypes() == 1 );
      assertTrue( (p).getLevel() == 2 );
      assertTrue( (p).getVersion() == 2 );
      m = null;
    }

    public void test_Model_createConstraint()
    {
      Model m = new  Model(2,2);
      Constraint p = m.createConstraint();
      assertTrue( m.getNumConstraints() == 1 );
      assertTrue( (p).getLevel() == 2 );
      assertTrue( (p).getVersion() == 2 );
      m = null;
    }

    public void test_Model_createEvent()
    {
      Model m = new  Model(2,2);
      Event p = m.createEvent();
      assertTrue( m.getNumEvents() == 1 );
      assertTrue( (p).getLevel() == 2 );
      assertTrue( (p).getVersion() == 2 );
      m = null;
    }

    public void test_Model_createEventAssignment()
    {
      Model m = new  Model(2,2);
      Event p = m.createEvent();
      EventAssignment ea = m.createEventAssignment();
      assertTrue( p.getNumEventAssignments() == 1 );
      assertTrue( (ea).getLevel() == 2 );
      assertTrue( (ea).getVersion() == 2 );
      m = null;
    }

    public void test_Model_createFunctionDefinition()
    {
      Model m = new  Model(2,2);
      FunctionDefinition p = m.createFunctionDefinition();
      assertTrue( m.getNumFunctionDefinitions() == 1 );
      assertTrue( (p).getLevel() == 2 );
      assertTrue( (p).getVersion() == 2 );
      m = null;
    }

    public void test_Model_createInitialAssignment()
    {
      Model m = new  Model(2,2);
      InitialAssignment p = m.createInitialAssignment();
      assertTrue( m.getNumInitialAssignments() == 1 );
      assertTrue( (p).getLevel() == 2 );
      assertTrue( (p).getVersion() == 2 );
      m = null;
    }

    public void test_Model_createKineticLaw()
    {
      Model m = new  Model(2,2);
      Reaction p = m.createReaction();
      KineticLaw kl = m.createKineticLaw();
      assertTrue( p.isSetKineticLaw() == true );
      assertTrue( (kl).getLevel() == 2 );
      assertTrue( (kl).getVersion() == 2 );
      m = null;
    }

    public void test_Model_createKineticLawParameters()
    {
      Model m = new  Model(2,2);
      Reaction r = m.createReaction();
      KineticLaw kl = m.createKineticLaw();
      Parameter p = m.createKineticLawParameter();
      assertTrue( r.isSetKineticLaw() == true );
      assertTrue( kl.getNumParameters() == 1 );
      assertTrue( (p).getLevel() == 2 );
      assertTrue( (p).getVersion() == 2 );
      m = null;
    }

    public void test_Model_createModifier()
    {
      Model m = new  Model(2,2);
      Reaction p = m.createReaction();
      ModifierSpeciesReference sr = m.createModifier();
      assertTrue( p.getNumModifiers() == 1 );
      assertTrue( (sr).getLevel() == 2 );
      assertTrue( (sr).getVersion() == 2 );
      m = null;
    }

    public void test_Model_createParameter()
    {
      Model m = new  Model(2,2);
      Parameter p = m.createParameter();
      assertTrue( m.getNumParameters() == 1 );
      assertTrue( (p).getLevel() == 2 );
      assertTrue( (p).getVersion() == 2 );
      m = null;
    }

    public void test_Model_createProduct()
    {
      Model m = new  Model(2,2);
      Reaction p = m.createReaction();
      SpeciesReference sr = m.createProduct();
      assertTrue( p.getNumProducts() == 1 );
      assertTrue( (sr).getLevel() == 2 );
      assertTrue( (sr).getVersion() == 2 );
      m = null;
    }

    public void test_Model_createReactant()
    {
      Model m = new  Model(2,2);
      Reaction p = m.createReaction();
      SpeciesReference sr = m.createReactant();
      assertTrue( p.getNumReactants() == 1 );
      assertTrue( (sr).getLevel() == 2 );
      assertTrue( (sr).getVersion() == 2 );
      m = null;
    }

    public void test_Model_createReaction()
    {
      Model m = new  Model(2,2);
      Reaction p = m.createReaction();
      assertTrue( m.getNumReactions() == 1 );
      assertTrue( (p).getLevel() == 2 );
      assertTrue( (p).getVersion() == 2 );
      m = null;
    }

    public void test_Model_createRule()
    {
      Model m = new  Model(2,2);
      Rule p = m.createAssignmentRule();
      assertTrue( m.getNumRules() == 1 );
      assertTrue( (p).getLevel() == 2 );
      assertTrue( (p).getVersion() == 2 );
      m = null;
    }

    public void test_Model_createSpecies()
    {
      Model m = new  Model(2,2);
      Species p = m.createSpecies();
      assertTrue( m.getNumSpecies() == 1 );
      assertTrue( (p).getLevel() == 2 );
      assertTrue( (p).getVersion() == 2 );
      m = null;
    }

    public void test_Model_createSpeciesType()
    {
      Model m = new  Model(2,2);
      SpeciesType p = m.createSpeciesType();
      assertTrue( m.getNumSpeciesTypes() == 1 );
      assertTrue( (p).getLevel() == 2 );
      assertTrue( (p).getVersion() == 2 );
      m = null;
    }

    public void test_Model_createUnit()
    {
      Model m = new  Model(2,2);
      UnitDefinition p = m.createUnitDefinition();
      Unit u = m.createUnit();
      assertTrue( p.getNumUnits() == 1 );
      assertTrue( (u).getLevel() == 2 );
      assertTrue( (u).getVersion() == 2 );
      m = null;
    }

    public void test_Model_createUnitDefinition()
    {
      Model m = new  Model(2,2);
      UnitDefinition p = m.createUnitDefinition();
      assertTrue( m.getNumUnitDefinitions() == 1 );
      assertTrue( (p).getLevel() == 2 );
      assertTrue( (p).getVersion() == 2 );
      m = null;
    }

    public void test_Model_setId1()
    {
      string id =  "1e1";
      int i = M.setId(id);
      assertTrue( i == libsbml.LIBSBML_INVALID_ATTRIBUTE_VALUE );
      assertEquals( false, M.isSetId() );
    }

    public void test_Model_setId2()
    {
      string id =  "e1";
      int i = M.setId(id);
      assertTrue( i == libsbml.LIBSBML_OPERATION_SUCCESS );
      assertTrue(( id == M.getId() ));
      assertEquals( true, M.isSetId() );
      i = M.setId("");
      assertTrue( i == libsbml.LIBSBML_OPERATION_SUCCESS );
      assertEquals( false, M.isSetId() );
    }

    public void test_Model_setId3()
    {
      string id =  "e1";
      int i = M.setId(id);
      assertTrue( i == libsbml.LIBSBML_OPERATION_SUCCESS );
      assertTrue(( id == M.getId() ));
      assertEquals( true, M.isSetId() );
      i = M.unsetId();
      assertTrue( i == libsbml.LIBSBML_OPERATION_SUCCESS );
      assertEquals( false, M.isSetId() );
    }

    public void test_Model_setModelHistory1()
    {
      ModelHistory mh = new  ModelHistory();
      int i = M.setModelHistory(mh);
      assertTrue( i == libsbml.LIBSBML_INVALID_OBJECT );
      assertEquals( false, M.isSetModelHistory() );
      i = M.unsetModelHistory();
      assertTrue( i == libsbml.LIBSBML_OPERATION_SUCCESS );
      assertEquals( false, M.isSetModelHistory() );
      mh = null;
    }

    public void test_Model_setModelHistory2()
    {
      int i = M.setModelHistory(null);
      assertTrue( i == libsbml.LIBSBML_OPERATION_SUCCESS );
      assertEquals( false, M.isSetModelHistory() );
      i = M.unsetModelHistory();
      assertTrue( i == libsbml.LIBSBML_OPERATION_SUCCESS );
      assertEquals( false, M.isSetModelHistory() );
    }

    public void test_Model_setName1()
    {
      string name =  "3Set_k2";
      int i = M.setName(name);
      assertTrue( i == libsbml.LIBSBML_OPERATION_SUCCESS );
      assertEquals( true, M.isSetName() );
    }

    public void test_Model_setName2()
    {
      string name =  "Set k2";
      int i = M.setName(name);
      assertTrue( i == libsbml.LIBSBML_OPERATION_SUCCESS );
      assertTrue(( name == M.getName() ));
      assertEquals( true, M.isSetName() );
      i = M.unsetName();
      assertTrue( i == libsbml.LIBSBML_OPERATION_SUCCESS );
      assertEquals( false, M.isSetName() );
    }

    public void test_Model_setName3()
    {
      int i = M.setName("");
      assertTrue( i == libsbml.LIBSBML_OPERATION_SUCCESS );
      assertEquals( false, M.isSetName() );
    }

    public void test_Model_setName4()
    {
      Model m = new  Model(1,2);
      int i = m.setName( "11dd");
      assertTrue( i == libsbml.LIBSBML_INVALID_ATTRIBUTE_VALUE );
      assertEquals( false, m.isSetName() );
    }

  }
}
