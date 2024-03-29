///  @file    TestRequiredElements.cs
///  @brief   Test hasRequiredElements unit tests
///  @author  Frank Bergmann (Csharp conversion)
///  @author  Akiya Jouraku (Csharp conversion)
///  @author  Sarah Keating 
/// 
///  $Id: TestRequiredElements.cs 11545 2010-07-23 02:19:10Z mhucka $
///  $HeadURL: https://sbml.svn.sourceforge.net/svnroot/sbml/trunk/libsbml/src/bindings/csharp/test/sbml/TestRequiredElements.cs $
/// 
///  ====== WARNING ===== WARNING ===== WARNING ===== WARNING ===== WARNING ======
/// 
///  DO NOT EDIT THIS FILE.
/// 
///  This file was generated automatically by converting the file located at
///  src/sbml/test/TestRequiredElements.cpp
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

  public class TestRequiredElements {
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


    public void test_AlgebraicRule()
    {
      AlgebraicRule ar = new AlgebraicRule(2,4);
      assertEquals( false, (ar.hasRequiredElements()) );
      ar.setMath(libsbml.parseFormula("ar"));
      assertEquals( true, ar.hasRequiredElements() );
      ar = null;
    }

    public void test_AssignmentRule()
    {
      AssignmentRule r = new AssignmentRule(2,4);
      assertEquals( false, (r.hasRequiredElements()) );
      r.setMath(libsbml.parseFormula("ar"));
      assertEquals( true, r.hasRequiredElements() );
      r = null;
    }

    public void test_Compartment()
    {
      Compartment c = new Compartment(2,4);
      assertEquals( true, c.hasRequiredElements() );
      c = null;
    }

    public void test_CompartmentType()
    {
      CompartmentType ct = new CompartmentType(2,4);
      assertEquals( true, ct.hasRequiredElements() );
      ct = null;
    }

    public void test_Constraint()
    {
      Constraint c = new Constraint(2,4);
      assertEquals( false, (c.hasRequiredElements()) );
      c.setMath(libsbml.parseFormula("a+b"));
      assertEquals( true, c.hasRequiredElements() );
      c = null;
    }

    public void test_Delay()
    {
      Delay d = new Delay(2,4);
      assertEquals( false, (d.hasRequiredElements()) );
      d.setMath(libsbml.parseFormula("a+b"));
      assertEquals( true, d.hasRequiredElements() );
      d = null;
    }

    public void test_Event()
    {
      Event e = new Event(2,4);
      assertEquals( false, (e.hasRequiredElements()) );
      Trigger t = new Trigger(2,4);
      e.setTrigger(t);
      assertEquals( false, (e.hasRequiredElements()) );
      e.createEventAssignment();
      assertEquals( true, e.hasRequiredElements() );
      e = null;
    }

    public void test_EventAssignment()
    {
      EventAssignment ea = new EventAssignment(2,4);
      assertEquals( false, (ea.hasRequiredElements()) );
      ea.setMath(libsbml.parseFormula("fd"));
      assertEquals( true, ea.hasRequiredElements() );
      ea = null;
    }

    public void test_FunctionDefinition()
    {
      FunctionDefinition fd = new FunctionDefinition(2,4);
      assertEquals( false, (fd.hasRequiredElements()) );
      fd.setMath(libsbml.parseFormula("fd"));
      assertEquals( true, fd.hasRequiredElements() );
      fd = null;
    }

    public void test_InitialAssignment()
    {
      InitialAssignment ia = new InitialAssignment(2,4);
      assertEquals( false, (ia.hasRequiredElements()) );
      ia.setMath(libsbml.parseFormula("ia"));
      assertEquals( true, ia.hasRequiredElements() );
      ia = null;
    }

    public void test_KineticLaw()
    {
      KineticLaw kl = new KineticLaw(2,4);
      assertEquals( false, (kl.hasRequiredElements()) );
      kl.setMath(libsbml.parseFormula("kl"));
      assertEquals( true, kl.hasRequiredElements() );
      kl = null;
    }

    public void test_Model()
    {
      Model m = new Model(2,4);
      assertEquals( true, m.hasRequiredElements() );
      m = null;
    }

    public void test_Model_L1V1()
    {
      Model m = new Model(1,1);
      assertEquals( false, (m.hasRequiredElements()) );
      m.createCompartment();
      assertEquals( false, (m.hasRequiredElements()) );
      m.createSpecies();
      assertEquals( false, (m.hasRequiredElements()) );
      m.createReaction();
      assertEquals( true, m.hasRequiredElements() );
      m = null;
    }

    public void test_Model_L1V2()
    {
      Model m = new Model(1,2);
      assertEquals( false, (m.hasRequiredElements()) );
      m.createCompartment();
      assertEquals( true, m.hasRequiredElements() );
      m = null;
    }

    public void test_ModifierSpeciesReference()
    {
      ModifierSpeciesReference msr = new ModifierSpeciesReference(2,4);
      assertEquals( true, msr.hasRequiredElements() );
      msr = null;
    }

    public void test_Parameter()
    {
      Parameter p = new Parameter(2,4);
      assertEquals( true, p.hasRequiredElements() );
      p = null;
    }

    public void test_RateRule()
    {
      RateRule r = new RateRule(2,4);
      assertEquals( false, (r.hasRequiredElements()) );
      r.setMath(libsbml.parseFormula("ar"));
      assertEquals( true, r.hasRequiredElements() );
      r = null;
    }

    public void test_Reaction()
    {
      Reaction r = new Reaction(2,4);
      assertEquals( true, r.hasRequiredElements() );
      r = null;
    }

    public void test_Species()
    {
      Species s = new Species(2,4);
      assertEquals( true, s.hasRequiredElements() );
      s = null;
    }

    public void test_SpeciesReference()
    {
      SpeciesReference sr = new SpeciesReference(2,4);
      assertEquals( true, sr.hasRequiredElements() );
      sr = null;
    }

    public void test_SpeciesType()
    {
      SpeciesType st = new SpeciesType(2,4);
      assertEquals( true, st.hasRequiredElements() );
      st = null;
    }

    public void test_StoichiometryMath()
    {
      StoichiometryMath sm = new StoichiometryMath(2,4);
      assertEquals( false, (sm.hasRequiredElements()) );
      sm.setMath(libsbml.parseFormula("ar"));
      assertEquals( true, sm.hasRequiredElements() );
      sm = null;
    }

    public void test_Trigger()
    {
      Trigger t = new Trigger(2,4);
      assertEquals( false, (t.hasRequiredElements()) );
      t.setMath(libsbml.parseFormula("ar"));
      assertEquals( true, t.hasRequiredElements() );
      t = null;
    }

    public void test_Unit()
    {
      Unit u = new Unit(2,4);
      assertEquals( true, u.hasRequiredElements() );
      u = null;
    }

    public void test_UnitDefinition()
    {
      UnitDefinition ud = new UnitDefinition(2,4);
      assertEquals( false, (ud.hasRequiredElements()) );
      ud.createUnit();
      assertEquals( true, ud.hasRequiredElements() );
      ud = null;
    }

    public void test_UnitDefinition_L1()
    {
      UnitDefinition ud = new UnitDefinition(1,2);
      assertEquals( true, ud.hasRequiredElements() );
      ud = null;
    }

  }
}
