///  @file    TestSBMLConvert.cs
///  @brief   SBMLConvert unit tests
///  @author  Frank Bergmann (Csharp conversion)
///  @author  Akiya Jouraku (Csharp conversion)
///  @author  Ben Bornstein 
/// 
///  $Id: TestSBMLConvert.cs 11545 2010-07-23 02:19:10Z mhucka $
///  $HeadURL: https://sbml.svn.sourceforge.net/svnroot/sbml/trunk/libsbml/src/bindings/csharp/test/sbml/TestSBMLConvert.cs $
/// 
///  ====== WARNING ===== WARNING ===== WARNING ===== WARNING ===== WARNING ======
/// 
///  DO NOT EDIT THIS FILE.
/// 
///  This file was generated automatically by converting the file located at
///  src/sbml/test/TestSBMLConvert.c
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

  public class TestSBMLConvert {
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


    public void test_SBMLConvert_addModifiersToReaction()
    {
      SBMLDocument d = new  SBMLDocument(1,2);
      Model m = d.createModel();
      Reaction r = m.createReaction();
      KineticLaw kl = r.createKineticLaw();
      kl.setFormula( "k1*S1*S2*S3*S4*S5");
      SimpleSpeciesReference ssr1;
      SimpleSpeciesReference ssr2;
      Species s1 = m.createSpecies();
      s1.setId( "S1" );
      Species s2 = m.createSpecies();
      s2.setId( "S2");
      Species s3 = m.createSpecies();
      s3.setId( "S3");
      Species s4 = m.createSpecies();
      s4.setId( "S4");
      Species s5 = m.createSpecies();
      s5.setId( "S5");
      SpeciesReference sr1 = r.createReactant();
      SpeciesReference sr2 = r.createReactant();
      SpeciesReference sr3 = r.createProduct();
      sr1.setSpecies( "S1");
      sr2.setSpecies( "S2");
      sr3.setSpecies( "S5");
      assertTrue( r.getNumModifiers() == 0 );
      assertTrue( d.setLevelAndVersion(2,1,false) == true );
      assertTrue( d.getLevel() == 2 );
      assertTrue( d.getVersion() == 1 );
      assertTrue( m.getReaction(0).getNumModifiers() == 2 );
      ssr1 = m.getReaction(0).getModifier(0);
      ssr2 = m.getReaction(0).getModifier(1);
      assertTrue((  "S3" == ssr1.getSpecies() ));
      assertTrue((  "S4" == ssr2.getSpecies() ));
      d = null;
    }

    public void test_SBMLConvert_convertFromL3()
    {
      SBMLDocument d = new  SBMLDocument(3,1);
      Model m = d.createModel();
      string sid =  "C";
      Compartment c = m.createCompartment();
      c.setId(sid);
      c.setSize(1.2);
      c.setUnits( "volume");
      assertTrue( d.setLevelAndVersion(1,1,true) == false );
      assertTrue( d.setLevelAndVersion(1,2,true) == false );
      assertTrue( d.setLevelAndVersion(2,1,true) == false );
      assertTrue( d.setLevelAndVersion(2,2,true) == false );
      assertTrue( d.setLevelAndVersion(2,3,true) == false );
      assertTrue( d.setLevelAndVersion(2,4,true) == false );
      assertTrue( d.setLevelAndVersion(3,1,true) == true );
    }

    public void test_SBMLConvert_convertToL1_SBMLDocument()
    {
      SBMLDocument d = new  SBMLDocument(2,1);
      assertTrue( d.setLevelAndVersion(1,2,false) == true );
      assertTrue( d.getLevel() == 1 );
      assertTrue( d.getVersion() == 2 );
      d = null;
    }

    public void test_SBMLConvert_convertToL1_Species_Amount()
    {
      SBMLDocument d = new  SBMLDocument(2,1);
      Model m = d.createModel();
      string sid =  "C";
      Compartment c = new  Compartment(2,4);
      Species s = new  Species(2,4);
      c.setId(sid);
      m.addCompartment(c);
      s.setCompartment(sid);
      s.setInitialAmount(2.34);
      m.addSpecies(s);
      assertTrue( d.setLevelAndVersion(1,2,true) == true );
      assertTrue( s.getInitialAmount() == 2.34 );
      d = null;
    }

    public void test_SBMLConvert_convertToL1_Species_Concentration()
    {
      SBMLDocument d = new  SBMLDocument(2,1);
      Model m = d.createModel();
      string sid =  "C";
      Compartment c = new  Compartment(2,1);
      Species s = new  Species(2,1);
      c.setId(sid);
      c.setSize(1.2);
      m.addCompartment(c);
      s.setId( "s"  );
      s.setCompartment(sid);
      s.setInitialConcentration(2.34);
      m.addSpecies(s);
      assertTrue( d.setLevelAndVersion(1,2,true) == true );
      Species s1 = m.getSpecies(0);
      assertTrue( s1 != null );
      assertTrue((  "C" == s1.getCompartment() ));
      assertTrue( m.getCompartment( "C").getSize() == 1.2 );
      assertTrue( s1.getInitialConcentration() == 2.34 );
      assertTrue( s1.isSetInitialConcentration() == true );
      d = null;
    }

    public void test_SBMLConvert_convertToL2_SBMLDocument()
    {
      SBMLDocument d = new  SBMLDocument(1,2);
      assertTrue( d.setLevelAndVersion(2,1,false) == true );
      assertTrue( d.getLevel() == 2 );
      assertTrue( d.getVersion() == 1 );
      assertTrue( d.setLevelAndVersion(2,2,false) == true );
      assertTrue( d.getLevel() == 2 );
      assertTrue( d.getVersion() == 2 );
      assertTrue( d.setLevelAndVersion(2,3,false) == true );
      assertTrue( d.getLevel() == 2 );
      assertTrue( d.getVersion() == 3 );
      d = null;
    }

    public void test_SBMLConvert_convertToL2v4_DuplicateAnnotations_doc()
    {
      SBMLDocument d = new  SBMLDocument(2,1);
      d.createModel();
      string annotation =  "<rdf/>\n<rdf/>";
      int i = (d).setAnnotation(annotation);
      assertTrue( d.getLevel() == 2 );
      assertTrue( d.getVersion() == 1 );
      assertTrue( (d).getAnnotation().getNumChildren() == 2 );
      assertTrue( d.setLevelAndVersion(2,4,true) == true );
      assertTrue( d.getLevel() == 2 );
      assertTrue( d.getVersion() == 4 );
      assertTrue( (d).getAnnotation().getNumChildren() == 1 );
      d = null;
    }

    public void test_SBMLConvert_convertToL2v4_DuplicateAnnotations_model()
    {
      SBMLDocument d = new  SBMLDocument(2,1);
      Model m = d.createModel();
      string annotation =  "<rdf/>\n<rdf/>";
      int i = (m).setAnnotation(annotation);
      assertTrue( d.getLevel() == 2 );
      assertTrue( d.getVersion() == 1 );
      assertTrue( (m).getAnnotation().getNumChildren() == 2 );
      assertTrue( d.setLevelAndVersion(2,4,true) == true );
      assertTrue( d.getLevel() == 2 );
      assertTrue( d.getVersion() == 4 );
      m = d.getModel();
      assertTrue( (m).getAnnotation().getNumChildren() == 1 );
      d = null;
    }

    public void test_SBMLConvert_convertToL3_defaultUnits()
    {
      SBMLDocument d = new  SBMLDocument(1,2);
      Model m = d.createModel();
      string sid =  "C";
      Compartment c = m.createCompartment();
      c.setId(sid);
      c.setSize(1.2);
      c.setUnits( "volume");
      assertTrue( m.getNumUnitDefinitions() == 0 );
      assertTrue( d.setLevelAndVersion(3,1,true) == true );
      assertTrue( m.getNumUnitDefinitions() == 2 );
      UnitDefinition ud = m.getUnitDefinition(0);
      assertTrue( ud != null );
      assertTrue((  "volume" == ud.getId() ));
      assertTrue( ud.getNumUnits() == 1 );
      Unit u = ud.getUnit(0);
      assertTrue( u.getKind() == libsbml.UNIT_KIND_LITRE );
      assertTrue( u.getExponent() == 1 );
      assertTrue( u.getMultiplier() == 1 );
      assertTrue( u.getScale() == 0 );
      ud = m.getUnitDefinition(1);
      assertTrue( ud != null );
      assertTrue((  "time" == ud.getId() ));
      assertTrue( ud.getNumUnits() == 1 );
      u = ud.getUnit(0);
      assertTrue( u.getKind() == libsbml.UNIT_KIND_SECOND );
      assertTrue( u.getExponent() == 1 );
      assertTrue( u.getMultiplier() == 1 );
      assertTrue( u.getScale() == 0 );
      assertTrue((  "time" == m.getTimeUnits() ));
      d = null;
    }

    public void test_SBMLConvert_convertToL3_localParameters()
    {
      SBMLDocument d = new  SBMLDocument(1,2);
      Model m = d.createModel();
      Compartment c = m.createCompartment();
      c.setId( "c" );
      Species s = m.createSpecies();
      s.setId( "s");
      s.setCompartment( "c");
      Reaction r = m.createReaction();
      SpeciesReference sr = r.createReactant();
      sr.setSpecies( "s");
      KineticLaw kl = r.createKineticLaw();
      kl.setFormula( "s*k");
      Parameter p = kl.createParameter();
      p.setId( "k");
      assertTrue( kl.getNumLocalParameters() == 0 );
      assertTrue( d.setLevelAndVersion(3,1,false) == true );
      m = d.getModel();
      r = m.getReaction(0);
      kl = r.getKineticLaw();
      assertTrue( kl.getNumLocalParameters() == 1 );
      LocalParameter lp = kl.getLocalParameter(0);
      d = null;
    }

    public void test_SBMLConvert_convertToL3_stoichiometryMath()
    {
      SBMLDocument d = new  SBMLDocument(2,1);
      Model m = d.createModel();
      Compartment c = m.createCompartment();
      c.setId( "c" );
      Species s = m.createSpecies();
      s.setId( "s");
      s.setCompartment( "c");
      Reaction r = m.createReaction();
      SpeciesReference sr = r.createReactant();
      sr.setSpecies( "s");
      StoichiometryMath sm = sr.createStoichiometryMath();
      ASTNode ast = libsbml.parseFormula("c*2");
      sm.setMath(ast);
      assertTrue( m.getNumRules() == 0 );
      assertTrue( sr.isSetId() == false );
      assertTrue( d.setLevelAndVersion(3,1,false) == true );
      m = d.getModel();
      r = m.getReaction(0);
      sr = r.getReactant(0);
      assertTrue( m.getNumRules() == 1 );
      assertTrue( sr.isSetId() == true );
      Rule rule = m.getRule(0);
      assertTrue( sr.getId() == rule.getVariable() );
      d = null;
    }

    public void test_SBMLConvert_invalidLevelVersion()
    {
      SBMLDocument d = new  SBMLDocument(2,1);
      Model m = d.createModel();
      string sid =  "C";
      Compartment c = m.createCompartment();
      c.setId(sid);
      c.setSize(1.2);
      c.setUnits( "volume");
      assertTrue( d.setLevelAndVersion(1,3,true) == false );
      assertTrue( d.setLevelAndVersion(2,5,true) == false );
      assertTrue( d.setLevelAndVersion(3,2,true) == false );
      assertTrue( d.setLevelAndVersion(4,1,true) == false );
    }

  }
}
