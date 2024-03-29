///  @file    TestL3Unit.cs
///  @brief   L3 Unit unit tests
///  @author  Frank Bergmann (Csharp conversion)
///  @author  Akiya Jouraku (Csharp conversion)
///  @author  Sarah Keating 
/// 
///  $Id: TestL3Unit.cs 11545 2010-07-23 02:19:10Z mhucka $
///  $HeadURL: https://sbml.svn.sourceforge.net/svnroot/sbml/trunk/libsbml/src/bindings/csharp/test/sbml/TestL3Unit.cs $
/// 
///  ====== WARNING ===== WARNING ===== WARNING ===== WARNING ===== WARNING ======
/// 
///  DO NOT EDIT THIS FILE.
/// 
///  This file was generated automatically by converting the file located at
///  src/sbml/test/TestL3Unit.c
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

  public class TestL3Unit {
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

    private Unit U;

    public bool isnan(double x)
    {
      return (x != x);
    }

  private const int SBML_INT_MAX = 2147483647;
    public void setUp()
    {
      U = new  Unit(3,1);
      if (U == null);
      {
      }
    }

    public void tearDown()
    {
      U = null;
    }

    public void test_L3_Unit_NS()
    {
      assertTrue( U.getNamespaces() != null );
      assertTrue( U.getNamespaces().getLength() == 1 );
      assertTrue((     "http://www.sbml.org/sbml/level3/version1/core" == U.getNamespaces().getURI(0) ));
    }

    public void test_L3_Unit_create()
    {
      assertTrue( U.getTypeCode() == libsbml.SBML_UNIT );
      assertTrue( U.getMetaId() == "" );
      assertTrue( U.getNotes() == null );
      assertTrue( U.getAnnotation() == null );
      assertTrue( U.getKind() == libsbml.UNIT_KIND_INVALID );
      assertEquals( true, isnan(U.getExponentAsDouble()) );
      assertEquals( true, isnan(U.getMultiplier()) );
      assertTrue( U.getScale() == SBML_INT_MAX );
      assertEquals( false, U.isSetKind() );
      assertEquals( false, U.isSetExponent() );
      assertEquals( false, U.isSetMultiplier() );
      assertEquals( false, U.isSetScale() );
    }

    public void test_L3_Unit_createWithNS()
    {
      XMLNamespaces xmlns = new  XMLNamespaces();
      xmlns.add( "http://www.sbml.org", "testsbml");
      SBMLNamespaces sbmlns = new  SBMLNamespaces(3,1);
      sbmlns.addNamespaces(xmlns);
      Unit u = new  Unit(sbmlns);
      assertTrue( u.getTypeCode() == libsbml.SBML_UNIT );
      assertTrue( u.getMetaId() == "" );
      assertTrue( u.getNotes() == null );
      assertTrue( u.getAnnotation() == null );
      assertTrue( u.getLevel() == 3 );
      assertTrue( u.getVersion() == 1 );
      assertTrue( u.getNamespaces() != null );
      assertTrue( u.getNamespaces().getLength() == 2 );
      assertTrue( u.getKind() == libsbml.UNIT_KIND_INVALID );
      assertEquals( true, isnan(u.getExponentAsDouble()) );
      assertEquals( true, isnan(u.getMultiplier()) );
      assertEquals( false, u.isSetKind() );
      assertEquals( false, u.isSetExponent() );
      assertEquals( false, u.isSetMultiplier() );
      assertEquals( false, u.isSetScale() );
      u = null;
    }

    public void test_L3_Unit_exponent()
    {
      double exponent = 0.2;
      assertEquals( false, U.isSetExponent() );
      assertEquals( true, isnan(U.getExponentAsDouble()) );
      U.setExponent(exponent);
      assertTrue( U.getExponentAsDouble() == exponent );
      assertEquals( true, U.isSetExponent() );
    }

    public void test_L3_Unit_free_NULL()
    {
    }

    public void test_L3_Unit_hasRequiredAttributes()
    {
      Unit u = new  Unit(3,1);
      assertEquals( false, u.hasRequiredAttributes() );
      u.setKind(libsbml.UNIT_KIND_MOLE);
      assertEquals( false, u.hasRequiredAttributes() );
      u.setExponent(0);
      assertEquals( false, u.hasRequiredAttributes() );
      u.setMultiplier(0.45);
      assertEquals( false, u.hasRequiredAttributes() );
      u.setScale(2);
      assertEquals( true, u.hasRequiredAttributes() );
      u = null;
    }

    public void test_L3_Unit_kind()
    {
      string kind =  "mole";
      assertEquals( false, U.isSetKind() );
      U.setKind(libsbml.UnitKind_forName(kind));
      assertTrue( U.getKind() == libsbml.UNIT_KIND_MOLE );
      assertEquals( true, U.isSetKind() );
    }

    public void test_L3_Unit_multiplier()
    {
      double multiplier = 0.2;
      assertEquals( false, U.isSetMultiplier() );
      assertEquals( true, isnan(U.getMultiplier()) );
      U.setMultiplier(multiplier);
      assertTrue( U.getMultiplier() == multiplier );
      assertEquals( true, U.isSetMultiplier() );
    }

    public void test_L3_Unit_scale()
    {
      int scale = 2;
      assertEquals( false, U.isSetScale() );
      assertTrue( U.getScale() == SBML_INT_MAX );
      U.setScale(scale);
      assertTrue( U.getScale() == scale );
      assertEquals( true, U.isSetScale() );
    }

  }
}
