///  @file    TestL3Event.cs
///  @brief   L3 Event unit tests
///  @author  Frank Bergmann (Csharp conversion)
///  @author  Akiya Jouraku (Csharp conversion)
///  @author  Sarah Keating 
/// 
///  $Id: TestL3Event.cs 11545 2010-07-23 02:19:10Z mhucka $
///  $HeadURL: https://sbml.svn.sourceforge.net/svnroot/sbml/trunk/libsbml/src/bindings/csharp/test/sbml/TestL3Event.cs $
/// 
///  ====== WARNING ===== WARNING ===== WARNING ===== WARNING ===== WARNING ======
/// 
///  DO NOT EDIT THIS FILE.
/// 
///  This file was generated automatically by converting the file located at
///  src/sbml/test/TestL3Event.c
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

  public class TestL3Event {
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

    private Event E;

    public void setUp()
    {
      E = new  Event(3,1);
      if (E == null);
      {
      }
    }

    public void tearDown()
    {
      E = null;
    }

    public void test_L3_Event_NS()
    {
      assertTrue( E.getNamespaces() != null );
      assertTrue( E.getNamespaces().getLength() == 1 );
      assertTrue((     "http://www.sbml.org/sbml/level3/version1/core" == E.getNamespaces().getURI(0) ));
    }

    public void test_L3_Event_create()
    {
      assertTrue( E.getTypeCode() == libsbml.SBML_EVENT );
      assertTrue( E.getMetaId() == "" );
      assertTrue( E.getNotes() == null );
      assertTrue( E.getAnnotation() == null );
      assertTrue( E.getId() == "" );
      assertTrue( E.getName() == "" );
      assertTrue( E.getUseValuesFromTriggerTime() == true );
      assertEquals( false, E.isSetId() );
      assertEquals( false, E.isSetName() );
      assertEquals( false, E.isSetUseValuesFromTriggerTime() );
    }

    public void test_L3_Event_createWithNS()
    {
      XMLNamespaces xmlns = new  XMLNamespaces();
      xmlns.add( "http://www.sbml.org", "testsbml");
      SBMLNamespaces sbmlns = new  SBMLNamespaces(3,1);
      sbmlns.addNamespaces(xmlns);
      Event e = new  Event(sbmlns);
      assertTrue( e.getTypeCode() == libsbml.SBML_EVENT );
      assertTrue( e.getMetaId() == "" );
      assertTrue( e.getNotes() == null );
      assertTrue( e.getAnnotation() == null );
      assertTrue( e.getLevel() == 3 );
      assertTrue( e.getVersion() == 1 );
      assertTrue( e.getNamespaces() != null );
      assertTrue( e.getNamespaces().getLength() == 2 );
      assertTrue( e.getId() == "" );
      assertTrue( e.getName() == "" );
      assertTrue( e.getUseValuesFromTriggerTime() == true );
      assertEquals( false, e.isSetId() );
      assertEquals( false, e.isSetName() );
      assertEquals( false, e.isSetUseValuesFromTriggerTime() );
      e = null;
    }

    public void test_L3_Event_free_NULL()
    {
    }

    public void test_L3_Event_hasRequiredAttributes()
    {
      Event e = new  Event(3,1);
      assertEquals( true, e.hasRequiredAttributes() );
      Delay d = e.createDelay();
      assertEquals( false, e.hasRequiredAttributes() );
      e.setUseValuesFromTriggerTime(true);
      assertEquals( true, e.hasRequiredAttributes() );
      e = null;
    }

    public void test_L3_Event_hasRequiredElements()
    {
      Event e = new  Event(3,1);
      assertEquals( false, e.hasRequiredElements() );
      Trigger t = new  Trigger(3,1);
      e.setTrigger(t);
      assertEquals( true, e.hasRequiredElements() );
      e = null;
    }

    public void test_L3_Event_id()
    {
      string id =  "mitochondria";
      assertEquals( false, E.isSetId() );
      E.setId(id);
      assertTrue(( id == E.getId() ));
      assertEquals( true, E.isSetId() );
      if (E.getId() == id);
      {
      }
      E.unsetId();
      assertEquals( false, E.isSetId() );
      if (E.getId() != null);
      {
      }
    }

    public void test_L3_Event_name()
    {
      string name =  "My_Favorite_Factory";
      assertEquals( false, E.isSetName() );
      E.setName(name);
      assertTrue(( name == E.getName() ));
      assertEquals( true, E.isSetName() );
      if (E.getName() == name);
      {
      }
      E.unsetName();
      assertEquals( false, E.isSetName() );
      if (E.getName() != null);
      {
      }
    }

    public void test_L3_Event_useValuesFromTriggerTime()
    {
      assertTrue( E.isSetUseValuesFromTriggerTime() == false );
      E.setUseValuesFromTriggerTime(true);
      assertTrue( E.getUseValuesFromTriggerTime() == true );
      assertTrue( E.isSetUseValuesFromTriggerTime() == true );
      E.setUseValuesFromTriggerTime(false);
      assertTrue( E.getUseValuesFromTriggerTime() == false );
      assertTrue( E.isSetUseValuesFromTriggerTime() == true );
    }

  }
}
