///  @file    TestListOf.cs
///  @brief   ListOf unit tests
///  @author  Frank Bergmann (Csharp conversion)
///  @author  Akiya Jouraku (Csharp conversion)
///  @author  Ben Bornstein 
/// 
///  $Id: TestListOf.cs 11545 2010-07-23 02:19:10Z mhucka $
///  $HeadURL: https://sbml.svn.sourceforge.net/svnroot/sbml/trunk/libsbml/src/bindings/csharp/test/sbml/TestListOf.cs $
/// 
///  ====== WARNING ===== WARNING ===== WARNING ===== WARNING ===== WARNING ======
/// 
///  DO NOT EDIT THIS FILE.
/// 
///  This file was generated automatically by converting the file located at
///  src/sbml/test/TestListOf.c
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

  public class TestListOf {
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
      assertTrue( lo.getMetaId() == "" );
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
      assertTrue( sp != elem );
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

  }
}
