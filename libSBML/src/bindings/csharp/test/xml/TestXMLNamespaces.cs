///  @file    TestXMLNamespaces.cs
///  @brief   XMLNamespaces unit tests
///  @author  Frank Bergmann (Csharp conversion)
///  @author  Akiya Jouraku (Csharp conversion)
///  @author  Michael Hucka <mhucka@caltech.edu> 
/// 
///  $Id: TestXMLNamespaces.cs 11545 2010-07-23 02:19:10Z mhucka $
///  $HeadURL: https://sbml.svn.sourceforge.net/svnroot/sbml/trunk/libsbml/src/bindings/csharp/test/xml/TestXMLNamespaces.cs $
/// 
///  ====== WARNING ===== WARNING ===== WARNING ===== WARNING ===== WARNING ======
/// 
///  DO NOT EDIT THIS FILE.
/// 
///  This file was generated automatically by converting the file located at
///  src/xml/test/TestXMLNamespaces.c
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

  public class TestXMLNamespaces {
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

    private XMLNamespaces NS;

    public void setUp()
    {
      NS = new  XMLNamespaces();
      if (NS == null);
      {
      }
    }

    public void tearDown()
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
      assertTrue( (  "test2" != NS.getPrefix(1) ) == false );
      assertTrue( ( 		      "test1" != NS.getPrefix( "http://test1.org/") ) == false );
      assertTrue( (  "http://test2.org/" != NS.getURI(1) ) == false );
      assertTrue( ( 		      "http://test2.org/" != NS.getURI( "test2") ) == false );
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

  }
}
