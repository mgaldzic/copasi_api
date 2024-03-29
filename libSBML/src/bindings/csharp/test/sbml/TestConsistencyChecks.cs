///  @file    TestConsistencyChecks.cs
///  @brief   Reads test-data/inconsistent.xml into memory and tests it.
///  @author  Frank Bergmann (Csharp conversion)
///  @author  Akiya Jouraku (Csharp conversion)
///  @author  Sarah Keating 
/// 
///  $Id: TestConsistencyChecks.cs 11545 2010-07-23 02:19:10Z mhucka $
///  $HeadURL: https://sbml.svn.sourceforge.net/svnroot/sbml/trunk/libsbml/src/bindings/csharp/test/sbml/TestConsistencyChecks.cs $
/// 
///  ====== WARNING ===== WARNING ===== WARNING ===== WARNING ===== WARNING ======
/// 
///  DO NOT EDIT THIS FILE.
/// 
///  This file was generated automatically by converting the file located at
///  src/sbml/test/TestConsistencyChecks.cpp
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

  public class TestConsistencyChecks {
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


    public void test_consistency_checks()
    {
      SBMLReader reader = new SBMLReader();
      SBMLDocument d;
      long errors;
      string filename =  "../../sbml/test/test-data/";
      filename += "inconsistent.xml";
      d = reader.readSBML(filename);
      if (d == null);
      {
      }
      errors = d.checkConsistency();
      assertTrue( errors == 1 );
      assertTrue( d.getError(0).getErrorId() == 10301 );
      d.getErrorLog().clearLog();
      d.setConsistencyChecks(libsbml.LIBSBML_CAT_IDENTIFIER_CONSISTENCY,false);
      errors = d.checkConsistency();
      assertTrue( errors == 1 );
      assertTrue( d.getError(0).getErrorId() == 20612 );
      d.getErrorLog().clearLog();
      d.setConsistencyChecks(libsbml.LIBSBML_CAT_GENERAL_CONSISTENCY,false);
      errors = d.checkConsistency();
      assertTrue( errors == 1 );
      assertTrue( d.getError(0).getErrorId() == 10701 );
      d.getErrorLog().clearLog();
      d.setConsistencyChecks(libsbml.LIBSBML_CAT_SBO_CONSISTENCY,false);
      errors = d.checkConsistency();
      assertTrue( errors == 1 );
      assertTrue( d.getError(0).getErrorId() == 10214 );
      d.getErrorLog().clearLog();
      d.setConsistencyChecks(libsbml.LIBSBML_CAT_MATHML_CONSISTENCY,false);
      errors = d.checkConsistency();
      assertTrue( errors == 3 );
      assertTrue( d.getError(0).getErrorId() == 99505 );
      assertTrue( d.getError(1).getErrorId() == 99505 );
      assertTrue( d.getError(2).getErrorId() == 80701 );
      d.getErrorLog().clearLog();
      d.setConsistencyChecks(libsbml.LIBSBML_CAT_UNITS_CONSISTENCY,false);
      errors = d.checkConsistency();
      assertTrue( errors == 0 );
      d = null;
    }

  }
}
