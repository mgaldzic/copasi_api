///  @file    TestFunctionDefinition.cs
///  @brief   SBML FunctionDefinition unit tests
///  @author  Frank Bergmann (Csharp conversion)
///  @author  Akiya Jouraku (Csharp conversion)
///  @author  Ben Bornstein 
/// 
///  $Id: TestFunctionDefinition.cs 11545 2010-07-23 02:19:10Z mhucka $
///  $HeadURL: https://sbml.svn.sourceforge.net/svnroot/sbml/trunk/libsbml/src/bindings/csharp/test/sbml/TestFunctionDefinition.cs $
/// 
///  ====== WARNING ===== WARNING ===== WARNING ===== WARNING ===== WARNING ======
/// 
///  DO NOT EDIT THIS FILE.
/// 
///  This file was generated automatically by converting the file located at
///  src/sbml/test/TestFunctionDefinition.c
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

  public class TestFunctionDefinition {
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

    private FunctionDefinition FD;

    public void setUp()
    {
      FD = new  FunctionDefinition(2,4);
      if (FD == null);
      {
      }
    }

    public void tearDown()
    {
      FD = null;
    }

    public void test_FunctionDefinition_create()
    {
      assertTrue( FD.getTypeCode() == libsbml.SBML_FUNCTION_DEFINITION );
      assertTrue( FD.getMetaId() == "" );
      assertTrue( FD.getNotes() == null );
      assertTrue( FD.getAnnotation() == null );
      assertTrue( FD.getId() == "" );
      assertTrue( FD.getName() == "" );
      assertTrue( FD.getMath() == null );
    }

    public void test_FunctionDefinition_createWith()
    {
      ASTNode math = libsbml.parseFormula("lambda(x, x^3)");
      FunctionDefinition fd = new  FunctionDefinition(2,4);
      fd.setId( "pow3");
      fd.setMath(math);
      ASTNode math1;
      string formula;
      assertTrue( fd.getTypeCode() == libsbml.SBML_FUNCTION_DEFINITION );
      assertTrue( fd.getMetaId() == "" );
      assertTrue( fd.getNotes() == null );
      assertTrue( fd.getAnnotation() == null );
      assertTrue( fd.getName() == "" );
      math1 = fd.getMath();
      assertTrue( math1 != null );
      formula = libsbml.formulaToString(math1);
      assertTrue( formula != null );
      assertTrue((  "lambda(x, x^3)" == formula ));
      assertTrue( fd.getMath() != math );
      assertEquals( true, fd.isSetMath() );
      assertTrue((  "pow3" == fd.getId() ));
      assertEquals( true, fd.isSetId() );
      math = null;
      fd = null;
    }

    public void test_FunctionDefinition_createWithNS()
    {
      XMLNamespaces xmlns = new  XMLNamespaces();
      xmlns.add( "http://www.sbml.org", "testsbml");
      SBMLNamespaces sbmlns = new  SBMLNamespaces(2,1);
      sbmlns.addNamespaces(xmlns);
      FunctionDefinition object1 = new  FunctionDefinition(sbmlns);
      assertTrue( object1.getTypeCode() == libsbml.SBML_FUNCTION_DEFINITION );
      assertTrue( object1.getMetaId() == "" );
      assertTrue( object1.getNotes() == null );
      assertTrue( object1.getAnnotation() == null );
      assertTrue( object1.getLevel() == 2 );
      assertTrue( object1.getVersion() == 1 );
      assertTrue( object1.getNamespaces() != null );
      assertTrue( object1.getNamespaces().getLength() == 2 );
      object1 = null;
    }

    public void test_FunctionDefinition_free_NULL()
    {
    }

    public void test_FunctionDefinition_getArguments()
    {
      ASTNode math;
      FD.setMath(libsbml.parseFormula("lambda(x, y, x^y)"));
      assertTrue( FD.getNumArguments() == 2 );
      math = FD.getArgument(0);
      assertTrue( math != null );
      assertEquals( true, math.isName() );
      assertTrue((  "x" == math.getName() ));
      assertTrue( math.getNumChildren() == 0 );
      math = FD.getArgument(1);
      assertTrue( math != null );
      assertEquals( true, math.isName() );
      assertTrue((  "y" == math.getName() ));
      assertTrue( math.getNumChildren() == 0 );
      assertTrue( FD.getArgument(0) == FD.getArgument( "x") );
      assertTrue( FD.getArgument(1) == FD.getArgument( "y") );
    }

    public void test_FunctionDefinition_getBody()
    {
      ASTNode math;
      ASTNode math1 = libsbml.parseFormula("lambda(x, x)");
      FD.setMath(math1);
      math = FD.getBody();
      assertTrue( math != null );
      assertEquals( true, math.isName() );
      assertTrue((  "x" == math.getName() ));
      assertTrue( math.getNumChildren() == 0 );
      math1 = null;
    }

    public void test_FunctionDefinition_setId()
    {
      string id =  "pow3";
      FD.setId(id);
      assertTrue(( id == FD.getId() ));
      assertEquals( true, FD.isSetId() );
      if (FD.getId() == id);
      {
      }
      FD.setId(FD.getId());
      assertTrue(( id == FD.getId() ));
      FD.setId("");
      assertEquals( false, FD.isSetId() );
      if (FD.getId() != null);
      {
      }
    }

    public void test_FunctionDefinition_setMath()
    {
      ASTNode math = libsbml.parseFormula("lambda(x, x^3)");
      ASTNode math1;
      string formula;
      FD.setMath(math);
      math1 = FD.getMath();
      assertTrue( math1 != null );
      formula = libsbml.formulaToString(math1);
      assertTrue( formula != null );
      assertTrue((  "lambda(x, x^3)" == formula ));
      assertTrue( FD.getMath() != math );
      assertEquals( true, FD.isSetMath() );
      FD.setMath(FD.getMath());
      math1 = FD.getMath();
      assertTrue( math1 != null );
      formula = libsbml.formulaToString(math1);
      assertTrue( formula != null );
      assertTrue((  "lambda(x, x^3)" == formula ));
      assertTrue( FD.getMath() != math );
      FD.setMath(null);
      assertEquals( false, FD.isSetMath() );
      if (FD.getMath() != null);
      {
      }
    }

    public void test_FunctionDefinition_setName()
    {
      string name =  "Cube_Me";
      FD.setName(name);
      assertTrue(( name == FD.getName() ));
      assertEquals( true, FD.isSetName() );
      if (FD.getName() == name);
      {
      }
      FD.setName(FD.getName());
      assertTrue(( name == FD.getName() ));
      FD.setName("");
      assertEquals( false, FD.isSetName() );
      if (FD.getName() != null);
      {
      }
    }

  }
}
