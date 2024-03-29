///  @file    TestReadFromFile2.cs
///  @brief   Tests for reading MathML from files into ASTNodes.
///  @author  Frank Bergmann (Csharp conversion)
///  @author  Akiya Jouraku (Csharp conversion)
///  @author  Sarah Keating 
/// 
///  $Id: TestReadFromFile2.cs 11545 2010-07-23 02:19:10Z mhucka $
///  $HeadURL: https://sbml.svn.sourceforge.net/svnroot/sbml/trunk/libsbml/src/bindings/csharp/test/math/TestReadFromFile2.cs $
/// 
///  ====== WARNING ===== WARNING ===== WARNING ===== WARNING ===== WARNING ======
/// 
///  DO NOT EDIT THIS FILE.
/// 
///  This file was generated automatically by converting the file located at
///  src/math/test/TestReadFromFile2.cpp
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

  public class TestReadFromFile2 {
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


    public void test_read_MathML_2()
    {
      SBMLReader reader = new SBMLReader();
      SBMLDocument d;
      Model m;
      FunctionDefinition fd;
      InitialAssignment ia;
      Rule r;
      string filename =  "../../math/test/test-data/";
      filename += "mathML_2.xml";
      d = reader.readSBML(filename);
      if (d == null);
      {
      }
      m = d.getModel();
      assertTrue( m != null );
      assertTrue( m.getNumFunctionDefinitions() == 2 );
      assertTrue( m.getNumInitialAssignments() == 1 );
      assertTrue( m.getNumRules() == 2 );
      fd = m.getFunctionDefinition(0);
      ASTNode fd_math = fd.getMath();
      assertTrue( fd_math.getType() == libsbml.AST_LAMBDA );
      assertTrue( fd_math.getNumChildren() == 1 );
      assertTrue((  "lambda()" == libsbml.formulaToString(fd_math) ));
      ASTNode child = fd_math.getChild(0);
      assertTrue( child.getType() == libsbml.AST_UNKNOWN );
      assertTrue( child.getNumChildren() == 0 );
      assertTrue((  "" == libsbml.formulaToString(child) ));
      fd = m.getFunctionDefinition(1);
      ASTNode fd1_math = fd.getMath();
      assertTrue( fd1_math.getType() == libsbml.AST_LAMBDA );
      assertTrue( fd1_math.getNumChildren() == 2 );
      assertTrue((                            "lambda(x, piecewise(p, leq(x, 4)))" == libsbml.formulaToString(fd1_math) ));
      ASTNode child1 = fd1_math.getRightChild();
      assertTrue( child1.getType() == libsbml.AST_FUNCTION_PIECEWISE );
      assertTrue( child1.getNumChildren() == 2 );
      assertTrue((                                      "piecewise(p, leq(x, 4))" == libsbml.formulaToString(child1) ));
      ASTNode c1 = child1.getChild(0);
      assertTrue( c1.getType() == libsbml.AST_NAME );
      assertTrue( c1.getNumChildren() == 0 );
      assertTrue((  "p" == libsbml.formulaToString(c1) ));
      ASTNode c2 = child1.getChild(1);
      assertTrue( c2.getType() == libsbml.AST_RELATIONAL_LEQ );
      assertTrue( c2.getNumChildren() == 2 );
      assertTrue((  "leq(x, 4)" == libsbml.formulaToString(c2) ));
      ia = m.getInitialAssignment(0);
      ASTNode ia_math = ia.getMath();
      assertTrue( ia_math.getType() == libsbml.AST_FUNCTION_PIECEWISE );
      assertTrue( ia_math.getNumChildren() == 4 );
      assertTrue((                      "piecewise(-x, lt(x, 0), 0, eq(x, 0))" == libsbml.formulaToString(ia_math) ));
      child1 = ia_math.getChild(0);
      ASTNode child2 = ia_math.getChild(1);
      ASTNode child3 = ia_math.getChild(2);
      ASTNode child4 = ia_math.getChild(3);
      assertTrue( child1.getType() == libsbml.AST_MINUS );
      assertTrue( child1.getNumChildren() == 1 );
      assertTrue((  "-x" == libsbml.formulaToString(child1) ));
      assertTrue( child2.getType() == libsbml.AST_RELATIONAL_LT );
      assertTrue( child2.getNumChildren() == 2 );
      assertTrue((  "lt(x, 0)" == libsbml.formulaToString(child2) ));
      assertTrue( child3.getType() == libsbml.AST_REAL );
      assertTrue( child3.getNumChildren() == 0 );
      assertTrue((  "0" == libsbml.formulaToString(child3) ));
      assertTrue( child4.getType() == libsbml.AST_RELATIONAL_EQ );
      assertTrue( child4.getNumChildren() == 2 );
      assertTrue((  "eq(x, 0)" == libsbml.formulaToString(child4) ));
      r = m.getRule(0);
      ASTNode r_math = r.getMath();
      assertTrue( r_math.getType() == libsbml.AST_CONSTANT_TRUE );
      assertTrue( r_math.getNumChildren() == 0 );
      assertTrue((  "true" == libsbml.formulaToString(r_math) ));
      r = m.getRule(1);
      ASTNode r1_math = r.getMath();
      assertTrue( r1_math.getType() == libsbml.AST_FUNCTION_LOG );
      assertTrue( r1_math.getNumChildren() == 2 );
      assertTrue((  "log(3, x)" == libsbml.formulaToString(r1_math) ));
      child1 = r1_math.getChild(0);
      child2 = r1_math.getChild(1);
      assertTrue( child1.getType() == libsbml.AST_REAL );
      assertTrue( child1.getNumChildren() == 0 );
      assertTrue((  "3" == libsbml.formulaToString(child1) ));
      assertTrue( child2.getType() == libsbml.AST_NAME );
      assertTrue( child2.getNumChildren() == 0 );
      assertTrue((  "x" == libsbml.formulaToString(child2) ));
      d = null;
    }

  }
}
