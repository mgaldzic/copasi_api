///  @file    TestReadFromFile1.cs
///  @brief   Tests for reading MathML from files into ASTNodes.
///  @author  Frank Bergmann (Csharp conversion)
///  @author  Akiya Jouraku (Csharp conversion)
///  @author  Sarah Keating 
/// 
///  $Id: TestReadFromFile1.cs 11545 2010-07-23 02:19:10Z mhucka $
///  $HeadURL: https://sbml.svn.sourceforge.net/svnroot/sbml/trunk/libsbml/src/bindings/csharp/test/math/TestReadFromFile1.cs $
/// 
///  ====== WARNING ===== WARNING ===== WARNING ===== WARNING ===== WARNING ======
/// 
///  DO NOT EDIT THIS FILE.
/// 
///  This file was generated automatically by converting the file located at
///  src/math/test/TestReadFromFile1.cpp
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

  public class TestReadFromFile1 {
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


    public void test_read_MathML_1()
    {
      SBMLReader reader = new SBMLReader();
      SBMLDocument d;
      Model m;
      FunctionDefinition fd;
      InitialAssignment ia;
      Rule r;
      KineticLaw kl;
      string filename =  "../../math/test/test-data/";
      filename += "mathML_1.xml";
      d = reader.readSBML(filename);
      if (d == null);
      {
      }
      m = d.getModel();
      assertTrue( m != null );
      assertTrue( m.getNumFunctionDefinitions() == 2 );
      assertTrue( m.getNumInitialAssignments() == 1 );
      assertTrue( m.getNumRules() == 2 );
      assertTrue( m.getNumReactions() == 1 );
      fd = m.getFunctionDefinition(0);
      ASTNode fd_math = fd.getMath();
      assertTrue( fd_math.getType() == libsbml.AST_LAMBDA );
      assertTrue( fd_math.getNumChildren() == 2 );
      assertTrue((  "lambda(x, )" == libsbml.formulaToString(fd_math) ));
      assertTrue( fd_math.getParentSBMLObject() == fd );
      ASTNode child = fd_math.getRightChild();
      assertTrue( child.getType() == libsbml.AST_UNKNOWN );
      assertTrue( child.getNumChildren() == 0 );
      assertTrue((  "" == libsbml.formulaToString(child) ));
      fd = m.getFunctionDefinition(1);
      ASTNode fd1_math = fd.getMath();
      assertTrue( fd1_math.getType() == libsbml.AST_LAMBDA );
      assertTrue( fd1_math.getNumChildren() == 2 );
      assertTrue((  "lambda(x, true)" == libsbml.formulaToString(fd1_math) ));
      assertTrue( fd1_math.getParentSBMLObject() == fd );
      ASTNode child1 = fd1_math.getRightChild();
      assertTrue( child1.getType() == libsbml.AST_CONSTANT_TRUE );
      assertTrue( child1.getNumChildren() == 0 );
      assertTrue((  "true" == libsbml.formulaToString(child1) ));
      ia = m.getInitialAssignment(0);
      ASTNode ia_math = ia.getMath();
      assertTrue( ia_math.getType() == libsbml.AST_UNKNOWN );
      assertTrue( ia_math.getNumChildren() == 0 );
      assertTrue((  "" == libsbml.formulaToString(ia_math) ));
      assertTrue( ia_math.getParentSBMLObject() == ia );
      r = m.getRule(0);
      ASTNode r_math = r.getMath();
      assertTrue( r_math.getType() == libsbml.AST_CONSTANT_TRUE );
      assertTrue( r_math.getNumChildren() == 0 );
      assertTrue((  "true" == libsbml.formulaToString(r_math) ));
      assertTrue( r_math.getParentSBMLObject() == r );
      r = m.getRule(1);
      ASTNode r1_math = r.getMath();
      assertTrue( r1_math.getType() == libsbml.AST_REAL );
      assertTrue( r1_math.getNumChildren() == 0 );
      assertTrue((  "INF" == libsbml.formulaToString(r1_math) ));
      assertTrue( r1_math.getParentSBMLObject() == r );
      kl = m.getReaction(0).getKineticLaw();
      ASTNode kl_math = kl.getMath();
      assertTrue( kl_math.getType() == libsbml.AST_REAL );
      assertTrue( kl_math.getNumChildren() == 0 );
      assertTrue((  "4.5" == libsbml.formulaToString(kl_math) ));
      assertTrue( kl_math.getParentSBMLObject() == kl );
      d = null;
    }

  }
}
