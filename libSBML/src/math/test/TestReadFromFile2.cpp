/**
 * @file    TestReadFromFile2.cpp
 * @brief   Tests for reading MathML from files into ASTNodes.
 * @author  Sarah Keating
 *
 * $Id: TestReadFromFile2.cpp 10866 2010-01-29 19:52:27Z mhucka $
 * $HeadURL: https://sbml.svn.sourceforge.net/svnroot/sbml/trunk/libsbml/src/math/test/TestReadFromFile2.cpp $
 *
 *<!---------------------------------------------------------------------------
 * This file is part of libSBML.  Please visit http://sbml.org for more
 * information about SBML, and the latest version of libSBML.
 *
 * Copyright 2005-2010 California Institute of Technology.
 * Copyright 2002-2005 California Institute of Technology and
 *                     Japan Science and Technology Corporation.
 * 
 * This library is free software; you can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation.  A copy of the license agreement is provided
 * in the file named "LICENSE.txt" included with this software distribution and
 * also available online as http://sbml.org/software/libsbml/license.html
 *----------------------------------------------------------------------- -->*/


#include <sbml/common/common.h>

#include <sbml/SBMLReader.h>
#include <sbml/SBMLTypes.h>

#include <sbml/math/ASTNode.h>



#include <string>

#include <check.h>

LIBSBML_CPP_NAMESPACE_USE

BEGIN_C_DECLS


extern char *TestDataDirectory;


START_TEST (test_read_MathML_2)
{
  SBMLReader         reader;
  SBMLDocument*      d;
  Model*             m;
  FunctionDefinition* fd;
  InitialAssignment* ia;
  Rule*              r;


  std::string filename(TestDataDirectory);
  filename += "mathML_2.xml";


  d = reader.readSBML(filename);

  if (d == NULL)
  {
    fail("readSBML(\"mathML_2.xml\") returned a NULL pointer.");
  }

  m = d->getModel();
  fail_unless( m != NULL, NULL );

  // check that whole model has been read in
  fail_unless( m->getNumFunctionDefinitions() == 2, NULL);
  fail_unless( m->getNumInitialAssignments() == 1, NULL);
  fail_unless( m->getNumRules() == 2, NULL );

  //<functionDefinition id="fd">
  //  <math xmlns="http://www.w3.org/1998/Math/MathML">
  //    <lambda>
  //      <apply/>
  //    </lambda>
  //  </math>
  //</functionDefinition>
  fd = m->getFunctionDefinition(0);
  const ASTNode *fd_math = fd->getMath();

  fail_unless (fd_math->getType() == AST_LAMBDA, NULL);
  fail_unless (fd_math->getNumChildren() == 1, NULL);
  fail_unless (!strcmp(SBML_formulaToString(fd_math), "lambda()"), NULL);

  ASTNode *child = fd_math->getChild(0);
  fail_unless (child->getType() == AST_UNKNOWN, NULL);
  fail_unless (child->getNumChildren() == 0, NULL);
  fail_unless (!strcmp(SBML_formulaToString(child), ""), NULL);

  //<functionDefinition id="fd1">
  //  <math xmlns="http://www.w3.org/1998/Math/MathML">
  //    <lambda>
  //      <bvar>
  //        <ci> x </ci>
  //      </bvar>
        //<piecewise>
        //  <piece>
        //    <ci> p </ci>
        //    <apply>
        //      <leq/>
        //      <ci> x </ci>
        //      <cn type="integer"> 4 </cn>
        //    </apply>
        //  </piece>
        //</piecewise>
  //    </lambda>
  //  </math>
  //</functionDefinition>
  fd = m->getFunctionDefinition(1);
  const ASTNode *fd1_math = fd->getMath();

  fail_unless (fd1_math->getType() == AST_LAMBDA, NULL);
  fail_unless (fd1_math->getNumChildren() == 2, NULL);
  fail_unless (!strcmp(SBML_formulaToString(fd1_math), 
                          "lambda(x, piecewise(p, leq(x, 4)))"), NULL);

  ASTNode *child1 = fd1_math->getRightChild();
  fail_unless (child1->getType() == AST_FUNCTION_PIECEWISE, NULL);
  fail_unless (child1->getNumChildren() == 2, NULL);
  fail_unless (!strcmp(SBML_formulaToString(child1), 
                                    "piecewise(p, leq(x, 4))"), NULL);

  ASTNode *c1 = child1->getChild(0);
  fail_unless (c1->getType() == AST_NAME, NULL);
  fail_unless (c1->getNumChildren() == 0, NULL);
  fail_unless (!strcmp(SBML_formulaToString(c1), "p"), NULL);

  ASTNode *c2 = child1->getChild(1);
  fail_unless (c2->getType() == AST_RELATIONAL_LEQ, NULL);
  fail_unless (c2->getNumChildren() == 2, NULL);
  fail_unless (!strcmp(SBML_formulaToString(c2), "leq(x, 4)"), NULL);


  
  //<initialAssignment symbol="p1">
    //<math xmlns="http://www.w3.org/1998/Math/MathML">
    //    <piecewise>
    //      <piece>
    //          <apply><minus/><ci> x </ci></apply>
    //          <apply><lt/><ci> x </ci> <cn> 0 </cn></apply>
    //      </piece>
    //      <piece>
    //          <cn> 0 </cn>
    //          <apply><eq/><ci> x </ci> <cn> 0 </cn></apply>
    //      </piece>
    //    </piecewise>
    //</math>
  //</initialAssignment>
  ia = m->getInitialAssignment(0);
  const ASTNode *ia_math = ia->getMath();

  fail_unless (ia_math->getType() == AST_FUNCTION_PIECEWISE, NULL);
  fail_unless (ia_math->getNumChildren() == 4, NULL);
  fail_unless (!strcmp(SBML_formulaToString(ia_math), 
                    "piecewise(-x, lt(x, 0), 0, eq(x, 0))"), NULL);

  child1 = ia_math->getChild(0);
  ASTNode *child2 = ia_math->getChild(1);
  ASTNode *child3 = ia_math->getChild(2);
  ASTNode *child4 = ia_math->getChild(3);

  fail_unless (child1->getType() == AST_MINUS, NULL);
  fail_unless (child1->getNumChildren() == 1, NULL);
  fail_unless (!strcmp(SBML_formulaToString(child1), "-x"), NULL);

  fail_unless (child2->getType() == AST_RELATIONAL_LT, NULL);
  fail_unless (child2->getNumChildren() == 2, NULL);
  fail_unless (!strcmp(SBML_formulaToString(child2), "lt(x, 0)"), NULL);

  fail_unless (child3->getType() == AST_REAL, NULL);
  fail_unless (child3->getNumChildren() == 0, NULL);
  fail_unless (!strcmp(SBML_formulaToString(child3), "0"), NULL);

  fail_unless (child4->getType() == AST_RELATIONAL_EQ, NULL);
  fail_unless (child4->getNumChildren() == 2, NULL);
  fail_unless (!strcmp(SBML_formulaToString(child4), "eq(x, 0)"), NULL);

  //<algebraicRule>
  //  <math xmlns="http://www.w3.org/1998/Math/MathML">
      //<apply>
      //<true/>
      //</apply>
  //  </math>
  //</algebraicRule>
  r = m->getRule(0);
  const ASTNode *r_math = r->getMath();

  fail_unless (r_math->getType() == AST_CONSTANT_TRUE, NULL);
  fail_unless (r_math->getNumChildren() == 0, NULL);
  fail_unless (!strcmp(SBML_formulaToString(r_math), "true"), NULL);

  //<assignmentRule variable="p2">
  //  <math xmlns="http://www.w3.org/1998/Math/MathML">
      //<apply>
      //  <log/>
      //  <logbase>
      //    <cn> 3 </cn>
      //  </logbase>
      //  <ci> x </ci>
      //</apply>
  //  </math>
  //</assignmentRule>
  r = m->getRule(1);
  const ASTNode *r1_math = r->getMath();

  fail_unless (r1_math->getType() == AST_FUNCTION_LOG, NULL);
  fail_unless (r1_math->getNumChildren() == 2, NULL);
  fail_unless (!strcmp(SBML_formulaToString(r1_math), "log(3, x)"), NULL);

  child1 = r1_math->getChild(0);
  child2 = r1_math->getChild(1);

  fail_unless (child1->getType() == AST_REAL, NULL);
  fail_unless (child1->getNumChildren() == 0, NULL);
  fail_unless (!strcmp(SBML_formulaToString(child1), "3"), NULL);

  fail_unless (child2->getType() == AST_NAME, NULL);
  fail_unless (child2->getNumChildren() == 0, NULL);
  fail_unless (!strcmp(SBML_formulaToString(child2), "x"), NULL);


  delete d;
}
END_TEST


Suite *
create_suite_TestReadFromFile2 (void)
{ 
  Suite *suite = suite_create("test-data/mathML_2.xml");
  TCase *tcase = tcase_create("test-data/mathML_2.xml");


  tcase_add_test(tcase, test_read_MathML_2);

  suite_add_tcase(suite, tcase);

  return suite;
}


END_C_DECLS
