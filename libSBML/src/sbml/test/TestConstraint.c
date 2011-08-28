/**
 * \file    TestConstraint.c
 * \brief   SBML Constraint unit tests
 * \author  Sarah Keating
 *
 * $Id: TestConstraint.c 10129 2009-08-28 12:23:22Z sarahkeating $
 * $HeadURL: https://sbml.svn.sourceforge.net/svnroot/sbml/trunk/libsbml/src/sbml/test/TestConstraint.c $
 */
/* Copyright 2003 California Institute of Technology and
 * Japan Science and Technology Corporation.
 *
 * This library is free software; you can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation; either version 2.1 of the License, or
 * any later version.
 *
 * This library is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY, WITHOUT EVEN THE IMPLIED WARRANTY OF
 * MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE.  The software and
 * documentation provided hereunder is on an "as is" basis, and the
 * California Institute of Technology and Japan Science and Technology
 * Corporation have no obligations to provide maintenance, support,
 * updates, enhancements or modifications.  In no event shall the
 * California Institute of Technology or the Japan Science and Technology
 * Corporation be liable to any party for direct, indirect, special,
 * incidental or consequential damages, including lost profits, arising
 * out of the use of this software and its documentation, even if the
 * California Institute of Technology and/or Japan Science and Technology
 * Corporation have been advised of the possibility of such damage.  See
 * the GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this library; if not, write to the Free Software Foundation,
 * Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA.
 *
 * The original code contained here was initially developed by:
 *
 *     Ben Bornstein
 *     The Systems Biology Markup Language Development Group
 *     ERATO Kitano Symbiotic Systems Project
 *     Control and Dynamical Systems, MC 107-81
 *     California Institute of Technology
 *     Pasadena, CA, 91125, USA
 *
 *     http://www.cds.caltech.edu/erato
 *     mailto:sbml-team@caltech.edu
 *
 * Contributor(s):
 */


#include <sbml/common/common.h>
#include <sbml/math/FormulaParser.h>

#include <sbml/SBase.h>
#include <sbml/Constraint.h>
#include <sbml/xml/XMLNamespaces.h>
#include <sbml/xml/XMLAttributes.h>
#include <sbml/xml/XMLTriple.h>
#include <sbml/xml/XMLNode.h>
#include <sbml/SBMLDocument.h>

#include <check.h>


static Constraint_t *C;


void
ConstraintTest_setup (void)
{
  C = Constraint_create(2, 4);

  if (C == NULL)
  {
    fail("Constraint_create() returned a NULL pointer.");
  }
}


void
ConstraintTest_teardown (void)
{
  Constraint_free(C);
}


START_TEST (test_Constraint_create)
{
  fail_unless( SBase_getTypeCode  ((SBase_t *) C) == SBML_CONSTRAINT );
  fail_unless( SBase_getMetaId    ((SBase_t *) C) == NULL );
  fail_unless( SBase_getNotes     ((SBase_t *) C) == NULL );
  fail_unless( SBase_getAnnotation((SBase_t *) C) == NULL );

  fail_unless( !Constraint_isSetMessage(C) );
  fail_unless( !Constraint_isSetMath    (C) );
}
END_TEST


START_TEST (test_Constraint_free_NULL)
{
  Constraint_free(NULL);
}
END_TEST


//START_TEST (test_Constraint_createWithMath)
//{
//  ASTNode_t       *math = SBML_parseFormula("1 + 1");
//  Constraint_t *c   = Constraint_createWithMath(math);
//
//
//
//  fail_unless( SBase_getTypeCode  ((SBase_t *) c) == SBML_CONSTRAINT );
//  fail_unless( SBase_getMetaId    ((SBase_t *) c) == NULL );
//
//  fail_unless( Constraint_getMath(c) != math );
//  fail_unless( !Constraint_isSetMessage(c) );
//  fail_unless( Constraint_isSetMath    (c) );
//  Constraint_free(c);
//}
//END_TEST


START_TEST (test_Constraint_setMath)
{
  ASTNode_t *math = SBML_parseFormula("2 * k");

  Constraint_setMath(C, math);

  fail_unless( Constraint_getMath(C) != math );
  fail_unless( Constraint_isSetMath(C) );

  /* Reflexive case (pathological) */
  Constraint_setMath(C, (ASTNode_t *) Constraint_getMath(C));

  fail_unless( Constraint_getMath(C) != math );

  Constraint_setMath(C, NULL);
  fail_unless( !Constraint_isSetMath(C) );

  if (Constraint_getMath(C) != NULL)
  {
    fail("Constraint_setMath(C, NULL) did not clear ASTNode.");
  }

  ASTNode_free(math);
}
END_TEST


START_TEST (test_Constraint_setMessage)
{
  const XMLNode_t *text = XMLNode_convertStringToXMLNode(" Some text ", NULL);
  XMLTriple_t *triple = XMLTriple_createWith("p", "http://www.w3.org/1999/xhtml", "");
  XMLAttributes_t *att = XMLAttributes_create();
  XMLNamespaces_t *xmlns = XMLNamespaces_create();
  XMLNamespaces_add(xmlns, "http://www.w3.org/1999/xhtml", "");
  
  XMLNode_t *p = XMLNode_createStartElementNS(triple, att, xmlns);
  XMLNode_addChild(p, text);
  
  XMLTriple_t *triple1 = XMLTriple_createWith("message", "", "");
  XMLAttributes_t *att1 = XMLAttributes_create();
  XMLNode_t *node = XMLNode_createStartElement(triple1, att1);

  XMLNode_addChild(node, p);

  Constraint_setMessage(C, node);

  fail_unless( Constraint_getMessage(C) != node );
  fail_unless( Constraint_isSetMessage(C) == 1);

  /* Reflexive case (pathological) */
  Constraint_setMessage(C, (XMLNode_t *) Constraint_getMessage(C));

  fail_unless( Constraint_getMessage(C) != node );

  fail_unless( Constraint_getMessageString(C) != NULL );

  Constraint_unsetMessage(C);
  fail_unless( !Constraint_isSetMessage(C) );

  if (Constraint_getMessage(C) != NULL)
  {
    fail("Constraint_unsetMessage(C) did not clear XMLNode.");
  }

  XMLNode_free(node);
}
END_TEST


START_TEST (test_Constraint_createWithNS )
{
  XMLNamespaces_t *xmlns = XMLNamespaces_create();
  XMLNamespaces_add(xmlns, "http://www.sbml.org", "testsbml");
  SBMLNamespaces_t *sbmlns = SBMLNamespaces_create(2,2);
  SBMLNamespaces_addNamespaces(sbmlns,xmlns);

  Constraint_t *object = 
    Constraint_createWithNS (sbmlns);


  fail_unless( SBase_getTypeCode  ((SBase_t *) object) == SBML_CONSTRAINT );
  fail_unless( SBase_getMetaId    ((SBase_t *) object) == NULL );
  fail_unless( SBase_getNotes     ((SBase_t *) object) == NULL );
  fail_unless( SBase_getAnnotation((SBase_t *) object) == NULL );

  fail_unless( SBase_getLevel       ((SBase_t *) object) == 2 );
  fail_unless( SBase_getVersion     ((SBase_t *) object) == 2 );

  fail_unless( Constraint_getNamespaces     (object) != NULL );
  fail_unless( XMLNamespaces_getLength(Constraint_getNamespaces(object)) == 2 );

  Constraint_free(object);
}
END_TEST


Suite *
create_suite_Constraint (void)
{
  Suite *suite = suite_create("Constraint");
  TCase *tcase = tcase_create("Constraint");


  tcase_add_checked_fixture( tcase,
                             ConstraintTest_setup,
                             ConstraintTest_teardown );

  tcase_add_test( tcase, test_Constraint_create      );
  //tcase_add_test( tcase, test_Constraint_createWithMath      );
  tcase_add_test( tcase, test_Constraint_free_NULL   );
  tcase_add_test( tcase, test_Constraint_setMath     );
  tcase_add_test( tcase, test_Constraint_setMessage  );
  tcase_add_test( tcase, test_Constraint_createWithNS         );

  suite_add_tcase(suite, tcase);

  return suite;
}
