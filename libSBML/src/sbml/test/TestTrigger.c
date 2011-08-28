/**
 * \file    TestTrigger.c
 * \brief   SBML Trigger unit tests
 * \author  Sarah Keating
 *
 * $Id: TestTrigger.c 10129 2009-08-28 12:23:22Z sarahkeating $
 * $HeadURL: https://sbml.svn.sourceforge.net/svnroot/sbml/trunk/libsbml/src/sbml/test/TestTrigger.c $
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
#include <sbml/math/FormulaFormatter.h>

#include <sbml/SBase.h>
#include <sbml/Trigger.h>
#include <sbml/xml/XMLNamespaces.h>
#include <sbml/SBMLDocument.h>

#include <check.h>


static Trigger_t *D;


void
TriggerTest_setup (void)
{
  D = Trigger_create(2, 4);

  if (D == NULL)
  {
    fail("Trigger_create() returned a NULL pointer.");
  }
}


void
TriggerTest_teardown (void)
{
  Trigger_free(D);
}


START_TEST (test_Trigger_create)
{
  fail_unless( SBase_getTypeCode((SBase_t *) D) == SBML_TRIGGER );
  fail_unless( SBase_getMetaId    ((SBase_t *) D) == NULL );
  fail_unless( SBase_getNotes     ((SBase_t *) D) == NULL );
  fail_unless( SBase_getAnnotation((SBase_t *) D) == NULL );

  fail_unless( Trigger_getMath(D) == NULL );
}
END_TEST


//START_TEST (test_Trigger_createWithMath)
//{
//  ASTNode_t            *math = SBML_parseFormula("x^3");
//  Trigger_t *fd   = Trigger_createWithMath(math);
//
//  const ASTNode_t * math1;
//  char * formula;
//
//  fail_unless( SBase_getTypeCode((SBase_t *) fd) == SBML_TRIGGER );
//  fail_unless( SBase_getMetaId    ((SBase_t *) fd) == NULL );
//  fail_unless( SBase_getNotes     ((SBase_t *) fd) == NULL );
//  fail_unless( SBase_getAnnotation((SBase_t *) fd) == NULL );
//
//
//  math1 = Trigger_getMath(fd);
//  fail_unless( math1 != NULL );
//
//  formula = SBML_formulaToString(math1);
//  fail_unless( formula != NULL );
//  fail_unless( !strcmp(formula, "x^3") );
//  fail_unless( Trigger_getMath(fd) != math );
//  fail_unless( Trigger_isSetMath(fd) );
//
//
//  Trigger_free(fd);
//}
//END_TEST


START_TEST (test_Trigger_free_NULL)
{
  Trigger_free(NULL);
}
END_TEST


START_TEST (test_Trigger_setMath)
{
  ASTNode_t *math = SBML_parseFormula("lambda(x, x^3)");

  const ASTNode_t * math1;
  char * formula;

  Trigger_setMath(D, math);

  math1 = Trigger_getMath(D);
  fail_unless( math1 != NULL );

  formula = SBML_formulaToString(math1);
  fail_unless( formula != NULL );
  fail_unless( !strcmp(formula, "lambda(x, x^3)") );
  fail_unless( Trigger_getMath(D) != math );
  fail_unless( Trigger_isSetMath(D) );

  /* Reflexive case (pathological) */
  Trigger_setMath(D, (ASTNode_t *) Trigger_getMath(D));
  math1 = Trigger_getMath(D);
  fail_unless( math1 != NULL );

  formula = SBML_formulaToString(math1);
  fail_unless( formula != NULL );
  fail_unless( !strcmp(formula, "lambda(x, x^3)") );

  Trigger_setMath(D, NULL);
  fail_unless( !Trigger_isSetMath(D) );

  if (Trigger_getMath(D) != NULL)
  {
    fail("Trigger_setMath(D, NULL) did not clear ASTNode.");
  }
}
END_TEST


START_TEST (test_Trigger_setMath1)
{
  ASTNode_t *math = SBML_parseFormula("2 * k");

  int i = Trigger_setMath(D, math);

  fail_unless( i == LIBSBML_OPERATION_SUCCESS);
  fail_unless( Trigger_getMath(D) != math );
  fail_unless( Trigger_isSetMath(D) );

  i = Trigger_setMath(D, NULL);
  
  fail_unless( i == LIBSBML_OPERATION_SUCCESS);
  fail_unless( Trigger_getMath(D) == NULL );
  fail_unless( !Trigger_isSetMath(D) );

  ASTNode_free(math);
}
END_TEST


START_TEST (test_Trigger_setMath2)
{
  ASTNode_t *math = ASTNode_createWithType(AST_TIMES);

  int i = Trigger_setMath(D, math);

  fail_unless( i == LIBSBML_INVALID_OBJECT);
  fail_unless( !Trigger_isSetMath(D) );

  ASTNode_free(math);
}
END_TEST


START_TEST (test_Trigger_createWithNS )
{
  XMLNamespaces_t *xmlns = XMLNamespaces_create();
  XMLNamespaces_add(xmlns, "http://www.sbml.org", "testsbml");
  SBMLNamespaces_t *sbmlns = SBMLNamespaces_create(2,1);
  SBMLNamespaces_addNamespaces(sbmlns,xmlns);

  Trigger_t *object = 
    Trigger_createWithNS (sbmlns);


  fail_unless( SBase_getTypeCode  ((SBase_t *) object) == SBML_TRIGGER );
  fail_unless( SBase_getMetaId    ((SBase_t *) object) == NULL );
  fail_unless( SBase_getNotes     ((SBase_t *) object) == NULL );
  fail_unless( SBase_getAnnotation((SBase_t *) object) == NULL );

  fail_unless( SBase_getLevel       ((SBase_t *) object) == 2 );
  fail_unless( SBase_getVersion     ((SBase_t *) object) == 1 );

  fail_unless( Trigger_getNamespaces     (object) != NULL );
  fail_unless( XMLNamespaces_getLength(Trigger_getNamespaces(object)) == 2 );

  Trigger_free(object);
}
END_TEST


Suite *
create_suite_Trigger (void)
{
  Suite *suite = suite_create("Trigger");
  TCase *tcase = tcase_create("Trigger");


  tcase_add_checked_fixture( tcase,
                             TriggerTest_setup,
                             TriggerTest_teardown );

  tcase_add_test( tcase, test_Trigger_create       );
  ////tcase_add_test( tcase, test_Trigger_createWithMath   );
  tcase_add_test( tcase, test_Trigger_setMath      );
  tcase_add_test( tcase, test_Trigger_setMath1     );
  tcase_add_test( tcase, test_Trigger_setMath2     );
  tcase_add_test( tcase, test_Trigger_free_NULL );
  tcase_add_test( tcase, test_Trigger_createWithNS         );

  suite_add_tcase(suite, tcase);

  return suite;
}
