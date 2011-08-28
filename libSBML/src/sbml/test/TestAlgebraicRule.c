/**
 * \file    TestAlgebraicRule.c
 * \brief   AlgebraicRule unit tests
 * \author  Ben Bornstein
 *
 * $Id: TestAlgebraicRule.c 10129 2009-08-28 12:23:22Z sarahkeating $
 * $HeadURL: https://sbml.svn.sourceforge.net/svnroot/sbml/trunk/libsbml/src/sbml/test/TestAlgebraicRule.c $
 */
/* Copyright 2002 California Institute of Technology and
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

#include <sbml/math/FormulaFormatter.h>
#include <sbml/math/FormulaParser.h>

#include <sbml/SBase.h>
#include <sbml/Rule.h>
#include <sbml/xml/XMLNamespaces.h>
#include <sbml/SBMLDocument.h>

#include <check.h>


static Rule_t *AR;


void
AlgebraicRuleTest_setup (void)
{
  AR = Rule_createAlgebraic(2, 4);

  if (AR == NULL)
  {
    fail("AlgebraicRule_createAlgebraic() returned a NULL pointer.");
  }
}


void
AlgebraicRuleTest_teardown (void)
{
  Rule_free(AR);
}


START_TEST (test_AlgebraicRule_create)
{
  fail_unless( SBase_getTypeCode  ((SBase_t *) AR) == SBML_ALGEBRAIC_RULE );
  fail_unless( SBase_getMetaId    ((SBase_t *) AR) == NULL );
  fail_unless( SBase_getNotes     ((SBase_t *) AR) == NULL );
  fail_unless( SBase_getAnnotation((SBase_t *) AR) == NULL );

  fail_unless( Rule_getFormula((Rule_t *) AR) == NULL );
  fail_unless( Rule_getMath   ((Rule_t *) AR) == NULL );
}
END_TEST


START_TEST (test_AlgebraicRule_createWithFormula)
{
  const ASTNode_t *math;
  char *formula;

  Rule_t *ar = Rule_createAlgebraic(2, 4);
  Rule_setFormula(ar, "1 + 1");


  fail_unless( SBase_getTypeCode  ((SBase_t *) ar) == SBML_ALGEBRAIC_RULE );
  fail_unless( SBase_getMetaId    ((SBase_t *) ar) == NULL );

  math = Rule_getMath((Rule_t *) ar);
  fail_unless(math != NULL);

  formula = SBML_formulaToString(math);
  fail_unless( formula != NULL );
  fail_unless( !strcmp(formula, "1 + 1") );

  fail_unless( !strcmp(Rule_getFormula((Rule_t *) ar), formula) );

  Rule_free(ar);
  safe_free(formula);
}
END_TEST


START_TEST (test_AlgebraicRule_createWithMath)
{
  ASTNode_t       *math = SBML_parseFormula("1 + 1");
  Rule_t *ar   = Rule_createAlgebraic(2, 4);
  Rule_setMath(ar, math);


  fail_unless( SBase_getTypeCode  ((SBase_t *) ar) == SBML_ALGEBRAIC_RULE );
  fail_unless( SBase_getMetaId    ((SBase_t *) ar) == NULL );

  fail_unless( !strcmp(Rule_getFormula((Rule_t *) ar), "1 + 1") );
  fail_unless( Rule_getMath((Rule_t *) ar) != math );

  Rule_free(ar);
}
END_TEST


START_TEST (test_AlgebraicRule_free_NULL)
{
  Rule_free(NULL);
}
END_TEST


START_TEST (test_AlgebraicRule_createWithNS )
{
  XMLNamespaces_t *xmlns = XMLNamespaces_create();
  XMLNamespaces_add(xmlns, "http://www.sbml.org", "testsbml");
  SBMLNamespaces_t *sbmlns = SBMLNamespaces_create(2,3);
  SBMLNamespaces_addNamespaces(sbmlns,xmlns);

  Rule_t *r = 
    Rule_createAlgebraicWithNS(sbmlns);


  fail_unless( SBase_getTypeCode  ((SBase_t *) r) == SBML_ALGEBRAIC_RULE );
  fail_unless( SBase_getMetaId    ((SBase_t *) r) == NULL );
  fail_unless( SBase_getNotes     ((SBase_t *) r) == NULL );
  fail_unless( SBase_getAnnotation((SBase_t *) r) == NULL );

  fail_unless( SBase_getLevel       ((SBase_t *) r) == 2 );
  fail_unless( SBase_getVersion     ((SBase_t *) r) == 3 );

  fail_unless( Rule_getNamespaces     (r) != NULL );
  fail_unless( XMLNamespaces_getLength(Rule_getNamespaces(r)) == 2 );


  Rule_free(r);
}
END_TEST


Suite *
create_suite_AlgebraicRule (void)
{
  Suite *suite = suite_create("AlgebraicRule");
  TCase *tcase = tcase_create("AlgebraicRule");


  tcase_add_checked_fixture( tcase,
                             AlgebraicRuleTest_setup,
                             AlgebraicRuleTest_teardown );

  tcase_add_test( tcase, test_AlgebraicRule_create         );
  tcase_add_test( tcase, test_AlgebraicRule_createWithFormula     );
  tcase_add_test( tcase, test_AlgebraicRule_createWithMath );
  tcase_add_test( tcase, test_AlgebraicRule_free_NULL      );
  tcase_add_test( tcase, test_AlgebraicRule_createWithNS         );

  suite_add_tcase(suite, tcase);

  return suite;
}
