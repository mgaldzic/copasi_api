/**
 * \file    TestModifierSpeciesReference.c
 * \brief   ModifierSpeciesReference unit tests
 * \author  Ben Bornstein
 *
 * $Id: TestModifierSpeciesReference.c 11381 2010-07-05 22:45:53Z mhucka $
 * $HeadURL: https://sbml.svn.sourceforge.net/svnroot/sbml/trunk/libsbml/src/sbml/test/TestModifierSpeciesReference.c $
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
#include <sbml/SpeciesReference.h>
#include <sbml/xml/XMLNamespaces.h>
#include <sbml/SBMLDocument.h>

#include <check.h>


static SpeciesReference_t *MSR;


void
ModifierSpeciesReferenceTest_setup (void)
{
  MSR = SpeciesReference_createModifier(2, 4);

  if (MSR == NULL)
  {
    fail( "SpeciesReference_createModifier() returned a NULL pointer." );
  }
}


void
ModifierSpeciesReferenceTest_teardown (void)
{
  SpeciesReference_free(MSR);
}


START_TEST (test_ModifierSpeciesReference_create)
{
  fail_unless(SBase_getTypeCode((SBase_t *) MSR) ==
              SBML_MODIFIER_SPECIES_REFERENCE);

  fail_unless( SBase_getMetaId    ((SBase_t *) MSR) == NULL );
  fail_unless( SBase_getNotes     ((SBase_t *) MSR) == NULL );
  fail_unless( SBase_getAnnotation((SBase_t *) MSR) == NULL );

  fail_unless( SpeciesReference_getSpecies(MSR) == NULL );

  fail_unless( !SpeciesReference_isSetSpecies(MSR) );

  fail_unless( SpeciesReference_isModifier(MSR) );
}
END_TEST


START_TEST (test_ModifierSpeciesReference_free_NULL)
{
  SpeciesReference_free(NULL);
}
END_TEST


START_TEST (test_ModifierSpeciesReference_setSpecies)
{
  const char *s;
  char *species = "s1";



  SpeciesReference_setSpecies(MSR, species);

  s = SpeciesReference_getSpecies(MSR);
  fail_unless( !strcmp(s, species) );
  fail_unless(SpeciesReference_isSetSpecies(MSR));

  if (SpeciesReference_getSpecies(MSR) == species)
  {
    fail( "ModifierSpeciesReference_setSpecies(...) "
          "did not make a copy of string." );
  }

  /* Reflexive case (pathological) */
  s = SpeciesReference_getSpecies(MSR);
  SpeciesReference_setSpecies(MSR, s);

  s = SpeciesReference_getSpecies(MSR);
  fail_unless( !strcmp(s, species) );

  SpeciesReference_setSpecies(MSR, NULL);
  fail_unless(!SpeciesReference_isSetSpecies(MSR));

  if (SpeciesReference_getSpecies(MSR) != NULL)
  {
    fail( "ModifierSpeciesReference_setSpecies(MSR, NULL) "
          "did not clear string." );
  }
}
END_TEST


START_TEST (test_ModifierSpeciesReference_createWithNS )
{
  XMLNamespaces_t *xmlns = XMLNamespaces_create();
  XMLNamespaces_add(xmlns, "http://www.sbml.org", "testsbml");
  SBMLNamespaces_t *sbmlns = SBMLNamespaces_create(2,1);
  SBMLNamespaces_addNamespaces(sbmlns,xmlns);

  SBase_t *object = (SBase_t *) SpeciesReference_createModifierWithNS(sbmlns);

  fail_unless( SBase_getTypeCode  ((SBase_t *) object) == SBML_MODIFIER_SPECIES_REFERENCE );
  fail_unless( SBase_getMetaId    ((SBase_t *) object) == NULL );
  fail_unless( SBase_getNotes     ((SBase_t *) object) == NULL );
  fail_unless( SBase_getAnnotation((SBase_t *) object) == NULL );

  fail_unless( SBase_getLevel       ((SBase_t *) object) == 2 );
  fail_unless( SBase_getVersion     ((SBase_t *) object) == 1 );

  fail_unless( SpeciesReference_getNamespaces ((SpeciesReference_t *) object) != NULL );
  const XMLNamespaces_t *n = SpeciesReference_getNamespaces((SpeciesReference_t *) object);
  fail_unless( XMLNamespaces_getLength( n ) == 2 );

  SpeciesReference_free((SpeciesReference_t *) object);
}
END_TEST


Suite *
create_suite_ModifierSpeciesReference (void)
{
  Suite *suite = suite_create("ModifierSpeciesReference");
  TCase *tcase = tcase_create("ModifierSpeciesReference");


  tcase_add_checked_fixture( tcase,
                             ModifierSpeciesReferenceTest_setup,
                             ModifierSpeciesReferenceTest_teardown );

  tcase_add_test( tcase, test_ModifierSpeciesReference_create     );
  tcase_add_test( tcase, test_ModifierSpeciesReference_free_NULL  );
  tcase_add_test( tcase, test_ModifierSpeciesReference_setSpecies );
  tcase_add_test( tcase, test_ModifierSpeciesReference_createWithNS         );

  suite_add_tcase(suite, tcase);

  return suite;
}
