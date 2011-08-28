/**
 * \file    TestListOf.c
 * \brief   ListOf unit tests
 * \author  Ben Bornstein
 *
 * $Id: TestListOf.c 11538 2010-07-22 23:24:40Z mhucka $
 * $HeadURL: https://sbml.svn.sourceforge.net/svnroot/sbml/trunk/libsbml/src/sbml/test/TestListOf.c $
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

#include <sbml/ListOf.h>
#include <sbml/SBase.h>
#include <sbml/Species.h>

#include <check.h>


START_TEST (test_ListOf_create)
{
  ListOf_t *lo = (ListOf_t*) ListOf_create();


  fail_unless( SBase_getTypeCode  ((SBase_t *) lo) == SBML_LIST_OF );
  fail_unless( SBase_getNotes     ((SBase_t *) lo) == NULL );
  fail_unless( SBase_getAnnotation((SBase_t *) lo) == NULL );
  fail_unless( SBase_getMetaId    ((SBase_t *) lo) == NULL );

  fail_unless( ListOf_size(lo) == 0 );

  ListOf_free(lo);
}
END_TEST


START_TEST (test_ListOf_free_NULL)
{
  ListOf_free(NULL);
}
END_TEST


START_TEST (test_ListOf_get)
{
  ListOf_t *lo = (ListOf_t*) ListOf_create();

  fail_unless( ListOf_size(lo) == 0 );

  SBase_t *sp = (SBase_t*)Species_create(2, 4);

  ListOf_append(lo, sp);

  fail_unless( ListOf_size(lo) == 1 );

  SBase_t *elem = ListOf_get(lo, 1);

  fail_unless( sp != elem );             /* ListOf_append makes a clone */

  Species_free((Species_t *) sp);

  ListOf_free(lo);
}
END_TEST


START_TEST (test_ListOf_remove)
{
  ListOf_t *lo = (ListOf_t*) ListOf_create();

  SBase_t *sp = (SBase_t*)Species_create(2, 4);

  fail_unless( ListOf_size(lo) == 0 );

  ListOf_append(lo, sp);
  ListOf_append(lo, sp);
  ListOf_append(lo, sp);
  ListOf_append(lo, sp);
  ListOf_append(lo, sp);

  fail_unless( ListOf_size(lo) == 5 );

  SBase_t *elem;
  elem = ListOf_remove(lo, 0);
  Species_free((Species_t*) elem);
  elem = ListOf_remove(lo, 0);
  Species_free((Species_t*) elem);
  elem = ListOf_remove(lo, 0);
  Species_free((Species_t*) elem);
  elem = ListOf_remove(lo, 0);
  Species_free((Species_t*) elem);
  elem = ListOf_remove(lo, 0);
  Species_free((Species_t*) elem);

  fail_unless( ListOf_size(lo) == 0 );

  ListOf_append(lo, sp);
  ListOf_append(lo, sp);
  ListOf_append(lo, sp);
  ListOf_append(lo, sp);
  ListOf_appendAndOwn(lo, sp);

  fail_unless( ListOf_size(lo) == 5 );

  ListOf_free(lo);
}
END_TEST


START_TEST (test_ListOf_clear)
{
  ListOf_t *lo = (ListOf_t*) ListOf_create();

  SBase_t *sp = (SBase_t*)Species_create(2, 4); 

  ListOf_append(lo, sp);
  ListOf_append(lo, sp);
  ListOf_append(lo, sp);
  ListOf_append(lo, sp);
  ListOf_append(lo, sp);

  fail_unless( ListOf_size(lo) == 5 );

  /* clear and delete */

  ListOf_clear(lo, 1);

  fail_unless( ListOf_size(lo) == 0 );

  ListOf_append(lo, sp);
  ListOf_append(lo, sp);
  ListOf_append(lo, sp);
  ListOf_append(lo, sp);
  ListOf_appendAndOwn(lo, sp);

  fail_unless( ListOf_size(lo) == 5 );

  /* delete each item */

  SBase_t *elem;
  elem = ListOf_get(lo, 0);
  Species_free((Species_t*) elem);
  elem = ListOf_get(lo, 1);
  Species_free((Species_t*) elem);
  elem = ListOf_get(lo, 2);
  Species_free((Species_t*) elem);
  elem = ListOf_get(lo, 3);
  Species_free((Species_t*) elem);
  elem = ListOf_get(lo, 4);
  Species_free((Species_t*) elem);

  /* clear only */

  ListOf_clear(lo, 0);

  fail_unless( ListOf_size(lo) == 0 );

  ListOf_free(lo);
  
}
END_TEST



Suite *
create_suite_ListOf (void) 
{ 
  Suite *suite = suite_create("ListOf");
  TCase *tcase = tcase_create("ListOf");
 

  tcase_add_test(tcase, test_ListOf_create    );
  tcase_add_test(tcase, test_ListOf_free_NULL );
  tcase_add_test(tcase, test_ListOf_get       );
  tcase_add_test(tcase, test_ListOf_remove    );
  tcase_add_test(tcase, test_ListOf_clear     );

  suite_add_tcase(suite, tcase);

  return suite;
}
