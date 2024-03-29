/**
 * \file    TestL3Event.c
 * \brief   L3 Event unit tests
 * \author  Sarah Keating
 *
 * $Id: TestL3Event.c 11430 2010-07-08 13:31:41Z sarahkeating $
 * $HeadURL: https://sbml.svn.sourceforge.net/svnroot/sbml/trunk/libsbml/src/sbml/test/TestL3Event.c $
 *
 *<!---------------------------------------------------------------------------
 * This file is part of libSBML.  Please visit http://sbml.org for more
 * information about SBML, and the latest version of libSBML.
 *
 * Copyright 2005-2007 California Institute of Technology.
 * Copyright 2002-2005 California Institute of Technology and
 *                     Japan Science and Technology Corporation.
 * 
 * This library is free software; you can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation.  A copy of the license agreement is provided
 * in the file named "LICENSE.txt" included with this software distribution
 * and also available online as http://sbml.org/software/libsbml/license.html
 *----------------------------------------------------------------------- -->*/


#include <sbml/common/common.h>

#include <sbml/SBase.h>
#include <sbml/Event.h>
#include <sbml/Trigger.h>
#include <sbml/xml/XMLNamespaces.h>
#include <sbml/SBMLDocument.h>

#include <check.h>


static Event_t *E;


void
L3EventTest_setup (void)
{
  E = Event_create(3, 1);

  if (E == NULL)
  {
    fail("Event_create(3, 1) returned a NULL pointer.");
  }
}


void
L3EventTest_teardown (void)
{
  Event_free(E);
}


START_TEST (test_L3_Event_create)
{
  fail_unless( SBase_getTypeCode  ((SBase_t *) E) == SBML_EVENT );
  fail_unless( SBase_getMetaId    ((SBase_t *) E) == NULL );
  fail_unless( SBase_getNotes     ((SBase_t *) E) == NULL );
  fail_unless( SBase_getAnnotation((SBase_t *) E) == NULL );

  fail_unless( Event_getId     (E) == NULL );
  fail_unless( Event_getName   (E) == NULL );
  fail_unless( Event_getUseValuesFromTriggerTime(E) == 1   );

  fail_unless( !Event_isSetId     (E) );
  fail_unless( !Event_isSetName   (E) );
  fail_unless( !Event_isSetUseValuesFromTriggerTime(E) );
}
END_TEST


START_TEST (test_L3_Event_free_NULL)
{
  Event_free(NULL);
}
END_TEST


START_TEST (test_L3_Event_id)
{
  char *id = "mitochondria";


  fail_unless( !Event_isSetId(E) );
  
  Event_setId(E, id);

  fail_unless( !strcmp(Event_getId(E), id) );
  fail_unless( Event_isSetId(E) );

  if (Event_getId(E) == id)
  {
    fail("Event_setId(...) did not make a copy of string.");
  }
 
  Event_unsetId(E);
  
  fail_unless( !Event_isSetId(E) );

  if (Event_getId(E) != NULL)
  {
    fail("Event_unsetId(E) did not clear string.");
  }
}
END_TEST


START_TEST (test_L3_Event_name)
{
  char *name = "My_Favorite_Factory";


  fail_unless( !Event_isSetName(E) );

  Event_setName(E, name);

  fail_unless( !strcmp(Event_getName(E), name) );
  fail_unless( Event_isSetName(E) );

  if (Event_getName(E) == name)
  {
    fail("Event_setName(...) did not make a copy of string.");
  }

  Event_unsetName(E);
  
  fail_unless( !Event_isSetName(E) );

  if (Event_getName(E) != NULL)
  {
    fail("Event_unsetName(E) did not clear string.");
  }
}
END_TEST


START_TEST (test_L3_Event_useValuesFromTriggerTime)
{
  fail_unless(Event_isSetUseValuesFromTriggerTime(E) == 0);

  Event_setUseValuesFromTriggerTime(E, 1);

  fail_unless(Event_getUseValuesFromTriggerTime(E) == 1);
  fail_unless(Event_isSetUseValuesFromTriggerTime(E) == 1);

  Event_setUseValuesFromTriggerTime(E, 0);

  fail_unless(Event_getUseValuesFromTriggerTime(E) == 0);
  fail_unless(Event_isSetUseValuesFromTriggerTime(E) == 1);

}
END_TEST


START_TEST (test_L3_Event_createWithNS )
{
  XMLNamespaces_t *xmlns = XMLNamespaces_create();
  XMLNamespaces_add(xmlns, "http://www.sbml.org", "testsbml");
  SBMLNamespaces_t *sbmlns = SBMLNamespaces_create(3,1);
  SBMLNamespaces_addNamespaces(sbmlns,xmlns);

  Event_t *e = 
    Event_createWithNS (sbmlns);


  fail_unless( SBase_getTypeCode  ((SBase_t *) e) == SBML_EVENT );
  fail_unless( SBase_getMetaId    ((SBase_t *) e) == NULL );
  fail_unless( SBase_getNotes     ((SBase_t *) e) == NULL );
  fail_unless( SBase_getAnnotation((SBase_t *) e) == NULL );

  fail_unless( SBase_getLevel       ((SBase_t *) e) == 3 );
  fail_unless( SBase_getVersion     ((SBase_t *) e) == 1 );

  fail_unless( Event_getNamespaces     (e) != NULL );
  fail_unless( XMLNamespaces_getLength(Event_getNamespaces(e)) == 2 );


  fail_unless( Event_getId     (e) == NULL );
  fail_unless( Event_getName   (e) == NULL );
  fail_unless( Event_getUseValuesFromTriggerTime(e) == 1   );

  fail_unless( !Event_isSetId     (e) );
  fail_unless( !Event_isSetName   (e) );
  fail_unless( !Event_isSetUseValuesFromTriggerTime(e) );

  Event_free(e);
}
END_TEST


START_TEST (test_L3_Event_hasRequiredAttributes )
{
  Event_t *e = Event_create (3, 1);

  fail_unless ( Event_hasRequiredAttributes(e));

  Delay_t *d = Event_createDelay(e);

  fail_unless ( !Event_hasRequiredAttributes(e));

  Event_setUseValuesFromTriggerTime(e, 1);

  fail_unless ( Event_hasRequiredAttributes(e));

  Event_free(e);
}
END_TEST


START_TEST (test_L3_Event_hasRequiredElements )
{
  Event_t *e = Event_create (3, 1);

  fail_unless ( !Event_hasRequiredElements(e));

  Trigger_t *t = Trigger_create(3, 1);
  Event_setTrigger(e, t);

  fail_unless ( Event_hasRequiredElements(e));

  Event_free(e);
}
END_TEST


START_TEST (test_L3_Event_NS)
{
  fail_unless( Event_getNamespaces     (E) != NULL );
  fail_unless( XMLNamespaces_getLength(Event_getNamespaces(E)) == 1 );
  fail_unless( !strcmp( XMLNamespaces_getURI(Event_getNamespaces(E), 0),
    "http://www.sbml.org/sbml/level3/version1/core"));
}
END_TEST


Suite *
create_suite_L3_Event (void)
{
  Suite *suite = suite_create("L3_Event");
  TCase *tcase = tcase_create("L3_Event");


  tcase_add_checked_fixture( tcase,
                             L3EventTest_setup,
                             L3EventTest_teardown );

  tcase_add_test( tcase, test_L3_Event_create              );
  tcase_add_test( tcase, test_L3_Event_free_NULL           );
  tcase_add_test( tcase, test_L3_Event_id               );
  tcase_add_test( tcase, test_L3_Event_name             );
  tcase_add_test( tcase, test_L3_Event_useValuesFromTriggerTime   );
  tcase_add_test( tcase, test_L3_Event_createWithNS         );
  tcase_add_test( tcase, test_L3_Event_hasRequiredAttributes        );
  tcase_add_test( tcase, test_L3_Event_hasRequiredElements        );
  tcase_add_test( tcase, test_L3_Event_NS              );

  suite_add_tcase(suite, tcase);

  return suite;
}
