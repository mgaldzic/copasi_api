/**
 * \file    TestL3Parameter.c
 * \brief   L3 Parameter unit tests
 * \author  Sarah Keating
 *
 * $Id: TestL3Parameter.c 10396 2009-12-03 13:36:00Z sarahkeating $
 * $HeadURL: https://sbml.svn.sourceforge.net/svnroot/sbml/trunk/libsbml/src/sbml/test/TestL3Parameter.c $
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
#include <sbml/Parameter.h>
#include <sbml/xml/XMLNamespaces.h>
#include <sbml/SBMLDocument.h>

#include <check.h>


static Parameter_t *P;


void
L3ParameterTest_setup (void)
{
  P = Parameter_create(3, 1);

  if (P == NULL)
  {
    fail("Parameter_create(3, 1) returned a NULL pointer.");
  }
}


void
L3ParameterTest_teardown (void)
{
  Parameter_free(P);
}


START_TEST (test_L3_Parameter_create)
{
  fail_unless( SBase_getTypeCode  ((SBase_t *) P) == SBML_PARAMETER );
  fail_unless( SBase_getMetaId    ((SBase_t *) P) == NULL );
  fail_unless( SBase_getNotes     ((SBase_t *) P) == NULL );
  fail_unless( SBase_getAnnotation((SBase_t *) P) == NULL );

  fail_unless( Parameter_getId     (P) == NULL );
  fail_unless( Parameter_getName   (P) == NULL );
  fail_unless( Parameter_getUnits  (P) == NULL );
  fail_unless( isnan(Parameter_getValue(P)));
  fail_unless( Parameter_getConstant(P) == 1   );

  fail_unless( !Parameter_isSetId     (P) );
  fail_unless( !Parameter_isSetName   (P) );
  fail_unless( !Parameter_isSetValue (P) );
  fail_unless( !Parameter_isSetUnits  (P) );
  fail_unless( !Parameter_isSetConstant(P) );
}
END_TEST


START_TEST (test_L3_Parameter_free_NULL)
{
  Parameter_free(NULL);
}
END_TEST


START_TEST (test_L3_Parameter_id)
{
  char *id = "mitochondria";


  fail_unless( !Parameter_isSetId(P) );
  
  Parameter_setId(P, id);

  fail_unless( !strcmp(Parameter_getId(P), id) );
  fail_unless( Parameter_isSetId(P) );

  if (Parameter_getId(P) == id)
  {
    fail("Parameter_setId(...) did not make a copy of string.");
  }
}
END_TEST


START_TEST (test_L3_Parameter_name)
{
  char *name = "My_Favorite_Factory";


  fail_unless( !Parameter_isSetName(P) );

  Parameter_setName(P, name);

  fail_unless( !strcmp(Parameter_getName(P), name) );
  fail_unless( Parameter_isSetName(P) );

  if (Parameter_getName(P) == name)
  {
    fail("Parameter_setName(...) did not make a copy of string.");
  }

  Parameter_unsetName(P);
  
  fail_unless( !Parameter_isSetName(P) );

  if (Parameter_getName(P) != NULL)
  {
    fail("Parameter_unsetName(P) did not clear string.");
  }
}
END_TEST


START_TEST (test_L3_Parameter_units)
{
  char *units = "volume";


  fail_unless( !Parameter_isSetUnits(P) );
  
  Parameter_setUnits(P, units);

  fail_unless( !strcmp(Parameter_getUnits(P), units) );
  fail_unless( Parameter_isSetUnits(P) );

  if (Parameter_getUnits(P) == units)
  {
    fail("Parameter_setUnits(...) did not make a copy of string.");
  }

  Parameter_unsetUnits(P);
  
  fail_unless( !Parameter_isSetUnits(P) );

  if (Parameter_getUnits(P) != NULL)
  {
    fail("Parameter_unsetUnits(P, NULL) did not clear string.");
  }
}
END_TEST


START_TEST (test_L3_Parameter_constant)
{
  fail_unless(Parameter_isSetConstant(P) == 0);

  Parameter_setConstant(P, 1);

  fail_unless(Parameter_getConstant(P) == 1);
  fail_unless(Parameter_isSetConstant(P) == 1);

  Parameter_setConstant(P, 0);

  fail_unless(Parameter_getConstant(P) == 0);
  fail_unless(Parameter_isSetConstant(P) == 1);

}
END_TEST


START_TEST (test_L3_Parameter_value)
{
  fail_unless( !Parameter_isSetValue(P));
  fail_unless( isnan(Parameter_getValue(P)));

  Parameter_setValue(P, 1.5);

  fail_unless( Parameter_isSetValue(P));
  fail_unless( Parameter_getValue(P) == 1.5);

  Parameter_unsetValue(P);

  fail_unless( !Parameter_isSetValue(P));
  fail_unless( isnan(Parameter_getValue(P)));
}
END_TEST


START_TEST (test_L3_Parameter_createWithNS )
{
  XMLNamespaces_t *xmlns = XMLNamespaces_create();
  XMLNamespaces_add(xmlns, "http://www.sbml.org", "testsbml");
  SBMLNamespaces_t *sbmlns = SBMLNamespaces_create(3,1);
  SBMLNamespaces_addNamespaces(sbmlns,xmlns);

  Parameter_t *p = 
    Parameter_createWithNS (sbmlns);


  fail_unless( SBase_getTypeCode  ((SBase_t *) p) == SBML_PARAMETER );
  fail_unless( SBase_getMetaId    ((SBase_t *) p) == NULL );
  fail_unless( SBase_getNotes     ((SBase_t *) p) == NULL );
  fail_unless( SBase_getAnnotation((SBase_t *) p) == NULL );

  fail_unless( SBase_getLevel       ((SBase_t *) p) == 3 );
  fail_unless( SBase_getVersion     ((SBase_t *) p) == 1 );

  fail_unless( Parameter_getNamespaces     (p) != NULL );
  fail_unless( XMLNamespaces_getLength(Parameter_getNamespaces(p)) == 2 );


  fail_unless( Parameter_getId     (p) == NULL );
  fail_unless( Parameter_getName   (p) == NULL );
  fail_unless( Parameter_getUnits  (p) == NULL );
  fail_unless( isnan(Parameter_getValue(p)));
  fail_unless( Parameter_getConstant(p) == 1   );

  fail_unless( !Parameter_isSetId     (p) );
  fail_unless( !Parameter_isSetName   (p) );
  fail_unless( !Parameter_isSetValue (p) );
  fail_unless( !Parameter_isSetUnits  (p) );
  fail_unless( !Parameter_isSetConstant(p) );

  Parameter_free(p);
}
END_TEST


START_TEST (test_L3_Parameter_hasRequiredAttributes )
{
  Parameter_t *p = Parameter_create (3, 1);

  fail_unless ( !Parameter_hasRequiredAttributes(p));

  Parameter_setId(p, "id");

  fail_unless ( !Parameter_hasRequiredAttributes(p));

  Parameter_setConstant(p, 0);

  fail_unless ( Parameter_hasRequiredAttributes(p));

  Parameter_free(p);
}
END_TEST


START_TEST (test_L3_Parameter_NS)
{
  fail_unless( Parameter_getNamespaces     (P) != NULL );
  fail_unless( XMLNamespaces_getLength(Parameter_getNamespaces(P)) == 1 );
  fail_unless( !strcmp( XMLNamespaces_getURI(Parameter_getNamespaces(P), 0),
    "http://www.sbml.org/sbml/level3/version1/core"));
}
END_TEST


Suite *
create_suite_L3_Parameter (void)
{
  Suite *suite = suite_create("L3_Parameter");
  TCase *tcase = tcase_create("L3_Parameter");


  tcase_add_checked_fixture( tcase,
                             L3ParameterTest_setup,
                             L3ParameterTest_teardown );

  tcase_add_test( tcase, test_L3_Parameter_create              );
  tcase_add_test( tcase, test_L3_Parameter_free_NULL           );
  tcase_add_test( tcase, test_L3_Parameter_id               );
  tcase_add_test( tcase, test_L3_Parameter_name             );
  tcase_add_test( tcase, test_L3_Parameter_units            );
  tcase_add_test( tcase, test_L3_Parameter_constant      );
  tcase_add_test( tcase, test_L3_Parameter_value);
  tcase_add_test( tcase, test_L3_Parameter_createWithNS         );
  tcase_add_test( tcase, test_L3_Parameter_hasRequiredAttributes        );
  tcase_add_test( tcase, test_L3_Parameter_NS              );

  suite_add_tcase(suite, tcase);

  return suite;
}
