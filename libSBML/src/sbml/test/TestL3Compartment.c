/**
 * \file    TestL3Compartment.c
 * \brief   L3 Compartment unit tests
 * \author  Sarah Keating
 *
 * $Id: TestL3Compartment.c 10396 2009-12-03 13:36:00Z sarahkeating $
 * $HeadURL: https://sbml.svn.sourceforge.net/svnroot/sbml/trunk/libsbml/src/sbml/test/TestL3Compartment.c $
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
#include <sbml/Compartment.h>
#include <sbml/xml/XMLNamespaces.h>
#include <sbml/SBMLDocument.h>

#include <check.h>


static Compartment_t *C;


void
L3CompartmentTest_setup (void)
{
  C = Compartment_create(3, 1);

  if (C == NULL)
  {
    fail("Compartment_create(3, 1) returned a NULL pointer.");
  }
}


void
L3CompartmentTest_teardown (void)
{
  Compartment_free(C);
}


START_TEST (test_L3_Compartment_create)
{
  fail_unless( SBase_getTypeCode  ((SBase_t *) C) == SBML_COMPARTMENT );
  fail_unless( SBase_getMetaId    ((SBase_t *) C) == NULL );
  fail_unless( SBase_getNotes     ((SBase_t *) C) == NULL );
  fail_unless( SBase_getAnnotation((SBase_t *) C) == NULL );

  fail_unless( Compartment_getId     (C) == NULL );
  fail_unless( Compartment_getName   (C) == NULL );
  fail_unless( Compartment_getUnits  (C) == NULL );
  fail_unless( Compartment_getOutside(C) == NULL );

  fail_unless( isnan(Compartment_getSpatialDimensionsAsDouble(C))   );
  fail_unless( isnan(Compartment_getVolume(C)));
  fail_unless( Compartment_getConstant(C) == 1   );

  fail_unless( !Compartment_isSetId     (C) );
  fail_unless( !Compartment_isSetSpatialDimensions (C) );
  fail_unless( !Compartment_isSetName   (C) );
  fail_unless( !Compartment_isSetSize   (C) );
  fail_unless( !Compartment_isSetVolume (C) );
  fail_unless( !Compartment_isSetUnits  (C) );
  fail_unless( !Compartment_isSetOutside(C) );
  fail_unless( !Compartment_isSetConstant(C) );
}
END_TEST


START_TEST (test_L3_Compartment_free_NULL)
{
  Compartment_free(NULL);
}
END_TEST


START_TEST (test_L3_Compartment_id)
{
  char *id = "mitochondria";


  fail_unless( !Compartment_isSetId(C) );
  
  Compartment_setId(C, id);

  fail_unless( !strcmp(Compartment_getId(C), id) );
  fail_unless( Compartment_isSetId(C) );

  if (Compartment_getId(C) == id)
  {
    fail("Compartment_setId(...) did not make a copy of string.");
  }
}
END_TEST


START_TEST (test_L3_Compartment_name)
{
  char *name = "My_Favorite_Factory";


  fail_unless( !Compartment_isSetName(C) );

  Compartment_setName(C, name);

  fail_unless( !strcmp(Compartment_getName(C), name) );
  fail_unless( Compartment_isSetName(C) );

  if (Compartment_getName(C) == name)
  {
    fail("Compartment_setName(...) did not make a copy of string.");
  }

  Compartment_unsetName(C);
  
  fail_unless( !Compartment_isSetName(C) );

  if (Compartment_getName(C) != NULL)
  {
    fail("Compartment_unsetName(C) did not clear string.");
  }
}
END_TEST


START_TEST (test_L3_Compartment_units)
{
  char *units = "volume";


  fail_unless( !Compartment_isSetUnits(C) );
  
  Compartment_setUnits(C, units);

  fail_unless( !strcmp(Compartment_getUnits(C), units) );
  fail_unless( Compartment_isSetUnits(C) );

  if (Compartment_getUnits(C) == units)
  {
    fail("Compartment_setUnits(...) did not make a copy of string.");
  }

  Compartment_unsetUnits(C);
  
  fail_unless( !Compartment_isSetUnits(C) );

  if (Compartment_getUnits(C) != NULL)
  {
    fail("Compartment_unsetUnits(C, NULL) did not clear string.");
  }
}
END_TEST


START_TEST (test_L3_Compartment_size)
{
  double size = 0.2;

  fail_unless( !Compartment_isSetSize(C));
  fail_unless( isnan(Compartment_getSize(C)));
  
  Compartment_setSize(C, size);

  fail_unless( Compartment_getSize(C) == size );
  fail_unless( Compartment_isSetSize(C) );

  Compartment_unsetSize(C);

  fail_unless( !Compartment_isSetSize(C) );
  fail_unless( isnan(Compartment_getSize(C)));
}
END_TEST


START_TEST (test_L3_Compartment_constant)
{
  fail_unless(Compartment_isSetConstant(C) == 0);

  Compartment_setConstant(C, 1);

  fail_unless(Compartment_getConstant(C) == 1);
  fail_unless(Compartment_isSetConstant(C) == 1);

  Compartment_setConstant(C, 0);

  fail_unless(Compartment_getConstant(C) == 0);
  fail_unless(Compartment_isSetConstant(C) == 1);

}
END_TEST


START_TEST (test_L3_Compartment_spatialDimensions)
{
  fail_unless( !Compartment_isSetSpatialDimensions(C));
  fail_unless( isnan(Compartment_getSpatialDimensionsAsDouble(C)));

  Compartment_setSpatialDimensionsAsDouble(C, 1.5);

  fail_unless( Compartment_isSetSpatialDimensions(C));
  fail_unless( Compartment_getSpatialDimensionsAsDouble(C) == 1.5);

  Compartment_unsetSpatialDimensions(C);

  fail_unless( !Compartment_isSetSpatialDimensions(C));
  fail_unless( isnan(Compartment_getSpatialDimensionsAsDouble(C)));
}
END_TEST


START_TEST (test_L3_Compartment_createWithNS )
{
  XMLNamespaces_t *xmlns = XMLNamespaces_create();
  XMLNamespaces_add(xmlns, "http://www.sbml.org", "testsbml");
  SBMLNamespaces_t *sbmlns = SBMLNamespaces_create(3,1);
  SBMLNamespaces_addNamespaces(sbmlns,xmlns);

  Compartment_t *c = 
    Compartment_createWithNS (sbmlns);


  fail_unless( SBase_getTypeCode  ((SBase_t *) c) == SBML_COMPARTMENT );
  fail_unless( SBase_getMetaId    ((SBase_t *) c) == NULL );
  fail_unless( SBase_getNotes     ((SBase_t *) c) == NULL );
  fail_unless( SBase_getAnnotation((SBase_t *) c) == NULL );

  fail_unless( SBase_getLevel       ((SBase_t *) c) == 3 );
  fail_unless( SBase_getVersion     ((SBase_t *) c) == 1 );

  fail_unless( Compartment_getNamespaces     (c) != NULL );
  fail_unless( XMLNamespaces_getLength(Compartment_getNamespaces(c)) == 2 );


  fail_unless( Compartment_getId     (c) == NULL );
  fail_unless( Compartment_getName   (c) == NULL );
  fail_unless( Compartment_getUnits  (c) == NULL );
  fail_unless( Compartment_getOutside(c) == NULL );

  fail_unless( isnan(Compartment_getSpatialDimensionsAsDouble(c))   );
  fail_unless( isnan(Compartment_getVolume(c)));
  fail_unless( Compartment_getConstant(c) == 1   );

  fail_unless( !Compartment_isSetId     (c) );
  fail_unless( !Compartment_isSetSpatialDimensions (c) );
  fail_unless( !Compartment_isSetName   (c) );
  fail_unless( !Compartment_isSetSize   (c) );
  fail_unless( !Compartment_isSetVolume (c) );
  fail_unless( !Compartment_isSetUnits  (c) );
  fail_unless( !Compartment_isSetOutside(c) );
  fail_unless( !Compartment_isSetConstant(c) );

  Compartment_free(c);
}
END_TEST


START_TEST (test_L3_Compartment_hasRequiredAttributes )
{
  Compartment_t *c = Compartment_create (3, 1);

  fail_unless ( !Compartment_hasRequiredAttributes(c));

  Compartment_setId(c, "id");

  fail_unless ( !Compartment_hasRequiredAttributes(c));

  Compartment_setConstant(c, 0);

  fail_unless ( Compartment_hasRequiredAttributes(c));

  Compartment_free(c);
}
END_TEST


START_TEST (test_L3_Compartment_NS)
{
  fail_unless( Compartment_getNamespaces     (C) != NULL );
  fail_unless( XMLNamespaces_getLength(Compartment_getNamespaces(C)) == 1 );
  fail_unless( !strcmp( XMLNamespaces_getURI(Compartment_getNamespaces(C), 0),
    "http://www.sbml.org/sbml/level3/version1/core"));
}
END_TEST


Suite *
create_suite_L3_Compartment (void)
{
  Suite *suite = suite_create("L3_Compartment");
  TCase *tcase = tcase_create("L3_Compartment");


  tcase_add_checked_fixture( tcase,
                             L3CompartmentTest_setup,
                             L3CompartmentTest_teardown );

  tcase_add_test( tcase, test_L3_Compartment_create              );
  tcase_add_test( tcase, test_L3_Compartment_free_NULL           );
  tcase_add_test( tcase, test_L3_Compartment_id               );
  tcase_add_test( tcase, test_L3_Compartment_name             );
  tcase_add_test( tcase, test_L3_Compartment_units            );
  tcase_add_test( tcase, test_L3_Compartment_size           );
  tcase_add_test( tcase, test_L3_Compartment_constant      );
  tcase_add_test( tcase, test_L3_Compartment_spatialDimensions);
  tcase_add_test( tcase, test_L3_Compartment_createWithNS         );
  tcase_add_test( tcase, test_L3_Compartment_hasRequiredAttributes        );
  tcase_add_test( tcase, test_L3_Compartment_NS              );

  suite_add_tcase(suite, tcase);

  return suite;
}
