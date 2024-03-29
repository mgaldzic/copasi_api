/**
 * \file    TestSpeciesReference_newSetters.c
 * \brief   SpeciesReference unit tests for new set function API
 * \author  Sarah Keating
 *
 * $Id: TestSpeciesReference_newSetters.c 11402 2010-07-07 01:43:53Z mhucka $
 * $HeadURL: https://sbml.svn.sourceforge.net/svnroot/sbml/trunk/libsbml/src/sbml/test/TestSpeciesReference_newSetters.c $
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
#include <sbml/SpeciesReference.h>
#include <sbml/xml/XMLNamespaces.h>
#include <sbml/SBMLDocument.h>
#include <sbml/StoichiometryMath.h>
#include <sbml/math/ASTNode.h>

#include <check.h>


static SpeciesReference_t *sr;


void
SpeciesReferenceTest1_setup (void)
{
  sr = SpeciesReference_create(2, 4);

  if (sr == NULL)
  {
    fail("SpeciesReference_create() returned a NULL pointer.");
  }
}


void
SpeciesReferenceTest1_teardown (void)
{
  SpeciesReference_free(sr);
}


START_TEST (test_SpeciesReference_setId1)
{
  int i = SpeciesReference_setId(sr, "cell");

  fail_unless( i == LIBSBML_OPERATION_SUCCESS );
  fail_unless( SpeciesReference_isSetId(sr) );
  fail_unless( !strcmp(SpeciesReference_getId(sr), "cell" ));
}
END_TEST


START_TEST (test_SpeciesReference_setId2)
{
  int i = SpeciesReference_setId(sr, "1cell");

  fail_unless( i == LIBSBML_INVALID_ATTRIBUTE_VALUE );
  fail_unless( !SpeciesReference_isSetId(sr) );
}
END_TEST


START_TEST (test_SpeciesReference_setId3)
{
  SpeciesReference_t *c = 
    SpeciesReference_create(2, 1);

  int i = SpeciesReference_setId(c, "cell");

  /* this is going to give different results dependenet on whether
   * layout is enabled but the test code
   */
#ifdef USE_LAYOUT
  fail_unless( i == LIBSBML_OPERATION_SUCCESS );
  fail_unless( SpeciesReference_isSetId(c) );
#else
  fail_unless( i == LIBSBML_UNEXPECTED_ATTRIBUTE );
  fail_unless( !SpeciesReference_isSetId(c) );
#endif
  SpeciesReference_free(c);
}
END_TEST


START_TEST (test_SpeciesReference_setId4)
{
  int i = SpeciesReference_setId(sr, "cell");

  fail_unless( i == LIBSBML_OPERATION_SUCCESS );
  fail_unless( SpeciesReference_isSetId(sr) );
  fail_unless( !strcmp(SpeciesReference_getId(sr), "cell" ));

  i = SpeciesReference_setId(sr, NULL);

  fail_unless( i == LIBSBML_OPERATION_SUCCESS );
  fail_unless( !SpeciesReference_isSetId(sr) );
}
END_TEST


START_TEST (test_SpeciesReference_setName1)
{
  int i = SpeciesReference_setName(sr, "cell");

  fail_unless( i == LIBSBML_OPERATION_SUCCESS );
  fail_unless( SpeciesReference_isSetName(sr) );

  i = SpeciesReference_unsetName(sr);

  fail_unless( i == LIBSBML_OPERATION_SUCCESS );
  fail_unless( !SpeciesReference_isSetName(sr) );
}
END_TEST


START_TEST (test_SpeciesReference_setName2)
{
  int i = SpeciesReference_setName(sr, "1cell");

  fail_unless( i == LIBSBML_INVALID_ATTRIBUTE_VALUE );
  fail_unless( !SpeciesReference_isSetName(sr) );

  i = SpeciesReference_unsetName(sr);

  fail_unless( i == LIBSBML_OPERATION_SUCCESS );
  fail_unless( !SpeciesReference_isSetName(sr) );
}
END_TEST


START_TEST (test_SpeciesReference_setName3)
{
  SpeciesReference_t *c = 
    SpeciesReference_create(2, 1);

  int i = SpeciesReference_setName(c, "cell");

  fail_unless( i == LIBSBML_UNEXPECTED_ATTRIBUTE );
  fail_unless( !SpeciesReference_isSetName(c) );

  SpeciesReference_free(c);
}
END_TEST


START_TEST (test_SpeciesReference_setName4)
{
  int i = SpeciesReference_setName(sr, "cell");

  fail_unless( i == LIBSBML_OPERATION_SUCCESS );
  fail_unless( SpeciesReference_isSetName(sr) );

  i = SpeciesReference_setName(sr, NULL);

  fail_unless( i == LIBSBML_OPERATION_SUCCESS );
  fail_unless( !SpeciesReference_isSetName(sr) );
}
END_TEST


START_TEST (test_SpeciesReference_setSpecies1)
{
  int i = SpeciesReference_setSpecies(sr, "mm");

  fail_unless( i == LIBSBML_OPERATION_SUCCESS );
  fail_unless( SpeciesReference_isSetSpecies(sr) );

}
END_TEST


START_TEST (test_SpeciesReference_setSpecies2)
{
  SpeciesReference_t *c = 
    SpeciesReference_create(2, 2);

  int i = SpeciesReference_setSpecies(c, "1cell");

  fail_unless( i == LIBSBML_INVALID_ATTRIBUTE_VALUE );
  fail_unless( !SpeciesReference_isSetSpecies(c) );

  SpeciesReference_free(c);
}
END_TEST


START_TEST (test_SpeciesReference_setSpecies3)
{
  SpeciesReference_t *c = 
    SpeciesReference_create(2, 2);

  int i = SpeciesReference_setSpecies(c, "mole");

  fail_unless( i == LIBSBML_OPERATION_SUCCESS );
  fail_unless( !strcmp(SpeciesReference_getSpecies(c), "mole") );
  fail_unless( SpeciesReference_isSetSpecies(c) );

  SpeciesReference_free(c);
}
END_TEST


START_TEST (test_SpeciesReference_setSpecies4)
{
  int i = SpeciesReference_setSpecies(sr, "mm");

  fail_unless( i == LIBSBML_OPERATION_SUCCESS );
  fail_unless( SpeciesReference_isSetSpecies(sr) );

  i = SpeciesReference_setSpecies(sr, NULL);

  fail_unless( i == LIBSBML_OPERATION_SUCCESS );
  fail_unless( !SpeciesReference_isSetSpecies(sr) );
}
END_TEST


START_TEST (test_SpeciesReference_setStoichiometry1)
{
  int i = SpeciesReference_setStoichiometry(sr, 2.0);

  fail_unless( i == LIBSBML_OPERATION_SUCCESS );
  fail_unless( SpeciesReference_getStoichiometry(sr) == 2.0 );
}
END_TEST


START_TEST (test_SpeciesReference_setStoichiometry2)
{
  SpeciesReference_t *c = 
    SpeciesReference_create(2, 2);

  int i = SpeciesReference_setStoichiometry(c, 4);

  fail_unless( i == LIBSBML_OPERATION_SUCCESS );
  fail_unless( SpeciesReference_getStoichiometry(c) == 4.0 );

  SpeciesReference_free(c);
}
END_TEST


START_TEST (test_SpeciesReference_setDenominator1)
{
  int i = SpeciesReference_setDenominator(sr, 2);

  fail_unless( i == LIBSBML_OPERATION_SUCCESS );
  fail_unless( SpeciesReference_getDenominator(sr) == 2 );
}
END_TEST


START_TEST (test_SpeciesReference_setDenominator2)
{
  SpeciesReference_t *c = 
    SpeciesReference_create(2, 2);

  int i = SpeciesReference_setDenominator(c, 4);

  fail_unless( i == LIBSBML_OPERATION_SUCCESS );
  fail_unless( SpeciesReference_getDenominator(c) == 4 );

  SpeciesReference_free(c);
}
END_TEST


START_TEST (test_SpeciesReference_setStoichiometryMath1)
{
  StoichiometryMath_t * sm = StoichiometryMath_create(2, 4);
  ASTNode_t *math = ASTNode_createWithType(AST_TIMES);
  ASTNode_t *a = ASTNode_create();
  ASTNode_t *b = ASTNode_create();
  ASTNode_setName(a, "a");
  ASTNode_setName(b, "b");
  ASTNode_addChild(math, a);
  ASTNode_addChild(math, b);
  StoichiometryMath_setMath(sm, math);

  int i = SpeciesReference_setStoichiometryMath(sr, sm);

  fail_unless( i == LIBSBML_OPERATION_SUCCESS);
  fail_unless( SpeciesReference_isSetStoichiometryMath(sr)   );
  fail_unless( SpeciesReference_getStoichiometry(sr) == 1 );

  i = SpeciesReference_unsetStoichiometryMath(sr);

  fail_unless( i == LIBSBML_OPERATION_SUCCESS);
  fail_unless( !SpeciesReference_isSetStoichiometryMath(sr)   );

  StoichiometryMath_free(sm);
}
END_TEST


START_TEST (test_SpeciesReference_setStoichiometryMath2)
{
  StoichiometryMath_t * sm = StoichiometryMath_create(2, 4);
  ASTNode_t *math = ASTNode_createWithType(AST_TIMES);
  ASTNode_t *a = ASTNode_create();
  ASTNode_setName(a, "a");
  ASTNode_addChild(math, a);
  StoichiometryMath_setMath(sm, math);

  int i = SpeciesReference_setStoichiometryMath(sr, sm);

  /* once the StoichiometryMath_setMath function does not set
   * an invalid ASTNode this changes to i == LIBSBML_OPERATION_SUCCESS
  fail_unless( i == LIBSBML_INVALID_OBJECT);
  */
  fail_unless ( i == LIBSBML_OPERATION_SUCCESS );
  fail_unless( SpeciesReference_isSetStoichiometryMath(sr)   );

  StoichiometryMath_free(sm);
}
END_TEST


START_TEST (test_SpeciesReference_setStoichiometryMath3)
{
  int i = SpeciesReference_setStoichiometryMath(sr, NULL);

  fail_unless( i == LIBSBML_OPERATION_SUCCESS);
  fail_unless( !SpeciesReference_isSetStoichiometryMath(sr)   );

  i = SpeciesReference_unsetStoichiometryMath(sr);

  fail_unless( i == LIBSBML_OPERATION_SUCCESS);
  fail_unless( !SpeciesReference_isSetStoichiometryMath(sr)   );

}
END_TEST


START_TEST (test_SpeciesReference_setStoichiometryMath4)
{
  StoichiometryMath_t * sm = StoichiometryMath_create(2, 4);
  ASTNode_t *math = NULL;
  StoichiometryMath_setMath(sm, math);

  int i = SpeciesReference_setStoichiometryMath(sr, sm);

  fail_unless( i == LIBSBML_OPERATION_SUCCESS);
  fail_unless( SpeciesReference_isSetStoichiometryMath(sr)   );
  fail_unless( SpeciesReference_getStoichiometry(sr) == 1 );

  i = SpeciesReference_unsetStoichiometryMath(sr);

  fail_unless( i == LIBSBML_OPERATION_SUCCESS);
  fail_unless( !SpeciesReference_isSetStoichiometryMath(sr)   );

  StoichiometryMath_free(sm);
}
END_TEST


START_TEST (test_SpeciesReference_setStoichiometryMath5)
{
  SpeciesReference_t *sr1 =
    SpeciesReference_create(1, 2);
  StoichiometryMath_t * sm = StoichiometryMath_create(2, 4);
  ASTNode_t *math = ASTNode_createWithType(AST_TIMES);
  ASTNode_t *a = ASTNode_create();
  ASTNode_t *b = ASTNode_create();
  ASTNode_setName(a, "a");
  ASTNode_setName(b, "b");
  ASTNode_addChild(math, a);
  ASTNode_addChild(math, b);
  StoichiometryMath_setMath(sm, math);

  int i = SpeciesReference_setStoichiometryMath(sr1, sm);

  fail_unless( i == LIBSBML_UNEXPECTED_ATTRIBUTE);
  fail_unless( !SpeciesReference_isSetStoichiometryMath(sr1)   );

  StoichiometryMath_free(sm);
  SpeciesReference_free(sr1);
}
END_TEST


START_TEST (test_SpeciesReference_setStoichiometryMath6)
{
  StoichiometryMath_t * sm = 
    StoichiometryMath_create(2, 1);

  int i = SpeciesReference_setStoichiometryMath(sr, sm);

  fail_unless( i == LIBSBML_VERSION_MISMATCH);
  fail_unless( !SpeciesReference_isSetStoichiometryMath(sr)   );

  StoichiometryMath_free(sm);
}
END_TEST


Suite *
create_suite_SpeciesReference_newSetters (void)
{
  Suite *suite = suite_create("SpeciesReference_newSetters");
  TCase *tcase = tcase_create("SpeciesReference_newSetters");


  tcase_add_checked_fixture( tcase,
                             SpeciesReferenceTest1_setup,
                             SpeciesReferenceTest1_teardown );

  tcase_add_test( tcase, test_SpeciesReference_setId1                );
  tcase_add_test( tcase, test_SpeciesReference_setId2                );
  tcase_add_test( tcase, test_SpeciesReference_setId3                );
  tcase_add_test( tcase, test_SpeciesReference_setId4                );
  tcase_add_test( tcase, test_SpeciesReference_setName1              );
  tcase_add_test( tcase, test_SpeciesReference_setName2              );
  tcase_add_test( tcase, test_SpeciesReference_setName3              );
  tcase_add_test( tcase, test_SpeciesReference_setName4              );
  tcase_add_test( tcase, test_SpeciesReference_setSpecies1           );
  tcase_add_test( tcase, test_SpeciesReference_setSpecies2           );
  tcase_add_test( tcase, test_SpeciesReference_setSpecies3           ); 
  tcase_add_test( tcase, test_SpeciesReference_setSpecies4           ); 
  tcase_add_test( tcase, test_SpeciesReference_setStoichiometry1     );
  tcase_add_test( tcase, test_SpeciesReference_setStoichiometry2     );
  tcase_add_test( tcase, test_SpeciesReference_setDenominator1       );
  tcase_add_test( tcase, test_SpeciesReference_setDenominator2       );
  tcase_add_test( tcase, test_SpeciesReference_setStoichiometryMath1 );
  tcase_add_test( tcase, test_SpeciesReference_setStoichiometryMath2 );
  tcase_add_test( tcase, test_SpeciesReference_setStoichiometryMath3 );
  tcase_add_test( tcase, test_SpeciesReference_setStoichiometryMath4 );
  tcase_add_test( tcase, test_SpeciesReference_setStoichiometryMath5 );
  tcase_add_test( tcase, test_SpeciesReference_setStoichiometryMath6 );


  suite_add_tcase(suite, tcase);

  return suite;
}
