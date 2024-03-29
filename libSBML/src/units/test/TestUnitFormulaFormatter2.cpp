/**
 * \file    TestUnitFormulaFormatter2.cpp
 * \brief   UnitFormulaFormatter unit tests
 * \author  Sarah Keating
 *
 * $Id: TestUnitFormulaFormatter2.cpp 11402 2010-07-07 01:43:53Z mhucka $
 * $HeadURL: https://sbml.svn.sourceforge.net/svnroot/sbml/trunk/libsbml/src/units/test/TestUnitFormulaFormatter2.cpp $
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
 *----------------------------------------------------------------------- -->
*/

#include <sbml/common/common.h>
#include <sbml/common/extern.h>

#include <sbml/SBMLReader.h>
#include <sbml/SBMLTypes.h>

#include <sbml/SBMLDocument.h>
#include <sbml/Model.h>

#include <sbml/units/UnitFormulaFormatter.h>

#include <check.h>

LIBSBML_CPP_NAMESPACE_USE

extern char *TestDataDirectory;

static UnitFormulaFormatter *uff;
static Model *m;
static SBMLDocument* d;

/* 
 * tests the results from different model
 * components that have units
 * e.g. compartment; species; parameter
 */

BEGIN_C_DECLS


void
UnitFormulaFormatter2Test_setup (void)
{
  d = new SBMLDocument();
 
  char *filename = safe_strcat(TestDataDirectory, "L3components.xml");


  d = readSBML(filename);
  m = d->getModel();

  uff = new UnitFormulaFormatter(m);

}


void
UnitFormulaFormatter2Test_teardown (void)
{
  delete uff;
  delete d;
}

START_TEST (test_UnitFormulaFormatter2_getUnitDefinition_compartment)
{
  UnitDefinition * ud = new UnitDefinition(3, 1);

  /* compartment with declared standard units */
  ud = uff->getUnitDefinitionFromCompartment(m->getCompartment(0));

  fail_unless(ud->getNumUnits() == 1);

  fail_unless(!strcmp(ud->getId().c_str(), ""), NULL);

  fail_unless(ud->getUnit(0)->getMultiplier() == 1);
  fail_unless(ud->getUnit(0)->getScale() == 0);
  fail_unless(ud->getUnit(0)->getExponent() == 1);
  fail_unless(ud->getUnit(0)->getOffset() == 0.0);
  fail_unless(ud->getUnit(0)->getKind() == UNIT_KIND_LITRE);

  /* compartment with declared units from unit definition */
  ud = uff->getUnitDefinitionFromCompartment(m->getCompartment(1));

  fail_unless(ud->getNumUnits() == 2);

  fail_unless(!strcmp(ud->getId().c_str(), ""), NULL);

  fail_unless(ud->getUnit(0)->getMultiplier() == 1);
  fail_unless(ud->getUnit(0)->getScale() == 0);
  fail_unless(ud->getUnit(0)->getExponent() == 1);
  fail_unless(ud->getUnit(0)->getOffset() == 0.0);
  fail_unless(ud->getUnit(0)->getKind() == UNIT_KIND_METRE);

  fail_unless(ud->getUnit(1)->getMultiplier() == 1);
  fail_unless(ud->getUnit(1)->getScale() == 0);
  fail_unless(ud->getUnit(1)->getExponent() == -1);
  fail_unless(ud->getUnit(1)->getOffset() == 0.0);
  fail_unless(ud->getUnit(1)->getKind() == UNIT_KIND_SECOND);

  /* compartment with spatial dimensions = 3 deduced from model */
  ud = uff->getUnitDefinitionFromCompartment(m->getCompartment(2));

  fail_unless(ud->getNumUnits() == 1);

  fail_unless(!strcmp(ud->getId().c_str(), ""), NULL);

  fail_unless(ud->getUnit(0)->getMultiplier() == 1);
  fail_unless(ud->getUnit(0)->getScale() == 0);
  fail_unless(ud->getUnit(0)->getExponent() == 1);
  fail_unless(ud->getUnit(0)->getOffset() == 0.0);
  fail_unless(ud->getUnit(0)->getKind() == UNIT_KIND_LITRE);


  /* compartment with spatial dimensions = 1 derived from model via ud */
  ud = uff->getUnitDefinitionFromCompartment(m->getCompartment(3));

  fail_unless(ud->getNumUnits() == 1);

  fail_unless(!strcmp(ud->getId().c_str(), ""), NULL);

  fail_unless(ud->getUnit(0)->getMultiplier() == 1);
  fail_unless(ud->getUnit(0)->getScale() == -2);
  fail_unless(ud->getUnit(0)->getExponent() == 1);
  fail_unless(ud->getUnit(0)->getOffset() == 0.0);
  fail_unless(ud->getUnit(0)->getKind() == UNIT_KIND_METRE);

  /* compartment with no declared units */
  ud = uff->getUnitDefinitionFromCompartment(m->getCompartment(4));

  fail_unless(ud->getNumUnits() == 0);


  /* check deals with invalid nodes */
  delete ud;
  UnitDefinition * ud1 = new UnitDefinition(m->getLevel(), m->getVersion());
  
  Compartment *c = new Compartment(m->getLevel(), m->getVersion());
  c->setId("c");
  c->setUnits("undefined");

  ud1 = uff->getUnitDefinitionFromCompartment(c);

  fail_unless (ud1->getNumUnits() == 0);

  delete c;
  delete ud1;
}
END_TEST

START_TEST (test_UnitFormulaFormatter2_getUnitDefinition_species)
{
  UnitDefinition * ud = new UnitDefinition(3, 1);
  
  /* species with declared standard units for substance and comp with units*/
  ud = uff->getUnitDefinitionFromSpecies(m->getSpecies(0));

  fail_unless(ud->getNumUnits() == 2);

  fail_unless(!strcmp(ud->getId().c_str(), ""), NULL);

  fail_unless(ud->getUnit(0)->getMultiplier() == 1);
  fail_unless(ud->getUnit(0)->getScale() == 0);
  fail_unless(ud->getUnit(0)->getExponent() == 1);
  fail_unless(ud->getUnit(0)->getOffset() == 0.0);
  fail_unless(ud->getUnit(0)->getKind() == UNIT_KIND_MOLE);

  fail_unless(ud->getUnit(1)->getMultiplier() == 1);
  fail_unless(ud->getUnit(1)->getScale() == 0);
  fail_unless(ud->getUnit(1)->getExponent() == -1);
  fail_unless(ud->getUnit(1)->getOffset() == 0.0);
  fail_unless(ud->getUnit(1)->getKind() == UNIT_KIND_LITRE);

  ud = uff->getSpeciesSubstanceUnitDefinition(m->getSpecies(0));

  fail_unless(ud->getNumUnits() == 1);

  fail_unless(!strcmp(ud->getId().c_str(), ""), NULL);

  fail_unless(ud->getUnit(0)->getMultiplier() == 1);
  fail_unless(ud->getUnit(0)->getScale() == 0);
  fail_unless(ud->getUnit(0)->getExponent() == 1);
  fail_unless(ud->getUnit(0)->getOffset() == 0.0);
  fail_unless(ud->getUnit(0)->getKind() == UNIT_KIND_MOLE);

  /* species with declared standard units for substance 
  and hasOnlySubstanceUnits = 1*/
  ud = uff->getUnitDefinitionFromSpecies(m->getSpecies(1));

  fail_unless(ud->getNumUnits() == 1);

  fail_unless(!strcmp(ud->getId().c_str(), ""), NULL);

  fail_unless(ud->getUnit(0)->getMultiplier() == 1);
  fail_unless(ud->getUnit(0)->getScale() == -2);
  fail_unless(ud->getUnit(0)->getExponent() == 1);
  fail_unless(ud->getUnit(0)->getOffset() == 0.0);
  fail_unless(ud->getUnit(0)->getKind() == UNIT_KIND_MOLE);

  ud = uff->getSpeciesSubstanceUnitDefinition(m->getSpecies(1));

  fail_unless(ud->getNumUnits() == 1);

  fail_unless(!strcmp(ud->getId().c_str(), ""), NULL);

  fail_unless(ud->getUnit(0)->getMultiplier() == 1);
  fail_unless(ud->getUnit(0)->getScale() == -2);
  fail_unless(ud->getUnit(0)->getExponent() == 1);
  fail_unless(ud->getUnit(0)->getOffset() == 0.0);
  fail_unless(ud->getUnit(0)->getKind() == UNIT_KIND_MOLE);

  /* species with units from model substance*/
  ud = uff->getUnitDefinitionFromSpecies(m->getSpecies(2));

  fail_unless(ud->getNumUnits() == 2);

  fail_unless(!strcmp(ud->getId().c_str(), ""), NULL);

  fail_unless(ud->getUnit(0)->getMultiplier() == 1);
  fail_unless(ud->getUnit(0)->getScale() == 0);
  fail_unless(ud->getUnit(0)->getExponent() == 1);
  fail_unless(ud->getUnit(0)->getOffset() == 0.0);
  fail_unless(ud->getUnit(0)->getKind() == UNIT_KIND_MOLE);

  fail_unless(ud->getUnit(1)->getMultiplier() == 1);
  fail_unless(ud->getUnit(1)->getScale() == 0);
  fail_unless(ud->getUnit(1)->getExponent() == -1);
  fail_unless(ud->getUnit(1)->getOffset() == 0.0);
  fail_unless(ud->getUnit(1)->getKind() == UNIT_KIND_LITRE);

  ud = uff->getSpeciesSubstanceUnitDefinition(m->getSpecies(2));
 
  fail_unless(ud->getNumUnits() == 1);

  fail_unless(!strcmp(ud->getId().c_str(), ""), NULL);

  fail_unless(ud->getUnit(0)->getMultiplier() == 1);
  fail_unless(ud->getUnit(0)->getScale() == 0);
  fail_unless(ud->getUnit(0)->getExponent() == 1);
  fail_unless(ud->getUnit(0)->getOffset() == 0.0);
  fail_unless(ud->getUnit(0)->getKind() == UNIT_KIND_MOLE);

  /* species with declared standard units for substance and no spatialSizeUnits*/
  ud = uff->getUnitDefinitionFromSpecies(m->getSpecies(3));

  fail_unless(ud->getNumUnits() == 0);

  ud = uff->getSpeciesSubstanceUnitDefinition(m->getSpecies(3));
 
  fail_unless(ud->getNumUnits() == 1);
 
  fail_unless(!strcmp(ud->getId().c_str(), ""), NULL);

  fail_unless(ud->getUnit(0)->getMultiplier() == 1);
  fail_unless(ud->getUnit(0)->getScale() == 0);
  fail_unless(ud->getUnit(0)->getExponent() == 1);
  fail_unless(ud->getUnit(0)->getOffset() == 0.0);
  fail_unless(ud->getUnit(0)->getKind() == UNIT_KIND_MOLE);


  ///* check deals with invalid nodes */
  delete ud;
  UnitDefinition * ud1 = new UnitDefinition(m->getLevel(), m->getVersion());
  
  Species *s = new Species(m->getLevel(), m->getVersion());
  s->setId("s");
  s->setUnits("undefined");

  ud1 = uff->getUnitDefinitionFromSpecies(s);

  fail_unless (ud1->getNumUnits() == 0);

  ud1 = uff->getSpeciesSubstanceUnitDefinition(s);

  fail_unless (ud1->getNumUnits() == 0);

  s->setUnits("mole"); // here the compartment size will be NULL

  ud1 = uff->getUnitDefinitionFromSpecies(s);

  fail_unless (ud1->getNumUnits() == 1);

  ud1 = uff->getSpeciesSubstanceUnitDefinition(s);

  fail_unless (ud1->getNumUnits() == 1);

  delete s;
  delete ud1;

}
END_TEST

START_TEST (test_UnitFormulaFormatter2_getUnitDefinition_parameter)
{
  UnitDefinition * ud = new UnitDefinition(3, 1);
 
  /* parameter with declared standard units */
  ud = uff->getUnitDefinitionFromParameter(m->getParameter(0));

  fail_unless(ud->getNumUnits() == 1);

  fail_unless(!strcmp(ud->getId().c_str(), ""), NULL);

  fail_unless(ud->getUnit(0)->getMultiplier() == 1);
  fail_unless(ud->getUnit(0)->getScale() == 0);
  fail_unless(ud->getUnit(0)->getExponent() == 1);
  fail_unless(ud->getUnit(0)->getOffset() == 0.0);
  fail_unless(ud->getUnit(0)->getKind() == UNIT_KIND_METRE);

  /* parameter with declared units from unit definition */
  ud = uff->getUnitDefinitionFromParameter(m->getParameter(1));

  fail_unless(ud->getNumUnits() == 2);

  fail_unless(!strcmp(ud->getId().c_str(), ""), NULL);

  fail_unless(ud->getUnit(0)->getMultiplier() == 1);
  fail_unless(ud->getUnit(0)->getScale() == 0);
  fail_unless(ud->getUnit(0)->getExponent() == 1);
  fail_unless(ud->getUnit(0)->getOffset() == 0.0);
  fail_unless(ud->getUnit(0)->getKind() == UNIT_KIND_METRE);

  fail_unless(ud->getUnit(1)->getMultiplier() == 1);
  fail_unless(ud->getUnit(1)->getScale() == 0);
  fail_unless(ud->getUnit(1)->getExponent() == -1);
  fail_unless(ud->getUnit(1)->getOffset() == 0.0);
  fail_unless(ud->getUnit(1)->getKind() == UNIT_KIND_SECOND);

  /* parameter with no declared units */
  ud = uff->getUnitDefinitionFromParameter(m->getParameter(2));

  fail_unless(ud->getNumUnits() == 0);

  fail_unless(!strcmp(ud->getId().c_str(), ""), NULL);
 
  /* check deals with invalid nodes */
  delete ud;
  UnitDefinition * ud1 = new UnitDefinition(m->getLevel(), m->getVersion());
  
  Parameter *p = new Parameter(m->getLevel(), m->getVersion());
  p->setId("p");
  p->setUnits("undefined");

  ud1 = uff->getUnitDefinitionFromParameter(p);

  fail_unless (ud1->getNumUnits() == 0);

  delete p;
  delete ud1;

}
END_TEST

START_TEST (test_UnitFormulaFormatter2_getUnitDefinition_function)
{
  UnitDefinition * ud;
 
  /* function applied to numbers only */
  ud = uff->getUnitDefinition(m->getRule(0)->getMath());

  fail_unless(ud->getNumUnits() == 0);

  fail_unless(!strcmp(ud->getId().c_str(), ""), NULL);


  /* function applied to components */
  uff->resetFlags();
  ud = uff->getUnitDefinition(m->getRule(1)->getMath());

  fail_unless(ud->getNumUnits() == 1);

  fail_unless(!strcmp(ud->getId().c_str(), ""), NULL);

  fail_unless(ud->getUnit(0)->getMultiplier() == 1);
  fail_unless(ud->getUnit(0)->getScale() == 0);
  fail_unless(ud->getUnit(0)->getExponent() == 1);
  fail_unless(ud->getUnit(0)->getOffset() == 0.0);
  fail_unless(ud->getUnit(0)->getKind() == UNIT_KIND_METRE);

  /* function with two arguments but only one bvar */
  uff->resetFlags();
  ud = uff->getUnitDefinition(m->getRule(2)->getMath());

  fail_unless(ud->getNumUnits() == 1);

  fail_unless(!strcmp(ud->getId().c_str(), ""), NULL);

  fail_unless(ud->getUnit(0)->getMultiplier() == 1);
  fail_unless(ud->getUnit(0)->getScale() == 0);
  fail_unless(ud->getUnit(0)->getExponent() == 3);
  fail_unless(ud->getUnit(0)->getOffset() == 0.0);
  fail_unless(ud->getUnit(0)->getKind() == UNIT_KIND_METRE);

  ///* function with two arguments but only one bvar */
  uff->resetFlags();
  ud = uff->getUnitDefinition(m->getRule(3)->getMath());

  fail_unless(ud->getNumUnits() == 1);

  fail_unless(!strcmp(ud->getId().c_str(), ""), NULL);

  fail_unless(ud->getUnit(0)->getMultiplier() == 1);
  fail_unless(ud->getUnit(0)->getScale() == 0);
  fail_unless(ud->getUnit(0)->getExponent() == -1);
  fail_unless(ud->getUnit(0)->getOffset() == 0.0);
  fail_unless(ud->getUnit(0)->getKind() == UNIT_KIND_SECOND);

  delete ud;

}
END_TEST


START_TEST (test_UnitFormulaFormatter2_getUnitDefinition_event)
{
  UnitDefinition * ud = new UnitDefinition(3, 1);

  /* event with no time units */
  ud = uff->getUnitDefinitionFromEventTime(m->getEvent(0));

  fail_unless(ud->getNumUnits() == 1);

  fail_unless(!strcmp(ud->getId().c_str(), ""), NULL);

  fail_unless(ud->getUnit(0)->getMultiplier() == 60);
  fail_unless(ud->getUnit(0)->getScale() == 0);
  fail_unless(ud->getUnit(0)->getExponent() == 1);
  fail_unless(ud->getUnit(0)->getOffset() == 0.0);
  fail_unless(ud->getUnit(0)->getKind() == UNIT_KIND_SECOND);


  /* check deals with invalid nodes */
  delete ud;
  UnitDefinition * ud1 = new UnitDefinition(m->getLevel(), m->getVersion());
  
  Event *e = new Event(m->getLevel(), m->getVersion());
  e->setId("p");
  m->unsetTimeUnits();

  ud1 = uff->getUnitDefinitionFromEventTime(e);

  fail_unless (ud1->getNumUnits() == 0);

  delete e;
  delete ud1;

}
END_TEST


START_TEST (test_UnitFormulaFormatter2_getUnitDefinition_model_extent)
{
  UnitDefinition * ud = new UnitDefinition(3, 1);

  /* model extent units */
  ud = uff->getExtentUnitDefinition();

  fail_unless(ud->getNumUnits() == 1);

  fail_unless(!strcmp(ud->getId().c_str(), ""), NULL);

  fail_unless(ud->getUnit(0)->getMultiplier() == 1);
  fail_unless(ud->getUnit(0)->getScale() == 0);
  fail_unless(ud->getUnit(0)->getExponent() == 1);
  fail_unless(ud->getUnit(0)->getOffset() == 0.0);
  fail_unless(ud->getUnit(0)->getKind() == UNIT_KIND_MOLE);


  /* check deals with invalid nodes */
  delete ud;
  UnitDefinition * ud1 = new UnitDefinition(m->getLevel(), m->getVersion());
  
  m->unsetExtentUnits();

  ud1 = uff->getExtentUnitDefinition();

  fail_unless (ud1->getNumUnits() == 0);

  delete ud1;

}
END_TEST


START_TEST (test_UnitFormulaFormatter2_getUnitDefinition_species_extent)
{
  UnitDefinition * ud = new UnitDefinition(3, 1);

  /* model extent units + conversion on species*/
  ud = uff->getSpeciesExtentUnitDefinition(m->getSpecies(0));

  fail_unless(ud->getNumUnits() == 2);
  
  fail_unless(!strcmp(ud->getId().c_str(), ""), NULL);

  fail_unless(ud->getUnit(0)->getMultiplier() == 1);
  fail_unless(ud->getUnit(0)->getScale() == 0);
  fail_unless(ud->getUnit(0)->getExponent() == 1);
  fail_unless(ud->getUnit(0)->getOffset() == 0.0);
  fail_unless(ud->getUnit(0)->getKind() == UNIT_KIND_MOLE);

  fail_unless(ud->getUnit(1)->getMultiplier() == 1);
  fail_unless(ud->getUnit(1)->getScale() == 0);
  fail_unless(ud->getUnit(1)->getExponent() == 1);
  fail_unless(ud->getUnit(1)->getOffset() == 0.0);
  fail_unless(ud->getUnit(1)->getKind() == UNIT_KIND_SECOND);


  /* model extent units + conversion on model*/
  ud = uff->getSpeciesExtentUnitDefinition(m->getSpecies(1));

  fail_unless(ud->getNumUnits() == 2);
  
  fail_unless(!strcmp(ud->getId().c_str(), ""), NULL);

  fail_unless(ud->getUnit(0)->getMultiplier() == 1);
  fail_unless(ud->getUnit(0)->getScale() == 0);
  fail_unless(ud->getUnit(0)->getExponent() == 1);
  fail_unless(ud->getUnit(0)->getOffset() == 0.0);
  fail_unless(ud->getUnit(0)->getKind() == UNIT_KIND_MOLE);

  fail_unless(ud->getUnit(1)->getMultiplier() == 1);
  fail_unless(ud->getUnit(1)->getScale() == 0);
  fail_unless(ud->getUnit(1)->getExponent() == 1);
  fail_unless(ud->getUnit(1)->getOffset() == 0.0);
  fail_unless(ud->getUnit(1)->getKind() == UNIT_KIND_SECOND);

  /* check deals with invalid nodes */
  delete ud;
  UnitDefinition * ud1 = new UnitDefinition(m->getLevel(), m->getVersion());
  
  m->unsetConversionFactor();

  ud1 = uff->getSpeciesExtentUnitDefinition(m->getSpecies(1));

  fail_unless (ud1->getNumUnits() == 0);
  
  m->getSpecies(0)->unsetConversionFactor();

  ud1 = uff->getSpeciesExtentUnitDefinition(m->getSpecies(0));

  fail_unless (ud1->getNumUnits() == 0);

  m->getSpecies(0)->setConversionFactor("cf");

  m->unsetExtentUnits();

  ud1 = uff->getSpeciesExtentUnitDefinition(m->getSpecies(0));

  fail_unless (ud1->getNumUnits() == 0);

  delete ud1;

}
END_TEST


Suite *
create_suite_UnitFormulaFormatter2 (void)
{
  Suite *suite = suite_create("UnitFormulaFormatter2");
  TCase *tcase = tcase_create("UnitFormulaFormatter2");

  tcase_add_checked_fixture(tcase,
                            UnitFormulaFormatter2Test_setup,
                            UnitFormulaFormatter2Test_teardown);

  tcase_add_test(tcase, test_UnitFormulaFormatter2_getUnitDefinition_compartment );
  tcase_add_test(tcase, test_UnitFormulaFormatter2_getUnitDefinition_species );
  tcase_add_test(tcase, test_UnitFormulaFormatter2_getUnitDefinition_parameter );
  tcase_add_test(tcase, test_UnitFormulaFormatter2_getUnitDefinition_function );
  tcase_add_test(tcase, test_UnitFormulaFormatter2_getUnitDefinition_event);
  tcase_add_test(tcase, test_UnitFormulaFormatter2_getUnitDefinition_model_extent);
  tcase_add_test(tcase, test_UnitFormulaFormatter2_getUnitDefinition_species_extent);

  suite_add_tcase(suite, tcase);

  return suite;
}


END_C_DECLS
