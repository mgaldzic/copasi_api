/**
 * \file    TestRequiredElements.cpp
 * \brief   Test hasRequiredElements unit tests
 * \author  Sarah Keating
 *
 * $Id: TestRequiredElements.cpp 11635 2010-08-03 04:05:22Z mhucka $
 * $HeadURL: https://sbml.svn.sourceforge.net/svnroot/sbml/trunk/libsbml/src/sbml/test/TestRequiredElements.cpp $
 */
/* Copyright 2002 California Institute of Technology and Japan Science and
 * Technology Corporation.
 *
 * This library is free software; you can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation.  A copy of the license agreement is
 * provided in the file named "LICENSE.txt" included with this software
 * distribution.  It is also available online at
 * http://sbml.org/software/libsbml/license.html
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this library; if not, write to the Free Software Foundation,
 * Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA.
 */


#include <sbml/common/common.h>

#include <sbml/SBase.h>
#include <sbml/Compartment.h>
#include <sbml/CompartmentType.h>
#include <sbml/Constraint.h>
#include <sbml/Delay.h>
#include <sbml/Event.h>
#include <sbml/EventAssignment.h>
#include <sbml/FunctionDefinition.h>
#include <sbml/InitialAssignment.h>
#include <sbml/KineticLaw.h>
#include <sbml/ListOf.h>
#include <sbml/Model.h>
#include <sbml/Parameter.h>
#include <sbml/Reaction.h>
#include <sbml/SBMLDocument.h>
#include <sbml/Species.h>
#include <sbml/SpeciesReference.h>
#include <sbml/SpeciesType.h>
#include <sbml/Unit.h>
#include <sbml/UnitDefinition.h>
#include <sbml/units/FormulaUnitsData.h>

#include <sbml/math/ASTNode.h>

#include <check.h>

/** @cond doxygen-ignored */

using namespace std;
LIBSBML_CPP_NAMESPACE_USE

/** @endcond */


CK_CPPSTART

START_TEST ( test_Compartment )
{
  Compartment* c = new Compartment(2, 4);
  
  fail_unless (c->hasRequiredElements());

  delete c;
}
END_TEST

START_TEST ( test_CompartmentType )
{
  CompartmentType* ct = new CompartmentType(2, 4);
  
  fail_unless (ct->hasRequiredElements());

  delete ct;
}
END_TEST

START_TEST ( test_Constraint )
{
  Constraint* c = new Constraint(2, 4);

  fail_unless (!(c->hasRequiredElements()));

  c->setMath(SBML_parseFormula("a+b"));
  
  fail_unless (c->hasRequiredElements());

  delete c;
}
END_TEST

START_TEST ( test_Delay )
{
  Delay* d = new Delay(2, 4);
  
  fail_unless (!(d->hasRequiredElements()));

  d->setMath(SBML_parseFormula("a+b"));

  fail_unless (d->hasRequiredElements());

  delete d;
}
END_TEST

START_TEST ( test_Event )
{
  Event* e = new Event(2, 4);
  
  fail_unless (!(e->hasRequiredElements()));

  Trigger *t = new Trigger(2, 4);

  e->setTrigger(t);

  fail_unless (!(e->hasRequiredElements()));

  e->createEventAssignment();

  fail_unless (e->hasRequiredElements());

  delete e;
}
END_TEST

START_TEST ( test_EventAssignment )
{
  EventAssignment* ea = new EventAssignment(2, 4);
  
  fail_unless (!(ea->hasRequiredElements()));

  ea->setMath(SBML_parseFormula("fd"));

  fail_unless (ea->hasRequiredElements());

  delete ea;
}
END_TEST

START_TEST ( test_FunctionDefinition )
{
  FunctionDefinition* fd = new FunctionDefinition(2, 4);
  
  fail_unless (!(fd->hasRequiredElements()));

  fd->setMath(SBML_parseFormula("fd"));

  fail_unless (fd->hasRequiredElements());

  delete fd;
}
END_TEST

START_TEST ( test_InitialAssignment )
{
  InitialAssignment* ia = new InitialAssignment(2, 4);
  
  fail_unless (!(ia->hasRequiredElements()));

  ia->setMath(SBML_parseFormula("ia"));

  fail_unless (ia->hasRequiredElements());

  delete ia;
}
END_TEST

START_TEST ( test_KineticLaw )
{
  KineticLaw* kl = new KineticLaw(2, 4);
   
  fail_unless (!(kl->hasRequiredElements()));

  kl->setMath(SBML_parseFormula("kl"));

  fail_unless (kl->hasRequiredElements());

  delete kl;
}
END_TEST

START_TEST ( test_Model )
{
  Model* m = new Model(2, 4);
  
  fail_unless (m->hasRequiredElements());

  delete m;
}
END_TEST

START_TEST ( test_Model_L1V2 )
{
  Model* m = new Model(1, 2);
  
  fail_unless (!(m->hasRequiredElements()));

  m->createCompartment();

  fail_unless (m->hasRequiredElements());

  delete m;
}
END_TEST

START_TEST ( test_Model_L1V1 )
{
  Model* m = new Model(1, 1);
  
  fail_unless (!(m->hasRequiredElements()));

  m->createCompartment();

  fail_unless (!(m->hasRequiredElements()));

  m->createSpecies();
  
  fail_unless (!(m->hasRequiredElements()));

  m->createReaction();

  fail_unless (m->hasRequiredElements());

  delete m;
}
END_TEST

START_TEST ( test_Parameter )
{
  Parameter* p = new Parameter(2, 4);
  
  fail_unless (p->hasRequiredElements());

  delete p;
}
END_TEST

START_TEST ( test_Reaction )
{
  Reaction* r = new Reaction(2, 4);
  
  fail_unless (r->hasRequiredElements());

  delete r;
}
END_TEST

START_TEST ( test_AlgebraicRule )
{
  AlgebraicRule* ar = new AlgebraicRule(2, 4);
  
  fail_unless (!(ar->hasRequiredElements()));

  ar->setMath(SBML_parseFormula("ar"));

  fail_unless (ar->hasRequiredElements());

  delete ar;
}
END_TEST

START_TEST ( test_AssignmentRule )
{
  AssignmentRule* r = new AssignmentRule(2, 4);
  
  fail_unless (!(r->hasRequiredElements()));

  r->setMath(SBML_parseFormula("ar"));

  fail_unless (r->hasRequiredElements());

  delete r;
}
END_TEST

START_TEST ( test_RateRule )
{
  RateRule* r = new RateRule(2, 4);
  
  fail_unless (!(r->hasRequiredElements()));

  r->setMath(SBML_parseFormula("ar"));

  fail_unless (r->hasRequiredElements());

  delete r;
}
END_TEST

START_TEST ( test_Species )
{
  Species* s = new Species(2, 4);
  
  fail_unless (s->hasRequiredElements());

  delete s;
}
END_TEST

START_TEST ( test_SpeciesReference )
{
  SpeciesReference* sr = new SpeciesReference(2, 4);
  
  fail_unless (sr->hasRequiredElements());

  delete sr;
}
END_TEST

START_TEST ( test_ModifierSpeciesReference )
{
  ModifierSpeciesReference* msr = new ModifierSpeciesReference(2, 4);
  
  fail_unless (msr->hasRequiredElements());

  delete msr;
}
END_TEST

START_TEST ( test_SpeciesType )
{
  SpeciesType* st = new SpeciesType(2, 4);
  
  fail_unless (st->hasRequiredElements());

  delete st;
}
END_TEST

START_TEST ( test_StoichiometryMath )
{
  StoichiometryMath* sm = new StoichiometryMath(2, 4);

  fail_unless (!(sm->hasRequiredElements()));

  sm->setMath(SBML_parseFormula("ar"));
  
  fail_unless (sm->hasRequiredElements());

  delete sm;
}
END_TEST

START_TEST ( test_Trigger )
{
  Trigger* t = new Trigger(2, 4);
  
  fail_unless (!(t->hasRequiredElements()));

  t->setMath(SBML_parseFormula("ar"));

  fail_unless (t->hasRequiredElements());

  delete t;
}
END_TEST

START_TEST ( test_Unit )
{
  Unit* u = new Unit(2, 4);
  
  fail_unless (u->hasRequiredElements());

  delete u;
}
END_TEST

START_TEST ( test_UnitDefinition )
{
  UnitDefinition* ud = new UnitDefinition(2, 4);
  
  fail_unless (!(ud->hasRequiredElements()));

  ud->createUnit();

  fail_unless (ud->hasRequiredElements());

  delete ud;
}
END_TEST

START_TEST ( test_UnitDefinition_L1 )
{
  UnitDefinition* ud = new UnitDefinition(1, 2);
  
  fail_unless (ud->hasRequiredElements());

  delete ud;
}
END_TEST

Suite *
create_suite_HasReqdElements (void)
{
  Suite *suite = suite_create("HasReqdElements");
  TCase *tcase = tcase_create("HasReqdElements");

  tcase_add_test( tcase, test_Compartment);
  tcase_add_test( tcase, test_CompartmentType);
  tcase_add_test( tcase, test_Constraint);
  tcase_add_test( tcase, test_Delay);
  tcase_add_test( tcase, test_Event);
  tcase_add_test( tcase, test_EventAssignment);
  tcase_add_test( tcase, test_FunctionDefinition);
  tcase_add_test( tcase, test_InitialAssignment);
  tcase_add_test( tcase, test_KineticLaw);
  tcase_add_test( tcase, test_Model);
  tcase_add_test( tcase, test_Model_L1V2);
  tcase_add_test( tcase, test_Model_L1V1);
  tcase_add_test( tcase, test_Parameter);
  tcase_add_test( tcase, test_Reaction);
  tcase_add_test( tcase, test_AlgebraicRule);
  tcase_add_test( tcase, test_AssignmentRule);
  tcase_add_test( tcase, test_RateRule);
  tcase_add_test( tcase, test_Species);
  tcase_add_test( tcase, test_SpeciesReference);
  tcase_add_test( tcase, test_ModifierSpeciesReference);
  tcase_add_test( tcase, test_StoichiometryMath);
  tcase_add_test( tcase, test_SpeciesType);
  tcase_add_test( tcase, test_Trigger);
  tcase_add_test( tcase, test_Unit);
  tcase_add_test( tcase, test_UnitDefinition);
  tcase_add_test( tcase, test_UnitDefinition_L1);


  suite_add_tcase(suite, tcase);

  return suite;
}
CK_CPPEND
