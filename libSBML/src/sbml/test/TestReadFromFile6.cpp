/**
 * \file    TestReadFromFile6.cpp
 * \brief   Reads test-data/l2v2-newComponents.xml into memory and tests it.
 * \author  Sarah Keating
 *
 * $Id: TestReadFromFile6.cpp 11333 2010-06-30 06:39:53Z mhucka $
 * $HeadURL: https://sbml.svn.sourceforge.net/svnroot/sbml/trunk/libsbml/src/sbml/test/TestReadFromFile6.cpp $
 */
/* Copyright 2004 California Institute of Technology and
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

#include <sbml/SBMLReader.h>
#include <sbml/SBMLWriter.h>
#include <sbml/SBMLTypes.h>

#include <string>

#include <check.h>

LIBSBML_CPP_NAMESPACE_USE

BEGIN_C_DECLS


extern char *TestDataDirectory;


START_TEST (test_read_l2v2_newComponents)
{
  SBMLReader        reader;
  SBMLDocument*     d;
  Model*            m;
  Compartment*      c;
  CompartmentType*  ct;
  Species*          s;
  Parameter*        p;
  Reaction*         r;
  SpeciesReference* sr;
  KineticLaw*       kl;
  Constraint*       con;
  InitialAssignment* ia;
  SpeciesType*      st;
  ListOfCompartmentTypes *loct;
  CompartmentType * ct1;
  ListOfConstraints *locon;
  Constraint *con1;
  ListOfInitialAssignments *loia;
  InitialAssignment *ia1;
  ListOfReactions *lor;
  Reaction *r1;
  ListOfSpeciesReferences *losr;
  SpeciesReference * sr1;
 
  const ASTNode*   ast;

  std::string filename(TestDataDirectory);
  filename += "l2v2-newComponents.xml";


  d = reader.readSBML(filename);

  if (d == NULL)
  {
    fail("readSBML(\"l2v2-newComponents.xml\") returned a NULL pointer.");
  }



  //
  // <sbml level="2" version="2" ...>
  //
  fail_unless( d->getLevel  () == 2, NULL );
  fail_unless( d->getVersion() == 2, NULL );


  //
  //<model id="l2v2_newComponents" sboTerm="SBO:0000004">
  //
  m = d->getModel();
  fail_unless( m != NULL, NULL );

  fail_unless(m->getId() == "l2v2_newComponents", NULL);
  fail_unless(m->getSBOTerm() == 4, NULL);
  fail_unless(m->getSBOTermID() == "SBO:0000004", NULL);


  //<listOfCompartments>
  //  <compartment id="cell" compartmentType="mitochondria" size="0.013" outside="m"/>
  //  <compartment id="m" compartmentType="mitochondria" size="0.013"/>
  //</listOfCompartments>
  fail_unless( m->getNumCompartments() == 2, NULL );

  c = m->getCompartment(0);
  fail_unless( c          != NULL  , NULL );
  fail_unless( c->getId() == "cell", NULL );
  fail_unless( c->getCompartmentType() == "mitochondria", NULL );
  fail_unless( c->getOutside() == "m", NULL );
  fail_unless( c->getSBOTerm() == -1, NULL);
  fail_unless( c->getSBOTermID() == "", NULL);

  c = m->getCompartment(1);
  fail_unless( c          != NULL  , NULL );
  fail_unless( c->getId() == "m", NULL );
  fail_unless( c->getCompartmentType() == "mitochondria", NULL );

  // <listOfCompartmentTypes>
  //  <compartmentType id="mitochondria"/>
  //</listOfCompartmentTypes>
  fail_unless( m->getNumCompartmentTypes() == 1, NULL );

  ct = m->getCompartmentType(0);
  fail_unless( ct         != NULL  , NULL );
  fail_unless( ct->getId() == "mitochondria", NULL );

  loct = m->getListOfCompartmentTypes();
  ct1 = loct->get(0);
  fail_unless(ct1 == ct);

  ct1 = loct->get("mitochondria");
  fail_unless(ct1 == ct);

  //<listOfSpeciesTypes>
  //  <speciesType id="Glucose"/>
  //</listOfSpeciesTypes>
  fail_unless( m->getNumSpeciesTypes() == 1, NULL );

  st = m->getSpeciesType(0);
  fail_unless( st         != NULL  , NULL );
  fail_unless( st->getId() == "Glucose", NULL );


    //<listOfSpecies>
  //  <species id="X0" speciesType="Glucose" compartment="cell"/>
  //  <species id="X1" compartment="cell" initialConcentration="0.013"/>
  //</listOfSpecies>
  fail_unless( m->getNumSpecies() == 2, NULL );

  s = m->getSpecies(0);
  fail_unless( s          != NULL  , NULL );
  fail_unless( s->getId() == "X0", NULL );
  fail_unless( s->getSpeciesType() == "Glucose", NULL );
  fail_unless( s->getCompartment() == "cell", NULL );
  fail_unless(!s->isSetInitialAmount(), NULL);
  fail_unless(!s->isSetInitialConcentration(), NULL);

  s = m->getSpecies(1);
  fail_unless( s          != NULL  , NULL );
  fail_unless( s->getId() == "X1", NULL );
  fail_unless( !s->isSetSpeciesType(), NULL);
  fail_unless( s->getCompartment() == "cell", NULL );
  fail_unless( s->getInitialConcentration() == 0.013, NULL);
  fail_unless(!s->isSetInitialAmount(), NULL);
  fail_unless(s->isSetInitialConcentration(), NULL);





  //<listOfParameters>
  //  <parameter id="y" value="2" units="dimensionless" sboTerm="SBO:0000002"/>
  //</listOfParameters>
  fail_unless( m->getNumParameters() == 1, NULL );

  p = m->getParameter(0);
  fail_unless( p         != NULL  , NULL );
  fail_unless( p->getId() == "y", NULL );
  fail_unless( p->getValue() == 2, NULL );
  fail_unless( p->getUnits() == "dimensionless", NULL );
  fail_unless( p->getId() == "y", NULL );
  fail_unless(p->getSBOTerm() == 2, NULL);
  fail_unless(p->getSBOTermID() == "SBO:0000002", NULL);

  //<listOfConstraints>
  //  <constraint sboTerm="SBO:0000064">
  //    <math xmlns="http://www.w3.org/1998/Math/MathML">
  //      <apply>
  //        <lt/>
  //        <cn type="integer"> 1 </cn>
  //        <ci> cell </ci>
  //      </apply>
  //    </math>
  //  </constraint>
  //</listOfConstraints>
  fail_unless( m->getNumConstraints() == 1, NULL );

  con = m->getConstraint(0);
  fail_unless( con         != NULL  , NULL );
  fail_unless(con->getSBOTerm() == 64, NULL);
  fail_unless(con->getSBOTermID() == "SBO:0000064", NULL);

  ast = con->getMath();
  fail_unless(!strcmp(SBML_formulaToString(ast), "lt(1, cell)"), NULL);

  locon = m->getListOfConstraints();
  con1 = locon->get(0);
  fail_unless(con1 == con);

  //<listOfInitialAssignments>
  //  <initialAssignment symbol="X0" sboTerm="SBO:0000064">
  //    <math xmlns="http://www.w3.org/1998/Math/MathML">
  //      <apply>
  //        <times/>
  //        <ci> y </ci>
  //        <ci> X1 </ci>
  //      </apply>
  //    </math>
  //  </initialAssignment>
  //</listOfInitialAssignments>
  fail_unless( m->getNumInitialAssignments() == 1, NULL );

  ia = m->getInitialAssignment(0);
  fail_unless( ia         != NULL  , NULL );
  fail_unless(ia->getSBOTerm() == 64, NULL);
  fail_unless(ia->getSBOTermID() == "SBO:0000064", NULL);
  fail_unless(ia->getSymbol() == "X0", NULL);

  ast = ia->getMath();
  fail_unless(!strcmp(SBML_formulaToString(ast), "y * X1"), NULL);

  loia = m->getListOfInitialAssignments();
  ia1 = loia->get(0);
  fail_unless( ia1 == ia);

  ia1 = loia->get("X0");
  fail_unless( ia1 == ia);

  //<listOfReactions>
  //  <reaction id="in" sboTerm="SBO:0000231">
  //    <listOfReactants>
  //      <speciesReference id="me" name="sarah" species="X0"/>
  //    </listOfReactants>
  //    <kineticLaw sboTerm="SBO:0000001">
  //      <math xmlns="http://www.w3.org/1998/Math/MathML">
  //        <apply>
  //          <divide/>
  //          <apply>
  //            <times/>
  //            <ci> v </ci>
  //            <ci> X0 </ci>
  //          </apply>
  //          <ci> t </ci>
  //        </apply>
  //      </math>
  //      <listOfParameters>
  //        <parameter id="v" units="litre"/>
  //        <parameter id="t" units="second"/>
  //      </listOfParameters>
  //    </kineticLaw>
  //  </reaction>
  //</listOfReactions>
  fail_unless( m->getNumReactions() == 1, NULL );

  r = m->getReaction(0);
  fail_unless( r         != NULL  , NULL );
  fail_unless(r->getSBOTerm() == 231, NULL);
  fail_unless(r->getSBOTermID() == "SBO:0000231", NULL);
  fail_unless(r->getId() == "in", NULL);

  lor = m->getListOfReactions();
  r1 = lor->get(0);
  fail_unless ( r1 == r);

  r1 = lor->get("in");
  fail_unless (r1 == r);

  fail_unless(r->isSetKineticLaw(), NULL);

  kl = r->getKineticLaw();
  fail_unless( kl         != NULL  , NULL );
  fail_unless(kl->getSBOTerm() == 1, NULL);
  fail_unless(kl->getSBOTermID() == "SBO:0000001", NULL);

  fail_unless(kl->isSetMath(), NULL);

  ast = kl->getMath();
  fail_unless(!strcmp(SBML_formulaToString(ast), "v * X0 / t"), NULL);

  fail_unless(kl->getNumParameters() == 2, NULL);

  p = kl->getParameter(0);
  fail_unless( p         != NULL  , NULL );
  fail_unless(p->getSBOTerm() == 2, NULL);
  fail_unless(p->getSBOTermID() == "SBO:0000002", NULL);
  fail_unless(p->getId() == "v", NULL);
  fail_unless(p->getUnits() == "litre", NULL);

  fail_unless(r->getNumReactants() == 1, NULL);
  fail_unless(r->getNumProducts() == 0, NULL);
  fail_unless(r->getNumModifiers() == 0, NULL);

  sr = r->getReactant(0);
  fail_unless( sr         != NULL  , NULL );
  fail_unless(sr->getName() == "sarah", NULL);
  fail_unless(sr->getId() == "me", NULL);
  fail_unless(sr->getSpecies() == "X0", NULL);

  losr = r->getListOfReactants();
  sr1 = static_cast <SpeciesReference *> (losr->get(0));
  fail_unless ( sr1 == sr);

  sr1 = static_cast <SpeciesReference *> (losr->get("me"));
  fail_unless ( sr1 == sr);

  delete d;
}
END_TEST


Suite *
create_suite_TestReadFromFile6 (void)
{ 
  Suite *suite = suite_create("test-data/l2v2-newComponents.xml");
  TCase *tcase = tcase_create("test-data/l2v2-newComponents.xml");


  tcase_add_test(tcase, test_read_l2v2_newComponents);

  suite_add_tcase(suite, tcase);

  return suite;
}


END_C_DECLS
