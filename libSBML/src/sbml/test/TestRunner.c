/**
 * \file    TestRunner.c
 * \brief   Runs all unit tests in the sbml module
 * \author  Ben Bornstein
 *
 * $Id: TestRunner.c 11197 2010-04-25 21:44:53Z mhucka $
 * $HeadURL: https://sbml.svn.sourceforge.net/svnroot/sbml/trunk/libsbml/src/sbml/test/TestRunner.c $
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


#include <string.h>
#include <stdlib.h>

#include <sbml/common/extern.h>
#include <sbml/util/memory.h>

#include <check.h>


/**
 * Test suite creation function prototypes.
 *
 * These functions are needed only for calls in main() below.  Therefore a
 * separate header file is not necessary and only adds a maintenance burden
 * to keep the two files synchronized.
 */
BEGIN_C_DECLS

Suite *create_suite_ReadSBML                      (void);
Suite *create_suite_WriteSBML                     (void);
Suite *create_suite_WriteL3SBML                   (void);

Suite *create_suite_AlgebraicRule                 (void);
Suite *create_suite_AssignmentRule                (void);
Suite *create_suite_Compartment                   (void);
Suite *create_suite_Compartment_newSetters        (void);
Suite *create_suite_L3_Compartment                (void);
Suite *create_suite_CompartmentType               (void);
Suite *create_suite_CompartmentType_newSetters    (void);
Suite *create_suite_Constraint                    (void);
Suite *create_suite_Constraint_newSetters         (void);
Suite *create_suite_CompartmentVolumeRule         (void);
Suite *create_suite_Delay                         (void);
Suite *create_suite_Event                         (void);
Suite *create_suite_L3_Event                      (void);
Suite *create_suite_Event_newSetters              (void);
Suite *create_suite_EventAssignment               (void);
Suite *create_suite_EventAssignment_newSetters    (void);
Suite *create_suite_FunctionDefinition            (void);
Suite *create_suite_FunctionDefinition_newSetters (void);
Suite *create_suite_InitialAssignment             (void);
Suite *create_suite_InitialAssignment_newSetters  (void);
Suite *create_suite_KineticLaw                    (void);
Suite *create_suite_KineticLaw_newSetters         (void);
Suite *create_suite_ListOf                        (void);
Suite *create_suite_Model                         (void);
Suite *create_suite_L3_Model                      (void);
Suite *create_suite_Model_newSetters              (void);
Suite *create_suite_ModifierSpeciesReference      (void);
Suite *create_suite_Parameter                     (void);
Suite *create_suite_L3_Parameter                  (void);
Suite *create_suite_Parameter_newSetters          (void);
Suite *create_suite_ParameterRule                 (void);
Suite *create_suite_RateRule                      (void);
Suite *create_suite_Reaction                      (void);
Suite *create_suite_L3_Reaction                   (void);
Suite *create_suite_Reaction_newSetters           (void);
Suite *create_suite_Rule                          (void);
Suite *create_suite_Rule_newSetters               (void);
Suite *create_suite_RuleType                      (void);
Suite *create_suite_SBase                         (void);
Suite *create_suite_SBase_newSetters              (void);
Suite *create_suite_SBMLConvert                   (void);
Suite *create_suite_SBMLConvertStrict             (void);
Suite *create_suite_SBMLDocument                  (void);
/* these files no longer exist in libSBML - SK 13/02/2007
Suite *create_suite_SBMLFormatter                 (void);
Suite *create_suite_SBMLHandler                   (void);
*/
Suite *create_suite_SBMLError                     (void);
Suite *create_suite_SBMLReader                    (void);
Suite *create_suite_SBMLWriter                    (void);
Suite *create_suite_SimpleSpeciesReference        (void);
Suite *create_suite_Species                       (void);
Suite *create_suite_L3_Species                    (void);
Suite *create_suite_Species_newSetters            (void);
Suite *create_suite_SpeciesConcentrationRule      (void);
Suite *create_suite_SpeciesReference              (void);
Suite *create_suite_L3_SpeciesReference           (void);
Suite *create_suite_SpeciesReference_newSetters   (void);
Suite *create_suite_SpeciesType                   (void);
Suite *create_suite_SpeciesType_newSetters        (void);
Suite *create_suite_StoichiometryMath             (void);
Suite *create_suite_Trigger                       (void);
Suite *create_suite_Unit                          (void);
Suite *create_suite_L3_Unit                       (void);
Suite *create_suite_Unit_newSetters               (void);
Suite *create_suite_UnitDefinition                (void);
Suite *create_suite_UnitDefinition_newSetters     (void);
Suite *create_suite_UnitKind                      (void);

Suite *create_suite_CopyAndClone                  (void);
Suite *create_suite_TestReadFromFile1             (void);
Suite *create_suite_TestReadFromFile2             (void);
Suite *create_suite_TestReadFromFile3             (void);
Suite *create_suite_TestReadFromFile4             (void);
Suite *create_suite_TestReadFromFile5             (void);
Suite *create_suite_TestReadFromFile6             (void);
Suite *create_suite_TestReadFromFile7             (void);
Suite *create_suite_TestReadFromFile8             (void);
Suite *create_suite_TestReadFromFile9             (void);

Suite *create_suite_TestConsistencyChecks         (void);
Suite *create_suite_ParentObject                  (void);
Suite *create_suite_SBMLNamespaces                (void);
Suite *create_suite_AncestorObject                (void);
Suite *create_suite_TestInternalConsistencyChecks (void);
Suite *create_suite_HasReqdAtt                    (void);
Suite *create_suite_HasReqdElements               (void);
Suite *create_suite_SyntaxChecker                 (void);
Suite *create_suite_SBMLConstructorException      (void);

Suite *create_suite_SBMLTransforms                (void);

END_C_DECLS


/**
 * Global.
 *
 * Declared extern in TestReadFromFileN suites.
 */
char *TestDataDirectory;


/**
 * Sets TestDataDirectory for the the TestReadFromFileN suites.
 *
 * For Automake's distcheck target to work properly, TestDataDirectory must
 * begin with the value of the environment variable SRCDIR.
 */
void
setTestDataDirectory (void)
{
  char *srcdir = getenv("srcdir");
  int  length  = (srcdir == NULL) ? 0 : strlen(srcdir);


  /**
   * strlen("/test-data/") = 11 + 1 (for NULL) = 12
   */
  TestDataDirectory = (char *) safe_calloc( length + 12, sizeof(char) );

  if (srcdir != NULL)
  {
    strcpy(TestDataDirectory, srcdir);
    strcat(TestDataDirectory, "/");
  }

  strcat(TestDataDirectory, "test-data/");
}


int
main (int argc, char* argv[]) 
{ 
  int num_failed;

  setTestDataDirectory();

  SRunner *runner = srunner_create( create_suite_ReadSBML               () );
  srunner_add_suite( runner, create_suite_WriteSBML                     () );
  srunner_add_suite( runner, create_suite_WriteL3SBML                   () );
  srunner_add_suite( runner, create_suite_AlgebraicRule                 () ); 
  srunner_add_suite( runner, create_suite_AssignmentRule                () );
  srunner_add_suite( runner, create_suite_Compartment                   () );
  srunner_add_suite( runner, create_suite_L3_Compartment                () );
  srunner_add_suite( runner, create_suite_Compartment_newSetters        () );
  srunner_add_suite( runner, create_suite_CompartmentType               () );
  srunner_add_suite( runner, create_suite_CompartmentType_newSetters    () );
  srunner_add_suite( runner, create_suite_CompartmentVolumeRule         () );
  srunner_add_suite( runner, create_suite_Constraint                    () );
  srunner_add_suite( runner, create_suite_Constraint_newSetters         () );
  srunner_add_suite( runner, create_suite_Delay                         () );
  srunner_add_suite( runner, create_suite_Event                         () );
  srunner_add_suite( runner, create_suite_L3_Event                      () );
  srunner_add_suite( runner, create_suite_Event_newSetters              () );
  srunner_add_suite( runner, create_suite_EventAssignment               () );
  srunner_add_suite( runner, create_suite_EventAssignment_newSetters    () );
  srunner_add_suite( runner, create_suite_FunctionDefinition            () );
  srunner_add_suite( runner, create_suite_FunctionDefinition_newSetters () );
  srunner_add_suite( runner, create_suite_InitialAssignment             () );
  srunner_add_suite( runner, create_suite_InitialAssignment_newSetters  () );
  srunner_add_suite( runner, create_suite_KineticLaw                    () );
  srunner_add_suite( runner, create_suite_KineticLaw_newSetters         () );
  srunner_add_suite( runner, create_suite_ListOf                        () );
  srunner_add_suite( runner, create_suite_Model                         () );
  srunner_add_suite( runner, create_suite_L3_Model                      () );
  srunner_add_suite( runner, create_suite_Model_newSetters              () );
  srunner_add_suite( runner, create_suite_ModifierSpeciesReference      () );
  srunner_add_suite( runner, create_suite_Parameter                     () );
  srunner_add_suite( runner, create_suite_L3_Parameter                  () );
  srunner_add_suite( runner, create_suite_Parameter_newSetters          () );
  srunner_add_suite( runner, create_suite_ParameterRule                 () );
  srunner_add_suite( runner, create_suite_RateRule                      () );
  srunner_add_suite( runner, create_suite_Reaction                      () );
  srunner_add_suite( runner, create_suite_L3_Reaction                   () );
  srunner_add_suite( runner, create_suite_Reaction_newSetters           () );
  srunner_add_suite( runner, create_suite_Rule                          () );
  srunner_add_suite( runner, create_suite_Rule_newSetters               () );
  srunner_add_suite( runner, create_suite_SBase                         () );
  srunner_add_suite( runner, create_suite_SBase_newSetters              () );
  srunner_add_suite( runner, create_suite_Species                       () );
  srunner_add_suite( runner, create_suite_L3_Species                    () );
  srunner_add_suite( runner, create_suite_Species_newSetters            () );
  srunner_add_suite( runner, create_suite_SpeciesReference              () );
  srunner_add_suite( runner, create_suite_L3_SpeciesReference           () );
  srunner_add_suite( runner, create_suite_SpeciesReference_newSetters   () );
  srunner_add_suite( runner, create_suite_SpeciesConcentrationRule      () );
  srunner_add_suite( runner, create_suite_SpeciesType                   () );
  srunner_add_suite( runner, create_suite_SpeciesType_newSetters        () );
  srunner_add_suite( runner, create_suite_StoichiometryMath             () );
  srunner_add_suite( runner, create_suite_Trigger                       () );
  srunner_add_suite( runner, create_suite_Unit                          () );
  srunner_add_suite( runner, create_suite_L3_Unit                       () );
  srunner_add_suite( runner, create_suite_Unit_newSetters               () );
  srunner_add_suite( runner, create_suite_UnitDefinition                () );
  srunner_add_suite( runner, create_suite_UnitDefinition_newSetters     () );
  srunner_add_suite( runner, create_suite_UnitKind                      () );
  srunner_add_suite( runner, create_suite_CopyAndClone                  () );
  srunner_add_suite( runner, create_suite_SBMLConvert                   () );
  srunner_add_suite( runner, create_suite_SBMLConvertStrict             () );
  srunner_add_suite( runner, create_suite_SBMLDocument                  () );
  srunner_add_suite( runner, create_suite_SBMLError                     () );
  srunner_add_suite( runner, create_suite_TestReadFromFile1             () );
  srunner_add_suite( runner, create_suite_TestReadFromFile2             () );
  srunner_add_suite( runner, create_suite_TestReadFromFile3             () );
  srunner_add_suite( runner, create_suite_TestReadFromFile4             () );
  srunner_add_suite( runner, create_suite_TestReadFromFile5             () );
  srunner_add_suite( runner, create_suite_TestReadFromFile6             () );
  srunner_add_suite( runner, create_suite_TestReadFromFile7             () );
  srunner_add_suite( runner, create_suite_TestReadFromFile8             () );
  srunner_add_suite( runner, create_suite_TestReadFromFile9             () );
  srunner_add_suite( runner, create_suite_TestConsistencyChecks         () );
  srunner_add_suite( runner, create_suite_ParentObject                  () );
  srunner_add_suite( runner, create_suite_AncestorObject                () );
  srunner_add_suite( runner, create_suite_TestInternalConsistencyChecks () );
  srunner_add_suite( runner, create_suite_HasReqdAtt                    () );
  srunner_add_suite( runner, create_suite_HasReqdElements               () );
  srunner_add_suite( runner, create_suite_SBMLNamespaces                () );
  srunner_add_suite( runner, create_suite_SyntaxChecker                 () );
  srunner_add_suite( runner, create_suite_SBMLConstructorException      () );
  srunner_add_suite( runner, create_suite_SBMLTransforms                () );


#ifdef TRACE_MEMORY
  srunner_set_fork_status(runner, CK_NOFORK);
#else
  if (argc > 1 && !strcmp(argv[1], "-nofork"))
  {
    srunner_set_fork_status( runner, CK_NOFORK );
  }
#endif

  srunner_run_all(runner, CK_NORMAL);
  num_failed = srunner_ntests_failed(runner);

#ifdef TRACE_MEMORY

  if (MemTrace_getNumLeaks() > 0)
  {
    MemTrace_printLeaks(stdout);
  }

  MemTrace_printStatistics(stdout);

#endif

  srunner_free(runner);

  return num_failed;
}
