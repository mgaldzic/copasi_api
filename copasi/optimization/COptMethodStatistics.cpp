// Begin CVS Header
//   $Source: /fs/turing/cvs/copasi_dev/copasi/optimization/COptMethodStatistics.cpp,v $
//   $Revision: 1.5 $
//   $Name: Build-33 $
//   $Author: shoops $
//   $Date: 2009/10/13 17:39:54 $
// End CVS Header

// Copyright (C) 2008 by Pedro Mendes, Virginia Tech Intellectual
// Properties, Inc., EML Research, gGmbH, University of Heidelberg,
// and The University of Manchester.
// All rights reserved.

// Copyright (C) 2001 - 2007 by Pedro Mendes, Virginia Tech Intellectual
// Properties, Inc. and EML Research, gGmbH.
// All rights reserved.

#include "copasi.h"

#include "COptMethod.h"
#include "COptProblem.h"
#include "COptItem.h"
#include "COptTask.h"

#include "COptMethodStatistics.h"

#include "report/CCopasiObjectReference.h"
#include "randomGenerator/CRandom.h"

COptMethodStatistics::COptMethodStatistics():
    COptMethod(CCopasiTask::optimization, CCopasiMethod::Statistics)
{
  initObjects();
}

COptMethodStatistics::COptMethodStatistics(const COptMethodStatistics & src):
    COptMethod(src)
{initObjects();}

/**
 * Destructor
 */
COptMethodStatistics::~COptMethodStatistics()
{//*** added similar to coptga
  cleanup();
}

/**
 * Optimizer Function
 * Returns: true if properly initialized
 * should return a boolean
 */
bool COptMethodStatistics::initialize()
{
  cleanup();

  if (!COptMethod::initialize()) return false;

  mBestValue = std::numeric_limits< C_FLOAT64 >::infinity();

  mVariableSize = mpOptItem->size();
  mIndividual.resize(mVariableSize);

  return true;
}

/**
 * Optimizer Function
 * Returns: nothing
 * should return a boolean
 */
//C_INT32 COptMethodStatistics::optimise()
bool COptMethodStatistics::optimise()
{
  bool Continue = true;

  if (!initialize()) return false;

  unsigned C_INT32 j;

  // Initialize the population
  // first individual is the initial guess
  for (j = 0; j < mVariableSize; j++)
    {
      C_FLOAT64 & mut = mIndividual[j];
      COptItem & OptItem = *(*mpOptItem)[j];

      mut = OptItem.getStartValue();

      // force it to be within the bounds
      switch (OptItem.checkConstraint(mut))
        {
          case - 1:
            mut = *OptItem.getLowerBoundValue();
            break;

          case 1:
            mut = *OptItem.getUpperBoundValue();
            break;
        }

      // We need to set the value here so that further checks take
      // account of the value.
      (*(*mpSetCalculateVariable)[j])(mut);
    }

  Continue = evaluate(mIndividual);

  mBestValue = mValue;
  Continue = mpOptProblem->setSolution(mBestValue, mIndividual);

  // We found a new best value lets report it.
  //if (mpReport) mpReport->printBody();
  mpParentTask->output(COutputInterface::DURING);

  return Continue;
}

// evaluate the fitness of one individual
bool COptMethodStatistics::evaluate(const CVector< C_FLOAT64 > & /* individual */)
{
  bool Continue = true;

  // We do not need to check whether the parametric constraints are fulfilled
  // since the parameters are created within the bounds.

  // evaluate the fitness
  Continue = mpOptProblem->calculate();

  // check whether the functional constraints are fulfilled
  if (!mpOptProblem->checkFunctionalConstraints())
    mValue = DBL_MAX;
  else
    mValue = mpOptProblem->getCalculateValue();

  return Continue;
}
