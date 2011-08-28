// Begin CVS Header
//   $Source: /fs/turing/cvs/copasi_dev/copasi/optimization/COptMethodSA.cpp,v $
//   $Revision: 1.20 $
//   $Name: Build-33 $
//   $Author: shoops $
//   $Date: 2010/09/02 14:30:57 $
// End CVS Header

// Copyright (C) 2010 by Pedro Mendes, Virginia Tech Intellectual
// Properties, Inc., University of Heidelberg, and The University
// of Manchester.
// All rights reserved.

// Copyright (C) 2008 by Pedro Mendes, Virginia Tech Intellectual
// Properties, Inc., EML Research, gGmbH, University of Heidelberg,
// and The University of Manchester.
// All rights reserved.

// Copyright (C) 2001 - 2007 by Pedro Mendes, Virginia Tech Intellectual
// Properties, Inc. and EML Research, gGmbH.
// All rights reserved.

// adapted by Stefan Hoops, August 2006, from the original Gepasi file
// nelmea.cpp :

#include "copasi.h"

#include "COptMethodSA.h"
#include "COptProblem.h"
#include "COptItem.h"
#include "COptTask.h"

#include "parameterFitting/CFitProblem.h"
#include "report/CCopasiObjectReference.h"
#include "randomGenerator/CRandom.h"

#define STORED 2
#define NS 10
#define K 1.0

COptMethodSA::COptMethodSA(const CCopasiContainer * pParent):
    COptMethod(CCopasiTask::optimization, CCopasiMethod::SimulatedAnnealing, pParent)
{
  addParameter("Start Temperature", CCopasiParameter::UDOUBLE, (C_FLOAT64) 1.0);
  addParameter("Cooling Factor", CCopasiParameter::UDOUBLE, (C_FLOAT64) 0.85);
  addParameter("Tolerance", CCopasiParameter::UDOUBLE, (C_FLOAT64) 1.e-006);
  addParameter("Random Number Generator", CCopasiParameter::UINT, (unsigned C_INT32) CRandom::mt19937);
  addParameter("Seed", CCopasiParameter::UINT, (unsigned C_INT32) 0);

  initObjects();
}

COptMethodSA::COptMethodSA(const COptMethodSA & src,
                           const CCopasiContainer * pParent):
    COptMethod(src, pParent)
{initObjects();}

COptMethodSA::~COptMethodSA()
{cleanup();}

void COptMethodSA::initObjects()
{
  addObjectReference("Current Temperature", mTemperature, CCopasiObject::ValueDbl);
}

bool COptMethodSA::optimise()
{
  if (!initialize())
    {
      if (mpCallBack)
        mpCallBack->finishItem(mhTemperature);

      return false;
    }

  unsigned C_INT32 i, j, k, m;

  unsigned C_INT32 h, a;
  C_FLOAT64 xc, p, c, nt, New;
  C_FLOAT64 fk[STORED];
  bool ready;

  // initial point is first guess but we have to make sure that we
  // are within the parameter domain
  for (i = 0; i < mVariableSize; i++)
    {
      const COptItem & OptItem = *(*mpOptItem)[i];

      switch (OptItem.checkConstraint())
        {
          case - 1:
            mCurrent[i] = *OptItem.getLowerBoundValue();
            break;

          case 1:
            mCurrent[i] = *OptItem.getUpperBoundValue();
            break;

          case 0:
            mCurrent[i] = OptItem.getStartValue();
            break;
        }

      (*(*mpSetCalculateVariable)[i])(mCurrent[i]);

      // The step must not contain any zeroes
      mStep[i] = std::max(fabs(mCurrent[i]), 1.0);
    }

  mCurrentValue = evaluate();

  if (!isnan(mEvaluationValue))
    {
      // and store that value
      mBestValue = mEvaluationValue;
      mContinue &= mpOptProblem->setSolution(mBestValue, mCurrent);

      // We found a new best value lets report it.
      mpParentTask->output(COutputInterface::DURING);
    }

  // store this in all positions
  for (a = 0; a < STORED; a++)
    fk[a] = mCurrentValue;

  mAccepted = 0;

  // set the number of steps at one single temperature
  nt = 5 * mVariableSize;

  if (nt < 100) nt = 100;

  // no temperature reductions yet
  k = 0;

  do // number of internal cycles: max(5*mVariableSize, 100) * NS * mVariableSize
    {
      for (m = 0; m < nt && mContinue; m++) // step adjustments
        {
          for (j = 0; j < NS && mContinue; j++) // adjustment in all directions
            {
              for (h = 0; h < mVariableSize && mContinue; h++) // adjustment in one direction
                {
                  // Calculate the step
                  xc = (2.0 * mpRandom->getRandomCC() - 1) * mStep[h];
                  New = mCurrent[h] + xc;

                  // Set the new parameter value
                  (*(*mpSetCalculateVariable)[h])(New);

                  // Check all parametric constraints
                  if (!mpOptProblem->checkParametricConstraints())
                    {
                      // Undo since not accepted
                      (*(*mpSetCalculateVariable)[h])(mCurrent[h]);
                      continue;
                    }

                  // evaluate the function
                  evaluate();

                  // Check all functional constraints
                  if (!mpOptProblem->checkFunctionalConstraints())
                    {
                      // Undo since not accepted
                      (*(*mpSetCalculateVariable)[h])(mCurrent[h]);
                      continue;
                    }

                  // here we have a completely feasible point
                  // keep if energy is reduced
                  if (mEvaluationValue <= mCurrentValue)
                    {
                      // only one value has changed...
                      mCurrent[h] = New;
                      mCurrentValue = mEvaluationValue;
                      i++;  // a new point
                      mAccepted[h]++; // a new point in this coordinate

                      if (mCurrentValue < mBestValue)
                        {
                          // and store that value
                          mBestValue = mEvaluationValue;
                          mContinue &= mpOptProblem->setSolution(mBestValue, mCurrent);

                          // We found a new best value lets report it.
                          mpParentTask->output(COutputInterface::DURING);
                        }
                    }
                  else
                    {
                      // keep with probability p, if energy is increased
                      p = exp((mCurrentValue - mEvaluationValue) / (K * mTemperature));

                      if (p > mpRandom->getRandomCO())
                        {
                          // only one value has changed...
                          mCurrent[h] = New;
                          mCurrentValue = mEvaluationValue;
                          i++;  // a new point
                          mAccepted[h]++; // a new point in this coordinate
                        }
                      else
                        // Undo since not accepted
                        (*(*mpSetCalculateVariable)[h])(mCurrent[h]);
                    }
                }
            }

          // update the step sizes
          for (a = 0; a < mVariableSize; a++)
            {
              c = (C_FLOAT64) mAccepted[a] / (C_FLOAT64) NS;

              if (c > 0.6)
                mStep[a] *= 1 + 5 * (c - 0.6);
              else if (c < 0.4)
                mStep[a] /= 1 + 5 * (0.4 - c);

              mAccepted[a] = 0;
            }
        }

      // the equilibrium energy

      k++;

      // if this is the first cycle ignore the convergence tests
      if (k == 1)
        ready = false;
      else
        {
          ready = true;

          // check termination criterion of not much change since last STORED times
          for (a = 0; a < STORED; a++)
            if (fabs(fk[a] - mCurrentValue) > mTolerance)
              {
                ready = false;
                break;
              }

          if (!ready)
            {
              for (a = 0; a < STORED - 1; a++)
                fk[a] = fk[a + 1];

              fk[STORED - 1] = mCurrentValue;
            }
          // check the termination criterion of not much larger than last optimal
          else if (fabs(mCurrentValue - mBestValue) > mTolerance)
            ready = false;
        }

      if (!ready)
        {
          i++;

          mCurrent = mpOptProblem->getSolutionVariables();

          for (a = 0; a < mVariableSize; a++)
            (*(*mpSetCalculateVariable)[a])(mCurrent[a]);

          mCurrentValue = mBestValue;
        }

      // update the temperature
      mTemperature *= mCoolingFactor;

      if (mpCallBack)
        mContinue &= mpCallBack->progressItem(mhTemperature);
    }
  while (!ready && mContinue);

  if (mpCallBack)
    mpCallBack->finishItem(mhTemperature);

  return true;
}

bool COptMethodSA::cleanup()
{
  return true;
}

const C_FLOAT64 & COptMethodSA::evaluate()
{
  // We do not need to check whether the parametric constraints are fulfilled
  // since the parameters are created within the bounds.

  mContinue &= mpOptProblem->calculate();
  mEvaluationValue = mpOptProblem->getCalculateValue();

  // When we leave the either functional domain
  // we set the objective value +Inf
  if (!mpOptProblem->checkFunctionalConstraints())
    mEvaluationValue = std::numeric_limits<C_FLOAT64>::infinity();

  return mEvaluationValue;
}

bool COptMethodSA::initialize()
{
  cleanup();

  if (!COptMethod::initialize()) return false;

  mTemperature = * getValue("Start Temperature").pUDOUBLE;
  mCoolingFactor = * getValue("Cooling Factor").pUDOUBLE;
  mTolerance = * getValue("Tolerance").pUDOUBLE;
  mpRandom =
    CRandom::createGenerator(* (CRandom::Type *) getValue("Random Number Generator").pUINT,
                             * getValue("Seed").pUINT);

  if (mpCallBack)
    mhTemperature =
      mpCallBack->addItem("Current Temperature",
                          CCopasiParameter::UDOUBLE,
                          & mTemperature,
                          NULL);

  mBestValue = std::numeric_limits<C_FLOAT64>::infinity();
  mContinue = true;

  mVariableSize = mpOptItem->size();

  mCurrent.resize(mVariableSize);
  mStep.resize(mVariableSize);
  mAccepted.resize(mVariableSize);

  return true;
}
