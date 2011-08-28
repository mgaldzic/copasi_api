// Begin CVS Header
//   $Source: /fs/turing/cvs/copasi_dev/copasi/optimization/COptMethodGA.cpp,v $
//   $Revision: 1.56 $
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

#include <float.h>

#include <limits.h>
#include <string.h>

#include "copasi.h"
#include "mathematics.h"
#include "COptMethodGA.h"
#include "COptProblem.h"
#include "COptItem.h"
#include "COptTask.h"

#include "randomGenerator/CRandom.h"
#include "utilities/CProcessReport.h"
#include "utilities/CSort.h"
#include "report/CCopasiObjectReference.h"

COptMethodGA::COptMethodGA(const CCopasiContainer * pParent):
    COptMethod(CCopasiTask::optimization, CCopasiMethod::GeneticAlgorithm, pParent),
    mGenerations(0),
    mPopulationSize(0),
    mpRandom(NULL),
    mVariableSize(0),
    mIndividual(0),
    mCrossOverFalse(0),
    mCrossOver(0),
    mEvaluationValue(DBL_MAX),
    mValue(0),
    mShuffle(0),
    mLosses(0),
    mPivot(0),
    mMutationVarians(0.1),
    mBestValue(DBL_MAX),
    mBestIndex(C_INVALID_INDEX),
    mGeneration(0)

{
  addParameter("Number of Generations", CCopasiParameter::UINT, (unsigned C_INT32) 200);
  addParameter("Population Size", CCopasiParameter::UINT, (unsigned C_INT32) 20);
  addParameter("Random Number Generator", CCopasiParameter::UINT, (unsigned C_INT32) CRandom::mt19937);
  addParameter("Seed", CCopasiParameter::UINT, (unsigned C_INT32) 0);

  initObjects();
}

COptMethodGA::COptMethodGA(const COptMethodGA & src,
                           const CCopasiContainer * pParent):
    COptMethod(src, pParent),
    mGenerations(0),
    mPopulationSize(0),
    mpRandom(NULL),
    mVariableSize(0),
    mIndividual(0),
    mCrossOverFalse(0),
    mCrossOver(0),
    mEvaluationValue(DBL_MAX),
    mValue(0),
    mShuffle(0),
    mLosses(0),
    mPivot(0),
    mMutationVarians(0.1),
    mBestValue(DBL_MAX),
    mBestIndex(C_INVALID_INDEX),
    mGeneration(0)
{initObjects();}

COptMethodGA::~COptMethodGA()
{cleanup();}

// evaluate the fitness of one individual
bool COptMethodGA::evaluate(const CVector< C_FLOAT64 > & /* individual */)
{
  bool Continue = true;

  // We do not need to check whether the parametric constraints are fulfilled
  // since the parameters are created within the bounds.

  // evaluate the fitness
  Continue &= mpOptProblem->calculate();

  // check whether the functional constraints are fulfilled
  if (!mpOptProblem->checkFunctionalConstraints())
    mEvaluationValue = std::numeric_limits<C_FLOAT64>::infinity();
  else
    mEvaluationValue = mpOptProblem->getCalculateValue();

  return Continue;
}

bool COptMethodGA::swap(unsigned C_INT32 from, unsigned C_INT32 to)
{
  CVector< C_FLOAT64 > * pTmp = mIndividual[to];
  mIndividual[to] = mIndividual[from];
  mIndividual[from] = pTmp;

  C_FLOAT64 dTmp = mValue[to];
  mValue[to] = mValue[from];
  mValue[from] = dTmp;

  C_INT32 iTmp = mLosses[to];
  mLosses[to] = mLosses[from];
  mLosses[from] = iTmp;

  return true;
}

//mutate one individual
bool COptMethodGA::mutate(CVector< C_FLOAT64 > & individual)
{
  unsigned C_INT32 j;

  // mutate the parameters
  for (j = 0; j < mVariableSize; j++)
    {
      COptItem & OptItem = *(*mpOptItem)[j];
      C_FLOAT64 & mut = individual[j];

      // calculate the mutatated parameter
      mut *= mpRandom->getRandomNormal(1, mMutationVarians);

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

  return true;
}

bool COptMethodGA::crossover(const CVector< C_FLOAT64 > & parent1,
                             const CVector< C_FLOAT64 > & parent2,
                             CVector< C_FLOAT64 > & child1,
                             CVector< C_FLOAT64 > & child2)
{
  unsigned C_INT32 i, crp;
  unsigned C_INT32 nCross = 0;

  mCrossOver = mCrossOverFalse;

  if (mVariableSize > 1)
    nCross = mpRandom->getRandomU(mVariableSize / 2);

  if (nCross == 0)
    {
      // if less than 0 just copy parent to child
      child1 = parent1;
      child2 = parent2;

      return true;
    }

  // choose cross over points;
  // We do not mind if a crossover point gets drawn twice
  for (i = 0; i < nCross; i++)
    {
      crp = mpRandom->getRandomU(mVariableSize - 1);
      mCrossOver[crp] = true;
    }

  const CVector< C_FLOAT64 > * pParent1 = & parent1;

  const CVector< C_FLOAT64 > * pParent2 = & parent2;

  const CVector< C_FLOAT64 > * pTmp;

  for (i = 0; i < mVariableSize; i++)
    {
      if (mCrossOver[i])
        {
          pTmp = pParent1;
          pParent1 = pParent2;
          pParent2 = pTmp;
        }

      child1[i] = (*pParent1)[i];
      child2[i] = (*pParent2)[i];
    }

  return true;
}

bool COptMethodGA::shuffle()
{
  unsigned C_INT32 i, from, to, tmp;

  //  Why sort first? We can just keep shuffling.
  //  for(i=0; i<mPopulationSize; i++) mShuffle[i] = i;

  for (i = 0; i < mPopulationSize / 2; i++)
    {
      from = mpRandom->getRandomU(mPopulationSize - 1);
      to = mpRandom->getRandomU(mPopulationSize - 1);

      tmp = mShuffle[to];
      mShuffle[to] = mShuffle[from];
      mShuffle[from] = tmp;
    }

  return true;
}

bool COptMethodGA::replicate()
{
  unsigned C_INT32 i;
  bool Continue = true;

  // generate a random order for the parents
  shuffle();

  // reproduce in consecutive pairs
  for (i = 0; i < mPopulationSize / 2; i++)
    crossover(*mIndividual[mShuffle[i * 2]],
              *mIndividual[mShuffle[i * 2 + 1]],
              *mIndividual[mPopulationSize + i * 2],
              *mIndividual[mPopulationSize + i * 2 + 1]);

  // check if there is one left over and just copy it
  if (mPopulationSize % 2 > 0)
    *mIndividual[2 * mPopulationSize - 1] = *mIndividual[mPopulationSize - 1];

  // mutate the offspring
  for (i = mPopulationSize; i < 2 * mPopulationSize && Continue; i++)
    {
      mutate(*mIndividual[i]);
      Continue &= evaluate(*mIndividual[i]);
      mValue[i] = mEvaluationValue;
    }

  return Continue;
}

// select mPopulationSize individuals
bool COptMethodGA::select()
{
  unsigned C_INT32 i, j, nopp, opp;
  unsigned C_INT32 TotalPopulation = 2 * mPopulationSize;

  // tournament competition
  mLosses = 0; // Set all wins to 0.

  // compete with ~ 20% of the TotalPopulation
  nopp = std::max<unsigned C_INT32>(1, mPopulationSize / 5);

  // parents and offspring are all in competition
  for (i = 0; i < TotalPopulation; i++)
    for (j = 0; j < nopp; j++)
      {
        // get random opponent
        opp = mpRandom->getRandomU(TotalPopulation - 1);

        if (mValue[i] < mValue[opp])
          mLosses[opp]++;
        else
          mLosses[i]++;
      }

  // selection of top mPopulationSize winners
  partialSortWithPivot(mLosses.array(),
                       mLosses.array() + mPopulationSize,
                       mLosses.array() + TotalPopulation,
                       mPivot);

  FSwapClass<COptMethodGA, unsigned C_INT32, bool> Swap(this, &COptMethodGA::swap);
  applyPartialPivot(mPivot, mPopulationSize, Swap);

  return true;
}

// check the best individual at this generation
unsigned C_INT32 COptMethodGA::fittest()
{
  unsigned C_INT32 i, BestIndex = C_INVALID_INDEX;
  C_FLOAT64 BestValue = DBL_MAX;

  for (i = 0; i < mPopulationSize && !mLosses[i]; i++)
    if (mValue[i] < BestValue)
      {
        BestIndex = i;
        BestValue = mValue[i];
      }

  return BestIndex;
}

// initialise the population
bool COptMethodGA::creation(unsigned C_INT32 first,
                            unsigned C_INT32 last)
{
  unsigned C_INT32 Last = std::min(last, mPopulationSize);

  unsigned C_INT32 i;
  unsigned C_INT32 j;

  C_FLOAT64 mn;
  C_FLOAT64 mx;
  C_FLOAT64 la;

  bool Continue = true;

  for (i = first; i < Last && Continue; i++)
    {
      for (j = 0; j < mVariableSize; j++)
        {
          // calculate lower and upper bounds
          COptItem & OptItem = *(*mpOptItem)[j];
          mn = *OptItem.getLowerBoundValue();
          mx = *OptItem.getUpperBoundValue();

          C_FLOAT64 & mut = (*mIndividual[i])[j];

          try
            {
              // determine if linear or log scale
              if ((mn < 0.0) || (mx <= 0.0))
                mut = mn + mpRandom->getRandomCC() * (mx - mn);
              else
                {
                  la = log10(mx) - log10(std::max(mn, DBL_MIN));

                  if (la < 1.8)
                    mut = mn + mpRandom->getRandomCC() * (mx - mn);
                  else
                    mut = pow(10.0, log10(std::max(mn, DBL_MIN)) + la * mpRandom->getRandomCC());
                }
            }

          catch (...)
            {
              mut = (mx + mn) * 0.5;
            }

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

      // calculate its fitness
      Continue &= evaluate(*mIndividual[i]);
      mValue[i] = mEvaluationValue;
    }

  return Continue;
}

void COptMethodGA::initObjects()
{
  addObjectReference("Current Generation", mGeneration, CCopasiObject::ValueInt);
}

bool COptMethodGA::initialize()
{
  cleanup();

  unsigned C_INT32 i;

  if (!COptMethod::initialize())
    {
      if (mpCallBack)
        mpCallBack->finishItem(mhGenerations);

      return false;
    }

  mGenerations = * getValue("Number of Generations").pUINT;
  mGeneration = 0;

  if (mpCallBack)
    mhGenerations =
      mpCallBack->addItem("Current Generation",
                          CCopasiParameter::UINT,
                          & mGeneration,
                          & mGenerations);

  mGeneration++;

  mPopulationSize = * getValue("Population Size").pUINT;
  mpRandom =
    CRandom::createGenerator(* (CRandom::Type *) getValue("Random Number Generator").pUINT,
                             * getValue("Seed").pUINT);

  mVariableSize = mpOptItem->size();

  mIndividual.resize(2*mPopulationSize);

  for (i = 0; i < 2*mPopulationSize; i++)
    mIndividual[i] = new CVector< C_FLOAT64 >(mVariableSize);

  mCrossOverFalse.resize(mVariableSize);
  mCrossOverFalse = false;
  mCrossOver.resize(mVariableSize);

  mValue.resize(2*mPopulationSize);
  mValue = std::numeric_limits<C_FLOAT64>::infinity();
  mBestValue = std::numeric_limits<C_FLOAT64>::infinity();

  mShuffle.resize(mPopulationSize);

  for (i = 0; i < mPopulationSize; i++)
    mShuffle[i] = i;

  mLosses.resize(2*mPopulationSize);

  // Initialize the variance for mutations
  mMutationVarians = 0.1;

  return true;
}

bool COptMethodGA::cleanup()
{
  unsigned C_INT32 i;

  pdelete(mpRandom);

  for (i = 0; i < mIndividual.size(); i++)
    pdelete(mIndividual[i]);

  return true;
}

bool COptMethodGA::optimise()
{
  bool Continue = true;

  if (!initialize())
    {
      if (mpCallBack)
        mpCallBack->finishItem(mhGenerations);

      return false;
    }

  // Counters to determine whether the optimization process has stalled
  // They count the number of generations without advances.
  unsigned C_INT32 Stalled, Stalled10, Stalled30, Stalled50;
  Stalled = Stalled10 = Stalled30 = Stalled50 = 0;

  unsigned C_INT32 i;

  // initialise the population
  // first individual is the initial guess
  for (i = 0; i < mVariableSize; i++)
    {
      C_FLOAT64 & mut = (*mIndividual[0])[i];
      COptItem & OptItem = *(*mpOptItem)[i];

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
      (*(*mpSetCalculateVariable)[i])(mut);
    }

  Continue &= evaluate(*mIndividual[0]);
  mValue[0] = mEvaluationValue;

  if (!isnan(mEvaluationValue))
    {
      // and store that value
      mBestValue = mValue[0];
      Continue &= mpOptProblem->setSolution(mBestValue, *mIndividual[0]);

      // We found a new best value lets report it.
      mpParentTask->output(COutputInterface::DURING);
    }

  // the others are random
  Continue &= creation(1, mPopulationSize);

  Continue &= select();
  mBestIndex = fittest();

  if (mBestIndex != C_INVALID_INDEX &&
      mValue[mBestIndex] < mBestValue)
    {
      // and store that value
      mBestValue = mValue[mBestIndex];
      Continue = mpOptProblem->setSolution(mBestValue, *mIndividual[mBestIndex]);

      // We found a new best value lets report it.
      mpParentTask->output(COutputInterface::DURING);
    }

  if (!Continue)
    {
      if (mpCallBack)
        mpCallBack->finishItem(mhGenerations);

      cleanup();
      return false;
    }

  // ITERATE FOR gener GENERATIONS
  for (mGeneration = 2;
       mGeneration <= mGenerations && Continue;
       mGeneration++, Stalled++, Stalled10++, Stalled30++, Stalled50++)
    {
      // perturb the population if we have stalled for a while
      if (Stalled > 50 && Stalled50 > 50)
        {
          Continue &= creation((unsigned C_INT32)(mPopulationSize / 2),
                               mPopulationSize);
          Stalled10 = Stalled30 = Stalled50 = 0;
        }
      else if (Stalled > 30 && Stalled30 > 30)
        {
          Continue &= creation((unsigned C_INT32)(mPopulationSize * 0.7),
                               mPopulationSize);
          Stalled10 = Stalled30 = 0;
        }
      else if (Stalled > 10 && Stalled10 > 10)
        {
          Continue &= creation((unsigned C_INT32)(mPopulationSize * 0.9),
                               mPopulationSize);
          Stalled10 = 0;
        }
      // replicate the individuals
      else
        Continue &= replicate();

      // select the most fit
      Continue &= select();

      // get the index of the fittest
      mBestIndex = fittest();

      if (mBestIndex != C_INVALID_INDEX &&
          mValue[mBestIndex] < mBestValue)
        {
          Stalled = Stalled10 = Stalled30 = Stalled50 = 0;
          mBestValue = mValue[mBestIndex];

          Continue &= mpOptProblem->setSolution(mBestValue, *mIndividual[mBestIndex]);

          // We found a new best value lets report it.
          //if (mpReport) mpReport->printBody();
          mpParentTask->output(COutputInterface::DURING);
        }

      if (mpCallBack)
        Continue &= mpCallBack->progressItem(mhGenerations);
    }

  if (mpCallBack)
    mpCallBack->finishItem(mhGenerations);

  cleanup();

  return Continue;
}
