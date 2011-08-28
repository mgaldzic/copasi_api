/* Begin CVS Header
$Source: /fs/turing/cvs/copasi_dev/copasi/elementaryFluxModes/CEFMAlgorithm.cpp,v $
$Revision: 1.28 $
$Name: Build-33 $
$Author: gauges $
$Date: 2010/09/07 09:16:21 $
End CVS Header */

// Copyright (C) 2008 by Pedro Mendes, Virginia Tech Intellectual
// Properties, Inc., EML Research, gGmbH, University of Heidelberg,
// and The University of Manchester.
// All rights reserved.

// Copyright (C) 2001 - 2007 by Pedro Mendes, Virginia Tech Intellectual
// Properties, Inc. and EML Research, gGmbH.
// All rights reserved.

/**
 *  CEFMAlgorithm class.
 *  Used to calculate elementary flux modes
 *
 *  Created for COPASI by Stefan Hoops 2002-05-08
 * (C) Stefan Hoops 2002
 */

#include "copasi.h"

#include "CEFMAlgorithm.h"
#include "CEFMProblem.h"
#include "CEFMTask.h"
#include "CFluxMode.h"
#include "CTableauMatrix.h"

#include "model/CModel.h"
#include "utilities/CProcessReport.h"
#include "report/CCopasiObjectReference.h"

CEFMAlgorithm::CSpeciesOrderNode::CSpeciesOrderNode():
    CCopasiNode< size_t >(),
    mTableauLines()
{}

CEFMAlgorithm::CSpeciesOrderNode::CSpeciesOrderNode(const CSpeciesOrderNode & src):
    CCopasiNode< size_t >(src),
    mTableauLines(src.mTableauLines)
{}

CEFMAlgorithm::CSpeciesOrderNode::CSpeciesOrderNode(const size_t & index,
    const CTableauMatrix & matrix):
    CCopasiNode< size_t >(index),
    mTableauLines()
{
  update(matrix);
}

CEFMAlgorithm::CSpeciesOrderNode::~CSpeciesOrderNode()
{}

void CEFMAlgorithm::CSpeciesOrderNode::update(const CTableauMatrix & matrix)
{
  mTableauLines.clear();

  std::list< const CTableauLine * >::const_iterator it = matrix.begin();
  std::list< const CTableauLine * >::const_iterator end = matrix.end();
  size_t TableauLineIndex = 0;

  for (; it != end; ++it, ++TableauLineIndex)
    {
      if ((*it)->getMultiplier(TableauLineIndex) != 0.0)
        {
          mTableauLines.push_back(TableauLineIndex);
        }
    }
}

CEFMAlgorithm::CEFMAlgorithm(const CCopasiContainer * pParent):
    CEFMMethod(CCopasiTask::fluxMode, CCopasiMethod::EFMAlgorithm, pParent),
    mpModel(NULL),
    mStoi(),
    mReversible(0),
    mpCurrentTableau(NULL),
    mpNextTableau(NULL),
    mIndexSet()
{initObjects();}

CEFMAlgorithm::CEFMAlgorithm(const CCopasiMethod::SubType subType, const CCopasiContainer * pParent):
    CEFMMethod(CCopasiTask::fluxMode, subType, pParent),
    mpModel(NULL),
    mStoi(),
    mReversible(0),
    mpCurrentTableau(NULL),
    mpNextTableau(NULL),
    mIndexSet()
{initObjects();}

CEFMAlgorithm::CEFMAlgorithm(const CEFMAlgorithm & src,
                             const CCopasiContainer * pParent):
    CEFMMethod(src, pParent),
    mpModel(NULL),
    mStoi(),
    mReversible(0),
    mpCurrentTableau(NULL),
    mpNextTableau(NULL),
    mIndexSet()
{initObjects();}

CEFMAlgorithm::~CEFMAlgorithm()
{
  DESTRUCTOR_TRACE;
  pdelete(mpCurrentTableau);
  pdelete(mpNextTableau);
}

void CEFMAlgorithm::initObjects()
{
  addObjectReference("Current Step", mStep, CCopasiObject::ValueInt);
}

bool CEFMAlgorithm::initialize()
{
  if (!CEFMMethod::initialize())
    {
      return false;
    }

  CEFMTask * pTask = dynamic_cast< CEFMTask *>(getObjectParent());

  if (pTask == NULL) return false;

  mpModel = pTask->getProblem()->getModel();

  if (mpModel == NULL) return false;

  mpFluxModes->clear();

  /* ModelStoi is the transpose of the models stoichiometry matrix */
  const CTransposeView< CMatrix< C_FLOAT64 > > ModelStoi(mpModel->getStoi());

  unsigned C_INT32 row, numRows = ModelStoi.numRows();
  unsigned C_INT32 col, numCols = ModelStoi.numCols();

  /* Size the stoichiometry matrix passed to the algorithm */
  mStoi.resize(numRows);
  std::vector< std::vector< C_FLOAT64 > >::iterator it = mStoi.begin();
  std::vector< std::vector< C_FLOAT64 > >::iterator end = mStoi.end();

  for (; it != end; ++it)
    it->resize(numCols);

  /* Get the reactions from the model */
  /* Note: We have as many reactions as we have rows in ModelStoi */
  const CCopasiVectorNS < CReaction > & Reaction = mpModel->getReactions();

  /* Vector to keep track of the rearrangements necessary to put the */
  /* reversible reactions to the top of stoichiometry matrix */
  mpReorderedReactions->resize(numRows);

  /* Reversible reaction counter */
  mReversible = 0;
  unsigned C_INT32 Insert;
  unsigned C_INT32 InsertReversible = 0;
  unsigned C_INT32 InsertIrreversible = numRows - 1;

  /* Build the transpose of the stoichiometry matrix, */

  /* sort reversible reactions to the top, count them */

  /* and keep track of the rearrangement */
  for (row = 0; row < numRows; row++)
    {
      if (Reaction[row]->isReversible())
        {
          Insert = InsertReversible++;
          mReversible++;
        }
      else
        Insert = InsertIrreversible--;

      (*mpReorderedReactions)[Insert] = Reaction[row];

      for (col = 0; col < numCols; col++)
        mStoi[Insert][col] = ModelStoi(row, col);
    }

  mStep = 0;
  mMaxStep = numCols;

  if (mpCallBack)
    mhSteps =
      mpCallBack->addItem("Current Step",
                          CCopasiParameter::UINT,
                          & mStep,
                          & mMaxStep);

  return true;
}

bool CEFMAlgorithm::calculate()
{
  if (!initialize())
    {
      if (mpCallBack)
        mpCallBack->finishItem(mhSteps);

      return false;
    }

  calculateFluxModes();

  return true;
}

void CEFMAlgorithm::calculateFluxModes()
{
  bool Continue = true;

  if (mStoi.size())
    {
      /* initialize the current tableau matrix */
      mpCurrentTableau = new CTableauMatrix(mStoi, mReversible);

      /* Do the iteration */
      mIndexSet.resize(mMaxStep);

      for (mStep = 0; mStep < mMaxStep; mStep++)
        mIndexSet[mStep] = mStep;

      while (findMinimalCombinationIndex() && Continue)
        {
          calculateNextTableau();

          if (mpCallBack)
            Continue &= mpCallBack->progressItem(mhSteps);

          static_cast<CCopasiTask *>(getObjectParent())->output(COutputInterface::DURING);
        }

      /* Build the elementary flux modes to be returned */
      if (Continue)
        buildFluxModes();

      /* Delete the current / final tableau matrix */
      pdelete(mpCurrentTableau);
    }

  if (mpCallBack)
    Continue &= mpCallBack->finishItem(mhSteps);
}

void CEFMAlgorithm::calculateNextTableau()
{
  std::list< const CTableauLine * >::iterator a;
  std::list< const CTableauLine * >::iterator b;
  C_FLOAT64 ma, mb;
#ifdef COPASI_DEBUG
  DebugFile << *mpCurrentTableau << std::endl;
#endif //COPASI_DEBUG

  mpNextTableau = new CTableauMatrix();

  /* Move all lines with zeros in the step column to the new tableau matrix */
  /* and remove them from the current tableau matrix */
  a = mpCurrentTableau->begin();

  bool Continue = true;

  unsigned C_INT32 Counter, MaxCounter, hCounter;

  Counter = 0;
  MaxCounter = mpCurrentTableau->size();

  if (mpCallBack)
    hCounter =
      mpCallBack->addItem("Current Line",
                          CCopasiParameter::UINT,
                          & Counter,
                          & MaxCounter);

  while (a != mpCurrentTableau->end() && Continue)
    if ((*a)->getMultiplier(mStep) == 0.0)
      {
        /* We have to make sure that "a" points to the next element in the */
        /* list after the removal of itself */

        if (a == mpCurrentTableau->begin())
          {
            mpNextTableau->addLine(*a, false);
            mpCurrentTableau->removeLine(a);
            a = mpCurrentTableau->begin();
          }
        else
          {
            /* We have to remember the previous element so that a++ points to */
            /* past the removed one */
            b = a--;
            mpNextTableau->addLine(*b, false);
            mpCurrentTableau->removeLine(b);
            a++;
          }

        Counter++;

        if (mpCallBack)
          Continue &= mpCallBack->progressItem(hCounter);
      }
    else
      a++;

  C_FLOAT64 Sign;

  /* Now we create all linear combinations of the remaining lines in the */
  /* current tableau */
  a = mpCurrentTableau->begin();

  //std::cout << "Tableau Marker A" << std::endl << std::endl;

  while (a != mpCurrentTableau->end() && Continue)
    {
      b = a;
      b++;

      /* We make sure that "mb" is positive */
      mb = (*a)->getMultiplier(mStep);

      if (mb < 0.0)
        {
          mb *= -1.0;
          Sign = 1.0;
        }
      else
        Sign = -1.0;

      //std::cout << "Tableau Marker B" << std::endl << std::endl;

      while (b != mpCurrentTableau->end() && Continue)
        {
          ma = Sign * (*b)->getMultiplier(mStep);

          /* The multiplier "ma" for irreversible reactions must be positive */
          if (ma > 0.0 || (*a)->isReversible())
            mpNextTableau->addLine(new CTableauLine(ma, **a, mb, **b));

          /*CTableauLine * debugLine = new CTableauLine(ma, **a, mb, **b);
          if(debugLine->isReversible())
            std::cout << "Reversible Rxn" << std::endl;
          else std::cout << "Irreversible Rxn" << std::endl;
          std::cout << "Flux Score: " << debugLine->getScore() << "Flux Mode Vector: "
            << debugLine->getFluxMode() << std::endl;*/

          b++;

          if (mpCallBack)
            Continue &= mpCallBack->proceed();
        }

      // We no longer need a since all linear combinations have been build;
      mpCurrentTableau->removeLine(a);
      a = mpCurrentTableau->begin();

      Counter++;

      if (mpCallBack)
        Continue &= mpCallBack->progressItem(hCounter);
    }

  if (mpCallBack)
    Continue &= mpCallBack->finishItem(hCounter);

  /* Assign the next tableau to the current tableau and cleanup */
  pdelete(mpCurrentTableau);

  mpCurrentTableau = mpNextTableau;

  mpNextTableau = NULL;
}

void CEFMAlgorithm::buildFluxModes()
{
  mpFluxModes->clear();

  std::list< const CTableauLine * >::iterator a = mpCurrentTableau->begin();
  std::list< const CTableauLine * >::iterator end = mpCurrentTableau->end();

  while (a != end)
    {
      mpFluxModes->push_back(CFluxMode(*a));
      a++;
    }
}

bool CEFMAlgorithm::findMinimalCombinationIndex()
{
  double minCombine = std::numeric_limits<double>::infinity();
  double combine = 0;
  unsigned C_INT32 minIndex = 0;
  unsigned C_INT32 counter;

  if (mIndexSet.size() == 0)
    return false;
  else if (mIndexSet.size() == 1)
    {
      mStep = mIndexSet[0];
      mIndexSet.pop_back();
      return true;
    }

  for (counter = 0; counter < mIndexSet.size(); counter++)
    {
      combine = calculateCombinations(mIndexSet[counter]);

      if (combine < minCombine)
        {
          minCombine = combine;
          minIndex = counter;
        }

      if (combine == 0)
        break;
    }

  mStep = mIndexSet[minIndex];
  mIndexSet.erase(mIndexSet.begin() + minIndex);

  return true;
}

double CEFMAlgorithm::calculateCombinations(unsigned C_INT32 index)
{
  double posIrr = 0;
  double negIrr = 0;
  double rev = 0;
  unsigned int row;

  //Reversible reactions
  std::list< const CTableauLine * >::const_iterator it = mpCurrentTableau->begin();
  std::list< const CTableauLine * >::const_iterator end = mpCurrentTableau->end();

  for (; it != end; ++it)
    {
      if ((*it)->isReversible() && (*it)->getMultiplier(index) != 0.0)
        {
          rev++;
        }
      else if ((*it)->getMultiplier(index) < 0.0)
        {
          negIrr++;
        }
      else if ((*it)->getMultiplier(index) > 0.0)
        {
          posIrr++;
        }
    }

  return (posIrr + rev) *(negIrr + rev);
}
