// Begin CVS Header
//   $Source: /fs/turing/cvs/copasi_dev/copasi/scan/CScanMethod.cpp,v $
//   $Revision: 1.59 $
//   $Name: Build-33 $
//   $Author: shoops $
//   $Date: 2010/02/15 18:18:35 $
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

/**
 *  CScanMethod class.
 *  This class describes the Scan method
 *
 *  Created for Copasi by Rohan Luktuke 2002
 */

#include <string>

#include "mathematics.h"
#include "copasi.h"
#include "model/CModel.h"
#include "model/CState.h"
#include "utilities/CReadConfig.h"
#include "randomGenerator/CRandom.h"
//#include "utilities/CWriteConfig.h"
#include "CScanProblem.h"
#include "CScanMethod.h"
#include "CScanTask.h"
#include "CopasiDataModel/CCopasiDataModel.h"
#include "report/CCopasiRootContainer.h"

// this will have to be defined somewhere else with the
// values of other distribution types
//#define SD_UNIFORM 0
//#define SD_GAUSS 1
//#define SD_BOLTZ 2
//#define SD_REGULAR 3

//**************** CScanItem classes ***************************

//static
CScanItem* CScanItem::createScanItemFromParameterGroup(CCopasiParameterGroup* si,
    CRandom* rg,
    CScanTask* /*st*/)
{
  if (!si) return NULL;

  CScanProblem::Type type = *(CScanProblem::Type*)(si->getValue("Type").pUINT);

  CScanItem* tmp = NULL;

  if (type == CScanProblem::SCAN_REPEAT)
    tmp = new CScanItemRepeat(si);

  if (type == CScanProblem::SCAN_LINEAR)
    tmp = new CScanItemLinear(si);

  if (type == CScanProblem::SCAN_RANDOM)
    tmp = new CScanItemRandom(si, rg);

  /*  if (type == CScanProblem::SCAN_BREAK)
      tmp = new CScanItemBreak(si, st);*/

  return tmp;
}

CScanItem::CScanItem(CCopasiParameterGroup* si)
    : mNumSteps(0),
    mpValue(NULL),
    mStoreValue(0.0),
    mIndex(0),
    mFlagFinished(false)
{
  //TODO: check if the parameters exist. Otherwise: crash
  // also in subclasses...

  mNumSteps = * si->getValue("Number of steps").pUINT;

  std::string tmpString = * si->getValue("Object").pCN;
  CCopasiDataModel* pDataModel = si->getObjectDataModel();
  assert(pDataModel != NULL);
  const CCopasiObject * tmpObject = pDataModel->getObject(tmpString);

  if (!tmpObject) {mpValue = NULL; return;}

  if (!tmpObject->isValueDbl()) {mpValue = NULL; return;}

  mpValue = const_cast<CCopasiObject *>(tmpObject);
}

unsigned C_INT32 CScanItem::getNumSteps() const {return mNumSteps;};

void CScanItem::restoreValue() const
{
  if (mpValue)
    mpValue->setObjectValue(mStoreValue);
};

void CScanItem::storeValue()
{
  if (mpValue)
    {
      assert(mpValue->isValueDbl());
      mStoreValue = * (C_FLOAT64 *) mpValue->getValuePointer();
    }
};

void CScanItem::reset()
{
  mIndex = 0;
  mFlagFinished = false;
  this->step(); //purely virtual
}

bool CScanItem::isFinished() const {return mFlagFinished;};

bool CScanItem::isValidScanItem()
{
  return true;
}

const CCopasiObject * CScanItem::getObject() const
{
  return mpValue;
}
//*******

CScanItemRepeat::CScanItemRepeat(CCopasiParameterGroup* si)
    : CScanItem(si)
{
  if (mNumSteps >= 1)
    --mNumSteps; // for the repeat item mNumSteps is the number of iterations, not of intervals
}

void CScanItemRepeat::step()
{
  //do something ...

  //the index
  if (mIndex > mNumSteps)
    mFlagFinished = true;

  ++mIndex;
}

//*******

CScanItemLinear::CScanItemLinear(CCopasiParameterGroup* si)
    : CScanItem(si),
    mLog(false)
{
  mLog = * si->getValue("log").pBOOL;
  mMin = * si->getValue("Minimum").pDOUBLE;
  mMax = * si->getValue("Maximum").pDOUBLE;

  if (mLog)
    {
      mMin = log(mMin);
      mMax = log(mMax);
    }

  mFaktor = (mMax - mMin) / mNumSteps;

  //TODO: log scanning of negative values?
}

void CScanItemLinear::step()
{
  //do something ...
  C_FLOAT64 Value = mMin + mIndex * mFaktor;

  if (mLog)
    Value = exp(Value);

  //the index
  if (mIndex > mNumSteps)
    mFlagFinished = true;

  if (mpValue) mpValue->setObjectValue(Value);

  ++mIndex;
}

bool CScanItemLinear::isValidScanItem()
{
  if (!CScanItem::isValidScanItem()) return false;

  if (!mpValue)
    {
      CCopasiMessage(CCopasiMessage::EXCEPTION, "Invalid or missing scan parameter.");
      return false;
    }

  if (mLog)
    {
      if (isnan(mFaktor) || mFaktor < - DBL_MAX || DBL_MAX < mFaktor)
        {
          //not a valid range for log
          CCopasiMessage(CCopasiMessage::EXCEPTION, "Only positive values for min and max are possible for a logarithmic scan.");
          return false;
        }
    }

  return true;
}

//*******

CScanItemRandom::CScanItemRandom(CCopasiParameterGroup* si, CRandom* rg)
    : CScanItem(si),
    mRg(rg),
    mRandomType(0),
    mLog(false)
{
  mRandomType = * si->getValue("Distribution type").pUINT;

  mLog = * si->getValue("log").pBOOL;
  mMin = * si->getValue("Minimum").pDOUBLE;
  mMax = * si->getValue("Maximum").pDOUBLE;

  if (mLog)
    {
      mMin = log(mMin);
      mMax = log(mMax);
    }

  mNumSteps = 0;
  mFaktor = (mMax - mMin);
}

void CScanItemRandom::step()
{
  C_FLOAT64 Value;

  //the index
  if (mIndex > mNumSteps)
    mFlagFinished = true;
  else
    {
      C_FLOAT64 tmpF;

      switch (mRandomType)
        {
          case 0:             //uniform
            Value = mMin + mRg->getRandomCC() * mFaktor;

            if (mLog)
              Value = exp(Value);

            break;

          case 1:             //normal
            tmpF = mRg->getRandomNormal01();
            Value = mMin + tmpF * mMax;

            if (mLog)
              Value = exp(Value);

            break;

          case 2:             //poisson
            Value = mRg->getRandomPoisson(mMin);
            //if (mLog)
            //  *mpValue = exp(*mpValue);
            break;
        }
    }

  if (mpValue) mpValue->setObjectValue(Value);

  ++mIndex;
}

bool CScanItemRandom::isValidScanItem()
{
  if (!CScanItem::isValidScanItem()) return false;

  if (!mpValue)
    {
      CCopasiMessage(CCopasiMessage::EXCEPTION, "Invalid or missing scan parameter.");
      return false;
    }

  //   if (mLog)
  //     {
  //       if ((mMin <= 0) || (mMax <= 0))
  //         {
  //           //not a valid range for log
  //           CCopasiMessage(CCopasiMessage::EXCEPTION, "Only positive values for min and max are possible\nfor a logarithmic scan.");
  //           return false;
  //}
  //}

  return true;
}

//*******

/*CScanItemBreak::CScanItemBreak(const CCopasiParameterGroup* si, CScanTask* st)
    : CScanItem(si),
    mPlotB(0),
    mReportB(0),
    mST(NULL)
{
  mReportB = *(unsigned C_INT32*)(si->getValue("Report break"));
  mPlotB = *(unsigned C_INT32*)(si->getValue("Plot break"));
  mST = st;
  mNumSteps = 0;
}

void CScanItemBreak::step()
{
  //the index
  if (mIndex > mNumSteps)
    mFlagFinished = true;
  else
    {
      //TODO: tell the task what exactly to do...
      mST->outputSeparatorCallback();
    }

  ++mIndex;
}*/

//**************** CScanMethod class ***************************

CScanMethod * CScanMethod::createMethod() {return new CScanMethod;}

CScanMethod::CScanMethod():
    CCopasiMethod(CCopasiTask::scan, CCopasiMethod::scanMethod),
    mpProblem(NULL),
    mpTask(NULL),
    mpRandomGenerator(NULL),
    mTotalSteps(1)
    //    mVariableSize(0),
    //    mpVariables(NULL)
{
  //  addParameter("Random Number Generator", CCopasiParameter::STRING,
  //               CRandom::TypeName[1]);
  //  addParameter("Random Number Seed", CCopasiParameter::INT, (C_INT32) 0);
  mpRandomGenerator = CRandom::createGenerator(CRandom::r250);
}

/*CScanMethod::CScanMethod(const CScanMethod & src,
                         const CCopasiContainer * pParent):
    CCopasiMethod(src, pParent),
    mpProblem(NULL),
    mpRandomGenerator(NULL)
    //    mVariableSize(0),
    //    mpVariables(NULL)
{}*/

CScanMethod::~CScanMethod()
{
  cleanupScanItems();
  delete mpRandomGenerator;
  mpRandomGenerator = NULL;
}

bool CScanMethod::cleanupScanItems()
{
  if (!mpProblem) return false;

  unsigned C_INT32 i, imax = mScanItems.size();

  for (i = 0; i < imax; ++i) if (mScanItems[i]) delete mScanItems[i];

  mScanItems.clear();
  return true;
}

bool CScanMethod::init()
{
  if (!mpProblem) return false;

  mpTask = dynamic_cast< CScanTask * >(getObjectParent());

  if (mpTask == NULL) return false;

  cleanupScanItems();
  mInitialRefreshes.clear();
  mTotalSteps = 1;
  std::set< const CCopasiObject * > ObjectSet;

  unsigned C_INT32 i, imax = mpProblem->getNumberOfScanItems();

  for (i = 0; i < imax; ++i)
    {
      mScanItems.push_back(CScanItem::createScanItemFromParameterGroup(mpProblem->getScanItem(i),
                           mpRandomGenerator,
                           (CScanTask*)(getObjectParent())));
      mTotalSteps *= mScanItems[i]->getNumSteps();
      ObjectSet.insert(mScanItems[i]->getObject());
    }

  ObjectSet.erase(NULL);
  mInitialRefreshes = mpProblem->getModel()->buildInitialRefreshSequence(ObjectSet);

  //set mLastNestingItem
  mLastNestingItem = -1;

  if (imax != 0)
    {
      //search from the end
      C_INT32 j;

      for (j = mScanItems.size() - 1; j >= 0; --j)
        {
          if (mScanItems[j]->isNesting())
            {
              mLastNestingItem = j;
              break;
            }
        }
    }

  return true;
}

bool CScanMethod::scan()
{
  if (!mpProblem) return false;

  //a hack to ensure that the first subtask is run with initial conditions
  //pDataModel->getModel()->setState(&mpProblem->getInitialState());

  bool success = true;

  unsigned C_INT32 i, imax = mScanItems.size();

  //store old parameter values
  for (i = 0; i < imax; ++i)
    mScanItems[i]->storeValue();

  //Do the scan...
  if (imax) //there are scan items
    success = loop(0);
  else
    success = calculate(); //nothing to scan, only one call to the subtask

  //restore old parameter values
  for (i = 0; i < imax; ++i)
    mScanItems[i]->restoreValue();

  return success;
}

bool CScanMethod::loop(unsigned C_INT32 level)
{
  bool isLastMasterItem = (level == (mScanItems.size() - 1)); //TODO

  CScanItem* currentSI = mScanItems[level];

  for (currentSI->reset(); !currentSI->isFinished(); currentSI->step())
    {
      //TODO: handle slave SIs

      if (isLastMasterItem)
        {
          if (!calculate()) return false;
        }
      else
        {
          if (!loop(level + 1)) return false;
        } //TODO

      //separator needs to be handled slightly differently if we are at the last item
      if (currentSI->isNesting())
        ((CScanTask*)(getObjectParent()))->outputSeparatorCallback((C_INT32)level == mLastNestingItem);
    }

  return true;
}

bool CScanMethod::calculate()
{
  std::vector< Refresh * >::iterator it = mInitialRefreshes.begin();
  std::vector< Refresh * >::iterator end = mInitialRefreshes.end();

  while (it != end)
    (**it++)();

  return mpTask->processCallback();
}

void CScanMethod::setProblem(CScanProblem * problem)
{mpProblem = problem;}

//virtual
bool CScanMethod::isValidProblem(const CCopasiProblem * pProblem)
{
  if (!CCopasiMethod::isValidProblem(pProblem)) return false;

  const CScanProblem * pP = dynamic_cast<const CScanProblem *>(pProblem);

  if (!pP)
    {
      //not a TrajectoryProblem
      CCopasiMessage(CCopasiMessage::EXCEPTION, "Problem is not a Scan problem.");
      return false;
    }

  unsigned C_INT32 i, imax = pP->getNumberOfScanItems();

  if (imax <= 0)
    {
      //no scan items
      CCopasiMessage(CCopasiMessage::WARNING, "There is nothing to scan.");
      return false;
    }

  for (i = 0; i < imax; ++i)
    {
      CScanItem * si = CScanItem::createScanItemFromParameterGroup(mpProblem->getScanItem(i),
                       mpRandomGenerator,
                       (CScanTask*)(getObjectParent()));

      if (!si)
        {
          //parameter group could not be interpreted
          CCopasiMessage(CCopasiMessage::ERROR, "Internal problem with scan definition.");
          return false;
        }

      if (!si->isValidScanItem())
        {
          //the self check of the scan item failed.
          //the message should be generated by the isValidScanItem() method.
          delete si;
          return false;
        }

      delete si;
      //mTotalSteps *= mScanItems[i]->getNumSteps();
    }

  return true;
}
