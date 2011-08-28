/* Begin CVS Header
$Source: /fs/turing/cvs/copasi_dev/copasi/parameterFitting/CExperimentSet.cpp,v $
$Revision: 1.33 $
$Name: Build-33 $
$Author: shoops $
$Date: 2010/07/16 19:01:59 $
End CVS Header */

// Copyright (C) 2008 by Pedro Mendes, Virginia Tech Intellectual
// Properties, Inc., EML Research, gGmbH, University of Heidelberg,
// and The University of Manchester.
// All rights reserved.

// Copyright (C) 2001 - 2007 by Pedro Mendes, Virginia Tech Intellectual
// Properties, Inc. and EML Research, gGmbH.
// All rights reserved.

#include <algorithm>
#include <limits>
#include <math.h>

#include "copasi.h"

#include "CExperimentSet.h"
#include "CExperiment.h"

#include "CopasiDataModel/CCopasiDataModel.h"
#include "report/CCopasiRootContainer.h"
#include "report/CKeyFactory.h"
#include "utilities/utility.h"

CExperimentSet::CExperimentSet(const CCopasiContainer * pParent,
                               const std::string & name):
    CCopasiParameterGroup(name, pParent, "CExperimentSet"),
    mpExperiments(NULL),
    mNonExperiments(0)
{initializeParameter();}

CExperimentSet::CExperimentSet(const CExperimentSet & src,
                               const CCopasiContainer * pParent):
    CCopasiParameterGroup(src, (pParent != NULL) ? pParent : src.getObjectDataModel()),
    mpExperiments(NULL),
    mNonExperiments(0)
{initializeParameter();}

CExperimentSet::CExperimentSet(const CCopasiParameterGroup & group,
                               const CCopasiContainer * pParent):
    CCopasiParameterGroup(group, (pParent != NULL) ? pParent : group.getObjectDataModel()),
    mpExperiments(NULL),
    mNonExperiments(0)
{initializeParameter();}

CExperimentSet::~CExperimentSet() {}

void CExperimentSet::initializeParameter()
{elevateChildren();}

bool CExperimentSet::elevateChildren()
{
  index_iterator it = mValue.pGROUP->begin();
  index_iterator end = mValue.pGROUP->end();

  for (; it != end; ++it)
    {
      if (dynamic_cast< CCopasiParameterGroup * >(*it) == NULL) continue;

      if (!elevate<CExperiment, CCopasiParameterGroup>(*it)) return false;
    }

  mpExperiments = static_cast<std::vector<CExperiment * > * >(mValue.pVOID);

  sort();

  return true;
}

bool CExperimentSet::compile(const std::vector< CCopasiContainer * > listOfContainer)
{
  bool success = true;

  // First we need to sort the experiments so that we can make use of continued
  // file reading.
  sort();

  std::set< CCopasiObject * > DependentObjects;

  std::ifstream in;
  std::string CurrentFileName("");
  unsigned C_INT32 CurrentLineNumber = 1;

  std::vector< CExperiment * >::iterator it = mpExperiments->begin() + mNonExperiments;
  std::vector< CExperiment * >::iterator end = mpExperiments->end();

  for (; it != end; ++it)
    {
      if (CurrentFileName != (*it)->getFileName())
        {
          CurrentFileName = (*it)->getFileName();
          CurrentLineNumber = 1;

          if (in.is_open())
            {
              in.close();
              in.clear();
            }

          in.open(utf8ToLocale(CurrentFileName).c_str(), std::ios::binary);

          if (in.fail())
            {
              CCopasiMessage(CCopasiMessage::ERROR, MCFitting + 8, CurrentFileName.c_str());
              return false; // File can not be opened.
            }
        }

      if (!(*it)->read(in, CurrentLineNumber)) return false;

      if (!(*it)->compile(listOfContainer)) return false;

      const std::map< CCopasiObject *, unsigned C_INT32 > & ExpDependentObjects
      = (*it)->getDependentObjects();
      std::map< CCopasiObject *, unsigned C_INT32 >::const_iterator itObject
      = ExpDependentObjects.begin();
      std::map< CCopasiObject *, unsigned C_INT32 >::const_iterator endObject
      = ExpDependentObjects.end();

      for (; itObject != endObject; ++itObject)
        DependentObjects.insert(itObject->first);
    }

  mDependentObjects.resize(DependentObjects.size());
  CCopasiObject ** ppInsert = mDependentObjects.array();
  std::set< CCopasiObject * >::const_iterator itObject = DependentObjects.begin();
  std::set< CCopasiObject * >::const_iterator endObject = DependentObjects.end();

  for (; itObject != endObject; ++itObject, ++ppInsert)
    *ppInsert = *itObject;

  // Allocation and initialization of statistical information
  mDependentObjectiveValues.resize(mDependentObjects.size());
  mDependentObjectiveValues = std::numeric_limits<C_FLOAT64>::quiet_NaN();

  mDependentRMS.resize(mDependentObjects.size());
  mDependentRMS = std::numeric_limits<C_FLOAT64>::quiet_NaN();

  mDependentErrorMean.resize(mDependentObjects.size());
  mDependentErrorMean = std::numeric_limits<C_FLOAT64>::quiet_NaN();

  mDependentErrorMeanSD.resize(mDependentObjects.size());
  mDependentErrorMeanSD = std::numeric_limits<C_FLOAT64>::quiet_NaN();

  mDependentDataCount.resize(mDependentObjects.size());
  mDependentDataCount = std::numeric_limits<unsigned C_INT32>::quiet_NaN();

  return success;
}

bool CExperimentSet::calculateStatistics()
{
  mDependentObjectiveValues.resize(mDependentObjects.size());
  mDependentObjectiveValues = 0.0;

  mDependentRMS.resize(mDependentObjects.size());
  mDependentRMS = 0.0;

  mDependentErrorMean.resize(mDependentObjects.size());
  mDependentErrorMean = 0.0;

  mDependentErrorMeanSD.resize(mDependentObjects.size());
  mDependentErrorMeanSD = 0.0;

  mDependentDataCount.resize(mDependentObjects.size());
  mDependentDataCount = 0;

  // calclate the per experiment and per dependent value statistics.
  std::vector< CExperiment * >::iterator it = mpExperiments->begin() + mNonExperiments;
  std::vector< CExperiment * >::iterator end = mpExperiments->end();

  unsigned C_INT32 i, Count;
  C_FLOAT64 Tmp;

  for (; it != end; ++it)
    {
      (*it)->calculateStatistics();

      CCopasiObject *const* ppObject = mDependentObjects.array();
      CCopasiObject *const* ppEnd = ppObject + mDependentObjects.size();

      for (i = 0; ppObject != ppEnd; ++ppObject, ++i)
        {
          Count = (*it)->getCount(*ppObject);

          if (Count)
            {
              mDependentObjectiveValues[i] += (*it)->getObjectiveValue(*ppObject);

              Tmp = (*it)->getRMS(*ppObject);
              mDependentRMS[i] += Tmp * Tmp * Count;

              mDependentErrorMean[i] += (*it)->getErrorMean(*ppObject);

              mDependentDataCount[i] += Count;
            }
        }
    }

  unsigned C_INT32 imax = mDependentObjects.size();

  for (i = 0; i != imax; i++)
    {
      Count = mDependentDataCount[i];

      if (Count)
        {
          mDependentRMS[i] = sqrt(mDependentRMS[i] / Count);
          mDependentErrorMean[i] /= Count;
        }
      else
        {
          mDependentRMS[i] = std::numeric_limits<C_FLOAT64>::quiet_NaN();
          mDependentErrorMean[i] = std::numeric_limits<C_FLOAT64>::quiet_NaN();
        }
    }

  it = mpExperiments->begin() + mNonExperiments;

  // We need to loop again to calculate the std. deviation.
  for (; it != end; ++it)
    {
      CCopasiObject *const* ppObject = mDependentObjects.array();
      CCopasiObject *const* ppEnd = ppObject + mDependentObjects.size();

      for (i = 0; ppObject != ppEnd; ++ppObject, ++i)
        {
          Count = (*it)->getCount(*ppObject);

          if (Count)
            mDependentErrorMeanSD[i] +=
              (*it)->getErrorMeanSD(*ppObject, mDependentErrorMean[i]);
        }
    }

  for (i = 0; i != imax; i++)
    {
      Count = mDependentDataCount[i];

      if (Count)
        mDependentErrorMeanSD[i] = sqrt(mDependentErrorMeanSD[i] / Count);
      else
        mDependentErrorMeanSD[i] = std::numeric_limits<C_FLOAT64>::quiet_NaN();
    }

  // This is the time to call the output handler to plot the fitted points.
  for (it = mpExperiments->begin() + mNonExperiments, imax = 0; it != end; ++it)
    imax = std::max(imax, (*it)->getDependentData().numRows());

  CCopasiTask * pParentTask = dynamic_cast< CCopasiTask *>(getObjectAncestor("Task"));
  assert(pParentTask != NULL);

  for (i = 0; i < imax; i++)
    {
      for (it = mpExperiments->begin() + mNonExperiments; it != end; ++it)
        (*it)->updateFittedPointValues(i);

      pParentTask->output(COutputInterface::AFTER);
    }

  return true;
}

const CVector< CCopasiObject * > & CExperimentSet::getDependentObjects() const
{return mDependentObjects;}

const CVector< C_FLOAT64 > & CExperimentSet::getDependentObjectiveValues() const
{return mDependentObjectiveValues;}

const CVector< C_FLOAT64 > & CExperimentSet::getDependentRMS() const
{return mDependentRMS;}

const CVector< C_FLOAT64 > & CExperimentSet::getDependentErrorMean() const
{return mDependentErrorMean;}

const CVector< C_FLOAT64 > & CExperimentSet::getDependentErrorMeanSD() const
{return mDependentErrorMeanSD;}

unsigned C_INT32 CExperimentSet::getExperimentCount() const
{return size() - mNonExperiments;}

CExperiment * CExperimentSet::addExperiment(const CExperiment & experiment)
{
  // We need to make sure that the experiment name is unique.
  std::string name = experiment.getObjectName();

  int i = 0;

  while (getParameter(name))
    {
      i++;
      name = StringPrint("%s_%d", experiment.getObjectName().c_str(), i);
    }

  CExperiment * pExperiment = new CExperiment(experiment);
  pExperiment->setObjectName(name);
  addParameter(pExperiment);

  sort();

  return pExperiment;
}

void CExperimentSet::removeExperiment(const unsigned C_INT32 & index)
{removeParameter(index + mNonExperiments);}

CExperiment * CExperimentSet::getExperiment(const unsigned C_INT32 & index)
{return (*mpExperiments)[index + mNonExperiments];}

const CExperiment * CExperimentSet::getExperiment(const unsigned C_INT32 & index) const
{return (*mpExperiments)[index + mNonExperiments];}

CExperiment * CExperimentSet::getExperiment(const std::string & name)
{return static_cast<CExperiment *>(getGroup(name));}

const CExperiment * CExperimentSet::getExperiment(const std::string & name) const
{return static_cast<const CExperiment *>(getGroup(name));}

bool CExperimentSet::hasDataForTaskType(const CCopasiTask::Type & type) const
{
  std::vector< CExperiment * >::const_iterator it = mpExperiments->begin() + mNonExperiments;
  std::vector< CExperiment * >::const_iterator end = mpExperiments->end();

  for (; it != end; ++it)
    {
      if ((*it)->getExperimentType() == type)
        {
          return true;
        }
    }

  return false;
}

const CCopasiTask::Type & CExperimentSet::getExperimentType(const unsigned C_INT32 & index) const
{return getExperiment(index)->getExperimentType();}

const CMatrix< C_FLOAT64 > & CExperimentSet::getIndependentData(const unsigned C_INT32 & index) const
{return getExperiment(index)->getIndependentData();}

const CMatrix< C_FLOAT64 > & CExperimentSet::getDependentData(const unsigned C_INT32 & index) const
{return getExperiment(index)->getDependentData();}

unsigned C_INT32 CExperimentSet::keyToIndex(const std::string & key) const
{
  const CExperiment * pExp = dynamic_cast<const CExperiment *>(CCopasiRootContainer::getKeyFactory()->get(key));

  if (!pExp) return C_INVALID_INDEX;

  unsigned C_INT32 i, imax = size();

  for (i = 0; i < imax; i++)
    if (pExp == getExperiment(i)) return i;

  return C_INVALID_INDEX;
}

void CExperimentSet::sort()
{
  // First we make sure that all experiments are at the end of the group
  index_iterator it = beginIndex();
  index_iterator end = endIndex();

  index_iterator swapTarget = beginIndex();
  mNonExperiments = 0;

  for (; it != end; ++it)
    if (dynamic_cast< CExperiment * >(*it) == NULL)
      {
        if (it != swapTarget)
          swap(it, swapTarget);

        swapTarget++;
        mNonExperiments++;
      }

  // Now sort the experiments
  std::vector< CExperiment * >::iterator startSort = mpExperiments->begin() + mNonExperiments;
  std::vector< CExperiment * >::iterator endSort = mpExperiments->end();

  std::sort(startSort, endSort, &CExperiment::compare);

  return;
}

std::vector< std::string > CExperimentSet::getFileNames() const
{
  std::vector< std::string > List;
  std::string currentFile = "";

  std::vector< CExperiment * >::iterator it = mpExperiments->begin() + mNonExperiments;
  std::vector< CExperiment * >::iterator end = mpExperiments->end();

  for (; it != end; ++it)
    if (currentFile != (*it)->getFileName())
      {
        currentFile = (*it)->getFileName();
        List.push_back(currentFile);
      }

  return List;
}

unsigned C_INT32 CExperimentSet::getDataPointCount() const
{
  unsigned C_INT32 Count = 0;
  std::vector< CExperiment * >::iterator it = mpExperiments->begin() + mNonExperiments;
  std::vector< CExperiment * >::iterator end = mpExperiments->end();

  for (; it != end; ++it)
    Count += (*it)->getDependentData().numRows() * (*it)->getDependentData().numCols();

  return Count;
}

