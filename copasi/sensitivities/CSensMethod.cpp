// Begin CVS Header
//   $Source: /fs/turing/cvs/copasi_dev/copasi/sensitivities/CSensMethod.cpp,v $
//   $Revision: 1.34 $
//   $Name: Build-33 $
//   $Author: shoops $
//   $Date: 2010/09/02 14:30:56 $
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
 *  CSensMethod class.
 */

#include "copasi.h"

//#include "utilities/CCopasiVector.h"
#include "CopasiDataModel/CCopasiDataModel.h"
#include "report/CCopasiRootContainer.h"
#include "CSensMethod.h"
#include "CSensProblem.h"

#include "model/CModel.h"
#include "model/CState.h"
#include "utilities/CProcessReport.h"

CSensMethod *
CSensMethod::createSensMethod(CCopasiMethod::SubType subType)
{
  CSensMethod * pMethod = NULL;

  switch (subType)
    {
      case unset:
      case sensMethod:
        pMethod = new CSensMethod(subType);
        break;

      default:
        fatalError();
    }

  return pMethod;
}

/**
 *  Default constructor.
 */
CSensMethod::CSensMethod(CCopasiMethod::SubType subType,
                         const CCopasiContainer * pParent):
    CCopasiMethod(CCopasiTask::sens, subType, pParent),
    mpProblem(NULL)
{
  addParameter("Delta factor",
               CCopasiParameter::UDOUBLE, (C_FLOAT64) 1e-3);
  mpDeltaFactor = (C_FLOAT64*)getValue("Delta factor").pUDOUBLE;

  addParameter("Delta minimum",
               CCopasiParameter::UDOUBLE, (C_FLOAT64) 1e-12);
  mpMinDelta = (C_FLOAT64*)getValue("Delta minimum").pUDOUBLE;

  CONSTRUCTOR_TRACE;
}

/**
 *  Copy constructor.
 *  @param "const CSensMethod &" src
 */
CSensMethod::CSensMethod(const CSensMethod & src,
                         const CCopasiContainer * pParent):
    CCopasiMethod(src, pParent),
    mpProblem(src.mpProblem)
{CONSTRUCTOR_TRACE;}

/**
 *  Destructor.
 */
CSensMethod::~CSensMethod()
{DESTRUCTOR_TRACE;}

//***********************************************************************************

bool CSensMethod::do_target_calculation(CCopasiArray & result, bool first)
{
  //perform the necessary updates
  std::vector< Refresh * >::iterator it = mInitialRefreshes.begin();
  std::vector< Refresh * >::iterator end = mInitialRefreshes.end();

  while (it != end)
    (**it++)();

  //****** do subtask ******************
  if (mpSubTask)
    {
      if (mpProblem->getSubTaskType() == CSensProblem::SteadyState)
        mpSubTask->process(/*first*/true);
      else
        mpSubTask->process(true);

      // for steady state calculation only the first calculation is done from
      //initial state, the remaining from the current state
    }
  else
    {
      //mpProblem->getModel()
      //mpProblem->getModel()->updateSimulatedValues();
      //mpProblem->getModel()->updateNonSimulatedValues();
    }

  mpProblem->getModel()->updateSimulatedValues(true);
  mpProblem->getModel()->updateNonSimulatedValues();

  //****** retrieve results ************

  //resize results array
  CCopasiArray::index_type resultindex;
  C_INT32 i, imax = mTargetfunctionPointers.size();

  if (imax > 1)
    resultindex.push_back(imax);

  result.resize(resultindex);

  //copy result
  for (i = 0; i < imax; ++i)
    {
      if (imax > 1)
        resultindex[0] = i;

      result[resultindex] = *(C_FLOAT64*)mTargetfunctionPointers[i]->getValuePointer();
    }

  //progress bar
  ++mProgress;

  if (mpProgressBar)
    {
      bool tmp = mpProgressBar->progressItem(mProgressHandler);
      return tmp;
    }

  return true;
}

C_FLOAT64 CSensMethod::do_variation(CCopasiObject* variable)
{
  C_FLOAT64 value;
  value = *(C_FLOAT64*)variable->getValuePointer();
  C_FLOAT64 delta;
  delta = fabs(value) * *mpDeltaFactor;

  if (delta < *mpMinDelta) delta = *mpMinDelta;

  setValue(variable, delta + value);

  return delta;
}

void CSensMethod::setValue(CCopasiObject* variable, C_FLOAT64 value)
{
  variable->setObjectValue(value);

  if (variable->getObjectName() == "Concentration")
    {
      CMetab* pMetab = dynamic_cast<CMetab*>(variable->getObjectAncestor("Metabolite"));

      if (pMetab)
        {
          pMetab->setConcentration(value);
          pMetab->refreshNumber();
        }
    }
}

void CSensMethod::calculate_difference(unsigned C_INT32 level, const C_FLOAT64 & delta,
                                       CCopasiArray & result, CCopasiArray::index_type & resultindex)
{
  assert(delta != 0.0);
  assert(mLocalData[level].tmp1.size() == mLocalData[level].tmp2.size());
  unsigned C_INT32 dim = mLocalData[level].tmp1.dimensionality();
  assert(resultindex.size() >= dim);

  CCopasiArray::index_type indexmax = mLocalData[level].tmp1.size();

  //init index with zero
  CCopasiArray::index_type indexit; indexit.resize(dim);
  unsigned C_INT32 i;

  for (i = 0; i < dim; ++i)
    indexit[i] = 0;

  //handle scalars separately
  if (dim == 0)
    {
      result[resultindex] = (mLocalData[level].tmp2[indexit] - mLocalData[level].tmp1[indexit]) / delta;
      return;
    }

  //now for all higher dimensionalities
  for (;;)
    {
      //do difference calculation
      for (i = 0; i < dim; ++i)
        resultindex[i] = indexit[i];  //TODO: use stl algorithm

      result[resultindex] = (mLocalData[level].tmp2[indexit] - mLocalData[level].tmp1[indexit]) / delta;

      //increase index
      ++indexit[dim - 1];

      //check overflow
      C_INT32 j;

      for (j = dim - 1; j >= 0; --j)
        {
          if (indexit[j] >= indexmax[j])
            {
              indexit[j] = 0;

              if (j > 0)
                ++indexit[j - 1];
              else
                return;
            }
          else
            break;
        }
    }
}

bool CSensMethod::calculate_one_level(unsigned C_INT32 level, CCopasiArray & result)
{
  //do first calculation
  if (level == 0)
    {
      if (!do_target_calculation(mLocalData[level].tmp1, true)) return false;
    }
  else
    {
      if (!calculate_one_level(level - 1, mLocalData[level].tmp1)) return false;
    }

  //resize results array
  CCopasiArray::index_type resultindex; resultindex = mLocalData[level].tmp1.size();

  if (mLocalData[level].variables.size() > 1)
    resultindex.push_back(mLocalData[level].variables.size());

  result.resize(resultindex);

  //loop over all variables
  C_INT32 i, imax = mLocalData[level].variables.size();

  for (i = 0; i < imax; ++i)
    {
      //store variable value
      C_FLOAT64 store = *(C_FLOAT64 *)mLocalData[level].variables[i]->getValuePointer();

      //change variable
      C_FLOAT64 delta = do_variation(mLocalData[level].variables[i]);

      //do second calculation
      if (level == 0)
        {
          if (!do_target_calculation(mLocalData[level].tmp2, false)) return false;
        }
      else
        {
          if (!calculate_one_level(level - 1, mLocalData[level].tmp2)) return false;
        }

      //restore variable
      //mLocalData[level].variables[i]->setObjectValue(store);
      setValue(mLocalData[level].variables[i], store);

      //calculate derivative
      if (imax > 1)
        resultindex[resultindex.size() - 1] = i;

      calculate_difference(level, delta, result, resultindex);
    }

  return true;
}

//********** SCALING *************************************************************

void CSensMethod::scaling_targetfunction(const C_FLOAT64 & factor,
    CCopasiArray::index_type & resultindex)
{
  unsigned C_INT32 dim = mLocalData[0].tmp1.dimensionality();
  assert(resultindex.size() >= dim);

  CCopasiArray::index_type indexmax = mLocalData[0].tmp1.size();
  //init index with zero
  CCopasiArray::index_type indexit; indexit.resize(dim);
  unsigned C_INT32 i;

  for (i = 0; i < dim; ++i)
    indexit[i] = 0;

  //handle scalars separately
  if (dim == 0)
    {
      mpProblem->getScaledResult()[resultindex] = mpProblem->getResult()[resultindex] * factor / mLocalData[0].tmp1[indexit];
      return;
    }

  //now for all higher dimensionalities
  for (;;)
    {
      for (i = 0; i < dim; ++i)
        resultindex[i] = indexit[i];  //TODO: use stl algorithm

      mpProblem->getScaledResult()[resultindex] = mpProblem->getResult()[resultindex] * factor / mLocalData[0].tmp1[indexit];

      //increase index
      ++indexit[dim - 1];

      //check overflow
      C_INT32 j;

      for (j = dim - 1; j >= 0; --j)
        {
          if (indexit[j] >= indexmax[j])
            {
              indexit[j] = 0;

              if (j > 0)
                ++indexit[j - 1];
              else
                return;
            }
          else
            break;
        }
    }
}

void CSensMethod::scaling_variables(C_INT32 level, const C_FLOAT64 & factor,
                                    CCopasiArray::index_type & resultindex)
{
  //loop over all variables
  C_INT32 i, imax = mLocalData[level].variables.size();

  for (i = 0; i < imax; ++i)
    {
      //get Value
      C_FLOAT64 value = *(C_FLOAT64 *)mLocalData[level].variables[i]->getValuePointer() * factor;

      //do recursive calculation
      if (imax > 1)
        resultindex[mLocalData[level].index] = i;

      if (level == 0)
        {
          scaling_targetfunction(value, resultindex);
        }
      else
        {
          scaling_variables(level - 1, value, resultindex);
        }
    }
}

void CSensMethod::do_scaling()
{
  CCopasiArray::index_type index;
  index.resize(mpProblem->getResult().dimensionality());
  scaling_variables(mLocalData.size() - 1, 1.0, index);
}

//****************************************************************************

C_FLOAT64 CSensMethod::do_collapsing_innerloop(CCopasiArray::index_type & fullindex)
{
  //fullindex[0]=0;
  //return mpProblem->getScaledResult()[fullindex];

  //assumes the sum is to be taken over the first dim of the scaled result array
  C_FLOAT64 tmpFloat, tmpSum = 0;
  unsigned C_INT32 i, imax = mpProblem->getScaledResult().size()[0];

  for (i = 0; i < imax; ++i)
    {
      fullindex[0] = i;
      tmpFloat = mpProblem->getScaledResult()[fullindex];

      if (tmpFloat != tmpFloat) continue;

      if (fabs(tmpFloat) >= DBL_MAX) continue;

      tmpSum += tmpFloat * tmpFloat;
    }

  return sqrt(tmpSum);
}

void CSensMethod::do_collapsing()
{
  if (mpProblem->collapsRequested())
    {
      CCopasiArray::index_type fullresultindex = mpProblem->getScaledResult().size();
      CCopasiArray::index_type collapsedresultindex = mpProblem->getCollapsedResult().size();

      C_INT32 shift = fullresultindex.size() - collapsedresultindex.size();

      if (shift != 1) return; //only supported if target functions list is 1D

      //***** skalar ********
      if (collapsedresultindex.size() == 0)
        {
          mpProblem->getCollapsedResult()[collapsedresultindex] =
            do_collapsing_innerloop(fullresultindex);
          return;
        }

      //***** higher dimensions *********
      unsigned C_INT32 i, dim = collapsedresultindex.size();
      CCopasiArray::index_type indexmax = mpProblem->getCollapsedResult().size();

      //set index to zero
      for (i = 0; i < dim; ++i) collapsedresultindex[i] = 0;

      for (;;)
        {
          fullresultindex[0] = 0;

          for (i = 0; i < dim; ++i)
            fullresultindex[i + shift] = collapsedresultindex[i];

          mpProblem->getCollapsedResult()[collapsedresultindex] =
            do_collapsing_innerloop(fullresultindex);

          //increase index
          ++collapsedresultindex[dim - 1];

          //check overflow
          C_INT32 j;

          for (j = dim - 1; j >= 0; --j)
            {
              if (collapsedresultindex[j] >= indexmax[j])
                {
                  collapsedresultindex[j] = 0;

                  if (j > 0)
                    ++collapsedresultindex[j - 1];
                  else
                    return;
                }
              else
                break;
            }
        }
    }
  else
    {}}

//****************************************************************************

bool CSensMethod::initialize(CSensProblem* problem)
{
  bool success = true;

  mpProblem = problem;
  assert(mpProblem);

  //initialize the target calculation
  mpSubTask = NULL;
  CCopasiDataModel* pDataModel = getObjectDataModel();
  assert(pDataModel != NULL);

  switch (mpProblem->getSubTaskType())
    {
      case CSensProblem::Evaluation:
        mpSubTask = NULL;
        break;

      case CSensProblem::SteadyState:
        mpSubTask = dynamic_cast<CCopasiTask*>
                    ((*pDataModel->getTaskList())["Steady-State"]);
        break;

      case CSensProblem::TimeSeries:
        mpSubTask = dynamic_cast<CCopasiTask*>
                    ((*pDataModel->getTaskList())["Time-Course"]);
        break;

        /*    case CSensProblem::LyapunovExp:
              mpSubTask = dynamic_cast<CCopasiTask*>
                          ((*pDataModel->getTaskList())["Lyapunov Exponents"]);
              break;*/
    }

  if (mpSubTask)
    {
      mpSubTask->getProblem()->setModel(mpProblem->getModel());
      mpSubTask->setCallBack(NULL);
      success &= mpSubTask->initialize(CCopasiTask::NO_OUTPUT, NULL, NULL);
    }

  //initialize the variables pointers
  std::set< const CCopasiObject * > ObjectSet;
  C_INT32 i, imax = mpProblem->getNumberOfVariables();
  mLocalData.resize(imax);

  for (i = 0; i < imax; ++i)
    {
      mLocalData[i].variables = mpProblem->getVariables(i).getVariablesPointerList(pDataModel);

      ObjectSet.insert(mLocalData[i].variables.begin(), mLocalData[i].variables.end());
    }

  //determine which refreshes need to be called when the variables are changed
  ObjectSet.erase(NULL);
  mInitialRefreshes.clear();
  mInitialRefreshes = mpProblem->getModel()->buildInitialRefreshSequence(ObjectSet);

  //initialize the target function pointers
  mTargetfunctionPointers = mpProblem->getTargetFunctions().getVariablesPointerList(pDataModel);

  //****** initialize result annotations ****************

  //determine dimensions of result
  CCopasiArray::index_type s;
  CCopasiArray::index_type sc; //size of collapsed result

  if (mTargetfunctionPointers.size() > 1)
    {
      s.push_back(mTargetfunctionPointers.size());
    }

  for (i = 0; i < imax; ++i)
    {
      if (mLocalData[i].variables.size() > 1)
        {
          mLocalData[i].index = s.size();
          s.push_back(mLocalData[i].variables.size());
          sc.push_back(mLocalData[i].variables.size());
        }
      else
        mLocalData[i].index = -1;
    }

  //resize result & annotations
  mpProblem->getResult().resize(s);
  mpProblem->getResultAnnotated()->resize();
  //mpProblem->getResultAnnotated()->setMode(CArrayAnnotation::OBJECTS);

  mpProblem->getScaledResult().resize(s);
  mpProblem->getScaledResultAnnotated()->resize();
  //mpProblem->getScaledResultAnnotated()->setMode(CArrayAnnotation::OBJECTS);

  if (mpProblem->collapsRequested())
    {
      mpProblem->getCollapsedResult().resize(sc);
      mpProblem->getCollapsedResultAnnotated()->resize();
      //mpProblem->getCollapsedResultAnnotated()->setMode(CArrayAnnotation::OBJECTS);
    }

  unsigned C_INT32 dim = 0;
  unsigned C_INT32 j;

  //target function annotations //TODO: only implemented for scalar and vector
  if (mTargetfunctionPointers.size() > 1)
    {
      std::ostringstream tmp;
      tmp << "Target functions, " << mpProblem->getTargetFunctions().getListTypeDisplayName();
      mpProblem->getResultAnnotated()->setDimensionDescription(dim, tmp.str());
      mpProblem->getScaledResultAnnotated()->setDimensionDescription(dim, tmp.str());

      for (j = 0; j < mTargetfunctionPointers.size(); ++j)
        {
          mpProblem->getResultAnnotated()->setAnnotationCN(dim, j, mTargetfunctionPointers[j]->getCN());
          mpProblem->getScaledResultAnnotated()->setAnnotationCN(dim, j, mTargetfunctionPointers[j]->getCN());
        }

      ++dim;
    }

  //variables annotiation
  unsigned C_INT32 dim2 = 0; //for collapsed result

  for (i = 0; i < imax; ++i)
    {
      if (mLocalData[i].variables.size() > 1)
        {
          std::ostringstream tmp;
          tmp << "Variables " << i + 1 << ", " << mpProblem->getVariables(i).getListTypeDisplayName();
          mpProblem->getResultAnnotated()->setDimensionDescription(dim, tmp.str());
          mpProblem->getScaledResultAnnotated()->setDimensionDescription(dim, tmp.str());

          if (mpProblem->collapsRequested())
            mpProblem->getCollapsedResultAnnotated()->setDimensionDescription(dim2, tmp.str());

          for (j = 0; j < mLocalData[i].variables.size(); ++j)
            {
              mpProblem->getResultAnnotated()->setAnnotationCN(dim, j, mLocalData[i].variables[j]->getCN());
              mpProblem->getScaledResultAnnotated()->setAnnotationCN(dim, j, mLocalData[i].variables[j]->getCN());

              if (mpProblem->collapsRequested())
                mpProblem->getCollapsedResultAnnotated()->setAnnotationCN(dim2, j, mLocalData[i].variables[j]->getCN());
            }

          ++dim;
          ++dim2;
        }
    }

  return success;
}

C_INT32 CSensMethod::getNumberOfSubtaskCalculations()
{
  C_INT32 ret = 1;
  unsigned C_INT32 i;

  for (i = 0; i < mLocalData.size(); ++i)
    {
      ret *= mLocalData[i].variables.size() + 1;
    }

  return ret;
}

bool CSensMethod::process(CProcessReport * handler)
{
  if (!mLocalData.size()) return false;

  //initialize progress bar
  mpProgressBar = handler;

  if (mpProgressBar)
    {
      mpProgressBar->setName("performing sensitivities calculation...");
      C_INT32 max = getNumberOfSubtaskCalculations();
      mProgress = 0;
      mProgressHandler = mpProgressBar->addItem("Completion",
                         CCopasiParameter::INT,
                         &mProgress, &max);

      if (mpSubTask)
        mpSubTask->setCallBack(mpProgressBar);
    }

  if (!calculate_one_level(mLocalData.size() - 1, mpProblem->getResult())) return false;

  do_scaling();

  do_collapsing();

  if (mpProgressBar) mpProgressBar->finishItem(mProgressHandler);

  return true;
}

//virtual
bool CSensMethod::isValidProblem(const CCopasiProblem * pProblem)
{
  if (!CCopasiMethod::isValidProblem(pProblem)) return false;

  const CSensProblem * pP = dynamic_cast<const CSensProblem *>(pProblem);

  if (!pP)
    {
      CCopasiMessage(CCopasiMessage::EXCEPTION, "Problem is not a sensitivities problem.");
      return false;
    }

  return true;

  //all sizes at least one

  //dimension of variables 0 or 1

  //if target is scan make sure the scan subtask is not sens.
}
