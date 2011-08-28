// Begin CVS Header
//   $Source: /fs/turing/cvs/copasi_dev/copasi/scan/CScanTask.cpp,v $
//   $Revision: 1.81 $
//   $Name: Build-33 $
//   $Author: shoops $
//   $Date: 2010/09/02 14:30:55 $
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
 * CScanTask class.
 *
 * This class implements a scan task which is comprised of a
 * of a problem and a method.
 *
 */
#include "copasi.h"
#include "CScanTask.h"
#include "CScanProblem.h"
#include "CScanMethod.h"
#include "utilities/CReadConfig.h"
#include "report/CKeyFactory.h"
#include "report/CReport.h"

#include "model/CModel.h"
#include "model/CState.h"

#include "trajectory/CTrajectoryTask.h"
#include "trajectory/CTrajectoryProblem.h"
#include "steadystate/CSteadyStateTask.h"
#include "steadystate/CSteadyStateProblem.h"
#include "utilities/COutputHandler.h"
#include "utilities/CProcessReport.h"
#include "CopasiDataModel/CCopasiDataModel.h"
#include "report/CCopasiRootContainer.h"

CScanTask::CScanTask(const CCopasiContainer * pParent):
    CCopasiTask(CCopasiTask::scan, pParent)
{
  mpProblem = new CScanProblem(this);
  mpMethod = CScanMethod::createMethod();
  this->add(mpMethod, true);
  static_cast< CScanMethod * >(mpMethod)->setProblem(static_cast< CScanProblem * >(mpProblem));
}

CScanTask::CScanTask(const CScanTask & src,
                     const CCopasiContainer * pParent):
    CCopasiTask(src, pParent)
{
  mpProblem = new CScanProblem(*(CScanProblem *) src.mpProblem, this);
  mpMethod = CScanMethod::createMethod();
  this->add(mpMethod, true);
  static_cast< CScanMethod * >(mpMethod)->setProblem(static_cast< CScanProblem * >(mpProblem));
}

CScanTask::~CScanTask()
{cleanup();}

void CScanTask::cleanup()
{}

bool CScanTask::initialize(const OutputFlag & of,
                           COutputHandler * pOutputHandler,
                           std::ostream * pOstream)
{
  assert(mpProblem && mpMethod);

  mpMethod->isValidProblem(mpProblem);

  bool success = true;

  if ((of & REPORT) &&
      pOutputHandler != NULL)
    {
      if (mReport.open(getObjectDataModel(), pOstream))
        pOutputHandler->addInterface(&mReport);
      else
        CCopasiMessage(CCopasiMessage::COMMANDLINE, MCCopasiTask + 5, getObjectName().c_str());
    }

  success &= initSubtask(of, pOutputHandler, mReport.getStream());
  success &= CCopasiTask::initialize(of, pOutputHandler, mReport.getStream());

  return success;
}

void CScanTask::load(CReadConfig & C_UNUSED(configBuffer))
{}

bool CScanTask::process(const bool & useInitialValues)
{
  if (!mpProblem) fatalError();

  if (!mpMethod) fatalError();

  //mpMethod->isValidProblem(mpProblem);

  CScanProblem * pProblem = dynamic_cast<CScanProblem *>(mpProblem);

  if (!pProblem) fatalError();

  CScanMethod * pMethod = dynamic_cast<CScanMethod *>(mpMethod);

  if (!pMethod) fatalError();

  bool success = true;

  //initSubtask();

  if (useInitialValues)
    {
      mpProblem->getModel()->applyInitialValues();
    }

  //TODO: reports

  //initialize the method (parsing the ScanItems)
  pMethod->setProblem(pProblem);

  if (!pMethod->init()) return false;

  // init progress bar
  mProgress = 0;

  if (mpCallBack)
    {
      mpCallBack->setName("performing parameter scan...");

      unsigned C_INT32 totalSteps = pMethod->getTotalNumberOfSteps();
      mhProgress = mpCallBack->addItem("Number of Steps",
                                       CCopasiParameter::UINT,
                                       &mProgress,
                                       &totalSteps);

      if (mpSubtask)
        mpSubtask->setCallBack(mpCallBack);
    }

  // init output handler (plotting)
  output(COutputInterface::BEFORE);

  //calling the scanner, output is done in the callback
  if (!pMethod->scan()) success = false;

  //finishing progress bar and output
  //if (mpCallBack) mpCallBack->finish();
  //if (mpOutputHandler) mpOutputHandler->finish();
  output(COutputInterface::AFTER);

  if (mpSubtask)
    mpSubtask->setCallBack(NULL);

  return success;
}

bool CScanTask::processCallback()
{
  bool success = mpSubtask->process(!mAdjustInitialConditions);

  //do output
  if (success && !mOutputInSubtask)
    output(COutputInterface::DURING);

  //do progress bar
  ++mProgress;

  if (mpCallBack) return mpCallBack->progressItem(mhProgress);

  return true;
}

bool CScanTask::outputSeparatorCallback(bool isLast)
{
  if ((!isLast) || mOutputInSubtask)
    separate(COutputInterface::DURING);

  return true;
}

bool CScanTask::initSubtask(const OutputFlag & /* of */,
                            COutputHandler * pOutputHandler,
                            std::ostream * pOstream)
{
  if (!mpProblem) fatalError();

  CScanProblem * pProblem = dynamic_cast<CScanProblem *>(mpProblem);

  if (!pProblem) fatalError();

  //get the parameters from the problem
  CCopasiTask::Type type = *(CCopasiTask::Type*) pProblem->getValue("Subtask").pUINT;
  CCopasiDataModel* pDataModel = getObjectDataModel();
  assert(pDataModel != NULL);

  switch (type)
    {
      case CCopasiTask::steadyState:
        mpSubtask = dynamic_cast<CCopasiTask*>
                    ((*pDataModel->getTaskList())["Steady-State"]);
        break;

      case CCopasiTask::timeCourse:
        mpSubtask = dynamic_cast<CCopasiTask*>
                    ((*pDataModel->getTaskList())["Time-Course"]);
        break;

      case CCopasiTask::mca:
        mpSubtask = dynamic_cast<CCopasiTask*>
                    ((*pDataModel->getTaskList())["Metabolic Control Analysis"]);
        break;

      case CCopasiTask::lyap:
        mpSubtask = dynamic_cast<CCopasiTask*>
                    ((*pDataModel->getTaskList())["Lyapunov Exponents"]);
        break;

      case CCopasiTask::optimization:
        mpSubtask = dynamic_cast<CCopasiTask*>
                    ((*pDataModel->getTaskList())["Optimization"]);
        break;

      case CCopasiTask::parameterFitting:
        mpSubtask = dynamic_cast<CCopasiTask*>
                    ((*pDataModel->getTaskList())["Parameter Estimation"]);
        break;

      case CCopasiTask::sens:
        mpSubtask = dynamic_cast<CCopasiTask*>
                    ((*pDataModel->getTaskList())["Sensitivities"]);
        break;

      default:
        mpSubtask = NULL;
    }

  mOutputInSubtask = * pProblem->getValue("Output in subtask").pBOOL;
  //if (type != CCopasiTask::timeCourse)
  //  mOutputInSubtask = false;

  mAdjustInitialConditions = * pProblem->getValue("Adjust initial conditions").pBOOL;

  if (!mpSubtask) return false;

  mpSubtask->getProblem()->setModel(pDataModel->getModel()); //TODO
  mpSubtask->setCallBack(NULL);

  if (mOutputInSubtask)
    return mpSubtask->initialize(OUTPUT, pOutputHandler, pOstream);
  else
    return mpSubtask->initialize(NO_OUTPUT, pOutputHandler, pOstream);

  return true;
}
