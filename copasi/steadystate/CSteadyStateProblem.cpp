/* Begin CVS Header
 $Source: /fs/turing/cvs/copasi_dev/copasi/steadystate/CSteadyStateProblem.cpp,v $
 $Revision: 1.28 $
 $Name: Build-33 $
 $Author: shoops $
 $Date: 2010/02/11 19:42:49 $
 End CVS Header */

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
 *  CSteadyStateProblem class.
 *  This class describes the steady state problem, i.e., it allows to specify
 *  for example initial conditions.
 *
 *  Created for Copasi by Stefan Hoops 2002
 */

#include <string>

#include "copasi.h"
#include "CSteadyStateProblem.h"

#include "CopasiDataModel/CCopasiDataModel.h"
#include "report/CCopasiRootContainer.h"
#include "model/CModel.h"

/**
 *  Default constructor.
 */
CSteadyStateProblem::CSteadyStateProblem(const CCopasiContainer * pParent):
    CCopasiProblem(CCopasiTask::steadyState, pParent)
    //mInitialState(),
    //mHasInitialState(false)
{
  addParameter("JacobianRequested", CCopasiParameter::BOOL, true);
  addParameter("StabilityAnalysisRequested", CCopasiParameter::BOOL, true);
  CONSTRUCTOR_TRACE;
}

/**
 *  Copy constructor.
 *  @param "const CSteadyStateProblem &" src
 */
CSteadyStateProblem::CSteadyStateProblem(const CSteadyStateProblem & src,
    const CCopasiContainer * pParent):
    CCopasiProblem(src, pParent)
    //mInitialState(src.mInitialState),
    //mHasInitialState(src.mHasInitialState)
{CONSTRUCTOR_TRACE;}

/**
 *  Destructor.
 */
CSteadyStateProblem::~CSteadyStateProblem()
{DESTRUCTOR_TRACE;}

bool CSteadyStateProblem::initialize()
{
  if (!mpModel) return false;

  return true;
}

/**
 * Set whether the jacobian is requested.
 * @param bool * jacobianRequested
 */
void CSteadyStateProblem::setJacobianRequested(bool & jacobianRequested)
{setValue("JacobianRequested", jacobianRequested);}

/**
 * Retrieve whether the jacobian is requested.
 * @return bool jacobianRequested
 */
bool CSteadyStateProblem::isJacobianRequested() const
{return * getValue("JacobianRequested").pBOOL;}

/**
 * Set whether stabilty analysis is requested.
 * @param bool * stabilityAnalysisRequested
 */
void CSteadyStateProblem::setStabilityAnalysisRequested(bool & stabilityAnalysisRequested)
{setValue("StabilityAnalysisRequested", stabilityAnalysisRequested);}

/**
 * Retrieve whether the stabilty analysis is requested.
 * @return bool stabilityAnalysisRequested
 */
bool CSteadyStateProblem::isStabilityAnalysisRequested() const
{return * getValue("StabilityAnalysisRequested").pBOOL;}

/**
 * Load a steadystate problem
 * @param "CReadConfig &" configBuffer
 */
void CSteadyStateProblem::load(CReadConfig & configBuffer,
                               CReadConfig::Mode C_UNUSED(mode))
{
  if (configBuffer.getVersion() < "4.0")
    {
      CCopasiDataModel* pDataModel = getObjectDataModel();
      assert(pDataModel != NULL);
      mpModel = pDataModel->getModel();
      //mInitialState = mpModel->getInitialState();
      //mHasInitialState = false;
      configBuffer.getVariable("RepStabilityAnalysis", "bool" ,
                               getValue("StabilityAnalysisRequested").pBOOL,
                               CReadConfig::LOOP);
      setValue("JacobianRequested",
               * getValue("StabilityAnalysisRequested").pBOOL);
    }
}
