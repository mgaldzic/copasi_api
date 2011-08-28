// Begin CVS Header
//   $Source: /fs/turing/cvs/copasi_dev/copasi/steadystate/CSteadyStateMethod.cpp,v $
//   $Revision: 1.33 $
//   $Name: Build-33 $
//   $Author: shoops $
//   $Date: 2010/07/16 19:03:26 $
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
 *  CSteadyStateMethod class.
 *  This class describes the interface to all steady state methods.
 *  The various method like Newton have to be derived from
 *  this class.
 *
 *  Created for COPASI by Stefan Hoops 2002
 */

#include "copasi.h"

#include "utilities/CCopasiVector.h"
#include "CSteadyStateMethod.h"
#include "CSteadyStateProblem.h"
#include "CSteadyStateTask.h"
#include "CEigen.h"

#include "model/CModel.h"
#include "model/CState.h"
#include "model/CCompartment.h"

CSteadyStateMethod *
CSteadyStateMethod::createSteadyStateMethod(CCopasiMethod::SubType subType)
{
  CSteadyStateMethod * pMethod = NULL;

  switch (subType)
    {
      case unset:
      case Newton:
        pMethod = new CNewtonMethod();
        break;

      default:
        fatalError();
    }

  return pMethod;
}

/**
 *  Default constructor.
 */
CSteadyStateMethod::CSteadyStateMethod(CCopasiMethod::SubType subType,
                                       const CCopasiContainer * pParent):
    CCopasiMethod(CCopasiTask::steadyState, subType, pParent),
    mpProblem(NULL)
{
  initializeParameter();
  CONSTRUCTOR_TRACE;
}

/**
 *  Copy constructor.
 *  @param "const CSteadyStateMethod &" src
 */
CSteadyStateMethod::CSteadyStateMethod(const CSteadyStateMethod & src,
                                       const CCopasiContainer * pParent):
    CCopasiMethod(src, pParent),
    mpProblem(src.mpProblem)
{
  initializeParameter();
  CONSTRUCTOR_TRACE;
}

/**
 *  Destructor.
 */
CSteadyStateMethod::~CSteadyStateMethod()
{DESTRUCTOR_TRACE;}

void CSteadyStateMethod::initializeParameter()
{
  CCopasiParameter *pParm;

  assertParameter("Resolution", CCopasiParameter::UDOUBLE, (C_FLOAT64) 1.0e-009);
  mpSSResolution = (C_FLOAT64*)getValue("Resolution").pUDOUBLE;
  mpDerivationResolution = (C_FLOAT64*)getValue("Resolution").pUDOUBLE;

  assertParameter("Derivation Factor", CCopasiParameter::UDOUBLE, (C_FLOAT64) 1.0e-003);
  mpDerivationFactor = (C_FLOAT64*)getValue("Derivation Factor").pUDOUBLE;

  // Check whether we have a method with the old parameter names
  if ((pParm = getParameter("Newton.DerivationFactor")) != NULL)
    {
      setValue("Derivation Factor", *pParm->getValue().pUDOUBLE);
      removeParameter("Newton.DerivationFactor");
    }

  if ((pParm = getParameter("Newton.Resolution")) != NULL)
    {
      setValue("Resolution", *pParm->getValue().pUDOUBLE);
      removeParameter("Newton.Resolution");
    }
}

bool CSteadyStateMethod::elevateChildren()
{
  initializeParameter();
  return true;
}

/**
 * This instructs the method to calculate a the steady state
 * starting with the initialState given.
 * The steady state is returned in the object pointed to by steadyState.
 * @param CState * steadyState
 * @param const CState * initialState
 * @param C_FLOAT64 * jacobian
 * @param CEigen * eigenValues
 * @return CSteadyStateMethod::ReturnCode returnCode
 */
CSteadyStateMethod::ReturnCode
CSteadyStateMethod::process(CState * pState,
                            CMatrix< C_FLOAT64 > & jacobianX,
                            CProcessReport * handler)
{
  mpParentTask = dynamic_cast<CSteadyStateTask *>(getObjectParent());
  assert(mpParentTask);

  mpSteadyState = pState;
  mpJacobianX = & jacobianX;
  mpProgressHandler = handler;

  return processInternal();
}

/**
 * This function has to be called at the return of any implementation
 * of the protected function process
 * @param bool success
 * @param const C_FLOAT64 & factor
 * @param const C_FLOAT64 & resolution
 * @return CSteadyStateMethod::ReturnCode returnCode
 */
CSteadyStateMethod::ReturnCode
CSteadyStateMethod::returnProcess(bool steadyStateFound)
{
  if (!steadyStateFound)
    return CSteadyStateMethod::notFound;

  if (!allPositive())
    return CSteadyStateMethod::foundNegative;

  if (isEquilibrium(*mpSSResolution))
    return CSteadyStateMethod::foundEquilibrium;

  return CSteadyStateMethod::found;
}

/**
 * This instructs the method to calculate a the steady state
 * @return CSteadyStateMethod::ReturnCode returnCode
 */
CSteadyStateMethod::ReturnCode
CSteadyStateMethod::processInternal()
{return CSteadyStateMethod::notFound;}

bool CSteadyStateMethod::isEquilibrium(const C_FLOAT64 & resolution) const
{
  const CCopasiVectorNS < CReaction > & Reaction =
    mpProblem->getModel()->getReactions();
  unsigned C_INT32 i, imax = Reaction.size();

  for (i = 0; i < imax; i++)
    if (Reaction[i]->getFlux() / Reaction[i]->getLargestCompartment().getValue() > resolution)
      return false; //TODO: smallest or largest ?

  return true;
}

bool CSteadyStateMethod::allPositive()
{
  // Assure that all values are updated.
  mpModel->updateSimulatedValues(true);

  CModelEntity *const* ppEntity = mpModel->getStateTemplate().beginIndependent();
  const C_FLOAT64 * pIt = mpModel->getState().beginIndependent();
  const C_FLOAT64 * pEnd = mpModel->getState().endDependent();

  // Skip Model quantities of type ODE
  for (; pIt != pEnd; ++pIt, ++ppEntity)
    if (dynamic_cast< const CCompartment *>(*ppEntity) != NULL ||
        dynamic_cast< const CMetab *>(*ppEntity) != NULL)
      break;

  // For all compartments of type ODE we check that the volume is positive
  for (; pIt != pEnd; ++pIt, ++ppEntity)
    {
      if (dynamic_cast< const CCompartment *>(*ppEntity) == NULL)
        break;

      if (*pIt < - *mpDerivationResolution)
        return false;
    }

  // We need to check that all metabolites have positive particle numbers
  // with respect to the given resolution.
  C_FLOAT64 ParticleResolution =
    - *mpDerivationResolution * mpModel->getQuantity2NumberFactor();

  for (; pIt != pEnd; ++pIt, ++ppEntity)
    {
      if (dynamic_cast< const CMetab *>(*ppEntity) == NULL)
        break;

      if (*pIt < ParticleResolution)
        return false;
    }

  // For all compartments of type ASSIGNMENT we check that the volume is positive
  for (; pIt != pEnd; ++pIt, ++ppEntity)
    {
      if (dynamic_cast< CCompartment *>(*ppEntity) == NULL)
        break;

      if (*pIt < - *mpDerivationResolution)
        return false;
    }

  return true;
}

//virtual
bool CSteadyStateMethod::isValidProblem(const CCopasiProblem * pProblem)
{
  if (!CCopasiMethod::isValidProblem(pProblem)) return false;

  const CSteadyStateProblem * pP = dynamic_cast<const CSteadyStateProblem *>(pProblem);

  if (!pP)
    {
      //not a TrajectoryProblem
      CCopasiMessage(CCopasiMessage::EXCEPTION, "Problem is not a steady state problem.");
      return false;
    }

  return true;
}

bool CSteadyStateMethod::initialize(const CSteadyStateProblem * pProblem)
{
  mpProblem = pProblem;
  mpModel = mpProblem->getModel();

  return true;
}

void CSteadyStateMethod::doJacobian(CMatrix< C_FLOAT64 > & jacobian,
                                    CMatrix< C_FLOAT64 > & jacobianX)
{
  mpModel->setState(*mpSteadyState);
  mpModel->updateSimulatedValues(true);

  mpModel->calculateJacobian(jacobian, *mpDerivationResolution, *mpDerivationResolution);
  mpModel->calculateJacobianX(jacobianX, *mpDerivationResolution, *mpDerivationResolution);
}

C_FLOAT64 CSteadyStateMethod::getStabilityResolution()
{
  //C_FLOAT64* pTmp = (C_FLOAT64*)getValue("Stability Resolution").pUDOUBLE;
  C_FLOAT64* pTmp = (C_FLOAT64*)getValue("Resolution").pUDOUBLE;
  assert(pTmp);
  return *pTmp;
}

void CSteadyStateMethod::calculateJacobianX(const C_FLOAT64 & oldMaxRate)
{
  mpModel->setState(*mpSteadyState);
  mpModel->updateSimulatedValues(true);
  mpModel->calculateJacobianX(*mpJacobianX,
                              std::min(*mpDerivationFactor, oldMaxRate),
                              *mpDerivationResolution);
}

std::string CSteadyStateMethod::getMethodLog() const
{
  return mMethodLog.str();
}
