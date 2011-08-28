// Begin CVS Header
//   $Source: /fs/turing/cvs/copasi_dev/copasi/trajectory/CTrajectoryMethod.cpp,v $
//   $Revision: 1.45 $
//   $Name: Build-33 $
//   $Author: shoops $
//   $Date: 2010/09/13 15:06:38 $
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
 *  CTrajectoryMethod class.
 *  This class describes the interface to all integration methods.
 *  The various method like LSODA or Gillespie have to be derived from
 *  this class.
 *
 *  Created for COPASI by Stefan Hoops 2002
 */

#include "copasi.h"

#include "CTrajectoryMethod.h"
#include "CLsodaMethod.h"
#include "CStochDirectMethod.h"
#include "CStochNextReactionMethod.h"
#include "CTrajAdaptiveSA.h"
#include "CHybridMethod.h"
#include "CHybridMethodLSODA.h"
#include "CTauLeapMethod.h"
#include "CHybridMethodLSODA.h"
#include "CTrajectoryMethodDsaLsodar.h"
#include "CTrajectoryProblem.h"
#include "model/CState.h"
#include "model/CCompartment.h"

CTrajectoryMethod *
CTrajectoryMethod::createTrajectoryMethod(CCopasiMethod::SubType subType,
    CTrajectoryProblem * /* pProblem */)
{
  CTrajectoryMethod * pMethod = NULL;

  switch (subType)
    {
      case unset:
      case deterministic:
        pMethod = new CLsodaMethod();
        break;

      case stochastic:
        pMethod = new CStochNextReactionMethod();
        break;

      case directMethod:
        pMethod = new CStochDirectMethod();
        break;

      case tauLeap:
        pMethod = new CTauLeapMethod();
        break;

      case adaptiveSA:
        pMethod = new CTrajAdaptiveSA();
        break;

      case hybrid:
        pMethod = CHybridMethod::createHybridMethod();
        break;

      case hybridLSODA:
        pMethod = CHybridMethodLSODA::createHybridMethodLSODA();
        break;

      case DsaLsodar:
        pMethod = new CTrajectoryMethodDsaLsodar();
        break;

      default:
        fatalError();
    }

  return pMethod;
}

/**
 *  Default constructor.
 */
CTrajectoryMethod::CTrajectoryMethod(const CCopasiMethod::SubType & subType,
                                     const CCopasiContainer * pParent) :
    CCopasiMethod(CCopasiTask::timeCourse, subType, pParent),
    mpCurrentState(NULL),
    mpProblem(NULL),
    mRoots(0)
{CONSTRUCTOR_TRACE;}

/**
 *  Copy constructor.
 *  @param "const CTrajectoryMethod &" src
 */
CTrajectoryMethod::CTrajectoryMethod(const CTrajectoryMethod & src,
                                     const CCopasiContainer * pParent):
    CCopasiMethod(src, pParent),
    mpCurrentState(src.mpCurrentState),
    mpProblem(src.mpProblem),
    mRoots(src.mRoots)
{CONSTRUCTOR_TRACE;}

/**
 *  Destructor.
 */
CTrajectoryMethod::~CTrajectoryMethod()
{DESTRUCTOR_TRACE;}

void CTrajectoryMethod::setCurrentState(CState * currentState)
{
  mpCurrentState = currentState;
}

/**
 *  Set a pointer to the problem.
 *  This method is used by CTrajectory
 *  @param "CTrajectoryProblem *" problem
 */
void CTrajectoryMethod::setProblem(CTrajectoryProblem * problem)
{mpProblem = problem;}

// virtual
void CTrajectoryMethod::stateChanged()
{}

/**
 *  This instructs the method to calculate a a time step of deltaT
 *  starting with the current state, i.e., the result of the previous
 *  step.
 *  The new state (after deltaT) is expected in the current state.
 *  The return value is the actual timestep taken.
 *  @param "const double &" deltaT
 */
CTrajectoryMethod::Status CTrajectoryMethod::step(const double & C_UNUSED(deltaT))
{return FAILURE;}

/**
 *  This instructs the method to calculate a a time step of deltaT
 *  starting with the initialState given.
 *  The new state (after deltaT) is expected in the current state.
 *  The return value is the actual timestep taken.
 *  @param "double &" deltaT
 *  @param "const CState *" initialState
 *  @return "const double &" actualDeltaT
 */
void CTrajectoryMethod::start(const CState * C_UNUSED(initialState))
{return;}

//virtual
bool CTrajectoryMethod::isValidProblem(const CCopasiProblem * pProblem)
{
  if (!CCopasiMethod::isValidProblem(pProblem)) return false;

  const CTrajectoryProblem * pTP = dynamic_cast<const CTrajectoryProblem *>(pProblem);

  if (!pTP)
    {
      //not a TrajectoryProblem
      CCopasiMessage(CCopasiMessage::EXCEPTION, MCTrajectoryMethod + 8);
      return false;
    }

  return true;
}

const CVector< C_INT > & CTrajectoryMethod::getRoots() const
{
  return mRoots;
}
