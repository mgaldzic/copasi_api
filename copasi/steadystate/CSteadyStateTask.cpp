// Begin CVS Header
//   $Source: /fs/turing/cvs/copasi_dev/copasi/steadystate/CSteadyStateTask.cpp,v $
//   $Revision: 1.84 $
//   $Name: Build-33 $
//   $Author: shoops $
//   $Date: 2010/07/16 19:03:27 $
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
 * CSteadyStateTask class.
 *
 * This class implements a steady state task which is comprised of a
 * of a problem and a method. Additionally calls to the reporting
 * methods are done when initialized.
 *
 * Created for COPASI by Stefan Hoops 2002
 */

#include "copasi.h"

#include "CSteadyStateTask.h"
#include "CSteadyStateProblem.h"
#include "CSteadyStateMethod.h"
#include "model/CModel.h"
#include "model/CState.h"
#include "model/CMetabNameInterface.h"
#include "report/CKeyFactory.h"
#include "report/CReport.h"

#include <sstream>

#define XXXX_Reporting

CSteadyStateTask::CSteadyStateTask(const CCopasiContainer * pParent):
    CCopasiTask(CCopasiTask::steadyState, pParent),
    mpSteadyState(NULL),
    mJacobian(),
    mJacobianX(),
    mpJacobianAnn(NULL),
    mpJacobianXAnn(NULL),
    mEigenValues("Eigenvalues of Jacobian", this),
    mEigenValuesX("Eigenvalues of reduced system Jacobian", this)
{
  mpProblem = new CSteadyStateProblem(this);
  mpMethod =
    CSteadyStateMethod::createSteadyStateMethod(CCopasiMethod::Newton);
  this->add(mpMethod, true);
  //mpMethod->setObjectParent(this);
  //((CSteadyStateMethod *) mpMethod)->setProblem((CSteadyStateProblem *) mpProblem);
  initObjects();
}

CSteadyStateTask::CSteadyStateTask(const CSteadyStateTask & src,
                                   const CCopasiContainer * pParent):
    CCopasiTask(src, pParent),
    mpSteadyState(src.mpSteadyState),
    mJacobian(src.mJacobian),
    mJacobianX(src.mJacobianX),
    mpJacobianAnn(NULL),
    mpJacobianXAnn(NULL),
    mEigenValues(src.mEigenValues, this),
    mEigenValuesX(src.mEigenValuesX, this)
{
  mpProblem =
    new CSteadyStateProblem(*(CSteadyStateProblem *) src.mpProblem, this);
  mpMethod =
    CSteadyStateMethod::createSteadyStateMethod(src.mpMethod->getSubType());
  this->add(mpMethod, true);
  //mpMethod->setObjectParent(this);
  //((CSteadyStateMethod *) mpMethod)->setProblem((CSteadyStateProblem *) mpProblem);
  initObjects();
}

CSteadyStateTask::~CSteadyStateTask()
{
  pdelete(mpSteadyState);
}

void CSteadyStateTask::cleanup()
{}

void CSteadyStateTask::initObjects()
{
  mpJacobianAnn = new CArrayAnnotation("Jacobian (complete system)", this,
                                       new CCopasiMatrixInterface<CMatrix<C_FLOAT64> >(&mJacobian), true);
  mpJacobianAnn->setMode(CArrayAnnotation::OBJECTS);
  mpJacobianAnn->setDescription("");
  mpJacobianAnn->setDimensionDescription(0, "Variables of the system, including dependent species");
  mpJacobianAnn->setDimensionDescription(1, "Variables of the system, including dependent species");

  mpJacobianXAnn = new CArrayAnnotation("Jacobian (reduced system)", this,
                                        new CCopasiMatrixInterface<CMatrix<C_FLOAT64> >(&mJacobianX), true);
  mpJacobianXAnn->setMode(CArrayAnnotation::OBJECTS);
  mpJacobianXAnn->setDescription("");
  mpJacobianXAnn->setDimensionDescription(0, "Independent variables of the system");
  mpJacobianXAnn->setDimensionDescription(1, "Independent variables of the system");

  mpEigenvaluesJacobianAnn = new CArrayAnnotation("Eigenvalues of Jacobian", this,
      new CCopasiMatrixInterface<CMatrix<C_FLOAT64> >(&mEigenvaluesMatrix), true);
  mpEigenvaluesJacobianAnn->setMode(CArrayAnnotation::VECTOR);
  mpEigenvaluesJacobianAnn->setDescription("");
  mpEigenvaluesJacobianAnn->setDimensionDescription(0, "n-th value");
  mpEigenvaluesJacobianAnn->setDimensionDescription(1, "Real/Imaginary part");

  mpEigenvaluesJacobianXAnn = new CArrayAnnotation("Eigenvalues of reduced system Jacobian", this,
      new CCopasiMatrixInterface<CMatrix<C_FLOAT64> >(&mEigenvaluesXMatrix), true);
//  mpEigenvaluesJacobianXAnn->setMode(CArrayAnnotation::VECTOR);
  mpEigenvaluesJacobianXAnn->setMode(CArrayAnnotation::OBJECTS);
  mpEigenvaluesJacobianXAnn->setDescription("");
  mpEigenvaluesJacobianXAnn->setDimensionDescription(0, "n-th value");
  mpEigenvaluesJacobianXAnn->setDimensionDescription(1, "Real/Imaginary part");
}

void CSteadyStateTask::print(std::ostream * ostream) const {(*ostream) << (*this);}

void CSteadyStateTask::load(CReadConfig & configBuffer)
{
  configBuffer.getVariable("SteadyState", "bool", &mScheduled,
                           CReadConfig::LOOP);

  ((CSteadyStateProblem *) mpProblem)->load(configBuffer);

  ((CSteadyStateMethod *) mpMethod)->load(configBuffer);
}

const CState * CSteadyStateTask::getState() const
{return mpSteadyState;}

const CMatrix< C_FLOAT64 > & CSteadyStateTask::getJacobian() const
{return mJacobian;}
const CMatrix< C_FLOAT64 > & CSteadyStateTask::getJacobianReduced() const
{return mJacobianX;}

const CArrayAnnotation * CSteadyStateTask::getJacobianAnnotated() const
{
  return mpJacobianAnn;
}

const CArrayAnnotation * CSteadyStateTask::getJacobianXAnnotated() const
{
  return mpJacobianXAnn;
}

const CEigen & CSteadyStateTask::getEigenValues() const
{
  return mEigenValues;
}
const CEigen & CSteadyStateTask::getEigenValuesReduced() const
{
  return mEigenValuesX;
}

bool CSteadyStateTask::updateMatrices()
{
  if (!mpProblem->getModel()) return false;

  const CStateTemplate & stateTemplate = mpProblem->getModel()->getStateTemplate();

  // init Jacobians
  unsigned C_INT32 sizeX = stateTemplate.getNumIndependent();
  mJacobianX.resize(sizeX, sizeX);
  unsigned C_INT32 size = sizeX + stateTemplate.getNumDependent();
  mJacobian.resize(size, size);

  // Jacobian Annotations

  mpJacobianAnn->resize();
  CModelEntity *const* ppEntities = stateTemplate.getEntities();
  const unsigned C_INT32 * pUserOrder = stateTemplate.getUserOrder().array();
  const unsigned C_INT32 * pUserOrderEnd = pUserOrder + stateTemplate.getUserOrder().size();

  pUserOrder++; // We skip the time which is the first.

  unsigned C_INT32 i, imax = size;

  for (i = 0; i < imax && pUserOrder != pUserOrderEnd; pUserOrder++)
    {
      const CModelEntity::Status & Status = ppEntities[*pUserOrder]->getStatus();

      if (Status == CModelEntity::ODE ||
          (Status == CModelEntity::REACTIONS && ppEntities[*pUserOrder]->isUsed()))
        {
          mpJacobianAnn->setAnnotationCN(0 , i, ppEntities[*pUserOrder]->getCN());
          mpJacobianAnn->setAnnotationCN(1 , i, ppEntities[*pUserOrder]->getCN());

          i++;
        }
    }

  mpJacobianXAnn->resize();

  ppEntities = stateTemplate.beginIndependent();
  imax = sizeX;

  for (i = 0; i < imax; ++i, ++ppEntities)
    {
      mpJacobianXAnn->setAnnotationCN(0 , i, (*ppEntities)->getCN());
      mpJacobianXAnn->setAnnotationCN(1 , i, (*ppEntities)->getCN());
    }

  // initial dimension of Eigenvalues of Jacobian
  mEigenvaluesMatrix.resize(size, 2);
  mEigenvaluesXMatrix.resize(sizeX, 2);

  return true;
}

bool CSteadyStateTask::initialize(const OutputFlag & of,
                                  COutputHandler * pOutputHandler,
                                  std::ostream * pOstream)
{
  assert(mpProblem && mpMethod);

  if (!mpMethod->isValidProblem(mpProblem)) return false;

  bool success = true;

  if (!updateMatrices())
    return false;

  success &= CCopasiTask::initialize(of, pOutputHandler, pOstream);

  pdelete(mpSteadyState);
  mpSteadyState = new CState(mpProblem->getModel()->getInitialState());

  mCalculateReducedSystem = (mpProblem->getModel()->getNumDependentReactionMetabs() != 0);

#ifdef xxxx
  // init Jacobians
  unsigned C_INT32 sizeX = mpSteadyState->getNumIndependent();
  mJacobianX.resize(sizeX, sizeX);
  unsigned C_INT32 size = sizeX + mpSteadyState->getNumDependent();
  mJacobian.resize(size, size);

  //jacobian annotations
  CStateTemplate & StateTemplate = mpProblem->getModel()->getStateTemplate();

  mpJacobianAnn->resize();
  CModelEntity **ppEntities = StateTemplate.getEntities();
  const unsigned C_INT32 * pUserOrder = StateTemplate.getUserOrder().array();
  const unsigned C_INT32 * pUserOrderEnd = pUserOrder + StateTemplate.getUserOrder().size();

  pUserOrder++; // We skip the time which is the first.

  unsigned C_INT32 i, imax = size;

  for (i = 0; i < imax && pUserOrder != pUserOrderEnd; pUserOrder++)
    {
      const CModelEntity::Status & Status = ppEntities[*pUserOrder]->getStatus();

      if (Status == CModelEntity::ODE ||
          (Status == CModelEntity::REACTIONS && ppEntities[*pUserOrder]->isUsed()))
        {
          mpJacobianAnn->setAnnotationCN(0 , i, ppEntities[*pUserOrder]->getCN());
          mpJacobianAnn->setAnnotationCN(1 , i, ppEntities[*pUserOrder]->getCN());

          i++;
        }
    }

  mpJacobianXAnn->resize();

  ppEntities = StateTemplate.beginIndependent();
  imax = sizeX;

  for (i = 0; i < imax; ++i, ++ppEntities)
    {
      mpJacobianXAnn->setAnnotationCN(0 , i, (*ppEntities)->getCN());
      mpJacobianXAnn->setAnnotationCN(1 , i, (*ppEntities)->getCN());
    }

#endif

  CSteadyStateProblem* pProblem =
    dynamic_cast<CSteadyStateProblem *>(mpProblem);
  assert(pProblem);

  success &= pProblem->initialize();

  CSteadyStateMethod* pMethod =
    dynamic_cast<CSteadyStateMethod *>(mpMethod);
  assert(pMethod);

  success &= pMethod->initialize(pProblem);

  return success;
}

bool CSteadyStateTask::process(const bool & useInitialValues)
{
  if (mpInitialState != NULL)
    {
      *mpInitialState = mpProblem->getModel()->getInitialState();
    }

  if (useInitialValues)
    {
      mpProblem->getModel()->applyInitialValues();
    }

  *mpSteadyState = mpProblem->getModel()->getState();

  // A steady-state makes only sense in an autonomous model,
  // i.e., the time of the steady-state must not be changed
  // during simulation.
  C_FLOAT64 InitialTime = mpSteadyState->getTime();

  CSteadyStateMethod* pMethod =
    dynamic_cast<CSteadyStateMethod *>(mpMethod);
  assert(pMethod);

  CSteadyStateProblem* pProblem =
    dynamic_cast<CSteadyStateProblem *>(mpProblem);
  assert(pMethod);

  output(COutputInterface::BEFORE);

  //call the method
  mResult = pMethod->process(mpSteadyState,
                             mJacobianX,
                             mpCallBack);

  if (mResult == CSteadyStateMethod::notFound)
    restore();

  //update jacobian
  if (pProblem->isJacobianRequested() ||
      pProblem->isStabilityAnalysisRequested())
    {
      pMethod->doJacobian(mJacobian, mJacobianX);
    }

  //mpProblem->getModel()->setState(mpSteadyState);
  //mpProblem->getModel()->updateRates();

  //calculate eigenvalues
  if (pProblem->isStabilityAnalysisRequested())
    {
      mEigenValues.calcEigenValues(mJacobian);
      mEigenValuesX.calcEigenValues(mJacobianX);

      mEigenValues.stabilityAnalysis(pMethod->getStabilityResolution());
      mEigenValuesX.stabilityAnalysis(pMethod->getStabilityResolution());
    }

  // Reset the time.
  mpSteadyState->setTime(InitialTime);

  C_FLOAT64 * pTo;
  size_t i;

  // construct Eigenvalues of Jacobian
  CVector< C_FLOAT64 > vectorEigen_R = mEigenValues.getR();
  CVector< C_FLOAT64 > vectorEigen_I = mEigenValues.getI();

#ifdef DEBUG_UI
  C_INT32 size = vectorEigen_R.size() + vectorEigen_I.size();

  std::cout << "vectorEigen_R.size() = " << vectorEigen_R.size() << " + vectorEigen_I.size() = " << vectorEigen_I.size() << " == " << size << std::endl;
  std::cout << "size = " << mEigenvaluesXMatrix.size() << std::endl;
#endif
  assert(vectorEigen_R.size() == vectorEigen_I.size());

  pTo = mEigenvaluesMatrix.array();

  for (i = 0; i < vectorEigen_R.size(); ++i)
    {
      *pTo = vectorEigen_R[i]; ++pTo;
      *pTo = vectorEigen_I[i]; ++pTo;
    }

#ifdef DEBUG_UI
  std::cout << mEigenvaluesMatrix << std::endl;
#endif

  // construct Eigenvalues of Jacobian of reduced system
  CVector< C_FLOAT64 > vectorEigenX_R = mEigenValuesX.getR();
  CVector< C_FLOAT64 > vectorEigenX_I = mEigenValuesX.getI();

#ifdef DEBUG_UI
  C_INT32 sizeX = vectorEigenX_R.size() + vectorEigenX_I.size();

  std::cout << "vectorEigenX_R.size() = " << vectorEigenX_R.size() << " + vectorEigenX_I.size() = " << vectorEigenX_I.size() << " == " << sizeX << std::endl;
  std::cout << "size = " << mEigenvaluesXMatrix.size() << std::endl;
#endif

  assert(vectorEigenX_R.size() == vectorEigenX_I.size());

  pTo = mEigenvaluesXMatrix.array();

  for (i = 0; i < vectorEigenX_R.size(); ++i)
    {
      *pTo = vectorEigenX_R[i]; ++pTo;
      *pTo = vectorEigenX_I[i]; ++pTo;
    }

#ifdef DEBUG_UI
  std::cout << mEigenvaluesXMatrix << std::endl;
#endif

  output(COutputInterface::AFTER);

  return (mResult != CSteadyStateMethod::notFound);
}

void CSteadyStateTask::setInitialState()
{
  CModel * pModel = mpProblem->getModel();
  pModel->setInitialState(pModel->getState());
}

bool CSteadyStateTask::restore()
{
  bool success = CCopasiTask::restore();

  if (mUpdateModel)
    {
      CModel * pModel = mpProblem->getModel();

      pModel->setState(*mpSteadyState);
      pModel->updateSimulatedValues(true);
      pModel->setInitialState(pModel->getState());
      pModel->updateInitialValues();
    }

  return success;
}

std::ostream &operator<<(std::ostream &os, const CSteadyStateTask &A)
{
  switch (A.getResult())
    {
      case CSteadyStateMethod::found:
        os << "A steady state with given resolution was found." << std::endl;
        break;

      case CSteadyStateMethod::notFound:
        os << "No steady state with given resolution was found!" << std::endl;
        os << "(below are the last unsuccessful trial values)" << std::endl;
        break;

      case CSteadyStateMethod::foundEquilibrium:
        os << "An equilibrium steady state (zero fluxes) was found." << std::endl;
        break;

      case CSteadyStateMethod::foundNegative:
        os << "An invalid steady state (negative concentrations) was found." << std::endl;
    }

  os << std::endl;

  // Update all necessary values.
  CState * pState = const_cast<CState *>(A.getState());

  if (!pState) return os;

  CModel * pModel = A.mpProblem->getModel();

  if (!pModel) return os;

  pModel->setState(*pState);
  pModel->updateSimulatedValues(true);
  pModel->updateNonSimulatedValues();

  // Metabolite Info: Name, Concentration, Concentration Rate, Particle Number, Particle Rate, Transition Time
  const CCopasiVector<CMetab> & Metabolites = pModel->getMetabolites();
  const CMetab * pMetab;

  unsigned C_INT32 i, imax = Metabolites.size();

  os << "Species" << "\t";
  os << "Concentration";

  std::string Units = pModel->getConcentrationUnitsDisplayString();

  if (Units != "")
    os << " (" << Units << ")";

  os << "\t";

  os << "Concentration Rate";
  Units = pModel->getConcentrationRateUnitsDisplayString();

  if (Units != "")
    os << " (" << Units << ")";

  os << "\t";

  os << "Particle Number" << "\t";

  os << "Particle Number Rate";
  Units = pModel->getFrequencyUnitsDisplayString();

  if (Units != "")
    os << " (" << Units << ")";

  os << "\t";

  os << "Transition Time";
  Units = pModel->getTimeUnitsDisplayString();

  if (Units != "")
    os << " (" << Units << ")";

  os << std::endl;

  for (i = 0; i < imax; ++i)
    {
      pMetab = Metabolites[i];
      os << CMetabNameInterface::getDisplayName(pModel, *pMetab) << "\t";
      os << pMetab->getConcentration() << "\t";
      os << pMetab->getConcentrationRate() << "\t";
      os << pMetab->getValue() << "\t";
      os << pMetab->getRate() << "\t";
      os << pMetab->getTransitionTime() << std::endl;
    }

  os << std::endl;

  // Reaction Info: Name, Flux, Particle Flux
  const CCopasiVector<CReaction>& Reactions = pModel->getReactions();
  const CReaction * pReaction;

  imax = Reactions.size();

  os << "Reaction" << "\t";

  os << "Flux";
  Units = pModel->getQuantityRateUnitsDisplayString();

  if (Units != "")
    os << " (" << Units << ")";

  os << "\t";

  os << "Particle Flux";
  Units = pModel->getFrequencyUnitsDisplayString();;

  if (Units != "")
    os << " (" << Units << ")";

  os << std::endl;

  for (i = 0; i < imax; ++i)
    {
      pReaction = Reactions[i];
      os << pReaction->getObjectName() << "\t";
      os << pReaction->getFlux() << "\t";
      os << pReaction->getParticleFlux() << std::endl;
    }

  os << std::endl;

  if (static_cast<CSteadyStateProblem *>(A.mpProblem)->isJacobianRequested())
    {
      os << *A.mpJacobianAnn << std::endl;

      if (static_cast<CSteadyStateProblem *>(A.mpProblem)->isStabilityAnalysisRequested())
        {
          os << "Eigenvalues\treal\timaginary" << std::endl;
          imax = A.mEigenValues.getR().size();

          for (i = 0; i < imax; i++)
            os << "\t" << A.mEigenValues.getR()[i] << "\t" << A.mEigenValues.getI()[i] << std::endl;

          os << std::endl;
        }

      os << *A.mpJacobianXAnn << std::endl;

      if (static_cast<CSteadyStateProblem *>(A.mpProblem)->isStabilityAnalysisRequested())
        {
          os << "Eigenvalues\treal\timaginary" << std::endl;
          imax = A.mEigenValuesX.getR().size();

          for (i = 0; i < imax; i++)
            os << "\t" << A.mEigenValuesX.getR()[i] << "\t" << A.mEigenValuesX.getI()[i] << std::endl;

          os << std::endl;
        }
    }

  if (static_cast<CSteadyStateProblem *>(A.mpProblem)->isStabilityAnalysisRequested())
    {
      os << "Stability Analysis of the Reduced System" << std::endl;
      os << A.mEigenValuesX << std::endl;
    }

  return os;
}
