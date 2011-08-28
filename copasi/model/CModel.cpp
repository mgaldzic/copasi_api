// Begin CVS Header
//   $Source: /fs/turing/cvs/copasi_dev/copasi/model/CModel.cpp,v $
//   $Revision: 1.395 $
//   $Name: Build-33 $
//   $Author: shoops $
//   $Date: 2010/09/21 16:48:02 $
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
//

#include <string.h>
#include <limits.h>
#include <string>
#include <vector>
#include <limits>
#include <float.h>
#include <cmath>
#include <algorithm>

#ifdef SunOS
# include <ieeefp.h>
#else
# include <float.h>
#endif

#include "copasi.h"

#include "CCompartment.h"
#include "CMetab.h"
#include "CModel.h"
#include "CState.h"
#include "CModelValue.h"
#include "function/CFunctionDB.h"
#include "report/CCopasiObjectReference.h"
#include "report/CKeyFactory.h"
#include "utilities/CCopasiException.h"
#include "utilities/CCopasiMessage.h"
#include "utilities/CCopasiVector.h"
#include "CopasiDataModel/CCopasiDataModel.h"
#include "report/CCopasiRootContainer.h"
#include "utilities/CVector.h"
#include "utilities/CluX.h"
#include "utilities/utility.h"
#include "utilities/CProcessReport.h"
#include "CReactionInterface.h"
#include "utilities/CAnnotatedMatrix.h"
#include "CMetabNameInterface.h"
#include "CMathModel.h"

#include "blaswrap.h"
#include "clapackwrap.h"

#ifdef COPASI_DEBUG
#define CCHECK {check();}
#else
#define CCHECK
#endif

#define MNumMetabolitesReactionDependent (mNumMetabolitesReaction - mNumMetabolitesReactionIndependent)

const char * CModel::VolumeUnitNames[] =
  {"dimensionless", "m\xc2\xb3", "l", "ml", "\xc2\xb5l", "nl", "pl", "fl", NULL};

const char * CModel::AreaUnitNames[] =
  {"dimensionless", "m\xc2\xb2", "dm\xc2\xb2", "cm\xc2\xb2", "mm\xc2\xb2", "\xc2\xb5m\xc2\xb2", "nm\xc2\xb2", "pm\xc2\xb2", "fm\xc2\xb2", NULL};

const char * CModel::LengthUnitNames[] =
  {"dimensionless", "m", "dm", "cm", "mm", "\xc2\xb5m", "nm", "pm", "fm", NULL};

const char * CModel::TimeUnitNames[] =
  {"dimensionless", "d", "h", "min", "s", "ms", "\xc2\xb5s", "ns", "ps", "fs", NULL};

// "mol" is the correct name, however in the COPASI XML files "Mol" is used
// up to build 18

const char * CModel::QuantityUnitOldXMLNames[] =
  {"dimensionless", "Mol", "mMol", "\xc2\xb5Mol", "nMol", "pMol", "fMol", "#", NULL};

const char * CModel::QuantityUnitNames[] =
  {"dimensionless", "mol", "mmol", "\xc2\xb5mol", "nmol", "pmol", "fmol", "#", NULL};

const char * CModel::ModelTypeNames[] =
  {"deterministic", "stochastic", NULL};

CModel::CModel(CCopasiContainer* pParent):
    CModelEntity("New Model", pParent, "Model"),
    mInitialState(),
    mCurrentState(),
    mStateTemplate(*this, this->mInitialState, this->mCurrentState),
    mVolumeUnit(ml),
    mAreaUnit(m2),
    mLengthUnit(m),
    mTimeUnit(s),
    mQuantityUnit(mMol),
    mType(deterministic),
    mCompartments("Compartments", this),
    mMetabolites("Metabolites", this),
    mMetabolitesX("Reduced Model Metabolites", this),
    mSteps("Reactions", this),
    mEvents("Events", this),
    mParticleFluxes(),
    mValues("Values", this),
    mMoieties("Moieties", this),
    mStoi(),
    mpStoiAnnotation(NULL),
    mStoiReordered(),
    mRedStoi(),
    mpRedStoiAnnotation(NULL),
    mNumMetabolitesUnused(0),
    mNumMetabolitesODE(0),
    mNumMetabolitesReaction(0),
    mNumMetabolitesAssignment(0),
    mNumMetabolitesReactionIndependent(0),
    mL(),
    mpLinkMatrixAnnotation(NULL),
    mLView(mL, mNumMetabolitesReactionIndependent),
    mAvogadro(6.02214179e23),
    mQuantity2NumberFactor(1.0),
    mNumber2QuantityFactor(1.0),
    mpCompileHandler(NULL),
    mInitialRefreshes(),
    mSimulatedRefreshes(),
    mApplyInitialValuesRefreshes(),
    mNonSimulatedRefreshes(),
    mReorderNeeded(false),
    mIsAutonomous(true),
    mBuildInitialSequence(true),
    mpMathModel(NULL)
{
  initObjects();

  setStatus(TIME);
  setUsed(true);

  *mpIValue = 0.0;
  *mpValue = std::numeric_limits<C_FLOAT64>::quiet_NaN();

  unsigned C_INT32 i, imax = mSteps.size();

  for (i = 0; i < imax; i++)
    mSteps[i]->compile(/*mCompartments*/);

  initializeMetabolites();

  forceCompile(NULL);

  /* This following 2 lines added by Liang Xu
  Because of the failure to initialize the parameter when creating a new models
  */
  setQuantityUnit(mQuantityUnit); // set the factors
  setVolumeUnit(mVolumeUnit); // set the factors

  CONSTRUCTOR_TRACE;
}

// CModel::CModel(const CModel & src):
//     CModelEntity(src),
//     mInitialState(),
//     mCurrentState(),
//     mStateTemplate(*this, this->mInitialState, this->mCurrentState),
//     mComments(src.mComments),
//     mVolumeUnit(src.mVolumeUnit),
//     mTimeUnit(src.mTimeUnit),
//     mQuantityUnit(src.mQuantityUnit),
//     mType(src.mType),
//     mCompartments(src.mCompartments, this),
//     mMetabolites(src.mMetabolites, this),
//     mMetabolitesX(src.mMetabolitesX, this),
//     mSteps(src.mSteps, this),
//     mEvents(src.mEvents,this),
//     mParticleFluxes(src.mParticleFluxes),
//     mValues(src.mValues, this),
//     mMoieties(src.mMoieties, this),
//     mStoi(src.mStoi),
//     mpStoiAnnotation(NULL),
//     mStoiReordered(src.mStoiReordered),
//     mRedStoi(src.mRedStoi),
//     mpRedStoiAnnotation(NULL),
//     mNumMetabolitesUnused(src.mNumMetabolitesUnused),
//     mNumMetabolitesODE(src.mNumMetabolitesODE),
//     mNumMetabolitesReaction(src.mNumMetabolitesReaction),
//     mNumMetabolitesAssignment(src.mNumMetabolitesAssignment),
//     mNumMetabolitesIndependent(src.mNumMetabolitesIndependent),
//     mL(src.mL),
//     mpLinkMatrixAnnotation(NULL),
//     mLView(mL, mNumMetabolitesIndependent),
//     mQuantity2NumberFactor(src.mQuantity2NumberFactor),
//     mNumber2QuantityFactor(src.mNumber2QuantityFactor),
//     mpCompileHandler(NULL),
//     mInitialRefreshes(),
//     mSimulatedRefreshes(),
//     mConstantRefreshes(),
//     mNonSimulatedRefreshes(),
//     mReorderNeeded(false),
//     mIsAutonomous(false)
// {
//   CONSTRUCTOR_TRACE;
//   initObjects();
//
//   unsigned C_INT32 i, imax = mSteps.size();
//
//   for (i = 0; i < imax; i++)
//     mSteps[i]->compile(/*mCompartments*/);
//
//   initializeMetabolites();
//
//   forceCompile(NULL);
//}

CModel::~CModel()
{
  mpIValue = NULL;
  mpValue = NULL;

  pdelete(mpStoiAnnotation);
  pdelete(mpRedStoiAnnotation);
  pdelete(mpLinkMatrixAnnotation);

  CCopasiRootContainer::getKeyFactory()->remove(mKey);

  DESTRUCTOR_TRACE;
}

// virtual
std::string CModel::getChildObjectUnits(const CCopasiObject * pObject) const
{
  if (pObject->getObjectName() == "Initial Time" ||
      pObject->getObjectName() == "Time")
    return getTimeUnitName();

  return "";
}

void CModel::cleanup()
{
  /* The real objects */
  mCompartments.cleanup();
  mSteps.cleanup();
  mMoieties.cleanup();

  /* The references */
  //mStepsInd.resize(0);
  mMetabolites.clear();
  mMetabolitesX.clear();
  mParticleFluxes.resize(0);
}

C_INT32 CModel::load(CReadConfig & configBuffer)
{
  C_INT32 Size = 0;
  C_INT32 Fail = 0;
  unsigned C_INT32 i;
  std::string tmp;

  // For old Versions we must read the list of Metabolites beforehand
  if ((Fail = configBuffer.getVariable("TotalMetabolites", "C_INT32",
                                       &Size, CReadConfig::LOOP)))
    return Fail;

  // :TODO: Remove OldMetabolites as part of the data model.
  CCopasiDataModel* pDataModel = getObjectDataModel();
  assert(pDataModel != NULL);
  pDataModel->pOldMetabolites->load(configBuffer, Size);

  if ((Fail = configBuffer.getVariable("Title", "string", &tmp,
                                       CReadConfig::LOOP)))
    return Fail;

  setObjectName(tmp);

  std::string Notes;

  try
    {
      Fail = configBuffer.getVariable("Comments", "multiline", &Notes,
                                      CReadConfig::SEARCH);
    }

  catch (CCopasiException Exception)
    {
      if ((MCReadConfig + 1) == Exception.getMessage().getNumber())
        Notes = "";
      else
        throw Exception;
    }

  setNotes(Notes);

  try
    {
      Fail = configBuffer.getVariable("TimeUnit", "string", &tmp,
                                      CReadConfig::LOOP);
    }
  catch (CCopasiException Exception)
    {
      if ((MCReadConfig + 1) == Exception.getMessage().getNumber())
        tmp = ""; //unknown?
      else
        throw Exception;
    }

  setTimeUnit(tmp); // set the factors

  try
    {
      Fail = configBuffer.getVariable("ConcentrationUnit", "string", &tmp,
                                      CReadConfig::LOOP);
    }
  catch (CCopasiException Exception)
    {
      if ((MCReadConfig + 1) == Exception.getMessage().getNumber())
        tmp = "";  //unknown?
      else
        throw Exception;
    }

  setQuantityUnit(tmp); // set the factors

  try
    {
      Fail = configBuffer.getVariable("VolumeUnit", "string", &tmp,
                                      CReadConfig::LOOP);
    }
  catch (CCopasiException Exception)
    {
      if ((MCReadConfig + 1) == Exception.getMessage().getNumber())
        tmp = ""; //unknown?
      else
        throw Exception;
    }

  setVolumeUnit(tmp); // set the factors

  *mpIValue = 0;

  if ((Fail = configBuffer.getVariable("TotalCompartments", "C_INT32", &Size,
                                       CReadConfig::LOOP)))
    return Fail;

  mCompartments.load(configBuffer, Size);

  // Create the correct compartment / metabolite relationships
  CMetab *pMetabolite;

  for (i = 0; i < pDataModel->pOldMetabolites->size(); i++)
    {
      pMetabolite = new CMetab;
      mCompartments[(*pDataModel->pOldMetabolites)[i]->getIndex()]->
      addMetabolite(pMetabolite);

      (*pMetabolite) = *(*pDataModel->pOldMetabolites)[i];
    }

  initializeMetabolites();

  if ((Fail = CCopasiRootContainer::getFunctionList()->load(configBuffer))) // slow
    return Fail;

  if ((Fail = configBuffer.getVariable("TotalSteps", "C_INT32", &Size,
                                       CReadConfig::LOOP)))
    return Fail;

  mSteps.load(configBuffer, Size); // slow

  for (i = 0; i < mSteps.size(); i++)
    mSteps[i]->compile(/*mCompartments*/);

  pDataModel->pOldMetabolites->cleanup();

  setCompileFlag();
  return Fail;
}

bool CModel::compile()
{
  mpValueReference->addDirectDependency(this);

  CMatrix< C_FLOAT64 > LU;

  unsigned C_INT32 CompileStep = 0;
  unsigned C_INT32 hCompileStep;

  if (mpCompileHandler)
    {
      mpCompileHandler->setName("Compiling model...");
      unsigned C_INT32 totalSteps = 7;
      hCompileStep = mpCompileHandler->addItem("Compile Process",
                     CCopasiParameter::UINT,
                     & CompileStep,
                     &totalSteps);
    }

  // To assure that we do not produce access violations we clear the refresh sequences
  // first
  mInitialRefreshes.clear();
  mSimulatedRefreshes.clear();
  mApplyInitialValuesRefreshes.clear();
  mNonSimulatedRefreshes.clear();

  CompileStep = 0;

  if (mpCompileHandler && !mpCompileHandler->progressItem(hCompileStep)) return false;

  buildStoi();
  CompileStep = 1;
  
  if (mpCompileHandler && !mpCompileHandler->progressItem(hCompileStep)) return false;
 
  buildLinkZero();
  CompileStep = 2;
  
  if (mpCompileHandler && !mpCompileHandler->progressItem(hCompileStep)) return false;
  
  buildRedStoi();
  CompileStep = 3;

  if (mpCompileHandler && !mpCompileHandler->progressItem(hCompileStep)) return false;

  buildMoieties();
  CompileStep = 4;

  if (mpCompileHandler && !mpCompileHandler->progressItem(hCompileStep)) return false;

  buildStateTemplate();
  CompileStep = 5;

  if (mpCompileHandler && !mpCompileHandler->progressItem(hCompileStep)) return false;

  bool success = true;

  try
    {
      success &= buildInitialSequence();
      success &= buildApplyInitialValuesSequence();
      success &= buildSimulatedSequence();
      success &= buildNonSimulatedSequence();
    }
  catch (...)
    {
      success = false;
    }

  CompileStep = 6;
  
  if (mpCompileHandler && !mpCompileHandler->progressItem(hCompileStep)) return false;

  buildUserOrder();

  if (mpCompileHandler) mpCompileHandler->finishItem(hCompileStep);

  //update annotations
  updateMatrixAnnotations();

  success &= compileEvents();
  success &= mpMathModel->compile(this);

  if (!success)
    {
      mIsAutonomous = false;
      return false;
    }

  mCompileIsNecessary = false;
  determineIsAutonomous();

  //writeDependenciesToDotFile();

  return true;
}

void CModel::compileDefaultMetabInitialValueDependencies()
{
  CCopasiVector< CMetab >::iterator it = mMetabolites.begin();
  CCopasiVector< CMetab >::iterator end = mMetabolites.end();

  for (; it != end; ++it)
  {
    (*it)->compileInitialValueDependencies(true);
   }

  mBuildInitialSequence = true;

  return;
}

void CModel::setCompileFlag(bool flag)
{
  mCompileIsNecessary = flag;
}

bool CModel::compileIfNecessary(CProcessReport* pProcessReport)
{
  if (!mCompileIsNecessary)
    {
      compileDefaultMetabInitialValueDependencies();

      return true;
    }

  mpCompileHandler = pProcessReport;
  
  try
    {
      compile();
    }

  catch (...)
    {
    }
    
  mpCompileHandler = NULL;

  return true;
}

bool CModel::forceCompile(CProcessReport* pProcessReport)
{
  setCompileFlag();
  return compileIfNecessary(pProcessReport);
}

void CModel::buildStoi()
{
  unsigned C_INT32 i;

  initializeMetabolites();

  unsigned C_INT32 numRows, numCols;
  numRows = mNumMetabolitesReaction;
  numCols = mSteps.size();

  mParticleFluxes.resize(numCols);
  mStoi.resize(numRows, numCols);
  mStoi = 0.0;

  unsigned C_INT32 hProcess;

  if (mpCompileHandler)
    {
      i = 0;
      hProcess = mpCompileHandler->addItem("Building Stoichiometry",
                                           CCopasiParameter::UINT,
                                           &i,
                                           &numCols);
    }

  C_FLOAT64 * pCol, *pColEnd;
  pCol = mStoi.array();
  pColEnd = mStoi.array() + numCols;

  C_FLOAT64 * pRow, *pRowEnd;
  pRowEnd = mStoi.array() + numRows * numCols;

  CCopasiVector< CReaction >::iterator itStep = mSteps.begin();
  CCopasiVector< CMetab >::const_iterator itMetab;

  for (; pCol < pColEnd; ++pCol, ++itStep, ++i)
    {
      if (mpCompileHandler && !mpCompileHandler->progressItem(hProcess)) return;

      // Since we are stepping through the reactions we can check whether
      // the kinetic functions are usable.
      if (!(*itStep)->getFunction()->isUsable())
        CCopasiMessage(CCopasiMessage::ERROR, MCReaction + 11,
                       (*itStep)->getObjectName().c_str(),
                       (*itStep)->getFunction()->getObjectName().c_str());

      const CCopasiVector< CChemEqElement > & Balance =
        (*itStep)->getChemEq().getBalances();
      CCopasiVector< CChemEqElement >::const_iterator itBalance = Balance.begin();
      CCopasiVector< CChemEqElement >::const_iterator endBalance = Balance.end();

      for (; itBalance != endBalance; ++itBalance)
        {
          const std::string & key = (*itBalance)->getMetaboliteKey();

          for (pRow = pCol, itMetab = mMetabolitesX.begin() + mNumMetabolitesODE;
               pRow < pRowEnd;
               pRow += numCols, ++itMetab)
            if ((*itMetab)->getKey() == key)
              {
                *pRow = (*itBalance)->getMultiplicity();
                break;
              }
        }
    }

  // We need to have all unused and fixed metabolites at the end of mMetabolites.
  // However we can only detect unused metabolites after building the
  // stoichiometry matrix.
  handleUnusedMetabolites();

#ifdef DEBUG_MATRIX
  DebugFile << "Stoichiometry Matrix" << std::endl;
  DebugFile << mStoi << std::endl;
#endif

  if (mpCompileHandler)
    mpCompileHandler->finishItem(hProcess);

  return;
}

bool CModel::handleUnusedMetabolites()
{
  unsigned C_INT32 numRows, numCols;
  numRows = mStoi.numRows();
  numCols = mStoi.numCols();

  C_FLOAT64 * pStoi, *pStoiEnd, *pRowEnd;
  pStoi = mStoi.array();
  pStoiEnd = mStoi.array() + numRows * numCols;

  unsigned C_INT32 i, NumUnused;
  C_FLOAT64 tmp;
  std::vector< unsigned C_INT32 > Unused;

  for (i = 0; i < numRows; i++)
    {
      tmp = 0;

      for (pRowEnd = pStoi + numCols; pStoi < pRowEnd; ++pStoi)
        tmp += fabs(*pStoi);

      if (tmp < DBL_MIN) Unused.push_back(i);
    }

  NumUnused = Unused.size();

  if (NumUnused == 0) return false;

  // We treat unused variables in the same way as fixed, i.e.
  // they will be sorted to the end of the metabolite list.

  numRows -= NumUnused;

  CMatrix< C_FLOAT64 > NewStoi(numRows, numCols);
  C_FLOAT64 * pNewStoi = NewStoi.array();
  std::vector< CMetab * > UsedMetabolites(numRows);
  std::vector< CMetab * >::iterator itUsedMetabolites = UsedMetabolites.begin();
  std::vector< CMetab * > UnusedMetabolites(NumUnused);
  std::vector< CMetab * >::iterator itUnusedMetabolites = UnusedMetabolites.begin();
  std::vector< unsigned C_INT32 >::const_iterator itUnused = Unused.begin();
  std::vector< unsigned C_INT32 >::const_iterator endUnused = Unused.end();

  CCopasiVector< CMetab >::iterator itMetab = mMetabolitesX.begin() + mNumMetabolitesODE;
  CCopasiVector< CMetab >::iterator endMetab = itMetab + mNumMetabolitesReaction;

  // Build new stoichiometry Matrix
  pStoi = mStoi.array();

  for (i = 0; itMetab != endMetab; ++itMetab, i++, pStoi += numCols)
    {
      if (itUnused != endUnused && i == *itUnused)
        {
          (*itMetab)->setUsed(false);
          *itUnusedMetabolites = (*itMetab);

          ++itUnusedMetabolites;
          ++itUnused;
        }
      else
        {
          (*itMetab)->setUsed(true);
          *itUsedMetabolites = (*itMetab);
          ++itUsedMetabolites;

          // The row needs to be copied to the new stoichiometry matrix
          memcpy(pNewStoi, pStoi, sizeof(C_FLOAT64) * numCols);
          pNewStoi += numCols;
        }
    }

  // Reorder metabolites
  // Skip the metabolites determined by ODE
  itMetab = mMetabolitesX.begin() + mNumMetabolitesODE;

  // Handle the metabolites determined by actually determined by reactions
  itUsedMetabolites = UsedMetabolites.begin();
  std::vector< CMetab * >::iterator itMetabolitesEnd = UsedMetabolites.end();

  for (; itUsedMetabolites != itMetabolitesEnd; ++itUsedMetabolites, ++itMetab)
    *itMetab = *itUsedMetabolites;

  // Handle metabolites determined by assignment and marked as fixed
  // This is just a shift of NumUnused.
  endMetab = itMetab + mNumMetabolitesAssignment + mNumMetabolitesUnused;

  for (; itMetab != endMetab; ++itMetab)
    *itMetab = *(itMetab + NumUnused);

  // Handle newly marked unused metabolites
  itUnusedMetabolites = UnusedMetabolites.begin();
  itMetabolitesEnd = UnusedMetabolites.end();

  for (; itUnusedMetabolites != itMetabolitesEnd; ++itUnusedMetabolites, ++itMetab)
    *itMetab = *itUnusedMetabolites;

  // Now its time to update the number of metabolites determined by reactions
  // and the one being unused.
  mNumMetabolitesReaction -= NumUnused;
  mNumMetabolitesUnused += NumUnused;

  // Update stoichiometry matrix
  mStoi = NewStoi;

  return true;
}

void CModel::buildRedStoi()
{
  unsigned C_INT32 i;
  unsigned C_INT32 numCols = mStoi.numCols();

  mRedStoi.resize(mNumMetabolitesReactionIndependent, numCols);
  mStoiReordered.resize(mStoi.numRows(), numCols);

  C_FLOAT64 * pRedStoi = mRedStoi.array();
  C_FLOAT64 * pStoiReordered = mStoiReordered.array();
  unsigned C_INT32 * pRow = mRowLU.array();

  // Create a temporary copy of the metabolites determined by reactions to reorder them
  // accordingly.
  CCopasiVector< CMetab >::iterator itMetabX = mMetabolitesX.begin() + mNumMetabolitesODE;
  CCopasiVector< CMetab >::iterator endMetabX = itMetabX + mNumMetabolitesReaction;
  std::vector< CMetab * > MetabolitesReaction;
  MetabolitesReaction.resize(mNumMetabolitesReaction);
  std::vector< CMetab * >::iterator itMetabolitesReaction = MetabolitesReaction.begin();

  for (; itMetabX != endMetabX; ++itMetabX, ++itMetabolitesReaction)
    *itMetabolitesReaction = *itMetabX;

  /* just have to swap rows */
  itMetabX = mMetabolitesX.begin() + mNumMetabolitesODE;

  for (i = 0; i < mNumMetabolitesReactionIndependent; i++, pRow++, pRedStoi += numCols, pStoiReordered += numCols, itMetabX++)
    {
      memcpy(pRedStoi, mStoi[*pRow], sizeof(C_FLOAT64) * numCols);
      memcpy(pStoiReordered, mStoi[*pRow], sizeof(C_FLOAT64) * numCols);
      *itMetabX = MetabolitesReaction[*pRow];
    }

  for (; i < mNumMetabolitesReaction; i++, pRow++, pStoiReordered += numCols, itMetabX++)
    {
      memcpy(pStoiReordered, mStoi[*pRow], sizeof(C_FLOAT64) * numCols);
      *itMetabX = MetabolitesReaction[*pRow];
    }

#ifdef DEBUG_MATRIX
  DebugFile << "Reduced Stoichiometry Matrix" << std::endl;
  DebugFile << mRedStoi << std::endl;
#endif

  return;
}

void CModel::updateMatrixAnnotations()
{
  mpLinkMatrixAnnotation->resize();
  mpStoiAnnotation->resize();
  mpRedStoiAnnotation->resize();

  CCopasiVector< CMetab >::const_iterator it = mMetabolitesX.begin() + mNumMetabolitesODE;
  CCopasiVector< CMetab >::const_iterator end = it + mNumMetabolitesReactionIndependent;

  CCopasiObjectName CN;
  unsigned C_INT32 j;

  for (j = 0; it != end; ++it, j++)
    {
      CN = (*it)->getCN();

      mpStoiAnnotation->setAnnotationCN(0, j, CN);
      mpLinkMatrixAnnotation->setAnnotationCN(0, j, CN);
      mpLinkMatrixAnnotation->setAnnotationCN(1, j, CN);
      mpRedStoiAnnotation->setAnnotationCN(0, j, CN);
    }

  end += MNumMetabolitesReactionDependent;

  for (; it != end; ++it, j++)
    {
      CN = (*it)->getCN();

      mpStoiAnnotation->setAnnotationCN(0, j, CN);
      mpLinkMatrixAnnotation->setAnnotationCN(0, j, CN);
    }
}

void CModel::updateMoietyValues()
{
  CCopasiVector< CMoiety >::iterator it = mMoieties.begin();
  CCopasiVector< CMoiety >::iterator end = mMoieties.end();

  for (; it != end; ++it)
    (*it)->refreshInitialValue();
}

void CModel::buildMoieties()
{
  unsigned C_INT32 i, imax = MNumMetabolitesReactionDependent;
  unsigned C_INT32 j;

  CCopasiVector< CMetab >::iterator it =
    mMetabolitesX.begin() + mNumMetabolitesODE + mNumMetabolitesReactionIndependent; //begin of dependent metabs
  C_FLOAT64 * pFactor = mL.array();

  CMoiety *pMoiety;

  mMoieties.cleanup();

  for (i = 0; i < imax; i++, ++it)
    {
      pMoiety = new CMoiety((*it)->getObjectName());
      pMoiety->add(1.0, *it);

      for (j = 0; j < mNumMetabolitesReactionIndependent; j++, pFactor++)
        if (fabs(*pFactor) > std::numeric_limits< C_FLOAT64 >::epsilon())
          pMoiety->add(- *pFactor, mMetabolitesX[j + mNumMetabolitesODE]);

      mMoieties.add(pMoiety, true);
    }

  updateMoietyValues();
  return;
}

//this is supposed to be so fast it can be called often to be kept up to date
//all the time. At the moment it creates the mMetabolites and sorts the fixed
//metabs to the end
void CModel::initializeMetabolites()
{
  // Create a vector of pointers to all metabolites.
  // Note, the metabolites physically exist in the compartments.
  mMetabolites.clear();

  CCopasiVector< CCompartment >::iterator itCompartment = mCompartments.begin();
  CCopasiVector< CCompartment >::iterator endCompartment = mCompartments.end();
  CCopasiVector< CMetab >::iterator itMetab;
  CCopasiVector< CMetab >::iterator endMetab;

  for (; itCompartment != endCompartment; ++itCompartment)
    {
      itMetab = (*itCompartment)->getMetabolites().begin();
      endMetab = (*itCompartment)->getMetabolites().end();

      for (; itMetab != endMetab; ++itMetab)
        {
          // Reset all moieties
          (*itMetab)->setDependentOn(NULL);
          mMetabolites.add(*itMetab);
        }
    }

  // We sort the metabolites by type
  itMetab = mMetabolites.begin();
  endMetab = mMetabolites.end();

  std::vector< CMetab *> ODEMetabs;
  std::vector< CMetab *> ReactionMetabs;
  std::vector< CMetab *> AssignmentMetabs;
  std::vector< CMetab *> FixedMetabs;

  for (; itMetab != endMetab; ++itMetab)
    switch ((*itMetab)->getStatus())
      {
        case FIXED:
          FixedMetabs.push_back(*itMetab);
          (*itMetab)->setUsed(false);
          break;

        case ASSIGNMENT:
          AssignmentMetabs.push_back(*itMetab);
          (*itMetab)->setUsed(true);
          break;

        case ODE:
          ODEMetabs.push_back(*itMetab);
          (*itMetab)->setUsed(true);
          break;

        case REACTIONS:
          ReactionMetabs.push_back(*itMetab);
          (*itMetab)->setUsed(true);
          break;

        default:
          fatalError();
      }

  // Update mMetabolitesX to reflect the reordering.
  // We need to to this to allow the use of the full model for simulation.
  mMetabolitesX.resize(mMetabolites.size());

  mNumMetabolitesODE = ODEMetabs.size();
  itMetab = mMetabolitesX.begin();
  std::vector< CMetab *>::const_iterator itSorted = ODEMetabs.begin();
  std::vector< CMetab *>::const_iterator endSorted = ODEMetabs.end();

  for (; itSorted != endSorted; ++itSorted, ++itMetab)
    *itMetab = *itSorted;

  mNumMetabolitesReaction = ReactionMetabs.size();
  itSorted = ReactionMetabs.begin();
  endSorted = ReactionMetabs.end();

  for (; itSorted != endSorted; ++itSorted, ++itMetab)
    *itMetab = *itSorted;

  mNumMetabolitesAssignment = AssignmentMetabs.size();
  itSorted = AssignmentMetabs.begin();
  endSorted = AssignmentMetabs.end();

  for (; itSorted != endSorted; ++itSorted, ++itMetab)
    *itMetab = *itSorted;

  mNumMetabolitesUnused = FixedMetabs.size();
  itSorted = FixedMetabs.begin();
  endSorted = FixedMetabs.end();

  for (; itSorted != endSorted; ++itSorted, ++itMetab)
    *itMetab = *itSorted;
}

//**********************************************************************

CCopasiVectorNS < CReaction > & CModel::getReactions()
{return mSteps;}

const CCopasiVectorNS < CReaction > & CModel::getReactions() const
{return mSteps;}

const CVector< C_FLOAT64 > & CModel::getParticleFlux() const
{return mParticleFluxes;}

CCopasiVector< CMetab > & CModel::getMetabolites()
{return mMetabolites;}

const CCopasiVector< CMetab > & CModel::getMetabolites() const
{return mMetabolites;}

CCopasiVector< CMetab > & CModel::getMetabolitesX()
{CCHECK return mMetabolitesX;}

const CCopasiVector< CMetab > & CModel::getMetabolitesX() const
{CCHECK return mMetabolitesX;}

const CCopasiVectorN< CModelValue > & CModel::getModelValues() const
{return mValues;}

CCopasiVectorN< CModelValue > & CModel::getModelValues()
{return mValues;}

CCopasiVectorN < CEvent > & CModel::getEvents()
{return mEvents;}

const CCopasiVectorN < CEvent > & CModel::getEvents() const
{return mEvents;}

//********

unsigned C_INT32 CModel::getNumMetabs() const
{return mMetabolites.size();}

unsigned C_INT32 CModel::getNumVariableMetabs() const
{return mNumMetabolitesODE + mNumMetabolitesReaction + mNumMetabolitesAssignment;}

unsigned C_INT32 CModel::getNumODEMetabs() const
{CCHECK return mNumMetabolitesODE;}

unsigned C_INT32 CModel::getNumAssignmentMetabs() const
{CCHECK return mNumMetabolitesAssignment;}

unsigned C_INT32 CModel::getNumIndependentReactionMetabs() const
{CCHECK return mNumMetabolitesReactionIndependent;}

unsigned C_INT32 CModel::getNumDependentReactionMetabs() const
{CCHECK return mNumMetabolitesReaction - mNumMetabolitesReactionIndependent;}

unsigned C_INT32 CModel::getTotSteps() const
{return mSteps.size();}

unsigned C_INT32 CModel::getNumModelValues() const
{return mValues.size();}

const std::string & CModel::getKey() const
{return mKey;}

CCopasiVectorNS < CCompartment > & CModel::getCompartments()
{return mCompartments;}

const CCopasiVectorNS < CCompartment > & CModel::getCompartments() const
{return mCompartments;}

/**
 *  Get the Reduced Stoichiometry Matrix of this Model
 */
const CMatrix < C_FLOAT64 >& CModel::getRedStoi() const
{CCHECK return mRedStoi;}

/**
 *  Get the Stoichiometry Matrix of this Model
 */
const CMatrix < C_FLOAT64 >& CModel::getStoi() const
{CCHECK return mStoi;}

/**
 *  Get the reordered stoichiometry matrix of this model
 */
const CMatrix < C_FLOAT64 >& CModel::getStoiReordered() const
{CCHECK return mStoiReordered;}

const CCopasiVector < CMoiety > & CModel::getMoieties() const
{return mMoieties;}

const CModel::CLinkMatrixView & CModel::getL() const
{CCHECK return mLView;}

const CMatrix< C_FLOAT64 > & CModel::getL0() const
{return mL;}

const CStateTemplate & CModel::getStateTemplate() const
{CCHECK return mStateTemplate;}

CStateTemplate & CModel::getStateTemplate()
{CCHECK return mStateTemplate;}

const std::set< const CCopasiObject * > & CModel::getUptoDateObjects() const
{return mSimulatedUpToDateObjects;}

bool CModel::setTitle(const std::string &title)
{
  if (title == "")
    return setObjectName("NoTitle");

  return setObjectName(title);
}

void CModel::setInitialTime(const C_FLOAT64 & time)
{*mpIValue = time;}

const C_FLOAT64 & CModel::getInitialTime() const
{return *mpIValue;}

void CModel::setTime(const C_FLOAT64 & time)
{*mpValue = time;}

const C_FLOAT64 & CModel::getTime() const
{return *mpValue;}

const CVector<unsigned C_INT32> & CModel::getMetabolitePermutation() const
{CCHECK return mRowLU;}

//**********************************************************************

/**
 *        Returns the index of the metab
 */
unsigned C_INT32 CModel::findMetabByName(const std::string & Target) const
{
  unsigned C_INT32 i, s;
  std::string name;

  s = mMetabolites.size();

  for (i = 0; i < s; i++)
    {
      name = mMetabolites[i]->getObjectName();

      if (name == Target)
        return i;
    }

  return C_INVALID_INDEX;
}

/**
 *        Returns the index of the Moiety
 */
unsigned C_INT32 CModel::findMoiety(const std::string &Target) const
{
  unsigned C_INT32 i, s;
  std::string name;

  s = mMoieties.size();

  for (i = 0; i < s; i++)
    {
      name = mMoieties[i]->getObjectName();

      if (name == Target)
        return i;
    }

  return C_INVALID_INDEX;
}

//**********************************************************************

void CModel::applyInitialValues()
{
  // Copy the initial state to the current state,
  setState(mInitialState);

  // Since the initial state is in itself consistent we should not need to
  // do anything further. However, for species of type ODE and ASSIGNMENT
  // the effective state variable is the concentration, i.e., we need to update
  // their concentration here.
  std::vector< Refresh * >::const_iterator itRefresh = mApplyInitialValuesRefreshes.begin();
  std::vector< Refresh * >::const_iterator endRefresh = mApplyInitialValuesRefreshes.end();

  while (itRefresh != endRefresh)
    (**itRefresh++)();

  // Update all dependent objects needed for simulation.
  updateSimulatedValues(false);

  mpMathModel->applyInitialValues();
}

void CModel::clearMoieties()
{
  mMoieties.clear();
}

bool CModel::buildStateTemplate()
{
  CVector<CModelEntity *> Entities(mMetabolitesX.size() + mCompartments.size() + mValues.size());
  CModelEntity ** ppEntity = Entities.array();
  
  CCopasiVector< CModelValue >::iterator itValue = mValues.begin();
  CCopasiVector< CModelValue >::iterator endValue = mValues.end();

  for (; itValue != endValue; ++itValue)
    if ((*itValue)->getStatus() == ODE)
      {
        *ppEntity = *itValue;
        (*ppEntity++)->setUsed(true);
      }

  CCopasiVector< CCompartment >::iterator itCompartment = mCompartments.begin();
  CCopasiVector< CCompartment >::iterator endCompartment = mCompartments.end();

  for (; itCompartment != endCompartment; ++itCompartment)
    if ((*itCompartment)->getStatus() == ODE)
      {
        *ppEntity = *itCompartment;
        (*ppEntity++)->setUsed(true);
      }

  CCopasiVector< CMetab >::iterator itMetab = mMetabolitesX.begin();
  CCopasiVector< CMetab >::iterator endMetab = mMetabolitesX.end();
  
  for (; itMetab != endMetab; ++itMetab)
  if (*itMetab)
  {
       if (!(*itMetab)->isUsed()) break;
      *ppEntity++ = *itMetab;
   }
  
  itCompartment = mCompartments.begin();

  for (; itCompartment != endCompartment; ++itCompartment)
    if ((*itCompartment)->getStatus() == ASSIGNMENT)
      {
        *ppEntity = *itCompartment;
        (*ppEntity++)->setUsed(true);
      }

  itValue = mValues.begin();

  for (; itValue != endValue; ++itValue)
    if ((*itValue)->getStatus() == ASSIGNMENT)
      {
        *ppEntity = *itValue;
        (*ppEntity++)->setUsed(true);
      }

  for (; itMetab != endMetab; ++itMetab)
    *ppEntity++ = *itMetab;

  itValue = mValues.begin();

  for (; itValue != endValue; ++itValue)
    if ((*itValue)->isFixed())
      *ppEntity++ = *itValue;

  itCompartment = mCompartments.begin();

  for (; itCompartment != endCompartment; ++itCompartment)
    if ((*itCompartment)->isFixed())
      *ppEntity++ = *itCompartment;

  mStateTemplate.reorder(Entities);
  mReorderNeeded = false;

  // Now all entities and reactions can be compiled
  ppEntity = mStateTemplate.beginIndependent();
  CModelEntity ** ppEntityEnd = mStateTemplate.endFixed();

  for (; ppEntity != ppEntityEnd; ++ppEntity)
    (*ppEntity)->compile();

  CCopasiVector< CReaction >::iterator itReaction = mSteps.begin();
  CCopasiVector< CReaction >::iterator endReaction = mSteps.end();

  for (; itReaction != endReaction; ++itReaction)
    (*itReaction)->compile();

  return true;
}

bool CModel::buildUserOrder()
{
  CVector<CModelEntity *> Entities(mMetabolites.size() + mCompartments.size() + mValues.size());
  CModelEntity ** ppEntity = Entities.array();

  CCopasiVector< CMetab >::iterator itMetab = mMetabolites.begin();
  CCopasiVector< CMetab >::iterator endMetab = mMetabolites.end();

  for (; itMetab != endMetab; ++itMetab)
    *ppEntity++ = *itMetab;;

  CCopasiVector< CCompartment >::iterator itCompartment = mCompartments.begin();
  CCopasiVector< CCompartment >::iterator endCompartment = mCompartments.end();

  for (; itCompartment != endCompartment; ++itCompartment)
    *ppEntity++ = *itCompartment;

  CCopasiVector< CModelValue >::iterator itValue = mValues.begin();
  CCopasiVector< CModelValue >::iterator endValue = mValues.end();

  for (; itValue != endValue; ++itValue)
    *ppEntity++ = *itValue;

  mStateTemplate.setUserOrder(Entities);

  mJacobianPivot.resize(mStateTemplate.getNumIndependent() + MNumMetabolitesReactionDependent);
  //now sized to the number of entities with ODEs + all metabolites dependent on reactions

  const unsigned C_INT32 * pUserOrder = mStateTemplate.getUserOrder().array();
  const unsigned C_INT32 * pUserOrderEnd = pUserOrder + mStateTemplate.getUserOrder().size();
  ppEntity = mStateTemplate.getEntities();

  unsigned C_INT32 i;

  for (i = 0; pUserOrder != pUserOrderEnd; ++pUserOrder)
    {
      const Status & Status = ppEntity[*pUserOrder]->getStatus();

      if (Status == ODE ||
          (Status == REACTIONS && ppEntity[*pUserOrder]->isUsed()))
        mJacobianPivot[i++] = *pUserOrder - 1;
    }

  return true;
}

bool CModel::buildInitialSequence()
{
  bool success = true;

  // The objects which are changed are all initial values of of all model entities including
  // fixed and unused once. Additionally, all kinetic parameters are possibly changed.
  // This is basically all the parameters in the parameter overview whose value is editable.

  // Issue 1170: We need to add elements of the stoichiometry, reduced stoichiometry,
  // and link matrix.

  std::set< const CCopasiObject * > Objects;

  // The initial values of the model entities
  CModelEntity **ppEntity = mStateTemplate.beginIndependent() - 1; // Offset for time
  CModelEntity **ppEntityEnd = mStateTemplate.endFixed();

  for (; ppEntity != ppEntityEnd; ++ppEntity)
    {
      // Assignments have no initial values
      if ((*ppEntity)->getStatus() != ASSIGNMENT ||
          (*ppEntity)->getInitialValueReference()->getDirectDependencies().size() == 0)
        Objects.insert((*ppEntity)->getInitialValueReference());
    }

  // The reaction parameters
  CCopasiVector< CReaction >::const_iterator itReaction = mSteps.begin();
  CCopasiVector< CReaction >::const_iterator endReaction = mSteps.end();
  unsigned C_INT32 i, imax;

  for (; itReaction != endReaction; ++itReaction)
    {
      const CCopasiParameterGroup & Group = (*itReaction)->getParameters();

      for (i = 0, imax = Group.size(); i < imax; i++)
        Objects.insert(Group.getParameter(i)->getObject(CCopasiObjectName("Reference=Value")));
    }

  // Fix for Issue 1170: We need to add elements of the stoichiometry, reduced stoichiometry,
  // and link matrices.
  if (mpStoiAnnotation != NULL)
    mpStoiAnnotation->appendElementReferences(Objects);

  if (mpRedStoiAnnotation != NULL)
    mpRedStoiAnnotation->appendElementReferences(Objects);

  if (mpLinkMatrixAnnotation != NULL)
    mpLinkMatrixAnnotation->appendElementReferences(Objects);

  try
    {
      mInitialRefreshes = buildInitialRefreshSequence(Objects);
    }
  catch (...)
    {
      mInitialRefreshes.clear();
      success = false;
    }

  mBuildInitialSequence = false;

  return success;
}

bool CModel::updateInitialValues()
{
  if (mCompileIsNecessary)
    {
      compileIfNecessary(NULL);
    }

  if (mBuildInitialSequence)
    {
      buildInitialSequence();
    }

  // Update all initial values.
  std::vector< Refresh * >::const_iterator itRefresh = mInitialRefreshes.begin();
  std::vector< Refresh * >::const_iterator endRefresh = mInitialRefreshes.end();

  while (itRefresh != endRefresh)
    (**itRefresh++)();

  return true;
}

bool CModel::buildSimulatedSequence()
{
  bool success = true;

  // We need to add each used model entity to the objects which need to be updated.
  mSimulatedUpToDateObjects.clear();

  // For CModelValues and CCompartment ODEs we need to add the Rate
  // For CMetab ODEs we need to add the Particle Rate
  CModelEntity **ppEntity = mStateTemplate.beginIndependent();
  CModelEntity **ppEntityEnd = mStateTemplate.endIndependent() - mNumMetabolitesReactionIndependent;

  for (; ppEntity != ppEntityEnd; ++ppEntity) //loop over all entities with ODEs
    mSimulatedUpToDateObjects.insert((*ppEntity)->getRateReference());

  // We do not add the rates for metabolites of type REACTION. These are automatically calculated
  // with dgemm in calculate derivatives based on the reaction fluxes added below.
  // In the case that other simulated values depend on such a rate this is taken care by
  // calculating all dependencies.
  // This mechanism may lead occasionally to multiple calculations of rates of metabolites when used
  // in assignments or ODEs. However this is acceptable and more than compensated by the performance
  // gains of dgemm.

  // Furthermore all reaction fluxes have to be calculated too (see CMetab REACTION above)
  CCopasiVector< CReaction >::iterator itReaction = mSteps.begin();
  CCopasiVector< CReaction >::iterator endReaction = mSteps.end();

  for (; itReaction != endReaction; ++itReaction)
    mSimulatedUpToDateObjects.insert((*itReaction)->getParticleFluxReference());

  // We now detect unused assignments, i.e., the result of an assignment is not
  // used during updateSimulatedValues except for itself or another unused assignment.
  bool UnusedFound = true;

  std::set<const CCopasiObject * > Candidate;
  std::set< const CCopasiObject * >::iterator it;
  std::set< const CCopasiObject * >::iterator end = mSimulatedUpToDateObjects.end();
  CCopasiObject * pObject;
  CMetab * pMetab;
  ppEntityEnd = mStateTemplate.endDependent();

  while (UnusedFound)
    {
      UnusedFound = false;
      ppEntity = mStateTemplate.beginDependent() + MNumMetabolitesReactionDependent;
      //now points to the first entity with assignment

      for (; ppEntity != ppEntityEnd; ++ppEntity) //over all entities with assignments
        if ((*ppEntity)->isUsed())
          {
            if ((*ppEntity)->getStatus() != ASSIGNMENT)
              pObject = *ppEntity;
            else if ((pMetab = dynamic_cast< CMetab *>(*ppEntity)) != NULL)
              pObject = pMetab->getConcentrationReference();
            else
              pObject = (*ppEntity)->getValueReference();

            Candidate.insert(pObject);

            for (it = mSimulatedUpToDateObjects.begin(), end = mSimulatedUpToDateObjects.end(); it != end; ++it)
              if (*it != pObject &&
                  (*it)->dependsOn(Candidate))
                break;

            if (it == end)
              {
                UnusedFound = true;
                mReorderNeeded = true;
                (*ppEntity)->setUsed(false);
              }

            Candidate.erase(pObject);
          }
    }

  if (mReorderNeeded)
    {
      CVector< CModelEntity * > Reorder(mStateTemplate.size() - 1);
      CModelEntity ** ppReorder = Reorder.array();

      ppEntity = mStateTemplate.beginIndependent();
      ppEntityEnd = mStateTemplate.beginDependent() + MNumMetabolitesReactionDependent;

      //range is for all entities with ODEs + all metabs dependent on Reactions
      for (; ppEntity != ppEntityEnd; ++ppEntity, ++ppReorder)
        *ppReorder = *ppEntity;

      // :TODO: This must be enhanced as the mMetaboliteX and the state template may get out of sync
      // when we use assignments for metabolites.

      // the entities with assignments are reordered according to the isUsed() flag
      ppEntityEnd = mStateTemplate.endDependent();

      for (; ppEntity != ppEntityEnd; ++ppEntity) //over all entities with assignments
        if ((*ppEntity)->isUsed())
          *ppReorder++ = *ppEntity;

      ppEntity = mStateTemplate.beginDependent() + MNumMetabolitesReactionDependent;

      for (; ppEntity != ppEntityEnd; ++ppEntity) //again over all entities with assignments
        if (!(*ppEntity)->isUsed())
          *ppReorder++ = *ppEntity;

      ppEntityEnd = mStateTemplate.endFixed();

      for (; ppEntity != ppEntityEnd; ++ppEntity, ++ppReorder) //over all fixed entities
        *ppReorder = *ppEntity;

      mStateTemplate.reorder(Reorder);
      mReorderNeeded = false;

      // We need to recompile as pointers to values may have changed
      ppEntity = mStateTemplate.beginIndependent();
      ppEntityEnd = mStateTemplate.endFixed();

      for (; ppEntity != ppEntityEnd; ++ppEntity)
        (*ppEntity)->compile();

      itReaction = mSteps.begin();
      endReaction = mSteps.end();

      for (; itReaction != endReaction; ++itReaction)
        (*itReaction)->compile();

      // The compile might have broken some refresh pointers we need to rebuild the constant sequence
      buildApplyInitialValuesSequence();
    }

  std::set< const CCopasiObject * > UpToDate;

  try
    {
      mSimulatedRefreshes = CCopasiObject::buildUpdateSequence(mSimulatedUpToDateObjects, UpToDate);
    }
  catch (...)
    {
      mSimulatedRefreshes.clear();
      success = false;
    }

  return success;
}

bool CModel::buildApplyInitialValuesSequence()
{
  bool success = true;

  mApplyInitialValuesRefreshes.clear();

  std::set< const CCopasiObject * > Objects;

  const CMetab * pMetab;

  CModelEntity ** ppEntity =  mStateTemplate.beginIndependent();
  CModelEntity ** ppEntityEnd = mStateTemplate.endFixed();

  for (; ppEntity != ppEntityEnd; ++ppEntity)
    {
      if ((*ppEntity)->getStatus() == ODE)
        {
          // TODO We need only to calculate rates which are constant since the other will
          // be updated by the simulation request.
          Objects.insert((*ppEntity)->getRateReference());
        }

      // Species of type assignment have a second pseudo state value the concentration,
      // which always can be directly calculated.
      if ((*ppEntity)->getStatus() == ASSIGNMENT &&
          (pMetab = dynamic_cast< const CMetab * >(*ppEntity)) != NULL)
        {
          mApplyInitialValuesRefreshes.push_back(pMetab->getConcentrationReference()->getApplyInitialValueRefresh());
        }
    }

  std::set< const CCopasiObject * > UpToDate;

  try
    {
      std::vector< Refresh * > RateRefreshes =
        CCopasiObject::buildUpdateSequence(Objects, UpToDate);

      mApplyInitialValuesRefreshes.insert(mApplyInitialValuesRefreshes.end(),
                                          RateRefreshes.begin(),
                                          RateRefreshes.end());
    }
  catch (...)
    {
      mApplyInitialValuesRefreshes.clear();
      success = false;
    }

  return success;
}

bool CModel::buildNonSimulatedSequence()
{
  bool success = true;

  std::set< const CCopasiObject * > Objects;

  // Compartments
  CCopasiVector< CCompartment >::iterator itComp = mCompartments.begin();
  CCopasiVector< CCompartment >::iterator endComp = mCompartments.end();

  for (; itComp != endComp; ++itComp)
    {
      Objects.insert((*itComp)->getValueReference());

      switch ((*itComp)->getStatus())
        {
          case ODE:
            Objects.insert((*itComp)->getRateReference());
            break;

          default:
            break;
        }
    }

  // Metabolites
  CCopasiVector< CMetab >::iterator itMetab = mMetabolites.begin();
  CCopasiVector< CMetab >::iterator endMetab = mMetabolites.end();

  for (; itMetab != endMetab; ++itMetab)
    {
      Objects.insert((*itMetab)->getConcentrationReference());
      Objects.insert((*itMetab)->getValueReference());

      switch ((*itMetab)->getStatus())
        {
          case REACTIONS:
          case ODE:
            Objects.insert((*itMetab)->getObject(CCopasiObjectName("Reference=TransitionTime")));
            Objects.insert((*itMetab)->getObject(CCopasiObjectName("Reference=Rate")));
            Objects.insert((*itMetab)->getRateReference());
            break;

          default:
            break;
        }
    }

  // Reactions
  CCopasiVector< CReaction >::iterator itStep = mSteps.begin();
  CCopasiVector< CReaction >::iterator endStep = mSteps.end();

  for (; itStep != endStep; ++itStep)
    {
      Objects.insert((*itStep)->getObject(CCopasiObjectName("Reference=Flux")));
      Objects.insert((*itStep)->getParticleFluxReference());
    }

  // Model Values
  CCopasiVector< CModelValue >::iterator itValue = mValues.begin();
  CCopasiVector< CModelValue >::iterator endValue = mValues.end();

  for (; itValue != endValue; ++itValue)
    {
      Objects.insert((*itValue)->getValueReference());

      switch ((*itValue)->getStatus())
        {
          case ODE:
            Objects.insert((*itValue)->getRateReference());
            break;

          default:
            break;
        }
    }

  try
    {
      mNonSimulatedRefreshes = CCopasiObject::buildUpdateSequence(Objects, mSimulatedUpToDateObjects);
    }
  catch (...)
    {
      mNonSimulatedRefreshes.clear();
      success = false;
    }

  return success;
}

const CState & CModel::getInitialState() const
{return mInitialState;}

const CState & CModel::getState() const
{return mCurrentState;}

void CModel::setInitialState(const CState & state)
{
  mInitialState = state;

  if (mIsAutonomous &&
      !mCompileIsNecessary)
    mInitialState.setTime(0.0);

  return;
}

void CModel::setState(const CState & state)
{
  mCurrentState = state;

  return;
}

void CModel::updateSimulatedValues(const bool & updateMoieties)
{
  // Depending on which model we are using we need to update
  // the particle numbers for the dependent metabolites.
  if (updateMoieties)
    {
      C_FLOAT64 * pDependent = mCurrentState.beginDependent();
      CCopasiVector< CMoiety >::iterator itMoiety = mMoieties.begin();
      CCopasiVector< CMoiety >::iterator endMoiety = mMoieties.end();

      for (; itMoiety != endMoiety; ++itMoiety, ++pDependent)
        *pDependent = (*itMoiety)->dependentNumber();
    }

  std::vector< Refresh * >::const_iterator itRefresh = mSimulatedRefreshes.begin();
  std::vector< Refresh * >::const_iterator endRefresh = mSimulatedRefreshes.end();

  while (itRefresh != endRefresh)
    (**itRefresh++)();

  // Store the particle fluxes for further calculations
  CCopasiVector< CReaction >::iterator it = mSteps.begin();
  CCopasiVector< CReaction >::iterator end = mSteps.end();

  C_FLOAT64 * pFlux = mParticleFluxes.array();

  for (; it != end; ++it, ++pFlux)
    *pFlux = (*it)->getParticleFlux();
}

void CModel::updateNonSimulatedValues(void)
{
  std::vector< Refresh * >::const_iterator itRefresh = mNonSimulatedRefreshes.begin();
  std::vector< Refresh * >::const_iterator endRefresh = mNonSimulatedRefreshes.end();

  while (itRefresh != endRefresh)
    (**itRefresh++)();

  return;
}

void CModel::calculateDerivatives(C_FLOAT64 * derivatives)
{
  C_FLOAT64 * pTmp = derivatives;

  // First retrieve derivatives of quantities determined by ODE
  // The offset 1 is for the model time which is always the first
  // state variable.
  CModelEntity ** ppIt = mStateTemplate.getEntities() + 1;
  CModelEntity ** ppEnd =
    ppIt + mStateTemplate.getNumIndependent() - mNumMetabolitesReactionIndependent;

  for (; ppIt != ppEnd; ++ppIt, ++pTmp)
    *pTmp = (*ppIt)->getRate();

  // Now calculate derivatives of all metabolites determined by reactions
  char T = 'N';
  C_INT M = 1;
  C_INT N = mNumMetabolitesReaction;
  C_INT K = mSteps.size();
  C_FLOAT64 Alpha = 1.0;
  C_FLOAT64 Beta = 0.0;

  if (K != 0)
    dgemm_(&T, &T, &M, &N, &K, &Alpha, mParticleFluxes.array(), &M,
           mStoiReordered.array(), &K, &Beta, pTmp, &M);
}

void CModel::calculateDerivativesX(C_FLOAT64 * derivativesX)
{
  C_FLOAT64 * pTmp = derivativesX;

  // First retrieve derivatives of quantities determined by ODE
  // The offset 1 is for the model time which is always the first
  // state variable.
  CModelEntity ** ppIt = mStateTemplate.getEntities() + 1;
  CModelEntity ** ppEnd =
    ppIt + mStateTemplate.getNumIndependent() - mNumMetabolitesReactionIndependent;

  for (; ppIt != ppEnd; ++ppIt, ++pTmp)
    *pTmp = (*ppIt)->getRate();

  // Now calculate derivatives of the independent metabolites determined by reactions
  char T = 'N';
  C_INT M = 1;
  C_INT N = mNumMetabolitesReactionIndependent;
  C_INT K = mSteps.size();
  C_FLOAT64 Alpha = 1.0;
  C_FLOAT64 Beta = 0.0;

  if (K != 0)
    dgemm_(&T, &T, &M, &N, &K, &Alpha, mParticleFluxes.array(), &M,
           mRedStoi.array(), &K, &Beta, pTmp, &M);
}

void CModel::calculateElasticityMatrix(const C_FLOAT64 & factor,
                                       const C_FLOAT64 & resolution)
{
  unsigned C_INT32 Col;
  unsigned C_INT32 nCol = mElasticities.numCols();

  C_FLOAT64 * itE;
  C_FLOAT64 * beginE = mElasticities.array();

  CCopasiVector< CReaction >::const_iterator itReaction;
  CCopasiVector< CReaction >::const_iterator beginReaction = mSteps.begin();
  CCopasiVector< CReaction >::const_iterator endReaction = mSteps.end();

  CModelEntity ** itEntity;
  CModelEntity ** beginEntity = mStateTemplate.beginIndependent();
  CModelEntity ** endEntity = mStateTemplate.endDependent();

  for (itEntity = beginEntity, Col = 0; itEntity != endEntity; ++itEntity, ++Col)
    {
      // :TODO: This only works for entities of type metabolites.
      //        The scaling factor for other entities should be 1.
      const C_FLOAT64 invVolume =
        1.0 / static_cast<CMetab *>(*itEntity)->getCompartment()->getValue();
      C_FLOAT64 * pX =
        const_cast<C_FLOAT64 *>(&static_cast<CMetab *>(*itEntity)->getConcentration());

      for (itReaction = beginReaction, itE = beginE + Col;
           itReaction != endReaction;
           ++itReaction, itE += nCol)
        * itE = invVolume * (*itReaction)->calculatePartialDerivative(pX, factor, resolution);
    }

  // DebugFile << "CModel::calculateElasticityMatrix()" << std::endl;
  // DebugFile << mElasticities << std::endl;

  return;
}

void CModel::calculateJacobian(CMatrix< C_FLOAT64 > & jacobian,
                               const C_FLOAT64 & derivationFactor,
                               const C_FLOAT64 & /* resolution */)
{
  unsigned C_INT32 Dim =
    mCurrentState.getNumIndependent() + MNumMetabolitesReactionDependent;
  //Dim now contains the number of entities with ODEs + number of metabs depending on reactions.

  unsigned C_INT32 Col;

  jacobian.resize(Dim, Dim);
  CMatrix< C_FLOAT64 > Jacobian(Dim, Dim);

  C_FLOAT64 Store;
  C_FLOAT64 X1;
  C_FLOAT64 X2;
  C_FLOAT64 InvDelta;

  CVector< C_FLOAT64 > Y1(Dim);
  CVector< C_FLOAT64 > Y2(Dim);

  C_FLOAT64 * pY1;
  C_FLOAT64 * pY2;

  C_FLOAT64 * pX = mCurrentState.beginIndependent();
  C_FLOAT64 * pXEnd = pX + Dim;

  C_FLOAT64 * pJacobian;
  C_FLOAT64 * pJacobianEnd = Jacobian.array() + Dim * Dim;

  for (Col = 0; pX != pXEnd; ++pX, ++Col)
    {
      Store = *pX;

      // We only need to make sure that we do not have an underflow problem
      if (fabs(Store) < 100 * DBL_MIN)
        {
          X1 = 0.0;

          if (Store < 0.0)
            X2 = -200.0 * DBL_MIN;
          else
            X2 = 200.0 * DBL_MIN;;
        }
      else
        {
          X1 = Store * (1.0 + derivationFactor);
          X2 = Store * (1.0 - derivationFactor);
        }

      InvDelta = 1.0 / (X2 - X1);

      *pX = X1;
      updateSimulatedValues(false);
      calculateDerivatives(Y1.array());

      *pX = X2;
      updateSimulatedValues(false);
      calculateDerivatives(Y2.array());

      *pX = Store;

      pJacobian = Jacobian.array() + Col;
      pY1 = Y1.array();
      pY2 = Y2.array();

      for (; pJacobian < pJacobianEnd; pJacobian += Dim, ++pY1, ++pY2)
        * pJacobian = (*pY2 - *pY1) * InvDelta;
    }

  updateSimulatedValues(false);

  //  jacobian = Jacobian;
  //  return;

  // :TODO: this can be incorporated into the above avoiding a temporary matrix.

  // We need to bring the jacobian into the expected order, i.e.,
  // convert it to the user defined order
  unsigned C_INT32 * pPermRow = mJacobianPivot.array();
  unsigned C_INT32 * pPermEnd = pPermRow + mJacobianPivot.size();
  unsigned C_INT32 * pPermCol;

  C_FLOAT64 * pTo;
  pTo = jacobian.array();

  for (; pPermRow < pPermEnd; ++pPermRow)
    {
      pJacobian = Jacobian.array() + *pPermRow * Dim;

      for (pPermCol = mJacobianPivot.array(); pPermCol < pPermEnd; ++pPermCol, ++pTo)
        *pTo = *(pJacobian + *pPermCol);
    }

  // DebugFile << jacobian << std::endl;
}

void CModel::calculateJacobianX(CMatrix< C_FLOAT64 > & jacobianX,
                                const C_FLOAT64 & derivationFactor,
                                const C_FLOAT64 & /* resolution */)
{
  C_FLOAT64 DerivationFactor = std::max(derivationFactor, 100.0 * std::numeric_limits< C_FLOAT64 >::epsilon());

  unsigned C_INT32 Dim = mCurrentState.getNumIndependent();
  unsigned C_INT32 Col;

  jacobianX.resize(Dim, Dim);

  C_FLOAT64 Store;
  C_FLOAT64 X1;
  C_FLOAT64 X2;
  C_FLOAT64 InvDelta;

  CVector< C_FLOAT64 > Y1(Dim);
  CVector< C_FLOAT64 > Y2(Dim);

  C_FLOAT64 * pY1;
  C_FLOAT64 * pY2;

  C_FLOAT64 * pX = mCurrentState.beginIndependent();
  C_FLOAT64 * pXEnd = pX + Dim;

  C_FLOAT64 * pJacobian;
  C_FLOAT64 * pJacobianEnd = jacobianX.array() + Dim * Dim;

  for (Col = 0; pX != pXEnd; ++pX, ++Col)
    {
      Store = *pX;

      // We only need to make sure that we do not have an underflow problem
      if (fabs(Store) < 100 * DBL_MIN)
        {
          X1 = 0.0;

          if (Store < 0.0)
            X2 = -200.0 * DBL_MIN;
          else
            X2 = 200.0 * DBL_MIN;;
        }
      else
        {
          X1 = Store * (1.0 + DerivationFactor);
          X2 = Store * (1.0 - DerivationFactor);
        }

      InvDelta = 1.0 / (X2 - X1);

      *pX = X1;
      updateSimulatedValues(true);
      calculateDerivativesX(Y1.array());

      *pX = X2;
      updateSimulatedValues(true);
      calculateDerivativesX(Y2.array());

      *pX = Store;

      pJacobian = jacobianX.array() + Col;
      pY1 = Y1.array();
      pY2 = Y2.array();

      for (; pJacobian < pJacobianEnd; pJacobian += Dim, ++pY1, ++pY2)
        * pJacobian = (*pY2 - *pY1) * InvDelta;
    }

  updateSimulatedValues(true);
}

C_FLOAT64 CModel::calculateDivergence() const
{
  fatalError(); //not yet implemented
  return 0.0;
}

bool CModel::setVolumeUnit(const std::string & name)
{
  return setVolumeUnit(toEnum(name.c_str(), VolumeUnitNames, ml));
}

bool CModel::setVolumeUnit(const CModel::VolumeUnit & unit)
{
  mVolumeUnit = unit;
  return true;
}

std::string CModel::getVolumeUnitName() const
{
  return VolumeUnitNames[mVolumeUnit];
}

CModel::VolumeUnit CModel::getVolumeUnitEnum() const
{
  return mVolumeUnit;
}

//****

bool CModel::setAreaUnit(const std::string & name)
{
  return setAreaUnit(toEnum(name.c_str(), AreaUnitNames, m2));
}

bool CModel::setAreaUnit(const CModel::AreaUnit & unit)
{
  mAreaUnit = unit;
  return true;
}

std::string CModel::getAreaUnitName() const
{
  return AreaUnitNames[mAreaUnit];
}

CModel::AreaUnit CModel::getAreaUnitEnum() const
{
  return mAreaUnit;
}

//****
bool CModel::setLengthUnit(const std::string & name)
{
  return setLengthUnit(toEnum(name.c_str(), LengthUnitNames, m));
}

bool CModel::setLengthUnit(const CModel::LengthUnit & unit)
{
  mLengthUnit = unit;
  return true;
}

std::string CModel::getLengthUnitName() const
{
  return LengthUnitNames[mLengthUnit];
}

CModel::LengthUnit CModel::getLengthUnitEnum() const
{
  return mLengthUnit;
}

//****

bool CModel::setTimeUnit(const std::string & name)
{
  return setTimeUnit(toEnum(name.c_str(), TimeUnitNames, s));
}

bool CModel::setTimeUnit(const CModel::TimeUnit & unit)
{
  mTimeUnit = unit;
  return true;
}

std::string CModel::getTimeUnitName() const
{
  return TimeUnitNames[mTimeUnit];
}

CModel::TimeUnit CModel::getTimeUnitEnum() const
{
  return mTimeUnit;
}

//****

bool CModel::setQuantityUnit(const std::string & name)
{
  QuantityUnit unit = toEnum(name.c_str(), QuantityUnitNames, OldXML);

  if (unit == OldXML)
    unit = toEnum(name.c_str(), QuantityUnitOldXMLNames, mMol);

  return setQuantityUnit(unit);
}

bool CModel::setQuantityUnit(const CModel::QuantityUnit & unit)
{
  bool success = true;
  mQuantityUnit = unit;

  switch (unit)
    {
      case Mol:
        mQuantity2NumberFactor = mAvogadro;
        break;

      case mMol:
        mQuantity2NumberFactor = mAvogadro * 1E-3;
        break;

      case microMol:
        mQuantity2NumberFactor = mAvogadro * 1E-6;
        break;

      case nMol:
        mQuantity2NumberFactor = mAvogadro * 1E-9;
        break;

      case pMol:
        mQuantity2NumberFactor = mAvogadro * 1E-12;
        break;

      case fMol:
        mQuantity2NumberFactor = mAvogadro * 1E-15;
        break;

      case number:
        mQuantity2NumberFactor = 1.0;
        break;

      case dimensionlessQuantity:
        mQuantity2NumberFactor = 1.0;
        break;

      default:
        mQuantityUnit = number;
        mQuantity2NumberFactor = 1.0;
        success = false;
    }

  mNumber2QuantityFactor = 1.0 / mQuantity2NumberFactor;

  //adapt particle numbers
  C_INT32 i, imax = mMetabolites.size();

  for (i = 0; i < imax; ++i)
    {
      //update particle numbers
      mMetabolites[i]->setInitialConcentration(mMetabolites[i]->getInitialConcentration());
      mMetabolites[i]->setConcentration(mMetabolites[i]->getConcentration());
    }

  return success;
}

std::string CModel::getQuantityUnitName() const
{
  return QuantityUnitNames[mQuantityUnit];
}

std::string CModel::getQuantityUnitOldXMLName() const
{
  return QuantityUnitOldXMLNames[mQuantityUnit];
}

CModel::QuantityUnit CModel::getQuantityUnitEnum() const
{
  return mQuantityUnit;
}

void CModel::setModelType(const CModel::ModelType & modelType)
{mType = modelType;}

const CModel::ModelType & CModel::getModelType() const
{return mType;}

void CModel::setAvogadro(const C_FLOAT64 & avogadro)
{
  mAvogadro = avogadro;
}

const C_FLOAT64 & CModel::getAvogadro() const
{
  return mAvogadro;
}

const C_FLOAT64 & CModel::getQuantity2NumberFactor() const
{return mQuantity2NumberFactor;}

const C_FLOAT64 & CModel::getNumber2QuantityFactor() const
{return mNumber2QuantityFactor;}

//*****

//**********************************************************************

bool CModel::appendDependentModelObjects(const std::set< const CCopasiObject * > & deletedObjects,
    std::set< const CCopasiObject * > & dependentReactions,
    std::set< const CCopasiObject * > & dependentMetabolites,
    std::set< const CCopasiObject * > & dependentCompartments,
    std::set< const CCopasiObject * > & dependentModelValues,
    std::set< const CCopasiObject * > & dependentEvents) const
{
  // We need a local copy since we recursively add deleted objects.
  std::set< const CCopasiObject * > DeletedObjects = deletedObjects;

  bool ObjectsAppended = false;
  bool DeleteObjects = DeletedObjects.size() > 0;

  // This is this implemented recursively. Since deleting a container may result
  // in the deletion of objects not dependent on the original set of deleted objects.

  while (DeleteObjects)
    {
      DeleteObjects = false;

      DeleteObjects |= appendDependentReactions(DeletedObjects, dependentReactions);

      if (dependentReactions.size() > 0)
        {
          std::set< const CCopasiObject * >::const_iterator it, itEnd = dependentReactions.end();

          for (it = dependentReactions.begin(); it != itEnd; ++it)
            if (DeletedObjects.find(*it) == DeletedObjects.end())
              {
                std::set< const CCopasiObject * > AdditionalObjects =
                  static_cast< const CReaction * >(*it)->getDeletedObjects();

                std::set< const CCopasiObject * >::const_iterator itDeleted = AdditionalObjects.begin();
                std::set< const CCopasiObject * >::const_iterator endDeleted = AdditionalObjects.end();

                for (; itDeleted != endDeleted; ++itDeleted)
                  DeletedObjects.insert(*itDeleted);
              }
        }

      DeleteObjects |= appendDependentMetabolites(DeletedObjects, dependentMetabolites);

      if (dependentMetabolites.size() > 0)
        {
          std::set< const CCopasiObject * >::const_iterator it, itEnd = dependentMetabolites.end();

          for (it = dependentMetabolites.begin(); it != itEnd; ++it)
            if (DeletedObjects.find(*it) == DeletedObjects.end())
              {
                std::set< const CCopasiObject * > AdditionalObjects =
                  static_cast< const CMetab * >(*it)->getDeletedObjects();

                std::set< const CCopasiObject * >::const_iterator itDeleted = AdditionalObjects.begin();
                std::set< const CCopasiObject * >::const_iterator endDeleted = AdditionalObjects.end();

                for (; itDeleted != endDeleted; ++itDeleted)
                  DeletedObjects.insert(*itDeleted);
              }
        }

      DeleteObjects |= appendDependentModelValues(DeletedObjects, dependentModelValues);

      if (dependentModelValues.size() > 0)
        {
          std::set< const CCopasiObject * >::const_iterator it, itEnd = dependentModelValues.end();

          for (it = dependentModelValues.begin(); it != itEnd; ++it)
            if (DeletedObjects.find(*it) == DeletedObjects.end())
              {
                std::set< const CCopasiObject * > AdditionalObjects =
                  static_cast< const CModelValue * >(*it)->getDeletedObjects();

                std::set< const CCopasiObject * >::const_iterator itDeleted = AdditionalObjects.begin();
                std::set< const CCopasiObject * >::const_iterator endDeleted = AdditionalObjects.end();

                for (; itDeleted != endDeleted; ++itDeleted)
                  DeletedObjects.insert(*itDeleted);
              }
        }

      DeleteObjects |= appendDependentCompartments(DeletedObjects, dependentCompartments);

      if (dependentCompartments.size() > 0)
        {
          std::set< const CCopasiObject * >::const_iterator it, itEnd = dependentCompartments.end();

          for (it = dependentCompartments.begin(); it != itEnd; ++it)
            if (DeletedObjects.find(*it) == DeletedObjects.end())
              {
                std::set< const CCopasiObject * > AdditionalObjects =
                  static_cast< const CCompartment * >(*it)->getDeletedObjects();

                std::set< const CCopasiObject * >::const_iterator itDeleted = AdditionalObjects.begin();
                std::set< const CCopasiObject * >::const_iterator endDeleted = AdditionalObjects.end();

                for (; itDeleted != endDeleted; ++itDeleted)
                  DeletedObjects.insert(*itDeleted);
              }
        }

      DeleteObjects |= appendDependentEvents(DeletedObjects, dependentEvents);

      ObjectsAppended |= DeleteObjects;
    }

  return ObjectsAppended;
}

bool CModel::appendDependentReactions(std::set< const CCopasiObject * > candidates,
                                      std::set< const CCopasiObject * > & dependents) const
{
  const_cast< CModel * >(this)->compileIfNecessary(NULL);

  size_t Size = dependents.size();

  CCopasiVectorN< CReaction >::const_iterator it = mSteps.begin();
  CCopasiVectorN< CReaction >::const_iterator end = mSteps.end();

  std::set< const CCopasiObject * >::const_iterator itSet;
  std::set< const CCopasiObject * >::const_iterator endSet;

  for (; it != end; ++it)
    if (candidates.find(*it) == candidates.end())
      {
        std::set< const CCopasiObject * > Ignored;

        // We need to ignore our own local reaction parameters.
        CCopasiParameterGroup::index_iterator itParameter = (*it)->getParameters().beginIndex();
        CCopasiParameterGroup::index_iterator endParameter = (*it)->getParameters().endIndex();
        std::set< const CCopasiObject * >::iterator itIgnored;

        for (; itParameter != endParameter; ++itParameter)
          {
            if ((itIgnored = candidates.find((*itParameter)->getObject(CCopasiObjectName("Reference=Value")))) != candidates.end())
              {
                Ignored.insert(*itIgnored);
                candidates.erase(itIgnored);
              }
          }

        if ((*it)->dependsOn(candidates))
          {
            dependents.insert((*it));
            continue;
          }

        std::set< const CCopasiObject * > DeletedObjects = (*it)->getDeletedObjects();
        itSet = DeletedObjects.begin();
        endSet = DeletedObjects.end();

        for (; itSet != endSet; ++itSet)
          if (candidates.find(*itSet) == candidates.end() &&
              (*itSet)->dependsOn(candidates))
            {
              dependents.insert((*it));
              break;
            }

        // Add the ignored parameters back.
        std::set< const CCopasiObject * >::iterator endIgnored = Ignored.end();

        for (itIgnored = Ignored.begin(); itIgnored != endIgnored; ++itIgnored)
          candidates.insert(*itIgnored);
      }

  return Size < dependents.size();
}

bool CModel::appendDependentMetabolites(std::set< const CCopasiObject * > candidates,
                                        std::set< const CCopasiObject * > & dependents) const
{
  const_cast< CModel * >(this)->compileIfNecessary(NULL);

  size_t Size = dependents.size();

  CCopasiVectorN< CCompartment >::const_iterator itComp = mCompartments.begin();
  CCopasiVectorN< CCompartment >::const_iterator endComp = mCompartments.end();

  CCopasiVectorN< CMetab >::const_iterator it;
  CCopasiVectorN< CMetab >::const_iterator end;

  std::set< const CCopasiObject * >::const_iterator itSet;
  std::set< const CCopasiObject * >::const_iterator endSet;

  for (; itComp != endComp; ++itComp)
    {
      it = (*itComp)->getMetabolites().begin();
      end = (*itComp)->getMetabolites().end();

      for (; it != end; ++it)
        if (candidates.find((*it)->getCompartment()) != candidates.end())
          dependents.insert((*it));
        else if (candidates.find(*it) == candidates.end())
          {
            if (candidates.find((*it)->getCompartment()->getObject(CCopasiObjectName("Reference=Volume"))) != candidates.end() ||
                (*it)->dependsOn(candidates))
              {
                dependents.insert((*it));
                continue;
              }

            std::set< const CCopasiObject * > DeletedObjects = (*it)->getDeletedObjects();

            if ((*it)->getStatus() == REACTIONS)
              {
                DeletedObjects.erase((*it)->getObject(CCopasiObjectName("Reference=Rate")));
                DeletedObjects.erase((*it)->getObject(CCopasiObjectName("Reference=ParticleNumberRate")));
                DeletedObjects.erase((*it)->getObject(CCopasiObjectName("Reference=TransitionTime")));
              }

            itSet = DeletedObjects.begin();
            endSet = DeletedObjects.end();

            for (; itSet != endSet; ++itSet)
              if (candidates.find(*itSet) == candidates.end() &&
                  (*itSet)->dependsOn(candidates))
                {
                  dependents.insert((*it));
                  break;
                }
          }
    }

  return Size < dependents.size();
}

bool CModel::appendDependentCompartments(std::set< const CCopasiObject * > candidates,
    std::set< const CCopasiObject * > & dependents) const
{
  const_cast< CModel * >(this)->compileIfNecessary(NULL);

  size_t Size = dependents.size();

  CCopasiVectorN< CCompartment >::const_iterator it = mCompartments.begin();
  CCopasiVectorN< CCompartment >::const_iterator end = mCompartments.end();

  std::set< const CCopasiObject * >::const_iterator itSet;
  std::set< const CCopasiObject * >::const_iterator endSet;

  for (; it != end; ++it)
    if (candidates.find(*it) == candidates.end())
      {
        if ((*it)->dependsOn(candidates))
          {
            dependents.insert((*it));
            continue;
          }

        std::set< const CCopasiObject * > DeletedObjects = (*it)->getDeletedObjects();
        itSet = DeletedObjects.begin();
        endSet = DeletedObjects.end();

        for (; itSet != endSet; ++itSet)
          if (candidates.find(*itSet) == candidates.end() &&
              (*itSet)->dependsOn(candidates))
            {
              dependents.insert((*it));
              break;
            }
      }

  return Size < dependents.size();
}

bool CModel::appendDependentModelValues(std::set< const CCopasiObject * > candidates,
                                        std::set< const CCopasiObject * > & dependents) const
{
  const_cast< CModel * >(this)->compileIfNecessary(NULL);

  size_t Size = dependents.size();

  CCopasiVectorN< CModelValue >::const_iterator it = mValues.begin();
  CCopasiVectorN< CModelValue >::const_iterator end = mValues.end();

  std::set< const CCopasiObject * >::const_iterator itSet;
  std::set< const CCopasiObject * >::const_iterator endSet;

  for (; it != end; ++it)
    if (candidates.find(*it) == candidates.end())
      {
        if ((*it)->dependsOn(candidates))
          {
            dependents.insert((*it));
            continue;
          }

        std::set< const CCopasiObject * > DeletedObjects = (*it)->getDeletedObjects();
        itSet = DeletedObjects.begin();
        endSet = DeletedObjects.end();

        for (; itSet != endSet; ++itSet)
          if (candidates.find(*itSet) == candidates.end() &&
              (*itSet)->dependsOn(candidates))
            {
              dependents.insert((*it));
              break;
            }
      }

  return Size < dependents.size();
}

bool CModel::appendDependentEvents(std::set< const CCopasiObject * > candidates,
                                   std::set< const CCopasiObject * > & dependents) const
{
  const_cast< CModel * >(this)->compileIfNecessary(NULL);

  size_t Size = dependents.size();

  CCopasiVectorN< CEvent >::const_iterator it = mEvents.begin();
  CCopasiVectorN< CEvent >::const_iterator end = mEvents.end();

  std::set< const CCopasiObject * >::const_iterator itSet;
  std::set< const CCopasiObject * >::const_iterator endSet;

  for (; it != end; ++it)
    if (candidates.find(*it) == candidates.end())
      {
        if ((*it)->dependsOn(candidates))
          {
            dependents.insert((*it));
          }
      }

  return Size < dependents.size();
}

//**********************************************************************

CMetab* CModel::createMetabolite(const std::string & name,
                                 const std::string & compartment,
                                 const C_FLOAT64 & iconc,
                                 const CMetab::Status & status)
{
  unsigned C_INT32 Index;

  if (mCompartments.size() == 0)
    return NULL;

  if (compartment == "")
    Index = 0;
  else if ((Index = mCompartments.getIndex(compartment)) == C_INVALID_INDEX)
    return NULL;

  if (mCompartments[Index]->getMetabolites().getIndex(name) != C_INVALID_INDEX)
    return NULL;

  CMetab * pMetab = new CMetab(name);

  if (!mCompartments[Index]->addMetabolite(pMetab))
    {
      delete pMetab;
      return NULL;
    }

  pMetab->setStatus(status);
  pMetab->setInitialConcentration(iconc);
  pMetab->refreshInitialValue();

  if (!mMetabolites.add(pMetab))
    return NULL;

  mCompileIsNecessary = true;

  return pMetab;
}

bool CModel::removeMetabolite(const unsigned C_INT32 index,
                              const bool & recursive)
{
  const CMetab* pMetabolite = getMetabolites()[index];
  return removeMetabolite(pMetabolite, recursive);
}

bool CModel::removeMetabolite(const std::string & key,
                              const bool & recursive)
{
  CMetab* pMetabolite =
    dynamic_cast< CMetab * >(CCopasiRootContainer::getKeyFactory()->get(key));
  return removeMetabolite(pMetabolite, recursive);
}

bool CModel::removeMetabolite(const CMetab* pMetabolite,
                              const bool & recursive)
{
  if (!pMetabolite)
    return false;

  if (recursive)
    {
      removeDependentModelObjects(pMetabolite->getDeletedObjects());
    }

  /* Assure that all references are removed */
  mMetabolites.remove((CMetab *)pMetabolite);
  mMetabolitesX.remove((CMetab *)pMetabolite);

  pdelete(pMetabolite);

  clearMoieties();
  mCompileIsNecessary = true;

  return true;
}

CCompartment* CModel::createCompartment(const std::string & name,
                                        const C_FLOAT64 & volume)
{
  // check if there is already a volume with this name
  if (mCompartments.getIndex(name) != C_INVALID_INDEX)
    return NULL;

  CCompartment * cpt = new CCompartment(name);

  cpt->setInitialValue(volume);
  //cpt->setVolume(volume);

  if (!mCompartments.add(cpt, true))
    {
      delete cpt;
      return NULL;
    }

  mCompileIsNecessary = true;
  return cpt;
}

bool CModel::removeCompartment(const unsigned C_INT32 index,
                               const bool & recursive)
{
  const CCompartment * pCompartment = getCompartments()[index];
  return removeCompartment(pCompartment, recursive);
}

bool CModel::removeCompartment(const std::string & key,
                               const bool & recursive)
{
  CCompartment *pCompartment =
    dynamic_cast< CCompartment * >(CCopasiRootContainer::getKeyFactory()->get(key));
  return removeCompartment(pCompartment, recursive);
}

bool CModel::removeCompartment(const CCompartment * pCompartment,
                               const bool & recursive)
{
  if (!pCompartment)
    return false;

  if (recursive)
    {
      removeDependentModelObjects(pCompartment->getDeletedObjects());
    }

  //Check if Compartment with that name exists
  unsigned C_INT32 index =
    mCompartments.CCopasiVector< CCompartment >::getIndex(pCompartment);

  if (index == C_INVALID_INDEX)
    return false;

  mCompartments.CCopasiVector< CCompartment >::remove(index);

  mCompileIsNecessary = true;

  return true;
}

CReaction* CModel::createReaction(const std::string & name)
{
  if (mSteps.getIndex(name) != C_INVALID_INDEX)
    return NULL;

  CReaction * pReaction = new CReaction(name);

  if (!mSteps.add(pReaction, true))
    {
      delete pReaction;
      return NULL;
    }

  mCompileIsNecessary = true;
  return pReaction;
}

bool CModel::removeReaction(const std::string & key,
                            const bool & recursive)
{
  CReaction * pReaction =
    dynamic_cast< CReaction * >(CCopasiRootContainer::getKeyFactory()->get(key));
  return removeReaction(pReaction, recursive);
}

bool CModel::removeReaction(const unsigned C_INT32 index,
                            const bool & recursive)
{
  const CReaction * pReaction = getReactions()[index];
  return removeReaction(pReaction, recursive);
}

bool CModel::removeReaction(const CReaction * pReaction,
                            const bool & recursive)
{
  if (!pReaction)
    return false;

  if (recursive)
    {
      removeDependentModelObjects(pReaction->getDeletedObjects());
    }

  //Check if Reaction exists
  unsigned C_INT32 index =
    mSteps.CCopasiVector< CReaction >::getIndex(pReaction);

  if (index == C_INVALID_INDEX)
    return false;

  mSteps.CCopasiVector< CReaction >::remove(index);

  clearMoieties();
  mCompileIsNecessary = true;

  return true;
}

bool CModel::removeLocalReactionParameter(const std::string & key,
    const bool & recursive)
{
  CCopasiParameter * pParameter =
    dynamic_cast< CCopasiParameter * >(CCopasiRootContainer::getKeyFactory()->get(key));

  if (pParameter == NULL)
    return false;

  if (recursive)
    {
      std::set< const CCopasiObject * > DeletedObjects;
      DeletedObjects.insert(pParameter->getObject(CCopasiObjectName("Reference=Value")));

      removeDependentModelObjects(DeletedObjects);
    }

  return true;
}

CModelValue* CModel::createModelValue(const std::string & name,
                                      const C_FLOAT64 & value)
{
  // check if there is already a value with this name
  if (mValues.getIndex(name) != C_INVALID_INDEX)
    return NULL;

  CModelValue * cmv = new CModelValue(name);

  cmv->setInitialValue(value);
  cmv->setValue(value);

  if (!mValues.add(cmv, true))
    {
      delete cmv;
      return NULL;
    }

  mCompileIsNecessary = true;
  return cmv;
}

void CModel::removeDependentModelObjects(const std::set<const CCopasiObject*> & deletedObjects)
{
  std::set<const CCopasiObject*> Reactions;
  std::set<const CCopasiObject*> Metabolites;
  std::set<const CCopasiObject*> Values;
  std::set<const CCopasiObject*> Compartments;
  std::set<const CCopasiObject*> Events;

  appendDependentModelObjects(deletedObjects, Reactions, Metabolites, Compartments, Values, Events);

  std::set<const CCopasiObject*>::const_iterator it, end;

  for (it = Reactions.begin(), end = Reactions.end(); it != end; ++it)
    removeReaction((*it)->getKey(), false);

  for (it = Metabolites.begin(), end = Metabolites.end(); it != end; ++it)
    removeMetabolite((*it)->getKey(), false);

  for (it = Compartments.begin(), end = Compartments.end(); it != end; ++it)
    removeCompartment((*it)->getKey(), false);

  for (it = Values.begin(), end = Values.end(); it != end; ++it)
    removeModelValue((*it)->getKey(), false);

  for (it = Events.begin(), end = Events.end(); it != end; ++it)
    removeEvent((*it)->getKey(), false);

  return;
}

bool CModel::removeModelValue(const unsigned C_INT32 index,
                              const bool & recursive)
{
  const CModelValue * pMV = getModelValues()[index];
  return removeModelValue(pMV, recursive);
}
bool CModel::removeModelValue(const std::string & key,
                              const bool & recursive)
{
  CModelValue *pModelValue =
    dynamic_cast< CModelValue * >(CCopasiRootContainer::getKeyFactory()->get(key));
  return removeModelValue(pModelValue, recursive);
}

bool CModel::removeModelValue(const CModelValue * pModelValue,
                              const bool & recursive)
{
  if (!pModelValue)
    return false;

  if (recursive)
    {
      removeDependentModelObjects(pModelValue->getDeletedObjects());
    }

  //Check if Value with that name exists
  unsigned C_INT32 index =
    mValues.CCopasiVector< CModelValue >::getIndex(pModelValue);

  if (index == C_INVALID_INDEX)
    return false;

  mValues.CCopasiVector< CModelValue >::remove(index);

  mCompileIsNecessary = true;

  return true;
}

CEvent* CModel::createEvent(const std::string & name)
{
  if (mEvents.getIndex(name) != C_INVALID_INDEX)
    return NULL;

  CEvent * pEvent = new CEvent(name, this);

  // Assure that the event order is unique.
  // We assume that the existing events are ordered consecutively without
  // gaps.
  unsigned C_INT32 Order = 0;

  CCopasiVectorN< CEvent >::const_iterator it = mEvents.begin();
  CCopasiVectorN< CEvent >::const_iterator end = mEvents.end();

  for (; it != end; ++it)
    {
      if (Order < (*it)->getOrder())
        {
          Order = (*it)->getOrder();
        }
    }

  pEvent->setOrder(Order + 1, false);

  if (!mEvents.add(pEvent, true))
    {
      delete pEvent;
      return NULL;
    }

  mCompileIsNecessary = true;
  return pEvent;
}

bool CModel::removeEvent(const unsigned C_INT32 index,
                         const bool & recursive)
{
  const CEvent * pEvent = mEvents[index];

  return removeEvent(pEvent, recursive);
}

bool CModel::removeEvent(const std::string & key,
                         const bool & recursive)
{
  CEvent * pEvent = dynamic_cast< CEvent * >(CCopasiRootContainer::getKeyFactory()->get(key));

  return removeEvent(pEvent, recursive);
}

bool CModel::removeEvent(const CEvent * pEvent,
                         const bool & /* recursive */)
{
  if (!pEvent)
    return false;

  //Check if Event exists
  unsigned C_INT32 index =
    mEvents.CCopasiVector< CEvent >::getIndex(pEvent);

  if (index == C_INVALID_INDEX)
    return false;

  mEvents.CCopasiVector< CEvent >::remove(index);

  clearMoieties();

  mCompileIsNecessary = true;

  return true;
}

void CModel::synchronizeEventOrder(const CEvent * pEvent,
                                   const unsigned C_INT32 newOrder)
{
  const unsigned C_INT32 & OldOrder = pEvent->getOrder();

  // If the OldOrder is the default for newly created events
  // we do nothing. This assumes that whenever CEvent::setOrder is called
  // for newly created event we assure that the order is unique.
  if (OldOrder == C_INVALID_INDEX)
    {
      return;
    }

  CCopasiVectorN< CEvent >::iterator it = mEvents.begin();
  CCopasiVectorN< CEvent >::iterator end = mEvents.end();

  if (newOrder < OldOrder)
    {
      // We need to increase the order of the events which are in the
      // interval [newOrder, OldOrder)
      for (; it != end; ++it)
        {
          const unsigned C_INT32 & Order = (*it)->getOrder();

          if (newOrder <= Order && Order < OldOrder)
            {
              (*it)->setOrder(Order + 1, false);
            }
        }
    }
  else if (newOrder > OldOrder)
    {
      // We need to decrease the order of the events which are in the
      // interval (OldOrder, newOrder]
      for (; it != end; ++it)
        {
          const unsigned C_INT32 & Order = (*it)->getOrder();

          if (OldOrder < Order && Order <= newOrder)
            {
              (*it)->setOrder(Order - 1, false);
            }
        }
    }
}

//*****************************************************************

bool CModel::convert2NonReversible()
{
  // TODO check if there are any reversible reactions
  // TODO warn the user
  // TODO tell the GUI about changes -> not from here
  // TODO generate report ?
  // TODO check if newly generated reaction names are valid
  // TODO map, so that the same function is split only once

  bool success = true;

  std::vector<std::string> reactionsToDelete;

  CReaction *reac0, *reac1, *reac2;
  CReactionInterface ri1(this), ri2(this);
  std::string fn, rn1, rn2;

  //CModel* model = dynamic_cast< CModel * >(CCopasiRootContainer::getKeyFactory()->get(objKey));
  //if (!model) return false;

  CCopasiVectorN< CReaction > & steps = this->getReactions();

  unsigned C_INT32 i, imax = steps.size();

  for (i = 0; i < imax; ++i)
    if (steps[i]->isReversible())
      {
        reac0 = steps[i];
        rn1 = reac0->getObjectName() + " (forward)";
        rn2 = reac0->getObjectName() + " (backward)";

        fn = reac0->getFunction()->getObjectName();

        const CFunction* pFunc = reac0->getFunction();

        bool massaction = (fn == "Mass action (reversible)");

        std::pair<CFunction *, CFunction *> tmp;

        if (massaction)
          {
            //set functions to mass action (irrev)
            tmp.first = dynamic_cast<CFunction*>
                        (CCopasiRootContainer::getFunctionList()-> findFunction("Mass action (irreversible)"));
            assert(tmp.first);
            tmp.second = tmp.first;
          }
        else //not mass action
          {
            //try splitting
            tmp = pFunc->splitFunction(NULL, pFunc->getObjectName() + " (forward part)",
                                       pFunc->getObjectName() + " (backward part)");

            if ((tmp.first == NULL) || (tmp.second == NULL))
              {
                // Create a message that the conversion for this reaction failed.
                CCopasiMessage(CCopasiMessage::ERROR, MCReaction + 12,
                               reac0->getObjectName().c_str(), fn.c_str());
                success = false;

                pdelete(tmp.first);
                pdelete(tmp.second);
                continue;
              }

            if (tmp.first) CCopasiRootContainer::getFunctionList()->addAndAdaptName(tmp.first);

            if (tmp.second) CCopasiRootContainer::getFunctionList()->addAndAdaptName(tmp.second);
          }

        C_INT32 i, imax;

        //**** create 1st reaction.
        reac1 = createReaction(rn1);
        reac1->setReversible(false);
        //substrates
        imax = reac0->getChemEq().getSubstrates().size();

        for (i = 0; i < imax; ++i)
          reac1->addSubstrate(reac0->getChemEq().getSubstrates()[i]->getMetaboliteKey(),
                              reac0->getChemEq().getSubstrates()[i]->getMultiplicity());

        //products
        imax = reac0->getChemEq().getProducts().size();

        for (i = 0; i < imax; ++i)
          reac1->addProduct(reac0->getChemEq().getProducts()[i]->getMetaboliteKey(),
                            reac0->getChemEq().getProducts()[i]->getMultiplicity());

        //function
        reac1->setFunction(tmp.first);

        //**** create 2nd reaction.
        reac2 = createReaction(rn2);
        reac2->setReversible(false);
        //substrates -> products
        imax = reac0->getChemEq().getSubstrates().size();

        for (i = 0; i < imax; ++i)
          reac2->addProduct(reac0->getChemEq().getSubstrates()[i]->getMetaboliteKey(),
                            reac0->getChemEq().getSubstrates()[i]->getMultiplicity());

        //products -> substrates
        imax = reac0->getChemEq().getProducts().size();

        for (i = 0; i < imax; ++i)
          reac2->addSubstrate(reac0->getChemEq().getProducts()[i]->getMetaboliteKey(),
                              reac0->getChemEq().getProducts()[i]->getMultiplicity());

        //function
        reac2->setFunction(tmp.second);

        //mapping for both reactions
        if (massaction)
          {
            // the parameter names of the massaction kinetics are hardcoded here.
            if (reac0->isLocalParameter("k1"))
              reac1->setParameterValue("k1", reac0->getParameterValue("k1"));
            else
              reac1->setParameterMapping("k1", reac0->getParameterMapping("k1")[0]);

            reac1->setParameterMappingVector("substrate", reac0->getParameterMapping("substrate"));

            if (reac0->isLocalParameter("k2"))
              reac2->setParameterValue("k1", reac0->getParameterValue("k2"));
            else
              reac2->setParameterMapping("k1", reac0->getParameterMapping("k2")[0]);

            reac2->setParameterMappingVector("substrate", reac0->getParameterMapping("product"));
          }
        else //not mass action
          {
            const CFunctionParameters & fps = reac0->getFunctionParameters();
            imax = fps.size();

            for (i = 0; i < imax; ++i)
              {
                const CFunctionParameter * fp = fps[i];
                assert(fp);
                assert(fp->getType() == CFunctionParameter::FLOAT64);

                switch (fp->getUsage())
                  {
                    case CFunctionParameter::SUBSTRATE:
                    case CFunctionParameter::PRODUCT:
                    case CFunctionParameter::MODIFIER:
                      reac1->setParameterMapping(fp->getObjectName(),
                                                 reac0->getParameterMapping(fp->getObjectName())[0]);
                      reac2->setParameterMapping(fp->getObjectName(),
                                                 reac0->getParameterMapping(fp->getObjectName())[0]);
                      break;

                    case CFunctionParameter::PARAMETER:

                      if (reac0->isLocalParameter(fp->getObjectName()))
                        {
                          reac1->setParameterValue(fp->getObjectName(),
                                                   reac0->getParameterValue(fp->getObjectName()));
                          reac2->setParameterValue(fp->getObjectName(),
                                                   reac0->getParameterValue(fp->getObjectName()));
                        }
                      else
                        {
                          reac1->setParameterMapping(fp->getObjectName(),
                                                     reac0->getParameterMapping(fp->getObjectName())[0]);
                          reac2->setParameterMapping(fp->getObjectName(),
                                                     reac0->getParameterMapping(fp->getObjectName())[0]);
                        }

                      break;

                    default:
                      reac1->setParameterMapping(fp->getObjectName(),
                                                 reac0->getParameterMapping(fp->getObjectName())[0]);
                      reac2->setParameterMapping(fp->getObjectName(),
                                                 reac0->getParameterMapping(fp->getObjectName())[0]);
                      break;
                  }
              }
          }

        reac1->compile();
        reac2->compile();

        //remove the old reaction
        reactionsToDelete.push_back(reac0->getObjectName());
      }

  imax = reactionsToDelete.size();

  for (i = 0; i < imax; ++i)
    steps.remove(reactionsToDelete[i]);

  return success;
}

//**********************************************************************

void CModel::initObjects()
{
  mKey = CCopasiRootContainer::getKeyFactory()->add("Model", this);

  // The regular CModelEntity mechanism does not work since
  // CModel is created before mStateTemplate :(
  C_FLOAT64 InitialValue = *mpIValue;
  C_FLOAT64 Value = *mpValue;
  pdelete(mpIValue);
  pdelete(mpValue);
  mStateTemplate.add(this);
  *mpIValue = InitialValue;
  *mpValue = Value;

  mpIValueReference->setObjectName("Initial Time");
  mpValueReference->setObjectName("Time");

  mRate = 1.0;

  addObjectReference("Comments", *const_cast<std::string *>(&getNotes()));

  // These are broken since they contain pointers to values :(
  //  addVectorReference("Fluxes", mFluxes, CCopasiObject::ValueDbl);
  //  addVectorReference("Particle Fluxes", mParticleFluxes, CCopasiObject::ValueDbl);

  addMatrixReference("Stoichiometry", mStoi, CCopasiObject::ValueDbl);
  addMatrixReference("Reduced Model Stoichiometry", mRedStoi, CCopasiObject::ValueDbl);

  addMatrixReference("Link Matrix"   , mLView, CCopasiObject::ValueDbl);
  addObjectReference("Quantity Unit", mQuantityUnit);
  addObjectReference("Quantity Conversion Factor", mQuantity2NumberFactor, CCopasiObject::ValueDbl);
  addObjectReference("Avogadro Constant", mAvogadro, CCopasiObject::ValueDbl);

  mpStoiAnnotation = new CArrayAnnotation("Stoichiometry(ann)", this, new CCopasiMatrixInterface<CMatrix<C_FLOAT64> >(&mStoiReordered), true);
  mpStoiAnnotation->setDescription("Stoichiometry Matrix");
  mpStoiAnnotation->setMode(0, CArrayAnnotation::OBJECTS);
  mpStoiAnnotation->setDimensionDescription(0, "Species that are controlled by reactions");
  mpStoiAnnotation->setMode(1, CArrayAnnotation::VECTOR_ON_THE_FLY);
  mpStoiAnnotation->setDimensionDescription(1, "Reactions");
  mpStoiAnnotation->setCopasiVector(1, &mSteps);

  mpRedStoiAnnotation = new CArrayAnnotation("Reduced stoichiometry(ann)", this, new CCopasiMatrixInterface<CMatrix<C_FLOAT64> >(&mRedStoi), true);
  mpRedStoiAnnotation->setDescription("Reduced stoichiometry Matrix");
  mpRedStoiAnnotation->setMode(0, CArrayAnnotation::OBJECTS);
  mpRedStoiAnnotation->setDimensionDescription(0, "Species (reduced system)");
  mpRedStoiAnnotation->setMode(1, CArrayAnnotation::VECTOR_ON_THE_FLY);
  mpRedStoiAnnotation->setDimensionDescription(1, "Reactions");
  mpRedStoiAnnotation->setCopasiVector(1, &mSteps);

  mpLinkMatrixAnnotation = new CArrayAnnotation("Link matrix(ann)", this, new CCopasiMatrixInterface<CLinkMatrixView>(&mLView), true);
  mpLinkMatrixAnnotation->setDescription("Link matrix");
  mpLinkMatrixAnnotation->setMode(0, CArrayAnnotation::OBJECTS);
  mpLinkMatrixAnnotation->setDimensionDescription(0, "Species that are controlled by reactions (full system)");
  mpLinkMatrixAnnotation->setMode(1, CArrayAnnotation::OBJECTS);
  mpLinkMatrixAnnotation->setDimensionDescription(1, "Species (reduced system)");

  mpMathModel = new CMathModel(this);
}

bool CModel::hasReversibleReaction() const
{
  unsigned C_INT32 i, imax = mSteps.size();

  for (i = 0; i < imax; ++i) if (mSteps[i]->isReversible()) return true;

  return false;
}

//**********************************************************************
//                   CLinkMatrixView
//**********************************************************************

const CModel::CLinkMatrixView::elementType CModel::CLinkMatrixView::mZero = 0.0;
const CModel::CLinkMatrixView::elementType CModel::CLinkMatrixView::mUnit = 1.0;

CModel::CLinkMatrixView::CLinkMatrixView(const CMatrix< C_FLOAT64 > & A,
    const unsigned C_INT32 & numIndependent):
    mA(A),
    mNumIndependent(numIndependent)
{CONSTRUCTOR_TRACE;}

CModel::CLinkMatrixView::~CLinkMatrixView()
{DESTRUCTOR_TRACE;}

CModel::CLinkMatrixView &
CModel::CLinkMatrixView::operator = (const CModel::CLinkMatrixView & rhs)
{
  const_cast< CMatrix< C_FLOAT64 > &>(mA) = rhs.mA;
  const_cast< unsigned C_INT32 & >(mNumIndependent) = rhs.mNumIndependent;

  return *this;
}

unsigned C_INT32 CModel::CLinkMatrixView::numRows() const
{return mNumIndependent + mA.numRows();}

unsigned C_INT32 CModel::CLinkMatrixView::numCols() const
{return mA.numCols();}

std::ostream &operator<<(std::ostream &os,
                         const CModel::CLinkMatrixView & A)
{
  unsigned C_INT32 i, imax = A.numRows();
  unsigned C_INT32 j, jmax = A.numCols();
  os << "Matrix(" << imax << "x" << jmax << ")" << std::endl;

  for (i = 0; i < imax; i++)
    {
      for (j = 0; j < jmax; j++)
        os << "\t" << A(i, j);

      os << std::endl;
    }

  return os;
}

std::string CModel::suitableForStochasticSimulation() const
{
  unsigned C_INT32 i, reactSize = mSteps.size();
  C_INT32 multInt;
  unsigned C_INT32 j;
  C_FLOAT64 multFloat;
  //  C_INT32 metabSize = mMetabolites->size();

  for (i = 0; i < reactSize; i++) // for every reaction
    {
      // TEST getCompartmentNumber() == 1
      //if (mSteps[i]->getCompartmentNumber() != 1) return - 1;

      // TEST isReversible() == 0
      if (mSteps[i]->isReversible() != 0)
        return "At least one reaction is reversible. That means stochastic simulation is not possible. \nYou can use \"Tools|Convert to irreversible\" which will split the reversible reactions \n into two irreversible reactions. However you should check the kinetics afterwards.";

      // TEST integer stoichiometry
      // Iterate through each the metabolites
      // Juergen: the number of rows of mStoi equals the number of non-fixed metabs!
      //  for (j=0; i<metabSize; j++)
      for (j = 0; j < mStoi.numRows(); j++)
        {
          multFloat = mStoi(j, i);
          multInt = static_cast<C_INT32>(floor(multFloat + 0.5)); // +0.5 to get a rounding out of the static_cast to int!

          if ((multFloat - multInt) > 0.01)
            return "Not all stoichiometries are integer numbers. \nThat means that discrete simulation is not possible.";
        }
    }

  for (i = 0; i < mMetabolites.size(); ++i)
    {
      if (mMetabolites[i]->getInitialValue() > LLONG_MAX)
        return "At least one particle number in the initial state is too big.";
    }

  return ""; // Model is appropriate for hybrid simulation
}

#ifdef COPASI_DEBUG
void CModel::check() const
{}
#endif

void CModel::buildLinkZero()
{
  // Prior to a call to buildLinkZero the stoichiometry matrix mStoi must
  // have been constructed.

  mRedStoi = mStoi;
  
  C_INT NumReactions = mRedStoi.numCols();
  C_INT NumSpecies = mRedStoi.numRows();
  C_INT LDA = std::max<C_INT>(1, NumReactions);

  CVector< C_INT > JPVT(NumSpecies);
  JPVT = 0;

  C_INT32 Dim = std::min(NumReactions, NumSpecies);

  if (Dim == 0)
    {
      C_INT32 i;
      mRowLU.resize(NumSpecies);

      for (i = 0; i < NumSpecies; i++)
        mRowLU[i] = i;

      mNumMetabolitesReactionIndependent = 0;
      mL.resize(NumSpecies - 0, 0);

      return;
    }

  CVector< C_FLOAT64 > TAU(Dim);

  CVector< C_FLOAT64 > WORK(1);
  C_INT LWORK = -1;
  C_INT INFO;

  // QR factorization of the stoichiometry matrix
  /*
   *  -- LAPACK routine (version 3.0) --
   *     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
   *     Courant Institute, Argonne National Lab, and Rice University
   *     June 30, 1999
   *
   *  Purpose
   *  =======
   *
   *  DGEQP3 computes a QR factorization with column pivoting of a
   *  matrix A:  A*P = Q*R  using Level 3 BLAS.
   *
   *  Arguments
   *  =========
   *
   *  M       (input) INTEGER
   *          The number of rows of the matrix A. M >= 0.
   *
   *  N       (input) INTEGER
   *          The number of columns of the matrix A.  N >= 0.
   *
   *  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
   *          On entry, the M-by-N matrix A.
   *          On exit, the upper triangle of the array contains the
   *          min(M,N)-by-N upper trapezoidal matrix R; the elements below
   *          the diagonal, together with the array TAU, represent the
   *          orthogonal matrix Q as a product of min(M,N) elementary
   *          reflectors.
   *
   *  LDA     (input) INTEGER
   *          The leading dimension of the array A. LDA >= max(1,M).
   *
   *  JPVT    (input/output) INTEGER array, dimension (N)
   *          On entry, if JPVT(J).ne.0, the J-th column of A is permuted
   *          to the front of A*P (a leading column); if JPVT(J)=0,
   *          the J-th column of A is a free column.
   *          On exit, if JPVT(J)=K, then the J-th column of A*P was the
   *          the K-th column of A.
   *
   *  TAU     (output) DOUBLE PRECISION array, dimension (min(M,N))
   *          The scalar factors of the elementary reflectors.
   *
   *  WORK    (workspace/output) DOUBLE PRECISION array, dimension (LWORK)
   *          On exit, if INFO=0, WORK(1) returns the optimal LWORK.
   *
   *  LWORK   (input) INTEGER
   *          The dimension of the array WORK. LWORK >= 3*N+1.
   *          For optimal performance LWORK >= 2*N+(N+1)*NB, where NB
   *          is the optimal blocksize.
   *
   *          If LWORK = -1, then a workspace query is assumed; the routine
   *          only calculates the optimal size of the WORK array, returns
   *          this value as the first entry of the WORK array, and no error
   *          message related to LWORK is issued by XERBLA.
   *
   *  INFO    (output) INTEGER
   *          = 0: successful exit.
   *          < 0: if INFO = -i, the i-th argument had an illegal value.
   *
   *  Further Details
   *  ===============
   *
   *  The matrix Q is represented as a product of elementary reflectors
   *
   *     Q = H(1) H(2) . . . H(k), where k = min(m,n).
   *
   *  Each H(i) has the form
   *
   *     H(i) = I - tau * v * v'
   *
   *  where tau is a real/complex scalar, and v is a real/complex vector
   *  with v(1:i-1) = 0 and v(i) = 1; v(i+1:m) is stored on exit in
   *  A(i+1:m,i), and tau in TAU(i).
   *
   *  Based on contributions by
   *    G. Quintana-Orti, Depto. de Informatica, Universidad Jaime I, Spain
   *    X. Sun, Computer Science Dept., Duke University, USA
   *
   */

#ifdef DEBUG_MATRIX
  DebugFile << CTransposeView< CMatrix< C_FLOAT64 > >(mRedStoi) << std::endl;
#endif

  dgeqp3_(&NumReactions, &NumSpecies, mRedStoi.array(), &LDA,
          JPVT.array(), TAU.array(), WORK.array(), &LWORK, &INFO);

  if (INFO < 0) fatalError();

  LWORK = (C_INT) WORK[0];
  WORK.resize(LWORK);
  
  dgeqp3_(&NumReactions, &NumSpecies, mRedStoi.array(), &LDA,
          JPVT.array(), TAU.array(), WORK.array(), &LWORK, &INFO);

  if (INFO < 0) fatalError();

  C_INT32 i;
  mRowLU.resize(NumSpecies);

  for (i = 0; i < NumSpecies; i++)
    mRowLU[i] = JPVT[i] - 1;

#ifdef DEBUG_MATRIX
  DebugFile << "QR Factorization:" << std::endl;
  DebugFile << "Row permutation:\t" << mRowLU << std::endl;
  DebugFile << CTransposeView< CMatrix< C_FLOAT64 > >(mRedStoi) << std::endl;
#endif

  C_INT independent = 0;

  while (independent < Dim &&
         fabs(mRedStoi(independent, independent)) > 100.0 * std::numeric_limits< C_FLOAT64 >::epsilon()) independent++;

  // Resize mL
  mNumMetabolitesReactionIndependent = independent;
  mL.resize(NumSpecies - independent, independent);

  if (NumSpecies == independent || independent == 0) return;

  /* to take care of differences between fortran's and c's memory  access,
     we need to take the transpose, i.e.,the upper triangular */
  char cL = 'U';
  char cU = 'N'; /* values in the diagonal of R */
  // Calculate Row Echelon form of R.
  // First invert R_1,1
  /* int dtrtri_(char *uplo,
   *             char *diag,
   *             integer *n,
   *             doublereal * A,
   *             integer *lda,
   *             integer *info);
   *  -- LAPACK routine (version 3.0) --
   *     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
   *     Courant Institute, Argonne National Lab, and Rice University
   *     March 31, 1993
   *
   *  Purpose
   *  =======
   *
   *  DTRTRI computes the inverse of a real upper or lower triangular
   *  matrix A.
   *
   *  This is the Level 3 BLAS version of the algorithm.
   *
   *  Arguments
   *  =========
   *
   *  uplo    (input) CHARACTER*1
   *          = 'U':  A is upper triangular;
   *          = 'L':  A is lower triangular.
   *
   *  diag    (input) CHARACTER*1
   *          = 'N':  A is non-unit triangular;
   *          = 'U':  A is unit triangular.
   *
   *  n       (input) INTEGERstd::cout << "bool \n"; 
   *          The order of the matrix A.  n >= 0.
   *
   *  A       (input/output) DOUBLE PRECISION array, dimension (lda,n)
   *          On entry, the triangular matrix A.  If uplo = 'U', the
   *          leading n-by-n upper triangular part of the array A contains
   *          the upper triangular matrix, and the strictly lower
   *          triangular part of A is not referenced.  If uplo = 'L', the
   *          leading n-by-n lower triangular part of the array A contains
   *          the lower triangular matrix, and the strictly upper
   *          triangular part of A is not referenced.  If diag = 'U', the
   *          diagonal elements of A are also not referenced and are
   *          assumed to be 1.
   *          On exit, the (triangular) inverse of the original matrix, in
   *          the same storage format.
   *
   *  lda     (input) INTEGER
   *          The leading dimension of the array A.  lda >= max(1,n).
   *
   *  info    (output) INTEGER
   *          = 0: successful exit
   *          < 0: if info = -i, the i-th argument had an illegal value
   *          > 0: if info = i, A(i,i) is exactly zero.  The triangular
   *               matrix is singular and its inverse can not be computed.
   */
  dtrtri_(&cL, &cU, &independent, mRedStoi.array(), &LDA, &INFO);

  if (INFO < 0) fatalError();

#ifdef DEBUG_MATRIX
  DebugFile << "Invert R_1,1:" << std::endl;
  DebugFile << CTransposeView< CMatrix< C_FLOAT64 > >(mRedStoi) << std::endl;
#endif

  C_INT32 j, k;

  // Compute Link_0 = inverse(R_1,1) * R_1,2
  // :TODO: Use dgemm
  C_FLOAT64 * pTmp1 = &mL(0, 0);
  C_FLOAT64 * pTmp2;
  C_FLOAT64 * pTmp3;

  for (j = 0; j < NumSpecies - independent; j++)
    for (i = 0; i < independent; i++, pTmp1++)
      {
        pTmp2 = &mRedStoi(j + independent, i);
        pTmp3 = &mRedStoi(i, i);

        // assert(&mL(j, i) == pTmp3);
        *pTmp1 = 0.0;

        for (k = i; k < independent; k++, pTmp2++, pTmp3 += NumReactions)
          {
            // assert(&mRedStoi(j + independent, k) == pTmp2);
            // assert(&mRedStoi(k, i) == pTmp3);

            *pTmp1 += *pTmp3 * *pTmp2;
          }

        if (fabs(*pTmp1) < 100.0 * std::numeric_limits< C_FLOAT64 >::epsilon()) *pTmp1 = 0.0;
      }

#ifdef DEBUG_MATRIX
  DebugFile << "Link Zero Matrix:" << std::endl;
  DebugFile << mL << std::endl;
#endif // DEBUG_MATRIX

  return;
}

const bool & CModel::isAutonomous() const
{return mIsAutonomous;}

void CModel::determineIsAutonomous()
{
  mIsAutonomous = true;

  // If the model is not empty we check whether anything depends on time
  if (mCompartments.size() != 0 ||
      mValues.size() != 0)
    {
      std::set< const CCopasiObject * > TimeDependent;

      appendDependentReactions(getDeletedObjects(), TimeDependent);
      appendDependentMetabolites(getDeletedObjects(), TimeDependent);
      appendDependentCompartments(getDeletedObjects(), TimeDependent);
      appendDependentModelValues(getDeletedObjects(), TimeDependent);
      appendDependentEvents(getDeletedObjects(), TimeDependent);

      mIsAutonomous = (TimeDependent.begin() == TimeDependent.end());
    }

  // An autonomous models always start simulation at T = 0
  if (mIsAutonomous)
    setInitialValue(0.0);
}

bool CModel::compileEvents()
{
  bool success = true;

  std::vector< CCopasiContainer * > ListOfContainer;

  CCopasiVectorN< CEvent >::iterator it = mEvents.begin();
  CCopasiVectorN< CEvent >::iterator end = mEvents.end();

  for (; it != end; ++ it)
    {
      success &= (*it)->compile(ListOfContainer);
    }

  return success;
}

const std::vector< Refresh * > & CModel::getListOfInitialRefreshes() const
{return mInitialRefreshes;}

const std::vector< Refresh * > & CModel::getListOfSimulatedRefreshes() const
{return mSimulatedRefreshes;}

const std::vector< Refresh * > & CModel::getListOfConstantRefreshes() const
{return mApplyInitialValuesRefreshes;}

const std::vector< Refresh * > & CModel::getListOfNonSimulatedRefreshes() const
{return mNonSimulatedRefreshes;}

std::vector< Refresh * >
CModel::buildInitialRefreshSequence(std::set< const CCopasiObject * > & changedObjects)
{
  // First we remove all objects which are of type assignment from the changed objects
  // since this may not be changed as they are under control of the assignment.
  std::set< const CCopasiObject * >::iterator itSet;
  std::set< const CCopasiObject * >::iterator endSet;
  std::set< const CCopasiObject * > Objects;
  std::set< const CCopasiObject * > Context;

  CModelEntity **ppEntity;
  CModelEntity **ppEntityEnd = mStateTemplate.endFixed();

  const CModelEntity * pEntity;
  CMetab * pMetab;

  // If the changed objects are empty we assume that all changeable objects have been changed
  if (changedObjects.size() == 0)
    {
      // The objects which are changed are all initial values of of all model entities including
      // fixed and unused once. Additionally, all kinetic parameters are possibly changed.
      // This is basically all the parameters in the parameter overview whose value is editable.

      // :TODO: Theoretically, it is possible that also task parameters influence the initial
      // state of a model but that is currently not handled.

      // The initial values of the model entities
      ppEntity = mStateTemplate.beginIndependent() - 1; // Offset for time

      for (; ppEntity != ppEntityEnd; ++ppEntity)
        {
          // If we have an initial expression we have no initial values
          if (((*ppEntity)->getInitialExpression() != "" ||
               (*ppEntity)->getStatus() == ASSIGNMENT) &&
              (*ppEntity)->getInitialValueReference()->getDirectDependencies().size() > 0)
            continue;

          // Metabolites have two initial values
          if ((pMetab = dynamic_cast< CMetab * >(*ppEntity)) != NULL)
            {
              // The concentration is assumed to be fix accept when this would lead to circular dependencies,
              // for the parent's compartment's initial volume.
              if (pMetab->isInitialConcentrationChangeAllowed())
                changedObjects.insert(pMetab->getInitialConcentrationReference());
              else
                changedObjects.insert(pMetab->getInitialValueReference());
            }
          else
            changedObjects.insert((*ppEntity)->getInitialValueReference());
        }

      // The reaction parameters
      CCopasiVector< CReaction >::const_iterator itReaction = mSteps.begin();
      CCopasiVector< CReaction >::const_iterator endReaction = mSteps.end();
      unsigned C_INT32 i, imax;

      for (; itReaction != endReaction; ++itReaction)
        {
          const CCopasiParameterGroup & Group = (*itReaction)->getParameters();

          for (i = 0, imax = Group.size(); i < imax; i++)
            changedObjects.insert(Group.getParameter(i)->getObject(CCopasiObjectName("Reference=Value")));
        }

      // Fix for Issue 1170: We need to add elements of the stoichiometry, reduced stoichiometry,
      // and link matrices.
      if (mpStoiAnnotation != NULL)
        mpStoiAnnotation->appendElementReferences(changedObjects);

      if (mpRedStoiAnnotation != NULL)
        mpRedStoiAnnotation->appendElementReferences(changedObjects);

      if (mpLinkMatrixAnnotation != NULL)
        mpLinkMatrixAnnotation->appendElementReferences(changedObjects);
    }
  else
    {
      // Remove all objects with initial assignments
      itSet = changedObjects.begin();
      endSet = changedObjects.end();

      for (; itSet != endSet; ++itSet)
        if ((pEntity = dynamic_cast< const CModelEntity * >((*itSet)->getObjectParent())) != NULL &&
            (pEntity->getInitialExpression() != "" ||
             pEntity->getStatus() == ASSIGNMENT) &&
            pEntity->getInitialValueReference()->getDirectDependencies().size() > 0)
          Objects.insert(*itSet);

      for (itSet = Objects.begin(), endSet = Objects.end(); itSet != endSet; ++itSet)
        changedObjects.erase(*itSet);
    }

  // We need to add all initial values which are dynamically calculated.
  // These are initial assignments and either a metabolite's initial particle number
  // or concentration.
  ppEntity = mStateTemplate.beginIndependent() - 1; // Offset for time
  Objects.clear();

  for (; ppEntity != ppEntityEnd; ++ppEntity)
    {
      if (((*ppEntity)->getInitialExpression() != "" ||
           (*ppEntity)->getStatus() == ASSIGNMENT))
        {
          Objects.insert((*ppEntity)->getInitialValueReference());
          continue;
        }

      // For metabolites we need to add the initial concentration or the initial
      // particle number..
      if ((pMetab = dynamic_cast< CMetab * >(*ppEntity)) != NULL)
        {
          if (changedObjects.count(pMetab->getInitialConcentrationReference()) != 0)
            {
              pMetab->compileInitialValueDependencies(false);
              Objects.insert(pMetab->getInitialValueReference());
            }
          else
            {
              pMetab->compileInitialValueDependencies(true);
              Objects.insert(pMetab->getInitialConcentrationReference());
            }
        }
    }

  // We need to add the total particle number of moieties.
  CCopasiVector< CMoiety >::iterator itMoiety = mMoieties.begin();
  CCopasiVector< CMoiety >::iterator endMoiety = mMoieties.end();

  for (; itMoiety != endMoiety; ++itMoiety)
    Objects.insert((*itMoiety)->getInitialValueReference());

  std::set< const CCopasiObject * > DependencySet;
  std::set< const CCopasiObject * > VerifiedSet;
  std::pair<std::set< const CCopasiObject * >::iterator, bool> InsertedObject;

  assert(Objects.count(NULL) == 0);

  // Check whether we have any circular dependencies
  for (itSet = Objects.begin(), endSet = Objects.end(); itSet != endSet; ++itSet)
    if ((*itSet)->hasCircularDependencies(DependencySet, VerifiedSet, Context))
      CCopasiMessage(CCopasiMessage::EXCEPTION, MCObject + 1, (*itSet)->getCN().c_str());

  // Build the complete set of dependencies
  for (itSet = Objects.begin(); itSet != endSet; ++itSet)
    {
      // At least the object itself needs to be up to date.
      InsertedObject = DependencySet.insert(*itSet);

      // Add all its dependencies
      if (InsertedObject.second)
        (*itSet)->getAllDependencies(DependencySet, Context);
    }

  // Remove all objects which do not depend on the changed objects, or do not have a
  // refresh method.
  Objects.clear();

  // We now check which objects we need to refresh
  for (itSet = DependencySet.begin(), endSet = DependencySet.end(); itSet != endSet; ++itSet)
    {
      // No refresh method
      if ((*itSet)->getRefresh() == NULL)
        Objects.insert(*itSet);
      // Is a changed object
      else if (changedObjects.count(*itSet) != 0)
        Objects.insert(*itSet);
      // Not dependent on the changed objects.
      else if (!(*itSet)->dependsOn(changedObjects))
        Objects.insert(*itSet);
    }

  for (itSet = Objects.begin(), endSet = Objects.end(); itSet != endSet; ++itSet)
    DependencySet.erase(*itSet);

  // Create a properly sorted list.
  std::list< const CCopasiObject * > SortedList =
    sortObjectsByDependency(DependencySet.begin(), DependencySet.end(), Context);

  std::list< const CCopasiObject * >::iterator itList;
  std::list< const CCopasiObject * >::iterator endList;

  // Build the vector of pointers to refresh methods
  Refresh * pRefresh;

  std::vector< Refresh * > UpdateVector;
  std::vector< Refresh * >::const_iterator itUpdate;
  std::vector< Refresh * >::const_iterator endUpdate;

  itList = SortedList.begin();
  endList = SortedList.end();

  for (; itList != endList; ++itList)
    {
      pRefresh = (*itList)->getRefresh();
      itUpdate = UpdateVector.begin();
      endUpdate = UpdateVector.end();

      while (itUpdate != endUpdate && !(*itUpdate)->isEqual(pRefresh)) ++itUpdate;

      if (itUpdate == endUpdate)
        UpdateVector.push_back(pRefresh);
    }

  return UpdateVector;
}

CVector< C_FLOAT64 > CModel::initializeAtolVector(const C_FLOAT64 & atol, const bool & reducedModel) const
{
  CVector< C_FLOAT64 > Atol;

  if (reducedModel)
    Atol.resize(mStateTemplate.getNumIndependent());
  else
    Atol.resize(mStateTemplate.getNumIndependent() + getNumDependentReactionMetabs());

  C_FLOAT64 * pAtol = Atol.array();
  C_FLOAT64 * pEnd = pAtol + Atol.size();

  C_FLOAT64 InitialValue;
  C_FLOAT64 Limit;

  CModelEntity *const* ppEntity = getStateTemplate().beginIndependent();
  const CMetab * pMetab;

  for (; pAtol != pEnd; ++pAtol, ++ppEntity)
    {
      *pAtol = atol;

      InitialValue = fabs((*ppEntity)->getInitialValue());

      if ((pMetab = dynamic_cast< const CMetab * >(*ppEntity)) != NULL)
        {
          Limit =
            fabs(pMetab->getCompartment()->getInitialValue()) * mQuantity2NumberFactor;

          if (InitialValue != 0.0)
            *pAtol *= std::min(Limit, InitialValue);
          else
            *pAtol *= std::max(1.0, Limit);
        }
      else if (InitialValue != 0.0)
        *pAtol *= std::min(1.0, InitialValue);
    }

  return Atol;
}

#include "utilities/CDimension.h"

std::string CModel::printParameterOverview()
{
  std::ostringstream oss;
  CModel* model = this;

  oss << "Initial time: " << model->getInitialTime() << " " << model->getTimeUnitName() << std::endl;

  oss << std::endl;

  unsigned C_INT32 i, imax, j, jmax;

  //Compartments
  const CCopasiVector< CCompartment > & comps = model->getCompartments();
  imax = comps.size();

  if (imax)
    {
      oss << "Initial volumes:\n\n";

      for (i = 0; i < imax; ++i)
        oss << comps[i]->getObjectName() << " \t" << comps[i]->getInitialValue()
        << " " << model->getVolumeUnitsDisplayString() << "\n";

      oss << "\n";
    }

  //Species
  const CCopasiVector< CMetab > & metabs = model->getMetabolites();
  imax = metabs.size();

  if (imax)
    {
      oss << "Initial concentrations:\n\n";

      for (i = 0; i < imax; ++i)
        oss << CMetabNameInterface::getDisplayName(model, *metabs[i]) << " \t"
        << metabs[i]->getInitialConcentration() << " "
        << model->getConcentrationUnitsDisplayString() << "\n";

      oss << "\n";
    }

  //global Parameters
  const CCopasiVector< CModelValue > & params = model->getModelValues();
  imax = params.size();

  if (imax)
    {
      oss << "Initial values of global quantities:\n\n";

      for (i = 0; i < imax; ++i)
        oss << params[i]->getObjectName() << " \t"
        << params[i]->getInitialValue() << "\n";

      oss << "\n";
    }

  //Reactions
  const CCopasiVector< CReaction > & reacs = model->getReactions();
  imax = reacs.size();

  if (imax)
    {
      oss << "Reaction parameters:\n\n";
      CReaction* reac;

      for (i = 0; i < imax; ++i)
        {
          reac = reacs[i];
          oss << reac->getObjectName() << "\n";

          //calculate units
          CFindDimensions units(reac->getFunction(), getQuantityUnitEnum() == dimensionlessQuantity,
                                getVolumeUnitEnum() == dimensionlessVolume,
                                getTimeUnitEnum() == dimensionlessTime,
                                getAreaUnitEnum() == dimensionlessArea,
                                getLengthUnitEnum() == dimensionlessLength);
          units.setUseHeuristics(true);
          units.setChemicalEquation(&reac->getChemEq());
          units.findDimensions(reac->getCompartmentNumber() > 1);

          const CFunctionParameters & params = reac->getFunctionParameters();
          jmax = params.size();

          for (j = 0; j < jmax; ++j)
            if (params[j]->getUsage() == CFunctionParameter::PARAMETER)
              {
                CCopasiObject * obj = CCopasiRootContainer::getKeyFactory()->get(reac->getParameterMappings()[j][0]);

                if (!obj) continue;

                if (reac->isLocalParameter(j))
                  {
                    CCopasiParameter * par = dynamic_cast<CCopasiParameter*>(obj); //must be a CCopasiParameter

                    if (!par) continue; //or rather fatal error?

                    oss << "    " << params[j]->getObjectName() << " \t"
                    << *par->getValue().pDOUBLE << " "
                    << units.getDimensions()[j].getDisplayString(getObjectDataModel()) << "\n";
                  }
                else
                  {
                    CModelValue * par = dynamic_cast<CModelValue*>(obj); //must be a CModelValue

                    if (!par) continue; //or rather fatal error?

                    oss << "    " << params[j]->getObjectName() << " \t"
                    << "-> " + par->getObjectName()
                    << " (" << units.getDimensions()[j].getDisplayString(getObjectDataModel()) << ")\n";
                  }
              }

          oss << "\n";
        }
    }

  return oss.str();
}

std::string CModel::getTimeUnitsDisplayString() const
{
  if (mTimeUnit == dimensionlessTime)
    return "";

  return TimeUnitNames[mTimeUnit];
}

std::string CModel::getFrequencyUnitsDisplayString() const
{
  if (mTimeUnit == dimensionlessTime)
    return "";

  return std::string("1/") + TimeUnitNames[mTimeUnit];
}

std::string CModel::getVolumeUnitsDisplayString() const
{
  if (mVolumeUnit == dimensionlessVolume)
    return "";

  return VolumeUnitNames[mVolumeUnit];
}

std::string CModel::getVolumeRateUnitsDisplayString() const
{
  if (mVolumeUnit == dimensionlessVolume)
    {
      if (mTimeUnit == dimensionlessTime)
        return "";

      return std::string("1/") + TimeUnitNames[mTimeUnit];
    }

  if (mTimeUnit == dimensionlessTime)
    return VolumeUnitNames[mVolumeUnit];

  return std::string(VolumeUnitNames[mVolumeUnit]) + "/" + TimeUnitNames[mTimeUnit];
}

std::string CModel::getConcentrationUnitsDisplayString() const
{
  std::string Units;

  if (mQuantityUnit == dimensionlessQuantity)
    {
      if (mVolumeUnit == dimensionlessVolume)
        return "";

      return std::string("1/") + VolumeUnitNames[mVolumeUnit];
    }

  Units = QuantityUnitNames[mQuantityUnit];

  if (mVolumeUnit == dimensionlessVolume)
    return Units;

  return Units + "/" + VolumeUnitNames[mVolumeUnit];
}

std::string CModel::getConcentrationRateUnitsDisplayString() const
{
  std::string Units;

  if (mQuantityUnit == dimensionlessQuantity)
    {
      Units = "1";

      if (mVolumeUnit == dimensionlessVolume)
        {
          if (mTimeUnit == dimensionlessTime)
            return "";

          return Units + "/" + TimeUnitNames[mTimeUnit];
        }
      else
        {
          if (mTimeUnit == dimensionlessTime)
            return Units + "/" + VolumeUnitNames[mVolumeUnit];

          return Units + "/(" + VolumeUnitNames[mVolumeUnit] + "*" + TimeUnitNames[mTimeUnit] + ")";
        }
    }

  Units = QuantityUnitNames[mQuantityUnit];

  if (mVolumeUnit == dimensionlessVolume)
    {
      if (mTimeUnit == dimensionlessTime)
        return Units;

      return Units + "/" + TimeUnitNames[mTimeUnit];
    }

  if (mTimeUnit == dimensionlessTime)
    return Units + "/" + VolumeUnitNames[mVolumeUnit];

  return Units + "/(" + VolumeUnitNames[mVolumeUnit] + "*" + TimeUnitNames[mTimeUnit] + ")";
}

std::string CModel::getQuantityRateUnitsDisplayString() const
{
  std::string Units;

  if (mQuantityUnit == dimensionlessQuantity)
    {
      if (mTimeUnit == dimensionlessTime)
        return "";

      return std::string("1/") + TimeUnitNames[mTimeUnit];
    }

  Units = QuantityUnitNames[mQuantityUnit];

  if (mTimeUnit == dimensionlessTime)
    return Units;

  return Units + "/" + TimeUnitNames[mTimeUnit];
}

/****** Below will be removed when the math model completed ******/

void CModel::evaluateRoots(CVectorCore< C_FLOAT64 > & rootValues,
                           const bool & ignoreDiscrete)
{
  return mpMathModel->evaluateRoots(rootValues, ignoreDiscrete);
}

bool CModel::processQueue(const C_FLOAT64 & time,
                          const bool & equality,
                          CProcessQueue::resolveSimultaneousAssignments pResolveSimultaneousAssignments)
{
  return mpMathModel->processQueue(time, equality, pResolveSimultaneousAssignments);
}

void CModel::processRoots(const C_FLOAT64 & time,
                          const bool & equality,
                          const CVector< C_INT > & roots)
{
  return mpMathModel->processRoots(time, equality, roots);
}

const C_FLOAT64 & CModel::getProcessQueueExecutionTime() const
{
  return mpMathModel->getProcessQueueExecutionTime();
}

size_t CModel::getNumRoots() const
{
  return mpMathModel->getNumRoots();
}

void CModel::calculateRootDerivatives(CVector< C_FLOAT64 > & rootDerivatives)
{
  return mpMathModel->calculateRootDerivatives(rootDerivatives);
}

const CVector< CMathTrigger::CRootFinder * > & CModel::getRootFinders() const
{
  return mpMathModel->getRootFinders();
}

