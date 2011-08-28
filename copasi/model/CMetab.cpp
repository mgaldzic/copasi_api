// Begin CVS Header
//   $Source: /fs/turing/cvs/copasi_dev/copasi/model/CMetab.cpp,v $
//   $Revision: 1.153 $
//   $Name: Build-33 $
//   $Author: shoops $
//   $Date: 2010/08/02 19:41:21 $
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

#include <iostream>
#include <string>
#include <vector>
#include <limits>

#include "copasi.h"

#include "CopasiDataModel/CCopasiDataModel.h"
#include "utilities/utility.h"
#include "report/CKeyFactory.h"
#include "CCompartment.h"
#include "CModel.h"
#include "CMetab.h"
#include "CMetabNameInterface.h"
#include "function/CExpression.h"
#include "report/CCopasiRootContainer.h"

#define METAB_MOIETY 7

//static
C_FLOAT64 CMetab::convertToNumber(const C_FLOAT64 & concentration,
                                  const CCompartment & compartment,
                                  const CModel & model)
{return concentration * compartment.getInitialValue() * model.getQuantity2NumberFactor();}

//static
C_FLOAT64 CMetab::convertToConcentration(const C_FLOAT64 & number,
    const CCompartment & compartment,
    const CModel & model)
{return number / compartment.getInitialValue() * model.getNumber2QuantityFactor();}

CMetab::CMetab(const std::string & name,
               const CCopasiContainer * pParent):
    CModelEntity(name, pParent, "Metabolite",
                 CCopasiObject::NonUniqueName),
    mConc(std::numeric_limits<C_FLOAT64>::quiet_NaN()),
    mIConc(0.0),
    mConcRate(0.0),
    mTT(0.0),
    mpCompartment(NULL),
    mpMoiety(NULL)
{
  mKey = CCopasiRootContainer::getKeyFactory()->add("Metabolite", this);
  initObjects();

  setStatus(REACTIONS);

  if (getObjectParent())
    {
      initCompartment(NULL);

      setInitialConcentration(1.0);
      setConcentration(1.0);
    }

  CONSTRUCTOR_TRACE;
}

CMetab::CMetab(const CMetab & src,
               const CCopasiContainer * pParent):
    CModelEntity(src, pParent),
    mConc(src.mConc),
    mIConc(src.mIConc),
    mConcRate(src.mConcRate),
    mTT(src.mTT),
    mpCompartment(NULL),
    mpMoiety(src.mpMoiety)
{
  mKey = CCopasiRootContainer::getKeyFactory()->add("Metabolite", this);

  initObjects();

  initCompartment(src.mpCompartment);
  CONSTRUCTOR_TRACE;
}

CMetab &CMetab::operator=(const CMetabOld & RHS)
{
  setObjectName(RHS.getObjectName());

  setStatus(RHS.mStatus);

  // We need to set the initial particle number since that is the expected for the initial state
  C_FLOAT64 InitialParticleNumber =
    CMetab::convertToNumber(RHS.mIConc, *mpCompartment, *mpModel);
  setInitialValue(InitialParticleNumber);

  mRate = 0.0;
  mConcRate = 0.0;
  mTT = 0.0;

  return *this;  // Assignment operator returns left side.
}

CMetab::~CMetab()
{
  CCopasiRootContainer::getKeyFactory()->remove(mKey);
  DESTRUCTOR_TRACE;
}

// virtual
std::string CMetab::getChildObjectUnits(const CCopasiObject * pObject) const
{
  if (mpModel == NULL) return "";

  const std::string & Name = pObject->getObjectName();

  if (Name == "InitialParticleNumber" ||
      Name == "ParticleNumber")
    return "";
  else if (Name == "ParticleNumberRate")
    return mpModel->getFrequencyUnitsDisplayString();
  else if (Name == "InitialConcentration" ||
           Name == "Concentration")
    return mpModel->getConcentrationUnitsDisplayString();
  else if (Name == "Rate")
    return mpModel->getConcentrationRateUnitsDisplayString();
  else if (Name == "TransitionTime")
    return mpModel->getTimeUnitsDisplayString();

  return "";
}

void CMetab::cleanup() {}

void CMetab::initCompartment(const CCompartment * pCompartment)
{
  mpCompartment =
    dynamic_cast< const CCompartment * >(getObjectAncestor("Compartment"));

  if (!mpCompartment)
    mpCompartment = pCompartment;
}

const C_FLOAT64 & CMetab::getConcentration() const {return mConc;}

const C_FLOAT64 & CMetab::getInitialConcentration() const {return mIConc;}

const CCompartment * CMetab::getCompartment() const {return mpCompartment;}

const CModel * CMetab::getModel() const {return mpModel;}

void CMetab::setTransitionTime(const C_FLOAT64 & transitionTime)
{mTT = transitionTime;}

const C_FLOAT64 & CMetab::getTransitionTime() const {return mTT;}

bool CMetab::setObjectParent(const CCopasiContainer * pParent)
{
  CModelEntity::setObjectParent(pParent);
  initCompartment(NULL);

  Status CurrentStatus = getStatus();

  if (CurrentStatus != FIXED)
    setStatus(FIXED);
  else
    setStatus(REACTIONS);

  setStatus(CurrentStatus);

  return true;
}

// ***** set quantities ********

void CMetab::setConcentration(const C_FLOAT64 concentration)
{
  if (isFixed()) return;

  mConc = concentration;
}

void CMetab::setInitialConcentration(const C_FLOAT64 & initialConcentration)
{
  mIConc = initialConcentration;

  return;
}

void CMetab::refreshInitialValue()
{
  *mpIValue = mIConc * mpCompartment->getInitialValue() * mpModel->getQuantity2NumberFactor();
}

void CMetab::refreshInitialConcentration()
{
  if (mpInitialExpression != NULL &&
      mpInitialExpression->getInfix() != "")
    mIConc = mpInitialExpression->calcValue();
  else
    mIConc =
      *mpIValue / mpCompartment->getInitialValue() * mpModel->getNumber2QuantityFactor();

  if (isFixed()) mConc = mIConc;
}

void CMetab::refreshConcentration()
{
  mConc = *mpValue / mpCompartment->getValue() * mpModel->getNumber2QuantityFactor();
}

void CMetab::refreshNumber()
{
  *mpValue = mConc * mpCompartment->getValue() * mpModel->getQuantity2NumberFactor();
}

void CMetab::setStatus(const CModelEntity::Status & status)
{
  Status OldStatus = getStatus();

  CModelEntity::setStatus(status);

  if (status == OldStatus) return;

  std::set< const CCopasiObject * > Dependencies;

  const CCopasiObject * pVolumeReference = NULL;

  if (mpCompartment)
    pVolumeReference = mpCompartment->getObject(CCopasiObjectName("Reference=Volume"));

  switch (getStatus())
    {
      case FIXED:
        break;

      case ASSIGNMENT:
        Dependencies.insert(mpConcReference);

        if (pVolumeReference)
          Dependencies.insert(pVolumeReference);

        mpValueReference->setDirectDependencies(Dependencies);
        mpValueReference->setRefresh(this, &CMetab::refreshNumber);

        mpConcReference->setDirectDependencies(mpExpression->getDirectDependencies());
        mpConcReference->setRefresh(this, &CMetab::calculate);

        // The dependencies and refresh of the rate are correct (see CModelEntity::setStatus).

        mpConcRateReference->setDirectDependencies(mpRateReference->getDirectDependencies());
        mpConcRateReference->clearRefresh();
        break;

      case ODE:
        mpValueReference->setDirectDependencies(Dependencies);
        mpValueReference->clearRefresh();

        Dependencies.insert(mpValueReference);

        if (pVolumeReference)
          Dependencies.insert(pVolumeReference);

        mpConcReference->setDirectDependencies(Dependencies);
        mpConcReference->setRefresh(this, &CMetab::refreshConcentration);

        Dependencies.clear();
        Dependencies.insert(mpConcRateReference);

        if (pVolumeReference)
          Dependencies.insert(pVolumeReference);

        mpRateReference->setDirectDependencies(Dependencies);
        mpRateReference->setRefresh(this, &CMetab::refreshRate);

        mpConcRateReference->setDirectDependencies(mpExpression->getDirectDependencies());
        mpConcRateReference->setRefresh(this, &CMetab::calculate);

        break;

      case REACTIONS:
        mpValueReference->setDirectDependencies(Dependencies);
        mpValueReference->clearRefresh();

        if (pVolumeReference)
          Dependencies.insert(pVolumeReference);

        Dependencies.insert(mpValueReference);
        mpConcReference->setDirectDependencies(Dependencies);
        mpConcReference->setRefresh(this, &CMetab::refreshConcentration);

        mpRateReference->setRefresh(this, &CMetab::refreshRate);

        Dependencies.clear();

        if (pVolumeReference)
          Dependencies.insert(pVolumeReference);

        Dependencies.insert(mpRateReference);
        mpConcRateReference->setDirectDependencies(Dependencies);
        mpConcRateReference->setRefresh(this, &CMetab::refreshConcentrationRate);
        break;

      default:
        break;
    }

  if (mpModel && mpCompartment) refreshConcentration();
}

bool CMetab::compile()
{
  bool success = true;

  // We first clear all dependencies and refreshes
  // Particle Number
  mpValueReference->clearDirectDependencies();
  mpValueReference->clearRefresh();

  // Rate (particle number rate)
  mRateVector.clear();
  mpRateReference->clearDirectDependencies();
  mpRateReference->clearRefresh();

  // Concentration
  mpConcReference->clearDirectDependencies();
  mpConcReference->clearRefresh();

  // Concentration Rate
  mpConcRateReference->clearDirectDependencies();
  mpConcRateReference->clearRefresh();

  // Transition Time
  mpTTReference->clearDirectDependencies();
  mpTTReference->clearRefresh();

  // Prepare the compilation
  std::set<const CCopasiObject *> Dependencies;
  std::vector< CCopasiContainer * > listOfContainer;
  listOfContainer.push_back(getObjectAncestor("Model"));

  CCopasiDataModel* pDataModel = NULL;
  const CCopasiObject * pVolumeReference = NULL;

  if (mpCompartment)
    pVolumeReference = mpCompartment->getValueReference();

  // Compile the value (particle number)
  Dependencies.insert(mpConcReference);

  if (pVolumeReference)
    Dependencies.insert(pVolumeReference);

  // We no longer need to distinguish the cases since the reference are now context sensitive
  mpValueReference->setDirectDependencies(Dependencies);
  mpValueReference->setRefresh(this, &CMetab::refreshNumber);

  Dependencies.clear();

  // Compiling of the rest.
  switch (getStatus())
    {
      case FIXED:
        // Concentration
        Dependencies.insert(mpValueReference);

        if (pVolumeReference)
          Dependencies.insert(pVolumeReference);

        mpConcReference->setDirectDependencies(Dependencies);
        mpConcReference->setRefresh(this, &CMetab::refreshConcentration);

        // Fixed values
        mRate = 0.0;
        mConcRate = 0.0;
        mTT = std::numeric_limits<C_FLOAT64>::infinity();
        break;

      case ASSIGNMENT:
        // Concentration
        success = mpExpression->compile(listOfContainer);
        mpConcReference->setDirectDependencies(mpExpression->getDirectDependencies());
        mpConcReference->setRefresh(this, &CMetab::calculate);

        // Implicit initial expression
        pdelete(mpInitialExpression);
        pDataModel = getObjectDataModel();
        mpInitialExpression = CExpression::createInitialExpression(*mpExpression, pDataModel);
        mpInitialExpression->setObjectName("InitialExpression");

        // Fixed values
        mRate = std::numeric_limits< C_FLOAT64 >::quiet_NaN();
        mConcRate = std::numeric_limits< C_FLOAT64 >::quiet_NaN();
        mTT = std::numeric_limits< C_FLOAT64 >::quiet_NaN();
        break;

      case ODE:
        // Concentration
        Dependencies.insert(mpValueReference);

        if (pVolumeReference)
          Dependencies.insert(pVolumeReference);

        mpConcReference->setDirectDependencies(Dependencies);
        mpConcReference->setRefresh(this, &CMetab::refreshConcentration);

        // Rate (particle number rate)
        success = mpExpression->compile(listOfContainer);
        Dependencies = mpExpression->getDirectDependencies();

        if (pVolumeReference)
          Dependencies.insert(pVolumeReference);

        mpRateReference->setDirectDependencies(Dependencies);
        mpRateReference->setRefresh(this, &CMetab::calculate);
        Dependencies.clear();

        // Concentration Rate
        Dependencies.insert(mpRateReference);
        Dependencies.insert(mpConcReference);

        if (pVolumeReference)
          Dependencies.insert(pVolumeReference);

        if (mpCompartment)
          Dependencies.insert(mpCompartment->getRateReference());

        mpConcRateReference->setDirectDependencies(Dependencies);
        mpConcRateReference->setRefresh(this, &CMetab::refreshConcentrationRate);
        Dependencies.clear();

        // Transition Time
        Dependencies.insert(mpValueReference);
        Dependencies.insert(mpRateReference);
        mpTTReference->setDirectDependencies(Dependencies);
        mpTTReference->setRefresh(this, &CMetab::refreshTransitionTime);
        Dependencies.clear();
        break;

      case REACTIONS:
        // Concentration
        Dependencies.insert(mpValueReference);

        if (pVolumeReference)
          Dependencies.insert(pVolumeReference);

        mpConcReference->setDirectDependencies(Dependencies);
        mpConcReference->setRefresh(this, &CMetab::refreshConcentration);

        Dependencies.clear();

        // Create the rate vector
        {
          CCopasiVectorN< CReaction >::const_iterator it = mpModel->getReactions().begin();
          CCopasiVectorN< CReaction >::const_iterator end = mpModel->getReactions().end();

          for (; it != end; ++it)
            {
              const CCopasiVector< CChemEqElement > &Balances =
                (*it)->getChemEq().getBalances();
              CCopasiVector< CChemEqElement >::const_iterator itChem = Balances.begin();
              CCopasiVector< CChemEqElement >::const_iterator endChem = Balances.end();

              for (; itChem != endChem; ++itChem)
                if ((*itChem)->getMetaboliteKey() == mKey)
                  break;

              if (itChem != endChem)
                {
                  Dependencies.insert((*it)->getParticleFluxReference());

                  std::pair< C_FLOAT64, const C_FLOAT64 * > Insert;
                  Insert.first = (*itChem)->getMultiplicity();
                  Insert.second = &(*it)->getParticleFlux();

                  mRateVector.push_back(Insert);
                }
            }
        }

        // Rate (particle number rate)
        mpRateReference->setRefresh(this, &CMetab::refreshRate);
        mpRateReference->setDirectDependencies(Dependencies);

        // Transition Time
        mpTTReference->setRefresh(this, &CMetab::refreshTransitionTime);
        mpTTReference->setDirectDependencies(Dependencies);
        Dependencies.clear();

        // Concentration Rate
        Dependencies.insert(mpRateReference);
        Dependencies.insert(mpConcReference);

        if (pVolumeReference)
          Dependencies.insert(pVolumeReference);

        if (mpCompartment)
          Dependencies.insert(mpCompartment->getRateReference());

        mpConcRateReference->setDirectDependencies(Dependencies);
        mpConcRateReference->setRefresh(this, &CMetab::refreshConcentrationRate);
        Dependencies.clear();

        break;
      default:
        break;
    }

  // The initial values
  success &= compileInitialValueDependencies(true);

  return success;
}

bool CMetab::compileInitialValueDependencies(const bool & updateConcentration)
{
  bool success = true;
  std::set<const CCopasiObject *> Dependencies;
  std::vector< CCopasiContainer * > listOfContainer;
  listOfContainer.push_back(getObjectAncestor("Model"));

  // If we have an initial expression we must update both
  if (mpInitialExpression != NULL &&
      mpInitialExpression->getInfix() != "")
    {
      // Initial concentration
      success &= mpInitialExpression->compile(listOfContainer);
      mpIConcReference->setDirectDependencies(mpInitialExpression->getDirectDependencies());

      // Initial particle number
      Dependencies.insert(mpIConcReference);

      if (mpCompartment)
        Dependencies.insert(mpCompartment->getInitialValueReference());

      mpIValueReference->setDirectDependencies(Dependencies);
      Dependencies.clear();

      // If we have constant initial expression, we update the initial concentration.
      if (mpInitialExpression->isUsable() &&
          mpIConcReference->getDirectDependencies().size() == 0)
        mIConc = mpInitialExpression->calcValue();

      return success;
    }

  if (updateConcentration)
    {
      // Initial particle number
      mpIValueReference->setDirectDependencies(Dependencies);

      // Initial concentration
      Dependencies.insert(mpIValueReference);

      if (mpCompartment)
        Dependencies.insert(mpCompartment->getInitialValueReference());

      mpIConcReference->setDirectDependencies(Dependencies);
    }
  else
    {
      // Initial concentration
      mpIConcReference->setDirectDependencies(Dependencies);

      // Initial particle number
      Dependencies.insert(mpIConcReference);

      if (mpCompartment)
        Dependencies.insert(mpCompartment->getInitialValueReference());

      mpIValueReference->setDirectDependencies(Dependencies);
    }

  return success;
}

bool CMetab::isInitialConcentrationChangeAllowed()
{
  compileInitialValueDependencies(false);

  std::set< const CCopasiObject * > Candidates;
  std::set< const CCopasiObject * > Verified;
  std::set< const CCopasiObject * > Context;

  bool Allowed = !mpIValueReference->hasCircularDependencies(Candidates, Verified, Context);

  if (!Allowed)
    compileInitialValueDependencies(true);

  return Allowed;
}

void CMetab::calculate()
{
  switch (getStatus())
    {
      case FIXED:
        break;

      case ASSIGNMENT:
        mConc = mpExpression->calcValue();
        break;

      case ODE:
        mRate = mpCompartment->getValue() * mpExpression->calcValue() * mpModel->getQuantity2NumberFactor();
        break;

      case REACTIONS:

        if (isDependent())
          *mpValue = mpMoiety->getDependentNumber();

        break;

      default:
        break;
    }
}

void CMetab::refreshRate()
{
  switch (getStatus())
    {
      case FIXED:
        break;

      case ASSIGNMENT:
        break;

      case ODE:
        break;

      case REACTIONS:
      {
        mRate = 0.0;

        std::vector< std::pair< C_FLOAT64, const C_FLOAT64 * > >::const_iterator it =
          mRateVector.begin();
        std::vector< std::pair< C_FLOAT64, const C_FLOAT64 * > >::const_iterator end =
          mRateVector.end();

        for (; it != end; ++it)
          mRate += it->first * *it->second;
      }
      break;

      default:
        break;
    }
}

void CMetab::refreshTransitionTime()
{
  switch (getStatus())
    {
      case FIXED:
        break;

      case ASSIGNMENT:
        break;

      case ODE:
        mTT = *mpValue / fabs(mRate);
        break;

      case REACTIONS:
      {
        C_FLOAT64 PositiveFlux = 0;
        C_FLOAT64 NegativeFlux = 0;
        C_FLOAT64 Flux;

        std::vector< std::pair< C_FLOAT64, const C_FLOAT64 * > >::const_iterator it =
          mRateVector.begin();
        std::vector< std::pair< C_FLOAT64, const C_FLOAT64 * > >::const_iterator end =
          mRateVector.end();

        for (; it != end; ++it)
          {
            Flux = it->first * *it->second;

            if (Flux > 0.0)
              PositiveFlux += Flux;
            else
              NegativeFlux -= Flux;
          }

        Flux = std::min(PositiveFlux, NegativeFlux);

        if (Flux == 0.0)
          mTT = std::numeric_limits<C_FLOAT64>::infinity();
        else
          mTT = *mpValue / Flux;
      }
      break;

      default:
        break;
    }
}

void CMetab::initObjects()
{
  mpIValueReference->setObjectName("InitialParticleNumber");
  mpIValueReference->setRefresh(this, &CMetab::refreshInitialValue);

  // We need to have mpValueRefernce point to a CParticleReference object.
  C_FLOAT64 * pValue = static_cast< C_FLOAT64 * >(mpValueReference->getValuePointer());
  assert(pValue != NULL);
  pdelete(mpValueReference);
  mpValueReference = new CParticleReference(this, *pValue);

  mpRateReference->setObjectName("ParticleNumberRate");

  mpIConcReference =
    static_cast<CCopasiObjectReference<C_FLOAT64> *>(addObjectReference("InitialConcentration", mIConc, CCopasiObject::ValueDbl));
  mpIConcReference->setRefresh(this, &CMetab::refreshInitialConcentration);

  mpConcReference = new CConcentrationReference(this, mConc);

  mpConcRateReference =
    static_cast<CCopasiObjectReference<C_FLOAT64> *>(addObjectReference("Rate", mConcRate, CCopasiObject::ValueDbl));

  mpTTReference =
    static_cast<CCopasiObjectReference<C_FLOAT64> *>(addObjectReference("TransitionTime", mTT, CCopasiObject::ValueDbl));
}

std::set< const CCopasiObject * > CMetab::getDeletedObjects() const
{
  std::set< const CCopasiObject * > Deleted = CModelEntity::getDeletedObjects();

  Deleted.insert(mpIConcReference);
  Deleted.insert(mpConcReference);
  Deleted.insert(mpConcRateReference);
  Deleted.insert(mpTTReference);

  return Deleted;
}

CCopasiObject * CMetab::getInitialConcentrationReference() const
{return mpIConcReference;}

CConcentrationReference * CMetab::getConcentrationReference() const
{return mpConcReference;}

C_FLOAT64 CMetab::getConcentrationRate() const
{
  const_cast<CMetab *>(this)->refreshConcentrationRate();
  return mConcRate;
}

void CMetab::refreshConcentrationRate()
{
  // d(Particle Number*conversion/Volume)/dt
  // = (Rate / Volume - Particle Number / Volume^2 * Volume Rate) * conversion
  // = (Rate - Particle Number / Volume * Volume Rate) / Volume * conversion
  // = (Rate  * conversion - Concentration * Volume Rate) / Volume

  mConcRate =
    (mRate * mpModel->getNumber2QuantityFactor() - mConc * mpCompartment->getRate())
    / mpCompartment->getValue();
}

// virtual
const CCopasiObject * CMetab::getValueObject() const
{
  return mpConcReference;
}

// virtual
void * CMetab::getValuePointer() const
{
  return const_cast<C_FLOAT64 *>(&mConc);
}

std::ostream & operator<<(std::ostream &os, const CMetab & d)
{
  os << "    ++++CMetab: " << d.getObjectName() << std::endl;
  os << "        mConc " << d.mConc << " mIConc " << d.mIConc << std::endl;
  os << "        mValue (particle number) " << *d.mpValue << " mIValue " << *d.mpIValue << std::endl;
  os << "        mRate " << d.mRate << " mTT " << d.mTT << " mStatus " << d.getStatus() << std::endl;

  if (d.mpCompartment)
    os << "        mpCompartment == " << d.mpCompartment << std::endl;
  else
    os << "        mpCompartment == 0 " << std::endl;

  if (d.mpModel)
    os << "        mpModel == " << d.mpModel << std::endl;
  else
    os << "        mpModel == 0 " << std::endl;

  os << "    ----CMetab " << std::endl;

  return os;
}

C_INT32 CMetab::load(CReadConfig &configbuffer)
{
  C_INT32 Fail = 0;

  std::string tmp;
  Fail = configbuffer.getVariable("Metabolite", "string",
                                  (void *) & tmp,
                                  CReadConfig::SEARCH);

  if (Fail)
    return Fail;

  setObjectName(tmp);

  Fail = configbuffer.getVariable("InitialConcentration", "C_FLOAT64",
                                  (void *) & mIConc);

  setInitialConcentration(mIConc);
  setConcentration(mIConc);

  Status GepasiStatus;
  Fail = configbuffer.getVariable("Type", "C_INT16",
                                  (void *) & GepasiStatus);

  if (Fail)
    return Fail;

  setStatus(GepasiStatus);

  // sanity check
  if ((GepasiStatus < 0) || (GepasiStatus > 7))
    {
      CCopasiMessage(CCopasiMessage::WARNING,
                     "The file specifies a non-existing type "
                     "for '%s'.\nReset to internal species.",
                     getObjectName().c_str());
      setStatus(REACTIONS);
    }

  // sanity check
  if ((GepasiStatus != METAB_MOIETY) && (mIConc < 0.0))
    {
      CCopasiMessage(CCopasiMessage::WARNING,
                     "The file specifies a negative concentration "
                     "for '%s'.\nReset to default.",
                     getObjectName().c_str());
      mIConc = 1.0;
    }

  return Fail;
}

std::string CMetab::getObjectDisplayName(bool regular, bool richtext) const
{
  CModel* tmp = dynamic_cast<CModel*>(this->getObjectAncestor("Model"));

  if (tmp)
    {
      return CMetabNameInterface::getDisplayName(tmp, *this);
    }

  return CCopasiObject::getObjectDisplayName(regular, richtext);
}

void CMetab::setDependentOn(const CMoiety * pMoiety)
{mpMoiety = pMoiety;}

bool CMetab::isDependent() const
{return mpMoiety != NULL;}

//******************* CMetabOld ***************************************************

CMetabOld::CMetabOld(const std::string & name,
                     const CCopasiContainer * pParent):
    CCopasiContainer(name, pParent, "Old Metabolite"),
    mIConc(1.0),
    mStatus(CModelEntity::REACTIONS),
    mCompartment()
{CONSTRUCTOR_TRACE;}

CMetabOld::CMetabOld(const CMetabOld & src,
                     const CCopasiContainer * pParent):
    CCopasiContainer(src, pParent),
    mIConc(src.mIConc),
    mStatus(src.mStatus),
    mCompartment(src.mCompartment)
{CONSTRUCTOR_TRACE;}

CMetabOld::~CMetabOld() {DESTRUCTOR_TRACE;}

void CMetabOld::cleanup() {}

C_INT32 CMetabOld::load(CReadConfig &configbuffer)
{
  C_INT32 Fail = 0;
  std::string tmp;
  Fail = configbuffer.getVariable("Metabolite", "string",
                                  (void *) & tmp,
                                  CReadConfig::SEARCH);

  if (Fail)
    return Fail;

  setObjectName(tmp);

  Fail = configbuffer.getVariable("Concentration", "C_FLOAT64",
                                  (void *) & mIConc);

  if (Fail)
    return Fail;

  Fail = configbuffer.getVariable("Compartment", "C_INT32",
                                  (void *) & mCompartment);

  if (Fail)
    return Fail;

  C_INT32 Status;

  Fail = configbuffer.getVariable("Type", "C_INT32",
                                  (void *) & Status);

  if (Status == 0)
    mStatus = CModelEntity::FIXED;
  else
    mStatus = CModelEntity::REACTIONS;

  // sanity check
  if ((mStatus < 0) || (mStatus > 7))
    {
      CCopasiMessage(CCopasiMessage::WARNING,
                     "The file specifies a non-existing type "
                     "for '%s'.\nReset to internal species.",
                     getObjectName().c_str());
      mStatus = CModelEntity::REACTIONS;
    }

  // sanity check
  if ((mStatus != METAB_MOIETY) && (mIConc < 0.0))
    {
      CCopasiMessage(CCopasiMessage::WARNING,
                     "The file specifies a negative concentration "
                     "for '%s'.\nReset to default.",
                     getObjectName().c_str());
      mIConc = 1.0;
    }

  return Fail;
}

C_INT32 CMetabOld::getIndex() const {return mCompartment;}

// static
std::set< const CCopasiObject * > CConcentrationReference::EmptyDependencies;

CConcentrationReference::CConcentrationReference(const CCopasiContainer * pParent,
    C_FLOAT64 & reference) :
    CCopasiObjectReference< C_FLOAT64 >("Concentration", pParent, reference),
    mpApplyInitialValuesRefresh(NULL)
{
  const CMetab * pMetab = static_cast< const CMetab * >(pParent);

  mpApplyInitialValuesRefresh =
    new RefreshTemplate< CMetab >(const_cast< CMetab * >(pMetab),
                                  &CMetab::refreshConcentration);
}

CConcentrationReference::CConcentrationReference(const CConcentrationReference & src,
    const CCopasiContainer * pParent) :
    CCopasiObjectReference< C_FLOAT64 >(src, pParent),
    mpApplyInitialValuesRefresh(NULL)
{
  const CMetab * pMetab = static_cast< const CMetab * >(pParent);

  mpApplyInitialValuesRefresh =
    new RefreshTemplate< CMetab >(const_cast< CMetab * >(pMetab),
                                  &CMetab::refreshConcentration);
}

CConcentrationReference::~CConcentrationReference()
{
  pdelete(mpApplyInitialValuesRefresh);
}

// virtual
const std::set< const CCopasiObject * > &
CConcentrationReference::getDirectDependencies(const std::set< const CCopasiObject * > & context) const
{
  // If the concentration was changed in the context it has no further dependencies
  if (context.count(this) > 0)
    {
      return EmptyDependencies;
    }

  return CCopasiObjectReference< C_FLOAT64 >::getDirectDependencies();
}

Refresh * CConcentrationReference::getApplyInitialValueRefresh() const
{
  return mpApplyInitialValuesRefresh;
}

// static
std::set< const CCopasiObject * > CParticleReference::EmptyDependencies;

CParticleReference::CParticleReference(const CCopasiContainer * pParent,
                                       C_FLOAT64 & reference) :
    CCopasiObjectReference< C_FLOAT64 >("ParticleNumber", pParent, reference)
{}

CParticleReference::CParticleReference(const CParticleReference & src,
                                       const CCopasiContainer * pParent) :
    CCopasiObjectReference< C_FLOAT64 >(src, pParent)
{}

CParticleReference::~CParticleReference()
{}

// virtual
const std::set< const CCopasiObject * > &
CParticleReference::getDirectDependencies(const std::set< const CCopasiObject * > & context) const
{
  // If the particle number was changed in the context it has no further dependencies
  if (context.count(this) > 0)
    {
      return EmptyDependencies;
    }

  const CMetab * pSpecies = static_cast< const CMetab * >(getObjectParent());

  // In an assignment the particles have always dependencies
  if (pSpecies == NULL ||
      pSpecies->getStatus() == CModelEntity::ASSIGNMENT)
    {
      return CCopasiObjectReference< C_FLOAT64 >::getDirectDependencies();
    }

  // In other cases we need to find out whether the concentration was changed in
  // the context.
  const CCopasiObject * pConcentrationReference = pSpecies->getConcentrationReference();

  if (pConcentrationReference != NULL &&
      context.count(pConcentrationReference) > 0)
    {
      return CCopasiObjectReference< C_FLOAT64 >::getDirectDependencies();
    }

  return EmptyDependencies;
}
