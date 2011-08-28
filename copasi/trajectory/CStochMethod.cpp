// Begin CVS Header
//   $Source: /fs/turing/cvs/copasi_dev/copasi/trajectory/CStochMethod.cpp,v $
//   $Revision: 1.78.2.1 $
//   $Name: Build-33 $
//   $Author: shoops $
//   $Date: 2010/10/01 13:26:32 $
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

#ifdef WIN32
# pragma warning (disable: 4786)
# pragma warning (disable: 4243)
// warning C4355: 'this' : used in base member initializer list
# pragma warning (disable: 4355)
#endif  // WIN32

#include <limits.h>

#include <mathematics.h>
#include <vector>
#include <numeric>
#include <limits>
#include <set>

#include "copasi.h"
#include "CStochMethod.h"
#include "CStochDirectMethod.h"
#include "CStochNextReactionMethod.h"
#include "utilities/CCopasiVector.h"
#include "function/CFunction.h"
#include "randomGenerator/CRandom.h"
#include "CTrajectoryMethod.h"
#include "CTrajectoryProblem.h"
#include "model/CState.h"
#include "model/CCompartment.h"
#include "model/CModel.h"

C_INT32 CStochMethod::checkModel(CModel * C_UNUSED(pmodel))
{
  // Here several checks will be performed to validate the model
  return 2; // suggest next reaction method
}

CStochMethod::CStochMethod(const CCopasiContainer * pParent):
    CTrajectoryMethod(CCopasiMethod::stochastic, pParent)
{
  initializeParameter();
  mpRandomGenerator = CRandom::createGenerator(CRandom::mt19937);
}

CStochMethod::CStochMethod(const CStochMethod & src,
                           const CCopasiContainer * pParent):
    CTrajectoryMethod(src, pParent)
{
  initializeParameter();
  mpRandomGenerator = CRandom::createGenerator(CRandom::mt19937);
}

CStochMethod::~CStochMethod()
{
  delete mpRandomGenerator;
  mpRandomGenerator = NULL;
}

void CStochMethod::initializeParameter()
{
  CCopasiParameter *pParm;

  assertParameter("Max Internal Steps", CCopasiParameter::INT, (C_INT32) 1000000);
  assertParameter("Subtype", CCopasiParameter::UINT, (unsigned C_INT32) 2);
  assertParameter("Use Random Seed", CCopasiParameter::BOOL, false);
  assertParameter("Random Seed", CCopasiParameter::UINT, (unsigned C_INT32) 1);

  // Check whether we have a method with the old parameter names
  if ((pParm = getParameter("STOCH.MaxSteps")) != NULL)
    {
      setValue("Max Internal Steps", *pParm->getValue().pINT);
      removeParameter("STOCH.MaxSteps");

      if ((pParm = getParameter("STOCH.Subtype")) != NULL)
        {
          setValue("Subtype", *pParm->getValue().pUINT);
          removeParameter("STOCH.Subtype");
        }

      if ((pParm = getParameter("STOCH.UseRandomSeed")) != NULL)
        {
          setValue("Use Random Seed", *pParm->getValue().pBOOL);
          removeParameter("STOCH.UseRandomSeed");
        }

      if ((pParm = getParameter("STOCH.RandomSeed")) != NULL)
        {
          setValue("Random Seed", *pParm->getValue().pUINT);
          removeParameter("STOCH.RandomSeed");
        }
    }
}

bool CStochMethod::elevateChildren()
{
  initializeParameter();
  return true;
}

CTrajectoryMethod::Status CStochMethod::step(const double & deltaT)
{
  // write the current state to the model:
  //mpProblem->getModel()->setState(mpCurrentState); //?

  // check for possible overflows:
  unsigned C_INT32 i;
  unsigned C_INT32 imax;

  // :TODO: Bug 774: This assumes that the number of variable metabs is the number
  // of metabs determined by reaction. In addition they are expected at the beginning of the
  // MetabolitesX which is not the case if we have metabolites of type ODE.
  for (i = 0, imax = mpProblem->getModel()->getNumVariableMetabs(); i < imax; i++)
    if (mpProblem->getModel()->getMetabolitesX()[i]->getValue() >= mMaxIntBeforeStep)
      {
        CCopasiMessage(CCopasiMessage::EXCEPTION, "at least one particle number got to big.");
        // TODO:throw exception or something like that
      }

  // do several steps:
  C_FLOAT64 time = mpCurrentState->getTime();
  C_FLOAT64 endtime = time + deltaT;

  size_t Steps = 0;

  while (time < endtime)
    {
      time = doSingleStep(time, endtime);

      if (++Steps > mMaxSteps)
        {
          CCopasiMessage(CCopasiMessage::EXCEPTION, MCTrajectoryMethod + 12);
        }
    }

  mpCurrentState->setTime(time);

  // get back the particle numbers:

  /* Set the variable Metabolites */
  C_FLOAT64 * Dbl = mpCurrentState->beginIndependent() + mFirstMetabIndex - 1;

  for (i = 0, imax = mpProblem->getModel()->getNumVariableMetabs(); i < imax; i++, Dbl++)
    *Dbl = mpProblem->getModel()->getMetabolitesX()[i]->getValue();

  return NORMAL;
}

void CStochMethod::start(const CState * initialState)
{
  /* get configuration data */
  mMaxSteps = * getValue("Max Internal Steps").pINT;

  bool useRandomSeed = * getValue("Use Random Seed").pBOOL;
  unsigned C_INT32 randomSeed = * getValue("Random Seed").pUINT;

  if (useRandomSeed) mpRandomGenerator->initialize(randomSeed);

  //mpCurrentState is initialized. This state is not used internally in the
  //stochastic solver, but it is used for returning the result after each step.
  *mpCurrentState = *initialState;

  mpModel = mpProblem->getModel();
  assert(mpModel);

  if (mpModel->getModelType() == CModel::deterministic)
    mDoCorrection = true;
  else
    mDoCorrection = false;

  mHasAssignments = modelHasAssignments(mpModel);

  unsigned C_INT32 i;

  //initialize the vector of ints that contains the particle numbers
  //for the discrete simulation. This also floors all particle numbers in the model.
  mNumbers.resize(mpModel->getMetabolitesX().size());

  for (i = 0; i < mNumbers.size(); ++i)
    {
      mNumbers[i] = (C_INT64) mpModel->getMetabolitesX()[i]->getValue();
      mpModel->getMetabolitesX()[i]->setValue(mNumbers[i]);
      mpModel->getMetabolitesX()[i]->refreshConcentration();
    }

  mFirstMetabIndex = mpModel->getStateTemplate().getIndex(mpModel->getMetabolitesX()[0]);

  mpModel->updateSimulatedValues(false); //for assignments
  //mpModel->updateNonSimulatedValues(); //for assignments

  mNumReactions = mpModel->getReactions().size();

  mAmu.clear(); mAmuOld.clear();

  for (i = 0; i < mNumReactions; i++)
    {
      mAmu.push_back(0);
      mAmuOld.push_back(0);
    }

  setupDependencyGraphAndBalances();
  updatePropensities();

  // call init of the specific simulation method
  initMethod(mpCurrentState->getTime());

  mMaxStepsReached = false;
  return;
}

C_INT32 CStochMethod::updatePropensities()
{
  //mA0Old = mA0;
  mA0 = 0;

  for (unsigned C_INT32 i = 0; i < mNumReactions; i++)
    {
      mAmuOld[i] = mAmu[i];
      calculateAmu(i);
      mA0 += mAmu[i];
    }

  return 0;
}

C_INT32 CStochMethod::calculateAmu(C_INT32 index)
{
  if (!mDoCorrection)
    {
      mAmu[index] = mpModel->getReactions()[index]->calculateParticleFlux();
      return 0;
    }

  // We need the product of the cmu and hmu for this step.
  // We calculate this in one go, as there are fewer steps to
  // perform and we eliminate some possible rounding errors.
  C_FLOAT64 amu = 1; // initially
  //C_INT32 total_substrates = 0;
  C_INT32 num_ident = 0;
  C_INT64 number = 0;
  C_INT64 lower_bound;
  // substrate_factor - The substrates, raised to their multiplicities,
  // multiplied with one another. If there are, e.g. m substrates of type m,
  // and n of type N, then substrate_factor = M^m * N^n.
  C_FLOAT64 substrate_factor = 1;
  // First, find the reaction associated with this index.
  // Keep a pointer to this.
  // Iterate through each substrate in the reaction
  const std::vector<CStochBalance> & substrates = mLocalSubstrates[index];

  int flag = 0;

  for (unsigned C_INT32 i = 0; i < substrates.size(); i++)
    {
      num_ident = substrates[i].mMultiplicity;

      if (num_ident > 1)
        {
          flag = 1;
          number = mNumbers[substrates[i].mIndex];
          lower_bound = number - num_ident;
          substrate_factor = substrate_factor * pow((double) number, (int)(num_ident - 1));  //optimization

          number--; //optimization

          while (number > lower_bound)
            {
              amu *= number;
              number--;
            }
        }
    }

  if ((amu == 0) || (substrate_factor == 0))  // at least one substrate particle number is zero
    {
      mAmu[index] = 0;
      return 0;
    }

  // rate_factor is the rate function divided by substrate_factor.
  // It would be more efficient if this was generated directly, since in effect we
  // are multiplying and then dividing by the same thing (substrate_factor)!
  C_FLOAT64 rate_factor = mpModel->getReactions()[index]->calculateParticleFlux();

  if (flag)
    {
      //cout << "Rate factor = " << rate_factor << endl;
      amu *= rate_factor / substrate_factor;
      mAmu[index] = amu;
    }
  else
    {
      mAmu[index] = rate_factor;
    }

  return 0;
}

C_INT32 CStochMethod::updateSystemState(C_INT32 rxn)
{
  // Change the particle numbers according to which step took place.
  // First, get the vector of balances in the reaction we've got.
  // (This vector expresses the number change of each metabolite
  // in the reaction.) Then step through each balance, using its
  // multiplicity to calculate a new value for the associated
  // metabolite. Finally, update the metabolite.

  const std::vector<CStochBalance> & bals = mLocalBalances[rxn];

  std::vector<CStochBalance>::const_iterator bi;
  CMetab* pTmpMetab;

  for (bi = bals.begin(); bi != bals.end(); bi++)
    {
      mNumbers[bi->mIndex] = mNumbers[bi->mIndex] + bi->mMultiplicity;
      pTmpMetab = mpModel->getMetabolitesX()[bi->mIndex];
      pTmpMetab->setValue(mNumbers[bi->mIndex]);
      pTmpMetab->refreshConcentration();
    }

  if (mHasAssignments)
    {
      // this is less efficient but can deal with assignments.
      //TODO: handle dependencies for assignments also.
      mpModel->updateSimulatedValues(false);

      //now potentially species with assignments have non integer
      //particle numbers. This needs to be rounded. Also the updated
      //particle numbers need to be copied to the vector of integers.
      //(the integer values may be used to calculate the propensities for
      //higher order kinetics).
      unsigned C_INT32 i, imax = mNumbers.size();

      for (i = 0; i < imax; ++i)
        {
          if (mpModel->getMetabolitesX()[i]->getStatus() == CModelEntity::ASSIGNMENT)
            {
              mNumbers[i] = (C_INT64) mpModel->getMetabolitesX()[i]->getValue();
              mpModel->getMetabolitesX()[i]->setValue(mNumbers[i]);
              mpModel->getMetabolitesX()[i]->refreshConcentration();
            }
        }

      //now the propensities can be updated
      updatePropensities();
    }
  else
    {
      const std::set<unsigned C_INT32> & dep_nodes = mDG.getDependents(rxn);

      std::set<unsigned C_INT32>::const_iterator it;
      unsigned int ii;

      for (it = dep_nodes.begin(); it != dep_nodes.end(); it++)
        {
          ii = *it;
          mAmuOld[ii] = mAmu[ii];
          calculateAmu(ii);
        }

      //mA0Old = mA0;

      mA0 = 0;
      mA0 = std::accumulate(mAmu.begin(), mAmu.end(), mA0);
    }

  return 0;
}

C_INT32 CStochMethod::generateReactionIndex()
{
  C_FLOAT64 rand1 = mpRandomGenerator->getRandomOO() * mA0;
  C_FLOAT64 sum = 0;
  unsigned C_INT32 index = 0;

  while (index < (mpModel->getReactions().size() - 1))
    {
      sum += mAmu[index] /* /mA0 */;

      if (rand1 <= sum)
        {return index;}

      index++;
    }

  return index;
}

C_FLOAT64 CStochMethod::generateReactionTime()
{
  if (mA0 == 0) return std::numeric_limits<C_FLOAT64>::infinity();

  C_FLOAT64 rand2 = mpRandomGenerator->getRandomOO();
  return - 1 * log(rand2) / mA0;
}

C_FLOAT64 CStochMethod::generateReactionTime(C_INT32 reaction_index)
{
  if (mAmu[reaction_index] == 0) return std::numeric_limits<C_FLOAT64>::infinity();

  C_FLOAT64 rand2 = mpRandomGenerator->getRandomOO();
  return - 1 * log(rand2) / mAmu[reaction_index];
}

void CStochMethod::setupDependencyGraphAndBalances()
{
  mDG.clear();
  std::vector< std::set<std::string>* > DependsOn;
  std::vector< std::set<std::string>* > Affects;
  unsigned C_INT32 i, j;
  // Do for each reaction:

  for (i = 0; i < mNumReactions; i++)
    {
      // Get the set of metabolites  which affect the value of amu for this
      // reaction i.e. the set on which amu depends. This may be  more than
      // the set of substrates, since the kinetics can involve other
      // reactants, e.g. catalysts. We thus need to step through the
      // rate function and pick out every reactant which can vary.
      DependsOn.push_back(getDependsOn(i));
      // Get the set of metabolites which are affected when this reaction takes place
      Affects.push_back(getAffects(i));
    }

  // For each possible pair of reactions i and j, if the intersection of
  // Affects(i) with DependsOn(j) is non-empty, add a dependency edge from i to j.
  mDG.resize(mNumReactions);

  for (i = 0; i < mNumReactions; i++)
    {
      for (j = 0; j < mNumReactions; j++)
        {
          // Determine whether the intersection of these two sets is non-empty
          // Could also do this with set_intersection generic algorithm, but that
          // would require operator<() to be defined on the set elements.

          std::set<std::string>::iterator iter = Affects[i]->begin();

          for (; iter != Affects[i]->end(); iter++)
            {
              if (DependsOn[j]->count(*iter))
                {
                  // The set intersection is non-empty
                  mDG.addDependent(i, j);
                  break;
                }
            }
        }

      // Ensure that self edges are included
      //mDG.addDependent(i, i);
    }

  // Create local copy of balances and substrates list
  CStochBalance bb;
  C_INT32 maxBalance = 0;

  mLocalBalances.clear();
  mLocalBalances.resize(mNumReactions);
  mLocalSubstrates.clear();
  mLocalSubstrates.resize(mNumReactions);

  for (i = 0; i < mNumReactions; i++)
    {
      const CCopasiVector<CChemEqElement> * bbb;

      bbb = &mpModel->getReactions()[i]->getChemEq().getBalances();

      //TODO clear old local balances and substrates
      for (j = 0; j < bbb->size(); j++)
        {
          //bb.mIndex = mpModel->getMetabolites().getIndex((*bbb)[j]->getMetabolite().getObjectName(());
          //bb.mIndex = mpModel->findMetabByKey((*bbb)[j]->getMetaboliteKey());
          assert((*bbb)[j]->getMetabolite());
          bb.mIndex = mpModel->getMetabolitesX().getIndex((*bbb)[j]->getMetabolite());
          bb.mMultiplicity = static_cast<C_INT32>(floor((*bbb)[j]->getMultiplicity() + 0.5));

          if (((*bbb)[j]->getMetabolite()->getStatus()) != CModelEntity::FIXED)
            {
              if (bb.mMultiplicity > maxBalance) maxBalance = bb.mMultiplicity;

              mLocalBalances[i].push_back(bb);
            }
        }

      bbb = &mpModel->getReactions()[i]->getChemEq().getSubstrates();

      for (j = 0; j < bbb->size(); j++)
        {
          //bb.mIndex = mpModel->getMetabolites().getIndex((*bbb)[j]->getMetabolite().getObjectName(());
          //bb.mIndex = mpModel->findMetabByKey((*bbb)[j]->getMetaboliteKey());
          assert((*bbb)[j]->getMetabolite());
          bb.mIndex = mpModel->getMetabolitesX().getIndex((*bbb)[j]->getMetabolite());
          bb.mMultiplicity = static_cast<C_INT32>(floor((*bbb)[j]->getMultiplicity() + 0.5));

          if (1)
            {
              mLocalSubstrates[i].push_back(bb);
            }
        }
    }

  mMaxBalance = maxBalance;
  mMaxIntBeforeStep =        /*INT_MAX*/ LLONG_MAX - 1 - mMaxSteps * mMaxBalance;

  // Delete the memory allocated in getDependsOn() and getAffects()
  // since this is allocated in other functions.
  for (i = 0; i < mNumReactions; i++)
    {
      delete DependsOn[i];
      delete Affects[i];
    }
}

std::set<std::string> *CStochMethod::getDependsOn(C_INT32 reaction_index)
{
  std::set<std::string> *retset = new std::set<std::string>;

  unsigned C_INT32 i, imax = mpModel->getReactions()[reaction_index]->getFunctionParameters().size();
  unsigned C_INT32 j, jmax;

  for (i = 0; i < imax; ++i)
    {
      if (mpModel->getReactions()[reaction_index]->getFunctionParameters()[i]->getUsage() == CFunctionParameter::PARAMETER)
        continue;

      //metablist = mpModel->getReactions()[reaction_index]->getParameterMappingMetab(i);
      const std::vector <std::string> & metabKeylist =
        mpModel->getReactions()[reaction_index]->getParameterMappings()[i];
      jmax = metabKeylist.size();

      for (j = 0; j < jmax; ++j)
        {
          retset->insert(metabKeylist[j]);
        }
    }

  return retset;
}

std::set<std::string> *CStochMethod::getAffects(C_INT32 reaction_index)
{
  std::set<std::string> *retset = new std::set<std::string>;

  // Get the balances  associated with the reaction at this index
  // XXX We first get the chemical equation, then the balances, since the getBalances method in CReaction is unimplemented!
  const CCopasiVector<CChemEqElement> & balances = mpModel->getReactions()[reaction_index]->getChemEq().getBalances();

  for (unsigned C_INT32 i = 0; i < balances.size(); i++)
    {
      if (!balances[i]->getMetabolite()) continue;

      if (fabs(balances[i]->getMultiplicity()) >= 0.1)
        if (balances[i]->getMetabolite()->getStatus() != CModelEntity::FIXED)
          {
            retset->insert(balances[i]->getMetabolite()->getKey());
          }
    }

  return retset;
}

//virtual
bool CStochMethod::isValidProblem(const CCopasiProblem * pProblem)
{
  if (!CTrajectoryMethod::isValidProblem(pProblem)) return false;

  const CTrajectoryProblem * pTP = dynamic_cast<const CTrajectoryProblem *>(pProblem);

  if (pTP->getDuration() < 0.0)
    {
      //back integration not possible
      CCopasiMessage(CCopasiMessage::ERROR, MCTrajectoryMethod + 9);
      return false;
    }

  if (pTP->getModel()->getTotSteps() < 1)
    {
      //at least one reaction necessary
      CCopasiMessage(CCopasiMessage::ERROR, MCTrajectoryMethod + 17);
      return false;
    }

  //check for rules
  C_INT32 i, imax = pTP->getModel()->getNumModelValues();

  for (i = 0; i < imax; ++i)
    {
      if (pTP->getModel()->getModelValues()[i]->getStatus() == CModelEntity::ODE)
        {
          //ode rule found
          CCopasiMessage(CCopasiMessage::ERROR, MCTrajectoryMethod + 18);
          return false;
        }
    }

  imax = pTP->getModel()->getNumMetabs();

  for (i = 0; i < imax; ++i)
    {
      if (pTP->getModel()->getMetabolites()[i]->getStatus() == CModelEntity::ODE)
        {
          //ode rule found
          CCopasiMessage(CCopasiMessage::ERROR, MCTrajectoryMethod + 20);
          return false;
        }
    }

  imax = pTP->getModel()->getCompartments().size();

  for (i = 0; i < imax; ++i)
    {
      if (pTP->getModel()->getCompartments()[i]->getStatus() == CModelEntity::ODE)
        {
          //ode rule found
          CCopasiMessage(CCopasiMessage::ERROR, MCTrajectoryMethod + 21);
          return false;
        }
    }

  //   if (modelHasAssignments(pTP->getModel()))
  //     {
  //       CCopasiMessage(CCopasiMessage::EXCEPTION, MCTrajectoryMethod + 19);
  //       return false;
  //}
  //  this test is disabled since we have preliminary support for assignments

  //TODO: rewrite CModel::suitableForStochasticSimulation() to use
  //      CCopasiMessage
  std::string message = pTP->getModel()->suitableForStochasticSimulation();

  if (message != "")
    {
      //model not suitable, message describes the problem
      CCopasiMessage(CCopasiMessage::ERROR, message.c_str());
      return false;
    }

  if (* getValue("Max Internal Steps").pINT <= 0)
    {
      //max steps should be at least 1
      CCopasiMessage(CCopasiMessage::ERROR, MCTrajectoryMethod + 15);
      return false;
    }

  //events are not supported at the moment
  if (pTP->getModel()->getEvents().size() > 0)
    {
      CCopasiMessage(CCopasiMessage::ERROR, MCTrajectoryMethod + 23);
      return false;
    }

  return true;
}

//static
bool CStochMethod::modelHasAssignments(const CModel* pModel)
{
  C_INT32 i, imax = pModel->getNumModelValues();

  for (i = 0; i < imax; ++i)
    {
      if (pModel->getModelValues()[i]->getStatus() == CModelEntity::ASSIGNMENT)
        if (pModel->getModelValues()[i]->isUsed())
          {
            //used assignment found
            return true;
          }
    }

  imax = pModel->getNumMetabs();

  for (i = 0; i < imax; ++i)
    {
      if (pModel->getMetabolites()[i]->getStatus() == CModelEntity::ASSIGNMENT)
        if (pModel->getMetabolites()[i]->isUsed())
          {
            //used assignment found
            return true;
          }
    }

  imax = pModel->getCompartments().size();

  for (i = 0; i < imax; ++i)
    {
      if (pModel->getCompartments()[i]->getStatus() == CModelEntity::ASSIGNMENT)
        if (pModel->getCompartments()[i]->isUsed())
          {
            //used assignment found
            return true;
          }
    }

  return false;
}
