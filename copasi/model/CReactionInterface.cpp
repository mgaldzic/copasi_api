// Begin CVS Header
//   $Source: /fs/turing/cvs/copasi_dev/copasi/model/CReactionInterface.cpp,v $
//   $Revision: 1.39 $
//   $Name: Build-33 $
//   $Author: shoops $
//   $Date: 2009/10/27 16:52:47 $
// End CVS Header

// Copyright (C) 2008 by Pedro Mendes, Virginia Tech Intellectual
// Properties, Inc., EML Research, gGmbH, University of Heidelberg,
// and The University of Manchester.
// All rights reserved.

// Copyright (C) 2001 - 2007 by Pedro Mendes, Virginia Tech Intellectual
// Properties, Inc. and EML Research, gGmbH.
// All rights reserved.

#include <string>

#include "CReactionInterface.h"
#include "CReaction.h"
#include "CModel.h"
#include "CChemEqElement.h"
#include "CopasiDataModel/CCopasiDataModel.h"
#include "function/CFunctionDB.h"
#include "report/CKeyFactory.h"
#include "model/CMetabNameInterface.h"
#include "report/CCopasiRootContainer.h"

CReactionInterface::CReactionInterface(CModel * pModel):
    mpModel(pModel),
    mChemEqI(pModel),
    mpFunction(NULL),
    mpParameters(NULL)
{
  assert(mpModel != NULL);
  emptyString = "";
}

CReactionInterface::~CReactionInterface()
{
  pdelete(mpParameters);
  //pdelete(mpChemEq);
}

const std::vector<std::string> & CReactionInterface::getListOfMetabs(CFunctionParameter::Role role) const
{
  return mChemEqI.getListOfDisplayNames(role);
}

std::vector< std::string > CReactionInterface::getListOfPossibleFunctions() const
{
  TriLogic reversible;

  if (isReversible() == false)
    reversible = TriFalse;
  else
    reversible = TriTrue;

  std::vector<CFunction*> functionVector =
    CCopasiRootContainer::getFunctionList()->suitableFunctions(
      mChemEqI.getMolecularity(CFunctionParameter::SUBSTRATE),
      mChemEqI.getMolecularity(CFunctionParameter::PRODUCT),
      reversible);

  std::vector<std::string> ret;
  unsigned C_INT32 i, imax = functionVector.size();

  for (i = 0; i < imax; ++i)
    ret.push_back(functionVector[i]->getObjectName());

  return ret;
}

void CReactionInterface::initFromReaction(const std::string & key)
{
  mReactionReferenceKey = key;

  const CReaction *rea;
  rea = dynamic_cast< CReaction *>(CCopasiRootContainer::getKeyFactory()->get(key));
  assert(rea);
  initFromReaction(rea);
}

void CReactionInterface::initFromReaction(const C_INT32 index)
{
  const CReaction *rea =
    (*CCopasiRootContainer::getDatamodelList())[0]->getModel()->getReactions()[index];

  mReactionReferenceKey = rea->getKey();
  assert(rea);
  initFromReaction(rea);
}

void CReactionInterface::initFromReaction(const CReaction *rea)
{
  //name
  mReactionName = rea->getObjectName();

  //chemical equation
  mChemEqI.loadFromChemEq(rea->getChemEq());

  if (rea->getFunction() && (rea->getFunction() != CCopasiRootContainer::getUndefinedFunction()))
    {
      //function
      mpFunction = rea->getFunction();
      pdelete(mpParameters)
      mpParameters = new CFunctionParameters(mpFunction->getVariables());

      //mapping
      loadMappingAndValues(*rea);
    }
  else
    {
      setFunctionWithEmptyMapping("");
    }

#ifdef COPASI_DEBUG
  printDebug();
#endif // COPASI_DEBUG
}

bool CReactionInterface::loadMappingAndValues(const CReaction & rea)
{
  bool success = true;
  std::vector< std::vector<std::string> >::const_iterator it;
  std::vector< std::vector<std::string> >::const_iterator iEnd;
  std::vector<std::string>::const_iterator jt;
  std::vector<std::string>::const_iterator jEnd;
  unsigned C_INT32 i;

  std::string metabName;
  const CModelEntity* pObj;
  std::vector<std::string> SubList;

  mNameMap.clear();
  mValues.resize(size(), 0.1);
  mIsLocal.resize(size(), false);

  it = rea.getParameterMappings().begin();
  iEnd = rea.getParameterMappings().end();

  for (i = 0; it != iEnd; ++it, ++i)
    {
      if (isVector(i))
        {
          assert((getUsage(i) == CFunctionParameter::SUBSTRATE)
                 || (getUsage(i) == CFunctionParameter::PRODUCT)
                 || (getUsage(i) == CFunctionParameter::MODIFIER));

          SubList.clear();

          for (jt = it->begin(), jEnd = it->end(); jt != jEnd; ++jt)
            {
              metabName = CMetabNameInterface::getDisplayName(mpModel, *jt);
              assert(metabName != "");
              SubList.push_back(metabName);
            }
        }
      else
        {
          assert(it->size() == 1);
          SubList.resize(1); SubList[0] = "unknown";

          switch (getUsage(i))
            {
              case CFunctionParameter::SUBSTRATE:
              case CFunctionParameter::PRODUCT:
              case CFunctionParameter::MODIFIER:
                metabName = CMetabNameInterface::getDisplayName(mpModel, *(it->begin()));
                assert(metabName != "");
                SubList[0] = metabName;
                //TODO: check if the metabolite is in the chemical equation with the correct rule
                break;

              case CFunctionParameter::VOLUME:
                pObj = dynamic_cast<const CCompartment*>(CCopasiRootContainer::getKeyFactory()->get(*(it->begin())));
                assert(pObj);
                SubList[0] = pObj->getObjectName();
                break;

              case CFunctionParameter::TIME:
                pObj = dynamic_cast<const CModel*>(CCopasiRootContainer::getKeyFactory()->get(*(it->begin())));
                assert(pObj);
                SubList[0] = pObj->getObjectName();
                break;

              case CFunctionParameter::PARAMETER:
                mIsLocal[i] = rea.isLocalParameter(i);
                mValues[i] = rea.getParameterValue(getParameterName(i));
                pObj = dynamic_cast<const CModelValue*>(CCopasiRootContainer::getKeyFactory()->get(*(it->begin())));

                //assert(pObj);
                if (pObj) SubList[0] = pObj->getObjectName();

                break;

              default:
                break;
            }
        }

      mNameMap.push_back(SubList);
    }

  return success;
}

bool CReactionInterface::writeBackToReaction(CReaction * rea)
{
  if (!isValid()) return false; // do nothing

  if (!(*mpParameters == mpFunction->getVariables())) return false; // do nothing

  bool success = true;

  //CReaction *rea;
  if (rea == NULL)
    rea = dynamic_cast< CReaction *>(CCopasiRootContainer::getKeyFactory()->get(mReactionReferenceKey));

  if (rea == NULL) return false;

  if (!rea->setObjectName(mReactionName))
    success = false;

  // Now we can safely write to the equation as we are sure that only unique metabolites
  // may have the empty string as compartments
  mChemEqI.writeToChemEq(rea->getChemEq());

  // TODO. check if function has changed since it was set in the R.I.
  rea->setFunction(mpFunction->getObjectName());

  unsigned C_INT32 j, jmax;
  unsigned C_INT32 i, imax = size();
  std::pair< std::string, std::string > Names;

  for (i = 0; i < imax; ++i)
    {
      switch (getUsage(i))
        {
          case CFunctionParameter::PARAMETER:

            if (mIsLocal[i])
              rea->setParameterValue(getParameterName(i), mValues[i]);
            else
              {
                rea->setParameterValue(getParameterName(i), mValues[i], false);
                rea->setParameterMapping(i, mpModel->getModelValues()[mNameMap[i][0]]->getKey());
              }

            break;

          case CFunctionParameter::VOLUME:
            rea->setParameterMapping(i, mpModel->getCompartments()[mNameMap[i][0]]->getKey());
            break;

          case CFunctionParameter::TIME:
            rea->setParameterMapping(i, mpModel->getKey()); //time is the value of the model
            break;

          case CFunctionParameter::SUBSTRATE:
          case CFunctionParameter::PRODUCT:
          case CFunctionParameter::MODIFIER:

            if (isVector(i))
              {
                rea->clearParameterMapping(i);
                jmax = mNameMap[i].size();

                for (j = 0; j < jmax; ++j)
                  {
                    Names = CMetabNameInterface::splitDisplayName(mNameMap[i][j]);
                    rea->addParameterMapping(i, CMetabNameInterface::getMetaboliteKey(mpModel, Names.first, Names.second));
                  }
              }
            else
              {
                Names = CMetabNameInterface::splitDisplayName(mNameMap[i][0]);
                rea->setParameterMapping(i, CMetabNameInterface::getMetaboliteKey(mpModel, Names.first, Names.second));
              }

            break;

          default:
            break;
        }
    }

  rea->compile();
  mpModel->setCompileFlag(); //TODO: check if really necessary

  return success;
}

void CReactionInterface::clearFunction()
{
  mpFunction = NULL;
  pdelete(mpParameters);
  //mValid = false;

  mValues.clear();
  mNameMap.clear();
}

void CReactionInterface::setChemEqString(const std::string & eq, const std::string & newFunction)
{
  mChemEqI.setChemEqString(eq);
  findAndSetFunction(newFunction);
}

void CReactionInterface::setReversibility(bool rev, const std::string & newFunction)
{
  mChemEqI.setReversibility(rev);
  findAndSetFunction(newFunction);
}

void CReactionInterface::reverse(bool rev, const std::string & newFunction)
{
  mChemEqI.setReversibility(rev);
  mChemEqI.reverse();
  findAndSetFunction(newFunction);
}

void CReactionInterface::findAndSetFunction(const std::string & newFunction)
{
  std::vector<std::string> fl = getListOfPossibleFunctions();
  unsigned C_INT32 i, imax = fl.size();

  //no valid function?
  if (imax == 0)
    {
      setFunctionAndDoMapping("");
      return;
    }

  //first try the function provided as argument
  if (newFunction != "")
    {
      for (i = 0; i < imax; ++i)
        if (fl[i] == newFunction)
          {
            setFunctionAndDoMapping(fl[i]);
            return;
          }
    }

  //next try if the current function works
  std::string currentFunctionName = getFunctionName();

  if ("" != currentFunctionName)
    {
      for (i = 0; i < imax; ++i)
        if (fl[i] == currentFunctionName)
          {
            setFunctionAndDoMapping(fl[i]);
            return;
          }
    }

  // if not found then see if there is a best match in the list (i.e. a corresponding rev/irrev function).
  //if there is a current function try to find a related new function
  std::string s;

  if (currentFunctionName != "")
    {
      s = currentFunctionName.substr(0, currentFunctionName.find('(') - 1);

      //'-1' so as to strip off the white space before '('
      //This is purely heuristics
      for (i = 0; i < imax; i++)
        {
          if (fl[i].find(s) != std::string::npos)    // if find succeeds, the return value is likely to be 0
            //if (fl[i].find(s) >= 0) - for some reason this doesn't work
            {
              setFunctionAndDoMapping(fl[i]);
              return;
            }
        }
    }

  //try mass action next
  s = "Mass action";

  for (i = 0; i < imax; i++)
    {
      if (fl[i].find(s) != std::string::npos)    // if find succeeds, the return value is likely to be 0
        //if (fl[i].find(s) >= 0) - for some reason this doesn't work
        {
          setFunctionAndDoMapping(fl[i]);
          return;
        }
    }

  //try constant flux next
  s = "Constant flux";

  for (i = 0; i < imax; i++)
    {
      if (fl[i].find(s) != std::string::npos)    // if find succeeds, the return value is likely to be 0
        //if (fl[i].find(s) >= 0) - for some reason this doesn't work
        {
          setFunctionAndDoMapping(fl[i]);

          // If we have a reaction of the type X + Y = and the function
          // is Constant flux (reversible) we need to set the parameter v to
          // be negative to avoid negative concentrations during time course simulations.
          // Note, this fix is only done when we are assigning a default rate law.
          if (mChemEqI.getReversibility() == true &&
              mChemEqI.getListOfDisplayNames(CFunctionParameter::PRODUCT).size() == 0)
            {
              C_FLOAT64 v = - fabs(getLocalValue(0));
              setLocalValue(0, v);
            }

          return;
        }
    }

  //if everything else fails just take the first function from the list
  //this should not be reached since constant flux is a valid kinetic function for every reaction
  setFunctionAndDoMapping(fl[0]);
}

void CReactionInterface::connectFromScratch(CFunctionParameter::Role role)
{
  unsigned C_INT32 i, imax = mpParameters->getNumberOfParametersByUsage(role);

  if (!imax) return;

  // get the list of chem eq elements
  std::vector<std::string> el = getExpandedMetabList(role);

  // get the first parameter with the respective role
  CFunctionParameter::DataType Type;
  unsigned C_INT32 pos = 0;
  Type = mpParameters->getParameterByUsage(role, pos)->getType();

  if (Type == CFunctionParameter::VFLOAT64)
    {
      mNameMap[pos - 1] = el;
    }
  else if (Type == CFunctionParameter::FLOAT64)
    {
      if (el.size() > 0)
        mNameMap[pos - 1][0] = el[0];
      else
        {mNameMap[pos - 1][0] = "unknown"; /*mValid = false;*/}

      for (i = 1; i < imax; ++i)
        {
          Type = mpParameters->getParameterByUsage(role, pos)->getType();

          if (Type != CFunctionParameter::FLOAT64) fatalError();

          if (el.size() > i)
            mNameMap[pos - 1][0] = el[i];
          else
            {mNameMap[pos - 1][0] = "unknown"; /*mValid = false;*/}
        }
    }
  else fatalError();
}

bool CReactionInterface::isLocked(unsigned C_INT32 index) const
{return isLocked(getUsage(index));}

bool CReactionInterface::isLocked(CFunctionParameter::Role usage) const
{
  switch (usage)
    {
      case CFunctionParameter::MODIFIER:
        return false;
        break;

      case CFunctionParameter::TIME:
        return true;
        break;

      case CFunctionParameter::SUBSTRATE:
      case CFunctionParameter::PRODUCT:
      {
        // get number of parameters
        unsigned C_INT32 paramSize = mpParameters->getNumberOfParametersByUsage(usage);

        if (paramSize == 0)
          return true;

        // get index of first parameter
        unsigned C_INT32 pos = 0;
        mpParameters->getParameterByUsage(usage, pos); --pos;

        if (isVector(pos))
          {
            assert(paramSize == 1);
            return true;
          }
        else
          {
            return (mChemEqI.getListOfDisplayNames(usage).size() == 1);
          }
      }
      break;

      case CFunctionParameter::PARAMETER:
        return mpModel->getModelValues().size() <= 1;
        break;

      case CFunctionParameter::VOLUME:
        return mpModel->getCompartments().size() <= 1;
        break;

      default:
        break;
    }

  return false;
}

std::set< const CCopasiObject * > CReactionInterface::getDeletedParameters() const
{
  std::set< const CCopasiObject * > ToBeDeleted;

  // We need to compare the current visible local parameter with the one stored
  // in the reaction.
  const CReaction * pReaction
  = dynamic_cast< CReaction *>(CCopasiRootContainer::getKeyFactory()->get(mReactionReferenceKey));

  if (pReaction == NULL)
    return ToBeDeleted;

  if (pReaction->getFunction() == NULL)
    return ToBeDeleted;

  const CFunctionParameters & OriginalParameters
  = pReaction->getFunction()->getVariables();

  C_INT32 j, jmax = size();
  C_INT32 i, imax = OriginalParameters.size();
  const CFunctionParameter * pParameter;

  for (i = 0; i < imax; ++i)
    {
      pParameter = OriginalParameters[i];

      if (pParameter->getUsage() == CFunctionParameter::PARAMETER &&
          pReaction->isLocalParameter(i))
        {
          const std::string & Name = pParameter->getObjectName();

          //find parameter with same name in current parameters
          for (j = 0; j < jmax; ++j)
            if (Name == getParameterName(j)) break;

          if (j < jmax && mIsLocal[j])
            continue;

          // The old parameter is either not found or is no longer local, i.e., it needs to
          // be added to values to be deleted.
          ToBeDeleted.insert(pReaction->getParameters().getParameter(Name));
        }
    }

  return ToBeDeleted;
}

void CReactionInterface::initMapping()
{
  mpParameters = new CFunctionParameters(mpFunction->getVariables());
  //make sure mpParameters is deleted! (e.g. in copyMapping())
  mNameMap.resize(size());
  mValues.resize(size());
  mIsLocal.resize(size());
  C_INT32 i, imax = size();

  for (i = 0; i < imax; ++i)
    {
      if (isVector(i))
        mNameMap[i].resize(0);
      else
        {
          mNameMap[i].resize(1);
          mNameMap[i][0] = "unknown";
        }

      if (getUsage(i) == CFunctionParameter::PARAMETER)
        mIsLocal[i] = true;
      else
        mIsLocal[i] = false;

      mValues[i] = 0.1;
    }
}

void CReactionInterface::copyMapping()
{
  if (!mpParameters) //nothing to copy
    {
      initMapping();
      return;
    }

  // save the old information
  CFunctionParameters *oldParameters = mpParameters;
  std::vector<std::vector<std::string> > oldMap = mNameMap;
  std::vector<C_FLOAT64> oldValues = mValues;
  std::vector<bool> oldIsLocal = mIsLocal;

  //create new mapping
  initMapping();

  //try to preserve as much information from the old mapping as possible
  C_INT32 j, jmax = oldParameters->size();
  C_INT32 i, imax = size();

  for (i = 0; i < imax; ++i)
    {
      //find parameter with same name in old parameters
      for (j = 0; j < jmax; ++j)
        if ((*oldParameters)[j]->getObjectName() == getParameterName(i)) break;

      if (j == jmax) continue;

      //see if usage matches
      if (getUsage(i) != (*oldParameters)[j]->getUsage()) continue;

      //see if vector property matches
      if (isVector(i) != ((*oldParameters)[j]->getType() == CFunctionParameter::VFLOAT64))
        continue;

      mIsLocal[i] = oldIsLocal[j];
      mValues[i] = oldValues[j];

      switch (getUsage(i))
        {
          case CFunctionParameter::SUBSTRATE:
          case CFunctionParameter::PRODUCT:
            //TODO: check with chemeq
            mNameMap[i] = oldMap[j];
            break;

          case CFunctionParameter::MODIFIER:
            //TODO: check existence?
            mNameMap[i] = oldMap[j];
            //TODO: add to chemeq
            break;

          case CFunctionParameter::PARAMETER:
          case CFunctionParameter::VOLUME:
          case CFunctionParameter::TIME:
            //TODO: check existence?
            mNameMap[i] = oldMap[j];
            break;

          default:
            break;
        }
    }

  pdelete(oldParameters);
}

void CReactionInterface::connectNonMetabolites()
{
  C_INT32 i, imax = size();

  for (i = 0; i < imax; ++i)
    {
      //only do something if the current mapping is "unknown"
      if (mNameMap[i].size())
        if (mNameMap[i][0] != "unknown")
          continue;

      if (isLocalValue(i))
        continue;

      switch (getUsage(i))
        {
          case CFunctionParameter::SUBSTRATE:
          case CFunctionParameter::PRODUCT:
          case CFunctionParameter::MODIFIER:
            break;

          case CFunctionParameter::PARAMETER:

            if (mpModel->getModelValues().size() == 1)
              mNameMap[i][0] = mpModel->getModelValues()[0]->getObjectName();

            break;

          case CFunctionParameter::VOLUME:
            //if (mpModel->getCompartments().size()==1)
            //  mNameMap[i][0] = mpModel->getCompartments()[0]->getObjectName();
          {
            const CCompartment* comp = mChemEqI.getCompartment();

            if (comp)
              mNameMap[i][0] = comp->getObjectName();
          }
          break;

          case CFunctionParameter::TIME:
            mNameMap[i][0] = mpModel->getObjectName();
            break;

          default:
            break;
        }
    }
}

void CReactionInterface::setFunctionWithEmptyMapping(const std::string & fn)
{
  if ((fn == "") || (fn == "undefined"))
    {clearFunction(); return;}

  //get the function
  mpFunction = dynamic_cast<CFunction *>
               (CCopasiRootContainer::getFunctionList()->findLoadFunction(fn));

  if (!mpFunction) fatalError();

  pdelete(mpParameters);
  initMapping(); //empty mapping
}

void CReactionInterface::setFunctionAndDoMapping(const std::string & fn)
{
  if ((fn == "") || (fn == "undefined"))
    {clearFunction(); return;}

  //get the function
  mpFunction = dynamic_cast<CFunction *>
               (CCopasiRootContainer::getFunctionList()->findLoadFunction(fn));

  if (!mpFunction) fatalError();

  copyMapping();
  connectNonMetabolites();

  //guess initial connections between metabs and function parameters
  connectFromScratch(CFunctionParameter::SUBSTRATE);
  connectFromScratch(CFunctionParameter::PRODUCT);
  connectFromScratch(CFunctionParameter::MODIFIER); // we can not be pedantic about modifiers
  // because modifiers are not taken into acount
  // when looking for suitable functions
}

void CReactionInterface::updateModifiersInChemEq()
{
  mChemEqI.clearModifiers();
  unsigned C_INT32 j, jmax = size();

  for (j = 0; j < jmax; ++j)
    if (getUsage(j) == CFunctionParameter::MODIFIER) //all the modifiers in the table
      if (getMapping(j) != "unknown")
        mChemEqI.addModifier(getMapping(j));
}

void CReactionInterface::setMapping(unsigned C_INT32 index, std::string mn)
{
  mIsLocal[index] = false;

  switch (getUsage(index))
    {
      case CFunctionParameter::VOLUME:
      case CFunctionParameter::PARAMETER:
      case CFunctionParameter::TIME:
        assert(!isVector(index));
        mNameMap[index][0] = mn;
        break;

      case CFunctionParameter::SUBSTRATE:
      case CFunctionParameter::PRODUCT:

        if (isVector(index))
          mNameMap[index].push_back(mn);
        else
          {
            mNameMap[index][0] = mn;

            //TODO: check the following
            // if we have two parameters of this usage change the other one.
            unsigned C_INT32 listSize = mChemEqI.getListOfDisplayNames(getUsage(index)).size();

            if ((listSize == 2) && (mpParameters->getNumberOfParametersByUsage(getUsage(index)) == 2))
              {
                // get index of other parameter
                unsigned C_INT32 pos = 0;
                mpParameters->getParameterByUsage(getUsage(index), pos);

                if ((pos - 1) == index) mpParameters->getParameterByUsage(getUsage(index), pos);

                --pos;

                // get name if other metab
                std::vector<std::string> ml = getListOfMetabs(getUsage(index));
                std::string otherMetab;
                if (ml[0] == mn) otherMetab = ml[1]; else otherMetab = ml[0];

                // set other parameter
                mNameMap[pos][0] = otherMetab;
              }
          }

        break;

      case CFunctionParameter::MODIFIER:
        assert(!isVector(index));
        mNameMap[index][0] = mn;
        updateModifiersInChemEq();
        break;

      default:
        break;
    }
}

std::vector<std::string> CReactionInterface::getExpandedMetabList(CFunctionParameter::Role role) const
{
  const std::vector<std::string> & names = mChemEqI.getListOfDisplayNames(role);
  const std::vector<C_FLOAT64> & mults = mChemEqI.getListOfMultiplicities(role);

  unsigned C_INT32 j, jmax;
  unsigned C_INT32 i, imax = names.size();

  std::vector<std::string> ret;

  for (i = 0; i < imax; ++i)
    {
      if (role == CFunctionParameter::MODIFIER) jmax = 1;
      else jmax = (C_INT32)mults[i];

      for (j = 0; j < jmax; ++j)
        ret.push_back(names[i]);
    }

  return ret;
}

bool CReactionInterface::createMetabolites()
{
  bool created = mChemEqI.createNonExistingMetabs();

  // Update the parameter mapping to assure that the new names match.
  if (created)
    setFunctionAndDoMapping(getFunctionName());

  return created;
}

bool CReactionInterface::createOtherObjects() const
{
  bool ret = false;

  unsigned C_INT32 i, imax = size();

  for (i = 0; i < imax; ++i)
    {
      switch (getUsage(i))
        {
          case CFunctionParameter::VOLUME:

            if (mNameMap[i][0] == "unknown" || mNameMap[i][0] == "") break;

            if (mpModel->createCompartment(mNameMap[i][0], 1.0))
              ret = true;

            break;

          case CFunctionParameter::PARAMETER:

            if (mNameMap[i][0] == "unknown" || mNameMap[i][0] == "") break;

            if (!isLocalValue(i))
              if (mpModel->createModelValue(mNameMap[i][0], 1.0))
                ret = true;

            break;

          default:
            break;
        }
    }

  return ret;
}

bool CReactionInterface::isMulticompartment() const
{
  return mChemEqI.isMulticompartment();
}

bool CReactionInterface::isValid() const
{
  if (!mpFunction) return false;

  if (mpFunction->getObjectName() == "undefined") return false;

  //A reaction is invalid if it has a metab, a global parameter, or a compartment "unknown"
  unsigned C_INT j, jmax = size();

  for (j = 0; j < jmax; ++j)
    if ((mNameMap[j][0] == "unknown") && (!mIsLocal[j]))
      return false;

  return true;
}

#ifdef COPASI_DEBUG
void CReactionInterface::printDebug() const
{
  std::cout << "Reaction interface   " << std::endl;
  std::cout << " Reaction:   " << getReactionName() << std::endl;
  std::cout << "  Function: " << getFunctionName() << std::endl;
  std::cout << "  ChemEq:   " << getChemEqString() << std::endl;

  unsigned C_INT32 i, imax = size();

  for (i = 0; i < imax; ++i)
    {
      std::cout << "    ---  " << i << ": " << getParameterName(i)
                << ", vector: " << isVector(i) << " local: " << isLocalValue(i)
                << ", value: " << mValues[i] << std::endl;

      unsigned C_INT32 j, jmax = mNameMap[i].size();

      for (j = 0; j < jmax; ++j)
        std::cout << "            " << mNameMap[i][j] << std::endl;
    }

  std::cout << std::endl;
}
#endif // COPASI_DEBUG
