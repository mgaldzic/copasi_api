// Begin CVS Header
//   $Source: /fs/turing/cvs/copasi_dev/copasi/tss/CODEExporterC.cpp,v $
//   $Revision: 1.11 $
//   $Name: Build-33 $
//   $Author: shoops $
//   $Date: 2010/07/16 19:03:28 $
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

#include <locale>
#include <math.h>
#include "copasi.h"

#include "CopasiDataModel/CCopasiDataModel.h"
#include "report/CCopasiRootContainer.h"

#include "CODEExporterC.h"

#include "model/CModel.h"
#include "model/CMetab.h"
#include "model/CMetabNameInterface.h"
#include "utilities/CCopasiVector.h"
#include "model/CReaction.h"
#include "model/CMoiety.h"
#include "model/CChemEqElement.h"
#include "function/CExpression.h"
#include "function/CFunction.h"
#include "function/CFunctionDB.h"
#include "report/CKeyFactory.h"
#include "function/CEvaluationTree.h"
#include "function/CEvaluationNode.h"
#include "function/CEvaluationNodeObject.h"
#include "function/CEvaluationNodeOperator.h"
#include "function/CEvaluationNodeFunction.h"
#include "utilities/CCopasiTree.h"

#include "trajectory/CTrajectoryTask.h"
#include "trajectory/CTrajectoryProblem.h"

#include <iostream>
#include <fstream>
#include <ctype.h>

/**
 ** Constructor for the exporter.
 */
CODEExporterC::CODEExporterC()
{}

bool CODEExporterC::exportTitleData(const CModel* copasiModel, std::ofstream & outFile)
{

  unsigned C_INT32 metab_size = copasiModel->getMetabolitesX().size();
  unsigned C_INT32 indep_size = copasiModel->getNumIndependentReactionMetabs();
  unsigned C_INT32 ode_size = copasiModel->getNumODEMetabs();
  unsigned C_INT32 comps_size = copasiModel->getCompartments().size();
  unsigned C_INT32 modvals_size = copasiModel->getModelValues().size();
  unsigned C_INT32 reacs_size = copasiModel->getReactions().size();

  unsigned C_INT32 i, count;
  const CCopasiVector< CReaction > & reacs = copasiModel->getReactions();
  CReaction* reac;

  count = 0;

  for (i = 0; i < reacs_size; ++i)
    {
      reac = reacs[i];
      count = count + reac->getParameters().size();
    }

  outFile << "#ifdef SIZE_DEFINITIONS" << std::endl;
  outFile << "#define N_METABS " << metab_size << std::endl;
  outFile << "#define N_ODE_METABS " << ode_size << std::endl;
  outFile << "#define N_INDEP_METABS " << indep_size << std::endl;
  outFile << "#define N_COMPARTMENTS " << comps_size << std::endl;
  outFile << "#define N_GLOBAL_PARAMS " << modvals_size << std::endl;
  outFile << "#define N_KIN_PARAMS " << count << std::endl;
  outFile << "#define N_REACTIONS " << reacs_size << std::endl;

  outFile << "#endif // SIZE_DEFINITIONS" << std::endl;
  outFile << std::endl;

  outFile << "#ifdef TIME" << std::endl;
  outFile << "#define T  <set here a user name for the time variable> " << std::endl;
  outFile << "#endif // TIME" << std::endl;

  return true;
}

void CODEExporterC::setReservedNames()
{return;}  // TODO

std::string CODEExporterC::translateTimeVariableName()
{return "T";}

std::string CODEExporterC::setExportName(const CModelEntity::Status & status, unsigned C_INT32 n[], unsigned C_INT32 dependent)
{
  std::ostringstream name;

  switch (status)
    {
      case CModelEntity::FIXED:
        name << "p[" << n[0] << "]";
        n[0] ++;
        break;
      case CModelEntity::REACTIONS:

        if (!dependent)
          {
            name << "x[" << n[1] << "]";
            n[1] ++;
          }
        else
          {
            name << "y[" << n[2] << "]";
            n[2] ++;
          }

        break;
      case CModelEntity::ODE:
        name << "x[" << n[1] << "]";
        n[1] ++;
        break;
      case CModelEntity::ASSIGNMENT:
        name << "y[" << n[2] << "]";
        n[2] ++;
        break;
      default:
        return " ";
        break;
    }

  return name.str();
}

std::string CODEExporterC::setConcentrationName(const CModelEntity::Status & status, unsigned C_INT32 n[], unsigned C_INT32 dependent)
{
  std::ostringstream name;

  switch (status)
    {
      case CModelEntity::FIXED:
        name << "p_c[" << n[0] << "]";
        n[0] ++;
        break;
      case CModelEntity::REACTIONS:

        if (!dependent)
          {
            name << "x_c[" << n[1] << "]";
            n[1] ++;
          }
        else
          {
            name << "y_c[" << n[2] << "]";
            n[2] ++;
          }

        break;
      case CModelEntity::ODE:
        name << "x_c[" << n[1] << "]";
        n[1] ++;
        break;
      case CModelEntity::ASSIGNMENT:
        name << "y_c[" << n[2] << "]";
        n[2] ++;
        break;
      default:
        return " ";
        break;
    }

  return name.str();
}

/**
 **      This method adapt a Copasi name for C syntax:
 **      Names can not start with a number.
 **      Any other combination of letters and numbers is valid as is the underscore.
 **/
std::string CODEExporterC::translateObjectName(const std::string & realName)
{
  std::locale C("C");
  char ch;

  std::ostringstream tmpName;

  unsigned C_INT32 realName_size = realName.size();
  unsigned C_INT32 i;

  ch = realName[0];

  if (!std::isalpha(ch, C))
    {
      tmpName << "_";

      if (std::isdigit(ch, C)) tmpName << ch;
    }
  else tmpName << ch;

  for (i = 1; i < realName_size; i++)
    {
      ch = realName[i];

      if (std::isalpha(ch, C))
        {
          if (std::isspace(realName[i - 1], C) && std::islower(ch, C))
            tmpName << (char) toupper(ch);
          else
            tmpName << ch;
        }

      if (std::isdigit(ch, C)) tmpName << ch;

      if (std::ispunct(ch, C))
        switch (ch)
          {
            case '_':
              tmpName << ch;
              break;
            case '-':
            case '{':
            case '}':
            case '(':
            case ')':
            case '[':
            case ']':
              tmpName << "_";
              break;
            default:
              break;
          }
    }

  return testName(tmpName.str());
}
/**
 **      This method tests whether the given C name already assigned,
 **      put the new name (in cappital letters:
 **      all names can be upper or lower case)
 **      in the set of assigned names
 **      or  modify the name
 **/

std::string CODEExporterC::testName(const std::string & name)
{
  std::locale C("C");
  char ch;

  std::ostringstream newname, tmp;

  unsigned C_INT32 name_size = name.size();
  unsigned C_INT32 i;

  for (i = 0; i < name_size; i++)
    {
      ch = name[i];

      if (std::isalpha(ch, C) && std::islower(ch, C))
        tmp << (char) toupper(ch);
      else
        tmp << ch;
    }

  if (NameSet.find(tmp.str()) == NameSet.end())
    {
      NameSet.insert(tmp.str());
      Frequancy[tmp.str()] = 0;

      return name;
    }
  else
    {
      Frequancy[tmp.str()]++;
      newname << name << "_" << Frequancy[tmp.str()];

      return testName(newname.str());
    }
}

std::string CODEExporterC::setODEName(const std::string & objName)
{
  return "d" + objName;
}

bool CODEExporterC::preprocess(const CModel* copasiModel)

{
  unsigned C_INT32 n[3] = {0, 0, 0};
  unsigned C_INT32 n_c[3] = {0, 0, 0};
  unsigned C_INT32 i, j;
  unsigned C_INT32 dependent;

  setReservedNames();

  NameMap[timeKey] = translateTimeVariableName();

  const CCopasiVector< CMetab > & metabs = copasiModel->getMetabolitesX();
  unsigned C_INT32 metabs_size = metabs.size();

  for (i = 0; i < metabs_size; i++)
    {
      CMetab * metab = metabs[i];

      //if (metab->isUsed())
      {
        std::string smname;
        std::string name;
        dependent = metab->isDependent();

        smname = setExportName(metab->getStatus(), n, dependent);
        name = setConcentrationName(metab->getStatus(), n_c, dependent);

        NameMap[metab->getKey()] = name;

        std::ostringstream smKey;
        smKey << "sm_" << metab->getKey();

        NameMap[smKey.str()] = smname;

        if ((metab->getStatus() == CModelEntity::REACTIONS && !metab->isDependent()) || metab->getStatus() == CModelEntity::ODE)
          {
            std::ostringstream odeKey;
            odeKey << "ode_" << metab->getKey();

            NameMap[odeKey.str()] = setODEName(smname);
          }
      }
    }

  unsigned C_INT32 comps_size = copasiModel->getCompartments().size();
  const CCopasiVector< CCompartment > & comps = copasiModel->getCompartments();

  for (i = 0; i < comps_size; i++)
    {
      CCompartment * comp = comps[i];

      std::string name;
      dependent = 0;

      name = setExportName(comp->getStatus(), n, dependent);
      NameMap[comp->getKey()] = name;

      if (comp->getStatus() == CModelEntity::ODE)
        {
          std::ostringstream odeKey;
          odeKey << "ode_" << comp->getKey();

          NameMap[odeKey.str()] = setODEName(name);
        }
    }

  unsigned C_INT32 modvals_size = copasiModel->getModelValues().size();
  const CCopasiVector< CModelValue > & modvals = copasiModel->getModelValues();

  for (i = 0; i < modvals_size; i++)
    {
      CModelValue* modval = modvals[i];
      std::string name = setExportName(modval->getStatus(), n, 0);
      NameMap[modval->getKey()] = name;

      if (modval->getStatus() == CModelEntity::ODE)
        {
          std::ostringstream odeKey;
          odeKey << "ode_" << modval->getKey();
          NameMap[odeKey.str()] = setODEName(name);
        }
    }

  unsigned C_INT32 reacs_size = copasiModel->getReactions().size();

  const CCopasiVector< CReaction > & reacs = copasiModel->getReactions();

  std::set<std::string> tmpset;

  for (i = 0; i < reacs_size; ++i)
    {
      unsigned C_INT32 params_size;

      params_size = reacs[i]->getParameters().size();

      for (j = 0; j < params_size; ++j)
        {
          std::ostringstream name;

          name << "p[" << n[0] << "]";
          n[0] ++;

          NameMap[reacs[i]->getParameters().getParameter(j)->getKey()] = name.str();
        }

      const CFunction* func = reacs[i]->getFunction();

      std::string name = func->getObjectName();

      if (func->getRoot())
        setExportNameOfFunction(func->getRoot(), tmpset);

      if (func->getType() != CEvaluationTree::MassAction)
        if (tmpset.find(name) == tmpset.end())
          {
            NameMap[func->getKey()] = translateObjectName(name);
            tmpset.insert(name);
          }
    }

  return true;
}

void CODEExporterC::setExportNameOfFunction(const CEvaluationNode* pNode, std::set<std::string> & tmpset)
{
  if (pNode)
    {
      CFunctionDB* pFunctionDB = CCopasiRootContainer::getFunctionList();
      CCopasiTree<CEvaluationNode>::const_iterator treeIt = pNode;

      while (treeIt != NULL)
        {
          if (CEvaluationNode::type(treeIt->getType()) == CEvaluationNode::CALL)
            {
              const CFunction* ifunc;
              ifunc = static_cast<CFunction*>(pFunctionDB->findFunction((*treeIt).getData()));

              setExportNameOfFunction(ifunc->getRoot(), tmpset);

              if (ifunc->getType() != CEvaluationTree::MassAction)
                if (tmpset.find(ifunc->getObjectName()) == tmpset.end())
                  {
                    NameMap[ifunc->getKey()] = translateObjectName(ifunc->getObjectName());
                    tmpset.insert(ifunc->getObjectName());
                  }
            }

          ++treeIt;
        }
    }
}

bool CODEExporterC::exportSingleObject(std::ostringstream & which, std::string & name, std::string & expression, std::string & comments)
{
  which << name << " = " << expression << ";"
  << '\t' << "//" << comments << std::endl;

  return true;
}

bool CODEExporterC::exportSingleMetabolite(const CMetab* metab, std::string & expression, std::string & comments)
{

  std::string name;

  std::ostringstream smKey;
  smKey << "sm_" << metab->getKey();
  name = NameMap[smKey.str()];

  switch (metab->getStatus())
    {

      case CModelEntity::FIXED:

        if (!exportSingleObject(fixed, name, expression, comments))
          return false;

        break;
      case CModelEntity::REACTIONS:
      case CModelEntity::ODE:
      {
        if (metab->isDependent())
          {
            if (!exportSingleObject(assignment, name, expression, comments))
              return false;
          }
        else
          {
            if (!exportSingleObject(initial, name, expression, comments))
              return false;
          }

        break;
      }
      case CModelEntity::ASSIGNMENT:
      {
        if (!exportSingleObject(assignment, name, expression, comments))
          return false;

        break;
      }
      default:
        return false;
        break;
    }

  return true;
}

bool CODEExporterC::exportSingleCompartment(const CCompartment* comp, std::string & expression, std::string & comments)
{
  switch (comp->getStatus())
    {
      case CModelEntity::FIXED:
      {
        if (!exportSingleObject(fixed, NameMap[comp->getKey()], expression, comments))
          return false;

        break;
      }
      case CModelEntity::ODE:
      {
        if (!exportSingleObject(initial, NameMap[comp->getKey()], expression, comments))
          return false;

        break;
      }
      case CModelEntity::ASSIGNMENT:
      {
        if (!exportSingleObject(assignment, NameMap[comp->getKey()], expression, comments))
          return false;

        break;
      }
      default:
        return false;
        break;
    }

  return true;
}

bool CODEExporterC::exportSingleModVal(const CModelValue* modval, std::string & expression, std::string & comments)
{
  switch (modval->getStatus())
    {
      case CModelEntity::FIXED:
      {
        if (!exportSingleObject(fixed, NameMap[modval->getKey()], expression, comments))
          return false;

        break;
      }
      case CModelEntity::ODE:
      {
        if (!exportSingleObject(initial, NameMap[modval->getKey()], expression, comments))
          return false;

        break;
      }
      case CModelEntity::ASSIGNMENT:
      {
        if (!exportSingleObject(assignment, NameMap[modval->getKey()], expression, comments))
          return false;

        break;
      }
      default:
        return false;
        break;
    }

  return true;
}

bool CODEExporterC::exportSingleModelEntity(const CModelEntity* tmp, std::string & expression, std::string & comments)
{

  std::string name;

  const CMetab* metab;
  metab = dynamic_cast< const CMetab * >(tmp);

  if (metab)
    {
      std::ostringstream smKey;
      smKey << "sm_" << metab->getKey();
      name = NameMap[smKey.str()];
    }
  else
    name = NameMap[tmp->getKey()];

  switch (tmp->getStatus())
    {
      case CModelEntity::FIXED:
      {
        if (!exportSingleObject(fixed, name, expression, comments))
          return false;

        break;
      }
      case CModelEntity::ODE:
      {
        if (!exportSingleObject(initial, name, expression, comments))
          return false;

        break;
      }
      case CModelEntity::ASSIGNMENT:
      {
        if (!exportSingleObject(assignment, name, expression, comments))
          return false;

        break;
      }
      default:
        return false;
        break;
    }

  return true;
}

bool CODEExporterC::exportSingleParameter(const CCopasiParameter* param, std::string & expression, std::string & comments)
{
  if (!exportSingleObject(fixed, NameMap[param->getKey()], expression, comments)) return false;

  return true;
}

bool CODEExporterC::exportKineticFunction(CReaction* /* reac */)
{return true;}

bool CODEExporterC::exportKineticFunctionGroup(const CModel* copasiModel)
{
  const CCopasiVector< CReaction > & reacs = copasiModel->getReactions();
  unsigned C_INT32 reacs_size = reacs.size();
  CReaction* reac;

  std::set<std::string> isExported;

  unsigned C_INT32 i;

  for (i = 0; i < reacs_size; ++i)
    {
      reac = reacs[i];
      const CFunction* func = reac->getFunction();

      if (func->getRoot())
        findFunctionsCalls(func->getRoot(), isExported);

      if (func->getType() != CEvaluationTree::MassAction)
        {
          if (!exportSingleFunction(func, isExported)) return false;
        }
    }

  return true;
}

void CODEExporterC::findFunctionsCalls(const CEvaluationNode* pNode, std::set<std::string>& isExported)
{
  if (pNode)
    {
      CFunctionDB* pFunctionDB = CCopasiRootContainer::getFunctionList();
      CCopasiTree<CEvaluationNode>::const_iterator treeIt = pNode;

      while (treeIt != NULL)
        {
          if (CEvaluationNode::type(treeIt->getType()) == CEvaluationNode::CALL)
            {
              const CFunction* ifunc;
              ifunc = static_cast<CFunction*>(pFunctionDB->findFunction((*treeIt).getData()));

              findFunctionsCalls(ifunc->getRoot(), isExported);

              if (ifunc->getType() != CEvaluationTree::MassAction)
                {
                  if (!exportSingleFunction(ifunc, isExported)) return;
                }
            }

          ++treeIt;
        }
    }
}

bool CODEExporterC::exportSingleFunction(const CFunction *func, std::set<std::string>& isExported)
{
  CFunctionDB* pFunctionDB = CCopasiRootContainer::getFunctionList();

  CFunction* tmpfunc = NULL;
  tmpfunc = new CFunction(*func);

  if (func->getType() != CEvaluationTree::MassAction)
    {

      CCopasiTree< CEvaluationNode>::iterator treeIt = tmpfunc->getRoot();
      CCopasiTree< CEvaluationNode>::iterator newIt = treeIt;

      unsigned C_INT32 j, varbs_size = tmpfunc->getVariables().size();
      std::map< std::string, std::string > parameterNameMap;
      std::set<std::string> parameterNameSet;

      std::map< CFunctionParameter::Role, std::string > constName;
      std::map< CFunctionParameter::Role, unsigned C_INT32 > tmpIndex;

      constName[CFunctionParameter::SUBSTRATE] = "sub_"; tmpIndex[CFunctionParameter::SUBSTRATE] = 0;
      constName[CFunctionParameter::PRODUCT] = "prod_"; tmpIndex[CFunctionParameter::PRODUCT] = 0;
      constName[CFunctionParameter::PARAMETER] = "param_"; tmpIndex[CFunctionParameter::PARAMETER] = 0;
      constName[CFunctionParameter::MODIFIER] = "modif_"; tmpIndex[CFunctionParameter::MODIFIER] = 0;
      constName[CFunctionParameter::VOLUME] = "volume_"; tmpIndex[CFunctionParameter::VOLUME] = 0;
      constName[CFunctionParameter::VARIABLE] = "varb_"; tmpIndex[CFunctionParameter::VARIABLE] = 0;
      constName[CFunctionParameter::TIME] = "time_"; tmpIndex[CFunctionParameter::VARIABLE] = 0;

      for (j = 0; j < varbs_size; ++j)
        {
          if (parameterNameSet.find(tmpfunc->getVariables()[j]->getObjectName()) == parameterNameSet.end())
            {
              std::ostringstream tmpName;
              CFunctionParameter::Role role = tmpfunc->getVariables()[j]->getUsage();

              tmpName << constName[role] << tmpIndex[role];
              parameterNameMap[ tmpfunc->getVariables()[j]->getObjectName()] = tmpName.str();
              parameterNameSet.insert(tmpfunc->getVariables()[j]->getObjectName());
              tmpIndex[role]++;
            }
        }

      CODEExporter::modifyTreeForMassAction(tmpfunc);

      while (newIt != NULL)
        {
          if (CEvaluationNode::type(newIt->getType()) == CEvaluationNode::VARIABLE)
            {
              newIt->setData(parameterNameMap[ tmpfunc->getVariables()[newIt->getData()]->getObjectName()]);
            }

          if (CEvaluationNode::type(newIt->getType()) == CEvaluationNode::CALL)
            {
              const CFunction* callfunc;
              callfunc = static_cast<CFunction*>(pFunctionDB->findFunction((*newIt).getData()));

              if (callfunc->getType() != CEvaluationTree::MassAction)
                newIt->setData(NameMap[callfunc->getKey()]);
            }

          ++newIt;
        }

      std::string name = func->getObjectName();

      if (isExported.find(name) == isExported.end())
        {
          unsigned C_INT32 j, varbs_size = tmpfunc->getVariables().size();

          functions << "double " << NameMap[func->getKey()] << "(";
          headers << "double " << NameMap[func->getKey()] << "(";

          for (j = 0; j < varbs_size; ++j)
            {
              functions << "double " << parameterNameMap[ tmpfunc->getVariables()[j]->getObjectName().c_str()];

              if (j != varbs_size - 1) functions << ", ";

              headers << "double " << parameterNameMap[ tmpfunc->getVariables()[j]->getObjectName().c_str()];

              if (j != varbs_size - 1) headers << ", ";
            }

          functions << ") ";
          functions << '\t' << "//" << name << std::endl;
          functions << "{return  " << tmpfunc->getRoot()->getDisplay_C_String(tmpfunc).c_str() << ";} " << std::endl;

          headers << "); " << std::endl;

          isExported.insert(name);
        }
    }

  return true;
}

std::string CODEExporterC::KineticFunction2ODEmember(const CReaction *reac)
{
  std::ostringstream equation;

  if (reac->getFunction()->getType() != CEvaluationTree::MassAction)
    {
      const CFunctionParameters & params = reac->getFunctionParameters();
      unsigned C_INT32 k, params_size = params.size();
      const std::vector<std::vector<std::string> > & keyMap = reac->getParameterMappings();

      equation << NameMap[reac->getFunction()->getKey()] << "(";

      for (k = 0; k < params_size; ++k)
        {
          CFunctionParameter::Role role = params[k]->getUsage();

          CCopasiObject * obj = CCopasiRootContainer::getKeyFactory()->get(keyMap[k][0]);

          if ((role == CFunctionParameter::SUBSTRATE)
              || (role == CFunctionParameter::PRODUCT)
              || (role == CFunctionParameter::MODIFIER))
            {
              CMetab* metab;
              metab = dynamic_cast< CMetab * >(obj);
              equation << NameMap[metab->getKey()];
            }

          if (role == CFunctionParameter::PARAMETER)
            {
              if (!(reac->isLocalParameter(k)))
                {
                  CModelValue* modval;
                  modval = dynamic_cast< CModelValue * >(obj);
                  equation << NameMap[modval->getKey()];
                }
              else
                {
                  CCopasiParameter* param;
                  param = dynamic_cast< CCopasiParameter * >(obj);
                  equation << NameMap[param->getKey()];
                }
            }

          if (role == CFunctionParameter::VOLUME)
            {
              CCompartment* comp;
              comp = dynamic_cast< CCompartment * >(obj);
              equation << NameMap[comp->getKey()];
            }

          if (role == CFunctionParameter::TIME)
            {
              equation << "T";
            }

          if (k != params_size - 1) equation << ", ";
        }

      equation << ")";
    }
  else
    {
      const CCopasiVector<CChemEqElement> & substrs = reac->getChemEq().getSubstrates();
      const CCopasiVector<CChemEqElement> & prods = reac->getChemEq().getProducts();
      const std::vector<std::vector<std::string> > & keyMap = reac->getParameterMappings();
      CCopasiObject * obj;

      unsigned C_INT32 substrs_size = substrs.size(), prods_size = prods.size();
      unsigned C_INT32 k, m, mult;

      CChemEqElement* substr;
      CChemEqElement* prod;

      const CMassAction cMassAction = *static_cast<const CMassAction*>(reac->getFunction());

      equation << "(";

      obj = CCopasiRootContainer::getKeyFactory()->get(keyMap[0][0]);

      if (!(reac->isLocalParameter(0)))
        {
          CModelValue* modval;
          modval = dynamic_cast< CModelValue * >(obj);

          equation << NameMap[modval->getKey()];
        }
      else
        {
          CCopasiParameter* param;
          param = dynamic_cast< CCopasiParameter * >(obj);

          equation << NameMap[param->getKey()];
        }

      for (k = 0; k < substrs_size; ++k)
        {
          substr = substrs[k];
          mult = (unsigned C_INT32) substr->getMultiplicity();

          assert(substr->getMetabolite());
          equation << " * " << NameMap[substr->getMetabolite()->getKey()];

          if (mult > 1)
            for (m = 1; m < mult; ++m)
              equation << " * " << NameMap[substr->getMetabolite()->getKey()];
        }

      if (cMassAction.isReversible() == TriTrue)
        {
          equation << " - ";

          obj = CCopasiRootContainer::getKeyFactory()->get(keyMap[2][0]);

          if (!(reac->isLocalParameter(2)))
            {
              CModelValue* modval;
              modval = dynamic_cast< CModelValue * >(obj);

              equation << NameMap[modval->getKey()];
            }
          else
            {
              CCopasiParameter* param;
              param = dynamic_cast< CCopasiParameter * >(obj);

              equation << NameMap[param->getKey()];
            }

          for (k = 0; k < prods_size; ++k)
            {
              prod = prods[k];
              mult = (unsigned C_INT32) prod->getMultiplicity();

              assert(prod->getMetabolite());
              equation << " * " << NameMap[prod->getMetabolite()->getKey()];

              if (mult > 1)
                for (m = 1; m < mult; ++m)
                  equation << " * " << NameMap[prod->getMetabolite()->getKey()];
            }
        }

      equation << ") ";
    }

  return equation.str();
}

bool CODEExporterC::exportSingleODE(const CModelEntity* mentity, std::string & equation, std::string & comments)
{
  std::ostringstream odeKey;
  odeKey << "ode_" << mentity->getKey();

  if (!exportSingleObject(ode, NameMap[odeKey.str()], equation, comments)) return false;

  return true;
}

std::string CODEExporterC::exportTitleString(const unsigned C_INT32 tmp)
{
  switch (tmp)
    {
      case INITIAL:
        return "#ifdef INITIAL";
      case FIXED:
        return "#ifdef FIXED";
      case ASSIGNMENT:
        return "#ifdef ASSIGNMENT";
      case HEADERS:
        return "#ifdef  FUNCTIONS_HEADERS";
      case FUNCTIONS:
        return "#ifdef  FUNCTIONS";
      case ODEs:
        return "#ifdef ODEs";
      default:
        return " ";
    }
}

std::string CODEExporterC::getDisplayExpressionString(CExpression * tmp)
{
  std::string str1;

  str1 = tmp->getRoot()->getDisplay_C_String(tmp).c_str();

  return str1;
}

std::string CODEExporterC::exportClosingString(const unsigned C_INT32 tmp)
{
  switch (tmp)
    {
      case INITIAL:
        return "#endif INITIAL\n";
      case FIXED:
        return "#endif FIXED\n";
      case ASSIGNMENT:
        return "#endif ASSIGNMENT\n";
      case HEADERS:
        return "#endif FUNCTIONS_HEADERS\n";
      case FUNCTIONS:
        return "#endif FUNCTIONS\n";
      case ODEs:
        return "#endif ODEs\n";
      default:
        return " ";
    }
}
