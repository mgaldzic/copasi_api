// Begin CVS Header
//   $Source: /fs/turing/cvs/copasi_dev/copasi/model/CMMLOutput.cpp,v $
//   $Revision: 1.8 $
//   $Name: Build-33 $
//   $Author: shoops $
//   $Date: 2010/07/16 19:00:59 $
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

#include <sstream>

#include "CMMLOutput.h"
#include "CopasiDataModel/CCopasiDataModel.h"
#include "model/CModel.h"
#include "function/CExpression.h"
#include "report/CKeyFactory.h"
#include "copasi/report/CCopasiRootContainer.h"

CMMLOutput::CMMLOutput()
{}

#include "utilities/copasimathml.h"

void CMMLOutput::writeLHS(std::ostream & out,
                          const std::string & metabName,
                          const std::string & compName, unsigned C_INT32 l)
{
  out << SPC(l + 0) << "<mfrac>" << std::endl;
  out << SPC(l + 1) << "<mrow>" << std::endl;
  out << SPC(l + 2) << "<mo>d</mo>" << std::endl;
  out << SPC(l + 2) << "<mfenced>" << std::endl;
  out << SPC(l + 3) << "<mrow>" << std::endl;
  //out << SPC(l + 4) << "<mi>" << CMathMl::fixName(metabName) << "</mi>" << std::endl;
  out << SPC(l + 4) << "<mi>[" << metabName << "]</mi>" << std::endl;
  out << SPC(l + 4) << "<mo>" << "&CenterDot;" << "</mo>" << std::endl;
  out << SPC(l + 4) << "<msub><mi>V</mi><mi>" << CMathMl::fixName(compName) << "</mi></msub>" << std::endl;

  out << SPC(l + 3) << "</mrow>" << std::endl;

  out << SPC(l + 2) << "</mfenced>" << std::endl;
  out << SPC(l + 1) << "</mrow>" << std::endl;
  out << SPC(l + 1) << "<mrow>" << std::endl;
  out << SPC(l + 2) << "<mo>d</mo><mi>t</mi>" << std::endl;
  out << SPC(l + 1) << "</mrow>" << std::endl;
  out << SPC(l + 0) << "</mfrac>" << std::endl;
}

void CMMLOutput::writeRHS(std::ostream & out,
                          const CMetab* pMetab, const CReaction* pReac,
                          bool numbers, bool expand, bool expandFull,
                          unsigned C_INT32 l)
{
  if (!pMetab)
    {
      out << SPC(l + 0) << "Error: invalid metabolite" << std::endl;
      return;
    }

  if (!pReac)
    {
      out << SPC(l + 0) << "Error: invalid reaction" << std::endl;
      return;
    }

  const CCopasiVector < CChemEqElement > & balances = pReac->getChemEq().getBalances();

  C_FLOAT64 balance = 0;

  unsigned C_INT32 i, imax = balances.size();

  for (i = 0; i < imax; ++i)
    {
      if (balances[i]->getMetaboliteKey() == pMetab->getKey())
        {
          balance = balances[i]->getMultiplicity();
          break;
        }
    }

  if (0 == balance) return;

  out << SPC(l + 0) << "<mrow>" << std::endl;

  //balance
  if (balance == 1.0)
    {
      out << SPC(l + 1) << "<mo>" << "+" << "</mo>" << std::endl;
    }
  else if (balance == -1.0)
    {
      out << SPC(l + 1) << "<mo>" << "-" << "</mo>" << std::endl;
    }
  else if (balance < 0.0)
    {
      out << SPC(l + 1) << "<mo>" << "-" << "</mo><mn>" << -balance << "</mn>"
      << "<mo>" << "&CenterDot;" << "</mo>" << std::endl;
    }
  else // >0
    {
      out << SPC(l + 1) << "<mo>" << "+" << "</mo><mn>" << balance << "</mn>"
      << "<mo>" << "&CenterDot;" << "</mo>" << std::endl;
    }

  //compartment volume?
  if (pReac->getCompartmentNumber() == 1)
    {
      std::string compName = pMetab->getCompartment()->getObjectName();
      out << SPC(l + 1) << "<msub><mi>V</mi><mi>" << CMathMl::fixName(compName)
      << "</mi></msub>" << std::endl;
      out << SPC(l + 1) << "<mo>" << "&CenterDot;" << "</mo>" << std::endl;
    }

  //kinetic function
  if (pReac->getFunction())
    {
      std::vector<std::vector<std::string> > params;
      createParameterMapping(pReac, params, numbers);

      if (expand) out << SPC(l + 1) << "<mfenced>" << std::endl;

      pReac->getFunction()->writeMathML(out, params, expand, expandFull, l + 1);

      if (expand) out << SPC(l + 1) << "</mfenced>" << std::endl;
    }

  out << SPC(l + 0) << "</mrow>" << std::endl;
}

void CMMLOutput::createParameterMapping(const CReaction* pReac,
                                        std::vector<std::vector<std::string> > & params,
                                        bool numbers)
{
  assert(pReac);
  assert(pReac->getFunction());

  const CFunctionParameters & functionParams = pReac->getFunctionParameters();
  unsigned C_INT32 j, jmax;
  unsigned C_INT32 i, imax = functionParams.size();
  params.resize(imax);

  for (i = 0; i < imax; ++i)
    {
      params[i].resize(1);

      std::string name;

      //std::ostringstream number;
      switch (functionParams[i]->getUsage())
        {
          case CFunctionParameter::SUBSTRATE:
          case CFunctionParameter::PRODUCT:
          case CFunctionParameter::MODIFIER:

            if (functionParams[i]->getType() == CFunctionParameter::FLOAT64)
              {
                name = CCopasiRootContainer::getKeyFactory()->get(pReac->getParameterMappings()[i][0])->getObjectDisplayName();
                //params[i][0] = "<mi>"+ CMathMl::fixName(name)+"</mi>";
                params[i][0] = "<mi>[" + name + "]</mi>";
              }
            else if (functionParams[i]->getType() == CFunctionParameter::VFLOAT64)
              {
                jmax = pReac->getParameterMappings()[i].size();
                params[i].resize(jmax);

                for (j = 0; j < jmax; ++j)
                  {
                    name = CCopasiRootContainer::getKeyFactory()->get(pReac->getParameterMappings()[i][j])->getObjectDisplayName();
                    //params[i][j] = "<mi>"+ CMathMl::fixName(name)+"</mi>";
                    params[i][j] = "<mi>[" + name + "]</mi>";
                  }
              }
            else assert(false);

            break;

          case CFunctionParameter::PARAMETER:

            if (pReac->isLocalParameter(i))
              {
                if (numbers)
                  {
                    std::ostringstream number;
                    number << pReac->getParameterValue(functionParams[i]->getObjectName());
                    params[i][0] = "<mn>" + number.str() + "</mn>";
                  }
                else
                  {
                    name = CCopasiRootContainer::getKeyFactory()->get(pReac->getParameterMappings()[i][0])->getObjectName();
                    //params[i][0] = "<mi>" + CMathMl::fixName(name) + "</mi>";
                    params[i][0] = "<msub><mi>" + CMathMl::fixName(name) + "</mi><mi>("
                                   + CMathMl::fixName(pReac->getObjectName()) + ")</mi></msub>";
                  }
              }
            else
              {
                name = CCopasiRootContainer::getKeyFactory()->get(pReac->getParameterMappings()[i][0])->getObjectName();
                params[i][0] = "<mi>" + CMathMl::fixName(name) + "</mi>";
                //params[i][0] = "<mi>ggg</mi>";
              }

            break;

          case CFunctionParameter::VOLUME:
            name = CCopasiRootContainer::getKeyFactory()->get(pReac->getParameterMappings()[i][0])->getObjectName();
            params[i][0] = "<msub><mi>V</mi><mi>" + CMathMl::fixName(name)
                           + "</mi></msub>";
            break;

          case CFunctionParameter::TIME:
            params[i][0] = "<mi>time</mi>";
            break;

          default:
            break;
        }
    }
}

void CMMLOutput::writeLHS_ModelValue(std::ostream & out,
                                     const std::string & valueName, unsigned C_INT32 l)
{
  out << SPC(l + 0) << "<mfrac>" << std::endl;
  out << SPC(l + 1) << "<mrow>" << std::endl;
  out << SPC(l + 2) << "<mo>d</mo>" << std::endl;
  out << SPC(l + 2) << "<mi>" << CMathMl::fixName(valueName) << "</mi>" << std::endl;
  out << SPC(l + 1) << "</mrow>" << std::endl;
  out << SPC(l + 1) << "<mrow>" << std::endl;
  out << SPC(l + 2) << "<mo>d</mo><mi>t</mi>" << std::endl;
  out << SPC(l + 1) << "</mrow>" << std::endl;
  out << SPC(l + 0) << "</mfrac>" << std::endl;
}

void CMMLOutput::writeRHS_ModelEntity(std::ostream & out,
                                      const CModelEntity* pEntity,
                                      bool expandFull, unsigned C_INT32 l)
{
  if (!pEntity)
    {
      out << SPC(l + 0) << "Error: invalid model entity" << std::endl;
      return;
    }

  if (!pEntity->getExpressionPtr())
    {
      out << SPC(l + 0) << "Error: no expression" << std::endl;
      return;
    }

  out << SPC(l + 0) << "<mrow>" << std::endl;

  pEntity->getExpressionPtr()->writeMathML(out, expandFull, l + 1);

  out << SPC(l + 0) << "</mrow>" << std::endl;
}

void CMMLOutput::writeDifferentialEquations(std::ostream & mml, CModel * model, bool localParameterNumbers, bool expand, bool expandFull)
{
  bool hasContents = false;
  //mml.str("");

  unsigned C_INT32 l = 0;
  mml << SPC(l) << "<mtable>" << std::endl;

  //write equations for compartments
  C_INT32 i, imax = model->getCompartments().size();

  for (i = 0; i < imax; i++)
    {
      if (model->getCompartments()[i]->getStatus() == CModelEntity::ODE)
        {
          mml << SPC(l + 1) << "<mtr>" << std::endl;

          //first column (lhs)
          mml << SPC(l + 2) << "<mtd>" << std::endl;

          mml << SPC(l + 3) << "<mfrac>" << std::endl;
          mml << SPC(l + 4) << "<mrow>" << std::endl;
          mml << SPC(l + 5) << "<mo>d</mo>" << std::endl;

          mml << SPC(l + 5) << "<msub><mi>V</mi><mi>"
          << CMathMl::fixName(model->getCompartments()[i]->getObjectName())
          << "</mi></msub>" << std::endl;

          mml << SPC(l + 4) << "</mrow>" << std::endl;
          mml << SPC(l + 3) << "<mrow>" << std::endl;
          mml << SPC(l + 2) << "<mo>d</mo><mi>t</mi>" << std::endl;
          mml << SPC(l + 1) << "</mrow>" << std::endl;
          mml << SPC(l + 0) << "</mfrac>" << std::endl;

          mml << SPC(l + 2) << "</mtd>" << std::endl;

          //second column ("=")
          mml << SPC(l + 2) << "<mtd>" << std::endl;
          mml << SPC(l + 3) << "<mo>=</mo>" << std::endl;
          mml << SPC(l + 2) << "</mtd>" << std::endl;

          //third column (rhs)
          mml << SPC(l + 2) << "<mtd columnalign='left'>" << std::endl;

          writeRHS_ModelEntity(mml, model->getCompartments()[i],
                               expandFull, l + 3);
          mml << SPC(l + 2) << "</mtd>" << std::endl;

          mml << SPC(l + 1) << "</mtr>" << std::endl;
        }
      else if (model->getCompartments()[i]->getStatus() == CModelEntity::ASSIGNMENT)
        {
          mml << SPC(l + 1) << "<mtr>" << std::endl;

          //first column (lhs)
          mml << SPC(l + 2) << "<mtd>" << std::endl;
          mml << SPC(l + 3) << "<msub><mi>V</mi><mi>"
          << CMathMl::fixName(model->getCompartments()[i]->getObjectName())
          << "</mi></msub>" << std::endl;
          mml << SPC(l + 2) << "</mtd>" << std::endl;

          //second column ("=")
          mml << SPC(l + 2) << "<mtd>" << std::endl;
          mml << SPC(l + 3) << "<mo>=</mo>" << std::endl;
          mml << SPC(l + 2) << "</mtd>" << std::endl;

          //third column (rhs)
          mml << SPC(l + 2) << "<mtd columnalign='left'>" << std::endl;
          writeRHS_ModelEntity(mml, model->getCompartments()[i],
                               expandFull, l + 3);
          mml << SPC(l + 2) << "</mtd>" << std::endl;

          mml << SPC(l + 1) << "</mtr>" << std::endl;
        }
    }

  //write equations for metabs
  imax = model->getMetabolites().size();

  for (i = 0; i < imax; i++)
    {
      if (model->getMetabolites()[i]->getStatus() == CModelEntity::REACTIONS)
        {

          std::set<std::string> reacKeys = listReactionsForMetab(model, model->getMetabolites()[i]->getKey());
          std::set<std::string>::const_iterator it, itEnd = reacKeys.end();

          for (it = reacKeys.begin(); it != itEnd; ++it)
            {
              hasContents = true;

              mml << SPC(l + 1) << "<mtr>" << std::endl;

              //first column (lhs)
              mml << SPC(l + 2) << "<mtd>" << std::endl;

              if (it == reacKeys.begin())
                writeLHS(mml, model->getMetabolites()[i]->getObjectDisplayName(),
                         model->getMetabolites()[i]->getCompartment()->getObjectName(), l + 3);

              mml << SPC(l + 2) << "</mtd>" << std::endl;

              //second column ("=")
              mml << SPC(l + 2) << "<mtd>" << std::endl;

              if (it == reacKeys.begin())
                mml << SPC(l + 3) << "<mo>=</mo>" << std::endl;

              mml << SPC(l + 2) << "</mtd>" << std::endl;

              //third column (rhs)
              mml << SPC(l + 2) << "<mtd columnalign='left'>" << std::endl;
              writeRHS(mml, model->getMetabolites()[i],
                       dynamic_cast<CReaction*>(CCopasiRootContainer::getKeyFactory()->get(*it)) ,
                       localParameterNumbers, expand, expandFull, l + 3);
              mml << SPC(l + 2) << "</mtd>" << std::endl;

              mml << SPC(l + 1) << "</mtr>" << std::endl;
            }
        }
      else if (model->getMetabolites()[i]->getStatus() == CModelEntity::ODE)
        {
          mml << SPC(l + 1) << "<mtr>" << std::endl;

          //first column (lhs)
          mml << SPC(l + 2) << "<mtd>" << std::endl;
          writeLHS(mml, model->getMetabolites()[i]->getObjectDisplayName(),
                   model->getMetabolites()[i]->getCompartment()->getObjectName(), l + 3);
          mml << SPC(l + 2) << "</mtd>" << std::endl;

          //second column ("=")
          mml << SPC(l + 2) << "<mtd>" << std::endl;
          mml << SPC(l + 3) << "<mo>=</mo>" << std::endl;
          mml << SPC(l + 2) << "</mtd>" << std::endl;

          //third column (rhs)
          mml << SPC(l + 2) << "<mtd columnalign='left'>" << std::endl;

          std::string compName = model->getMetabolites()[i]->getCompartment()->getObjectName();
          mml << SPC(l + 3) << "<msub><mi>V</mi><mi>" << CMathMl::fixName(compName)
          << "</mi></msub>" << std::endl;
          mml << SPC(l + 3) << "<mo>" << "&CenterDot;" << "</mo>" << std::endl;

          writeRHS_ModelEntity(mml, model->getMetabolites()[i],
                               expandFull, l + 3);
          mml << SPC(l + 2) << "</mtd>" << std::endl;

          mml << SPC(l + 1) << "</mtr>" << std::endl;
        }
      else if (model->getMetabolites()[i]->getStatus() == CModelEntity::ASSIGNMENT)
        {
          mml << SPC(l + 1) << "<mtr>" << std::endl;

          //first column (lhs)
          mml << SPC(l + 2) << "<mtd>" << std::endl;
          mml << SPC(l + 3) << "<mi>" << CMathMl::fixName("[" + model->getMetabolites()[i]->getObjectName() + "]") << "</mi>" << std::endl;
          mml << SPC(l + 2) << "</mtd>" << std::endl;

          //second column ("=")
          mml << SPC(l + 2) << "<mtd>" << std::endl;
          mml << SPC(l + 3) << "<mo>=</mo>" << std::endl;
          mml << SPC(l + 2) << "</mtd>" << std::endl;

          //third column (rhs)
          mml << SPC(l + 2) << "<mtd columnalign='left'>" << std::endl;
          writeRHS_ModelEntity(mml, model->getMetabolites()[i],
                               expandFull, l + 3);
          mml << SPC(l + 2) << "</mtd>" << std::endl;

          mml << SPC(l + 1) << "</mtr>" << std::endl;
        }
    }

  //write differential equations for model values
  imax = model->getModelValues().size();

  for (i = 0; i < imax; ++i)
    if (model->getModelValues()[i]->getStatus() == CModelEntity::ODE)
      {
        hasContents = true;
        mml << SPC(l + 1) << "<mtr>" << std::endl;

        //first column (lhs)
        mml << SPC(l + 2) << "<mtd columnalign='right'>" << std::endl;
        writeLHS_ModelValue(mml, model->getModelValues()[i]->getObjectName(), l + 3);
        mml << SPC(l + 2) << "</mtd>" << std::endl;

        //second column ("=")
        mml << SPC(l + 2) << "<mtd>" << std::endl;
        mml << SPC(l + 3) << "<mo>=</mo>" << std::endl;
        mml << SPC(l + 2) << "</mtd>" << std::endl;

        //third column (rhs)
        mml << SPC(l + 2) << "<mtd columnalign='left'>" << std::endl;
        writeRHS_ModelEntity(mml, model->getModelValues()[i],
                             expandFull, l + 3);
        mml << SPC(l + 2) << "</mtd>" << std::endl;

        mml << SPC(l + 1) << "</mtr>" << std::endl;
      }

  //write assignment rules
  imax = model->getModelValues().size();

  for (i = 0; i < imax; ++i)
    if (model->getModelValues()[i]->getStatus() == CModelEntity::ASSIGNMENT)
      {
        hasContents = true;
        mml << SPC(l + 1) << "<mtr>" << std::endl;

        //first column (lhs)
        mml << SPC(l + 2) << "<mtd columnalign='right'>" << std::endl;
        mml << SPC(l + 3) << "<mi>" << CMathMl::fixName(model->getModelValues()[i]->getObjectName()) << "</mi>" << std::endl;
        //writeLHS_ModelValue(mml, model->getModelValues()[i]->getObjectName(), l + 3);
        mml << SPC(l + 2) << "</mtd>" << std::endl;

        //second column ("=")
        mml << SPC(l + 2) << "<mtd>" << std::endl;
        mml << SPC(l + 3) << "<mo>=</mo>" << std::endl;
        mml << SPC(l + 2) << "</mtd>" << std::endl;

        //third column (rhs)
        mml << SPC(l + 2) << "<mtd columnalign='left'>" << std::endl;
        writeRHS_ModelEntity(mml, model->getModelValues()[i],
                             expandFull, l + 3);
        mml << SPC(l + 2) << "</mtd>" << std::endl;

        mml << SPC(l + 1) << "</mtr>" << std::endl;
      }

  mml << SPC(l) << "</mtable>" << std::endl;
}

std::set<std::string> CMMLOutput::listReactionsForMetab(const CModel* model,
    const std::string & key)
{
  std::set<std::string> Keys;
  const CCopasiVectorN<CReaction> & Reactions = model->getReactions();
  C_INT32 j, jmax = Reactions.size();

  for (j = 0; j < jmax; j++)
    {
      const CCopasiVector <CChemEqElement> &Balances = Reactions[j]->getChemEq().getBalances();
      C_INT32 i, imax = Balances.size();

      for (i = 0; i < imax; i++)
        if (key == Balances[i]->getMetaboliteKey() && Balances[i]->getMultiplicity() != 0)
          {
            Keys.insert(Reactions[j]->getKey());
            break;
          }
    }

  return Keys;
}
