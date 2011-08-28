// Begin CVS Header
//   $Source: /fs/turing/cvs/copasi_dev/copasi/function/CEvaluationNodeCall.cpp,v $
//   $Revision: 1.34 $
//   $Name: Build-33 $
//   $Author: gauges $
//   $Date: 2010/02/19 15:15:28 $
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

#include <sbml/math/ASTNode.h>

#include "copasi.h"

#include "mathematics.h"
#include "CEvaluationNode.h"
#include "CEvaluationTree.h"
#include "CFunction.h"
#include "CExpression.h"
#include "CFunctionDB.h"
#include "CopasiDataModel/CCopasiDataModel.h"
#include "utilities/utility.h"
#include "copasi/report/CCopasiRootContainer.h"

CEvaluationNodeCall::CEvaluationNodeCall():
    CEvaluationNode(CEvaluationNode::INVALID, ""),
    mpFunction(NULL),
    mpExpression(NULL),
    mCallNodes(),
    mpCallParameters(NULL),
    mBooleanRequired(false)
{mPrecedence = PRECEDENCE_NUMBER;}

CEvaluationNodeCall::CEvaluationNodeCall(const SubType & subType,
    const Data & data):
    CEvaluationNode((Type)(CEvaluationNode::CALL | subType), data),
    mpFunction(NULL),
    mpExpression(NULL),
    mCallNodes(),
    mpCallParameters(NULL),
    mQuotesRequired(false),
    mBooleanRequired(false)
{
  std::string::size_type len = mData.length();

  if (len > 1 && mData[0] == '"' && mData[len - 1] == '"')
    {
      mQuotesRequired = true;
    }

  mData = unQuote(mData);

  switch (subType)
    {
      case FUNCTION:
      case EXPRESSION:
        break;

      default:
        fatalError();
        break;
    }

  mPrecedence = PRECEDENCE_FUNCTION;
}

CEvaluationNodeCall::CEvaluationNodeCall(const CEvaluationNodeCall & src):
    CEvaluationNode(src),
    mpFunction(src.mpFunction),
    mpExpression(src.mpExpression),
    mCallNodes(src.mCallNodes),
    mpCallParameters(NULL),
    mQuotesRequired(src.mQuotesRequired),
    mBooleanRequired(src.mBooleanRequired)
{mpCallParameters = buildParameters(mCallNodes);}

CEvaluationNodeCall::~CEvaluationNodeCall() {}

const C_FLOAT64 & CEvaluationNodeCall::value() const
{
  C_FLOAT64 &Value = *const_cast<C_FLOAT64 *>(&mValue);

  switch (mType & 0x00FFFFFF)
    {
      case FUNCTION:
      {
        std::vector< CEvaluationNode * >::const_iterator it = mCallNodes.begin();
        std::vector< CEvaluationNode * >::const_iterator end = mCallNodes.end();

        for (; it != end; ++it)(*it)->value();
      }
      return Value = mpFunction->calcValue(*mpCallParameters);
      break;

      case EXPRESSION:
        return Value = mpExpression->calcValue();
        break;
      default:
        return Value = std::numeric_limits<C_FLOAT64>::quiet_NaN();
        break;
    }
}

bool CEvaluationNodeCall::compile(const CEvaluationTree * pTree)
{
  bool success = true;
  clearParameters(mpCallParameters, mCallNodes);

  switch (mType & 0x00FFFFFF)
    {
      case FUNCTION:
        mpFunction =
          dynamic_cast<CFunction *>(CCopasiRootContainer::getFunctionList()->findFunction(mData));

        if (!mpFunction) return false;

        // We need to check whether the provided arguments match the on needed by the
        // function;
        if (!verifyParameters(mCallNodes, mpFunction->getVariables())) return false;

        mpCallParameters = buildParameters(mCallNodes);
        break;

      case EXPRESSION:
        mpExpression =
          dynamic_cast<CExpression *>(CCopasiRootContainer::getFunctionList()->findFunction(mData));

        if (!mpExpression)
          {
            // We may have a function with no arguments the parser is not able to distinguish
            // between that and an expression.
            mpFunction =
              dynamic_cast<CFunction *>(CCopasiRootContainer::getFunctionList()->findFunction(mData));

            if (!mpFunction) return false;

            mType = (CEvaluationNode::Type)(CEvaluationNode::CALL | FUNCTION);
            success = compile(pTree);
          }
        else
          {
            success = mpExpression->compile(static_cast<const CExpression *>(pTree)->getListOfContainer());
          }

        break;

      default:
        success = false;
        break;
    }

  return success;
}

bool CEvaluationNodeCall::calls(std::set< std::string > & list) const
{
  if (list.count(mData)) return true;

  CEvaluationTree * pTree =
    CCopasiRootContainer::getFunctionList()->findFunction(mData);

  if (pTree) return pTree->calls(list);

  return false;
}

std::string CEvaluationNodeCall::getInfix() const
{
  std::string Infix;

  if (mQuotesRequired)
    {
      Infix = "\"" + quote(mData, "-+^*/%(){},\t\r\n\"") + "\"(";
    }
  else
    {
      Infix = quote(mData, "-+^*/%(){},\t\r\n") + "(";
    }

  switch (mType & 0x00FFFFFF)
    {
      case FUNCTION:
      {
        std::vector< CEvaluationNode * >::const_iterator it = mCallNodes.begin();
        std::vector< CEvaluationNode * >::const_iterator end = mCallNodes.end();

        if (it != end) Infix += (*it++)->getInfix();

        for (; it != end; ++it)
          Infix += "," + (*it)->getInfix();
      }

      break;

      case EXPRESSION:
        break;

      default:
        return "@";
        break;
    }

  return Infix + ")";
}

std::string CEvaluationNodeCall::getDisplayString(const CEvaluationTree * pTree) const
{
  std::string DisplayString;

  if (mQuotesRequired)
    {
      DisplayString = "\"" + quote(mData, "-+^*/%(){},\t\r\n\"") + "\"(";
    }
  else
    {
      DisplayString = quote(mData, "-+^*/%(){},\t\r\n") + "(";
    }

  switch (mType & 0x00FFFFFF)
    {
      case FUNCTION:
      {
        std::vector< CEvaluationNode * >::const_iterator it = mCallNodes.begin();
        std::vector< CEvaluationNode * >::const_iterator end = mCallNodes.end();

        if (it != end) DisplayString += (*it++)->getDisplayString(pTree);

        for (; it != end; ++it)
          DisplayString += "," + (*it)->getDisplayString(pTree);
      }

      break;

      case EXPRESSION:
        break;

      default:
        return "@";
        break;
    }

  return DisplayString + ")";
}

std::string CEvaluationNodeCall::getDisplay_C_String(const CEvaluationTree * pTree) const
{
  std::string DisplayString;

  if (mQuotesRequired)
    {
      DisplayString = "\"" + quote(mData, "-+^*/%(){},\t\r\n\"") + "\"(";
    }
  else
    {
      DisplayString = quote(mData, "-+^*/%(){},\t\r\n") + "(";
    }

  switch (mType & 0x00FFFFFF)
    {
      case FUNCTION:
      {
        std::vector< CEvaluationNode * >::const_iterator it = mCallNodes.begin();
        std::vector< CEvaluationNode * >::const_iterator end = mCallNodes.end();

        if (it != end) DisplayString += (*it++)->getDisplay_C_String(pTree);

        for (; it != end; ++it)
          DisplayString += "," + (*it)->getDisplay_C_String(pTree);
      }

      break;

      case EXPRESSION:
        break;
      default:
        return "@";
        break;
    }

  return DisplayString + ")";
}

std::string CEvaluationNodeCall::getDisplay_MMD_String(const CEvaluationTree * /* pTree */) const
{
  std::string DisplayString;

  if (mQuotesRequired)
    {
      DisplayString = "\"" + quote(mData, "-+^*/%(){},\t\r\n\"") + "\"(";
    }
  else
    {
      DisplayString = quote(mData, "-+^*/%(){},\t\r\n") + "(";
    }

  return DisplayString;
}

std::string CEvaluationNodeCall::getDisplay_XPP_String(const CEvaluationTree * /* pTree */) const
{
  std::string DisplayString;

  if (mQuotesRequired)
    {
      DisplayString = "\"" + quote(mData, "-+^*/%(){},\t\r\n\"") + "\"(";
    }
  else
    {
      DisplayString = quote(mData, "-+^*/%(){},\t\r\n") + "(";
    }

  return DisplayString;
}

CEvaluationNode* CEvaluationNodeCall::createNodeFromASTTree(const ASTNode& node)
{
  SubType subType = CEvaluationNodeCall::FUNCTION;
  std::string data = node.getName();

  CEvaluationNodeCall* pConvertedNode = new CEvaluationNodeCall(subType, data);
  unsigned int i, iMax = node.getNumChildren();

  for (i = 0; i < iMax; ++i)
    {
      pConvertedNode->addChild(CEvaluationTree::convertASTNode(*node.getChild(i)));
    }

  return pConvertedNode;
}

ASTNode* CEvaluationNodeCall::toAST(const CCopasiDataModel* pDataModel) const
{
  ASTNode* pNode = NULL;

  pNode = new ASTNode(AST_FUNCTION);
  const std::string funName = this->getData();
  CEvaluationTree* pFun = CCopasiRootContainer::getFunctionList()->findFunction(funName);
  assert(pFun != NULL);

  if (pFun == NULL || pFun->getSBMLId().empty()) fatalError();

  pNode->setName(pFun->getSBMLId().c_str());

  const CEvaluationNode* child = static_cast<const CEvaluationNode*>(this->getChild());

  while (child)
    {
      pNode->addChild(child->toAST(pDataModel));
      child = static_cast<const CEvaluationNode*>(child->getSibling());
    }

  return pNode;
}

bool CEvaluationNodeCall::addChild(CCopasiNode< Data > * pChild,
                                   CCopasiNode< Data > * pAfter)
{
  CCopasiNode< Data >::addChild(pChild, pAfter);
  mCallNodes.push_back(static_cast<CEvaluationNode *>(pChild));

  return true;
}

bool CEvaluationNodeCall::removeChild(CCopasiNode< Data > * pChild)
{
  std::vector<CEvaluationNode *>::iterator it = mCallNodes.begin();
  std::vector<CEvaluationNode *>::iterator end = mCallNodes.end();

  while (it != end && *it != pChild) ++it;

  if (it != end) mCallNodes.erase(it);

  return CCopasiNode< Data >::removeChild(pChild);
}

CCallParameters< C_FLOAT64 > *
CEvaluationNodeCall::buildParameters(const std::vector<CEvaluationNode *> & vector)
{
  std::vector<CEvaluationNode *>::const_iterator it = vector.begin();
  std::vector<CEvaluationNode *>::const_iterator end = vector.end();

  CCallParameters< C_FLOAT64 > * pCallParameters =
    new CCallParameters< C_FLOAT64 >(vector.size());
  unsigned C_INT32 i;

  for (i = 0; it != end; ++it, i++)
    {
      if (type((*it)->getType()) == CEvaluationNode::VECTOR)
        (*pCallParameters)[i].vector = buildParameters(static_cast<const CEvaluationNodeVector *>(*it)->getVector());
      else
        (*pCallParameters)[i].value = (*it)->getValuePointer();
    }

  return pCallParameters;
}

void
CEvaluationNodeCall::clearParameters(CCallParameters< C_FLOAT64 > * pCallParameters,
                                     const std::vector<CEvaluationNode *> & vector)
{
  if (!pCallParameters) return;

  std::vector<CEvaluationNode *>::const_iterator it = vector.begin();
  std::vector<CEvaluationNode *>::const_iterator end = vector.end();

  unsigned C_INT32 i;

  for (i = 0; it != end; ++it, i++)
    {
      if (type((*it)->getType()) == CEvaluationNode::VECTOR)
        clearParameters((*pCallParameters)[i].vector,
                        static_cast<const CEvaluationNodeVector *>(*it)->getVector());
    }

  delete pCallParameters;
  return;
}

bool
CEvaluationNodeCall::verifyParameters(const std::vector<CEvaluationNode *> & vector,
                                      const CFunctionParameters & functionParameters)
{
  if (vector.size() != functionParameters.size()) return false;

  std::vector<CEvaluationNode *>::const_iterator it = vector.begin();
  std::vector<CEvaluationNode *>::const_iterator end = vector.end();

  unsigned C_INT32 i;

  for (i = 0; it != end; ++it, i++)
    {
      if ((type((*it)->getType()) == CEvaluationNode::VECTOR &&
           functionParameters[i]->getType() != CFunctionParameter::VFLOAT64) ||
          functionParameters[i]->getType() == CFunctionParameter::VFLOAT64)
        return false;
    }

  return true;
}

const CEvaluationTree * CEvaluationNodeCall::getCalledTree() const
{
  switch (mType & 0x00FFFFFF)
    {
      case FUNCTION:
      case EXPRESSION:
        return CCopasiRootContainer::getFunctionList()->findFunction(mData);

      default:
        return NULL;
    }
}

#include "utilities/copasimathml.h"

void CEvaluationNodeCall::writeMathML(std::ostream & out,
                                      const std::vector<std::vector<std::string> > & env,
                                      bool expand,
                                      unsigned C_INT32 l) const
{
  switch (mType & 0x00FFFFFF)
    {
      case FUNCTION:
      {

        if (!expand || !mpFunction)
          {
            out << SPC(l) << "<mrow>" << std::endl;

            out << SPC(l + 1) << "<mi>" << mData << "</mi>" << std::endl;
            out << SPC(l + 1) << "<mo> &ApplyFunction; </mo>" << std::endl;
            out << SPC(l + 1) << "<mrow>" << std::endl;
            out << SPC(l + 2) << "<mo> (</mo>" << std::endl;
            out << SPC(l + 2) << "<mrow>" << std::endl;

            std::vector< CEvaluationNode * >::const_iterator it = mCallNodes.begin();
            std::vector< CEvaluationNode * >::const_iterator end = mCallNodes.end();

            if (it != end)(*it++)->writeMathML(out, env, expand, l + 3);

            for (; it != end; ++it)
              {

                out << SPC(l + 3) << "<mo> , </mo>" << std::endl;
                (*it)->writeMathML(out, env, expand, l + 3);
              }

            out << SPC(l + 2) << "</mrow>" << std::endl;
            out << SPC(l + 2) << "<mo>) </mo>" << std::endl;

            out << SPC(l + 1) << "</mrow>" << std::endl;
            out << SPC(l) << "</mrow>" << std::endl;
          }
        else
          {
            //construct the environment for the nested function
            std::vector<std::vector<std::string> > env2;

            std::vector< CEvaluationNode * >::const_iterator it = mCallNodes.begin();
            std::vector< CEvaluationNode * >::const_iterator end = mCallNodes.end();

            for (; it != end; ++it)
              {
                std::ostringstream oss;
                (*it)->writeMathML(oss, env, expand, l + 3);
                std::vector<std::string> tmpvector; tmpvector.push_back(oss.str());
                env2.push_back(tmpvector);
              }

            out << SPC(l) << "<mfenced>" << std::endl;
            mpFunction->writeMathML(out, env2, expand, expand, l + 1);
            out << SPC(l) << "</mfenced>" << std::endl;
          }
      }
      break;

      case EXPRESSION:
        break;
      default:
        break;
    }

  return;
}

void CEvaluationNodeCall::setBooleanRequired(const bool & booleanRequired)
{mBooleanRequired = booleanRequired;}

const bool & CEvaluationNodeCall::isBooleanRequired() const
{return mBooleanRequired;}

// virtual
bool CEvaluationNodeCall::isBoolean() const
{
  const CEvaluationTree * pEvaluationTree = getCalledTree();

  if (pEvaluationTree != NULL)
    {
      return pEvaluationTree->isBoolean();
    }

  return false;
}
