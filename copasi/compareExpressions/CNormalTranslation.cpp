// Begin CVS Header
//   $Source: /fs/turing/cvs/copasi_dev/copasi/compareExpressions/CNormalTranslation.cpp,v $
//   $Revision: 1.45 $
//   $Name: Build-33 $
//   $Author: shoops $
//   $Date: 2009/10/27 16:50:08 $
// End CVS Header

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

#include <sstream>
#include <string>
#include <set>
#include <vector>
#include <iostream>

#include "copasi.h"

#include "CNormalBase.h"
#include "ConvertToCEvaluationNode.h"
#include "CNormalTranslation.h"
#include "CNormalFraction.h"
#include "CNormalSum.h"
#include "CNormalProduct.h"
#include "CNormalLogical.h"

#include "function/CEvaluationTree.h"
#include "function/CEvaluationNodeConstant.h"

const double CNormalTranslation::ZERO = 1e-100;
const unsigned int CNormalTranslation::RECURSION_LIMIT = 20;
const CEvaluationNode CNormalTranslation::NEUTRAL_ELEMENT_ADD = CEvaluationNodeNumber(CEvaluationNodeNumber::DOUBLE, "0.0");
const CEvaluationNode CNormalTranslation::NEUTRAL_ELEMENT_MULTIPLY = CEvaluationNodeNumber(CEvaluationNodeNumber::DOUBLE, "1.0");
const CEvaluationNode CNormalTranslation::NEUTRAL_ELEMENT_OR = CEvaluationNodeConstant(CEvaluationNodeConstant::FALSE, "FALSE");
const CEvaluationNode CNormalTranslation::NEUTRAL_ELEMENT_AND = CEvaluationNodeConstant(CEvaluationNodeConstant::TRUE, "TRUE");
const CEvaluationNode CNormalTranslation::ZERO_NODE = CNormalTranslation::NEUTRAL_ELEMENT_ADD;
const CEvaluationNode CNormalTranslation::ONE_NODE = CNormalTranslation::NEUTRAL_ELEMENT_MULTIPLY;
const CEvaluationNode CNormalTranslation::PLUS_NODE = CEvaluationNodeOperator(CEvaluationNodeOperator::PLUS, "+");
const CEvaluationNode CNormalTranslation::TIMES_NODE = CEvaluationNodeOperator(CEvaluationNodeOperator::MULTIPLY, "*");

recursion_limit_exception::recursion_limit_exception(recursion_limit_exception::LIMIT_TYPE type): std::exception(), mType(type)
{}

/**
 * Simplify an evaluation tree given by the root node by creating a new simplified tree from the original one.
 * The tree itself is actually not created!
 * @return CEvaluationNode*, pointer to root node of the newly created tree.
 */
CEvaluationNode* CNormalTranslation::simplifyTree(const CEvaluationNode* node)
{
  const CEvaluationNode * child = dynamic_cast<const CEvaluationNode*>(node->getChild());
  CEvaluationNode * newchild = NULL;
  std::vector<CEvaluationNode*> children;

  while (child != NULL)
    {
      newchild = simplifyTree(child);
      child = dynamic_cast<const CEvaluationNode*>(child->getSibling());
      children.push_back(newchild);
    }

  CEvaluationNode* newnode = node->simplifyNode(children);
  return newnode;
}

/**
 * Creating a simplified tree by calling simplifyTree repeatedly until it cannot be simplified further.
 * The tree itself is actually not created!
 * @return CEvaluationNode*, pointer to root node of the newly created tree.
 */
CEvaluationNode * CNormalTranslation::simplifyTreeReptdly(const CEvaluationNode* root0)
{
  CEvaluationNode * root1 = simplifyTree(root0);

  if (root1->getInfix() != root0->getInfix())
    {
      CEvaluationNode * root2 = simplifyTreeReptdly(root1);
      delete root1;
      return root2;
    }
  else
    {
      return root1;
    }
}

/**
 * Translate and simplify a tree given by the root node into CNormal structure
 * @return CNormalFraction*
 */
CNormalFraction* CNormalTranslation::normAndSimplify(const CEvaluationNode* root0)
{
  //CEvaluationNode * root1 = simplifyTreeReptdly(root0);
  CEvaluationNode* root1 = CNormalTranslation::simplify(root0);
  CEvaluationNode* root2 = CNormalTranslation::expandPowerExponents(root1);
  delete root1;
  CNormalFraction* base = createNormalRepresentation(root2);
  base->simplify();

  delete root2;

  return base;
}

/**
 * Translate and simplify a tree by calling normAndSimplify repeatedly until it cannot be simplified further
 * @return CNormalFraction*
 */
CNormalFraction* CNormalTranslation::normAndSimplifyReptdly(const CEvaluationTree* tree0, unsigned int depth)
{
  if (depth > RECURSION_LIMIT) throw recursion_limit_exception(recursion_limit_exception::NORM_AND_SIMPLIFY);

  const CEvaluationNode* root0 = tree0->getRoot();

  CNormalFraction * base0 = normAndSimplify(root0);

  std::stringstream tmp;
  tmp << base0->toString();

  CEvaluationTree * tree1 = new CEvaluationTree("second tree", NULL, CEvaluationTree::Function);

  tree1->setInfix(tmp.str());

  if (tree1->getInfix() != tree0->getInfix())
    {
      CNormalFraction * base1 = normAndSimplifyReptdly(tree1, ++depth);
      delete tree1;
      delete base0;
      return base1;
    }
  else
    {
      delete tree1;
      return base0;
    }
}

/**
 * Translate and simplify a tree by calling normAndSimplify repeatedly until it cannot be simplified further
 * @return CNormalFraction*
 */
CNormalFraction* CNormalTranslation::normAndSimplifyReptdly(const CEvaluationNode* root0, unsigned int depth)
{
  if (depth > RECURSION_LIMIT) throw recursion_limit_exception(recursion_limit_exception::NORM_AND_SIMPLIFY);

  CNormalFraction * base0 = normAndSimplify(root0);

  std::stringstream tmp;

  //CEvaluationTree * tree1 = new CEvaluationTree("second tree", NULL, CEvaluationTree::Function);
  CEvaluationNode* pTmpNode = convertToCEvaluationNode(*base0);
  assert(pTmpNode != NULL);
  CNormalFraction* pFraction = dynamic_cast<CNormalFraction*>(base0);
  assert(pFraction != NULL);

  if (pTmpNode->getInfix() != root0->getInfix())
    {
      CNormalFraction * base1 = normAndSimplifyReptdly(pTmpNode, ++depth);
      delete pTmpNode;
      delete base0;
      return base1;
    }
  else
    {
      delete pTmpNode;
      return base0;
    }
}

CEvaluationNode* CNormalTranslation::expandPowerExponents(const CEvaluationNode* pRoot)
{
  CEvaluationNode* pResult = NULL;
  const CEvaluationNode* pChild = dynamic_cast<const CEvaluationNode*>(pRoot->getChild());
  std::vector<CEvaluationNode*> children;

  while (pChild != NULL)
    {
      CEvaluationNode* pNewChild = CNormalTranslation::expandPowerExponents(pChild);
      children.push_back(pNewChild);
      pChild = dynamic_cast<const CEvaluationNode*>(pChild->getSibling());
    }

  if (CEvaluationNode::type(pRoot->getType()) == CEvaluationNode::OPERATOR && ((CEvaluationNodeOperator::SubType)CEvaluationNode::subType(pRoot->getType())) == CEvaluationNodeOperator::POWER)
    {
      assert(children.size() == 2);
      std::vector<const CEvaluationNode*> summands;
      CNormalTranslation::findSummands(children[1], summands);
      // for each summand create a power node with a copy of the first child
      // in children as child 1 and a copy of the summand as child 2
      std::vector<CEvaluationNode*> numeratorNodes;
      std::vector<CEvaluationNode*> denominatorNodes;
      std::vector<const CEvaluationNode*>::iterator it = summands.begin(), endit = summands.end();

      while (it != endit)
        {
          CEvaluationNodeOperator* pPowerNode = new CEvaluationNodeOperator(CEvaluationNodeOperator::POWER, "^");
          pPowerNode->addChild(children[0]->copyBranch());

          if (CEvaluationNode::type((*it)->getType()) == CEvaluationNode::FUNCTION && ((CEvaluationNodeFunction::SubType)CEvaluationNode::subType((*it)->getType())) == CEvaluationNodeFunction::MINUS)
            {
              pPowerNode->addChild(dynamic_cast<const CEvaluationNode*>((*it)->getChild())->copyBranch());
              denominatorNodes.push_back(pPowerNode);
            }
          else if ((CEvaluationNode::type((*it)->getType()) == CEvaluationNode::NUMBER && dynamic_cast<const CEvaluationNodeNumber*>((*it))->value() < 0.0))
            {
              std::ostringstream os;
              os.precision(18);
              os << fabs(dynamic_cast<const CEvaluationNodeNumber*>(*it)->value());
              pPowerNode->addChild(new CEvaluationNodeNumber((CEvaluationNodeNumber::SubType)CEvaluationNode::subType((*it)->getType()), os.str().c_str()));
              denominatorNodes.push_back(pPowerNode);
            }
          else
            {
              pPowerNode->addChild((*it)->copyBranch());
              numeratorNodes.push_back(pPowerNode);
            }

          ++it;
        }

      delete children[0];
      delete children[1];

      // create the numerator chain
      if (numeratorNodes.empty())
        {
          pResult = new CEvaluationNodeNumber(CEvaluationNodeNumber::DOUBLE, "1.0");
        }
      else
        {
          pResult = CNormalTranslation::createChain(&CNormalTranslation::TIMES_NODE, &CNormalTranslation::ONE_NODE, numeratorNodes);
        }

      assert(pResult != NULL);

      // if there are items in the denominator vector create the denominator
      // chain and divide the numerato chain by the denominator chain
      if (!denominatorNodes.empty())
        {
          CEvaluationNodeOperator* pDivision = new CEvaluationNodeOperator(CEvaluationNodeOperator::DIVIDE, "/");
          pDivision->addChild(pResult);
          pDivision->addChild(CNormalTranslation::createChain(&CNormalTranslation::TIMES_NODE, &CNormalTranslation::ONE_NODE, denominatorNodes));
          pResult = pDivision;
        }
    }
  else
    {
      // copy the node and add the children
      pResult = CEvaluationNode::create(pRoot->getType(), pRoot->getData());
      std::vector<CEvaluationNode*>::iterator it = children.begin(), endit = children.end();

      while (it != endit)
        {
          pResult->addChild(*it);
          ++it;
        }
    }

  return pResult;
}

CEvaluationNode* CNormalTranslation::createOperatorChain(CEvaluationNodeOperator::SubType type, const char* data, const std::vector<const CEvaluationNode*>& nodes)
{
  CEvaluationNode* pResult = NULL;

  if (nodes.size() == 0)
    {
      pResult = new CEvaluationNodeNumber(CEvaluationNodeNumber::DOUBLE, "1.0");
    }
  else if (nodes.size() == 1)
    {
      pResult = nodes[0]->copyBranch();
    }
  else
    {
      // start from the back to create the deepest nodes first
      std::vector<const CEvaluationNode*>::const_reverse_iterator it = nodes.rbegin(), endit = nodes.rend();
      CEvaluationNode* pOperator = new CEvaluationNodeOperator(type, data);
      CEvaluationNode* pChild2 = (*it)->copyBranch();
      ++it;
      CEvaluationNode* pChild1 = (*it)->copyBranch();
      pOperator->addChild(pChild1);
      pOperator->addChild(pChild2);
      ++it;
      pChild2 = pOperator;

      while (it != endit)
        {
          pOperator = new CEvaluationNodeOperator(type, data);
          pOperator->addChild((*it)->copyBranch());
          pOperator->addChild(pChild2);
          pChild2 = pOperator;
          ++it;
        }

      pResult = pOperator;
    }

  return pResult;
}

void CNormalTranslation::findSummands(const CEvaluationNode* pRoot, std::vector<const CEvaluationNode*>& summands)
{
  if (CEvaluationNode::type(pRoot->getType()) == CEvaluationNode::OPERATOR && ((CEvaluationNodeOperator::SubType)CEvaluationNode::subType(pRoot->getType())) == CEvaluationNodeOperator::PLUS)
    {
      const CEvaluationNode* pChild1 = dynamic_cast<const CEvaluationNode*>(pRoot->getChild());
      assert(pChild1 != NULL);

      if (pChild1 != NULL)
        {
          const CEvaluationNode* pChild2 = dynamic_cast<const CEvaluationNode*>(pChild1->getSibling());
          assert(pChild2 != NULL);

          if (pChild2 != NULL)
            {
              assert(pChild2->getSibling() == NULL);

              if (CEvaluationNode::type(pChild1->getType()) == CEvaluationNode::OPERATOR && ((CEvaluationNodeOperator::SubType)CEvaluationNode::subType(pChild1->getType())) == CEvaluationNodeOperator::PLUS)
                {
                  CNormalTranslation::findSummands(pChild1, summands);
                }
              else
                {
                  summands.push_back(pChild1);
                }

              if (CEvaluationNode::type(pChild2->getType()) == CEvaluationNode::OPERATOR && ((CEvaluationNodeOperator::SubType)CEvaluationNode::subType(pChild2->getType())) == CEvaluationNodeOperator::PLUS)
                {
                  CNormalTranslation::findSummands(pChild2, summands);
                }
              else
                {
                  summands.push_back(pChild2);
                }
            }
        }
    }
  else
    {
      summands.push_back(pRoot);
    }
}

// new routines

/**
 * This method elminates subexpressions from an expression
 */
CEvaluationNode* CNormalTranslation::eliminate(const CEvaluationNode* pOrig)
{
  CEvaluationNode* pResult = pOrig->copyBranch();
  CEvaluationNode* pTmp = NULL;
  std::string infix = pResult->getInfix(); //base->toString();
  //delete base;
  bool changed = true;

  while (changed)
    {
      // first make elementary eliminations
      pTmp = CNormalTranslation::elementaryElimination(pResult);

      if (pTmp != pResult) delete pResult;

      // now get rid of nested powers a^b^c
      pResult = CNormalTranslation::eliminateNestedPowers(pTmp);
      delete pTmp;
      pTmp = pResult;
      // eliminate fractions within powers
      // (a/b)^3 -> a^3 / b^3
      // now get rid of directly nested fractions
      pResult = CNormalTranslation::eliminatePowersOfFractions(pTmp);
      delete pResult;
      pResult = CNormalTranslation::eliminateDirectlyNestedFractions(pTmp);
      delete pTmp;
      // now cancel since cancelation can lead to new nodes for which
      // elementary elimination would be possible, we might have to run
      // this loop again
      pTmp = CNormalTranslation::cancel(pResult);
      delete pResult;

      // check if we are done
      // we are done if the infix has not changed over one loop run
      if (/*base->toString()*/pTmp->getInfix() == infix)
        {
          changed = false;
        }
      else
        {
          infix = pTmp->getInfix(); //base->toString();
        }

      pResult = pTmp;
      //delete base;
      //base = NULL;
    }

  return pResult;
}

/**
 * This routine is responsible for recursively simplifying a given
 * CEvaluationNode based tree.
 */
CEvaluationNode* CNormalTranslation::simplify(const CEvaluationNode* pOrig)
{
  CEvaluationNode* pResult = NULL;
  bool finished = false;
  //CNormalFraction* base = createNormalRepresentation(pOrig);
  //assert(base != NULL);
  std::string infix = pOrig->getInfix(); //base->toString();
  std::string infix2 = infix;
  //delete base;
  //base = NULL;
  CEvaluationNode* pTmp = pOrig->copyBranch();
  unsigned int counter = 0;

  while (!finished)
    {
      ++counter;

      if (counter > RECURSION_LIMIT) throw recursion_limit_exception(recursion_limit_exception::SIMPLIFY);

      pResult = CNormalTranslation::eliminate(pTmp);
      delete pTmp;
      // now we evaluate everything that can be evaluated, e.g. operations on
      // numbers
      pTmp = CNormalTranslation::evaluateNumbers(pResult);

      if (pTmp != pResult)
        {
          delete pResult;
        }

      // this method combines identical multiplicants and summands
      pResult = CNormalTranslation::cancel(pTmp);
      delete pTmp;
      // now expand products in bases to power operators
      pTmp = CNormalTranslation::expandPowerBases(pResult);
      delete pResult;
      // now expand the exponents in the power nodes and multiply products
      // expansions can lead to new cancelations being possible so we might
      // need to rerun the whole loop
      pResult = CNormalTranslation::expandPowerNodes(pTmp);
      delete pTmp;
      pTmp = CNormalTranslation::expandProducts(pResult);

      if (pTmp != pResult)
        {
          delete pResult;
        }

      // check if we are done
      // we are done, once the infix has not changed during one loop run
      //base = createNormalRepresentation(pResult);
      //assert(base != NULL);
      pResult = pTmp;

      if (/*base->toString()*/pResult->getInfix() == infix)
        {
          finished = true;
        }
      else
        {
          infix = pResult->getInfix(); //base->toString();
        }

      pTmp = pResult;
    }

  pTmp = CNormalTranslation::product2fraction(pResult);
  delete pResult;
  pResult = pTmp;
  return pResult;
}

/**
 * This routine is responsible for all elementary eliminations, e.g. addition
 * of 0.
 * These steps can not lead to new simplifications in the children of the node
 * being simplified, so it is not necessary to run this on the children again.
 */
CEvaluationNode* CNormalTranslation::elementaryElimination(CEvaluationNode* pOrig)
{
  // this is done depth first
  CEvaluationNode* pResult = pOrig;
  CEvaluationNode* pChild = dynamic_cast<CEvaluationNode*>(pOrig->getChild());
  CEvaluationNode* pLastChild = pOrig;

  while (pChild != NULL)
    {
      CEvaluationNode* pNewChild = elementaryElimination(pChild);
      assert(pNewChild != NULL);

      if (pNewChild != pChild)
        {
          // remove the old child and add the new one
          pOrig->removeChild(pChild);
          delete pChild;
          pChild = pNewChild;
          pOrig->addChild(pNewChild, pLastChild);
        }

      pLastChild = pChild;
      pChild = dynamic_cast<CEvaluationNode*>(pChild->getSibling());
    }

  if (CEvaluationNode::type(pOrig->getType()) == CEvaluationNode::OPERATOR)
    {
      // check if we can eliminate anything
      // check if one of the children is (-)0, (-)1, NaN or INFINITY
      switch ((CEvaluationNodeOperator::SubType)CEvaluationNode::subType(pOrig->getType()))
        {
          case CEvaluationNodeOperator::POWER:
            pResult = CNormalTranslation::elementaryEliminationPower(pOrig);
            break;
          case CEvaluationNodeOperator::MODULUS:
            pResult = CNormalTranslation::elementaryEliminationModulus(pOrig);
            break;
          case CEvaluationNodeOperator::MULTIPLY:
            pResult = CNormalTranslation::elementaryEliminationMultiply(pOrig);
            break;
          case CEvaluationNodeOperator::DIVIDE:
            pResult = CNormalTranslation::elementaryEliminationDivide(pOrig);
            break;
          case CEvaluationNodeOperator::PLUS:
            pResult = CNormalTranslation::elementaryEliminationPlus(pOrig);
            break;
          case CEvaluationNodeOperator::MINUS:
            pResult = CNormalTranslation::elementaryEliminationMinus(pOrig);
            break;
          default:
            // we should never end up here
            fatalError();
            break;
        }
    }
  else if (CEvaluationNode::type(pOrig->getType()) == CEvaluationNode::FUNCTION)
    {
      pResult = CNormalTranslation::elementaryEliminationFunction(pOrig);
    }

  return pResult;
}

/**
 * This method makes elementary eliminations on function nodes
 */
CEvaluationNode* CNormalTranslation::elementaryEliminationFunction(CEvaluationNode* pFunctionNode)
{
  // PLUS(X) -> X
  // X(NaN) -> NaN
  // MINUX(X) where X is a number -> -X
  CEvaluationNode* pResult = pFunctionNode;
  CEvaluationNode* pChild = NULL;

  switch ((CEvaluationNodeFunction::SubType)CEvaluationNode::subType(pFunctionNode->getType()))
    {
      case CEvaluationNodeFunction::INVALID:
        break;
      case CEvaluationNodeFunction::PLUS:
        pChild = dynamic_cast<CEvaluationNode*>(pFunctionNode->getChild());
        assert(pChild != NULL);
        assert(pChild->getSibling() == NULL);
        pResult = pChild->copyBranch();
        break;
      case CEvaluationNodeFunction::MINUS:
        pChild = dynamic_cast<CEvaluationNode*>(pFunctionNode->getChild());
        assert(pChild != NULL);
        assert(pChild->getSibling() == NULL);

        if (CEvaluationNode::type(pChild->getType()) == CEvaluationNode::NUMBER)
          {
            std::ostringstream os;
            os.precision(18);
            os << -1.0*dynamic_cast<const CEvaluationNodeNumber*>(pChild)->value();
            pResult = new CEvaluationNodeNumber(CEvaluationNodeNumber::DOUBLE, os.str().c_str());
          }
        else if (CEvaluationNode::type(pChild->getType()) == CEvaluationNode::CONSTANT &&
                 ((CEvaluationNodeConstant::SubType)CEvaluationNode::subType(pChild->getType())) == CEvaluationNodeConstant::_NaN)
          {
            pResult = new CEvaluationNodeConstant(CEvaluationNodeConstant::_NaN, "NaN");
          }

        if (pResult == pFunctionNode)
          {
            // MINUS(X) -> -1.0 * X
            pResult = new CEvaluationNodeOperator(CEvaluationNodeOperator::MULTIPLY, "*");
            pResult->addChild(new CEvaluationNodeNumber(CEvaluationNodeNumber::DOUBLE, "-1.0"));
            pResult->addChild(pChild->copyBranch());
          }

        break;
      default:
        pChild = dynamic_cast<CEvaluationNode*>(pFunctionNode->getChild());

        while (pChild != NULL)
          {
            if (CEvaluationNode::type(pChild->getType()) == CEvaluationNode::CONSTANT && ((CEvaluationNodeConstant::SubType)CEvaluationNode::subType(pChild->getType())) == CEvaluationNodeConstant::_NaN)
              {
                pResult = new CEvaluationNodeConstant(CEvaluationNodeConstant::_NaN, "NaN");
                break;
              }

            pChild = dynamic_cast<CEvaluationNode*>(pChild->getSibling());
          }

        break;
    }

  return pResult;
}

/**
 * This method makes the elementary elimination on a power node.
 */
CEvaluationNode* CNormalTranslation::elementaryEliminationPower(CEvaluationNode* pPowerNode)
{
  CEvaluationNode* pResult = pPowerNode;
  assert(CEvaluationNode::type(pPowerNode->getType()) == CEvaluationNode::OPERATOR);
  assert(((CEvaluationNodeOperator::SubType)CEvaluationNode::subType(pPowerNode->getType())) == CEvaluationNodeOperator::POWER);
  CEvaluationNode* pChild1 = dynamic_cast<CEvaluationNode*>(pPowerNode->getChild());
  assert(pChild1 != NULL);
  CEvaluationNode* pChild2 = dynamic_cast<CEvaluationNode*>(pChild1->getSibling());
  assert(pChild2 != NULL);
  assert(pChild2->getSibling() == NULL);

  if (CEvaluationNode::type(pChild1->getType()) == CEvaluationNode::NUMBER)
    {
      // 0 and 1
      CEvaluationNodeNumber* pNumberNode = dynamic_cast<CEvaluationNodeNumber*>(pChild1);
      assert(pNumberNode != NULL);
      double value = pNumberNode->value();

      if (fabs(value) < ZERO)
        {
          // 0^(NaN) -> NaN
          if (pChild2->getType() == CEvaluationNode::CONSTANT && ((CEvaluationNodeConstant::SubType)CEvaluationNode::subType(pChild2->getType())) == CEvaluationNodeConstant::_NaN)
            {
              pResult = new CEvaluationNodeConstant(CEvaluationNodeConstant::_NaN, "NaN");
            }
          else if (pChild2->getType() == CEvaluationNode::NUMBER)
            {
              CEvaluationNodeNumber* pNumberNode2 = dynamic_cast<CEvaluationNodeNumber*>(pChild2);
              double value = pNumberNode2->value();

              // 0^0 -> NaN
              // 0^(-x) -> NaN
              if (fabs(value) < ZERO || value < 0.0)
                {
                  pResult = new CEvaluationNodeConstant(CEvaluationNodeConstant::_NaN, "NaN");
                }
            }

          // 0^x -> 0
          if (pResult == pPowerNode)
            {
              pResult = new CEvaluationNodeNumber(CEvaluationNodeNumber::DOUBLE, "0.0");
            }
        }
      else if (fabs(value - 1.0) < ZERO)
        {
          // 1^NaN -> NaN
          // 1^x -> 1
          if (CEvaluationNode::type(pChild2->getType()) == CEvaluationNode::CONSTANT && ((CEvaluationNodeConstant::SubType)CEvaluationNode::subType(pChild2->getType())) == CEvaluationNodeConstant::_NaN)
            {
              pResult = new CEvaluationNodeConstant(CEvaluationNodeConstant::_NaN, "NaN");
            }

          if (pResult == NULL)
            {
              pResult = new CEvaluationNodeNumber(CEvaluationNodeNumber::DOUBLE, "1.0");
            }
        }

      /* ignore -1 for now
         else if(fabs(value + 1.0) < ZERO)
         {
         }
         */
    }
  else if (CEvaluationNode::type(pChild1->getType()) == CEvaluationNode::CONSTANT)
    {
      // infinity and NaN
      if (((CEvaluationNodeConstant::SubType)CEvaluationNode::subType(pChild1->getType())) == CEvaluationNodeConstant::_NaN)
        {
          // NaN^x -> NaN
          pResult = new CEvaluationNodeConstant(CEvaluationNodeConstant::_NaN, "NaN");
        }
      else if (((CEvaluationNodeConstant::SubType)CEvaluationNode::subType(pChild1->getType())) == CEvaluationNodeConstant::_INFINITY)
        {
          // INFINITY^(-NaN) -> NaN
          // INFINITY^-x -> 0.0 // x being a positive number
          // INFINITY^x -> INFINITY // x being a positive number
          // INFINITY^0 -> 1
          if (CEvaluationNode::type(pChild2->getType()) == CEvaluationNode::NUMBER)
            {
              CEvaluationNodeNumber* pNumberNode2 = dynamic_cast<CEvaluationNodeNumber*>(pChild2);
              assert(pNumberNode2 != NULL);
              double value = pNumberNode2->value();

              // INFINITY^0 -> 1
              if (fabs(value) < ZERO)
                {
                  pResult = new CEvaluationNodeNumber(CEvaluationNodeNumber::DOUBLE, "1.0");
                }
              // INFINITY^x -> INFINITY // x being a positive number
              else if (value > 0.0)
                {
                  pResult = new CEvaluationNodeConstant(CEvaluationNodeConstant::_INFINITY, "inf");
                }
              // INFINITY^-x -> 0.0 // x being a positive number
              else
                {
                  pResult = new CEvaluationNodeNumber(CEvaluationNodeNumber::DOUBLE, "0.0");
                }
            }
          else if (CEvaluationNode::type(pChild2->getType()) == CEvaluationNode::CONSTANT && ((CEvaluationNodeConstant::SubType)CEvaluationNode::subType(pChild2->getType())) == CEvaluationNodeConstant::_NaN)
            {
              // INFINITY^NaN    -> NaN
              pResult = new CEvaluationNodeConstant(CEvaluationNodeConstant::_NaN, "NaN");
            }
          /* the minus function is eliminated
          else if(CEvaluationNode::type(pChild2->getType())==CEvaluationNode::FUNCTION && ((CEvaluationNodeFunction::SubType)CEvaluationNode::subType(pChild2->getType()))==CEvaluationNodeFunction::MINUS)
          {
              CEvaluationNode* pChild=dynamic_cast<CEvaluationNode*>(pChild2->getChild());
              // INFINITY^(-CONSTANT) -> 0.0 // where CONSTANT != NaN
              if(CEvaluationNode::type(pChild->getType())==CEvaluationNode::CONSTANT)
              {
                  pResult=new CEvaluationNodeNumber(CEvaluationNodeNumber::DOUBLE,"0.0");
              }
          }
          */
          else if (CEvaluationNode::type(pChild2->getType()) == CEvaluationNode::NUMBER)
            {
              CEvaluationNodeNumber* pNumberNode2 = dynamic_cast<CEvaluationNodeNumber*>(pChild2);
              assert(pNumberNode2 != NULL);
              double value = pNumberNode2->value();

              // INFINITY^0 -> 1.0
              if (fabs(value) < ZERO)
                {
                  pResult = new CEvaluationNodeNumber(CEvaluationNodeNumber::DOUBLE, "1.0");
                }
              // INFINITY^(-x) -> 0.0
              else if (value > 0.0)
                {
                  pResult = new CEvaluationNodeNumber(CEvaluationNodeNumber::DOUBLE, "0.0");
                }
              // INFINITY^x -> INFINITY
              else
                {
                  pResult = new CEvaluationNodeConstant(CEvaluationNodeConstant::_INFINITY, "inf");
                }
            }

          // INFINITY ^ x -> INFINITY
          if (pResult == NULL)
            {
              pResult = new CEvaluationNodeConstant(CEvaluationNodeConstant::_INFINITY, "inf");
            }
        }
    }
  else if (CEvaluationNode::type(pChild2->getType()) == CEvaluationNode::NUMBER)
    {
      // 0 and 1
      CEvaluationNodeNumber* pNumberNode2 = dynamic_cast<CEvaluationNodeNumber*>(pChild2);
      assert(pNumberNode2 != NULL);
      double value = pNumberNode2->value();

      // x^0 -> 1.0
      if (fabs(value) < ZERO)
        {
          pResult = new CEvaluationNodeNumber(CEvaluationNodeNumber::DOUBLE, "1.0");
        }
      else if (fabs(value - 1.0) < ZERO)
        {
          // make a deep copy of the first child
          pResult = pChild1->copyBranch();
        }

      /* ignore -1 since this may interfere with other simplification
       * mechanisms.
       * Negative exponents will be eliminated in the end.
       else if(fabs(value + 1.0) < ZERO)
       {
       }
       */
    }
  else if (CEvaluationNode::type(pChild2->getType()) == CEvaluationNode::CONSTANT)
    {
      // infinity and NaN
      if (((CEvaluationNodeConstant::SubType)CEvaluationNode::subType(pChild2->getType())) == CEvaluationNodeConstant::_NaN)
        {
          pResult = pChild2->copyBranch();
        }
      else if (((CEvaluationNodeConstant::SubType)CEvaluationNode::subType(pChild2->getType())) == CEvaluationNodeConstant::_INFINITY)
        {
          pResult = pChild2->copyBranch();
        }
    }

  return pResult;
}

/**
 * This method makes the elementary elimination on a modulus node.
 */
CEvaluationNode* CNormalTranslation::elementaryEliminationModulus(CEvaluationNode* pModulusNode)
{
  CEvaluationNode* pResult = pModulusNode;
  assert(CEvaluationNode::type(pModulusNode->getType()) == CEvaluationNode::OPERATOR);
  assert(((CEvaluationNodeOperator::SubType)CEvaluationNode::subType(pModulusNode->getType())) == CEvaluationNodeOperator::MODULUS);
  CEvaluationNode* pChild1 = dynamic_cast<CEvaluationNode*>(pModulusNode->getChild());
  assert(pChild1 != NULL);
  CEvaluationNode* pChild2 = dynamic_cast<CEvaluationNode*>(pChild1->getSibling());
  assert(pChild2 != NULL);
  assert(pChild2->getSibling() == NULL);

  // if one child is NaN, the result is NaN
  if ((CEvaluationNode::type(pChild1->getType()) == CEvaluationNode::CONSTANT &&
       ((CEvaluationNodeConstant::SubType)(CEvaluationNode::subType(pChild1->getType()))) == CEvaluationNodeConstant::_NaN) ||
      (CEvaluationNode::type(pChild2->getType()) == CEvaluationNode::CONSTANT &&
       ((CEvaluationNodeConstant::SubType)(CEvaluationNode::subType(pChild2->getType()))) == CEvaluationNodeConstant::_NaN))
    {
      pResult = new CEvaluationNodeConstant(CEvaluationNodeConstant::_NaN, "NAN");
    }

  // X%X -> 0
  CNormalFraction* base1 = createNormalRepresentation(pChild1);
  CNormalFraction* base2 = createNormalRepresentation(pChild2);

  if (base1->toString() == base2->toString())
    {
      pResult = new CEvaluationNodeNumber(CEvaluationNodeNumber::DOUBLE, "0.0");
    }
  else if (CEvaluationNode::type(pChild1->getType()) == CEvaluationNode::NUMBER)
    {
      // 0 and 1
      CEvaluationNodeNumber* pNumberNode = dynamic_cast<CEvaluationNodeNumber*>(pChild1);
      assert(pNumberNode != NULL);
      double value = pNumberNode->value();

      if (fabs(value) < ZERO)
        {
          pResult = new CEvaluationNodeNumber(CEvaluationNodeNumber::DOUBLE, "0.0");
        }
      else if (fabs(value - 1.0) < ZERO)
        {
          // 1%X where X is any number other than 1 will give 1.0
          // the case where X is 1 is already covered above
          if (CEvaluationNode::type(pChild2->getType()) == CEvaluationNode::NUMBER)
            {
              pResult = new CEvaluationNodeNumber(CEvaluationNodeNumber::DOUBLE, "1.0");
            }
        }

      // ignore the rest
    }
  else if (CEvaluationNode::type(pChild2->getType()) == CEvaluationNode::NUMBER)
    {
      // 0 and 1
    }

  delete base1;
  delete base2;
  base1 = NULL;
  base2 = NULL;
  return pResult;
}

/**
 * This method makes the elementary elimination on a multiply node.
 */
CEvaluationNode* CNormalTranslation::elementaryEliminationMultiply(CEvaluationNode* pMultiplyNode)
{
  CEvaluationNode* pResult = pMultiplyNode;
  assert(CEvaluationNode::type(pMultiplyNode->getType()) == CEvaluationNode::OPERATOR);
  assert(((CEvaluationNodeOperator::SubType)CEvaluationNode::subType(pMultiplyNode->getType())) == CEvaluationNodeOperator::MULTIPLY);
  CEvaluationNode* pChild1 = dynamic_cast<CEvaluationNode*>(pMultiplyNode->getChild());
  assert(pChild1 != NULL);
  CEvaluationNode* pChild2 = dynamic_cast<CEvaluationNode*>(pChild1->getSibling());
  assert(pChild2 != NULL);
  assert(pChild2->getSibling() == NULL);

  // if one child is NaN, the result is NaN
  if ((CEvaluationNode::type(pChild1->getType()) == CEvaluationNode::CONSTANT &&
       ((CEvaluationNodeConstant::SubType)(CEvaluationNode::subType(pChild1->getType()))) == CEvaluationNodeConstant::_NaN) ||
      (CEvaluationNode::type(pChild2->getType()) == CEvaluationNode::CONSTANT &&
       ((CEvaluationNodeConstant::SubType)(CEvaluationNode::subType(pChild2->getType()))) == CEvaluationNodeConstant::_NaN))
    {
      pResult = new CEvaluationNodeConstant(CEvaluationNodeConstant::_NaN, "NAN");
    }
  // if one child is 0, the result is 0
  else if ((CEvaluationNode::type(pChild1->getType()) == CEvaluationNode::NUMBER &&
            fabs(dynamic_cast<const CEvaluationNodeNumber*>(pChild1)->value()) < ZERO) ||
           (CEvaluationNode::type(pChild2->getType()) == CEvaluationNode::NUMBER &&
            fabs(dynamic_cast<const CEvaluationNodeNumber*>(pChild2)->value()) < ZERO))
    {
      pResult = new CEvaluationNodeNumber(CEvaluationNodeNumber::DOUBLE, "0.0");
    }
  // if one child is 1, the result is the other child
  else if ((CEvaluationNode::type(pChild1->getType()) == CEvaluationNode::NUMBER &&
            fabs(dynamic_cast<const CEvaluationNodeNumber*>(pChild1)->value() - 1.0) < ZERO))
    {
      pResult = pChild2->copyBranch();
    }
  else if ((CEvaluationNode::type(pChild2->getType()) == CEvaluationNode::NUMBER &&
            fabs(dynamic_cast<const CEvaluationNodeNumber*>(pChild2)->value() - 1.0) < ZERO))
    {
      pResult = pChild1->copyBranch();
    }

  return pResult;
}

/**
 * This method makes the elementary elimination on a divide node.
 */
CEvaluationNode* CNormalTranslation::elementaryEliminationDivide(CEvaluationNode* pDivideNode)
{
  CEvaluationNode* pResult = pDivideNode;
  assert(CEvaluationNode::type(pDivideNode->getType()) == CEvaluationNode::OPERATOR);
  assert(((CEvaluationNodeOperator::SubType)CEvaluationNode::subType(pDivideNode->getType())) == CEvaluationNodeOperator::DIVIDE);
  CEvaluationNode* pChild1 = dynamic_cast<CEvaluationNode*>(pDivideNode->getChild());
  assert(pChild1 != NULL);
  CEvaluationNode* pChild2 = dynamic_cast<CEvaluationNode*>(pChild1->getSibling());
  assert(pChild2 != NULL);
  assert(pChild2->getSibling() == NULL);
  // if one of the children is NaN, the result is NaN
  CNormalFraction* base1 = createNormalRepresentation(pChild1);
  CNormalFraction* base2 = createNormalRepresentation(pChild2);

  if ((CEvaluationNode::type(pChild1->getType()) == CEvaluationNode::CONSTANT &&
       ((CEvaluationNodeConstant::SubType)(CEvaluationNode::subType(pChild1->getType()))) == CEvaluationNodeConstant::_NaN) ||
      (CEvaluationNode::type(pChild2->getType()) == CEvaluationNode::CONSTANT &&
       ((CEvaluationNodeConstant::SubType)(CEvaluationNode::subType(pChild2->getType()))) == CEvaluationNodeConstant::_NaN))
    {
      pResult = new CEvaluationNodeConstant(CEvaluationNodeConstant::_NaN, "NAN");
    }
  // the second child is 0, the result is NaN
  else if ((CEvaluationNode::type(pChild2->getType()) == CEvaluationNode::NUMBER &&
            fabs(dynamic_cast<const CEvaluationNodeNumber*>(pChild2)->value()) < ZERO))
    {
      pResult = new CEvaluationNodeConstant(CEvaluationNodeConstant::_NaN, "NaN");
    }
  // if the first child is 0, the result is 0
  else if ((CEvaluationNode::type(pChild1->getType()) == CEvaluationNode::NUMBER &&
            fabs(dynamic_cast<const CEvaluationNodeNumber*>(pChild1)->value()) < ZERO))
    {
      pResult = new CEvaluationNodeNumber(CEvaluationNodeNumber::DOUBLE, "0.0");
    }
  // if both children are the same, the result is 1
  else if (base1->toString() == base2->toString())
    {
      pResult = new CEvaluationNodeNumber(CEvaluationNodeNumber::DOUBLE, "1.0");
    }
  // if the second child is 1, the result is the first child
  else if ((CEvaluationNode::type(pChild2->getType()) == CEvaluationNode::NUMBER &&
            fabs(dynamic_cast<const CEvaluationNodeNumber*>(pChild2)->value() - 1.0) < ZERO))
    {
      pResult = pChild1->copyBranch();
    }

  delete base1;
  delete base2;
  base1 = NULL;
  base2 = NULL;
  return pResult;
}

/**
 * This method makes the elementary elimination on a plus node.
 */
CEvaluationNode* CNormalTranslation::elementaryEliminationPlus(CEvaluationNode* pPlusNode)
{
  CEvaluationNode* pResult = pPlusNode;
  assert(CEvaluationNode::type(pPlusNode->getType()) == CEvaluationNode::OPERATOR);
  assert(((CEvaluationNodeOperator::SubType)CEvaluationNode::subType(pPlusNode->getType())) == CEvaluationNodeOperator::PLUS);
  CEvaluationNode* pChild1 = dynamic_cast<CEvaluationNode*>(pPlusNode->getChild());
  assert(pChild1 != NULL);
  CEvaluationNode* pChild2 = dynamic_cast<CEvaluationNode*>(pChild1->getSibling());
  assert(pChild2 != NULL);
  assert(pChild2->getSibling() == NULL);

  // if one child is NaN, the result is NaN
  if ((CEvaluationNode::type(pChild1->getType()) == CEvaluationNode::CONSTANT &&
       ((CEvaluationNodeConstant::SubType)(CEvaluationNode::subType(pChild1->getType()))) == CEvaluationNodeConstant::_NaN) ||
      (CEvaluationNode::type(pChild2->getType()) == CEvaluationNode::CONSTANT &&
       ((CEvaluationNodeConstant::SubType)(CEvaluationNode::subType(pChild2->getType()))) == CEvaluationNodeConstant::_NaN))
    {
      pResult = new CEvaluationNodeConstant(CEvaluationNodeConstant::_NaN, "NAN");
    }
  // the second child is 0, the result is the first child
  else if ((CEvaluationNode::type(pChild2->getType()) == CEvaluationNode::NUMBER &&
            fabs(dynamic_cast<const CEvaluationNodeNumber*>(pChild2)->value()) < ZERO))
    {
      pResult = pChild1->copyBranch();
    }
  // if the first child is 0, the result is the second child
  else if ((CEvaluationNode::type(pChild1->getType()) == CEvaluationNode::NUMBER &&
            fabs(dynamic_cast<const CEvaluationNodeNumber*>(pChild1)->value()) < ZERO))
    {
      pResult = pChild2->copyBranch();
    }

  return pResult;
}

/**
 * This method makes the elementary elimination on a minus node.
 */
CEvaluationNode* CNormalTranslation::elementaryEliminationMinus(CEvaluationNode* pMinusNode)
{
  CEvaluationNode* pResult = pMinusNode;
  assert(CEvaluationNode::type(pMinusNode->getType()) == CEvaluationNode::OPERATOR);
  assert(((CEvaluationNodeOperator::SubType)CEvaluationNode::subType(pMinusNode->getType())) == CEvaluationNodeOperator::MINUS);
  CEvaluationNode* pChild1 = dynamic_cast<CEvaluationNode*>(pMinusNode->getChild());
  assert(pChild1 != NULL);
  CEvaluationNode* pChild2 = dynamic_cast<CEvaluationNode*>(pChild1->getSibling());
  assert(pChild2 != NULL);
  assert(pChild2->getSibling() == NULL);
  // if one child is NaN, the result is NaN (one could also consider to put
  // the second condition first so that to NaN would cancel each other out
  CNormalFraction* base1 = createNormalRepresentation(pChild1);
  CNormalFraction* base2 = createNormalRepresentation(pChild2);

  if ((CEvaluationNode::type(pChild1->getType()) == CEvaluationNode::CONSTANT &&
       ((CEvaluationNodeConstant::SubType)(CEvaluationNode::subType(pChild1->getType()))) == CEvaluationNodeConstant::_NaN) ||
      (CEvaluationNode::type(pChild2->getType()) == CEvaluationNode::CONSTANT &&
       ((CEvaluationNodeConstant::SubType)(CEvaluationNode::subType(pChild2->getType()))) == CEvaluationNodeConstant::_NaN))
    {
      pResult = new CEvaluationNodeConstant(CEvaluationNodeConstant::_NaN, "NAN");
    }
  // if both nodes are equal, the result is 0.0
  else if (base1->toString() == base2->toString())
    {
      pResult = new CEvaluationNodeNumber(CEvaluationNodeNumber::DOUBLE, "0.0");
    }
  // the second child is 0, the result is the first child
  else if ((CEvaluationNode::type(pChild2->getType()) == CEvaluationNode::NUMBER &&
            fabs(dynamic_cast<const CEvaluationNodeNumber*>(pChild2)->value()) < ZERO))
    {
      pResult = pChild1->copyBranch();
    }
  // if the first child is 0, the result is -1 times the second child
  else if ((CEvaluationNode::type(pChild1->getType()) == CEvaluationNode::NUMBER &&
            fabs(dynamic_cast<const CEvaluationNodeNumber*>(pChild1)->value()) < ZERO))
    {
      pResult = new CEvaluationNodeOperator(CEvaluationNodeOperator::MULTIPLY, "*");
      pResult->addChild(new CEvaluationNodeNumber(CEvaluationNodeNumber::DOUBLE, "-1.0"));
      pResult->addChild(pChild2->copyBranch());
    }

  delete base1;
  delete base2;
  base1 = NULL;
  base2 = NULL;
  return pResult;
}

CEvaluationNode* CNormalTranslation::evaluateNumbers(CEvaluationNode* pOrig)
{
  // if the node is unmodified, return the original node
  // else try to make the modifiactions in place instead of copying the whole
  // subtree
  CEvaluationNode* pResult = pOrig;
  CEvaluationNode* pChild = dynamic_cast<CEvaluationNode*>(pOrig->getChild());
  CEvaluationNode* pLastChild = pOrig;

  while (pChild != NULL)
    {
      CEvaluationNode* pNewChild = CNormalTranslation::evaluateNumbers(pChild);
      assert(pNewChild != NULL);

      if (pNewChild != pChild)
        {
          // remove the old child and add the new one
          pOrig->removeChild(pChild);
          delete pChild;
          pChild = pNewChild;
          pOrig->addChild(pNewChild, pLastChild);
        }

      pLastChild = pChild;
      pChild = dynamic_cast<CEvaluationNode*>(pChild->getSibling());
    }

  if (CEvaluationNode::type(pOrig->getType()) == CEvaluationNode::OPERATOR)
    {
      CEvaluationNode* pChild1 = dynamic_cast<CEvaluationNode*>(pOrig->getChild());
      assert(pChild1 != NULL);
      CEvaluationNode* pChild2 = dynamic_cast<CEvaluationNode*>(pChild1->getSibling());
      assert(pChild2 != NULL);

      switch ((CEvaluationNodeOperator::SubType)CEvaluationNode::subType(pOrig->getType()))
        {
          case CEvaluationNodeOperator::POWER:

            if (CEvaluationNode::type(pChild1->getType()) == CEvaluationNode::NUMBER && CEvaluationNode::type(pChild2->getType()) == CEvaluationNode::NUMBER)
              {
                std::ostringstream os;
                const CEvaluationNodeNumber* pNumberNode1 = dynamic_cast<const CEvaluationNodeNumber*>(pChild1);
                assert(pNumberNode1 != NULL);
                const CEvaluationNodeNumber* pNumberNode2 = dynamic_cast<const CEvaluationNodeNumber*>(pChild2);
                assert(pNumberNode1 != NULL);
                os << pow(pNumberNode1->value(), pNumberNode2->value());
                pResult = new CEvaluationNodeNumber(CEvaluationNodeNumber::DOUBLE, os.str().c_str());
              }

            break;
          case CEvaluationNodeOperator::MODULUS:

            if (CEvaluationNode::type(pChild1->getType()) == CEvaluationNode::NUMBER && CEvaluationNode::type(pChild2->getType()) == CEvaluationNode::NUMBER)
              {
                std::ostringstream os;
                const CEvaluationNodeNumber* pNumberNode1 = dynamic_cast<const CEvaluationNodeNumber*>(pChild1);
                assert(pNumberNode1 != NULL);
                const CEvaluationNodeNumber* pNumberNode2 = dynamic_cast<const CEvaluationNodeNumber*>(pChild2);
                assert(pNumberNode2 != NULL);
                os << ((int)pNumberNode1->value()) % ((int)pNumberNode2->value());
                pResult = new CEvaluationNodeNumber(CEvaluationNodeNumber::DOUBLE, os.str().c_str());
              }

            break;
          case CEvaluationNodeOperator::MULTIPLY:
          case CEvaluationNodeOperator::DIVIDE:
          {
            std::vector<const CEvaluationNode*> multiplications, divisions;
            // multiplications and divisions contain the original nodes,
            // splitProduct doesn't copy nodes
            CNormalTranslation::splitProduct(pResult, multiplications, divisions, false);
            std::set<const CEvaluationNode*> multiplicationNumberNodes;
            unsigned int i, iMax = multiplications.size();
            const CEvaluationNode* pNode = NULL;

            for (i = 0; i < iMax; ++i)
              {
                pNode = multiplications[i];

                if (CEvaluationNode::type(pNode->getType()) == CEvaluationNode::NUMBER)
                  {
                    multiplicationNumberNodes.insert(pNode);
                  }
              }

            std::set<const CEvaluationNode*> divisionNumberNodes;
            iMax = divisions.size();

            for (i = 0; i < iMax; ++i)
              {
                pNode = divisions[i];

                if (CEvaluationNode::type(pNode->getType()) == CEvaluationNode::NUMBER)
                  {
                    divisionNumberNodes.insert(pNode);
                  }
              }

            if ((multiplicationNumberNodes.size() + divisionNumberNodes.size()) > 1)
              {
                // there are at least two number nodes, so we have to evaluate
                // the numbers
                double value = 1.0;
                std::set<const CEvaluationNode*>::iterator it = multiplicationNumberNodes.begin(), endit = multiplicationNumberNodes.end();

                while (it != endit)
                  {
                    value *= (*it)->value();
                    ++it;
                  }

                it = divisionNumberNodes.begin();
                endit = divisionNumberNodes.end();

                while (it != endit)
                  {
                    value /= (*it)->value();
                    ++it;
                  }

                std::vector<CEvaluationNode*> newMultiplications, newDivisions;

                if (fabs((value - 1.0)) >= ZERO)
                  {
                    std::ostringstream os;
                    os.precision(18);

                    if (fabs(value) < 1.0)
                      {
                        os << 1.0 / value;
                        CEvaluationNodeNumber* pEvaluated = new CEvaluationNodeNumber(CEvaluationNodeNumber::DOUBLE, os.str().c_str());
                        newDivisions.push_back(pEvaluated);
                      }
                    else
                      {
                        os << value;
                        CEvaluationNodeNumber* pEvaluated = new CEvaluationNodeNumber(CEvaluationNodeNumber::DOUBLE, os.str().c_str());
                        newMultiplications.push_back(pEvaluated);
                      }
                  }

                // now we have to copy all nonnumber nodes
                std::vector<const CEvaluationNode*>::iterator it2 = multiplications.begin(), endit2 = multiplications.end();

                while (it2 != endit2)
                  {
                    // if the node is not in multiplicationNumberNodes, we copy
                    // it
                    it = multiplicationNumberNodes.find(*it2);

                    if (it == multiplicationNumberNodes.end())
                      {
                        newMultiplications.push_back((*it2)->copyBranch());
                      }

                    ++it2;
                  }

                it2 = divisions.begin();
                endit2 = divisions.end();

                while (it2 != endit2)
                  {
                    // if the node is not in multiplicationNumberNodes, we copy
                    // it
                    it = divisionNumberNodes.find(*it2);

                    if (it == divisionNumberNodes.end())
                      {
                        newDivisions.push_back((*it2)->copyBranch());
                      }

                    ++it2;
                  }

                // now we create a new result node from the newMultiplications
                // and newDivisions
                if (newMultiplications.empty())
                  {
                    pResult = new CEvaluationNodeNumber(CEvaluationNodeNumber::DOUBLE, "1.0");
                  }
                else
                  {
                    pResult = CNormalTranslation::createChain(&CNormalTranslation::TIMES_NODE, &CNormalTranslation::ONE_NODE, newMultiplications);
                  }

                if (!newDivisions.empty())
                  {
                    CEvaluationNode* pTmpNode = new CEvaluationNodeOperator(CEvaluationNodeOperator::DIVIDE, "/");
                    pTmpNode->addChild(pResult);
                    pResult = pTmpNode;
                    pTmpNode = CNormalTranslation::createChain(&CNormalTranslation::TIMES_NODE, &CNormalTranslation::ONE_NODE, newDivisions);
                    pResult->addChild(pTmpNode);
                  }
              }
          }
          break;
          case CEvaluationNodeOperator::PLUS:
          case CEvaluationNodeOperator::MINUS:
          {
            std::vector<CEvaluationNode*> additions, subtractions;
            // splitSum copies the nodes that are returned
            CNormalTranslation::splitSum(pResult, additions, subtractions, false);
            CNormalTranslation::swapNegativeNumbers(additions, subtractions);
            std::set<const CEvaluationNode*> additionNumberNodes;
            unsigned int i, iMax = additions.size();
            const CEvaluationNode* pNode = NULL;

            for (i = 0; i < iMax; ++i)
              {
                pNode = additions[i];

                if (CEvaluationNode::type(pNode->getType()) == CEvaluationNode::NUMBER)
                  {
                    additionNumberNodes.insert(pNode);
                  }
              }

            std::set<const CEvaluationNode*> subtractionNumberNodes;
            iMax = subtractions.size();

            for (i = 0; i < iMax; ++i)
              {
                pNode = subtractions[i];

                if (CEvaluationNode::type(pNode->getType()) == CEvaluationNode::NUMBER)
                  {
                    subtractionNumberNodes.insert(pNode);
                  }
              }

            if ((additionNumberNodes.size() + subtractionNumberNodes.size()) > 1)
              {
                // there are at least two number nodes, so we have to evaluate
                // the numbers
                double value = 0.0;
                std::set<const CEvaluationNode*>::const_iterator it = additionNumberNodes.begin(), endit = additionNumberNodes.end();

                while (it != endit)
                  {
                    value += (*it)->value();
                    ++it;
                  }

                it = subtractionNumberNodes.begin();
                endit = subtractionNumberNodes.end();

                while (it != endit)
                  {
                    value -= (*it)->value();
                    ++it;
                  }

                std::vector<CEvaluationNode*> newAdditions, newSubtractions;

                if (fabs(value) >= ZERO)
                  {
                    std::ostringstream os;
                    os.precision(18);

                    if (value < 0.0)
                      {
                        os << -1.0 * value;
                        CEvaluationNodeNumber* pEvaluated = new CEvaluationNodeNumber(CEvaluationNodeNumber::DOUBLE, os.str().c_str());
                        newSubtractions.push_back(pEvaluated);
                      }
                    else
                      {
                        os << value;
                        CEvaluationNodeNumber* pEvaluated = new CEvaluationNodeNumber(CEvaluationNodeNumber::DOUBLE, os.str().c_str());
                        newAdditions.push_back(pEvaluated);
                      }
                  }

                // now we have to copy all nonnumber nodes
                std::vector<CEvaluationNode*>::const_iterator it2 = additions.begin(), endit2 = additions.end();

                while (it2 != endit2)
                  {
                    // if the node is not in additionNumberNodes, we copy
                    // it
                    it = additionNumberNodes.find(*it2);

                    if (it == additionNumberNodes.end())
                      {
                        newAdditions.push_back(*it2);
                      }
                    else
                      {
                        // delete the original node that was created by splitSum
                        delete *it2;
                      }

                    ++it2;
                  }

                it2 = subtractions.begin();
                endit2 = subtractions.end();

                while (it2 != endit2)
                  {
                    // if the node is not in subtractionNumberNodes, we copy
                    // it
                    it = subtractionNumberNodes.find(*it2);

                    if (it == subtractionNumberNodes.end())
                      {
                        newSubtractions.push_back(*it2);
                      }
                    else
                      {
                        // delete the original node that was created by splitSum
                        delete *it2;
                      }

                    ++it2;
                  }

                // now we create a new result node from the newAdditions
                // and newSubtractions
                if (newAdditions.empty())
                  {
                    if (newSubtractions.empty())
                      {
                        pResult = new CEvaluationNodeNumber(CEvaluationNodeNumber::DOUBLE, "0.0");
                      }
                  }
                else
                  {
                    pResult = CNormalTranslation::createChain(&CNormalTranslation::PLUS_NODE, &CNormalTranslation::ZERO_NODE, newAdditions);
                  }

                if (!newSubtractions.empty())
                  {
                    if (newAdditions.empty())
                      {
                        if (newSubtractions.size() == 1 && CEvaluationNode::type(newSubtractions[0]->getType()) == CEvaluationNode::NUMBER)
                          {
                            std::ostringstream os;
                            os.precision(18);
                            os << -1.0 * newSubtractions[0]->value();
                            pResult = new CEvaluationNodeNumber(CEvaluationNodeNumber::DOUBLE, os.str().c_str());
                            delete newSubtractions[0];
                          }
                        else
                          {
                            pResult = new CEvaluationNodeFunction(CEvaluationNodeFunction::MINUS, "-");
                            CEvaluationNode* pTmpNode = CNormalTranslation::createChain(&CNormalTranslation::PLUS_NODE, &CNormalTranslation::ZERO_NODE, newSubtractions);
                            pResult->addChild(pTmpNode);
                          }
                      }
                    else
                      {
                        CEvaluationNode* pTmpNode = new CEvaluationNodeOperator(CEvaluationNodeOperator::MINUS, "-");
                        pTmpNode->addChild(pResult);
                        pResult = pTmpNode;
                        pTmpNode = CNormalTranslation::createChain(&CNormalTranslation::PLUS_NODE, &CNormalTranslation::ZERO_NODE, newSubtractions);
                        pResult->addChild(pTmpNode);
                      }
                  }
              }
            else
              {
                // delete all nodes in additions and subtractions
                unsigned int i, iMax = additions.size();

                for (i = 0; i < iMax; ++i)
                  {
                    delete additions[i];
                  }

                iMax = subtractions.size();

                for (i = 0; i < iMax; ++i)
                  {
                    delete subtractions[i];
                  }
              }
          }
          break;
          case CEvaluationNodeOperator::INVALID:
            break;
        }
    }

  return pResult;
}

/**
 * This method removes nested power nodes, e.g. (a^b)^c -> a^(b*c)
 */
CEvaluationNode* CNormalTranslation::eliminateNestedPowers(const CEvaluationNode* pOrig)
{
  CEvaluationNode* pResult = NULL;
  std::vector<CEvaluationNode*> children;
  const CEvaluationNode* pChild = dynamic_cast<const CEvaluationNode*>(pOrig->getChild());

  while (pChild != NULL)
    {
      children.push_back(CNormalTranslation::eliminateNestedPowers(pChild));
      pChild = dynamic_cast<const CEvaluationNode*>(pChild->getSibling());
    }

  if (CEvaluationNode::type(pOrig->getType()) == CEvaluationNode::OPERATOR &&
      ((CEvaluationNodeOperator::SubType)CEvaluationNode::subType(pOrig->getType())) == CEvaluationNodeOperator::POWER)
    {
      // check if the first child is also a power node
      assert(children.size() == 2);

      if (CEvaluationNode::type(children[0]->getType()) == CEvaluationNode::OPERATOR &&
          ((CEvaluationNodeOperator::SubType)CEvaluationNode::subType(children[0]->getType())) == CEvaluationNodeOperator::POWER)
        {
          pResult = new CEvaluationNodeOperator(CEvaluationNodeOperator::POWER, "^");
          const CEvaluationNode* pChild = dynamic_cast<CEvaluationNode*>(children[0]->getChild());
          assert(pChild != NULL);
          pResult->addChild(pChild->copyBranch());
          CEvaluationNode* pMult = new CEvaluationNodeOperator(CEvaluationNodeOperator::MULTIPLY, "*");
          pChild = dynamic_cast<const CEvaluationNode*>(pChild->getSibling());
          assert(pChild != NULL);
          pMult->addChild(pChild->copyBranch());
          pMult->addChild(children[1]);
          // delete the unused child
          delete children[0];
          pResult->addChild(pMult);
        }
    }

  if (pResult == NULL)
    {
      pResult = pOrig->copyNode(children);
    }

  return pResult;
}

/**
 * This method splits a product into the individual elements
 */
void CNormalTranslation::splitProduct(const CEvaluationNode* pRoot, std::vector<const CEvaluationNode*>& multiplications, std::vector<const CEvaluationNode*>& divisions, bool division)
{
  if (CEvaluationNode::type(pRoot->getType()) == CEvaluationNode::OPERATOR && ((CEvaluationNodeOperator::SubType)CEvaluationNode::subType(pRoot->getType()) == CEvaluationNodeOperator::MULTIPLY || (CEvaluationNodeOperator::SubType)CEvaluationNode::subType(pRoot->getType()) == CEvaluationNodeOperator::DIVIDE))
    {
      const CEvaluationNode* pChild1 = dynamic_cast<const CEvaluationNode*>(pRoot->getChild());
      assert(pChild1 != NULL);
      const CEvaluationNode* pChild2 = dynamic_cast<const CEvaluationNode*>(pChild1->getSibling());
      assert(pChild2 != NULL);
      assert(pChild2->getSibling() == NULL);

      if (((CEvaluationNodeOperator::SubType)CEvaluationNode::subType(pRoot->getType())) == CEvaluationNodeOperator::MULTIPLY)
        {
          if (CEvaluationNode::type(pChild1->getType()) == CEvaluationNode::OPERATOR &&
              (((CEvaluationNodeOperator::SubType)CEvaluationNode::subType(pChild1->getType())) == CEvaluationNodeOperator::MULTIPLY ||
               ((CEvaluationNodeOperator::SubType)CEvaluationNode::subType(pChild1->getType())) == CEvaluationNodeOperator::DIVIDE))
            {
              CNormalTranslation::splitProduct(pChild1, multiplications, divisions, division);
            }
          else
            {
              if (division == false)
                {
                  multiplications.push_back(pChild1);
                }
              else
                {
                  divisions.push_back(pChild1);
                }
            }

          if (CEvaluationNode::type(pChild2->getType()) == CEvaluationNode::OPERATOR &&
              (((CEvaluationNodeOperator::SubType)CEvaluationNode::subType(pChild2->getType())) == CEvaluationNodeOperator::MULTIPLY ||
               ((CEvaluationNodeOperator::SubType)CEvaluationNode::subType(pChild2->getType())) == CEvaluationNodeOperator::DIVIDE))
            {
              CNormalTranslation::splitProduct(pChild2, multiplications, divisions, division);
            }
          else
            {
              if (division == false)
                {
                  multiplications.push_back(pChild2);
                }
              else
                {
                  divisions.push_back(pChild2);
                }
            }
        }
      else if (((CEvaluationNodeOperator::SubType)CEvaluationNode::subType(pRoot->getType())) == CEvaluationNodeOperator::DIVIDE)
        {
          if (CEvaluationNode::type(pChild1->getType()) == CEvaluationNode::OPERATOR &&
              (((CEvaluationNodeOperator::SubType)CEvaluationNode::subType(pChild1->getType())) == CEvaluationNodeOperator::MULTIPLY ||
               ((CEvaluationNodeOperator::SubType)CEvaluationNode::subType(pChild1->getType())) == CEvaluationNodeOperator::DIVIDE))
            {
              CNormalTranslation::splitProduct(pChild1, multiplications, divisions, division);
            }
          else
            {
              if (division == false)
                {
                  multiplications.push_back(pChild1);
                }
              else
                {
                  divisions.push_back(pChild1);
                }
            }

          if (CEvaluationNode::type(pChild2->getType()) == CEvaluationNode::OPERATOR &&
              (((CEvaluationNodeOperator::SubType)CEvaluationNode::subType(pChild2->getType())) == CEvaluationNodeOperator::MULTIPLY ||
               ((CEvaluationNodeOperator::SubType)CEvaluationNode::subType(pChild2->getType())) == CEvaluationNodeOperator::DIVIDE))
            {
              CNormalTranslation::splitProduct(pChild2, multiplications, divisions, !division);
            }
          else
            {
              if (division == false)
                {
                  divisions.push_back(pChild2);
                }
              else
                {
                  multiplications.push_back(pChild2);
                }
            }
        }
    }
  else
    {
      multiplications.push_back(pRoot);
    }
}

/**
 * This method splits a sum into the individual elements
 * The returned nodes are part of the original node and not copies.
 */
void CNormalTranslation::splitSum(const CEvaluationNode* pRoot, std::vector<const CEvaluationNode*>& additions, std::vector<const CEvaluationNode*>& subtractions, bool minus)
{
  // TODO this method might save some copy/delete cycles if the test for
  // TODO negative number was done before making copies of children and
  // TODO inserting them
  // TODO this would also simplify the code since the test would be put
  // TODO into a separate function
  if (CEvaluationNode::type(pRoot->getType()) == CEvaluationNode::OPERATOR && ((CEvaluationNodeOperator::SubType)CEvaluationNode::subType(pRoot->getType()) == CEvaluationNodeOperator::PLUS || (CEvaluationNodeOperator::SubType)CEvaluationNode::subType(pRoot->getType()) == CEvaluationNodeOperator::MINUS))
    {
      const CEvaluationNode* pChild1 = dynamic_cast<const CEvaluationNode*>(pRoot->getChild());
      assert(pChild1 != NULL);
      const CEvaluationNode* pChild2 = dynamic_cast<const CEvaluationNode*>(pChild1->getSibling());
      assert(pChild2 != NULL);
      assert(pChild2->getSibling() == NULL);

      if (((CEvaluationNodeOperator::SubType)CEvaluationNode::subType(pRoot->getType())) == CEvaluationNodeOperator::PLUS)
        {
          if (CEvaluationNode::type(pChild1->getType()) == CEvaluationNode::OPERATOR &&
              (((CEvaluationNodeOperator::SubType)CEvaluationNode::subType(pChild1->getType())) == CEvaluationNodeOperator::PLUS ||
               ((CEvaluationNodeOperator::SubType)CEvaluationNode::subType(pChild1->getType())) == CEvaluationNodeOperator::MINUS))
            {
              CNormalTranslation::splitSum(pChild1, additions, subtractions, minus);
            }
          else
            {
              if (minus == false)
                {
                  additions.push_back(pChild1);
                }
              else
                {
                  subtractions.push_back(pChild1);
                }
            }

          if (CEvaluationNode::type(pChild2->getType()) == CEvaluationNode::OPERATOR &&
              (((CEvaluationNodeOperator::SubType)CEvaluationNode::subType(pChild2->getType())) == CEvaluationNodeOperator::PLUS ||
               ((CEvaluationNodeOperator::SubType)CEvaluationNode::subType(pChild2->getType())) == CEvaluationNodeOperator::MINUS))
            {
              CNormalTranslation::splitSum(pChild2, additions, subtractions, minus);
            }
          else
            {
              if (minus == false)
                {
                  additions.push_back(pChild2);
                }
              else
                {
                  subtractions.push_back(pChild2);
                }
            }
        }
      else if (((CEvaluationNodeOperator::SubType)CEvaluationNode::subType(pRoot->getType())) == CEvaluationNodeOperator::MINUS)
        {
          if (CEvaluationNode::type(pChild1->getType()) == CEvaluationNode::OPERATOR &&
              (((CEvaluationNodeOperator::SubType)CEvaluationNode::subType(pChild1->getType())) == CEvaluationNodeOperator::PLUS ||
               ((CEvaluationNodeOperator::SubType)CEvaluationNode::subType(pChild1->getType())) == CEvaluationNodeOperator::MINUS))
            {
              CNormalTranslation::splitSum(pChild1, additions, subtractions, minus);
            }
          else
            {
              if (minus == false)
                {
                  additions.push_back(pChild1);
                }
              else
                {
                  subtractions.push_back(pChild1);
                }
            }

          if (CEvaluationNode::type(pChild2->getType()) == CEvaluationNode::OPERATOR &&
              (((CEvaluationNodeOperator::SubType)CEvaluationNode::subType(pChild2->getType())) == CEvaluationNodeOperator::PLUS ||
               ((CEvaluationNodeOperator::SubType)CEvaluationNode::subType(pChild2->getType())) == CEvaluationNodeOperator::MINUS))
            {
              CNormalTranslation::splitSum(pChild2, additions, subtractions, !minus);
            }
          else
            {
              if (minus == false)
                {
                  subtractions.push_back(pChild2);
                }
              else
                {
                  additions.push_back(pChild2);
                }
            }
        }
    }
  else
    {
      additions.push_back(pRoot);
    }
}

/**
 * This method splits a sum into the individual elements
 * The returned nodes are copies of the original.
 */
void CNormalTranslation::splitSum(const CEvaluationNode* pRoot, std::vector<CEvaluationNode*>& additions, std::vector<CEvaluationNode*>& subtractions, bool minus)
{
  std::vector<const CEvaluationNode*> tmpAdditions, tmpSubtractions;
  CNormalTranslation::splitSum(pRoot, tmpAdditions, tmpSubtractions, minus);
  unsigned int i, iMax = tmpAdditions.size();
  additions.reserve(iMax);

  for (i = 0; i < iMax; ++i)
    {
      additions.push_back(tmpAdditions[i]->copyBranch());
    }

  iMax = tmpSubtractions.size();
  subtractions.reserve(iMax);

  for (i = 0; i < iMax; ++i)
    {
      subtractions.push_back(tmpSubtractions[i]->copyBranch());
    }

  // TODO the code below was part of the old splitSum method that has largely
  // TODO been replaced by the new method that doesn't copy the nodes and is
  // TODO therefor faster
  // TODO If this code below is removed, the expression comparison for ordered
  // TODO bi bi goes into an endless loop, so I know there is still a bug somewhere
  // TODO in another routine which has to be fixed.
  // TODO actually the code below should be obsolete since swapNegativeNumbers
  // does this now.

  // check for negative numbers in additions and add them to subtractions
  // likewise check for negative numbers in substractions and add them to
  // additions
  // do the same for multiplications with a negative number
  std::vector<CEvaluationNode*>::iterator it = additions.begin(), endit = additions.end();

  while (it != endit)
    {
      if (CEvaluationNode::type((*it)->getType()) == CEvaluationNode::NUMBER && (*it)->value() < 0.0)
        {
          std::ostringstream os;
          os.precision(18);
          os << (*it)->value() * -1.0;
          CEvaluationNode* pTmpNumber = new CEvaluationNodeNumber(CEvaluationNodeNumber::DOUBLE, os.str().c_str());
          subtractions.push_back(pTmpNumber);
          delete *it;
          it = additions.erase(it);
          endit = additions.end();
          continue;
        }
      else if (CEvaluationNode::type((*it)->getType()) == CEvaluationNode::OPERATOR && (CEvaluationNodeOperator::SubType)CEvaluationNode::subType((*it)->getType()) == CEvaluationNodeOperator::MULTIPLY)
        {
          // actually there should be code that tests if both are negative
          // numbers
          if ((CEvaluationNode::type(dynamic_cast<const CEvaluationNode*>((*it)->getChild())->getType()) == CEvaluationNode::NUMBER && dynamic_cast<const CEvaluationNodeNumber*>((*it)->getChild())->value() < 0.0))
            {
              if (fabs(dynamic_cast<const CEvaluationNodeNumber*>((*it)->getChild())->value()) - 1.0 < ZERO)
                {
                  subtractions.push_back(dynamic_cast<const CEvaluationNode*>((*it)->getChild()->getSibling())->copyBranch());
                  delete *it;
                  it = additions.erase(it);
                  endit = additions.end();
                  continue;
                }
              else
                {
                  CEvaluationNode* pTmp = new CEvaluationNodeOperator(CEvaluationNodeOperator::MULTIPLY, "*");
                  std::ostringstream os;
                  os.precision(18);
                  os << static_cast<const CEvaluationNodeNumber*>((*it)->getChild())->value() * -1.0;
                  pTmp->addChild(new CEvaluationNodeNumber(CEvaluationNodeNumber::DOUBLE, os.str().c_str()));
                  pTmp->addChild(dynamic_cast<const CEvaluationNode*>((*it)->getChild()->getSibling())->copyBranch());
                  subtractions.push_back(pTmp);
                  delete *it;
                  it = additions.erase(it);
                  endit = additions.end();
                  continue;
                }
            }
          else if (CEvaluationNode::type(dynamic_cast<const CEvaluationNode*>((*it)->getChild()->getSibling())->getType()) == CEvaluationNode::NUMBER && dynamic_cast<const CEvaluationNodeNumber*>((*it)->getChild()->getSibling())->value() < 0.0)
            {
              if (fabs(dynamic_cast<const CEvaluationNodeNumber*>((*it)->getChild()->getSibling())->value()) - 1.0 < ZERO)
                {
                  subtractions.push_back(dynamic_cast<const CEvaluationNode*>((*it)->getChild())->copyBranch());
                  delete *it;
                  it = additions.erase(it);
                  endit = additions.end();
                  continue;
                }
              else
                {
                  CEvaluationNode* pTmp = new CEvaluationNodeOperator(CEvaluationNodeOperator::MULTIPLY, "*");
                  pTmp->addChild(dynamic_cast<const CEvaluationNode*>((*it)->getChild())->copyBranch());
                  std::ostringstream os;
                  os.precision(18);
                  os << static_cast<const CEvaluationNodeNumber*>((*it)->getChild()->getSibling())->value() * -1.0;
                  pTmp->addChild(new CEvaluationNodeNumber(CEvaluationNodeNumber::DOUBLE, os.str().c_str()));
                  subtractions.push_back(pTmp);
                  delete *it;
                  it = additions.erase(it);
                  endit = additions.end();
                  continue;
                }
            }
          else
            {
              ++it;
            }
        }
      else
        {
          ++it;
        }
    }

  it = subtractions.begin();
  endit = subtractions.end();

  while (it != endit)
    {
      if (CEvaluationNode::type((*it)->getType()) == CEvaluationNode::NUMBER && dynamic_cast<const CEvaluationNodeNumber*>(*it)->value() < 0.0)
        {
          std::ostringstream os;
          os.precision(18);
          os << (*it)->value() * -1.0;
          CEvaluationNode* pTmpNumber = new CEvaluationNodeNumber(CEvaluationNodeNumber::DOUBLE, os.str().c_str());
          additions.push_back(pTmpNumber);
          delete *it;
          it = subtractions.erase(it);
          endit = subtractions.end();
          continue;
        }
      else if (CEvaluationNode::type((*it)->getType()) == CEvaluationNode::OPERATOR && (CEvaluationNodeOperator::SubType)CEvaluationNode::subType((*it)->getType()) == CEvaluationNodeOperator::MULTIPLY)
        {
          // actually there should be code that tests if both are negative
          // numbers
          if ((CEvaluationNode::type(dynamic_cast<const CEvaluationNode*>((*it)->getChild())->getType()) == CEvaluationNode::NUMBER && dynamic_cast<const CEvaluationNodeNumber*>((*it)->getChild())->value() < 0.0))
            {
              if (fabs(dynamic_cast<const CEvaluationNodeNumber*>((*it)->getChild())->value()) - 1.0 < ZERO)
                {
                  additions.push_back(dynamic_cast<const CEvaluationNode*>((*it)->getChild()->getSibling())->copyBranch());
                  delete *it;
                  it = subtractions.erase(it);
                  endit = subtractions.end();
                  continue;
                }
              else
                {
                  CEvaluationNode* pTmp = new CEvaluationNodeOperator(CEvaluationNodeOperator::MULTIPLY, "*");
                  std::ostringstream os;
                  os.precision(18);
                  os << static_cast<const CEvaluationNodeNumber*>((*it)->getChild())->value() * -1.0;
                  pTmp->addChild(new CEvaluationNodeNumber(CEvaluationNodeNumber::DOUBLE, os.str().c_str()));
                  pTmp->addChild(dynamic_cast<const CEvaluationNode*>((*it)->getChild()->getSibling())->copyBranch());
                  additions.push_back(pTmp);
                  delete *it;
                  it = subtractions.erase(it);
                  endit = subtractions.end();
                  continue;
                }
            }
          else if (CEvaluationNode::type(dynamic_cast<const CEvaluationNode*>((*it)->getChild()->getSibling())->getType()) == CEvaluationNode::NUMBER && dynamic_cast<const CEvaluationNodeNumber*>((*it)->getChild()->getSibling())->value() < 0.0)
            {
              if (fabs(dynamic_cast<const CEvaluationNodeNumber*>((*it)->getChild()->getSibling())->value()) - 1.0 < ZERO)
                {
                  additions.push_back(dynamic_cast<const CEvaluationNode*>((*it)->getChild())->copyBranch());
                  delete *it;
                  it = subtractions.erase(it);
                  endit = subtractions.end();
                  continue;
                }
              else
                {
                  CEvaluationNode* pTmp = new CEvaluationNodeOperator(CEvaluationNodeOperator::MULTIPLY, "*");
                  pTmp->addChild(dynamic_cast<const CEvaluationNode*>((*it)->getChild())->copyBranch());
                  std::ostringstream os;
                  os.precision(18);
                  os << static_cast<const CEvaluationNodeNumber*>((*it)->getChild()->getSibling())->value() * -1.0;
                  pTmp->addChild(new CEvaluationNodeNumber(CEvaluationNodeNumber::DOUBLE, os.str().c_str()));
                  additions.push_back(pTmp);
                  delete *it;
                  it = subtractions.erase(it);
                  endit = subtractions.end();
                  continue;
                }
            }
          else
            {
              ++it;
            }
        }
      else
        {
          ++it;
        }
    }
}

/**
 * This method splits a sum into the individual elements
void CNormalTranslation::splitSum(const CEvaluationNode* pRoot, std::vector<CEvaluationNode*>& additions, std::vector<CEvaluationNode*>& subtractions, bool minus)
{
  // TODO this method might save some copy/delete cycles if the test for
  // TODO negative number was done before making copies of children and
  // TODO inserting them
  // TODO this would also simplify the code since the test would be put
  // TODO into a separate function
  if (CEvaluationNode::type(pRoot->getType()) == CEvaluationNode::OPERATOR && ((CEvaluationNodeOperator::SubType)CEvaluationNode::subType(pRoot->getType()) == CEvaluationNodeOperator::PLUS || (CEvaluationNodeOperator::SubType)CEvaluationNode::subType(pRoot->getType()) == CEvaluationNodeOperator::MINUS))
    {
      const CEvaluationNode* pChild1 = dynamic_cast<const CEvaluationNode*>(pRoot->getChild());
      assert(pChild1 != NULL);
      const CEvaluationNode* pChild2 = dynamic_cast<const CEvaluationNode*>(pChild1->getSibling());
      assert(pChild2 != NULL);
      assert(pChild2->getSibling() == NULL);
      if (((CEvaluationNodeOperator::SubType)CEvaluationNode::subType(pRoot->getType())) == CEvaluationNodeOperator::PLUS)
        {
          if (CEvaluationNode::type(pChild1->getType()) == CEvaluationNode::OPERATOR &&
              (((CEvaluationNodeOperator::SubType)CEvaluationNode::subType(pChild1->getType())) == CEvaluationNodeOperator::PLUS ||
               ((CEvaluationNodeOperator::SubType)CEvaluationNode::subType(pChild1->getType())) == CEvaluationNodeOperator::MINUS))
            {
              CNormalTranslation::splitSum(pChild1, additions, subtractions, minus);
            }
          else
            {
              if (minus == false)
                {
                  additions.push_back(pChild1->copyBranch());
                }
              else
                {
                  subtractions.push_back(pChild1->copyBranch());
                }
            }
          if (CEvaluationNode::type(pChild2->getType()) == CEvaluationNode::OPERATOR &&
              (((CEvaluationNodeOperator::SubType)CEvaluationNode::subType(pChild2->getType())) == CEvaluationNodeOperator::PLUS ||
               ((CEvaluationNodeOperator::SubType)CEvaluationNode::subType(pChild2->getType())) == CEvaluationNodeOperator::MINUS))
            {
              CNormalTranslation::splitSum(pChild2, additions, subtractions, minus);
            }
          else
            {
              if (minus == false)
                {
                  additions.push_back(pChild2->copyBranch());
                }
              else
                {
                  subtractions.push_back(pChild2->copyBranch());
                }
            }
        }
      else if (((CEvaluationNodeOperator::SubType)CEvaluationNode::subType(pRoot->getType())) == CEvaluationNodeOperator::MINUS)
        {
          if (CEvaluationNode::type(pChild1->getType()) == CEvaluationNode::OPERATOR &&
              (((CEvaluationNodeOperator::SubType)CEvaluationNode::subType(pChild1->getType())) == CEvaluationNodeOperator::PLUS ||
               ((CEvaluationNodeOperator::SubType)CEvaluationNode::subType(pChild1->getType())) == CEvaluationNodeOperator::MINUS))
            {
              CNormalTranslation::splitSum(pChild1, additions, subtractions, minus);
            }
          else
            {
              if (minus == false)
                {
                  additions.push_back(pChild1->copyBranch());
                }
              else
                {
                  subtractions.push_back(pChild1->copyBranch());
                }
            }
          if (CEvaluationNode::type(pChild2->getType()) == CEvaluationNode::OPERATOR &&
              (((CEvaluationNodeOperator::SubType)CEvaluationNode::subType(pChild2->getType())) == CEvaluationNodeOperator::PLUS ||
               ((CEvaluationNodeOperator::SubType)CEvaluationNode::subType(pChild2->getType())) == CEvaluationNodeOperator::MINUS))
            {
              CNormalTranslation::splitSum(pChild2, additions, subtractions, !minus);
            }
          else
            {
              if (minus == false)
                {
                  subtractions.push_back(pChild2->copyBranch());
                }
              else
                {
                  additions.push_back(pChild2->copyBranch());
                }
            }
        }
    }
  else
    {
      additions.push_back(pRoot->copyBranch());
    }
  // check for negative numbers in additions and add them to subtractions
  // likewise check for negative numbers in substractions and add them to
  // additions
  // do the same for multiplications with a negative number
  std::vector<CEvaluationNode*>::iterator it = additions.begin();
  while (it != additions.end())
    {
      if (CEvaluationNode::type((*it)->getType()) == CEvaluationNode::NUMBER && dynamic_cast<const CEvaluationNodeNumber*>(*it)->value() < 0.0)
        {
          std::ostringstream os;
          os.precision(18);
          os << static_cast<CEvaluationNodeNumber*>(*it)->value() * -1.0;
          CEvaluationNode* pTmpNumber = new CEvaluationNodeNumber(CEvaluationNodeNumber::DOUBLE, os.str().c_str());
          subtractions.push_back(pTmpNumber);
          delete *it;
          it = additions.erase(it);
        }
      else if (CEvaluationNode::type((*it)->getType()) == CEvaluationNode::OPERATOR && (CEvaluationNodeOperator::SubType)CEvaluationNode::subType((*it)->getType()) == CEvaluationNodeOperator::MULTIPLY)
        {
          // actually there should be code that tests if both are negative
          // numbers
          if ((CEvaluationNode::type(dynamic_cast<const CEvaluationNode*>((*it)->getChild())->getType()) == CEvaluationNode::NUMBER && dynamic_cast<const CEvaluationNodeNumber*>((*it)->getChild())->value() < 0.0))
            {
              if(fabs(dynamic_cast<const CEvaluationNodeNumber*>((*it)->getChild())->value()) - 1.0 < ZERO)
                {
                  subtractions.push_back(dynamic_cast<const CEvaluationNode*>((*it)->getChild()->getSibling())->copyBranch());
                  delete *it;
                  it = additions.erase(it);
                }
              else
                {
                  CEvaluationNode* pTmp = new CEvaluationNodeOperator(CEvaluationNodeOperator::MULTIPLY, "*");
                  std::ostringstream os;
                  os.precision(18);
                  os << static_cast<const CEvaluationNodeNumber*>((*it)->getChild())->value() * -1.0;
                  pTmp->addChild(new CEvaluationNodeNumber(CEvaluationNodeNumber::DOUBLE, os.str().c_str()));
                  pTmp->addChild(dynamic_cast<const CEvaluationNode*>((*it)->getChild()->getSibling())->copyBranch());
                  subtractions.push_back(pTmp);
                  delete *it;
                  it = additions.erase(it);
                }
            }
          else if (CEvaluationNode::type(dynamic_cast<const CEvaluationNode*>((*it)->getChild()->getSibling())->getType()) == CEvaluationNode::NUMBER && dynamic_cast<const CEvaluationNodeNumber*>((*it)->getChild()->getSibling())->value() < 0.0)
            {
              if(fabs(dynamic_cast<const CEvaluationNodeNumber*>((*it)->getChild()->getSibling())->value()) - 1.0 < ZERO)
                {
                  subtractions.push_back(dynamic_cast<const CEvaluationNode*>((*it)->getChild())->copyBranch());
                  delete *it;
                  it = additions.erase(it);
                }
              else
                {
                  CEvaluationNode* pTmp = new CEvaluationNodeOperator(CEvaluationNodeOperator::MULTIPLY, "*");
                  pTmp->addChild(dynamic_cast<const CEvaluationNode*>((*it)->getChild())->copyBranch());
                  std::ostringstream os;
                  os.precision(18);
                  os << static_cast<const CEvaluationNodeNumber*>((*it)->getChild()->getSibling())->value() * -1.0;
                  pTmp->addChild(new CEvaluationNodeNumber(CEvaluationNodeNumber::DOUBLE, os.str().c_str()));
                  subtractions.push_back(pTmp);
                  delete *it;
                  it = additions.erase(it);
                }
            }
          else
            {
              ++it;
            }
        }
      else
        {
          ++it;
        }
    }
  it = subtractions.begin();
  while (it != subtractions.end())
    {
      if (CEvaluationNode::type((*it)->getType()) == CEvaluationNode::NUMBER && dynamic_cast<const CEvaluationNodeNumber*>(*it)->value() < 0.0)
        {
          std::ostringstream os;
          os.precision(18);
          os << static_cast<CEvaluationNodeNumber*>(*it)->value() * -1.0;
          CEvaluationNode* pTmpNumber = new CEvaluationNodeNumber(CEvaluationNodeNumber::DOUBLE, os.str().c_str());
          additions.push_back(pTmpNumber);
          delete *it;
          it = subtractions.erase(it);
        }
      else if (CEvaluationNode::type((*it)->getType()) == CEvaluationNode::OPERATOR && (CEvaluationNodeOperator::SubType)CEvaluationNode::subType((*it)->getType()) == CEvaluationNodeOperator::MULTIPLY)
        {
          // actually there should be code that tests if both are negative
          // numbers
          if ((CEvaluationNode::type(dynamic_cast<const CEvaluationNode*>((*it)->getChild())->getType()) == CEvaluationNode::NUMBER && dynamic_cast<const CEvaluationNodeNumber*>((*it)->getChild())->value() < 0.0))
            {
              if(fabs(dynamic_cast<const CEvaluationNodeNumber*>((*it)->getChild())->value()) - 1.0 < ZERO)
                {
                  additions.push_back(dynamic_cast<const CEvaluationNode*>((*it)->getChild()->getSibling())->copyBranch());
                  delete *it;
                  it = subtractions.erase(it);
                }
              else
                {
                  CEvaluationNode* pTmp = new CEvaluationNodeOperator(CEvaluationNodeOperator::MULTIPLY, "*");
                  std::ostringstream os;
                  os.precision(18);
                  os << static_cast<const CEvaluationNodeNumber*>((*it)->getChild())->value() * -1.0;
                  pTmp->addChild(new CEvaluationNodeNumber(CEvaluationNodeNumber::DOUBLE, os.str().c_str()));
                  pTmp->addChild(dynamic_cast<const CEvaluationNode*>((*it)->getChild()->getSibling())->copyBranch());
                  additions.push_back(pTmp);
                  delete *it;
                  it = subtractions.erase(it);
                }
            }
          else if (CEvaluationNode::type(dynamic_cast<const CEvaluationNode*>((*it)->getChild()->getSibling())->getType()) == CEvaluationNode::NUMBER && dynamic_cast<const CEvaluationNodeNumber*>((*it)->getChild()->getSibling())->value() < 0.0)
            {
              if(fabs(dynamic_cast<const CEvaluationNodeNumber*>((*it)->getChild()->getSibling())->value()) - 1.0 < ZERO)
                {
                  additions.push_back(dynamic_cast<const CEvaluationNode*>((*it)->getChild())->copyBranch());
                  delete *it;
                  it = subtractions.erase(it);
                }
              else
                {
                  CEvaluationNode* pTmp = new CEvaluationNodeOperator(CEvaluationNodeOperator::MULTIPLY, "*");
                  pTmp->addChild(dynamic_cast<const CEvaluationNode*>((*it)->getChild())->copyBranch());
                  std::ostringstream os;
                  os.precision(18);
                  os << static_cast<const CEvaluationNodeNumber*>((*it)->getChild()->getSibling())->value() * -1.0;
                  pTmp->addChild(new CEvaluationNodeNumber(CEvaluationNodeNumber::DOUBLE, os.str().c_str()));
                  additions.push_back(pTmp);
                  delete *it;
                  it = subtractions.erase(it);
                }
            }
          else
            {
              ++it;
            }
        }
      else
        {
          ++it;
        }
    }
}
 */

/**
 * This method expands the exponents of power nodes, e.g. A^(x+y) -> A^x * A^y
 */
CEvaluationNode* CNormalTranslation::expandPowerNodes(const CEvaluationNode* pOrig)
{
  CEvaluationNode* pResult = NULL;
  std::vector<CEvaluationNode*> children;
  const CEvaluationNode* pChild = dynamic_cast<const CEvaluationNode*>(pOrig->getChild());

  while (pChild != NULL)
    {
      children.push_back(CNormalTranslation::expandPowerNodes(pChild));
      pChild = dynamic_cast<const CEvaluationNode*>(pChild->getSibling());
    }

  if (CEvaluationNode::type(pOrig->getType()) == CEvaluationNode::OPERATOR &&
      ((CEvaluationNodeOperator::SubType)CEvaluationNode::subType(pOrig->getType())) == CEvaluationNodeOperator::POWER)
    {
      std::vector<CEvaluationNode*> additions, subtractions;
      CNormalTranslation::splitSum(children[1], additions, subtractions, false);
      CNormalTranslation::swapNegativeNumbers(additions, subtractions);

      // the root node is a fraction
      // the denominator is a product of all subtraction nodes
      // the numerator is a product of all addition nodes
      if (!additions.empty() || !subtractions.empty())
        {
          // replace all nodes in additions and subtractions by
          // children[0]^node so we can use the generic method to create the
          // multiplication chain
          unsigned int i, iMax = additions.size();

          for (i = 0; i < iMax; ++i)
            {
              CEvaluationNode* pTmpNode = NULL;

              if (CEvaluationNode::type(additions[i]->getType()) == CEvaluationNode::NUMBER && fabs(static_cast<const CEvaluationNodeNumber*>(additions[i])->value() - 1.0) < 1e-12)
                {
                  delete additions[i];
                  pTmpNode = children[0]->copyBranch();
                }
              else
                {
                  pTmpNode = new CEvaluationNodeOperator(CEvaluationNodeOperator::POWER, "^");
                  pTmpNode->addChild(children[0]->copyBranch());
                  // don't copy additions element since this has been created in
                  // splitSum
                  pTmpNode->addChild(additions[i]);
                  additions[i] = pTmpNode;
                }

              additions[i] = pTmpNode;
            }

          iMax = subtractions.size();

          for (i = 0; i < iMax; ++i)
            {
              CEvaluationNode* pTmpNode = NULL;

              if (CEvaluationNode::type(subtractions[i]->getType()) == CEvaluationNode::NUMBER && fabs(static_cast<const CEvaluationNodeNumber*>(subtractions[i])->value() - 1.0) < 1e-12)
                {
                  pTmpNode = children[0]->copyBranch();
                  delete subtractions[i];
                }
              else
                {
                  pTmpNode = new CEvaluationNodeOperator(CEvaluationNodeOperator::POWER, "^");
                  pTmpNode->addChild(children[0]->copyBranch());
                  // don't copy subtractions element since this has been created in
                  // splitSum
                  pTmpNode->addChild(subtractions[i]);
                }

              subtractions[i] = pTmpNode;
            }

          // if we have only subtractions, the numerator of the resulting
          // exponent has to be 1
          if (additions.empty())
            {
              pResult = CNormalTranslation::ONE_NODE.copyBranch();
            }
          else
            {
              pResult = CNormalTranslation::createChain(&CNormalTranslation::TIMES_NODE, &CNormalTranslation::ONE_NODE, additions);
              additions.clear();
            }

          assert(pResult != NULL);

          if (!subtractions.empty())
            {
              CEvaluationNode* pTmpNode = new CEvaluationNodeOperator(CEvaluationNodeOperator::DIVIDE, "/");
              pTmpNode->addChild(pResult);
              pResult = pTmpNode;
              pTmpNode = CNormalTranslation::createChain(&CNormalTranslation::TIMES_NODE, &CNormalTranslation::ONE_NODE, subtractions);
              assert(pTmpNode != NULL);
              pResult->addChild(pTmpNode);
              subtractions.clear();
            }

          // delete the children
          delete children[0];
          delete children[1];
        }
    }

  if (pResult == NULL)
    {
      pResult = pOrig->copyNode(children);
    }

  return pResult;
}

/**
 * The methods get a vector of multiplication elements and a vector of division
 * elements and tries to find elements with the same power base in those two vectors.
 */
std::vector<std::pair<CEvaluationNode*, CEvaluationNode*> > CNormalTranslation::matchPowerBases(const std::vector<const CEvaluationNode*>& multiplications, const std::vector<const CEvaluationNode*>& divisions)
{
  std::vector<std::pair<std::pair<const CEvaluationNode*, std::string>, std::vector<CEvaluationNode*> > > matchMap;
  std::vector<const CEvaluationNode*>::const_iterator vit = multiplications.begin(), vendit = multiplications.end();

  while (vit != vendit)
    {
      const CEvaluationNode* pBase = (*vit);
      CEvaluationNode* pExponent = NULL;

      if (CEvaluationNode::type(pBase->getType()) == CEvaluationNode::OPERATOR && ((CEvaluationNodeOperator::SubType)CEvaluationNode::subType(pBase->getType())) == CEvaluationNodeOperator::POWER)
        {
          pBase = dynamic_cast<const CEvaluationNode*>(pBase->getChild());
          pExponent = dynamic_cast<const CEvaluationNode*>(pBase->getSibling())->copyBranch();
        }
      else
        {
          pExponent = new CEvaluationNodeNumber(CEvaluationNodeNumber::DOUBLE, "1.0");
        }

      // check if a base with the same infix is already in the map.
      // if not, add the base
      // if yes, add the exponent to the vector associated with the base
      std::vector<std::pair<std::pair<const CEvaluationNode*, std::string>, std::vector<CEvaluationNode*> > >::iterator mapIt = matchMap.begin(), mapEndit = matchMap.end();
      CNormalFraction* pBase2 = createNormalRepresentation(pBase);
      std::string base2String = pBase2->toString();
      delete pBase2;

      while (mapIt != mapEndit)
        {
          if (mapIt->first.second == base2String)
            {
              mapIt->second.push_back(pExponent);
              break;
            }

          ++mapIt;
        }

      if (mapIt == mapEndit)
        {
          std::vector<CEvaluationNode*> v;
          v.push_back(pExponent);
          matchMap.push_back(std::make_pair(std::pair<const CEvaluationNode*, std::string>(pBase, base2String), v));
        }

      ++vit;
    }

  std::vector<std::pair<std::pair<const CEvaluationNode*, std::string>, std::vector<CEvaluationNode*> > > matchMap2;
  vit = divisions.begin(), vendit = divisions.end();

  while (vit != vendit)
    {
      const CEvaluationNode* pBase = (*vit);
      CEvaluationNode* pExponent = NULL;

      if (CEvaluationNode::type(pBase->getType()) == CEvaluationNode::OPERATOR && ((CEvaluationNodeOperator::SubType)CEvaluationNode::subType(pBase->getType())) == CEvaluationNodeOperator::POWER)
        {
          pBase = dynamic_cast<const CEvaluationNode*>(pBase->getChild());
          pExponent = dynamic_cast<const CEvaluationNode*>(pBase->getSibling())->copyBranch();
        }
      else
        {
          pExponent = new CEvaluationNodeNumber(CEvaluationNodeNumber::DOUBLE, "1.0");
        }

      // check if a base with the same infix is already in the map.
      // if not, add the base
      // if yes, add the exponent to the vector associated with the base
      std::vector<std::pair<std::pair<const CEvaluationNode*, std::string>, std::vector<CEvaluationNode*> > >::iterator mapIt = matchMap2.begin(), mapEndit = matchMap2.end();
      CNormalFraction* pBase2 = createNormalRepresentation(pBase);
      std::string base2String = pBase2->toString();
      delete pBase2;

      while (mapIt != mapEndit)
        {
          if (mapIt->first.second == base2String)
            {
              mapIt->second.push_back(pExponent);
              break;
            }

          ++mapIt;
        }

      if (mapIt == mapEndit)
        {
          std::vector<CEvaluationNode*> v;
          v.push_back(pExponent);
          matchMap2.push_back(std::make_pair(std::pair<const CEvaluationNode*, std::string>(pBase, base2String), v));
        }

      ++vit;
    }

  // now combine the two maps
  std::vector<std::pair<CEvaluationNode*, std::pair<CEvaluationNode*, std::string> > > result;
  std::vector<std::pair<std::pair<const CEvaluationNode*, std::string>, std::vector<CEvaluationNode*> > >::iterator mapIt = matchMap.begin(), mapEndit = matchMap.end();

  while (mapIt != mapEndit)
    {
      CEvaluationNode* pNode = CNormalTranslation::createChain(&CNormalTranslation::PLUS_NODE, &CNormalTranslation::ZERO_NODE, mapIt->second);
      assert(pNode != NULL);
      result.push_back(std::make_pair(pNode, std::pair<CEvaluationNode*, std::string>(mapIt->first.first->copyBranch(), mapIt->first.second)));
      ++mapIt;
    }

  mapIt = matchMap2.begin(), mapEndit = matchMap2.end();

  while (mapIt != mapEndit)
    {
      std::vector<CEvaluationNode*> constVect;
      constVect.insert(constVect.begin(), mapIt->second.begin(), mapIt->second.end());
      CEvaluationNode* pNode = CNormalTranslation::createOperatorChain(CEvaluationNodeOperator::PLUS, "+", constVect);
      // now check if we already have a base with the same infix in the
      // results
      unsigned int i, iMax = result.size();

      for (i = 0; i < iMax; ++i)
        {
          if (result[i].second.second == mapIt->first.second)
            {
              CEvaluationNode* pTmpNode = new CEvaluationNodeOperator(CEvaluationNodeOperator::MINUS, "-");
              pTmpNode->addChild(result[i].first);
              pTmpNode->addChild(pNode);
              result[i] = std::make_pair(pTmpNode, result[i].second);
              break;
            }
        }

      if (i == iMax)
        {
          if (CEvaluationNode::type(pNode->getType()) == CEvaluationNode::NUMBER)
            {
              std::ostringstream os;
              os.precision(18);
              os << static_cast<CEvaluationNodeNumber*>(pNode)->value() * -1.0;
              CEvaluationNode* pTmpNumber = new CEvaluationNodeNumber(CEvaluationNodeNumber::DOUBLE, os.str().c_str());
              delete pNode;
              result.push_back(std::make_pair(pTmpNumber, std::pair<CEvaluationNode*, std::string>(mapIt->first.first->copyBranch(), mapIt->first.second)));
            }
          else
            {
              CEvaluationNode* pTmpNode = new CEvaluationNodeOperator(CEvaluationNodeOperator::MULTIPLY, "*");
              pTmpNode->addChild(new CEvaluationNodeNumber(CEvaluationNodeNumber::DOUBLE, "-1.0"));
              pTmpNode->addChild(pNode);
              result.push_back(std::make_pair(pTmpNode, std::pair<CEvaluationNode*, std::string>(mapIt->first.first->copyBranch(), mapIt->first.second)));
            }
        }

      // delete the obsolete nodes
      iMax = mapIt->second.size();

      for (i = 0; i < iMax; ++i)
        {
          delete mapIt->second[i];
        }

      ++mapIt;
    }

  // copy the result vector into the return data structure
  std::vector<std::pair<CEvaluationNode*, CEvaluationNode*> > tmp;
  unsigned int i, iMax = result.size();
  // since we know how many elements will end up in tmp, we can already reserve
  // the space
  tmp.reserve(iMax);

  for (i = 0; i < iMax; ++i)
    {
      tmp.push_back(std::pair<CEvaluationNode*, CEvaluationNode*>(result[i].first, result[i].second.first));
    }

  return tmp;
}

/**
 * The methods get a vector of addition elements and a vector of subtractions
 * elements and tries to find equal elements in those two vectors.
 */
std::vector<std::pair<CEvaluationNode*, CEvaluationNode*> > CNormalTranslation::matchSummands(const std::vector<CEvaluationNode*>& additions, const std::vector<CEvaluationNode*>& subtractions)
{
  // the individual elements could be  multiplication chains and there could
  // be a common factor somewhere in the chain
  // Since I only want to get rid of numbers, it might be enough to
  // consider only those multiplication chains the contain a number node and
  // something else, everything else is ambiguous anyway and depends on
  // the order of the nodes in the chain
  std::vector<std::pair<std::pair<const CEvaluationNode*, std::string>, std::vector<CEvaluationNode*> > > matchMap;

  std::vector<CEvaluationNode*>::const_iterator vit = additions.begin(), vendit = additions.end();

  while (vit != vendit)
    {
      const CEvaluationNode* pNode = (*vit);
      CEvaluationNode* pFactor = NULL;

      if (CEvaluationNode::type(pNode->getType()) == CEvaluationNode::OPERATOR && ((CEvaluationNodeOperator::SubType)CEvaluationNode::subType(pNode->getType())) == CEvaluationNodeOperator::MULTIPLY)
        {
          const CEvaluationNode* pChild1 = dynamic_cast<const CEvaluationNode*>(pNode->getChild());
          assert(pChild1 != NULL);
          const CEvaluationNode* pChild2 = dynamic_cast<const CEvaluationNode*>(pChild1->getSibling());
          assert(pChild2 != NULL);
          assert(pChild2->getSibling() == NULL);

          if (CEvaluationNode::type(pChild1->getType()) == CEvaluationNode::NUMBER)
            {
              pNode = pChild2;
              pFactor = pChild1->copyBranch();
            }
          else if (CEvaluationNode::type(pChild2->getType()) == CEvaluationNode::NUMBER)
            {
              pNode = pChild1;
              pFactor = pChild2->copyBranch();
            }
          else
            {
              pFactor = new CEvaluationNodeNumber(CEvaluationNodeNumber::DOUBLE, "1.0");
            }
        }
      else
        {
          pFactor = new CEvaluationNodeNumber(CEvaluationNodeNumber::DOUBLE, "1.0");
        }

      // check if a node with the same infix is already in the map.
      // if not, add the base
      // if yes, add the exponent to the vector associated with the base
      CNormalFraction* pBase2 = createNormalRepresentation(pNode);
      std::string base2String = pBase2->toString();
      std::vector<std::pair<std::pair<const CEvaluationNode*, std::string>, std::vector<CEvaluationNode*> > >::iterator mapIt = matchMap.begin(), mapEndit = matchMap.end();

      while (mapIt != mapEndit)
        {
          if (mapIt->first.second == base2String)
            {
              mapIt->second.push_back(pFactor);
              break;
            }

          ++mapIt;
        }

      delete pBase2;

      if (mapIt == mapEndit)
        {
          std::vector<CEvaluationNode*> v;
          v.push_back(pFactor);
          matchMap.push_back(std::make_pair(std::pair<const CEvaluationNode*, std::string>(pNode, base2String), v));
        }

      ++vit;
    }

  std::vector<std::pair<std::pair<const CEvaluationNode*, std::string>, std::vector<CEvaluationNode*> > > matchMap2;
  vit = subtractions.begin(), vendit = subtractions.end();

  while (vit != vendit)
    {
      const CEvaluationNode* pNode = (*vit);
      CEvaluationNode* pFactor = NULL;

      if (CEvaluationNode::type(pNode->getType()) == CEvaluationNode::OPERATOR && ((CEvaluationNodeOperator::SubType)CEvaluationNode::subType(pNode->getType())) == CEvaluationNodeOperator::MULTIPLY)
        {
          const CEvaluationNode* pChild1 = dynamic_cast<const CEvaluationNode*>(pNode->getChild());
          assert(pChild1 != NULL);
          const CEvaluationNode* pChild2 = dynamic_cast<const CEvaluationNode*>(pChild1->getSibling());
          assert(pChild2 != NULL);
          assert(pChild2->getSibling() == NULL);

          if (CEvaluationNode::type(pChild1->getType()) == CEvaluationNode::NUMBER)
            {
              pNode = pChild2;
              pFactor = pChild1->copyBranch();
            }
          else if (CEvaluationNode::type(pChild2->getType()) == CEvaluationNode::NUMBER)
            {
              pNode = pChild1;
              pFactor = pChild2->copyBranch();
            }
          else
            {
              pFactor = new CEvaluationNodeNumber(CEvaluationNodeNumber::DOUBLE, "1.0");
            }
        }
      else
        {
          pFactor = new CEvaluationNodeNumber(CEvaluationNodeNumber::DOUBLE, "1.0");
        }

      // check if a node with the same infix is already in the map.
      // if not, add the node
      // if yes, add the 1 to the vector associated with the base
      std::vector<std::pair<std::pair<const CEvaluationNode*, std::string>, std::vector<CEvaluationNode*> > >::iterator mapIt = matchMap2.begin(), mapEndit = matchMap2.end();
      CNormalFraction* pBase2 = createNormalRepresentation(pNode);
      std::string base2String = pBase2->toString();

      while (mapIt != mapEndit)
        {
          if (mapIt->first.second == base2String)
            {
              mapIt->second.push_back(pFactor);
              break;
            }

          ++mapIt;
        }

      delete pBase2;

      if (mapIt == mapEndit)
        {
          std::vector<CEvaluationNode*> v;
          v.push_back(pFactor);
          matchMap2.push_back(std::make_pair(std::pair<const CEvaluationNode*, std::string>(pNode, base2String), v));
        }

      ++vit;
    }

  // now combine the two maps
  std::vector<std::pair<CEvaluationNode*, std::pair<CEvaluationNode*, std::string> > > result;
  std::vector<std::pair<std::pair<const CEvaluationNode*, std::string>, std::vector<CEvaluationNode*> > >::iterator mapIt = matchMap.begin(), mapEndit = matchMap.end();

  while (mapIt != mapEndit)
    {
      CEvaluationNode* pNode = CNormalTranslation::createChain(&CNormalTranslation::PLUS_NODE, &CNormalTranslation::ZERO_NODE, mapIt->second);
      assert(pNode != NULL);
      result.push_back(std::make_pair(pNode, std::pair<CEvaluationNode*, std::string>(mapIt->first.first->copyBranch(), mapIt->first.second)));
      ++mapIt;
    }

  mapIt = matchMap2.begin(), mapEndit = matchMap2.end();

  while (mapIt != mapEndit)
    {
      std::vector<CEvaluationNode*> constVect;
      constVect.insert(constVect.begin(), mapIt->second.begin(), mapIt->second.end());
      CEvaluationNode* pNode = CNormalTranslation::createOperatorChain(CEvaluationNodeOperator::PLUS, "+", constVect);
      // now check if we already have a base with the same infix in the
      // results
      unsigned int i, iMax = result.size();

      for (i = 0; i < iMax; ++i)
        {
          if (result[i].second.second == mapIt->first.second)
            {
              CEvaluationNode* pTmpNode = new CEvaluationNodeOperator(CEvaluationNodeOperator::MINUS, "-");
              pTmpNode->addChild(result[i].first);
              pTmpNode->addChild(pNode);
              result[i] = std::make_pair(pTmpNode, result[i].second);
              break;
            }
        }

      if (i == iMax)
        {
          CEvaluationNode* pTmpNode = NULL;

          if (CEvaluationNode::type(pNode->getType()) == CEvaluationNode::NUMBER)
            {
              std::ostringstream os;
              os.precision(18);
              os << dynamic_cast<const CEvaluationNodeNumber*>(pNode)->value() * -1.0;
              pTmpNode = new CEvaluationNodeNumber(CEvaluationNodeNumber::DOUBLE, os.str().c_str());
              delete pNode;
            }
          else
            {
              pTmpNode = new CEvaluationNodeOperator(CEvaluationNodeOperator::MULTIPLY, "*");
              pTmpNode->addChild(new CEvaluationNodeNumber(CEvaluationNodeNumber::DOUBLE, "-1.0"));
              pTmpNode->addChild(pNode);
            }

          result.push_back(std::make_pair(pTmpNode, std::pair<CEvaluationNode*, std::string>(mapIt->first.first->copyBranch(), mapIt->first.second)));
        }

      // delete the obsolete nodes
      iMax = mapIt->second.size();

      for (i = 0; i < iMax; ++i)
        {
          delete mapIt->second[i];
        }

      ++mapIt;
    }

  std::vector<std::pair<CEvaluationNode*, CEvaluationNode*> > tmp;
  // copy the result vector into the expected return data type
  unsigned int i, iMax = result.size();
  // since we know how many item will end up in the vector we can already
  // reserve the space
  tmp.reserve(iMax);

  for (i = 0; i < iMax; ++i)
    {
      tmp.push_back(std::pair<CEvaluationNode*, CEvaluationNode*>(result[i].first, result[i].second.first));
    }

  return tmp;
}

/**
 * This method expands products. (A+B)*(C+D) -> (A*C)+(A*D)+(B*C)+(B*D)
CEvaluationNode* CNormalTranslation::expandProducts(const CEvaluationNode* pOrig)
{
  CEvaluationNode* pResult = NULL;
  // we have to create operation chains and do the mutliplication
  // on the numerator and the denominator chain if the node is a multiplication
  // or a division
  if (CEvaluationNode::type(pOrig->getType()) == CEvaluationNode::OPERATOR && ((CEvaluationNodeOperator::SubType)CEvaluationNode::subType(pOrig->getType()) == CEvaluationNodeOperator::MULTIPLY || (CEvaluationNodeOperator::SubType)CEvaluationNode::subType(pOrig->getType()) == CEvaluationNodeOperator::DIVIDE))
    {
      std::vector<const CEvaluationNode*> multiplications, divisions;
      CNormalTranslation::splitProduct(pOrig, multiplications, divisions, false);
      unsigned int i, iMax = multiplications.size();
      CEvaluationNode* pTmpResult;
      for (i = 0;i < iMax;++i)
        {
          if (pResult == NULL)
            {
              pResult = CNormalTranslation::expandProducts(multiplications[i]);
            }
          else
            {
              CEvaluationNode* pTmpNode = CNormalTranslation::expandProducts(multiplications[i]);
              pTmpResult = CNormalTranslation::multiply(pResult, pTmpNode);
              delete pResult;
              delete pTmpNode;
              pResult = pTmpResult;
            }
        }
      if (!divisions.empty())
        {
          CEvaluationNode* pDenominator = NULL;
          iMax = divisions.size();
          for (i = 0;i < iMax;++i)
            {
              if (pDenominator == NULL)
                {
                  pDenominator = CNormalTranslation::expandProducts(divisions[i]);
                }
              else
                {
                  CEvaluationNode* pTmpNode = CNormalTranslation::expandProducts(divisions[i]);
                  pTmpResult = CNormalTranslation::multiply(pDenominator, pTmpNode);
                  delete pDenominator;
                  delete pTmpNode;
                  pDenominator = pTmpResult;
                }
              //delete divisions[i];
            }
          pTmpResult = new CEvaluationNodeOperator(CEvaluationNodeOperator::DIVIDE, "/");
          pTmpResult->addChild(pResult);
          pTmpResult->addChild(pDenominator);
          pResult = pTmpResult;
        }
    }
  else
    {
      const CEvaluationNode* pChild = dynamic_cast<const CEvaluationNode*>(pOrig->getChild());
      std::vector<CEvaluationNode*> children;
      while (pChild != NULL)
        {
          children.push_back(CNormalTranslation::expandProducts(pChild));
          pChild = dynamic_cast<const CEvaluationNode*>(pChild->getSibling());
        }
      if (CEvaluationNode::type(pOrig->getType()) == CEvaluationNode::OPERATOR &&
          ((CEvaluationNodeOperator::SubType)CEvaluationNode::subType(pOrig->getType())) == CEvaluationNodeOperator::MULTIPLY)
        {
          assert(children.size() == 2);
          pResult = CNormalTranslation::multiply(children[0], children[1]);
          // delete the children
          delete children[0];
          delete children[1];
        }
      if (pResult == NULL)
        {
          pResult = pOrig->copyNode(children);
        }
    }
  return pResult;
}
 */

/**
 * This method expands products. (A+B)*(C+D) -> (A*C)+(A*D)+(B*C)+(B*D)
 */
CEvaluationNode* CNormalTranslation::expandProducts(CEvaluationNode* pOrig)
{
  // this is done depth first
  CEvaluationNode* pResult = pOrig;
  CEvaluationNode* pChild = dynamic_cast<CEvaluationNode*>(pOrig->getChild());
  CEvaluationNode* pLastChild = pOrig;

  while (pChild != NULL)
    {
      CEvaluationNode* pNewChild = CNormalTranslation::expandProducts(pChild);
      assert(pNewChild != NULL);

      if (pNewChild != pChild)
        {
          // remove the old child and add the new one
          pOrig->removeChild(pChild);
          delete pChild;
          pChild = pNewChild;
          pOrig->addChild(pNewChild, pLastChild);
        }

      pLastChild = pChild;
      pChild = dynamic_cast<CEvaluationNode*>(pChild->getSibling());
    }

  // we have to create operation chains and do the multiplication
  // on the numerator and the denominator chain if the node is a multiplication
  // or a division
  if (CEvaluationNode::type(pOrig->getType()) == CEvaluationNode::OPERATOR && ((CEvaluationNodeOperator::SubType)CEvaluationNode::subType(pOrig->getType()) == CEvaluationNodeOperator::MULTIPLY || (CEvaluationNodeOperator::SubType)CEvaluationNode::subType(pOrig->getType()) == CEvaluationNodeOperator::DIVIDE))
    {
      std::vector<const CEvaluationNode*> multiplications, divisions;
      CNormalTranslation::splitProduct(pOrig, multiplications, divisions, false);
      unsigned int i, iMax = multiplications.size();
      CEvaluationNode* pTmpResult;

      for (i = 0; i < iMax; ++i)
        {
          if (pResult == pOrig)
            {
              pResult = multiplications[i]->copyBranch();
              assert(pResult != NULL);
            }
          else
            {
              pTmpResult = CNormalTranslation::multiply(pResult, multiplications[i]);
              delete pResult;
              pResult = pTmpResult;
              assert(pResult != NULL);
            }
        }

      if (!divisions.empty())
        {
          CEvaluationNode* pDenominator = NULL;
          iMax = divisions.size();

          for (i = 0; i < iMax; ++i)
            {
              if (pDenominator == NULL)
                {
                  pDenominator = divisions[i]->copyBranch();
                  assert(pDenominator != NULL);
                }
              else
                {
                  pTmpResult = CNormalTranslation::multiply(pDenominator, divisions[i]);
                  delete pDenominator;
                  pDenominator = pTmpResult;
                  assert(pDenominator != NULL);
                }
            }

          pTmpResult = new CEvaluationNodeOperator(CEvaluationNodeOperator::DIVIDE, "/");
          pTmpResult->addChild(pResult);
          pTmpResult->addChild(pDenominator);
          pResult = pTmpResult;
        }
    }

  return pResult;
}

/**
 * Multiplies the two given nodes and returns the result.
 */
CEvaluationNode* CNormalTranslation::multiply(const CEvaluationNode* pNode1, const CEvaluationNode* pNode2)
{
  CEvaluationNode* pResult = NULL;
  std::vector<const CEvaluationNode*> additions1, subtractions1;
  CNormalTranslation::splitSum(pNode1, additions1, subtractions1, false);
  std::vector<const CEvaluationNode*> additions2, subtractions2;
  CNormalTranslation::splitSum(pNode2, additions2, subtractions2, false);
  // multiply every element in additions1 with every element in additions2
  // and subtractions2 the results for the multiplication with the elements
  // of subtractions2 must be multiplied by -1
  // multiply every element in subtraction1 with every element in additions2
  // and subtractions2 the results for the multiplication with the elements
  // of additions2 must be multiplied by -1
  std::vector<CEvaluationNode*> tmp;
  unsigned int i, iMax = additions1.size();

  for (i = 0; i < iMax; ++i)
    {
      unsigned int j, jMax = additions2.size();

      for (j = 0; j < jMax; ++j)
        {
          CEvaluationNode* pMult = new CEvaluationNodeOperator(CEvaluationNodeOperator::MULTIPLY, "*");
          pMult->addChild(additions1[i]->copyBranch());
          pMult->addChild(additions2[j]->copyBranch());
          tmp.push_back(pMult);
        }
    }

  iMax = subtractions1.size();

  for (i = 0; i < iMax; ++i)
    {
      unsigned int j, jMax = subtractions2.size();

      for (j = 0; j < jMax; ++j)
        {
          CEvaluationNode* pMult = new CEvaluationNodeOperator(CEvaluationNodeOperator::MULTIPLY, "*");
          pMult->addChild(subtractions1[i]->copyBranch());
          pMult->addChild(subtractions2[j]->copyBranch());
          tmp.push_back(pMult);
        }
    }

  if (!tmp.empty())
    {
      pResult = CNormalTranslation::createChain(&CNormalTranslation::PLUS_NODE, &CNormalTranslation::ZERO_NODE, tmp);
      assert(pResult != NULL);
      tmp.clear();
    }

  iMax = additions1.size();

  for (i = 0; i < iMax; ++i)
    {
      unsigned int j, jMax = subtractions2.size();

      for (j = 0; j < jMax; ++j)
        {
          CEvaluationNode* pMult = new CEvaluationNodeOperator(CEvaluationNodeOperator::MULTIPLY, "*");
          pMult->addChild(additions1[i]->copyBranch());
          pMult->addChild(subtractions2[j]->copyBranch());
          tmp.push_back(pMult);
        }
    }

  iMax = subtractions1.size();

  for (i = 0; i < iMax; ++i)
    {
      unsigned int j, jMax = additions2.size();

      for (j = 0; j < jMax; ++j)
        {
          CEvaluationNode* pMult = new CEvaluationNodeOperator(CEvaluationNodeOperator::MULTIPLY, "*");
          pMult->addChild(subtractions1[i]->copyBranch());
          pMult->addChild(additions2[j]->copyBranch());
          tmp.push_back(pMult);
        }
    }

  if (!tmp.empty())
    {
      if (pResult != NULL)
        {
          CEvaluationNode* pTmpNode = new CEvaluationNodeOperator(CEvaluationNodeOperator::MINUS, "-");
          pTmpNode->addChild(pResult);
          pResult = pTmpNode;
          pTmpNode = CNormalTranslation::createChain(&CNormalTranslation::PLUS_NODE, &CNormalTranslation::ZERO_NODE, tmp);
          assert(pTmpNode != NULL);
          pResult->addChild(pTmpNode);
        }
      else
        {
          if (tmp.size() == 1 && CEvaluationNode::type(tmp[0]->getType()) == CEvaluationNode::NUMBER)
            {
              std::ostringstream os;
              os.precision(18);
              os << tmp[0]->value() * -1.0;
              pResult = new CEvaluationNodeNumber(CEvaluationNodeNumber::DOUBLE, os.str().c_str());
              delete tmp[0];
            }
          else
            {
              CEvaluationNode* pTmpNode = new CEvaluationNodeOperator(CEvaluationNodeOperator::MULTIPLY, "*");
              pTmpNode->addChild(new CEvaluationNodeNumber(CEvaluationNodeNumber::DOUBLE, "-1.0"));
              pResult = pTmpNode;
              pTmpNode = CNormalTranslation::createChain(&CNormalTranslation::PLUS_NODE, &CNormalTranslation::ZERO_NODE, tmp);
              assert(pTmpNode != NULL);
              pResult->addChild(pTmpNode);
            }
        }
    }

  return pResult;
}

/**
 * This method does all the canceling on a given node and its children.
 */
CEvaluationNode* CNormalTranslation::cancel(const CEvaluationNode* pOrig)
{
  // TODO I think this method has much potential for improvement
  // TODO since the comparison code seems to spend about 85% of the time
  // TODO here, this is where I should start making optimizations
  //
  // try to find multiplication chains where something is divided by itself
  // or multiplied by -1 times itself
  // also consider powers (it's the bases that have to match)
  //
  // try to find addition changes where there is a subtraction of two
  // identical nodes or an addition of one node and the same node times -1
  CEvaluationNode* pResult = NULL;
  std::vector<CEvaluationNode*> children;

  if (CEvaluationNode::type(pOrig->getType()) == CEvaluationNode::OPERATOR)
    {
      if (((CEvaluationNodeOperator::SubType)CEvaluationNode::subType(pOrig->getType())) == CEvaluationNodeOperator::PLUS ||
          ((CEvaluationNodeOperator::SubType)CEvaluationNode::subType(pOrig->getType())) == CEvaluationNodeOperator::MINUS)
        {
          // we are in a sum
          std::vector<CEvaluationNode*> additions, subtractions;
          CNormalTranslation::splitSum(pOrig, additions, subtractions, false);
          CNormalTranslation::swapNegativeNumbers(additions, subtractions);
          // collect all nodes in additions and subtractions
          unsigned int i, iMax = additions.size();

          for (i = 0; i < iMax; ++i)
            {
              CEvaluationNode* pChild = CNormalTranslation::cancel(additions[i]);
              delete additions[i];
              additions[i] = pChild;
            }

          iMax = subtractions.size();

          for (i = 0; i < iMax; ++i)
            {
              CEvaluationNode* pChild = CNormalTranslation::cancel(subtractions[i]);
              delete subtractions[i];
              subtractions[i] = pChild;
            }

          // find identical nodes in additions and subtractions
          // The first entry in the pair is the collected factor
          // the second entry is the original branch
          // make sure the collected factor is again simplified
          std::vector<std::pair<CEvaluationNode*, CEvaluationNode*> > collected = CNormalTranslation::matchSummands(additions, subtractions);
          iMax = additions.size();

          for (i = 0; i < iMax; ++i)
            {
              delete additions[i];
            }

          additions.clear();
          iMax = subtractions.size();

          for (i = 0; i < iMax; ++i)
            {
              delete subtractions[i];
            }

          subtractions.clear();
          std::vector<CEvaluationNode*> chain;
          iMax = collected.size();

          for (i = 0; i < iMax; ++i)
            {
              std::pair<CEvaluationNode*, CEvaluationNode*> pair = collected[i];

              //CEvaluationNode* pTmpNode = CNormalTranslation::eliminate(pair.first);
              //delete pair.first;
              // if simplified node is 0.0, we ignore this node
              if (CEvaluationNode::type(pair.first->getType()) == CEvaluationNode::NUMBER &&
                  fabs(dynamic_cast<CEvaluationNodeNumber*>(pair.first)->value()) < ZERO)
                {
                  delete pair.first;
                  delete pair.second;
                }
              else if (CEvaluationNode::type(pair.first->getType()) == CEvaluationNode::NUMBER &&
                       fabs(dynamic_cast<CEvaluationNodeNumber*>(pair.first)->value() - 1.0) < ZERO)
                {
                  delete pair.first;
                  chain.push_back(pair.second);
                }
              else
                {
                  CEvaluationNode* pMult = new CEvaluationNodeOperator(CEvaluationNodeOperator::MULTIPLY, "*");
                  pMult->addChild(pair.first);
                  pMult->addChild(pair.second);
                  chain.push_back(pMult);
                }
            }

          pResult = CNormalTranslation::createChain(&CNormalTranslation::PLUS_NODE, &CNormalTranslation::ZERO_NODE, chain);
          assert(pResult != NULL);
        }
      else if (((CEvaluationNodeOperator::SubType)CEvaluationNode::subType(pOrig->getType())) == CEvaluationNodeOperator::MULTIPLY ||
               ((CEvaluationNodeOperator::SubType)CEvaluationNode::subType(pOrig->getType())) == CEvaluationNodeOperator::DIVIDE)
        {
          // we are in a product
          std::vector<const CEvaluationNode*> multiplications, divisions;
          CNormalTranslation::splitProduct(pOrig, multiplications, divisions, false);
          // collect all nodes in multiplications and divisions
          unsigned int i, iMax = multiplications.size();

          for (i = 0; i < iMax; ++i)
            {
              multiplications[i] = CNormalTranslation::cancel(multiplications[i]);
            }

          iMax = divisions.size();

          for (i = 0; i < iMax; ++i)
            {
              divisions[i] = CNormalTranslation::cancel(divisions[i]);
            }

          // find identical nodes in multiplications and divisions
          // The first entry in the pair is the collected power exponent
          // the second entry is the original power base
          // make sure the collected factor is again simplified
          std::vector<std::pair<CEvaluationNode*, CEvaluationNode*> > collected = CNormalTranslation::matchPowerBases(multiplications, divisions);
          iMax = multiplications.size();

          for (i = 0; i < iMax; ++i)
            {
              delete multiplications[i];
            }

          multiplications.clear();
          iMax = divisions.size();

          for (i = 0; i < iMax; ++i)
            {
              delete divisions[i];
            }

          divisions.clear();
          std::vector<CEvaluationNode*> numeratorChain;
          std::vector<CEvaluationNode*> denominatorChain;
          iMax = collected.size();

          for (i = 0; i < iMax; ++i)
            {
              std::pair<CEvaluationNode*, CEvaluationNode*> pair = collected[i];

              //CEvaluationNode* pTmpNode = CNormalTranslation::eliminate(pair.first);
              //delete pair.first;
              // if simplified node is a 0.0, we ignore this node
              if (CEvaluationNode::type(pair.first->getType()) == CEvaluationNode::NUMBER)
                {
                  if (fabs(dynamic_cast<CEvaluationNodeNumber*>(pair.first)->value()) < ZERO)
                    {
                      delete pair.first;
                      delete pair.second;
                    }
                  else if (dynamic_cast<CEvaluationNodeNumber*>(pair.first)->value() > 0.0)
                    {
                      if (fabs(dynamic_cast<CEvaluationNodeNumber*>(pair.first)->value() - 1.0) < ZERO)
                        {
                          delete pair.first;
                          numeratorChain.push_back(pair.second);
                        }
                      else
                        {
                          CEvaluationNode* pPower = new CEvaluationNodeOperator(CEvaluationNodeOperator::POWER, "^");
                          pPower->addChild(pair.second);
                          pPower->addChild(pair.first);
                          numeratorChain.push_back(pPower);
                        }
                    }
                  else
                    {
                      if (fabs(dynamic_cast<CEvaluationNodeNumber*>(pair.first)->value() + 1.0) < ZERO)
                        {
                          delete pair.first;
                          denominatorChain.push_back(pair.second);
                        }
                      else
                        {
                          CEvaluationNode* pPower = new CEvaluationNodeOperator(CEvaluationNodeOperator::POWER, "^");
                          pPower->addChild(pair.second);
                          std::ostringstream os;
                          os.precision(18);
                          os << fabs(dynamic_cast<const CEvaluationNodeNumber*>(pair.first)->value());
                          pPower->addChild(new CEvaluationNodeNumber((CEvaluationNodeNumber::SubType)CEvaluationNode::subType(pair.first->getType()), os.str().c_str()));
                          delete pair.first;
                          denominatorChain.push_back(pPower);
                        }
                    }
                }
              else
                {
                  CEvaluationNode* pPower = new CEvaluationNodeOperator(CEvaluationNodeOperator::POWER, "^");

                  // check if the node is -1.0 * SOMETHING
                  if (CEvaluationNode::type(pair.first->getType()) == CEvaluationNode::OPERATOR && (CEvaluationNodeOperator::SubType)CEvaluationNode::subType(pair.first->getType()) == CEvaluationNodeOperator::MULTIPLY
                      && CEvaluationNode::type(dynamic_cast<const CEvaluationNode*>(pair.first->getChild())->getType()) == CEvaluationNode::NUMBER)
                    {
                      if (fabs(static_cast<const CEvaluationNodeNumber*>(pair.first->getChild())->value() + 1.0) < ZERO)
                        {
                          pPower->addChild(pair.second);
                          pPower->addChild(dynamic_cast<const CEvaluationNode*>(pair.first->getChild()->getSibling())->copyBranch());
                          delete pair.first;
                          denominatorChain.push_back(pPower);
                        }
                      else if (fabs(static_cast<const CEvaluationNodeNumber*>(pair.first->getChild())->value()) < ZERO)
                        {
                          // delete the power node and add
                          delete pPower;
                          delete pair.first;
                          numeratorChain.push_back(pair.second);
                        }
                      else if (static_cast<const CEvaluationNodeNumber*>(pair.first->getChild())->value() < 0.0)
                        {
                          pPower->addChild(pair.second);
                          pPower->addChild(pair.first);
                          denominatorChain.push_back(pPower);
                        }
                      else
                        {
                          pPower->addChild(pair.second);
                          pPower->addChild(pair.first);
                          numeratorChain.push_back(pPower);
                        }
                    }
                  else
                    {
                      pPower->addChild(pair.second);
                      pPower->addChild(pair.first);
                      numeratorChain.push_back(pPower);
                    }
                }
            }

          // if there are only divisions, we have an empty numerator chain
          if (numeratorChain.empty())
            {
              pResult = CNormalTranslation::ONE_NODE.copyBranch();
            }
          else
            {
              pResult = CNormalTranslation::createChain(&CNormalTranslation::TIMES_NODE, &CNormalTranslation::ONE_NODE, numeratorChain);
              assert(pResult != NULL);
            }

          if (!denominatorChain.empty())
            {
              CEvaluationNode* pTmpNode = new CEvaluationNodeOperator(CEvaluationNodeOperator::DIVIDE, "/");
              pTmpNode->addChild(pResult);
              pResult = pTmpNode;
              pTmpNode = CNormalTranslation::createChain(&CNormalTranslation::TIMES_NODE, &CNormalTranslation::ONE_NODE, denominatorChain);
              assert(pTmpNode != NULL);
              pResult->addChild(pTmpNode);
            }
        }
      else
        {
          const CEvaluationNode* pChild = dynamic_cast<const CEvaluationNode*>(pOrig->getChild());

          while (pChild != NULL)
            {
              children.push_back(CNormalTranslation::cancel(pChild));
              pChild = dynamic_cast<const CEvaluationNode*>(pChild->getSibling());
            }

          assert(children.size() == 2);
        }
    }
  else
    {
      const CEvaluationNode* pChild = dynamic_cast<const CEvaluationNode*>(pOrig->getChild());

      while (pChild != NULL)
        {
          children.push_back(CNormalTranslation::cancel(pChild));
          pChild = dynamic_cast<const CEvaluationNode*>(pChild->getSibling());
        }
    }

  if (pResult == NULL)
    {
      pResult = pOrig->copyNode(children);
    }

  return pResult;
}

/**
 * This method eliminates directly nested fractions.
 */
CEvaluationNode* CNormalTranslation::eliminateDirectlyNestedFractions(const CEvaluationNode* pOrig)
{
  if (pOrig == NULL) return NULL;

  CEvaluationNode* pResult = NULL;

  if (CEvaluationNode::type(pOrig->getType()) == CEvaluationNode::OPERATOR && (CEvaluationNodeOperator::SubType)CEvaluationNode::subType(pOrig->getType()) == CEvaluationNodeOperator::DIVIDE)
    {
      std::vector<CEvaluationNode*> children;
      const CEvaluationNode* pChild = dynamic_cast<const CEvaluationNode*>(pOrig->getChild());

      while (pChild != NULL)
        {
          children.push_back(CNormalTranslation::eliminateDirectlyNestedFractions(pChild));
          pChild = dynamic_cast<const CEvaluationNode*>(pChild->getSibling());
        }

      // check if one of the children (or both) are a division
      assert(children.size() == 2);

      if (CEvaluationNode::type(children[0]->getType()) == CEvaluationNode::OPERATOR && (CEvaluationNodeOperator::SubType)CEvaluationNode::subType(children[0]->getType()) == CEvaluationNodeOperator::DIVIDE)
        {
          if (CEvaluationNode::type(children[1]->getType()) == CEvaluationNode::OPERATOR && (CEvaluationNodeOperator::SubType)CEvaluationNode::subType(children[1]->getType()) == CEvaluationNodeOperator::DIVIDE)
            {
              // both children are division
              pResult = new CEvaluationNodeOperator(CEvaluationNodeOperator::DIVIDE, "/");
              CEvaluationNode* pTmp = new CEvaluationNodeOperator(CEvaluationNodeOperator::MULTIPLY, "*");
              pTmp->addChild(dynamic_cast<const CEvaluationNode*>(children[0]->getChild())->copyBranch());
              pTmp->addChild(dynamic_cast<const CEvaluationNode*>(children[1]->getChild()->getSibling())->copyBranch());
              pResult->addChild(pTmp);
              pTmp = new CEvaluationNodeOperator(CEvaluationNodeOperator::MULTIPLY, "*");
              pTmp->addChild(dynamic_cast<const CEvaluationNode*>(children[0]->getChild()->getSibling())->copyBranch());
              pTmp->addChild(dynamic_cast<const CEvaluationNode*>(children[1]->getChild())->copyBranch());
              pResult->addChild(pTmp);
              delete children[0];
              delete children[1];
            }
          else
            {
              // only the first child is a division
              pResult = new CEvaluationNodeOperator(CEvaluationNodeOperator::DIVIDE, "/");
              pResult->addChild(dynamic_cast<const CEvaluationNode*>(children[0]->getChild())->copyBranch());
              CEvaluationNode* pTmp = new CEvaluationNodeOperator(CEvaluationNodeOperator::MULTIPLY, "*");
              pTmp->addChild(dynamic_cast<const CEvaluationNode*>(children[0]->getChild()->getSibling())->copyBranch());
              pTmp->addChild(children[1]);
              pResult->addChild(pTmp);
              delete children[0];
            }
        }
      else if (CEvaluationNode::type(children[1]->getType()) == CEvaluationNode::OPERATOR && (CEvaluationNodeOperator::SubType)CEvaluationNode::subType(children[1]->getType()) == CEvaluationNodeOperator::DIVIDE)
        {
          // only the second child is a division
          pResult = new CEvaluationNodeOperator(CEvaluationNodeOperator::DIVIDE, "/");
          CEvaluationNode* pTmp = new CEvaluationNodeOperator(CEvaluationNodeOperator::MULTIPLY, "*");
          pTmp->addChild(children[0]);
          pTmp->addChild(dynamic_cast<const CEvaluationNode*>(children[1]->getChild()->getSibling())->copyBranch());
          pResult->addChild(pTmp);
          pResult->addChild(dynamic_cast<const CEvaluationNode*>(children[1]->getChild())->copyBranch());
          delete children[1];
        }
      else
        {
          pResult = pOrig->copyNode(children);
        }
    }
  else
    {
      std::vector<CEvaluationNode*> children;
      const CEvaluationNode* pChild = dynamic_cast<const CEvaluationNode*>(pOrig->getChild());

      while (pChild != NULL)
        {
          children.push_back(CNormalTranslation::eliminateDirectlyNestedFractions(pChild));
          pChild = dynamic_cast<const CEvaluationNode*>(pChild->getSibling());
        }

      pResult = pOrig->copyNode(children);
    }

  return pResult;
}

CEvaluationNode* CNormalTranslation::eliminatePowersOfFractions(const CEvaluationNode* pOrig)
{
  if (pOrig == NULL) return NULL;

  CEvaluationNode* pResult = NULL;

  if (CEvaluationNode::type(pOrig->getType()) == CEvaluationNode::OPERATOR && (CEvaluationNodeOperator::SubType)CEvaluationNode::subType(pOrig->getType()) == CEvaluationNodeOperator::POWER)
    {
      std::vector<CEvaluationNode*> children;
      const CEvaluationNode* pChild = dynamic_cast<const CEvaluationNode*>(pOrig->getChild());

      while (pChild != NULL)
        {
          children.push_back(CNormalTranslation::eliminatePowersOfFractions(pChild));
          pChild = dynamic_cast<const CEvaluationNode*>(pChild->getSibling());
        }

      // check if the first child is a fraction
      assert(children.size() == 2);

      if (CEvaluationNode::type(children[0]->getType()) == CEvaluationNode::OPERATOR && (CEvaluationNodeOperator::SubType)CEvaluationNode::subType(children[0]->getType()) == CEvaluationNodeOperator::DIVIDE)
        {
          // the first child is a division
          pResult = new CEvaluationNodeOperator(CEvaluationNodeOperator::DIVIDE, "/");

          CEvaluationNode* pTmp = new CEvaluationNodeOperator(CEvaluationNodeOperator::POWER, "^");
          pTmp->addChild(dynamic_cast<const CEvaluationNode*>(children[0]->getChild())->copyBranch());
          pTmp->addChild(children[1]->copyBranch());
          pResult->addChild(pTmp);
          pTmp = new CEvaluationNodeOperator(CEvaluationNodeOperator::POWER, "^");
          pTmp->addChild(dynamic_cast<const CEvaluationNode*>(children[0]->getChild()->getSibling())->copyBranch());
          pTmp->addChild(children[1]);
          pResult->addChild(pTmp);
          delete children[0];
        }
      else
        {
          pResult = pOrig->copyNode(children);
        }
    }
  else
    {
      std::vector<CEvaluationNode*> children;
      const CEvaluationNode* pChild = dynamic_cast<const CEvaluationNode*>(pOrig->getChild());

      while (pChild != NULL)
        {
          children.push_back(CNormalTranslation::eliminatePowersOfFractions(pChild));
          pChild = dynamic_cast<const CEvaluationNode*>(pChild->getSibling());
        }

      pResult = pOrig->copyNode(children);
    }

  return pResult;
}
// TODO for comparing in infix, it should be brought into a normalform first

/**
 * This methods converts a product of fractions into a fraction of products.
 */
CEvaluationNode* CNormalTranslation::product2fraction(const CEvaluationNode* pOrig)
{
  CEvaluationNode* pResult = NULL;
  std::vector<CEvaluationNode*> children;
  const CEvaluationNode* pChild = dynamic_cast<const CEvaluationNode*>(pOrig->getChild());

  while (pChild != NULL)
    {
      children.push_back(CNormalTranslation::product2fraction(pChild));
      pChild = dynamic_cast<const CEvaluationNode*>(pChild->getSibling());
    }

  if (CEvaluationNode::type(pOrig->getType()) == CEvaluationNode::OPERATOR && (CEvaluationNodeOperator::SubType)CEvaluationNode::subType(pOrig->getType()) == CEvaluationNodeOperator::MULTIPLY)
    {
      CEvaluationNode* pNumerator1 = NULL;
      CEvaluationNode* pNumerator2 = NULL;
      CEvaluationNode* pDenominator1 = NULL;
      CEvaluationNode* pDenominator2 = NULL;
      assert(children.size() == 2);

      if (CEvaluationNode::type(children[0]->getType()) == CEvaluationNode::OPERATOR && (CEvaluationNodeOperator::SubType)CEvaluationNode::subType(children[0]->getType()) == CEvaluationNodeOperator::DIVIDE)
        {
          pNumerator1 = dynamic_cast<CEvaluationNode*>(children[0]->getChild());
          pDenominator1 = dynamic_cast<CEvaluationNode*>(children[0]->getChild()->getSibling());
        }
      else
        {
          pNumerator1 = children[0];
        }

      if (CEvaluationNode::type(children[1]->getType()) == CEvaluationNode::OPERATOR && (CEvaluationNodeOperator::SubType)CEvaluationNode::subType(children[1]->getType()) == CEvaluationNodeOperator::DIVIDE)
        {
          pNumerator2 = dynamic_cast<CEvaluationNode*>(children[1]->getChild());
          pDenominator2 = dynamic_cast<CEvaluationNode*>(children[1]->getChild()->getSibling());
        }
      else
        {
          pNumerator2 = children[1];
        }

      if (pDenominator1 || pDenominator2)
        {
          pResult = new CEvaluationNodeOperator(CEvaluationNodeOperator::DIVIDE, "/");
          CEvaluationNodeOperator* pTmpNode = new CEvaluationNodeOperator(CEvaluationNodeOperator::MULTIPLY, "*");
          pTmpNode->addChild(pNumerator1->copyBranch());
          pTmpNode->addChild(pNumerator2->copyBranch());
          pResult->addChild(pTmpNode);

          if (pDenominator1 && pDenominator2)
            {
              pTmpNode = new CEvaluationNodeOperator(CEvaluationNodeOperator::MULTIPLY, "*");
              pTmpNode->addChild(pDenominator1->copyBranch());
              pTmpNode->addChild(pDenominator2->copyBranch());
              pResult->addChild(pTmpNode);
            }
          else if (pDenominator1)
            {
              pResult->addChild(pDenominator1->copyBranch());
            }
          else
            {
              pResult->addChild(pDenominator2->copyBranch());
            }

          delete children[0];
          delete children[1];
        }
      else
        {
          pResult = pOrig->copyNode(children);
        }
    }
  else
    {
      pResult = pOrig->copyNode(children);
    }

  return pResult;
}

/**
 * Given a root node, this method traverses the tree and expands produtcs in
 * power bases to multiplications of power items.
 * It is the responsibility of the caller to delete the returned node.
 */
CEvaluationNode* CNormalTranslation::expandPowerBases(const CEvaluationNode* pRoot)
{
  CEvaluationNode* pResult = NULL;
  CEvaluationNode::Type type = pRoot->getType();

  if (CEvaluationNode::type(type) == CEvaluationNode::OPERATOR && (CEvaluationNodeOperator::SubType)CEvaluationNode::subType(type) == CEvaluationNodeOperator::POWER)
    {
      const CEvaluationNode* pBase = dynamic_cast<const CEvaluationNode*>(pRoot->getChild());
      const CEvaluationNode* pExp = dynamic_cast<const CEvaluationNode*>(pBase->getSibling());
      type = pBase->getType();

      if (CEvaluationNode::type(type) == CEvaluationNode::OPERATOR && ((CEvaluationNodeOperator::SubType)CEvaluationNode::subType(type) == CEvaluationNodeOperator::MULTIPLY || (CEvaluationNodeOperator::SubType)CEvaluationNode::subType(type) == CEvaluationNodeOperator::DIVIDE))
        {
          std::vector<const CEvaluationNode*> multiplications, divisions;
          std::vector<CEvaluationNode*> numeratorNodes, denominatorNodes;
          CNormalTranslation::splitProduct(pBase, multiplications, divisions, false);
          std::vector<const CEvaluationNode*>::const_iterator it = multiplications.begin(), endit = multiplications.end();
          CEvaluationNode* pPower = NULL;

          while (it != endit)
            {
              pPower = new CEvaluationNodeOperator(CEvaluationNodeOperator::POWER, "^");
              pPower->addChild(CNormalTranslation::expandPowerBases(*it));
              pPower->addChild(pExp->copyBranch());
              numeratorNodes.push_back(pPower);
              ++it;
            }

          pResult = CNormalTranslation::createChain(&CNormalTranslation::TIMES_NODE, &CNormalTranslation::ONE_NODE, numeratorNodes);
          assert(pResult != NULL);

          if (!divisions.empty())
            {
              it = divisions.begin(), endit = divisions.end();

              while (it != endit)
                {
                  pPower = new CEvaluationNodeOperator(CEvaluationNodeOperator::POWER, "^");
                  pPower->addChild(CNormalTranslation::expandPowerBases(*it));
                  pPower->addChild(pExp->copyBranch());
                  denominatorNodes.push_back(pPower);
                  ++it;
                }

              CEvaluationNode* pTmpResult = new CEvaluationNodeOperator(CEvaluationNodeOperator::DIVIDE, "/");
              pTmpResult->addChild(pResult);
              pResult = CNormalTranslation::createChain(&CNormalTranslation::TIMES_NODE, &CNormalTranslation::ONE_NODE, denominatorNodes);
              assert(pResult != NULL);
              pTmpResult->addChild(pResult);
              pResult = pTmpResult;
            }
        }
      else if (CEvaluationNode::type(type) == CEvaluationNode::OPERATOR && ((CEvaluationNodeOperator::SubType)CEvaluationNode::subType(type) == CEvaluationNodeOperator::PLUS || (CEvaluationNodeOperator::SubType)CEvaluationNode::subType(type) == CEvaluationNodeOperator::MINUS))
        {
          std::vector<CEvaluationNode*> additions, subtractions;
          CNormalTranslation::splitSum(pBase, additions, subtractions, false);
          CNormalTranslation::swapNegativeNumbers(additions, subtractions);
          std::pair<CEvaluationNode*, CEvaluationNode*> resultPair = CNormalTranslation::factorize(additions, subtractions);
          unsigned int i, iMax = additions.size();

          for (i = 0; i < iMax; ++i)
            {
              delete additions[i];
            }

          additions.clear();
          iMax = subtractions.size();

          for (i = 0; i < iMax; ++i)
            {
              delete subtractions[i];
            }

          subtractions.clear();

          if (resultPair.first != NULL)
            {
              pResult = new CEvaluationNodeOperator(CEvaluationNodeOperator::MULTIPLY, "*");
              CEvaluationNode* pTmpNode = new CEvaluationNodeOperator(CEvaluationNodeOperator::POWER, "^");
              pTmpNode->addChild(resultPair.first);
              pTmpNode->addChild(pExp->copyBranch());
              pResult->addChild(pTmpNode);
              pTmpNode = new CEvaluationNodeOperator(CEvaluationNodeOperator::POWER, "^");
              pTmpNode->addChild(resultPair.second);
              pTmpNode->addChild(pExp->copyBranch());
              pResult->addChild(pTmpNode);
            }
          else
            {
              // there are no common factors
              std::vector<CEvaluationNode*> children;
              const CEvaluationNode* pChild = dynamic_cast<const CEvaluationNode*>(pRoot->getChild());

              while (pChild != NULL)
                {
                  children.push_back(CNormalTranslation::expandPowerBases(pChild));
                  pChild = dynamic_cast<const CEvaluationNode*>(pChild->getSibling());
                }

              pResult = pRoot->copyNode(children);
            }
        }
      else
        {
          std::vector<CEvaluationNode*> children;
          const CEvaluationNode* pChild = dynamic_cast<const CEvaluationNode*>(pRoot->getChild());

          while (pChild != NULL)
            {
              children.push_back(CNormalTranslation::expandPowerBases(pChild));
              pChild = dynamic_cast<const CEvaluationNode*>(pChild->getSibling());
            }

          pResult = pRoot->copyNode(children);
        }
    }
  else
    {
      std::vector<CEvaluationNode*> children;
      const CEvaluationNode* pChild = dynamic_cast<const CEvaluationNode*>(pRoot->getChild());

      while (pChild != NULL)
        {
          children.push_back(CNormalTranslation::expandPowerBases(pChild));
          pChild = dynamic_cast<const CEvaluationNode*>(pChild->getSibling());
        }

      pResult = pRoot->copyNode(children);
    }

  return pResult;
}

/**
 * This method takes two vectors and checks if the elements in the two vectors
 * can be split into multiplications and divisions and if there a common factors in all resulting subgroups.
 */
std::pair<CEvaluationNode*, CEvaluationNode*> CNormalTranslation::factorize(const std::vector<CEvaluationNode*>& additions, const std::vector<CEvaluationNode*>& subtractions)
{
  std::vector<const CEvaluationNode*> commonMultiplications;
  std::vector<const CEvaluationNode*> commonDivisions;
  // additions must have at least one entry
  assert(additions.size() > 0);
  // get all multipllications and divisions from the first entry in additions
  std::vector<const CEvaluationNode*> multiplications, divisions;
  unsigned int i, iMax = additions.size();
  unsigned int iiMax = iMax + subtractions.size();
  std::vector<std::vector<const CEvaluationNode*> > multiplicationVectors, divisionVectors;

  for (i = 0; i < iiMax; ++i)
    {
      const CEvaluationNode* pTmpNode = (i < iMax) ? additions[i] : subtractions[i - iMax];
      CEvaluationNode::Type type = pTmpNode->getType();

      if (CEvaluationNode::type(type) == CEvaluationNode::OPERATOR)
        {
          CEvaluationNodeOperator::SubType subType = (CEvaluationNodeOperator::SubType)CEvaluationNode::subType(type);

          if (subType == CEvaluationNodeOperator::MULTIPLY || subType == CEvaluationNodeOperator::DIVIDE)
            {
              CNormalTranslation::splitProduct(pTmpNode, multiplications, divisions, false);
            }
          else
            {
              multiplications.push_back(pTmpNode);
            }
        }
      else
        {
          multiplications.push_back(pTmpNode);
        }

      multiplicationVectors.push_back(multiplications);
      divisionVectors.push_back(divisions);
      multiplications.clear();
      divisions.clear();
    }

  // now we first search for common multiplications
  multiplications = multiplicationVectors[0];
  std::vector<const CEvaluationNode*>::const_iterator it = multiplications.begin(), endit = multiplications.end();

  while (it != endit)
    {
      bool everywhere = true;
      std::vector<std::vector<const CEvaluationNode*> >::iterator innerIt = multiplicationVectors.begin(), innerEndit = multiplicationVectors.end();
      // we can leave out the first one since the item comes from there anyway,
      // so we know it is in that vector
      std::string infix = (*it)->getInfix();
      ++innerIt;

      while (innerIt != innerEndit)
        {
          bool found = false;
          std::vector<const CEvaluationNode*>::iterator innerIt2 = (*innerIt).begin(), innerEndit2 = (*innerIt).end();

          while (innerIt2 != innerEndit2)
            {
              if ((*innerIt2)->getInfix() == infix)
                {
                  found = true;
                  break;
                }

              ++innerIt2;
            }

          if (!found)
            {
              everywhere = false;
              break;
            }

          ++innerIt;
        }

      // if the item was found as a factor in all other additions and
      // subtractions, we know it is a common factor, we add it to the
      // commonFactors and update the additions and subtractions
      if (everywhere)
        {
          commonMultiplications.push_back(*it);
          std::vector<std::vector<const CEvaluationNode*> >::iterator innerIt = multiplicationVectors.begin();
          std::vector<std::vector<const CEvaluationNode*> >::iterator innerEndit = multiplicationVectors.end();

          while (innerIt != innerEndit)
            {
              std::vector<const CEvaluationNode*>::iterator innerIt2 = (*innerIt).begin();
              std::vector<const CEvaluationNode*>::iterator innerEndit2 = (*innerIt).end();

              while (innerIt2 != innerEndit2)
                {
                  if ((*innerIt2)->getInfix() == infix)
                    {
                      innerIt->erase(innerIt2);
                      break;
                    }

                  ++innerIt2;
                }

              ++innerIt;
            }
        }

      ++it;
    }

  // now we search for common divisions
  divisions = divisionVectors[0];

  if (!divisions.empty())
    {
      it = divisions.begin(), endit = divisions.end();

      while (it != endit)
        {
          bool everywhere = true;
          std::vector<std::vector<const CEvaluationNode*> >::iterator innerIt = divisionVectors.begin(), innerEndit = divisionVectors.end();
          // we can leav out the first one since the item comes from there anyway,
          // so we know it is in that vector
          std::string infix = (*it)->getInfix();
          ++innerIt;

          while (innerIt != innerEndit)
            {
              bool found = false;
              std::vector<const CEvaluationNode*>::iterator innerIt2 = (*innerIt).begin(), innerEndit2 = (*innerIt).end();

              while (innerIt2 != innerEndit2)
                {
                  if ((*innerIt2)->getInfix() == infix)
                    {
                      found = true;
                      break;
                    }

                  ++innerIt2;
                }

              if (!found)
                {
                  everywhere = false;
                  break;
                }

              ++innerIt;
            }

          // if the item was found as a factor in all other additions and
          // subtractions, we know it is a common factor, we add it to the
          // commonFactors and update the additions and subtractions
          if (everywhere)
            {
              commonDivisions.push_back(*it);
              innerIt = divisionVectors.begin();
              innerEndit = divisionVectors.end();

              while (innerIt != innerEndit)
                {
                  std::vector<const CEvaluationNode*>::iterator innerIt2 = (*innerIt).begin();
                  std::vector<const CEvaluationNode*>::iterator innerEndit2 = (*innerIt).end();

                  while (innerIt2 != innerEndit2)
                    {
                      if ((*innerIt2)->getInfix() == infix)
                        {
                          innerIt->erase(innerIt2);
                          break;
                        }

                      ++innerIt2;
                    }

                  ++innerIt;
                }
            }

          ++it;
        }
    }

  // create the two resulting nodes
  // first we have to create new additions and subtraction vectors which we
  // then combine into a subtraction
  // then we combine all commonMultiplications and commonDivisions into a
  // division
  // those two nodes are then returned in a pair
  CEvaluationNode* pFirstNode = NULL;
  CEvaluationNode* pSecondNode = NULL;

  if (!(commonMultiplications.empty() && commonDivisions.empty()))
    {
      unsigned int i, iMax = additions.size();
      unsigned int iiMax = iMax + subtractions.size();
      std::vector<CEvaluationNode*> newAdditions, newSubtractions;

      for (i = 0; i < iiMax; ++i)
        {
          // since the createOperatorChain automatically returns 1 as the result
          // if an empty vector is given, we don't have to worry about havinf
          // removed all items from the vectors above.
          CEvaluationNode* pTmpNode = CNormalTranslation::createOperatorChain(CEvaluationNodeOperator::MULTIPLY, "*", multiplicationVectors[i]);

          if (!divisionVectors[i].empty())
            {
              CEvaluationNode* pTmpNode2 = new CEvaluationNodeOperator(CEvaluationNodeOperator::DIVIDE, "/");
              pTmpNode2->addChild(pTmpNode);
              pTmpNode2->addChild(CNormalTranslation::createOperatorChain(CEvaluationNodeOperator::MULTIPLY, "*", divisionVectors[i]));
              pTmpNode = pTmpNode2;
            }

          if (i < iMax)
            {
              newAdditions.push_back(pTmpNode);
            }
          else
            {
              newSubtractions.push_back(pTmpNode);
            }
        }

      pSecondNode = CNormalTranslation::createChain(&CNormalTranslation::PLUS_NODE, &CNormalTranslation::ZERO_NODE, newAdditions);
      assert(pSecondNode != NULL);

      if (!newSubtractions.empty())
        {
          CEvaluationNode* pTmpNode = new CEvaluationNodeOperator(CEvaluationNodeOperator::MINUS, "-");
          pTmpNode->addChild(pSecondNode);
          pTmpNode->addChild(CNormalTranslation::createChain(&CNormalTranslation::PLUS_NODE, &CNormalTranslation::ZERO_NODE, newSubtractions));
          assert(pTmpNode != NULL);
          pSecondNode = pTmpNode;
        }

      pFirstNode = CNormalTranslation::createOperatorChain(CEvaluationNodeOperator::MULTIPLY, "*", commonMultiplications);

      if (!commonDivisions.empty())
        {
          CEvaluationNode* pTmpNode = new CEvaluationNodeOperator(CEvaluationNodeOperator::DIVIDE, "/");
          pTmpNode->addChild(pFirstNode);
          pTmpNode->addChild(CNormalTranslation::createOperatorChain(CEvaluationNodeOperator::MULTIPLY, "*", commonDivisions));
          pFirstNode = pTmpNode;
        }
    }

  return std::pair<CEvaluationNode*, CEvaluationNode*>(pFirstNode, pSecondNode);
}

/**
 * Given a vector of nodes, this method creates a multiplication chain of
 * all the nodes. The chain contains the original nodes and not copies.
 */
CEvaluationNode* CNormalTranslation::createOperatorChain(CEvaluationNodeOperator::SubType type, const char* data, const std::vector<CEvaluationNode*>& nodes)
{
  std::vector<const CEvaluationNode*> tmpV;
  std::vector<CEvaluationNode*>::const_iterator it = nodes.begin(), endit = nodes.end();

  while (it != endit)
    {
      tmpV.push_back(*it);
      ++it;
    }

  return CNormalTranslation::createOperatorChain(type, data, tmpV);
}

/**
 * This version of create chain copies the given elements in the vector and
 *then calls the createChain method which does not need copying.
 */
CEvaluationNode* CNormalTranslation::createChain(const CEvaluationNode* pLink, const CEvaluationNode* pNeutralElement, const std::vector<const CEvaluationNode*>& elements)
{
  std::vector<CEvaluationNode*> tmpVector;
  tmpVector.reserve(elements.size());
  std::vector<const CEvaluationNode*>::const_iterator it = elements.begin(), endit = elements.end();

  while (it != endit)
    {
      tmpVector.push_back((*it)->copyBranch());
      ++it;
    }

  return CNormalTranslation::createChain(pLink, pNeutralElement, tmpVector);
}

/**
 * This method creates a chain of operations. The individual elements are
 * linked with copies of pLink.
 * NULL is returned if elements is empty.
 * So if this method is used to create a chanin of OR linked elements which
 * will be embedded in another and linked chain, the neutral element should be
 * a TRUE node since AND combining something with true does not change the result.
 * The neutral element is the element that does not change the result of the
 * operation represented be pLink. So if pLink represents a multiplication,
 * the neutral element is the number node 1.0.
 * This method does not copy the elements in the given vector, but uses them in
 *the chain directly.
 */
CEvaluationNode* CNormalTranslation::createChain(const CEvaluationNode* pLink, const CEvaluationNode* /*pNeutralElement*/, const std::vector<CEvaluationNode*>& elements)
{
  CEvaluationNode* pResult = NULL;

  if (elements.size() == 1)
    {
      pResult = elements[0];
    }
  else if (elements.size() > 1)
    {
      std::vector<CEvaluationNode*>::const_reverse_iterator it = elements.rbegin(), endit = elements.rend();
      CEvaluationNode* pOperator = pLink->copyBranch();
      CEvaluationNode* pChild = *it;
      ++it;
      pOperator->addChild(*it);
      pOperator->addChild(pChild);
      ++it;
      pChild = pOperator;

      while (it != endit)
        {
          pOperator = pLink->copyBranch();
          pOperator->addChild(*it);
          pOperator->addChild(pChild);
          pChild = pOperator;
          ++it;
        }

      pResult = pOperator;
    }

  return pResult;
}

/**
 * This routine moves all negative numbers from vector v1 to v2
 * and changes the number to a positive number.
 */
void CNormalTranslation::swapNegativeNumbers(std::vector<CEvaluationNode*>& v1, std::vector<CEvaluationNode*>& v2)
{
  CEvaluationNode* pNode;
  CEvaluationNode::Type type, type1, type2;
  std::vector<CEvaluationNode*>::iterator it = v1.begin(), endit = v1.end();
  std::ostringstream os;

  while (it != endit)
    {
      pNode = *it;
      type = pNode->getType();

      if (CEvaluationNode::type(type) == CEvaluationNode::NUMBER && pNode->value() < 0.0)
        {
          it = v1.erase(it);
          endit = v1.end();
          os.str("");
          os.precision(18);
          os << pNode->value() * -1.0;
          v2.push_back(new CEvaluationNodeNumber(CEvaluationNodeNumber::DOUBLE, os.str().c_str()));
          delete pNode;
          continue;
        }
      else if (CEvaluationNode::type(type) == CEvaluationNode::OPERATOR && (CEvaluationNodeOperator::SubType)CEvaluationNode::subType(type) == CEvaluationNodeOperator::MULTIPLY)
        {
          CEvaluationNode* pChild1 = dynamic_cast<CEvaluationNode*>(pNode->getChild());
          assert(pChild1 != NULL);
          CEvaluationNode* pChild2 = dynamic_cast<CEvaluationNode*>(pChild1->getSibling());
          assert(pChild2 != NULL);
          type1 = CEvaluationNode::type(pChild1->getType());
          type2 = CEvaluationNode::type(pChild2->getType());

          if (type1 == CEvaluationNode::NUMBER || type2 == CEvaluationNode::NUMBER)
            {
              if (type1 == CEvaluationNode::NUMBER && type2 == CEvaluationNode::NUMBER)
                {
                  if (pChild1->value() < 0.0 && pChild2->value() < 0.0)
                    {
                      os.str("");
                      os.precision(18);
                      os << pChild1->value() * -1.0;
                      CEvaluationNodeNumber* pTmpNode = new CEvaluationNodeNumber(CEvaluationNodeNumber::DOUBLE, os.str().c_str());
                      pNode->removeChild(pChild1);
                      pNode->removeChild(pChild2);
                      pNode->addChild(pTmpNode);
                      delete pChild1;
                      os.str("");
                      os.precision(18);
                      os << pChild2->value() * -1.0;
                      pTmpNode = new CEvaluationNodeNumber(CEvaluationNodeNumber::DOUBLE, os.str().c_str());
                      pNode->addChild(pTmpNode);
                      delete pChild2;
                    }
                  else if (pChild1->value() < 0.0)
                    {
                      os.str("");
                      os.precision(18);
                      os << pChild1->value() * -1.0;
                      CEvaluationNodeNumber* pTmpNode = new CEvaluationNodeNumber(CEvaluationNodeNumber::DOUBLE, os.str().c_str());
                      pNode->removeChild(pChild1);
                      pNode->addChild(pTmpNode, pNode);
                      delete pChild1;
                      v2.push_back(pNode);
                      it = v1.erase(it);
                      endit = v1.end();
                      continue;
                    }
                  else if (pChild2->value() < 0.0)
                    {
                      os.str("");
                      os.precision(18);
                      os << pChild2->value() * -1.0;
                      CEvaluationNodeNumber* pTmpNode = new CEvaluationNodeNumber(CEvaluationNodeNumber::DOUBLE, os.str().c_str());
                      pNode->removeChild(pChild2);
                      pNode->addChild(pTmpNode, pChild1);
                      delete pChild2;
                      v2.push_back(pNode);
                      it = v1.erase(it);
                      endit = v1.end();
                      continue;
                    }
                }
              else if (type1 == CEvaluationNode::NUMBER)
                {
                  if (pChild1->value() < 0.0)
                    {
                      os.str("");
                      os.precision(18);
                      os << pChild1->value() * -1.0;
                      CEvaluationNodeNumber* pTmpNode = new CEvaluationNodeNumber(CEvaluationNodeNumber::DOUBLE, os.str().c_str());
                      pNode->removeChild(pChild1);
                      pNode->addChild(pTmpNode, pNode);
                      delete pChild1;
                      v2.push_back(pNode);
                      it = v1.erase(it);
                      endit = v1.end();
                      continue;
                    }
                }
              else
                {
                  if (pChild2->value() < 0.0)
                    {
                      os.str("");
                      os.precision(18);
                      os << pChild2->value() * -1.0;
                      CEvaluationNodeNumber* pTmpNode = new CEvaluationNodeNumber(CEvaluationNodeNumber::DOUBLE, os.str().c_str());
                      pNode->removeChild(pChild2);
                      pNode->addChild(pTmpNode, pChild1);
                      delete pChild2;
                      v2.push_back(pNode);
                      it = v1.erase(it);
                      endit = v1.end();
                      continue;
                    }
                }
            }
        }

      ++it;
    }

  it = v2.begin();
  endit = v2.end();

  while (it != endit)
    {
      pNode = *it;
      type = pNode->getType();

      if (CEvaluationNode::type(type) == CEvaluationNode::NUMBER && pNode->value() < 0.0)
        {
          it = v2.erase(it);
          endit = v2.end();
          os.str("");
          os.precision(18);
          os << pNode->value() * -1.0;
          v1.push_back(new CEvaluationNodeNumber(CEvaluationNodeNumber::DOUBLE, os.str().c_str()));
          delete pNode;
          continue;
        }
      else if (CEvaluationNode::type(type) == CEvaluationNode::OPERATOR && (CEvaluationNodeOperator::SubType)CEvaluationNode::subType(type) == CEvaluationNodeOperator::MULTIPLY)
        {
          CEvaluationNode* pChild1 = dynamic_cast<CEvaluationNode*>(pNode->getChild());
          assert(pChild1 != NULL);
          CEvaluationNode* pChild2 = dynamic_cast<CEvaluationNode*>(pChild1->getSibling());
          assert(pChild2 != NULL);
          type1 = CEvaluationNode::type(pChild1->getType());
          type2 = CEvaluationNode::type(pChild2->getType());

          if (type1 == CEvaluationNode::NUMBER || type2 == CEvaluationNode::NUMBER)
            {
              if (type1 == CEvaluationNode::NUMBER && type2 == CEvaluationNode::NUMBER)
                {
                  if (pChild1->value() < 0.0 && pChild2->value() < 0.0)
                    {
                      os.str("");
                      os.precision(18);
                      os << pChild1->value() * -1.0;
                      CEvaluationNodeNumber* pTmpNode = new CEvaluationNodeNumber(CEvaluationNodeNumber::DOUBLE, os.str().c_str());
                      pNode->removeChild(pChild1);
                      pNode->removeChild(pChild2);
                      pNode->addChild(pTmpNode);
                      delete pChild1;
                      os.str("");
                      os.precision(18);
                      os << pChild2->value() * -1.0;
                      pTmpNode = new CEvaluationNodeNumber(CEvaluationNodeNumber::DOUBLE, os.str().c_str());
                      pNode->addChild(pTmpNode);
                      delete pChild2;
                    }
                  else if (pChild1->value() < 0.0)
                    {
                      os.str("");
                      os.precision(18);
                      os << pChild1->value() * -1.0;
                      CEvaluationNodeNumber* pTmpNode = new CEvaluationNodeNumber(CEvaluationNodeNumber::DOUBLE, os.str().c_str());
                      pNode->removeChild(pChild1);
                      pNode->addChild(pTmpNode, pNode);
                      delete pChild1;
                      v1.push_back(pNode);
                      it = v2.erase(it);
                      endit = v2.end();
                      continue;
                    }
                  else if (pChild2->value() < 0.0)
                    {
                      os.str("");
                      os.precision(18);
                      os << pChild2->value() * -1.0;
                      CEvaluationNodeNumber* pTmpNode = new CEvaluationNodeNumber(CEvaluationNodeNumber::DOUBLE, os.str().c_str());
                      pNode->removeChild(pChild2);
                      pNode->addChild(pTmpNode, pChild1);
                      delete pChild2;
                      v1.push_back(pNode);
                      it = v2.erase(it);
                      endit = v2.end();
                      continue;
                    }
                }
              else if (type1 == CEvaluationNode::NUMBER)
                {
                  if (pChild1->value() < 0.0)
                    {
                      os.str("");
                      os.precision(18);
                      os << pChild1->value() * -1.0;
                      CEvaluationNodeNumber* pTmpNode = new CEvaluationNodeNumber(CEvaluationNodeNumber::DOUBLE, os.str().c_str());
                      pNode->removeChild(pChild1);
                      pNode->addChild(pTmpNode, pNode);
                      delete pChild1;
                      v1.push_back(pNode);
                      it = v2.erase(it);
                      endit = v2.end();
                      continue;
                    }
                }
              else
                {
                  if (pChild2->value() < 0.0)
                    {
                      os.str("");
                      os.precision(18);
                      os << pChild2->value() * -1.0;
                      CEvaluationNodeNumber* pTmpNode = new CEvaluationNodeNumber(CEvaluationNodeNumber::DOUBLE, os.str().c_str());
                      pNode->removeChild(pChild2);
                      pNode->addChild(pTmpNode, pChild1);
                      delete pChild2;
                      v1.push_back(pNode);
                      it = v2.erase(it);
                      endit = v2.end();
                      continue;
                    }
                }
            }
        }

      ++it;
    }
}
