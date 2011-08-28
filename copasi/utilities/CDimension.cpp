/* Begin CVS Header
 $Source: /fs/turing/cvs/copasi_dev/copasi/utilities/CDimension.cpp,v $
 $Revision: 1.9 $
 $Name: Build-33 $
 $Author: ssahle $
 $Date: 2009/05/08 22:38:18 $
 End CVS Header */

// Copyright (C) 2008 by Pedro Mendes, Virginia Tech Intellectual
// Properties, Inc., EML Research, gGmbH, University of Heidelberg,
// and The University of Manchester.
// All rights reserved.

// Copyright (C) 2001 - 2007 by Pedro Mendes, Virginia Tech Intellectual
// Properties, Inc. and EML Research, gGmbH.
// All rights reserved.

#include <sstream>
#include "CDimension.h"
#include "CopasiDataModel/CCopasiDataModel.h"
#include "report/CCopasiRootContainer.h"
#include "model/CModel.h"
#include "model/CChemEq.h"

CDimension::CDimension()
    : mD1(0), mD2(0), mD3(0), mD4(0), mD5(0),
    mUnknown(true),
    mContradiction(false)
{}

void CDimension::setUnknown()
{
  mUnknown = true;
  mContradiction = false;
}
bool CDimension::isUnknown() const
{
  return mUnknown;
}

void CDimension::setContradiction()
{
  mUnknown = false;
  mContradiction = true;
}
bool CDimension::isContradiction() const
{
  return mContradiction;
}

void CDimension::setDimension(const C_FLOAT64 & d1, const C_FLOAT64 & d2, const C_FLOAT64 & d3,
                              const C_FLOAT64 & d4, const C_FLOAT64 & d5)
{
  mUnknown = false;
  mContradiction = false;
  mD1 = d1;
  mD2 = d2;
  mD3 = d3;
  mD4 = d4;
  mD5 = d5;
}

//static
std::string CDimension::constructDisplayElement(const std::string & base, C_FLOAT64 exponent)
{
  if (exponent <= 0) return "";

  if (exponent == 1.0) return base;

  std::ostringstream ss;
  ss << base << "^" << exponent;
  return ss.str();
}

std::string CDimension::getDisplayString(const CCopasiDataModel* pDataModel) const
{
  if (isUnknown())
    return "?";

  if (isContradiction())
    return "EEE";

  assert(pDataModel != NULL);
  std::string vol = /*FROM_UTF8*/(pDataModel->getModel()->getVolumeUnitName());
  std::string time = /*FROM_UTF8*/(pDataModel->getModel()->getTimeUnitName());
  std::string quan = /*FROM_UTF8*/(pDataModel->getModel()->getQuantityUnitName());
  std::string area = /*FROM_UTF8*/(pDataModel->getModel()->getAreaUnitName());
  std::string len = /*FROM_UTF8*/(pDataModel->getModel()->getLengthUnitName());

  std::string tmp;

  //positive exponents
  std::string s1;
  s1 = constructDisplayElement(quan, mD1);

  tmp = constructDisplayElement(vol, mD2);

  if ((s1 != "") && (tmp != "")) s1 += "*";

  s1 += tmp;

  tmp = constructDisplayElement(time, mD3);

  if ((s1 != "") && (tmp != "")) s1 += "*";

  s1 += tmp;

  tmp = constructDisplayElement(area, mD4);

  if ((s1 != "") && (tmp != "")) s1 += "*";

  s1 += tmp;

  tmp = constructDisplayElement(len, mD5);

  if ((s1 != "") && (tmp != "")) s1 += "*";

  s1 += tmp;

  //negative exponents
  std::string s2;
  bool parflag = false;
  s2 = constructDisplayElement(quan, -mD1);

  tmp = constructDisplayElement(vol, -mD2);

  if ((s2 != "") && (tmp != ""))
    {s2 += "*"; parflag = true;}

  s2 += tmp;

  tmp = constructDisplayElement(time, -mD3);

  if ((s2 != "") && (tmp != ""))
    {s2 += "*"; parflag = true;}

  s2 += tmp;

  tmp = constructDisplayElement(area, -mD4);

  if ((s2 != "") && (tmp != ""))
    {s2 += "*"; parflag = true;}

  s2 += tmp;

  tmp = constructDisplayElement(len, -mD5);

  if ((s2 != "") && (tmp != ""))
    {s2 += "*"; parflag = true;}

  s2 += tmp;

  if (parflag) s2 = "(" + s2 + ")";

  //put it together..
  if ((s1 == "") && (s2 == ""))
    return "1";

  if (s1 == "")
    return "1/" + s2;

  if (s2 == "")
    return s1;

  return s1 + "/" + s2;
}

bool CDimension::operator==(const CDimension & rhs) const
{
  return (mUnknown == rhs.mUnknown)
         && (mContradiction == rhs.mContradiction)
         && (mD1 == rhs.mD1)
         && (mD2 == rhs.mD2)
         && (mD3 == rhs.mD3)
         && (mD4 == rhs.mD4)
         && (mD5 == rhs.mD5);
}

CDimension CDimension::operator+(const CDimension & rhs) const
{
  CDimension result;

  if (isContradiction() || rhs.isContradiction())
    result.setContradiction();
  else if (isUnknown() || rhs.isUnknown())
    result.setUnknown();
  else
    result.setDimension(mD1 + rhs.mD1, mD2 + rhs.mD2, mD3 + rhs.mD3, mD4 + rhs.mD4, mD5 + rhs.mD5);

  return result;
}

CDimension CDimension::operator-(const CDimension & rhs) const
{
  CDimension result;

  if (isContradiction() || rhs.isContradiction())
    result.setContradiction();
  else if (isUnknown() || rhs.isUnknown())
    result.setUnknown();
  else
    result.setDimension(mD1 - rhs.mD1, mD2 - rhs.mD2, mD3 - rhs.mD3, mD4 - rhs.mD4, mD5 - rhs.mD5);

  return result;
}

CDimension CDimension::operator*(const C_FLOAT64 & rhs) const
{
  CDimension result;

  if (isContradiction())
    result.setContradiction();
  else if (isUnknown())
    result.setUnknown();
  else
    result.setDimension(mD1 * rhs, mD2 * rhs, mD3 * rhs, mD4 * rhs, mD5 * rhs);

  return result;
}

CDimension CDimension::compare(const CDimension & rhs) const
{
  CDimension result;

  if (this->isContradiction() || rhs.isContradiction())
    result.setContradiction();
  else if (*this == rhs)
    result = *this;
  else if (this->isUnknown())
    result = rhs;
  else if (rhs.isUnknown())
    result = *this;
  else
    result.setContradiction();

  return result;
}

/**
 * Disabled becuase the in order to generate output, the dimensions instance
 * needs the datamodel.
std::ostream & operator<<(std::ostream &os, const CDimension & d)
{
  if (d.mUnknown) os << "Dim: unknown";
  else if (d.mContradiction) os << "Dim: conctradiction";
  else os << "Dim: (" << d.mD1 << ", " << d.mD2 << ", " << d.mD3 << ")  " << d.getDisplayString();

  return os;
}
 */

std::string CDimension::print(const CCopasiDataModel* pDataModel) const
{
  std::ostringstream os;

  if (this->mUnknown) os << "Dim: unknown";
  else if (this->mContradiction) os << "Dim: contradiction";
  else os << "Dim: (" << this->mD1 << ", " << this->mD2 << ", " << this->mD3
    << ", " << this->mD4 << ", " << this->mD5 << ")  "
    << this->getDisplayString(pDataModel);

  return os.str();;
}

void CDimension::fixDimensionless(bool d1, bool d2, bool d3, bool d4, bool d5)
{
  if (d1)
    mD1 = 0;

  if (d2)
    mD2 = 0;

  if (d3)
    mD3 = 0;

  if (d4)
    mD4 = 0;

  if (d5)
    mD5 = 0;
}

//*************************************************************************

#include "function/CFunction.h"

CFindDimensions::CFindDimensions(const CFunction* function, bool d1, bool d2, bool d3, bool d4, bool d5)
    : mpFunction(function),
    mRootDimension(),
    mUseHeuristics(false),
    mM1(-1.0), mM2(-1.0),
    mD1(d1), mD2(d2), mD3(d3), mD4(d4), mD5(d5)
{
  setupDimensions();
}

void CFindDimensions::setupDimensions()
{
  if (!mpFunction) return;

  mDimensions.resize(mpFunction->getVariables().size());

  unsigned C_INT32 i, imax = mpFunction->getVariables().size();

  for (i = 0; i < imax; ++i)
    {
      switch (mpFunction->getVariables()[i]->getUsage())
        {
          case CFunctionParameter::SUBSTRATE:
          case CFunctionParameter::PRODUCT:
          case CFunctionParameter::MODIFIER:
            mDimensions[i].setDimension(1, -1, 0, 0, 0); //concentration
            break;

          case CFunctionParameter::VOLUME:
            mDimensions[i].setUnknown(); // TODO Dimension(0, 1, 0); //volume
            break;

          case CFunctionParameter::TIME:
            mDimensions[i].setDimension(0, 0, 1, 0, 0); //time
            break;

          default:
            mDimensions[i].setUnknown();
            break;
        }

      mDimensions[i].fixDimensionless(mD1, mD2, mD3, mD4, mD5);
    }
}

void CFindDimensions::setUseHeuristics(bool flag)
{
  mUseHeuristics = flag;
}

const std::vector<CDimension> & CFindDimensions::getDimensions() const
{
  return mDimensions;
}

void CFindDimensions::findDimensions(CDimension rootDim)
{
  mRootDimension = rootDim;
  findDimensions();
}

void CFindDimensions::findDimensions(bool isMulticompartment)
{
  if (isMulticompartment)
    mRootDimension.setDimension(1, 0, -1, 0, 0); //amount of subs/time
  else
    mRootDimension.setDimension(1, -1, -1, 0, 0); //TODO !!! conc/time

  mRootDimension.fixDimensionless(mD1, mD2, mD3, mD4, mD5);

  findDimensions();
}

std::vector<std::string> CFindDimensions::findDimensionsBoth(const CCopasiDataModel* pDataModel)
{
  //first for single compartment
  findDimensions(false);
  std::vector<CDimension> store = mDimensions;

  //next for multiple compartments
  setupDimensions();
  findDimensions(true);

  //compare...
  std::vector<std::string> ret;
  std::vector<CDimension>::const_iterator it1, it2, it1end = store.end();

  for (it1 = store.begin(), it2 = mDimensions.begin(); it1 != it1end; ++it1, ++it2)
    {
      if (*it1 == *it2)
        ret.push_back(it1->getDisplayString(pDataModel));
      else
        ret.push_back(it1->getDisplayString(pDataModel) + " or " + it2->getDisplayString(pDataModel));
    }

  return ret;
}

void CFindDimensions::findDimensions()
{
  if (!mpFunction) return;

  if (dynamic_cast<const CMassAction*>(mpFunction))
    {
      findDimensionsMassAction();
      return;
    }

  unsigned C_INT32 i, imax = mpFunction->getVariables().size();

  for (i = 0; i < imax; ++i)
    if (mDimensions[i].isUnknown())
      findDimension(i);

  for (i = 0; i < imax; ++i)
    if (mDimensions[i].isUnknown())
      findDimension(i);

  for (i = 0; i < imax; ++i)
    if (mDimensions[i].isUnknown())
      findDimension(i);

  //TODO: conistency check for known dimensions?
}

void CFindDimensions::setChemicalEquation(const CChemEq* eq)
{
  //mpChemEq = eq;
  if (!eq)
    {
      mM1 = 0.0; mM2 = 0.0;
      return;
    }

  mM1 = eq->getMolecularity(CChemEq::SUBSTRATE);
  mM2 = eq->getMolecularity(CChemEq::PRODUCT);
}

void CFindDimensions::setMolecularitiesForMassAction(const unsigned C_INT32 & m1,
    const unsigned C_INT32 & m2)
{
  mM1 = (m1 != C_INVALID_INDEX) ? m1 : -1.0;
  mM2 = (m2 != C_INVALID_INDEX) ? m2 : -1.0;
}

void CFindDimensions::findDimensionsMassAction()
{
  if (mM1 < 0) return;

  CDimension conc; conc.setDimension(1.0, -1.0, 0.0, 0, 0); //TODO!!!

  mRootDimension.fixDimensionless(mD1, mD2, mD3, mD4, mD5);
  conc.fixDimensionless(mD1, mD2, mD3, mD4, mD5);

  if (mDimensions[0].isUnknown())
    {
      mDimensions[0] = mRootDimension - conc * mM1;
    }

  if (mDimensions.size() == 2) return; //irreversible

  if (mDimensions[2].isUnknown())
    {
      mDimensions[2] = mRootDimension - conc * mM2;
    }
}

//find dim for one parameter
void CFindDimensions::findDimension(unsigned C_INT32 index)
{
  if (!mpFunction) return;

  if (index >= mDimensions.size()) return;

  CDimension result;

  //find all variable nodes with given index
  std::vector<const CEvaluationNode*> nodes;
  const std::vector< CEvaluationNode * > & allnodes = mpFunction->getNodeList();
  std::vector< CEvaluationNode * >::const_iterator it, itEnd = allnodes.end();

  for (it = allnodes.begin(); it != itEnd; ++it)
    {
      if ((*it)->getType() == CEvaluationNode::VARIABLE)
        if (dynamic_cast<const CEvaluationNodeVariable*>(*it)->getIndex() == index)
          {
            nodes.push_back(*it);
          }
    }

  //find dimension for all nodes and compare results
  std::vector<const CEvaluationNode * >::const_iterator it2, it2End = nodes.end();

  for (it2 = nodes.begin(); it2 != it2End; ++it2)
    {
      result = result.compare(findDimension(*it2));
    }

  mDimensions[index] = result;
}

CDimension CFindDimensions::findDimension(const CEvaluationNode * node,
    const CEvaluationNode * requestingNode)
{
  CDimension result;

  //variable node
  const CEvaluationNodeVariable* varnode = dynamic_cast<const CEvaluationNodeVariable*>(node);

  if (varnode)
    {
      //is dim known?
      if (!mDimensions[varnode->getIndex()].isUnknown())
        result = mDimensions[varnode->getIndex()];

      else if (requestingNode) //do not try to evaluate dim recursively if asked by parent node
        {
          if (node->getParent())
            result.setUnknown();
          else //no parent
            fatalError();
        }
      else // !requestingNode; method is called from outside, evaluate recursively
        {
          if (node->getParent())
            {
              const CEvaluationNode* parent = dynamic_cast<const CEvaluationNode*>(node->getParent());

              if (parent) result = findDimension(parent, node);
            }
          else //no parent, use root dimension
            {
              result = mRootDimension;
            }
        }
    }

  //operator node
  const CEvaluationNodeOperator* opnode = dynamic_cast<const CEvaluationNodeOperator*>(node);

  if (opnode)
    {
      assert(requestingNode); //should not be called from outside

      switch (opnode->getType() & 0x00FFFFFF)
        {
          case CEvaluationNodeOperator::PLUS:
          case CEvaluationNodeOperator::MINUS:
          {
            CDimension r1, r2;

            if (requestingNode == node->getParent()) //called by parent
              {
                r1 = findDimension(opnode->getLeft(), node);
                r2 = findDimension(opnode->getRight(), node);
              }
            else //called by one of the children
              {
                const CEvaluationNode* parent = dynamic_cast<const CEvaluationNode*>(node->getParent());

                if (parent) r1 = findDimension(parent, node);
                else r1 = mRootDimension;

                if (requestingNode == opnode->getLeft())
                  r2 = findDimension(opnode->getRight(), node);

                if (requestingNode == opnode->getRight())
                  r2 = findDimension(opnode->getLeft(), node);
              }

            //r1 and r2 are now the dimensions of the two nodes that are not the
            //requesting node.
            /*if (r1.isContradiction() || r2.isContradiction())
              result.setContradiction();
            else if (r1 == r2)
              result = r1;
            else if (r1.isUnknown())
              result = r2;
            else if (r2.isUnknown())
              result =r1;
            else
              result.setContradiction();*/
            result = r1.compare(r2);
          }
          break;

          case CEvaluationNodeOperator::MULTIPLY:
          {
            CDimension r1, r2;

            if (requestingNode == node->getParent()) //called by parent
              {
                r1 = findDimension(opnode->getLeft(), node);
                r2 = findDimension(opnode->getRight(), node);

                result = r1 + r2;
              }
            else //called by one of the children
              {
                //parent node
                const CEvaluationNode* parent = dynamic_cast<const CEvaluationNode*>(node->getParent());

                if (parent) r1 = findDimension(parent, node);
                else r1 = mRootDimension;

                //other child
                if (requestingNode == opnode->getLeft())
                  r2 = findDimension(opnode->getRight(), node);

                if (requestingNode == opnode->getRight())
                  r2 = findDimension(opnode->getLeft(), node);

                result = r1 - r2;
              }
          }
          break;

          case CEvaluationNodeOperator::DIVIDE:
          {
            CDimension r1, r2;

            if (requestingNode == node->getParent()) //called by parent
              {
                r1 = findDimension(opnode->getLeft(), node);
                r2 = findDimension(opnode->getRight(), node);

                result = r1 - r2;
              }
            else //called by one of the children
              {
                //parent node
                const CEvaluationNode* parent = dynamic_cast<const CEvaluationNode*>(node->getParent());

                if (parent) r1 = findDimension(parent, node);
                else r1 = mRootDimension;

                //other child
                if (requestingNode == opnode->getLeft())
                  {
                    r2 = findDimension(opnode->getRight(), node);
                    result = r1 + r2;
                  }

                if (requestingNode == opnode->getRight())
                  {
                    r2 = findDimension(opnode->getLeft(), node);
                    result = r2 - r1;
                  }
              }
          }
          break;

          default:
            break;
        }
    }

  //number node
  const CEvaluationNodeNumber* numnode = dynamic_cast<const CEvaluationNodeNumber*>(node);

  if (numnode)
    {
      //heuristics!
      if (mUseHeuristics && (numnode->value() == 1.0))
        result.setDimension(0, 0, 0, 0, 0);
    }

  return result;
}

#ifdef COPASI_DEBUG
void CFindDimensions::printDebugOutput(const CCopasiDataModel* pDataModel) const
{
  std::cout << "mDimensions " << mDimensions.size() << std::endl;
  unsigned C_INT32 i, imax = mDimensions.size();

  for (i = 0; i < imax; ++i)
    std::cout << i << ": " << mDimensions[i].print(pDataModel) << std::endl;
}
#endif // COPASI_DEBUG
