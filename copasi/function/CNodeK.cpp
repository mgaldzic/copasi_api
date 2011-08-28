/* Begin CVS Header
   $Source: /fs/turing/cvs/copasi_dev/copasi/function/CNodeK.cpp,v $
   $Revision: 1.30 $
   $Name: Build-33 $
   $Author: shoops $
   $Date: 2006/04/27 01:28:26 $
   End CVS Header */

// Copyright � 2005 by Pedro Mendes, Virginia Tech Intellectual
// Properties, Inc. and EML Research, gGmbH.
// All rights reserved.

// CNodeK.cpp : classes for function tree
//
/////////////////////////////////////////////////////////////////////////////

#include "mathematics.h"
#include <stdio.h>

#include "copasi.h"
#include "utilities/CReadConfig.h"
#include "utilities/CCopasiMessage.h"
#include "utilities/utility.h"

#include "CNodeK.h"

CNodeK::CNodeK()
{
  CONSTRUCTOR_TRACE;
  mType = N_NOP;
  mSubtype = N_NOP;
  mLeft = NULL;
  mRight = NULL;
  mConstant = 0.0;
  mIndex = -1;
  mOldIndex = -1;
}

CNodeK::CNodeK(const CNodeK & src)
{
  CONSTRUCTOR_TRACE;
  mType = src.mType;
  mSubtype = src.mSubtype;
  mLeft = src.mLeft;
  mRight = src.mRight;
  mConstant = src.mConstant;
  mName = src.mName;
  mIndex = src.mIndex;
  mOldIndex = src.mOldIndex;
}

CNodeK::CNodeK(char type, char subtype)
{
  CONSTRUCTOR_TRACE;
  mType = type;
  mSubtype = subtype;
  mLeft = NULL;
  mRight = NULL;
  mConstant = 0.0;
  mIndex = -1;
  mOldIndex = -1;
}

CNodeK::CNodeK(const std::string & name)
{
  CONSTRUCTOR_TRACE;
  mType = N_IDENTIFIER;
  mSubtype = N_NOP;
  mLeft = NULL;
  mRight = NULL;
  mConstant = 0.0;
  mName = name;
  mIndex = -1;
  mOldIndex = -1;
}

CNodeK::CNodeK(C_FLOAT64 constant)
{
  CONSTRUCTOR_TRACE;
  mType = N_NUMBER;
  mSubtype = N_NOP;
  mLeft = NULL;
  mRight = NULL;
  mConstant = constant;
  mIndex = -1;
  mOldIndex = -1;
}

void CNodeK::cleanup()
{}

CNodeK::~CNodeK()
{
  DESTRUCTOR_TRACE;
}

C_INT32 CNodeK::load(CReadConfig & configbuffer)
{
  C_INT32 Fail = 0;

  if ((Fail = configbuffer.getVariable("Node", "node", &mType, &mSubtype,
                                       CReadConfig::SEARCH)))
    return Fail;

  /* This COPASI treats all these as identifiers */
  if (mType == N_SUBSTRATE ||
      mType == N_PRODUCT ||
      mType == N_MODIFIER ||
      mType == N_KCONSTANT)
    {
      mSubtype = mType;
      mType = N_IDENTIFIER;
    }

  // leave the Left & Right pointers out
  // value of the constant if one
  if (mType == N_NUMBER)
    {
      if ((Fail = configbuffer.getVariable("Value", "C_FLOAT64", &mConstant)))
        return Fail;
    }
  else if (mType == N_IDENTIFIER)
    {
      if ((Fail = configbuffer.getVariable("Index", "C_INT32", &mIndex)))
        return Fail;
      if ((Fail = configbuffer.getVariable("Name", "string", &mName)))
        return Fail;
    }

  return Fail;
}

/*
std::string CNodeK::getExplicitFunctionString(const std::vector< std::vector< std::string > > & callParameterNames,
    const std::string &r)
{
  char fstr[256];
  switch (mType)
    {
    case N_ROOT:
      return mLeft->getExplicitFunctionString(callParameterNames, r);
    case N_NUMBER:
      sprintf(fstr, "%-20g", mConstant);
      mExplicitFunction = fstr;
      break;
    case N_IDENTIFIER:
      FixSName(callParameterNames[mIndex][0], mExplicitFunction);
      if (mSubtype == N_KCONSTANT)
        mExplicitFunction += r;
      break;
    case N_OPERATOR:
      switch (mSubtype)
        {
        case '+':
          mExplicitFunction = mLeft->getExplicitFunctionString(callParameterNames, r) + "+"
                              + mRight->getExplicitFunctionString(callParameterNames, r);
          break;
        case '-':
          mExplicitFunction = mLeft->getExplicitFunctionString(callParameterNames, r) + "-"
                              + mRight->getExplicitFunctionString(callParameterNames, r);
          break;
        case '*':
          mExplicitFunction = "(" + mLeft->getExplicitFunctionString(callParameterNames, r)
                              + ")*(" + mRight->getExplicitFunctionString(callParameterNames, r)
                              + ")";
          break;
        case '/':
          mExplicitFunction = "(" + mLeft->getExplicitFunctionString(callParameterNames, r)
                              + ")/(" + mRight->getExplicitFunctionString(callParameterNames, r)
                              + ")";
          break;
        case '^':
          mExplicitFunction = "(" + mLeft->getExplicitFunctionString(callParameterNames, r)
                              + ")^(" + mRight->getExplicitFunctionString(callParameterNames, r)
                              + ")";
          break;
        default:
          mExplicitFunction.empty();
        }
      break;
    case N_FUNCTION:
      switch (mSubtype)
        {
        case '+':
          mExplicitFunction = "+(" + mLeft->getExplicitFunctionString(callParameterNames, r)
                              + ")";
          break;
        case '-':
          mExplicitFunction = "-(" + mLeft->getExplicitFunctionString(callParameterNames, r)
                              + ")";
          break;
        case N_EXP:
          mExplicitFunction = "e^(" + mLeft->getExplicitFunctionString(callParameterNames, r)
                              + ")";
          break;
        case N_LOG:
          mExplicitFunction = "log10(" + mLeft->getExplicitFunctionString(callParameterNames, r)
                              + ")";
          break;
        case N_LOG10:
          mExplicitFunction = "log(" + mLeft->getExplicitFunctionString(callParameterNames, r)
                              + ")";
          break;
        case N_SIN:
          mExplicitFunction = "sin(" + mLeft->getExplicitFunctionString(callParameterNames, r)
                              + ")";
          break;
        case N_COS:
          mExplicitFunction = "cos(" + mLeft->getExplicitFunctionString(callParameterNames, r)
                              + ")";
          break;
        case N_TAN:
          mExplicitFunction = "tan(" + mLeft->getExplicitFunctionString(callParameterNames, r)
                              + ")";
          break;
        case N_SEC:
          mExplicitFunction = "sec(" + mLeft->getExplicitFunctionString(callParameterNames, r)
                              + ")";
          break;
        case N_CSC:
          mExplicitFunction = "csc(" + mLeft->getExplicitFunctionString(callParameterNames, r)
                              + ")";
          break;
        case N_COT:
          mExplicitFunction = "cot(" + mLeft->getExplicitFunctionString(callParameterNames, r)
                              + ")";
          break;
        case N_SINH:
          mExplicitFunction = "sinh(" + mLeft->getExplicitFunctionString(callParameterNames, r)
                              + ")";
          break;
        case N_COSH:
          mExplicitFunction = "cosh(" + mLeft->getExplicitFunctionString(callParameterNames, r)
                              + ")";
          break;
        case N_TANH:
          mExplicitFunction = "tanh(" + mLeft->getExplicitFunctionString(callParameterNames, r)
                              + ")";
          break;
        case N_SECH:
          mExplicitFunction = "sech(" + mLeft->getExplicitFunctionString(callParameterNames, r)
                              + ")";
          break;
        case N_CSCH:
          mExplicitFunction = "csch(" + mLeft->getExplicitFunctionString(callParameterNames, r)
                              + ")";
          break;
        case N_COTH:
          mExplicitFunction = "coth(" + mLeft->getExplicitFunctionString(callParameterNames, r)
                              + ")";
          break;
        case N_ARCSIN:
          mExplicitFunction = "arcsin(" + mLeft->getExplicitFunctionString(callParameterNames, r)
                              + ")";
          break;
        case N_ARCCOS:
          mExplicitFunction = "arccos(" + mLeft->getExplicitFunctionString(callParameterNames, r)
                              + ")";
          break;
        case N_ARCTAN:
          mExplicitFunction = "arctan(" + mLeft->getExplicitFunctionString(callParameterNames, r)
                              + ")";
          break;
        case N_ARCSEC:
          mExplicitFunction = "arcsec(" + mLeft->getExplicitFunctionString(callParameterNames, r)
                              + ")";
          break;
        case N_ARCCSC:
          mExplicitFunction = "arccsc(" + mLeft->getExplicitFunctionString(callParameterNames, r)
                              + ")";
          break;
        case N_ARCCOT:
          mExplicitFunction = "arccot(" + mLeft->getExplicitFunctionString(callParameterNames, r)
                              + ")";
          break;
        case N_ARCSINH:
          mExplicitFunction = "arcsinh(" + mLeft->getExplicitFunctionString(callParameterNames, r)
                              + ")";
          break;
        case N_ARCCOSH:
          mExplicitFunction = "arccosh(" + mLeft->getExplicitFunctionString(callParameterNames, r)
                              + ")";
          break;
        case N_ARCTANH:
          mExplicitFunction = "arctanh(" + mLeft->getExplicitFunctionString(callParameterNames, r)
                              + ")";
          break;
        case N_ARCSECH:
          mExplicitFunction = "arcsech(" + mLeft->getExplicitFunctionString(callParameterNames, r)
                              + ")";
          break;
        case N_ARCCSCH:
          mExplicitFunction = "arccsch(" + mLeft->getExplicitFunctionString(callParameterNames, r)
                              + ")";
          break;
        case N_ARCCOTH:
          mExplicitFunction = "arccoth(" + mLeft->getExplicitFunctionString(callParameterNames, r)
                              + ")";
          break;
        case N_ABS:
          mExplicitFunction = "abs(" + mLeft->getExplicitFunctionString(callParameterNames, r)
                              + ")";
          break;
        case N_SQRT:
          mExplicitFunction = "sqrt(" + mLeft->getExplicitFunctionString(callParameterNames, r)
                              + ")";
          break;
        case N_GAUSS:
        case N_BOLTZ:
        case N_RND:
          mExplicitFunction = "(" + mLeft->getExplicitFunctionString(callParameterNames, r)
                              + ")";
          break;
        default:
          mExplicitFunction.empty();
        }
      break;
    default:
      mExplicitFunction.empty();
    }
  return mExplicitFunction;
}
 */

char CNodeK::getType() const
  {
    return mType;
  }

char CNodeK::getSubtype() const
  {
    return mSubtype;
  }

CNodeK & CNodeK::getLeft() const
  {
    if (!mLeft)
      fatalError(); // Call LeftIsValid first to avoid this!
    return *mLeft;
  }

CNodeK & CNodeK::getRight() const
  {
    if (!mRight)
      fatalError(); // Call RightIsValid first to avoid this!
    return *mRight;
  }

std::string CNodeK::getName() const
  {
    return mName;
#ifdef XXXX

    static unsigned C_INT ctr = 0;
    char name[9];
    if (isIdentifier())
      return mName;
    else
      {
        sprintf(name, "%X", ctr++);
        return name;
      }
#endif // XXXX
  }

C_FLOAT64 CNodeK::getConstant() const
  {
    return mConstant;
  }

C_INT32 CNodeK::getIndex() const
  {
    return mIndex;
  }

void CNodeK::setType(char type)
{
  mType = type;
}

void CNodeK::setSubtype(char subtype)
{
  mSubtype = subtype;
}

void CNodeK::setLeft(CNodeK & left)
{
  mLeft = &left;
}

void CNodeK::setLeft(CNodeK * pleft)
{
  mLeft = pleft;
}

void CNodeK::setRight(CNodeK & right)
{
  mRight = &right;
}

void CNodeK::setRight(CNodeK * pright)
{
  mRight = pright;
}

void CNodeK::setName(const std::string & name)
{
  mName = name;
}

void CNodeK::setConstant(C_FLOAT64 & constant)
{
  mConstant = constant;
}

void CNodeK::setIndex(C_INT32 index)
{
  mIndex = index;
}

void CNodeK::setOldIndex(C_INT32 oldindex)
{
  mOldIndex = oldindex;
}

C_INT16 CNodeK::isLeftValid() const
  {
    return (mLeft != NULL);
  }

C_INT16 CNodeK::isRightValid() const
  {
    return (mRight != NULL);
  }

C_INT16 CNodeK::isNumber() const
  {
    return (mType == N_NUMBER);
  }

C_INT16 CNodeK::isIdentifier() const
  {
    switch (mType)
      {
      case N_OBJECT:
      case N_IDENTIFIER:
      case N_SUBSTRATE:
      case N_PRODUCT:
      case N_MODIFIER:
      case N_KCONSTANT:
      case N_VOLUME:
        return true;
      default:
        return false;
      }
  }

C_INT16 CNodeK::isOperator() const
  {
    return mType == N_OPERATOR;
  }

C_INT16 CNodeK::leftPrecedence() const
  {
    switch (mType)
      {
      case N_OBJECT:
      case N_NUMBER:
      case N_IDENTIFIER:
      case N_FUNCTION:
        return 5;
      }
    // if we got here then it is an operator
    switch (mSubtype)
      {
      case '+':
      case '-':
        return 1;
      case '*':
      case '/':
        return 3;
      case '(':
        return 6;
      case '^':
        return 5;
      case ')':
      case '%':
        return 0;
      }
    return 0;
  }

C_INT16 CNodeK::rightPrecedence() const
  {
    switch (mType)
      {
      case N_OBJECT:
      case N_NUMBER:
      case N_IDENTIFIER:
        return 6;
      case N_FUNCTION:
        return 4;
      }
    // if we got here then it is an operator
    switch (mSubtype)
      {
      case '+':
      case '-':
        return 2;
      case '*':
      case '/':
        return 4;
      case ')':
        return 6;
      case '^':
        return 4;
      case '(':
      case '%':
        return 0;
      }
    return 0;
  }

/*
C_FLOAT64 CNodeK::value(const CCallParameters<C_FLOAT64> & callParameters) const
  {
    // if it is a constant or an identifier just return its value
    if (isNumber())
      return mConstant;
    switch (mType)
      {
      case N_OBJECT:
        return *(double*)((CCopasiObject*)mLeft)->getValuePointer();
        break;
      case N_IDENTIFIER:
        return * callParameters[mIndex].value;
        break;
      case N_OPERATOR:
        switch (mSubtype)
          {
          case '+':
            return mLeft->value(callParameters) + mRight->value(callParameters);
          case '-':
            return mLeft->value(callParameters) - mRight->value(callParameters);
          case '*':
            return mLeft->value(callParameters) * mRight->value(callParameters);
          case '/':
            return mLeft->value(callParameters) / mRight->value(callParameters);
          case '^':
            return pow(mLeft->value(callParameters), mRight->value(callParameters));
          default:
            fatalError();   // THROW EXCEPTION
            return 0.0;
          }
        break;
      case N_FUNCTION:
        switch (mSubtype)
          {
          case '+':
            return mLeft->value(callParameters);
          case '-':
            return - mLeft->value(callParameters);
          case N_EXP:
            return exp(mLeft->value(callParameters));
          case N_LOG:
            return log(mLeft->value(callParameters));
          case N_LOG10:
            return log10(mLeft->value(callParameters));


          case N_SIN:
            return sin(mLeft->value(callParameters));
          case N_COS:
            return cos(mLeft->value(callParameters));
          case N_TAN:
            return tan(mLeft->value(callParameters));


          case N_SEC:
            return 1 / cos(mLeft->value(callParameters));
          case N_CSC:
            return 1 / sin(mLeft->value(callParameters));
          case N_COT:
            return 1 / tan(mLeft->value(callParameters));


          case N_SINH:
            return sinh(mLeft->value(callParameters));
          case N_COSH:
            return cosh(mLeft->value(callParameters));
          case N_TANH:
            return tanh(mLeft->value(callParameters));


          case N_SECH:
            return 1 / cosh(mLeft->value(callParameters));
          case N_CSCH:
            return 1 / sinh(mLeft->value(callParameters));
          case N_COTH:
            return 1 / tanh(mLeft->value(callParameters));


          case N_ARCSIN:
            return asin(mLeft->value(callParameters));
          case N_ARCCOS:
            return acos(mLeft->value(callParameters));
          case N_ARCTAN:
            return atan(mLeft->value(callParameters));


          case N_ARCSEC:   //TODO
            return acos(1 / mLeft->value(callParameters));
          case N_ARCCSC:   //TODO
            return asin(1 / mLeft->value(callParameters));
          case N_ARCCOT:   //TODO
            return atan(1 / mLeft->value(callParameters));


          case N_ARCSINH:
            return asinh(mLeft->value(callParameters));
          case N_ARCCOSH:
            return acosh(mLeft->value(callParameters));
          case N_ARCTANH:
            return atanh(mLeft->value(callParameters));


          case N_ARCSECH:
            return acosh(1 / mLeft->value(callParameters));
          case N_ARCCSCH:
            return asinh(1 / mLeft->value(callParameters));
          case N_ARCCOTH:
            return atanh(1 / mLeft->value(callParameters));


          case N_ABS:
            return fabs(mLeft->value(callParameters));
          case N_SQRT:
            return sqrt(mLeft->value(callParameters));


          default:
            fatalError();   // THROW EXCEPTION
            return 0.0;
          }
        break;
      default:
        fatalError();   // THROW EXCEPTION
        return 0.0;
      }
    fatalError();   // THROW EXCEPTION
    return 0.0;
  }
 */

/*
#define SPC(level) std::string(level, ' ')


void CNodeK::writeMathML(std::ostream & out, C_INT32 level) const
  {
    bool flag = false;


    switch (mType)
      {
      case N_NUMBER:
        out << SPC(level) << "<mn>" << mConstant << "</mn>" << std::endl;
        break;
        //    case N_OBJECT:
        //      return *(double*)((CCopasiObject*)mLeft)->getReference();
        //      break;
      case N_IDENTIFIER:       //do some heuristics for indentifiers starting with "K" or "V"
        out << SPC(level);
        if (mName.substr(0, 1) == "K")
          out << "<msub><mi>K</mi><mi>" << mName.substr(1) << "</mi></msub>" << std::endl;
        else if (mName.substr(0, 1) == "V")
          out << "<msub><mi>V</mi><mi>" << mName.substr(1) << "</mi></msub>" << std::endl;
        else
          out << "<mi>" << mName << "</mi>" << std::endl;
        break;
      case N_OPERATOR:
        switch (mSubtype)
          {
          case '+':
            out << SPC(level) << "<mrow>" << std::endl;
            mLeft->writeMathML(out, level + 1);
            out << SPC(level + 1) << "<mo>" << "+" << "</mo>" << std::endl;
            mRight->writeMathML(out, level + 1);
            out << SPC(level) << "</mrow>" << std::endl;
            break;
          case '-':
            out << SPC(level) << "<mrow>" << std::endl;
            mLeft->writeMathML(out, level + 1);
            out << SPC(level + 1) << "<mo>" << "-" << "</mo>" << std::endl;


            //do we need "()" ?
            flag = (mRight->mType == N_OPERATOR) && ((mRight->mSubtype == '-') || (mRight->mSubtype == '+'));
            if (flag)
              {
                out << SPC(level + 1) << "<mfenced>" << std::endl;
              }
            mRight->writeMathML(out, level + 1);
            if (flag)
              {
                out << SPC(level + 1) << "</mfenced>" << std::endl;
              }
            out << SPC(level) << "</mrow>" << std::endl;
            break;
          case '*':
            out << SPC(level) << "<mrow>" << std::endl;


            //do we need "()" ?
            flag = (mLeft->mType == N_OPERATOR) && ((mLeft->mSubtype == '-') || (mLeft->mSubtype == '+'));
            if (flag)
              {
                out << SPC(level + 1) << "<mfenced>" << std::endl;
              }
            mLeft->writeMathML(out, level + 1);
            if (flag)
              {
                out << SPC(level + 1) << "</mfenced>" << std::endl;
              }
            out << SPC(level + 1) << "<mo>" << "&CenterDot;" << "</mo>" << std::endl;
            flag = (mRight->mType == N_OPERATOR) && ((mRight->mSubtype == '-') || (mRight->mSubtype == '+'));
            if (flag)
              {
                out << SPC(level) << "<mfenced>" << std::endl;
              }
            mRight->writeMathML(out, level + 1);
            if (flag)
              {
                out << SPC(level + 1) << "</mfenced>" << std::endl;
              }
            out << SPC(level) << "</mrow>" << std::endl;
            break;
          case '^':
            out << SPC(level) << "<msup>" << std::endl;


            //do we need "()" ?
            flag = (mLeft->mType == N_OPERATOR) && ((mLeft->mSubtype == '-') || (mLeft->mSubtype == '+')
                                                    || (mLeft->mSubtype == '*') || (mLeft->mSubtype == '/')
                                                    || (mLeft->mSubtype == '^'));
            if (flag)
              {
                out << SPC(level + 1) << "<mfenced>" << std::endl;
              }
            mLeft->writeMathML(out, level + 2);
            if (flag)
              {
                out << SPC(level + 1) << "</mfenced>" << std::endl;
              }


            out << SPC(level + 1) << "<mrow>" << std::endl;
            mRight->writeMathML(out, level + 2);
            out << SPC(level + 1) << "</mrow>" << std::endl;


            out << SPC(level) << "</msup>" << std::endl;
            break;
          case '/':
            out << SPC(level) << "<mfrac>" << std::endl;


            out << SPC(level + 1) << "<mrow>" << std::endl;
            mLeft->writeMathML(out, level + 2);
            out << SPC(level + 1) << "</mrow>" << std::endl;


            out << SPC(level + 1) << "<mrow>" << std::endl;
            mRight->writeMathML(out, level + 2);
            out << SPC(level + 1) << "</mrow>" << std::endl;


            out << SPC(level) << "</mfrac>" << std::endl;
            break;
          }
        break;
      case N_FUNCTION:
        switch (mSubtype)
          {
          case '+':     //do nothing
            mLeft->writeMathML(out, level);
            break;
          case '-':
            out << SPC(level) << "<mrow>" << std::endl;
            out << SPC(level + 1) << "<mo>" << "-" << "</mo>" << std::endl;


            //do we need "()" ?
            flag = (mLeft->mType == N_OPERATOR) && ((mLeft->mSubtype == '-') || (mLeft->mSubtype == '+'));
            if (flag)
              {
                out << SPC(level + 1) << "<mfenced>" << std::endl;
              }
            mLeft->writeMathML(out, level + 1);
            if (flag)
              {
                out << SPC(level + 1) << "</mfenced>" << std::endl;
              }
            out << SPC(level) << "</mrow>" << std::endl;
            break;
          case N_EXP:
            //return exp(mLeft->value(callParameters));
          case N_LOG:
            //return log(mLeft->value(callParameters));
          case N_LOG10:
            //return log10(mLeft->value(callParameters));
          case N_SIN:
            //return sin(mLeft->value(callParameters));
          case N_COS:
            //return cos(mLeft->value(callParameters));
          case N_TAN:
            //return cos(mLeft->value(callParameters));
          default:
            //fatalError();   // THROW EXCEPTION
            //return 0.0;
            break;
          }
        break;
        //default:
        //fatalError();   // THROW EXCEPTION
        //return 0.0;
      }
  }


#undef SPC
 */
