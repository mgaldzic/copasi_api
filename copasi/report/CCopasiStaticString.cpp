/* Begin CVS Header
   $Source: /fs/turing/cvs/copasi_dev/copasi/report/CCopasiStaticString.cpp,v $
   $Revision: 1.13 $
   $Name: Build-33 $
   $Author: shoops $
   $Date: 2006/04/27 01:31:09 $
   End CVS Header */

// Copyright � 2005 by Pedro Mendes, Virginia Tech Intellectual
// Properties, Inc. and EML Research, gGmbH.
// All rights reserved.

#include "copasi.h"

#include <sstream>
#include "CCopasiStaticString.h"
#include "CCopasiObjectName.h"

CCopasiStaticString::CCopasiStaticString(const std::string & name,
    const CCopasiContainer * pParent,
    const std::string & type,
    const unsigned C_INT32 & flag):
    CCopasiObject(name, pParent, type, flag | CCopasiObject::StaticString),
    mStaticString(name)
{}

CCopasiStaticString::CCopasiStaticString(const CCopasiStaticString & src,
    const CCopasiContainer * pParent):
    CCopasiObject(src, pParent),
    mStaticString(src.mStaticString)
{}

CCopasiStaticString::~CCopasiStaticString() {}

CCopasiStaticString & CCopasiStaticString::operator = (const std::string & rhs)
{
  mStaticString = rhs;
  setObjectName(mStaticString);
  return *this;
}

void CCopasiStaticString::print(std::ostream * ostream) const
  {(*ostream) << mStaticString;}

std::string CCopasiStaticString::getObjectDisplayName(bool C_UNUSED(regular) /*=true*/,
    bool C_UNUSED(richtext) /*=false*/) const
  {
    if (mStaticString == "\n")
      return "<linebreak>";

    return "'" + mStaticString + "'";
  }

const std::string & CCopasiStaticString::getStaticString() const
{return mStaticString;}

CCopasiReportSeparator::CCopasiReportSeparator(const std::string & name,
    const CCopasiContainer * pParent):
    CCopasiStaticString(name, pParent, "Separator", CCopasiObject::Separator)
{}

CCopasiReportSeparator::CCopasiReportSeparator(const CCopasiStaticString & src,
    const CCopasiContainer * pParent):
    CCopasiStaticString(src, pParent)
{}

CCopasiReportSeparator::~CCopasiReportSeparator() {}

std::string CCopasiReportSeparator::getObjectDisplayName(bool C_UNUSED(regular) /*=true*/,
    bool C_UNUSED(richtext) /*=false*/) const
  {return getObjectType();}

CCopasiReportSeparator & CCopasiReportSeparator::operator = (const std::string & rhs)
{
  * (CCopasiStaticString *) this = rhs;
  return *this;
}
