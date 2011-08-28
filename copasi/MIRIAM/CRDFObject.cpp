// Begin CVS Header
//   $Source: /fs/turing/cvs/copasi_dev/copasi/MIRIAM/CRDFObject.cpp,v $
//   $Revision: 1.8 $
//   $Name: Build-33 $
//   $Author: shoops $
//   $Date: 2008/06/10 20:31:11 $
// End CVS Header

// Copyright (C) 2008 by Pedro Mendes, Virginia Tech Intellectual
// Properties, Inc., EML Research, gGmbH, University of Heidelberg,
// and The University of Manchester.
// All rights reserved.

// Copyright (C) 2001 - 2007 by Pedro Mendes, Virginia Tech Intellectual
// Properties, Inc. and EML Research, gGmbH.
// All rights reserved.

#include "copasi.h"

#include "CRDFObject.h"
#include "CRDFLiteral.h"

CRDFObject::CRDFObject():
    mType(RESOURCE),
    mResource(""),
    mBlankNodeId(""),
    mIsLocalResource(false),
    mpLiteral(NULL)
{}

CRDFObject::CRDFObject(const CRDFObject & src):
    mType(src.mType),
    mResource(src.mResource),
    mBlankNodeId(src.mBlankNodeId),
    mIsLocalResource(src.mIsLocalResource),
    mpLiteral(src.mpLiteral == NULL ? NULL : new CRDFLiteral(*src.mpLiteral))
{}

CRDFObject& CRDFObject::operator =(const CRDFObject& rhs)
{
  if (this != &rhs)
    {
      mType = rhs.mType;
      mResource = rhs.mResource;
      mBlankNodeId = rhs.mBlankNodeId;
      mIsLocalResource = rhs.mIsLocalResource;
      mpLiteral = (rhs.mpLiteral == NULL ? NULL : new CRDFLiteral(*rhs.mpLiteral));
    }
  return *this;
}

CRDFObject::~CRDFObject()
{pdelete(mpLiteral);}

bool CRDFObject::operator ==(const CRDFObject& rhs) const
  {
    if (mType == rhs.mType)
      {
        switch (mType)
          {
          case CRDFObject::BLANK_NODE:
            if (mBlankNodeId == rhs.mBlankNodeId)
              return true;
            break;

          case CRDFObject::RESOURCE:
            if (mResource == rhs.mResource && mIsLocalResource == rhs.mIsLocalResource)
              return true;
            break;

          case CRDFObject::LITERAL:
            if (mpLiteral != NULL && rhs.mpLiteral != NULL)
              return false; // TODO (*mpLiteral == *rhs.mpLiteral);
            else if (mpLiteral == NULL && rhs.mpLiteral == NULL)
              return true;
            break;
          }
      }
    return false;
  }

void CRDFObject::setType(const CRDFObject::eObjectType & type)
{mType = type;}

const CRDFObject::eObjectType & CRDFObject::getType() const
  {return mType;}

void CRDFObject::setResource(const std::string & resource, const bool & isLocal)
{
  mResource = resource;
  mIsLocalResource = isLocal;
}

const std::string & CRDFObject::getResource() const
  {return mResource;}

const bool & CRDFObject::isLocal() const
  {return mIsLocalResource;}

void CRDFObject::setBlankNodeId(const std::string & blankNodeId)
{mBlankNodeId = blankNodeId;}

const std::string & CRDFObject::getBlankNodeID() const
  {return mBlankNodeId;}

void CRDFObject::setLiteral(const CRDFLiteral & literal)
{
  pdelete(mpLiteral);
  mpLiteral = new CRDFLiteral(literal);
}

CRDFLiteral & CRDFObject::getLiteral()
{
  if (mpLiteral == NULL)
    mpLiteral = new CRDFLiteral;

  return *mpLiteral;
}

const CRDFLiteral & CRDFObject::getLiteral() const
  {
    if (mpLiteral == NULL)
      mpLiteral = new CRDFLiteral;

    return *mpLiteral;
  }

void CRDFObject::clearData()
{
  mType = RESOURCE;
  mResource = "";
  mBlankNodeId = "";
  mIsLocalResource = false;
  mpLiteral = NULL;
}

std::ostream & operator << (std::ostream & os, const CRDFObject & object)
{
  switch (object.getType())
    {
    case CRDFObject::BLANK_NODE:
      os << object.getBlankNodeID();
      break;

    case CRDFObject::RESOURCE:
      os << object.getResource();
      break;

    case CRDFObject::LITERAL:
      os << object.getLiteral();
      break;
    }

  return os;
}
