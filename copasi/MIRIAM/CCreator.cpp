// Begin CVS Header
//   $Source: /fs/turing/cvs/copasi_dev/copasi/MIRIAM/CCreator.cpp,v $
//   $Revision: 1.10 $
//   $Name: Build-33 $
//   $Author: shoops $
//   $Date: 2009/02/19 19:50:46 $
// End CVS Header

// Copyright (C) 2008 by Pedro Mendes, Virginia Tech Intellectual
// Properties, Inc., EML Research, gGmbH, University of Heidelberg,
// and The University of Manchester.
// All rights reserved.

// Copyright (C) 2001 - 2007 by Pedro Mendes, Virginia Tech Intellectual
// Properties, Inc. and EML Research, gGmbH.
// All rights reserved.

#include "copasi.h"

#include "CRDFGraph.h"
#include "CCreator.h"
#include "CRDFLiteral.h"
#include "CModelMIRIAMInfo.h"

#include "CopasiDataModel/CCopasiDataModel.h"
#include "model/CModel.h"
#include "report/CKeyFactory.h"
#include "copasi/report/CCopasiRootContainer.h"

CCreator::CCreator(const std::string & objectName,
                   const CCopasiContainer * pParent):
    CCopasiContainer(objectName, pParent, "Creator"),
    mTriplet(),
    mNodePath(),
    mKey(CCopasiRootContainer::getKeyFactory()->add("Creator", this))
{}

CCreator::CCreator(const CRDFTriplet & triplet,
                   const std::string & objectName,
                   const CCopasiContainer * pParent):
    CCopasiContainer(objectName, pParent, "Creator"),
    mTriplet(triplet),
    mNodePath(),
    mKey(CCopasiRootContainer::getKeyFactory()->add("Creator", this))
{
  if (!mTriplet)
    return;

  mNodePath = mTriplet.pObject->getPath();
}

CCreator::CCreator(const CCreator & src,
                   const CCopasiContainer * pParent):
    CCopasiContainer(src, pParent),
    mTriplet(src.mTriplet),
    mNodePath(src.mNodePath),
    mKey(CCopasiRootContainer::getKeyFactory()->add("Creator", this))
{}

CCreator::~CCreator()
{
  CCopasiRootContainer::getKeyFactory()->remove(mKey);
}

const CRDFTriplet & CCreator::getTriplet() const
  {return mTriplet;}

const std::string & CCreator::getKey() const
  {return mKey;}

const std::string & CCreator::getFamilyName() const
  {return mTriplet.pObject->getFieldValue(CRDFPredicate::vcard_Family);}

const std::string & CCreator::getGivenName() const
  {return mTriplet.pObject->getFieldValue(CRDFPredicate::vcard_Given);}

const std::string & CCreator::getEmail() const
  {return mTriplet.pObject->getFieldValue(CRDFPredicate::vcard_EMAIL);}

const std::string & CCreator::getORG() const
  {return mTriplet.pObject->getFieldValue(CRDFPredicate::vcard_Orgname);}

void CCreator::setFamilyName(const std::string & familyName)
{mTriplet.pObject->setFieldValue(familyName, CRDFPredicate::vcard_Family, mNodePath);}

void CCreator::setGivenName(const std::string & givenName)
{mTriplet.pObject->setFieldValue(givenName, CRDFPredicate::vcard_Given, mNodePath);}

void CCreator::setEmail(const std::string & Email)
{mTriplet.pObject->setFieldValue(Email, CRDFPredicate::vcard_EMAIL, mNodePath);}

void CCreator::setORG(const std::string & Orgname)
{mTriplet.pObject->setFieldValue(Orgname, CRDFPredicate::vcard_Orgname, mNodePath);}
