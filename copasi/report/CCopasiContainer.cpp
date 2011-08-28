// Begin CVS Header
//   $Source: /fs/turing/cvs/copasi_dev/copasi/report/CCopasiContainer.cpp,v $
//   $Revision: 1.54 $
//   $Name: Build-33 $
//   $Author: shoops $
//   $Date: 2009/10/27 16:52:48 $
// End CVS Header

// Copyright (C) 2008 by Pedro Mendes, Virginia Tech Intellectual
// Properties, Inc., EML Research, gGmbH, University of Heidelberg,
// and The University of Manchester.
// All rights reserved.

// Copyright (C) 2001 - 2007 by Pedro Mendes, Virginia Tech Intellectual
// Properties, Inc. and EML Research, gGmbH.
// All rights reserved.

/*!
    \file CCopasiContainer.cpp
    \brief Implementation file of class CCopasiContainer
 */

/**
 * Class CCopasiContainer
 *
 * This class is the is used to group CCopasiObjects logically. It inself is
 * an object. Contained objects are still globally accessible.
 *
 * Copyright Stefan Hoops 2002
 */

#include "copasi.h"
#include "CCopasiObjectName.h"
#include "CCopasiContainer.h"
#include "CCopasiObjectReference.h"
#include "CCopasiStaticString.h"
#include "CCopasiTimer.h"

#include "utilities/CCopasiVector.h"
#include "report/CCopasiRootContainer.h"

const std::vector< CCopasiContainer * > CCopasiContainer::EmptyList;

CCopasiContainer::CCopasiContainer() :
    CCopasiObject(),
    mObjects()
{addObjectReference("Name", *const_cast<std::string *>(&getObjectName()));}

CCopasiContainer::CCopasiContainer(const std::string & name,
                                   const CCopasiContainer * pParent,
                                   const std::string & type,
                                   const unsigned C_INT32 & flag):
    CCopasiObject(name, pParent, type, flag | CCopasiObject::Container),
    mObjects()
{addObjectReference("Name", *const_cast<std::string *>(&getObjectName()));}

CCopasiContainer::CCopasiContainer(const CCopasiContainer & src,
                                   const CCopasiContainer * pParent):
    CCopasiObject(src, pParent),
    mObjects()
{addObjectReference("Name", *const_cast<std::string *>(&getObjectName()));}

CCopasiContainer::~CCopasiContainer()
{
  objectMap::iterator it = mObjects.begin();
  objectMap::iterator end = mObjects.end();

  for (; it != end; it++)
    if (it->second != NULL &&
        it->second->getObjectParent() == this)
      {
        it->second->setObjectParent(NULL);
        pdelete(it->second);
      }
}

const CCopasiObject * CCopasiContainer::getObject(const CCopasiObjectName & cn) const
{
  if (cn == "")
    {
      if (isRoot())
        return NULL;
      else
        return this;
    }

  std::string Name = cn.getObjectName();
  std::string Type = cn.getObjectType();

  if (getObjectName() == Name && getObjectType() == Type)
    return getObject(cn.getRemainder());

  //check if the first part of the cn matches one of the children (by name and type)
  std::pair< objectMap::const_iterator, objectMap::const_iterator > range =
    mObjects.equal_range(Name);

  objectMap::const_iterator it = range.first;

  while (it != range.second && it->second->getObjectType() != Type) ++it;

  if (it == range.second) //not found in the list of children
    {
      if (Type == "String")
        return new CCopasiStaticString(Name, this);
      else if (Type == "Separator")
        return new CCopasiReportSeparator(Name, this);
      else
        return NULL;
    }

  const CCopasiObject * pObject = NULL;

  if (it->second->isNameVector() || it->second->isVector())
    {
      if (cn.getElementName(0, false) == "")
        return it->second;

      pObject = it->second->getObject("[" + cn.getElementName(0, false) + "]");

      if (it->second->getObjectType() == "Reference" ||
          !pObject ||
          cn.getRemainder() == "")
        return pObject;
      else
        return pObject->getObject(cn.getRemainder());
    }

  //handle objects where the array flag is set. Currently this applies to the
  //CArrayAnnotation object. Since this is also a container, we have to do this
  //before handling general containers.
  if (it->second->isArray())
    {
      //we need to call the getObject() method of the child array with the
      //remainder of the cn, with the indices in square brackets, or with an empty string

      //if there are no indices there could still be a remainder (since the array can also be
      //a container)
      if (cn.getElementName(0, false) == "") //no indices
        return it->second->getObject(cn.getRemainder());

      //get the indices from the CN
      std::string indices;
      std::string tmp;
      C_INT32 ii = 0;

      while ((tmp = cn.getElementName(ii, false)) != "")
        {
          indices += "[" + tmp + "]";
          ++ii;
        }

      //try to get the array element from the indices
      pObject = it->second->getObject(indices);

      //if the element could not be resolved, just return NULL. If there is no
      //remainder, just return the array element.
      //In all other cases we call getObject() on the array element with the remainder
      //of the CN. The special treatment of the empty remainder is probably necessary
      //since not all implementations of getObject(cn) handle the case of empty CN.
      if (!pObject)
        return NULL;

      if (cn.getRemainder() == "")
        return pObject;
      else
        return pObject->getObject(cn.getRemainder());
    }

  //handle generic containers.
  if (it->second->isContainer())
    return it->second->getObject(cn.getRemainder());

  if (it->second->isMatrix())
    {
      if (cn.getElementName(0, false) == "")
        return it->second;

      pObject = it->second->getObject("[" + cn.getElementName(0, false) + "]" +                   //TODO really?
                                      "[" + cn.getElementName(1, false) + "]");

      if (it->second->getObjectType() == "Reference" || !pObject)
        return pObject;
      else
        return pObject->getObject(cn.getRemainder());
    }

  if (it->second->isReference() || it->second->isStaticString() ||
      cn.getRemainder() == "")
    return it->second;

  return it->second->getObject(cn.getRemainder());
}

const CCopasiContainer::objectMap & CCopasiContainer::getObjects() const
{return mObjects;}

const CCopasiObject * CCopasiContainer::getValueObject() const
{
  void * ptr = getValuePointer();

  if (ptr == NULL) return NULL;

  std::multimap< const std::string, CCopasiObject * >::const_iterator it = mObjects.begin();
  std::multimap< const std::string, CCopasiObject * >::const_iterator end = mObjects.end();

  for (; it != end; ++it)
    if (ptr == it->second->getValuePointer()) return it->second;

  return NULL;
}

void CCopasiContainer::initObjects() {}

bool CCopasiContainer::add(CCopasiObject * pObject,
                           const bool & adopt)
{
  /* We check whether we are already containing that object. */
  std::pair< objectMap::iterator, objectMap::iterator > range =
    mObjects.equal_range(pObject->getObjectName());
  objectMap::iterator it;

  for (it = range.first; it != range.second; ++it)
    if (it->second == pObject) break;

  if (it != range.second) return false;

  /* This object is not contained, so we can add it. */
  mObjects.insert
  (std::pair<const std::string, CCopasiObject * >(pObject->getObjectName(),
      pObject));

  if (adopt) pObject->setObjectParent(this);

  return true;
}

bool CCopasiContainer::remove(CCopasiObject * pObject)
{
  objectMap::iterator it = mObjects.begin();
  objectMap::iterator end = mObjects.end();

  /*
    std::pair< objectMap::iterator, objectMap::iterator > range =
      mObjects.equal_range(pObject->getObjectName());
    objectMap::iterator it;
  */
  for (; it != end; ++it)
    if (it->second == pObject) break;

  if (it == end) return false;

  mObjects.erase(it);

  return true;
}

// virtual
std::string CCopasiContainer::getUnits() const
{return "";}

// virtual
std::string CCopasiContainer::getChildObjectUnits(const CCopasiObject * /* pObject */) const
{return "";}
