// Begin CVS Header
//   $Source: /fs/turing/cvs/copasi_dev/copasi/report/CArrayElementReference.cpp,v $
//   $Revision: 1.5 $
//   $Name: Build-33 $
//   $Author: shoops $
//   $Date: 2009/01/07 19:04:15 $
// End CVS Header

// Copyright (C) 2008 by Pedro Mendes, Virginia Tech Intellectual
// Properties, Inc., EML Research, gGmbH, University of Heidelberg,
// and The University of Manchester.
// All rights reserved.

#include "CArrayElementReference.h"
#include "CCopasiContainer.h"
#include "utilities/CAnnotatedMatrix.h"

CArrayElementReference::CArrayElementReference(const std::string & index, const CCopasiContainer * pParent)
    : CCopasiObject(index, pParent, "ElementReference",
                    CCopasiObject::Reference |
                    CCopasiObject::NonUniqueName |
                    CCopasiObject::ValueDbl),
    //    mpReference(NULL),
    mIndex(index)
{
  assert(pParent != NULL);
}

void CArrayElementReference::updateMethod(const CCopasiAbstractArray::data_type & /*value*/)
{
  //  if (mpReference)
  //    *mpReference = value;
}

void * CArrayElementReference::getValuePointer() const
  {
    CArrayAnnotation * tmpAA = dynamic_cast<CArrayAnnotation*>(getObjectParent());
    if (!tmpAA) return NULL;

    //now get the array indices. At the moment only numerical indices...
    //this could be done in the constructor, actually
    unsigned C_INT32 ii = 0;
    CCopasiArray::index_type index;
    std::string tmpIndexString;
    while ((tmpIndexString = mIndex.getElementName(ii, false)) != "")
      {
        if (tmpIndexString == ".")
          break; //"." indicates 0-dimensional array. This means index should stay empty
        index.push_back(mIndex.getElementIndex(ii));
        ++ii;
      }

    if (index.size() != tmpAA->dimensionality())  //wrong number of indices for this array
      return NULL;

    for (ii = 0; ii < tmpAA->dimensionality(); ++ii)
      if (index[ii] >= tmpAA->size()[ii]) //out of range
        return NULL;

    return &(*tmpAA->array())[index];

    //TODO
    //TODO perhaps we should cache the pointer. This would mean we need to invalidate the pointer
    //if somthing with the array changes

    //return mpReference;
  }

std::string CArrayElementReference::getObjectDisplayName(bool regular, bool richtext) const
  {
    if (getObjectParent())
      {
        //if the array has as task as ancestor, use the task (skip the problem/method)
        CCopasiContainer* pT = getObjectAncestor("Task");

        std::string part;
        if (pT)
          part = pT->getObjectDisplayName(regular, richtext) + ".";
        else if (getObjectParent()->getObjectParent() && getObjectParent()->getObjectParent()->getObjectType() != "Model")
          part = getObjectParent()->getObjectParent()->getObjectDisplayName(regular, richtext) + ".";

        //now part contains the display name of the task, or the parent of the parent
        return part + getObjectParent()->getObjectName() + mIndex;
      }
    else
      return "Array" + mIndex;
  }

CCopasiObjectName CArrayElementReference::getCN() const
  {
    if (getObjectParent())
      return getObjectParent()->getCN() + mIndex;
    else
      return "Array" + mIndex;
  }

void CArrayElementReference::print(std::ostream * ostream) const
  {
    //if (mpReference)
    //  (*ostream) << *mpReference;

    //TODO perhaps we should cache the
    CCopasiAbstractArray::data_type * tmp = (double*)getValuePointer();
    if (tmp)
      (*ostream) << *tmp;
  };
