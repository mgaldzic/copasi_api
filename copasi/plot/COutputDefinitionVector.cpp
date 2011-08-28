/* Begin CVS Header
 $Source: /fs/turing/cvs/copasi_dev/copasi/plot/COutputDefinitionVector.cpp,v $
 $Revision: 1.4 $
 $Name: Build-33 $
 $Author: shoops $
 $Date: 2009/02/19 19:51:20 $
 End CVS Header */

// Copyright (C) 2008 by Pedro Mendes, Virginia Tech Intellectual
// Properties, Inc., EML Research, gGmbH, University of Heidelberg,
// and The University of Manchester.
// All rights reserved.

// Copyright (C) 2001 - 2007 by Pedro Mendes, Virginia Tech Intellectual
// Properties, Inc. and EML Research, gGmbH.
// All rights reserved.

#include "COutputDefinitionVector.h"
#include "report/CKeyFactory.h"
#include "copasi/report/CCopasiRootContainer.h"

COutputDefinitionVector::COutputDefinitionVector(const std::string & name,
    const CCopasiContainer * pParent):
    CCopasiVectorN< CPlotSpecification >(name, pParent),
    mKey(CCopasiRootContainer::getKeyFactory()->add("COutputDefinitionVector", this))
{}

COutputDefinitionVector::~COutputDefinitionVector()
{
  cleanup();
}

void COutputDefinitionVector::cleanup()
{
  CCopasiRootContainer::getKeyFactory()->remove(mKey);
}

const std::string& COutputDefinitionVector::getKey()
{
  return mKey;
}

CPlotSpecification* COutputDefinitionVector::createPlotSpec(const std::string & name,
    CPlotItem::Type type)
{
  unsigned C_INT32 i;
  for (i = 0; i < size(); i++)
    if ((*this)[i]->getObjectName() == name)
      return NULL; // duplicate name

  CPlotSpecification* pNewPlotSpec = new CPlotSpecification(name, this, type);
  pNewPlotSpec->setObjectName(name);

  add(pNewPlotSpec);
  return pNewPlotSpec;
}

bool COutputDefinitionVector::removePlotSpec(const std::string & key)
{
  CPlotSpecification* pPl =
    dynamic_cast<CPlotSpecification*>(CCopasiRootContainer::getKeyFactory()->get(key));
  unsigned C_INT32 index = this->CCopasiVector<CPlotSpecification>::getIndex(pPl);
  if (index == C_INVALID_INDEX)
    return false;

  this->CCopasiVector<CPlotSpecification>::remove(index);

  return true;
}
