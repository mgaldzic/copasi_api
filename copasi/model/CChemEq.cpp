/* Begin CVS Header
 $Source: /fs/turing/cvs/copasi_dev/copasi/model/CChemEq.cpp,v $
 $Revision: 1.50 $
 $Name: Build-33 $
 $Author: shoops $
 $Date: 2009/01/07 19:00:14 $
 End CVS Header */

// Copyright (C) 2008 by Pedro Mendes, Virginia Tech Intellectual
// Properties, Inc., EML Research, gGmbH, University of Heidelberg,
// and The University of Manchester.
// All rights reserved.

// Copyright (C) 2001 - 2007 by Pedro Mendes, Virginia Tech Intellectual
// Properties, Inc. and EML Research, gGmbH.
// All rights reserved.

// CChemEqElement
//
// A class describing an element of a chemical equation
// (C) Stefan Hoops 2001
//

#include "mathematics.h"
#include "copasi.h"
#include "CChemEq.h"
#include "utilities/CReadConfig.h"
#include "utilities/CCopasiVector.h"
#include "CMetabNameInterface.h"
#include "CCompartment.h"

CChemEq::CChemEq(const std::string & name,
                 const CCopasiContainer * pParent):
    CCopasiContainer(name, pParent, "Chemical Equation"),
    mReversible(false),
    mSubstrates("Substrates", this),
    mProducts("Products", this),
    mModifiers("Modifiers", this),
    mBalances("Balances", this)
{CONSTRUCTOR_TRACE;}

CChemEq::CChemEq(const CChemEq & src,
                 const CCopasiContainer * pParent):
    CCopasiContainer(src, pParent),
    mReversible(src.mReversible),
    mSubstrates(src.mSubstrates, this),
    mProducts(src.mProducts, this),
    mModifiers(src.mModifiers, this),
    mBalances(src.mBalances, this)
{CONSTRUCTOR_TRACE;}

CChemEq::~CChemEq(){cleanup(); DESTRUCTOR_TRACE;}

void CChemEq::cleanup()
{
  mSubstrates.cleanup();
  mProducts.cleanup();
  mModifiers.cleanup();
  mBalances.cleanup();
}

const CCopasiVector < CChemEqElement > & CChemEq::getSubstrates() const
  {return mSubstrates;}

const CCopasiVector < CChemEqElement > & CChemEq::getProducts() const
  {return mProducts;}

const CCopasiVector < CChemEqElement > & CChemEq::getModifiers() const
  {return mModifiers;}

const CCopasiVector < CChemEqElement > & CChemEq::getBalances() const
  {return mBalances;}

bool CChemEq::addMetabolite(const std::string & key, const C_FLOAT64 multiplicity, const MetaboliteRole & role)
{
  CChemEqElement element;
  element.setMetabolite(key);
  element.setMultiplicity(multiplicity);

  switch (role)
    {
    case CChemEq::SUBSTRATE:
      addElement(mSubstrates, element);
      addElement(mBalances, element, CChemEq::SUBSTRATE);
      break;
    case CChemEq::PRODUCT:
      addElement(mProducts, element);
      addElement(mBalances, element);
      break;
    case CChemEq::MODIFIER:
      addElement(mModifiers, element);
      break;
    default:
      fatalError();
      break;
    }

  return true;
}

unsigned C_INT32 CChemEq::getCompartmentNumber() const
  {
    unsigned C_INT32 i, imax = mBalances.size();
    unsigned C_INT32 j, jmax;
    unsigned C_INT32 Number;
    std::vector<const CCompartment *> Compartments;

    for (i = 0, Number = 0; i < imax; i++)
      {
        if (!mBalances[i]->getMetabolite())
          continue;

        for (j = 0, jmax = Compartments.size(); j < jmax; j++)
          if (Compartments[j] == mBalances[i]->getMetabolite()->getCompartment())
            break;

        if (j == jmax)
          {
            Number ++;
            Compartments.push_back(mBalances[i]->getMetabolite()->getCompartment());
          }
      }

    return Number;
  }

const CCompartment & CChemEq::getLargestCompartment() const
  {
    unsigned C_INT32 indexSubstrates = C_INVALID_INDEX;
    unsigned C_INT32 indexProducts = C_INVALID_INDEX;
    unsigned C_INT32 i, imax;

    C_FLOAT64 tmp, maxVol = -1.0;

    for (i = 0, imax = mSubstrates.size(); i < imax; i++)
      {
        if (!mSubstrates[i]->getMetabolite()) continue;

        tmp = mSubstrates[i]->getMetabolite()->getCompartment()->getValue();

        if (tmp > maxVol)
          {
            maxVol = tmp;
            indexSubstrates = i;
          }
      }

    for (i = 0, imax = mProducts.size(); i < imax; i++)
      {
        if (!mProducts[i]->getMetabolite()) continue;

        tmp = mProducts[i]->getMetabolite()->getCompartment()->getValue();

        if (tmp > maxVol)
          {
            maxVol = tmp;
            indexProducts = i;
          }
      }

    if (indexProducts != C_INVALID_INDEX)
      return *mProducts[indexProducts]->getMetabolite()->getCompartment();

    if (indexSubstrates != C_INVALID_INDEX)
      return *mSubstrates[indexSubstrates]->getMetabolite()->getCompartment();

    fatalError();

    return *mSubstrates[indexSubstrates]->getMetabolite()->getCompartment();
  }

void CChemEq::addElement(CCopasiVector < CChemEqElement > & structure,
                         const CChemEqElement & element,
                         CChemEq::MetaboliteRole role)
{
  unsigned C_INT32 i;

  std::string key = element.getMetaboliteKey();

  if (key == "")
    return; // don�t add empty element

  for (i = 0; i < structure.size(); i++)
    if (key == structure[i]->getMetaboliteKey())
      break;

  if (i >= structure.size())
    {
      CChemEqElement * Element = new CChemEqElement(element);

      if (role == CChemEq::SUBSTRATE)
        Element->setMultiplicity(- Element->getMultiplicity());

      structure.add(Element, true);
    }
  else if (role == CChemEq::SUBSTRATE)
    structure[i]->addToMultiplicity(- element.getMultiplicity());
  else
    structure[i]->addToMultiplicity(element.getMultiplicity());
}

C_INT32 CChemEq::getMolecularity(const MetaboliteRole role) const
  {
    const CCopasiVector<CChemEqElement> * tmpVector = NULL;

    switch (role)
      {
      case CChemEq::SUBSTRATE:
        tmpVector = &mSubstrates;
        break;
      case CChemEq::PRODUCT:
        tmpVector = &mProducts;
        break;
      case CChemEq::MODIFIER:
        tmpVector = &mModifiers;
        break;
      default:
        fatalError();
        break;
      }

    C_INT32 ccc, i, imax = tmpVector->size();
    ccc = 0;
    for (i = 0; i < imax; ++i)
      ccc += (C_INT32)floor((*tmpVector)[i]->getMultiplicity());

    return ccc;
  }

std::ostream & operator<<(std::ostream &os, const CChemEq & d)
{
  os << "CChemEq:" << std::endl;
  //os << "   mChemicalEquation:          " << d.getChemicalEquation() << std::endl;
  //os << "   mChemicalEquationConverted: " << d.getChemicalEquationConverted() << std::endl;

  os << "   mSubstrates:" << std::endl;
  os << d.mSubstrates;
  os << "   mProducts:" << std::endl;
  os << d.mProducts;
  os << "   mBalances:" << std::endl;
  os << d.mBalances;

  os << "----CChemEq" << std::endl;
  return os;
}
