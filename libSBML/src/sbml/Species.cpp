/**
 * @file    Species.cpp
 * @brief   Implementations of Species and ListOfSpecies.
 * @author  Ben Bornstein
 *
 * $Id: Species.cpp 11635 2010-08-03 04:05:22Z mhucka $
 * $HeadURL: https://sbml.svn.sourceforge.net/svnroot/sbml/trunk/libsbml/src/sbml/Species.cpp $
 *
 *<!---------------------------------------------------------------------------
 * This file is part of libSBML.  Please visit http://sbml.org for more
 * information about SBML, and the latest version of libSBML.
 *
 * Copyright 2005-2010 California Institute of Technology.
 * Copyright 2002-2005 California Institute of Technology and
 *                     Japan Science and Technology Corporation.
 * 
 * This library is free software; you can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation.  A copy of the license agreement is provided
 * in the file named "LICENSE.txt" included with this software distribution
 * and also available online as http://sbml.org/software/libsbml/license.html
 *----------------------------------------------------------------------- -->*/

#include <limits>

#include <sbml/xml/XMLNode.h>
#include <sbml/xml/XMLAttributes.h>
#include <sbml/xml/XMLInputStream.h>
#include <sbml/xml/XMLOutputStream.h>

#include <sbml/SBMLVisitor.h>
#include <sbml/SBMLError.h>
#include <sbml/SBMLDocument.h>
#include <sbml/Model.h>
#include <sbml/Species.h>

/** @cond doxygen-ignored */

using namespace std;

/** @endcond */

LIBSBML_CPP_NAMESPACE_BEGIN

Species::Species (unsigned int level, unsigned int version) :
   SBase ( level, version )
  , mId                       ( ""    )
  , mName                     ( ""    )
  , mInitialAmount            ( 0.0   )
  , mInitialConcentration     ( 0.0   )
  , mHasOnlySubstanceUnits    ( false )
  , mBoundaryCondition        ( false )
  , mCharge                   ( 0     )
  , mConstant                 ( false )
  , mIsSetInitialAmount       ( false )
  , mIsSetInitialConcentration( false )
  , mIsSetCharge              ( false )
  , mConversionFactor         ( ""    )
  , mIsSetBoundaryCondition   ( false )
  , mIsSetHasOnlySubstanceUnits (false )
  , mIsSetConstant             ( false )
{
  if (!hasValidLevelVersionNamespaceCombination())
    throw SBMLConstructorException();

  // if level 3 values have no defaults
  if (level == 3)
  {
    mInitialAmount = numeric_limits<double>::quiet_NaN();
    mInitialConcentration = numeric_limits<double>::quiet_NaN();
  }
}


Species::Species (SBMLNamespaces *sbmlns) :
    SBase                     ( sbmlns    )
  , mId                       ( ""    )
  , mName                     ( ""    )
  , mInitialAmount            ( 0.0   )
  , mInitialConcentration     ( 0.0   )
  , mHasOnlySubstanceUnits    ( false )
  , mBoundaryCondition        ( false )
  , mCharge                   ( 0     )
  , mConstant                 ( false )
  , mIsSetInitialAmount       ( false )
  , mIsSetInitialConcentration( false )
  , mIsSetCharge              ( false )
  , mConversionFactor         ( ""    )
  , mIsSetBoundaryCondition   ( false )
  , mIsSetHasOnlySubstanceUnits (false )
  , mIsSetConstant             ( false )
{
  if (!hasValidLevelVersionNamespaceCombination())
    throw SBMLConstructorException();

  // if level 3 values have no defaults
  if (sbmlns->getLevel() == 3)
  {
    mInitialAmount = numeric_limits<double>::quiet_NaN();
    mInitialConcentration = numeric_limits<double>::quiet_NaN();
  }
}


/** @cond doxygen-libsbml-internal */

/* constructor for validators */
Species::Species() :
  SBase()
{
}

/** @endcond */
                          

/*
 * Destroys this Species.
 */
Species::~Species ()
{
}


/*
 * Copy constructor. Creates a copy of this Species.
 */
Species::Species(const Species& orig) :
   SBase             ( orig                    )
 , mId               ( orig.mId                )  
 , mName             ( orig.mName              )
{
    this->mSpeciesType = orig.mSpeciesType;
    this->mCompartment = orig.mCompartment;

    this->mInitialAmount = orig.mInitialAmount;
    this->mInitialConcentration = orig.mInitialConcentration;

    this->mSubstanceUnits = orig.mSubstanceUnits;
    this->mSpatialSizeUnits = orig.mSpatialSizeUnits;

    this->mHasOnlySubstanceUnits = orig.mHasOnlySubstanceUnits;
    this->mBoundaryCondition = orig.mBoundaryCondition;
    this->mCharge = orig.mCharge;
    this->mConstant = orig.mConstant;

    this->mIsSetInitialAmount = orig.mIsSetInitialAmount;
    this->mIsSetInitialConcentration = orig.mIsSetInitialConcentration;
    this->mIsSetCharge = orig.mIsSetCharge;
    this->mConversionFactor         = orig.mConversionFactor;
    this->mIsSetBoundaryCondition   = orig.mIsSetBoundaryCondition;
    this->mIsSetHasOnlySubstanceUnits = orig.mIsSetHasOnlySubstanceUnits;
    this->mIsSetConstant             = orig.mIsSetConstant;
}


/*
 * Assignment operator.
 */
Species& Species::operator=(const Species& orig)
{
  if(&orig!=this)
  {
    this->SBase::operator =(orig);
    this->mId = orig.mId;
    this->mName = orig.mName;
    this->mSpeciesType = orig.mSpeciesType;
    this->mCompartment = orig.mCompartment;

    this->mInitialAmount = orig.mInitialAmount;
    this->mInitialConcentration = orig.mInitialConcentration;

    this->mSubstanceUnits = orig.mSubstanceUnits;
    this->mSpatialSizeUnits = orig.mSpatialSizeUnits;

    this->mHasOnlySubstanceUnits = orig.mHasOnlySubstanceUnits;
    this->mBoundaryCondition = orig.mBoundaryCondition;
    this->mCharge = orig.mCharge;
    this->mConstant = orig.mConstant;

    this->mIsSetInitialAmount = orig.mIsSetInitialAmount;
    this->mIsSetInitialConcentration = orig.mIsSetInitialConcentration;
    this->mIsSetCharge = orig.mIsSetCharge;
    this->mConversionFactor         = orig.mConversionFactor;
    this->mIsSetBoundaryCondition   = orig.mIsSetBoundaryCondition;
    this->mIsSetHasOnlySubstanceUnits = orig.mIsSetHasOnlySubstanceUnits;
    this->mIsSetConstant             = orig.mIsSetConstant;
  }

  return *this;
}


/*
 * Accepts the given SBMLVisitor.
 *
 * @return the result of calling <code>v.visit()</code>, which indicates
 * whether or not the Visitor would like to visit the Model's next
 * Species (if available).
 */
bool
Species::accept (SBMLVisitor& v) const
{
  return v.visit(*this);
}


/*
 * @return a (deep) copy of this Species.
 */
Species*
Species::clone () const
{
  return new Species(*this);
}


/*
 * Initializes the fields of this Species to their defaults:
 *
 *   - boundaryCondition     = false
 *   - constant              = false  (L2 only)
 *   - hasOnlySubstanceUnits = false  (L2 only)
 */
void
Species::initDefaults ()
{
  ////level 3 has no defaults
  //if (getLevel() < 3) 
  //{
    setBoundaryCondition     (false);
    setConstant              (false);
    setHasOnlySubstanceUnits (false);
  //}
}


/*
 * @return the id of this SBML object.
 */
const string&
Species::getId () const
{
  return mId;
}


/*
 * @return the name of this SBML object.
 */
const string&
Species::getName () const
{
  return (getLevel() == 1) ? mId : mName;
}


/*
 * @return the speciesType of this Species.
 */
const string&
Species::getSpeciesType () const
{
  return mSpeciesType;
}


/*
 * @return the compartment of this Species.
 */
const string&
Species::getCompartment () const
{
  return mCompartment;
}


/*
 * @return the initialAmount of this Species.
 */
double
Species::getInitialAmount () const
{
  double initialAmount = mInitialAmount;
  
  // need to cover case where user has changed level 
  // and expects an initial amount where there was none
  if ( getLevel() == 1 && isSetInitialConcentration() )
  {
    const Compartment *c = getModel()->getCompartment(getCompartment());
    if (c)
    {
      initialAmount = mInitialConcentration * c->getSize();
    }
  }

  return initialAmount;
}


/*
 * @return the initialConcentration of this Species.
 */
double
Species::getInitialConcentration () const
{
  return mInitialConcentration;
}


/*
 * @return the substanceUnits of this Species.
 */
const string&
Species::getSubstanceUnits () const
{
  return mSubstanceUnits;
}


/*
 * @return the spatialSizeUnits of this Species.
 */
const string&
Species::getSpatialSizeUnits () const
{
  return mSpatialSizeUnits;
}


/*
 * @return the units of this Species (L1 only).
 */
const string&
Species::getUnits () const
{
  return mSubstanceUnits;
}


/*
 * @return true if this Species hasOnlySubstanceUnits, false otherwise.
 */
bool
Species::getHasOnlySubstanceUnits () const
{
  return mHasOnlySubstanceUnits;
}


/*
 * @return true if this Species has boundaryCondition
 * true, false otherwise.
 */
bool
Species::getBoundaryCondition () const
{
  return mBoundaryCondition;
}


/*
 * @return the charge of this Species.
 */
int
Species::getCharge () const
{
  return mCharge;
}


/*
 * @return true if this Species is constant, false otherwise.
 */
bool
Species::getConstant () const
{
  return mConstant;
}


/*
 * @return the conversionFactor of this Species, as a string.
 */
const std::string& 
Species::getConversionFactor () const
{
  return mConversionFactor;
}


/*
 * @return true if the id of this SBML object has been set, false
 * otherwise.
 */
bool
Species::isSetId () const
{
  return (mId.empty() == false);
}


/*
 * @return true if the name of this SBML object has been set, false
 * otherwise.
 */
bool
Species::isSetName () const
{
  return (getLevel() == 1) ? (mId.empty() == false) : 
                            (mName.empty() == false);
}


/*
 * @return true if the speciesType of this Species has been set, false
 * otherwise.
 */
bool
Species::isSetSpeciesType () const
{
  return (mSpeciesType.empty() == false);
}


/*
 * @return true if the compartment of this Species has been set, false
 * otherwise.
 */
bool
Species::isSetCompartment () const
{
  return (mCompartment.empty() == false);
}


/*
 * @return true if the initialAmount of this Species has been set, false
 * otherwise.
 *
 * In SBML L1, a Species initialAmount is required and therefore <b>should
 * always be set</b>.  In L2, initialAmount is optional and as such may or
 * may not be set.
 */
bool
Species::isSetInitialAmount () const
{
  return mIsSetInitialAmount;
}


/*
 * @return true if the initialConcentration of this Species has been set,
 * false otherwise.
 */
bool
Species::isSetInitialConcentration () const
{
  return mIsSetInitialConcentration;
}


/*
 * @return true if the substanceUnits of this Species has been set, false
 * otherwise.
 */
bool
Species::isSetSubstanceUnits () const
{
  return (mSubstanceUnits.empty() == false);
}


/*
 * @return true if the spatialSizeUnits of this Species has been set, false
 * otherwise.
 */
bool
Species::isSetSpatialSizeUnits () const
{
  return (mSpatialSizeUnits.empty() == false);
}


/*
 * @return true if the units of this Species has been set, false otherwise
 * (L1 only).
 */
bool
Species::isSetUnits () const
{
  return isSetSubstanceUnits();
}


/*
 * @return true if the charge of this Species has been set, false
 * otherwise.
 */
bool
Species::isSetCharge () const
{
  return mIsSetCharge;
}


/*
 * @return true if the conversionFactor of this Species has been set, false
 * otherwise.
 */
bool
Species::isSetConversionFactor () const
{
  return (mConversionFactor.empty() == false);
}


/*
 * Predicate returning @c true or @c false depending on whether this
 * Species's "boundaryCondition" attribute has been set.
 */
bool 
Species::isSetBoundaryCondition () const
{
  return mIsSetBoundaryCondition;
}


/*
 * Predicate returning @c true or @c false depending on whether this
 * Species's "hasOnlySubstanceUnits" attribute has been set.
 */
bool 
Species::isSetHasOnlySubstanceUnits () const
{
  return mIsSetHasOnlySubstanceUnits;
}


/*
 * Predicate returning @c true or @c false depending on whether this
 * Species's "constant" attribute has been set.
 */
bool 
Species::isSetConstant () const
{
  return mIsSetConstant;
}


/*
 * Sets the id of this SBML object to a copy of sid.
 */
int
Species::setId (const std::string& sid)
{
  /* since the setId function has been used as an
   * alias for setName we cant require it to only
   * be used on a L2 model
   */
/*  if (getLevel() == 1)
  {
    return LIBSBML_UNEXPECTED_ATTRIBUTE;
  }
*/
  if (!(SyntaxChecker::isValidSBMLSId(sid)))
  {
    return LIBSBML_INVALID_ATTRIBUTE_VALUE;
  }
  else
  {
    mId = sid;
    return LIBSBML_OPERATION_SUCCESS;
  }
}


/*
 * Sets the name of this SBML object to a copy of name.
 */
int
Species::setName (const std::string& name)
{
  /* if this is setting an L2 name the type is string
   * whereas if it is setting an L1 name its type is SId
   */
  if (getLevel() == 1)
  {
    if (!(SyntaxChecker::isValidSBMLSId(name)))
    {
      return LIBSBML_INVALID_ATTRIBUTE_VALUE;
    }
    else
    {
      mId = name;
      return LIBSBML_OPERATION_SUCCESS;
    }
  }
  else
  {
    mName = name;
    return LIBSBML_OPERATION_SUCCESS;
  }
}


/*
 * Sets the speciesType field of this Species to a copy of sid.
 */
int
Species::setSpeciesType (const std::string& sid)
{
  if ( (getLevel() < 2)
    || (getLevel() == 2 && getVersion() == 1))
  {
    return LIBSBML_UNEXPECTED_ATTRIBUTE;
  }
  else if (!(SyntaxChecker::isValidSBMLSId(sid)))
  {
    return LIBSBML_INVALID_ATTRIBUTE_VALUE;
  }
  else
  {
    mSpeciesType = sid;
    return LIBSBML_OPERATION_SUCCESS;
  }
}


/*
 * Sets the compartment of this Species to a copy of sid.
 */
int
Species::setCompartment (const std::string& sid)
{
  if (!(SyntaxChecker::isValidSBMLSId(sid)))
  {
    return LIBSBML_INVALID_ATTRIBUTE_VALUE;
  }
  else
  {
    mCompartment = sid;
    return LIBSBML_OPERATION_SUCCESS;
  }
}


/*
 * Sets the initialAmount of this Species to value and marks the field as
 * set.  This method also unsets the initialConcentration field.
 */
int
Species::setInitialAmount (double value)
{
  mInitialAmount      = value;
  mIsSetInitialAmount = true;

  unsetInitialConcentration();
  return LIBSBML_OPERATION_SUCCESS;
}


/*
 * Sets the initialConcentration of this Species to value and marks the
 * field as set.  This method also unsets the initialAmount field.
 */
int
Species::setInitialConcentration (double value)
{
  if ( getLevel() < 2)
  {
    return LIBSBML_UNEXPECTED_ATTRIBUTE;
  }
  else
  {
    mInitialConcentration      = value;
    mIsSetInitialConcentration = true;

    unsetInitialAmount();
    return LIBSBML_OPERATION_SUCCESS;
  }
}


/*
 * Sets the substanceUnits of this Species to a copy of sid.
 */
int
Species::setSubstanceUnits (const std::string& sid)
{
  if (!(SyntaxChecker::isValidSBMLSId(sid)))
  {
    return LIBSBML_INVALID_ATTRIBUTE_VALUE;
  }
  else
  {
    mSubstanceUnits = sid;
    return LIBSBML_OPERATION_SUCCESS;
  }
}


/*
 * Sets the spatialSizeUnits of this Species to a copy of sid.
 */
int
Species::setSpatialSizeUnits (const std::string& sid)
{
  if ( (getLevel() != 2)
    || (getLevel() == 2 && getVersion() > 2))
  {
    return LIBSBML_UNEXPECTED_ATTRIBUTE;
  }
  else if (!(SyntaxChecker::isValidSBMLSId(sid)))
  {
    return LIBSBML_INVALID_ATTRIBUTE_VALUE;
  }
  else
  {
    mSpatialSizeUnits = sid;
    return LIBSBML_OPERATION_SUCCESS;
  }
}


/*
 * Sets the units of this Species to a copy of sname (L1 only).
 */
int
Species::setUnits (const std::string& sname)
{
  return setSubstanceUnits(sname);
}


/*
 * Sets the hasOnlySubstanceUnits field of this Species to value.
 */
int
Species::setHasOnlySubstanceUnits (bool value)
{
  if (getLevel() < 2)
  {
    return LIBSBML_UNEXPECTED_ATTRIBUTE;
  }
  else
  {
    mHasOnlySubstanceUnits = value;
    mIsSetHasOnlySubstanceUnits = true;
    return LIBSBML_OPERATION_SUCCESS;
  }
}


/*
 * Sets the boundaryCondition of this Species to value.
 */
int
Species::setBoundaryCondition (bool value)
{
  mBoundaryCondition = value;
  mIsSetBoundaryCondition = true;
  return LIBSBML_OPERATION_SUCCESS;
}


/*
 * Sets the charge of this Species to value and marks the field as set.
 */
int
Species::setCharge (int value)
{
  if ( !((getLevel() == 1)
    || (getLevel() == 2 && getVersion() == 1)))
  {
    return LIBSBML_UNEXPECTED_ATTRIBUTE;
  }
  else
  {
    mCharge      = value;
    mIsSetCharge = true;
    return LIBSBML_OPERATION_SUCCESS;
  }
}


/*
 * Sets the constant field of this Species to value.
 */
int
Species::setConstant (bool value)
{
  if ( getLevel() < 2 )
  {
    mConstant = value;
    return LIBSBML_UNEXPECTED_ATTRIBUTE;
  }
  else
  {
    mConstant = value;
    mIsSetConstant = true;
    return LIBSBML_OPERATION_SUCCESS;
  }
}


/*
 * Sets the conversionFactor field of this Species to a copy of sid.
 */
int
Species::setConversionFactor (const std::string& sid)
{
  if (getLevel() < 3)
  {
    return LIBSBML_UNEXPECTED_ATTRIBUTE;
  }
  else if (!(SyntaxChecker::isValidSBMLSId(sid)))
  {
    return LIBSBML_INVALID_ATTRIBUTE_VALUE;
  }
  else
  {
    mConversionFactor = sid;
    return LIBSBML_OPERATION_SUCCESS;
  }
}


/*
 * Unsets the name of this SBML object.
 */
int
Species::unsetName ()
{
  if (getLevel() == 1) 
  {
    mId.erase();
  }
  else 
  {
    mName.erase();
  }

  if (getLevel() == 1 && mId.empty())
  {
    return LIBSBML_OPERATION_SUCCESS;
  }
  else if (mName.empty())
  {
    return LIBSBML_OPERATION_SUCCESS;
  }
  else
  {
    return LIBSBML_OPERATION_FAILED;
  }
}


/*
 * Unsets the speciesType of this Species.
 */
int
Species::unsetSpeciesType ()
{
  mSpeciesType.erase();

  if (mSpeciesType.empty()) 
  {
    return LIBSBML_OPERATION_SUCCESS;
  }
  else
  {
    return LIBSBML_OPERATION_FAILED;
  }
}


/*
 * Unsets the initialAmount of this Species.
 */
int
Species::unsetInitialAmount ()
{
  mInitialAmount      = numeric_limits<double>::quiet_NaN();
  mIsSetInitialAmount = false;
  
  if (!isSetInitialAmount())
  {
    return LIBSBML_OPERATION_SUCCESS;
  }
  else
  {
    return LIBSBML_OPERATION_FAILED;
  }
}


/*
 * Unsets the initialConcentration of this Species.
 */
int
Species::unsetInitialConcentration ()
{
  mInitialConcentration      = numeric_limits<double>::quiet_NaN();
  mIsSetInitialConcentration = false;

  if (!isSetInitialConcentration())
  {
    return LIBSBML_OPERATION_SUCCESS;
  }
  else
  {
    return LIBSBML_OPERATION_FAILED;
  }
}


/*
 * Unsets the substanceUnits of this Species.
 */
int
Species::unsetSubstanceUnits ()
{
  mSubstanceUnits.erase();
  
  if (mSubstanceUnits.empty()) 
  {
    return LIBSBML_OPERATION_SUCCESS;
  }
  else
  {
    return LIBSBML_OPERATION_FAILED;
  }
}


/*
 * Unsets the spatialSizeUnits of this Species.
 */
int
Species::unsetSpatialSizeUnits ()
{
  mSpatialSizeUnits.erase();

  if (mSpatialSizeUnits.empty()) 
  {
    return LIBSBML_OPERATION_SUCCESS;
  }
  else
  {
    return LIBSBML_OPERATION_FAILED;
  }
}


/*
 * Unsets the units of this Species (L1 only).
 */
int
Species::unsetUnits ()
{
  return unsetSubstanceUnits();
}


/*
 * Unsets the charge of this Species.
 */
int
Species::unsetCharge ()
{
  if ( !((getLevel() == 1)
    || (getLevel() == 2 && getVersion() == 1)))
  {
    return LIBSBML_UNEXPECTED_ATTRIBUTE;
  }

  mCharge      = 0;
  mIsSetCharge = false;

  if (!isSetCharge())
  {
    return LIBSBML_OPERATION_SUCCESS;
  }
  else
  {
    return LIBSBML_OPERATION_FAILED;
  }    
}


/*
 * Unsets the conversionFactor of this Species.
 */
int
Species::unsetConversionFactor ()
{
  mConversionFactor.erase();

  if (mConversionFactor.empty()) 
  {
    return LIBSBML_OPERATION_SUCCESS;
  }
  else
  {
    return LIBSBML_OPERATION_FAILED;
  }
}


/*
  * Constructs and returns a UnitDefinition that expresses the units of this 
  * Species.
  */
UnitDefinition *
Species::getDerivedUnitDefinition()
{
  /* if we have the whole model but it is not in a document
   * it is still possible to determine the units
   */
  Model * m = static_cast <Model *> (getAncestorOfType(SBML_MODEL));

  if (m)
  {
    if (!m->isPopulatedListFormulaUnitsData())
    {
      m->populateListFormulaUnitsData();
    }
    
    if (m->getFormulaUnitsData(getId(), getTypeCode()))
    {
      return m->getFormulaUnitsData(getId(), getTypeCode())
                                             ->getUnitDefinition();
    }
    else
    {
      return NULL;
    }  
  }
  else
  {
    return NULL;
  }
}


/*
  * Constructs and returns a UnitDefinition that expresses the units of this 
  * Compartment.
  */
const UnitDefinition *
Species::getDerivedUnitDefinition() const
{
  return const_cast <Species *> (this)->getDerivedUnitDefinition();
}


/*
 * @return the SBMLTypeCode_t of this SBML object or SBML_UNKNOWN
 * (default).
 *
 * @see getElementName()
 */
SBMLTypeCode_t
Species::getTypeCode () const
{
  return SBML_SPECIES;
}


/*
 * @return the name of this element ie "specie" (L1) or "species" (L2).
 */
const string&
Species::getElementName () const
{
  static const string specie  = "specie";
  static const string species = "species";

  return (getLevel() == 1 && getVersion() == 1) ? specie : species;
}


bool 
Species::hasRequiredAttributes() const
{
  bool allPresent = true;

  /* required attributes for species: 
   * @li id (name L1)
   * @li compartment
   * @li initialAmount (L1 only)
   * @li hasOnlySubstanceUnits (L3 on)
   * @li boundaryCondition (L3 on)
   * @li constant (L3 on)
   */

  if (!isSetId())
    allPresent = false;

  if (!isSetCompartment())
    allPresent = false;

  if (getLevel() == 1 && !isSetInitialAmount())
    allPresent = false;

  if (getLevel() > 2 && !isSetHasOnlySubstanceUnits())
    allPresent = false;

  if (getLevel() > 2 && !isSetBoundaryCondition())
    allPresent = false;

  if (getLevel() > 2 && !isSetConstant())
    allPresent = false;

  return allPresent;
}


/** @cond doxygen-libsbml-internal */
/*
 * Subclasses should override this method to read values from the given
 * XMLAttributes set into their specific fields.  Be sure to call your
 * parents implementation of this method as well.
 */
void
Species::readAttributes (const XMLAttributes& attributes)
{
  SBase::readAttributes(attributes);

  const unsigned int level   = getLevel  ();
  switch (level)
  {
  case 1:
    readL1Attributes(attributes);
    break;
  case 2:
    readL2Attributes(attributes);
    break;
  case 3:
  default:
    readL3Attributes(attributes);
    break;
  }
}
/** @endcond */


/** @cond doxygen-libsbml-internal */
/*
 * Subclasses should override this method to read values from the given
 * XMLAttributes set into their specific fields.  Be sure to call your
 * parents implementation of this method as well.
 */
void
Species::readL1Attributes (const XMLAttributes& attributes)
{
  SBase::readAttributes(attributes);

  const unsigned int level   = getLevel  ();
  const unsigned int version = getVersion();
 
  std::vector<std::string> expectedAttributes;
  expectedAttributes.clear();
  expectedAttributes.push_back("name");
  expectedAttributes.push_back("compartment");
  expectedAttributes.push_back("initialAmount");
  expectedAttributes.push_back("boundaryCondition");
  expectedAttributes.push_back("charge");
  expectedAttributes.push_back("units");

  // check that all attributes are expected
  for (int i = 0; i < attributes.getLength(); i++)
  {
    std::vector<std::string>::const_iterator end = expectedAttributes.end();
    std::vector<std::string>::const_iterator begin = expectedAttributes.begin();
    std::string name = attributes.getName(i);
    std::string prefix = attributes.getPrefix(i);
    // only check attributes in the sbml namespace   
    if (prefix.empty() || prefix == "sbml")
    {
      if (std::find(begin, end, name) == end)
      {
        logUnknownAttribute(name, level, version, "<species>");
      }
    }
  }


  //
  // name: SName   { use="required" }  (L1v1, L1v2)
  //
  bool assigned = attributes.readInto("name", mId, getErrorLog(), true);
  if (assigned && mId.size() == 0)
  {
    logEmptyString("name", level, version, "<species>");
  }
  if (!SyntaxChecker::isValidSBMLSId(mId)) logError(InvalidIdSyntax);

  //
  // compartment: SName  { use="required" }  (L1v1, L2v1)
  //
  attributes.readInto("compartment", mCompartment, getErrorLog(), true);

  //
  // initialAmount: double  { use="required" }  (L1v1, L1v2)
  //
  mIsSetInitialAmount = attributes.readInto("initialAmount", mInitialAmount,
                                                  getErrorLog(), true);

  //
  //          units: SName  { use="optional" }  (L1v1, L1v2)
  //
  assigned = attributes.readInto("units", mSubstanceUnits);
  if (assigned && mSubstanceUnits.size() == 0)
  {
    logEmptyString("units", level, version, "<species>");
  }
  if (!SyntaxChecker::isValidUnitSId(mSubstanceUnits))
  {
    logError(InvalidUnitIdSyntax);
  }

  //
  // boundaryCondition: boolean
  // { use="optional" default="false" }  (L1v1, L1v2, L2v1->)
  //
  attributes.readInto("boundaryCondition", mBoundaryCondition);

  //
  // charge: integer  { use="optional" }  (L1v1, L1v2, L2v1)
  //
  mIsSetCharge = attributes.readInto("charge", mCharge);
}
/** @endcond */


/** @cond doxygen-libsbml-internal */
/*
 * Subclasses should override this method to read values from the given
 * XMLAttributes set into their specific fields.  Be sure to call your
 * parents implementation of this method as well.
 */
void
Species::readL2Attributes (const XMLAttributes& attributes)
{
  SBase::readAttributes(attributes);

  const unsigned int level   = getLevel  ();
  const unsigned int version = getVersion();

  std::vector<std::string> expectedAttributes;
  expectedAttributes.clear();
  expectedAttributes.push_back("name");
  expectedAttributes.push_back("compartment");
  expectedAttributes.push_back("initialAmount");
  expectedAttributes.push_back("boundaryCondition");
  expectedAttributes.push_back("charge");
  expectedAttributes.push_back("metaid");
  expectedAttributes.push_back("id");
  expectedAttributes.push_back("initialConcentration");
  expectedAttributes.push_back("substanceUnits");
  expectedAttributes.push_back("hasOnlySubstanceUnits");
  expectedAttributes.push_back("constant");

  if (version > 1)
  {
    expectedAttributes.push_back("speciesType");
  }

  if (version < 3)
  {
    expectedAttributes.push_back("spatialSizeUnits");
  }

  if (version > 2)
  {
    expectedAttributes.push_back("sboTerm");
  }

  // check that all attributes are expected
  for (int i = 0; i < attributes.getLength(); i++)
  {
    std::vector<std::string>::const_iterator end = expectedAttributes.end();
    std::vector<std::string>::const_iterator begin = expectedAttributes.begin();
    std::string name = attributes.getName(i);
    std::string prefix = attributes.getPrefix(i);
    // only check attributes in the sbml namespace   
    if (prefix.empty() || prefix == "sbml")
    {
      if (std::find(begin, end, name) == end)
      {
        logUnknownAttribute(name, level, version, "<species>");
      }
    }
  }


  //
  //   id: SId     { use="required" }  (L2v1->)
  //
  bool assigned = attributes.readInto("id", mId, getErrorLog(), true);
  if (assigned && mId.size() == 0)
  {
    logEmptyString("id", level, version, "<species>");
  }
  if (!SyntaxChecker::isValidSBMLSId(mId)) logError(InvalidIdSyntax);

  //
  // compartment: SId    { use="required" }  (L2v1->)
  //
  attributes.readInto("compartment", mCompartment, getErrorLog(), true);

  //
  // initialAmount: double  { use="optional" }  (L2v1->)
  //
  mIsSetInitialAmount = attributes.readInto("initialAmount", mInitialAmount);

  //
  // substanceUntis: SId    { use="optional" }  (L2v1->)
  //
  assigned = attributes.readInto("substanceUnits", mSubstanceUnits);
  if (assigned && mSubstanceUnits.size() == 0)
  {
    logEmptyString("substanceUnits", level, version, "<species>");
  }
  if (!SyntaxChecker::isValidUnitSId(mSubstanceUnits))
  {
    logError(InvalidUnitIdSyntax);
  }

  //
  // boundaryCondition: boolean
  // { use="optional" default="false" }  (L1v1, L1v2, L2v1->)
  //
  attributes.readInto("boundaryCondition", mBoundaryCondition);

  //
  // charge: integer  { use="optional" }  (L1v1, L1v2, L2v1)
  // charge: integer  { use="optional" }  deprecated (L2v2)
  //
  mIsSetCharge = attributes.readInto("charge", mCharge);

  //
  // name: string  { use="optional" }  (L2v1->)
  //
  attributes.readInto("name", mName);

  //
  // speciesType: SId  { use="optional" }  (L2v2-> L2v4)
  //
  if (version > 1)
  {
    attributes.readInto("speciesType", mSpeciesType);
  }

  //
  // initialConcentration: double  { use="optional" }  (L2v1->)
  //
  mIsSetInitialConcentration =
      attributes.readInto("initialConcentration", mInitialConcentration);

  //
  // spatialSizeUnits: SId  { use="optional" }  (L2v1, L2v2) removed in l2v3
  //
  if (version < 3)
  {
    assigned = attributes.readInto("spatialSizeUnits", mSpatialSizeUnits);
    if (assigned && mSpatialSizeUnits.size() == 0)
    {
      logEmptyString("spatialSizeUnits", level, version, "<species>");
    }
    if (!SyntaxChecker::isValidUnitSId(mSpatialSizeUnits))
    {
      logError(InvalidUnitIdSyntax);
    }
  }

  //
  // hasOnlySubstanceUnits: boolean
  // { use="optional" default="false" }  (L2v1->)
  //
  attributes.readInto("hasOnlySubstanceUnits", mHasOnlySubstanceUnits);

  //
  // constant: boolean  { use="optional" default="false" }  (L2v1->)
  //
  attributes.readInto("constant", mConstant);

  //
  // sboTerm: SBOTerm { use="optional" }  (L2v3->)
  //
  if (version > 2) 
    mSBOTerm = SBO::readTerm(attributes, this->getErrorLog(), level, version);
}
/** @endcond */


/** @cond doxygen-libsbml-internal */
/*
 * Subclasses should override this method to read values from the given
 * XMLAttributes set into their specific fields.  Be sure to call your
 * parents implementation of this method as well.
 */
void
Species::readL3Attributes (const XMLAttributes& attributes)
{
  SBase::readAttributes(attributes);

  const unsigned int level   = getLevel  ();
  const unsigned int version = getVersion();
  std::vector<std::string> expectedAttributes;
  expectedAttributes.clear();
  expectedAttributes.push_back("name");
  expectedAttributes.push_back("compartment");
  expectedAttributes.push_back("initialAmount");
  expectedAttributes.push_back("boundaryCondition");
  expectedAttributes.push_back("metaid");
  expectedAttributes.push_back("id");
  expectedAttributes.push_back("initialConcentration");
  expectedAttributes.push_back("substanceUnits");
  expectedAttributes.push_back("hasOnlySubstanceUnits");
  expectedAttributes.push_back("constant");
  expectedAttributes.push_back("sboTerm");
  expectedAttributes.push_back("conversionFactor");

  // check that all attributes are expected
  for (int i = 0; i < attributes.getLength(); i++)
  {
    std::vector<std::string>::const_iterator end = expectedAttributes.end();
    std::vector<std::string>::const_iterator begin = expectedAttributes.begin();
    std::string name = attributes.getName(i);
    std::string prefix = attributes.getPrefix(i);
    // only check attributes in the sbml namespace   
    if (prefix.empty() || prefix == "sbml")
    {
      if (std::find(begin, end, name) == end)
      {
        logUnknownAttribute(name, level, version, "<species>");
      }
    }
  }


  //
  //   id: SId     { use="required" }  (L2v1->)
  //
  bool assigned = attributes.readInto("id", mId, getErrorLog());
  if (!assigned)
  {
    getErrorLog()->logError(AllowedAttributesOnSpecies, level, version);
  }
  if (assigned && mId.size() == 0)
  {
    logEmptyString("id", level, version, "<species>");
  }
  if (!SyntaxChecker::isValidSBMLSId(mId)) logError(InvalidIdSyntax);

  //
  // compartment: SId    { use="required" }  (L2v1->)
  //
  assigned = attributes.readInto("compartment", mCompartment, getErrorLog());
  if (!assigned)
  {
    logError(MissingSpeciesCompartment, level, version);
  }

  //
  // initialAmount: double  { use="optional" }  (L2v1->)
  //
  mIsSetInitialAmount = attributes.readInto("initialAmount", mInitialAmount);

  //
  // substanceUntis: SId    { use="optional" }  (L2v1->)
  //
  const string units = (level == 1) ? "units" : "substanceUnits";
  assigned = attributes.readInto(units, mSubstanceUnits);
  if (assigned && mSubstanceUnits.size() == 0)
  {
    logEmptyString("substanceUnits", level, version, "<species>");
  }
  if (!SyntaxChecker::isValidUnitSId(mSubstanceUnits))
  {
    logError(InvalidUnitIdSyntax);
  }

  //
  // boundaryCondition: boolean
  // { use="required" }  (L3v1->)
  //
  mIsSetBoundaryCondition = attributes.readInto("boundaryCondition", 
                               mBoundaryCondition, getErrorLog());
  if (!mIsSetBoundaryCondition)
  {
    logError(AllowedAttributesOnSpecies, level, version);
  }

  //
  // name: string  { use="optional" }  (L2v1->)
  //
  attributes.readInto("name", mName);

  //
  // initialConcentration: double  { use="optional" }  (L2v1->)
  //
  mIsSetInitialConcentration =
        attributes.readInto("initialConcentration", mInitialConcentration);

  //
  // hasOnlySubstanceUnits: boolean
  // { use="required" } (L3v1 -> )
  mIsSetHasOnlySubstanceUnits = attributes.readInto(
                         "hasOnlySubstanceUnits", mHasOnlySubstanceUnits,
                          getErrorLog());
  if (!mIsSetHasOnlySubstanceUnits)
  {
    logError(AllowedAttributesOnSpecies, level, version);
  }

  //
  // constant: boolean  { use="required" }  (L3v1->)
  //
  mIsSetConstant = attributes.readInto("constant", mConstant, getErrorLog());
  if (!mIsSetConstant)
  {
    logError(AllowedAttributesOnSpecies, level, version);
  }

  mSBOTerm = SBO::readTerm(attributes, this->getErrorLog(), level, version);

  //
  // conversionFactor: SIdRef {use="optional" } (L3v1 ->)
  //
  assigned = attributes.readInto("conversionFactor", mConversionFactor);
  if (assigned && mConversionFactor.size() == 0)
  {
    logEmptyString("conversionFactor", level, version, "<species>");
  }
  if (!SyntaxChecker::isValidSBMLSId(mConversionFactor))
  {
    logError(InvalidIdSyntax);
  }
}
/** @endcond */


/** @cond doxygen-libsbml-internal */
/*
 * Subclasses should override this method to write their XML attributes
 * to the XMLOutputStream.  Be sure to call your parents implementation
 * of this method as well.
 */
void
Species::writeAttributes (XMLOutputStream& stream) const
{
  SBase::writeAttributes(stream);

  const unsigned int level   = getLevel  ();
  const unsigned int version = getVersion();

  //
  // name: SName   { use="required" }  (L1v1, L1v2)
  //   id: SId     { use="required" }  (L2v1, L2v2)
  //
  const string id = (level == 1) ? "name" : "id";
  stream.writeAttribute(id, mId);

  if (level > 1)
  {
    //
    // name: string  { use="optional" }  (L2v1->)
    //
    stream.writeAttribute("name", mName);

    //
    // speciesType: SId  { use="optional" }  (L2v2->)
    //
    if (level == 2 && version > 1)
    {
      stream.writeAttribute("speciesType", mSpeciesType);
    }
  }

  //
  // compartment: SName  { use="required" }  (L1v1, L2v1)
  // compartment: SId    { use="required" }  (L2v1->)
  //
  stream.writeAttribute("compartment", mCompartment);

  //
  // initialAmount: double  { use="required" }  (L1v1, L1v2)
  // initialAmount: double  { use="optional" }  (L2v1->)
  //
  if ( isSetInitialAmount() )
  {
    stream.writeAttribute("initialAmount", mInitialAmount);
  }

  //
  // initialConcentration: double  { use="optional" }  (L2v1-> )
  //
  else if ( level > 1 && isSetInitialConcentration() )
  {
    stream.writeAttribute("initialConcentration", mInitialConcentration);
  }

  //
  // If user is converting a model from L2 to L1, two possiblities exist
  // that have not been covered:
  //
  //   1.  InitialConcentration has been used, so it must be converted to
  //       initialAmount
  //
  //   2.  No initialAmount/initialAmount has been set, but initialAmount
  //       is required in L1
  //
  else if (level == 1)
  {
    if ( isSetInitialConcentration() )
    {
      const Model*       m = getModel();
      const Compartment* c = m ? m->getCompartment( getCompartment() ) : 0;

      if (c)
      {
        double amount = mInitialConcentration * c->getSize();
        stream.writeAttribute("initialAmount", amount);
      }
      else
      {
        stream.writeAttribute("initialAmount", mInitialConcentration);
      }
    }
    else
    {
      stream.writeAttribute("initialAmount", mInitialAmount);
    }
  }

  //
  //          units: SName  { use="optional" }  (L1v1, L1v2)
  // substanceUntis: SId    { use="optional" }  (L2v1->)
  //
  const string units = (level == 1) ? "units" : "substanceUnits";
  stream.writeAttribute( units, getUnits() );

  if (level > 1)
  {
    if (level == 2 && version < 3)
    {
      //
      // spatialSizeUnits: SId  { use="optional" }  (L2v1, L2v2)
      //
      stream.writeAttribute("spatialSizeUnits", mSpatialSizeUnits);
    }

    //
    // hasOnlySubstanceUnits: boolean
    // { use="optional" default="false" }  (L2v1->)
    // { use="required" }  (L3v1->)
    //
    if (level == 2 && mHasOnlySubstanceUnits)
    {
      stream.writeAttribute( "hasOnlySubstanceUnits", 
                              mHasOnlySubstanceUnits );
    }
    else if (level > 2)
    {
      stream.writeAttribute( "hasOnlySubstanceUnits", 
                              mHasOnlySubstanceUnits );
    }
  }

  //
  // boundaryCondition: boolean
  // { use="optional" default="false" }  (L1v1, L1v2, L2v1->)
  // { use="required" }  (L3v1->)
  //
  if (level < 3)
  {
    if (mBoundaryCondition)
      stream.writeAttribute("boundaryCondition", mBoundaryCondition);
  }
  else
  {
    stream.writeAttribute("boundaryCondition", mBoundaryCondition);
  }


  //
  // charge: integer  { use="optional" }  (L1v1, L1v2, L2v1)
  // charge: integer  { use="optional" }  deprecated (L2v2)
  //
  if ( isSetCharge() )
  {
    stream.writeAttribute("charge", mCharge);
  }

  if (level > 1)
  {
    //
    // constant: boolean  { use="optional" default="false" }  (L2v1->)
    // constant: boolean  { use="required" }  (L3v1->)
    //
    if (level == 2 && mConstant != false)
    {
      stream.writeAttribute("constant", mConstant);
    }
    else if (level > 2)
    {
      stream.writeAttribute("constant", mConstant);
    }

    //
    // sboTerm: SBOTerm { use="optional" }  (L2v3->)
    //
    if (!(level == 2 && version < 3)) 
      SBO::writeTerm(stream, mSBOTerm);
  }

  if (level > 2)
  {
    stream.writeAttribute("conversionFactor", mConversionFactor);
  }
}
/** @endcond */



/*
 * @return a (deep) copy of this ListOfSpecies.
 */
ListOfSpecies*
ListOfSpecies::clone () const
{
  return new ListOfSpecies(*this);
}


/*
 * @return the SBMLTypeCode_t of SBML objects contained in this ListOf or
 * SBML_UNKNOWN (default).
 */
SBMLTypeCode_t
ListOfSpecies::getItemTypeCode () const
{
  return SBML_SPECIES;
}


/*
 * @return the name of this element ie "listOfSpecies".
 */
const string&
ListOfSpecies::getElementName () const
{
  static const string name = "listOfSpecies";
  return name;
}


/* return nth item in list */
Species *
ListOfSpecies::get(unsigned int n)
{
  return static_cast<Species*>(ListOf::get(n));
}


/* return nth item in list */
const Species *
ListOfSpecies::get(unsigned int n) const
{
  return static_cast<const Species*>(ListOf::get(n));
}


/**
 * Used by ListOf::get() to lookup an SBase based by its id.
 */
struct IdEqS : public unary_function<SBase*, bool>
{
  const string& id;

  IdEqS (const string& id) : id(id) { }
  bool operator() (SBase* sb) 
       { return static_cast <Species *> (sb)->getId() == id; }
};


/* return item by id */
Species*
ListOfSpecies::get (const std::string& sid)
{
  return const_cast<Species*>( 
    static_cast<const ListOfSpecies&>(*this).get(sid) );
}


/* return item by id */
const Species*
ListOfSpecies::get (const std::string& sid) const
{
  vector<SBase*>::const_iterator result;

  result = find_if( mItems.begin(), mItems.end(), IdEqS(sid) );
  return (result == mItems.end()) ? 0 : static_cast <Species*> (*result);
}


/* Removes the nth item from this list */
Species*
ListOfSpecies::remove (unsigned int n)
{
   return static_cast<Species*>(ListOf::remove(n));
}


/* Removes item in this list by id */
Species*
ListOfSpecies::remove (const std::string& sid)
{
  SBase* item = 0;
  vector<SBase*>::iterator result;

  result = find_if( mItems.begin(), mItems.end(), IdEqS(sid) );

  if (result != mItems.end())
  {
    item = *result;
    mItems.erase(result);
  }

  return static_cast <Species*> (item);
}


/** @cond doxygen-libsbml-internal */
/*
 * @return the ordinal position of the element with respect to its siblings
 * or -1 (default) to indicate the position is not significant.
 */
int
ListOfSpecies::getElementPosition () const
{
  return 6;
}
/** @endcond */


/** @cond doxygen-libsbml-internal */
/*
 * @return the SBML object corresponding to next XMLToken in the
 * XMLInputStream or NULL if the token was not recognized.
 */
SBase*
ListOfSpecies::createObject (XMLInputStream& stream)
{
  const string& name   = stream.peek().getName();
  SBase*        object = 0;


  if (name == "species" || name == "specie")
  {
    try
    {
      object = new Species(getSBMLNamespaces());
    }
    catch (SBMLConstructorException*)
    {
      object = new Species(SBMLDocument::getDefaultLevel(),
        SBMLDocument::getDefaultVersion());
    }
    catch ( ... )
    {
      object = new Species(SBMLDocument::getDefaultLevel(),
        SBMLDocument::getDefaultVersion());
    }
    
    if (object) mItems.push_back(object);
  }

  return object;
}
/** @endcond */


/** @cond doxygen-c-only */


/**
 * Creates a new Species_t structure using the given SBML @p level
 * and @p version values.
 *
 * @param level an unsigned int, the SBML Level to assign to this
 * Species
 *
 * @param version an unsigned int, the SBML Version to assign to this
 * Species
 *
 * @return a pointer to the newly created Species_t structure.
 *
 * @note Once a Species has been added to an SBMLDocument, the @p
 * level and @p version for the document @em override those used to create
 * the Species.  Despite this, the ability to supply the values at
 * creation time is an important aid to creating valid SBML.  Knowledge of
 * the intended SBML Level and Version  determine whether it is valid to
 * assign a particular value to an attribute, or whether it is valid to add
 * an object to an existing SBMLDocument.
 */
LIBSBML_EXTERN
Species_t *
Species_create (unsigned int level, unsigned int version)
{
  try
  {
    Species* obj = new Species(level,version);
    return obj;
  }
  catch (SBMLConstructorException)
  {
    return NULL;
  }
}


/**
 * Creates a new Species_t structure using the given
 * SBMLNamespaces_t structure.
 *
 * @param sbmlns SBMLNamespaces, a pointer to an SBMLNamespaces structure
 * to assign to this Species
 *
 * @return a pointer to the newly created Species_t structure.
 *
 * @note Once a Species has been added to an SBMLDocument, the
 * @p sbmlns namespaces for the document @em override those used to create
 * the Species.  Despite this, the ability to supply the values at creation time
 * is an important aid to creating valid SBML.  Knowledge of the intended SBML
 * Level and Version determine whether it is valid to assign a particular value
 * to an attribute, or whether it is valid to add an object to an existing
 * SBMLDocument.
 */
LIBSBML_EXTERN
Species_t *
Species_createWithNS (SBMLNamespaces_t* sbmlns)
{
  try
  {
    Species* obj = new Species(sbmlns);
    return obj;
  }
  catch (SBMLConstructorException)
  {
    return NULL;
  }
}


/**
 * Frees the given Species_t structure.
 *
 * @param s the Species_t structure to be freed.
 */
LIBSBML_EXTERN
void
Species_free (Species_t *s)
{
  delete s;
}


/**
 * Creates a deep copy of the given Species_t structure
 * 
 * @param p the Species_t structure to be copied
 * 
 * @return a (deep) copy of the given Species_t structure.
 */
LIBSBML_EXTERN
Species_t *
Species_clone (const Species_t *s)
{
  return static_cast<Species*>( s->clone() );
}


/**
 * Initializes the attributes of the given Species_t structure to the
 * defaults defined in the specification of the relevant Level/Version of
 * SBML.
 * 
 * @li sets "boundaryCondition" to @c 1 (true)
 * @li (Level 2 only) sets "constant" to @c 0 (false)
 * @li (Level 2 only) sets "hasOnlySubstanceUnits" to @c 0 (false)
 *
 * @param s the Species_t structure.
 */
LIBSBML_EXTERN
void
Species_initDefaults (Species_t *s)
{
  s->initDefaults();
}


/**
 * Returns a list of XMLNamespaces_t associated with this Species_t
 * structure.
 *
 * @param s the Species_t structure
 * 
 * @return pointer to the XMLNamespaces_t structure associated with 
 * this SBML object
 */
LIBSBML_EXTERN
const XMLNamespaces_t *
Species_getNamespaces(Species_t *s)
{
  return s->getNamespaces();
}

/**
 * Takes a Species_t structure and returns its identifier.
 *
 * @param s the Species_t structure whose identifier is sought
 * 
 * @return the identifier of the Species_t structure @p s, as a pointer
 * to a string.
 */
LIBSBML_EXTERN
const char *
Species_getId (const Species_t *s)
{
  return s->isSetId() ? s->getId().c_str() : NULL;
}


/**
 * Takes a Species_t structure and returns its name.
 *
 * @param s the Species_t structure whose name is sought.
 *
 * @return the name of the Species_t structure @p s, as a pointer to a
 * string.
 */
LIBSBML_EXTERN
const char *
Species_getName (const Species_t *s)
{
  return s->isSetName() ? s->getName().c_str() : NULL;
}


/**
 * Get the species type of this Species_t structure, as indicated by the
 * Species_t structure object's "speciesType" attribute value.
 *
 * @param s the Species_t structure
 * 
 * @return the value of the "speciesType" attribute of the
 * Species_t structure @p s as a string.
 */
LIBSBML_EXTERN
const char *
Species_getSpeciesType (const Species_t *s)
{
  return s->isSetSpeciesType() ? s->getSpeciesType().c_str() : NULL;
}


/**
 * Get the compartment in which this species is located.
 *
 * @param s the Species_t structure
 * 
 * @return the value of the "compartment" attribute of the given Species_t
 * structure, as a string.
 */
LIBSBML_EXTERN
const char *
Species_getCompartment (const Species_t *s)
{
  return s->isSetCompartment() ? s->getCompartment().c_str() : NULL;
}


/**
 * Get the value of the "initialAmount" attribute.
 *
 * @param s the Species_t structure
 * 
 * @return the "initialAmount" attribute of the given Species_t structure,
 * as a float-point number.
 */
LIBSBML_EXTERN
double
Species_getInitialAmount (const Species_t *s)
{
  return s->getInitialAmount();
}


/**
 * Get the value of the "initialConcentration" attribute.
 *
 * @param s the Species_t structure
 * 
 * @return the "initialConcentration" of the given Species_t structure, as
 * a float-point number.
 */
LIBSBML_EXTERN
double
Species_getInitialConcentration (const Species_t *s)
{
  return s->getInitialConcentration();
}


/**
 * Get the value of the "substanceUnit" attribute.
 *
 * @param s the Species_t structure
 * 
 * @return the "substanceUnits" attribute of the given Species_t structure,
 * as a string.
 */
LIBSBML_EXTERN
const char *
Species_getSubstanceUnits (const Species_t *s)
{
  return s->isSetSubstanceUnits() ? s->getSubstanceUnits().c_str() : NULL;
}


/**
 * Get the value of the "spatialSizeUnits" attribute.
 *
 * @param s the Species_t structure
 * 
 * @return the spatialSizeUnits of the given Species.
 * 
 * @warning In versions of SBML Level~2 before Version 3, the Species_t
 * structure included an attribute called "spatialSizeUnits", which allowed
 * explicitly setting the units of size for initial concentration.  This
 * attribute was removed in SBML Level 2 Version 3.  LibSBML retains this
 * attribute for compatibility with older definitions of Level 2, but its
 * use is strongly discouraged because it is incompatible with Level 2
 * Versions 3 and 4.
 */
LIBSBML_EXTERN
const char *
Species_getSpatialSizeUnits (const Species_t *s)
{
  return s->isSetSpatialSizeUnits() ? s->getSpatialSizeUnits().c_str() : NULL;
}


/**
 * (SBML Level 1 only) Get the value of the "units" attribute.
 *
 * @param s the Species_t structure
 * 
 * @return the units of the given Species_t structure.
 */
LIBSBML_EXTERN
const char *
Species_getUnits (const Species_t *s)
{
  return s->isSetUnits() ? s->getUnits().c_str() : NULL;
}


/**
 * Get the value of the "hasOnlySubstanceUnits" attribute.
 *
 * @param s the Species_t structure
 * 
 * @return nonzero (true) if the given Species_t structure's
 * "hasOnlySubstanceUnits" attribute value is nonzero, zero (0) otherwise.
 */
LIBSBML_EXTERN
int
Species_getHasOnlySubstanceUnits (const Species_t *s)
{
  return static_cast<int>( s->getHasOnlySubstanceUnits() );
}


/**
 * Get the value of the "boundaryCondition" attribute.
 *
 * @param s the Species_t structure
 * 
 * @return nonzero (true) if the given Species_t structure's
 * "boundaryCondition" attribute value is nonzero, zero (0) otherwise.
 */
LIBSBML_EXTERN
int
Species_getBoundaryCondition (const Species_t *s)
{
  return static_cast<int>( s->getBoundaryCondition() );
}


/**
 * Get the value of the "charge" attribute.
 *
 * @param s the Species_t structure
 * 
 * @return the charge of the given Species_t structure.
 *
 * @note Beginning in SBML Level 2 Version 2, the "charge" attribute on
 * Species_t is deprecated and its use strongly discouraged.  Its presence
 * is considered a misfeature in earlier definitions of SBML because its
 * implications for the mathematics of a model were never defined, and in
 * any case, no known modeling system ever used it.  Instead, models take
 * account of charge values directly in their definitions of species by
 * (for example) having separate species identities for the charged and
 * uncharged versions of the same species.  This allows the condition to
 * affect model mathematics directly.  LibSBML retains this method for
 * easier compatibility with SBML Level 1.
 */
LIBSBML_EXTERN
int
Species_getCharge (const Species_t *s)
{ 
  return s->getCharge();
}


/**
 * Get the value of the "constant" attribute.
 *
 * @param s the Species_t structure
 * 
 * @return nonzero (true) if the given Species_t structure's "constant"
 * attribute value is nonzero, @c false otherwise.
 */
LIBSBML_EXTERN
int
Species_getConstant (const Species_t *s)
{
  return static_cast<int>( s->getConstant() );
}


/**
 * Get the value of the "conversionFactor" attribute.
 *
 * @param s the Species_t structure
 * 
 * @return the "conversionFactor" attribute of the given Species_t structure,
 * as a string.
 */
LIBSBML_EXTERN
const char *
Species_getConversionFactor (const Species_t *s)
{
  return s->isSetConversionFactor() ? s->getConversionFactor().c_str() : NULL;
}


/**
 * Predicate returning true or false depending on whether the attribute
 * "id" of the given Species_t structure has been set.
 *
 * @param s the Species_t structure
 * 
 * @return true (nonzero) if the "id" attribute of the given Species_t
 * structure has been set, false (0) otherwise.
 */
LIBSBML_EXTERN
int
Species_isSetId (const Species_t *s)
{
  return static_cast<int>( s->isSetId() );
}


/**
 * Predicate returning true or false depending on whether the attribute
 * "name" of the given Species_t structure has been set.
 *
 * @param s the Species_t structure
 * 
 * @return true (nonzero) if the "name" attribute of the given Species_t
 * structure has been set, false (0) otherwise.
 */
LIBSBML_EXTERN
int
Species_isSetName (const Species_t *s)
{
  return static_cast<int>( s->isSetName() );
}


/**
 * Predicate returning true or false depending on whether the attribute
 * "speciesType" of the given Species_t structure has been set.
 *
 * @param s the Species_t structure
 * 
 * @return true (nonzero) if the "speciesType" attribute of the given
 * Species_t structure has been set, false (0) otherwise.
 */
LIBSBML_EXTERN
int
Species_isSetSpeciesType (const Species_t *s)
{
  return static_cast<int>( s->isSetSpeciesType() );
}


/**
 * Predicate returning true or false depending on whether the attribute
 * "compartment" of the given Species_t structure has been set.
 *
 * @param s the Species_t structure
 * 
 * @return true (nonzero) if the "compartment" attribute of the given
 * Species_t structure has been set, false (0) otherwise.
 */
LIBSBML_EXTERN
int
Species_isSetCompartment (const Species_t *s)
{
  return static_cast<int>( s->isSetCompartment() );
}


/**
 * Predicate returning true or false depending on whether the attribute
 * "initialAmount" of the given Species_t structure has been set.
 *
 * @param s the Species_t structure
 * 
 * @return true (nonzero) if the "initialAmount" attribute of the given
 * Species_t structure has been set, false (0) otherwise.
 *
 * @note In SBML Level 1, Species_t' "initialAmount" is required and
 * therefore <em>should always be set</em>.  (However, in Level 1, the
 * attribute has no default value either, so this method will not return
 * nonzero until a value has been assigned.)  In SBML Level 2,
 * "initialAmount" is optional and as such may or may not be set.
 */
LIBSBML_EXTERN
int
Species_isSetInitialAmount (const Species_t *s)
{
  return static_cast<int>( s->isSetInitialAmount() );
}


/**
 * Predicate returning true or false depending on whether the attribute
 * "compartment" of the given Species_t structure has been set.
 *
 * @param s the Species_t structure
 * 
 * @return true (nonzero) if the "compartment" attribute of the given
 * Species_t structure has been set, false (0) otherwise.
 */
LIBSBML_EXTERN
int
Species_isSetInitialConcentration (const Species_t *s)
{
  return static_cast<int>( s->isSetInitialConcentration() );
}


/**
 * Predicate returning true or false depending on whether the attribute
 * "substanceUnits" of the given Species_t structure has been set.
 *
 * @param s the Species_t structure
 * 
 * @return true (nonzero) if the "substanceUnits" attribute of the given
 * Species_t structure has been set, false (0) otherwise.
 */
LIBSBML_EXTERN
int
Species_isSetSubstanceUnits (const Species_t *s)
{
  return static_cast<int>( s->isSetSubstanceUnits() );
}


/**
 * Predicate returning true or false depending on whether the attribute
 * "spatialSizeUnits" of the given Species_t structure has been set.
 *
 * @param s the Species_t structure
 * 
 * @return true (nonzero) if the "spatialSizeUnits" attribute of the given
 * Species_t structure has been set, false (0) otherwise.
 * 
 * @warning In versions of SBML Level~2 before Version 3, the class
 * Species included an attribute called "spatialSizeUnits", which allowed
 * explicitly setting the units of size for initial concentration.  This
 * attribute was removed in SBML Level 2 Version 3.  LibSBML retains this
 * attribute for compatibility with older definitions of Level 2, but its
 * use is strongly discouraged because it is incompatible with Level 2
 * Versions 3 and 4.
 */
LIBSBML_EXTERN
int
Species_isSetSpatialSizeUnits (const Species_t *s)
{
  return static_cast<int>( s->isSetSpatialSizeUnits() );
}


/**
 * (SBML Level 1 only) Predicate returning true or false depending on
 * whether the attribute "units" of the given Species_t structure has been
 * set.
 *
 * @param s the Species_t structure
 * 
 * @return true (nonzero) if the "units" attribute of the given
 * Species_t structure has been set, false (0) otherwise.
 */
LIBSBML_EXTERN
int
Species_isSetUnits (const Species_t *s)
{
  return static_cast<int>( s->isSetUnits() );
}


/**
 * Predicate returning true or false depending on whether the attribute
 * "charge" of the given Species_t structure has been set.
 *
 * @param s the Species_t structure
 * 
 * @return true (nonzero) if the "charge" attribute of the given
 * Species_t structure has been set, false (0) otherwise.
 *
 * @note Beginning in SBML Level 2 Version 2, the "charge" attribute on
 * Species in SBML is deprecated and its use strongly discouraged.  Its
 * presence is considered a misfeature in earlier definitions of SBML
 * because its implications for the mathematics of a model were never
 * defined, and in any case, no known modeling system ever used it.
 * Instead, models take account of charge values directly in their
 * definitions of species by (for example) having separate species
 * identities for the charged and uncharged versions of the same species.
 * This allows the condition to affect model mathematics directly.
 * LibSBML retains this method for easier compatibility with SBML Level 1.
 */
LIBSBML_EXTERN
int
Species_isSetCharge (const Species_t *s)
{
  return static_cast<int>( s->isSetCharge() );
}


/**
 * Predicate returning true or false depending on whether the attribute
 * "conversionFactor" of the given Species_t structure has been set.
 *
 * @param s the Species_t structure
 * 
 * @return true (nonzero) if the "conversionFactor" attribute of the given Species_t
 * structure has been set, false (0) otherwise.
 */
LIBSBML_EXTERN
int
Species_isSetConversionFactor (const Species_t *s)
{
  return static_cast<int>( s->isSetConversionFactor() );
}


/**
 * Predicate returning true or false depending on whether the attribute
 * "constant" of the given Species_t structure has been set.
 *
 * @param s the Species_t structure
 * 
 * @return true (nonzero) if the "constant" attribute of the given Species_t
 * structure has been set, false (0) otherwise.
 */
LIBSBML_EXTERN
int
Species_isSetConstant (const Species_t *s)
{
  return static_cast<int>( s->isSetConstant() );
}


/**
 * Predicate returning true or false depending on whether the attribute
 * "boundaryCondition" of the given Species_t structure has been set.
 *
 * @param s the Species_t structure
 * 
 * @return true (nonzero) if the "boundaryCondition" attribute of the given Species_t
 * structure has been set, false (0) otherwise.
 */
LIBSBML_EXTERN
int
Species_isSetBoundaryCondition (const Species_t *s)
{
  return static_cast<int>( s->isSetBoundaryCondition() );
}


/**
 * Predicate returning true or false depending on whether the attribute
 * "hasOnlySubstanceUnits" of the given Species_t structure has been set.
 *
 * @param s the Species_t structure
 * 
 * @return true (nonzero) if the "hasOnlySubstanceUnits" attribute of the given Species_t
 * structure has been set, false (0) otherwise.
 */
LIBSBML_EXTERN
int
Species_isSetHasOnlySubstanceUnits (const Species_t *s)
{
  return static_cast<int>( s->isSetHasOnlySubstanceUnits() );
}


/**
 * Sets the "id" attribute of the given Species_t structure.
 *
 * This function copies the string given in @p sid.  If the string is
 * NULL, this function performs unsetId() instead.
 *
 * @param s the Species_t structure
 * 
 * @param sid the identifier string to which the "id" attribute should be
 * set.
 *
 * @return integer value indicating success/failure of the
 * function.  @if clike The value is drawn from the
 * enumeration #OperationReturnValues_t. @endif The possible values
 * returned by this function are:
 * @li LIBSBML_OPERATION_SUCCESS
 * @li LIBSBML_INVALID_ATTRIBUTE_VALUE
 *
 * @note Using this function with an id of NULL is equivalent to
 * unsetting the "id" attribute.
 */
LIBSBML_EXTERN
int
Species_setId (Species_t *s, const char *sid)
{
  return (sid == NULL) ? s->setId("") : s->setId(sid);
}


/**
 * Sets the "name" attribute of the given Species_t structure.
 *
 * This function copies the string given in @p name.  If the string is
 * NULL, this function performs unsetName() instead.
 *
 * @param s the Species_t structure
 * 
 * @param name the name string to which the "name" attribute should be set.
 *
 * @return integer value indicating success/failure of the
 * function.  @if clike The value is drawn from the
 * enumeration #OperationReturnValues_t. @endif The possible values
 * returned by this function are:
 * @li LIBSBML_OPERATION_SUCCESS
 * @li LIBSBML_INVALID_ATTRIBUTE_VALUE
 *
 * @note Using this function with the name set to NULL is equivalent to
 * unsetting the "name" attribute.
 */
LIBSBML_EXTERN
int
Species_setName (Species_t *s, const char *name)
{
  return (name == NULL) ? s->unsetName() : s->setName(name);
}


/**
 * Sets the "speciesType" attribute of the given Species_t structure.
 *
 * This function copies the string given in @p sid.  If the string
 * is NULL, this function performs unsetSpeciesType() instead.
 *
 * @param s the Species_t structure
 * 
 * @param speciesType the identifer to which the "speciesType" attribute
 * should be set.
 *
 * @return integer value indicating success/failure of the
 * function.  @if clike The value is drawn from the
 * enumeration #OperationReturnValues_t. @endif The possible values
 * returned by this function are:
 * @li LIBSBML_OPERATION_SUCCESS
 * @li LIBSBML_INVALID_ATTRIBUTE_VALUE
 * @li LIBSBML_UNEXPECTED_ATTRIBUTE
 *
 * @note Using this function with an id of NULL is equivalent to
 * unsetting the "speciesType" attribute.
 */
LIBSBML_EXTERN
int
Species_setSpeciesType (Species_t *s, const char *sid)
{
  return (sid == NULL) ? s->unsetSpeciesType() : s->setSpeciesType(sid);
}


/**
 * Sets the "compartment" attribute of the given Species_t structure.
 *
 * This function copies the string given in @p compartment.  If the string
 * is NULL, this function performs unsetCompartment() instead.
 *
 * @param s the Species_t structure
 * 
 * @param sid the identifer to which the "compartment" attribute
 * should be set.
 *
 * @return integer value indicating success/failure of the
 * function.  @if clike The value is drawn from the
 * enumeration #OperationReturnValues_t. @endif The possible values
 * returned by this function are:
 * @li LIBSBML_OPERATION_SUCCESS
 * @li LIBSBML_INVALID_ATTRIBUTE_VALUE
 *
 * @note Using this function with an id of NULL is equivalent to
 * unsetting the "compartment" attribute.
 */
LIBSBML_EXTERN
int
Species_setCompartment (Species_t *s, const char *sid)
{
  return (sid == NULL) ? s->setCompartment("") : s->setCompartment(sid);
}


/**
 * Sets the "initialAmount" attribute value of the given Species_t
 * structure.
 *
 * As a side-effect, calling this function also unsets the
 * "initialConcentration" attribute.
 *
 * @param s the Species_t structure
 *
 * @param value the numerical value for the "initialAmount" attribute
 *
 * @return integer value indicating success/failure of the
 * function.  @if clike The value is drawn from the
 * enumeration #OperationReturnValues_t. @endif The possible values
 * returned by this function are:
 * @li LIBSBML_OPERATION_SUCCESS
 */
LIBSBML_EXTERN
int
Species_setInitialAmount (Species_t *s, double value)
{
  return s->setInitialAmount(value);
}


/**
 * Sets the "initialConcentration" attribute value of the given Species_t
 * structure.
 *
 * As a side-effect, calling this function also unsets the "initialAmount"
 * attribute.
 *
 * @param s the Species_t structure
 *
 * @param value the numerical value for the "initialConcentration"
 * attribute
 *
 * @return integer value indicating success/failure of the
 * function.  @if clike The value is drawn from the
 * enumeration #OperationReturnValues_t. @endif The possible values
 * returned by this function are:
 * @li LIBSBML_OPERATION_SUCCESS
 * @li LIBSBML_UNEXPECTED_ATTRIBUTE
 */
LIBSBML_EXTERN
int
Species_setInitialConcentration (Species_t *s, double value)
{
  return s->setInitialConcentration(value);
}


/**
 * Sets the "substanceUnits" attribute of the given Species_t structure.
 *
 * This function copies the string given in @p sid.  If the string
 * is NULL, this function performs unsetSubstanceUnits() instead.
 *
 * @param s the Species_t structure
 * 
 * @param substanceUnits the identifer to which the "substanceUnits"
 * attribute should be set.
 *
 * @return integer value indicating success/failure of the
 * function.  @if clike The value is drawn from the
 * enumeration #OperationReturnValues_t. @endif The possible values
 * returned by this function are:
 * @li LIBSBML_OPERATION_SUCCESS
 * @li LIBSBML_INVALID_ATTRIBUTE_VALUE
 *
 * @note Using this function with an id of NULL is equivalent to
 * unsetting the "substanceUnits" attribute.
 */
LIBSBML_EXTERN
int
Species_setSubstanceUnits (Species_t *s, const char *sid)
{
  return (sid == NULL) ? s->unsetSubstanceUnits() : s->setSubstanceUnits(sid);
}


/**
 * Sets the "spatialSizeUnits" attribute of the given Species_t structure.
 *
 * This function copies the string given in @p sid.  If the string is NULL,
 * this function performs unsetSpatialSizeUnits() instead.
 *
 * @param s the Species_t structure
 * 
 * @param spatialSizeUnits the identifer to which the "spatialSizeUnits"
 * attribute should be set.
 *
 * @return integer value indicating success/failure of the
 * function.  @if clike The value is drawn from the
 * enumeration #OperationReturnValues_t. @endif The possible values
 * returned by this function are:
 * @li LIBSBML_OPERATION_SUCCESS
 * @li LIBSBML_INVALID_ATTRIBUTE_VALUE
 * @li LIBSBML_UNEXPECTED_ATTRIBUTE
 *
 * @note Using this function with an id of NULL is equivalent to
 * unsetting the "spatialSizeUnits" attribute.
 * 
 * @warning In versions of SBML Level~2 before Version 3, the class
 * Species included an attribute called "spatialSizeUnits", which allowed
 * explicitly setting the units of size for initial concentration.  This
 * attribute was removed in SBML Level 2 Version 3.  LibSBML retains this
 * attribute for compatibility with older definitions of Level 2, but its
 * use is strongly discouraged because it is incompatible with Level 2
 * Versions 3 and 4.
 */
LIBSBML_EXTERN
int
Species_setSpatialSizeUnits (Species_t *s, const char *sid)
{
  return (sid == NULL) ? s->unsetSpatialSizeUnits() : s->setSpatialSizeUnits(sid);
}


/**
 * (SBML Level 1 only) Sets the "units" attribute of the given Species_t
 * structure.
 *
 * This function copies the string given in @p sid.  If the string is NULL,
 * this function performs unsetUnits() instead.
 *
 * @param s the Species_t structure
 * 
 * @param units the identifer to which the "units" attribute
 * should be set.
 *
 * @return integer value indicating success/failure of the
 * function.  @if clike The value is drawn from the
 * enumeration #OperationReturnValues_t. @endif The possible values
 * returned by this function are:
 * @li LIBSBML_OPERATION_SUCCESS
 * @li LIBSBML_INVALID_ATTRIBUTE_VALUE
 *
 * @note Using this function with an id of NULL is equivalent to
 * unsetting the "units" attribute.
 */
LIBSBML_EXTERN
int
Species_setUnits (Species_t *s, const char *sname)
{
  return (sname == NULL) ? s->unsetUnits() : s->setUnits(sname);
}


/**
 * Sets the "hasOnlySubstanceUnits" attribute of the given Species_t
 * structure.
 *
 * @param s the Species_t structure
 * 
 * @param value nonzero to indicate true, zero to indicate false.
 *
 * @return integer value indicating success/failure of the
 * function.  @if clike The value is drawn from the
 * enumeration #OperationReturnValues_t. @endif The possible values
 * returned by this function are:
 * @li LIBSBML_OPERATION_SUCCESS
 * @li LIBSBML_UNEXPECTED_ATTRIBUTE
 */
LIBSBML_EXTERN
int
Species_setHasOnlySubstanceUnits (Species_t *s, int value)
{
  return s->setHasOnlySubstanceUnits( static_cast<bool>(value) );
}


/**
 * Sets the "boundaryCondition" attribute of the given Species_t
 * structure.
 *
 * @param s the Species_t structure
 * 
 * @param value nonzero to indicate true, zero to indicate false.
 *
 * @return integer value indicating success/failure of the
 * function.  @if clike The value is drawn from the
 * enumeration #OperationReturnValues_t. @endif The possible values
 * returned by this function are:
 * @li LIBSBML_OPERATION_SUCCESS
 */
LIBSBML_EXTERN
int
Species_setBoundaryCondition (Species_t *s, int value)
{
  return s->setBoundaryCondition( static_cast<bool>(value) );
}


/**
 * Sets the "charge" attribute of the given Species_t
 * structure.
 *
 * @param s the Species_t structure
 * 
 * @param value the value of charge to assign to the "charge" attribute
 *
 * @return integer value indicating success/failure of the
 * function.  @if clike The value is drawn from the
 * enumeration #OperationReturnValues_t. @endif The possible values
 * returned by this function are:
 * @li LIBSBML_OPERATION_SUCCESS
 * @li LIBSBML_UNEXPECTED_ATTRIBUTE
 *
 * @note Beginning in SBML Level 2 Version 2, the "charge" attribute on
 * Species in SBML is deprecated and its use strongly discouraged.  Its
 * presence is considered a misfeature in earlier definitions of SBML
 * because its implications for the mathematics of a model were never
 * defined, and in any case, no known modeling system ever used it.
 * Instead, models take account of charge values directly in their
 * definitions of species by (for example) having separate species
 * identities for the charged and uncharged versions of the same species.
 * This allows the condition to affect model mathematics directly.
 * LibSBML retains this method for easier compatibility with SBML Level 1.
 */
LIBSBML_EXTERN
int
Species_setCharge (Species_t *s, int value)
{
  return s->setCharge(value);
}


/**
 * Sets the "constant" attribute of the given Species_t
 * structure.
 *
 * @param s the Species_t structure
 * 
 * @param value nonzero to indicate true, zero to indicate false.
 *
 * @return integer value indicating success/failure of the
 * function.  @if clike The value is drawn from the
 * enumeration #OperationReturnValues_t. @endif The possible values
 * returned by this function are:
 * @li LIBSBML_OPERATION_SUCCESS
 * @li LIBSBML_UNEXPECTED_ATTRIBUTE
 */
LIBSBML_EXTERN
int
Species_setConstant (Species_t *s, int value)
{
  return s->setConstant( static_cast<bool>(value) );
}


/**
 * Sets the "conversionFactor" attribute of the given Species_t structure.
 *
 * This function copies the string given in @p sid.  If the string
 * is NULL, this function performs unsetConversionFactor() instead.
 *
 * @param s the Species_t structure
 * 
 * @param conversionFactor the identifer to which the "conversionFactor" attribute
 * should be set.
 *
 * @return integer value indicating success/failure of the
 * function.  @if clike The value is drawn from the
 * enumeration #OperationReturnValues_t. @endif The possible values
 * returned by this function are:
 * @li LIBSBML_OPERATION_SUCCESS
 * @li LIBSBML_INVALID_ATTRIBUTE_VALUE
 * @li LIBSBML_UNEXPECTED_ATTRIBUTE
 *
 * @note Using this function with an id of NULL is equivalent to
 * unsetting the "conversionFactor" attribute.
 */
LIBSBML_EXTERN
int
Species_setConversionFactor (Species_t *s, const char *sid)
{
  return (sid == NULL) ? s->unsetConversionFactor() : 
                         s->setConversionFactor(sid);
}


/**
 * Unsets the "name" attribute of the given Species_t structure.
 *
 * @return integer value indicating success/failure of the
 * function.  @if clike The value is drawn from the
 * enumeration #OperationReturnValues_t. @endif The possible values
 * returned by this function are:
 * @li LIBSBML_OPERATION_SUCCESS
 * @li LIBSBML_OPERATION_FAILED
 */
LIBSBML_EXTERN
int
Species_unsetName (Species_t *s)
{
  return s->unsetName();
}


/**
 * Unsets the "speciesType" attribute of the given Species_t structure.
 *
 * @param s the Species_t structure whose attribute is to be unset.
 *
 * @return integer value indicating success/failure of the
 * function.  @if clike The value is drawn from the
 * enumeration #OperationReturnValues_t. @endif The possible values
 * returned by this function are:
 * @li LIBSBML_OPERATION_SUCCESS
 * @li LIBSBML_OPERATION_FAILED
 */
LIBSBML_EXTERN
int
Species_unsetSpeciesType (Species_t *s)
{
  return s->unsetSpeciesType();
}


/**
 * Unsets the "initialAmount" attribute of the given Species_t structure.
 *
 * @param s the Species_t structure whose attribute is to be unset.
 *
 * @return integer value indicating success/failure of the
 * function.  @if clike The value is drawn from the
 * enumeration #OperationReturnValues_t. @endif The possible values
 * returned by this function are:
 * @li LIBSBML_OPERATION_SUCCESS
 */
LIBSBML_EXTERN
int
Species_unsetInitialAmount (Species_t *s)
{
  return s->unsetInitialAmount();
}


/**
 * Unsets the "initialConcentration" attribute of the given
 * Species_t structure.
 *
 * @param s the Species_t structure whose attribute is to be unset.
 *
 * @return integer value indicating success/failure of the
 * function.  @if clike The value is drawn from the
 * enumeration #OperationReturnValues_t. @endif The possible values
 * returned by this function are:
 * @li LIBSBML_OPERATION_SUCCESS
 */
LIBSBML_EXTERN
int
Species_unsetInitialConcentration (Species_t *s)
{
  return s->unsetInitialConcentration();
}


/**
 * Unsets the "substanceUnits" attribute of the given
 * Species_t structure.
 *
 * @param s the Species_t structure whose attribute is to be unset.
 *
 * @return integer value indicating success/failure of the
 * function.  @if clike The value is drawn from the
 * enumeration #OperationReturnValues_t. @endif The possible values
 * returned by this function are:
 * @li LIBSBML_OPERATION_SUCCESS
 * @li LIBSBML_OPERATION_FAILED
 */
LIBSBML_EXTERN
int
Species_unsetSubstanceUnits (Species_t *s)
{
  return s->unsetSubstanceUnits();
}


/**
 * Unsets the "spatialSizeUnits" attribute of the given
 * Species_t structure.
 *
 * @param s the Species_t structure whose attribute is to be unset.
 *
 * @return integer value indicating success/failure of the
 * function.  @if clike The value is drawn from the
 * enumeration #OperationReturnValues_t. @endif The possible values
 * returned by this function are:
 * @li LIBSBML_OPERATION_SUCCESS
 * @li LIBSBML_OPERATION_FAILED
 * 
 * @warning In versions of SBML Level~2 before Version 3, the class
 * Species included an attribute called "spatialSizeUnits", which allowed
 * explicitly setting the units of size for initial concentration.  This
 * attribute was removed in SBML Level 2 Version 3.  LibSBML retains this
 * attribute for compatibility with older definitions of Level 2, but its
 * use is strongly discouraged because it is incompatible with Level 2
 * Versions 3 and 4.
 */
LIBSBML_EXTERN
int
Species_unsetSpatialSizeUnits (Species_t *s)
{
  return s->unsetSpatialSizeUnits();
}


/**
 * Unsets the "units" attribute of the given
 * Species_t structure.
 *
 * @param s the Species_t structure whose attribute is to be unset.
 *
 * @return integer value indicating success/failure of the
 * function.  @if clike The value is drawn from the
 * enumeration #OperationReturnValues_t. @endif The possible values
 * returned by this function are:
 * @li LIBSBML_OPERATION_SUCCESS
 * @li LIBSBML_OPERATION_FAILED
 */
LIBSBML_EXTERN
int
Species_unsetUnits (Species_t *s)
{
  return s->unsetUnits();
}


/**
 * Unsets the "charge" attribute of the given
 * Species_t structure.
 *
 * @param s the Species_t structure whose attribute is to be unset.
 *
 * @return integer value indicating success/failure of the
 * function.  @if clike The value is drawn from the
 * enumeration #OperationReturnValues_t. @endif The possible values
 * returned by this function are:
 * @li LIBSBML_OPERATION_SUCCESS
 * @li LIBSBML_OPERATION_FAILED
 *
 * @note Beginning in SBML Level 2 Version 2, the "charge" attribute on
 * Species in SBML is deprecated and its use strongly discouraged.  Its
 * presence is considered a misfeature in earlier definitions of SBML
 * because its implications for the mathematics of a model were never
 * defined, and in any case, no known modeling system ever used it.
 * Instead, models take account of charge values directly in their
 * definitions of species by (for example) having separate species
 * identities for the charged and uncharged versions of the same species.
 * This allows the condition to affect model mathematics directly.
 * LibSBML retains this method for easier compatibility with SBML Level 1.
 */
LIBSBML_EXTERN
int
Species_unsetCharge (Species_t *s)
{
  return s->unsetCharge();
}


/**
 * Unsets the "conversionFactor" attribute of the given Species_t structure.
 *
 * @param s the Species_t structure whose attribute is to be unset.
 *
 * @return integer value indicating success/failure of the
 * function.  @if clike The value is drawn from the
 * enumeration #OperationReturnValues_t. @endif The possible values
 * returned by this function are:
 * @li LIBSBML_OPERATION_SUCCESS
 * @li LIBSBML_OPERATION_FAILED
 */
LIBSBML_EXTERN
int
Species_unsetConversionFactor (Species_t *s)
{
  return s->unsetConversionFactor();
}


/**
 * Constructs and returns a UnitDefinition_t structure that expresses 
 * the units of this Species_t structure.
 *
 * @param s the Species_t structure whose units are to be returned.
 *
 * @return a UnitDefinition_t structure that expresses the units 
 * of this Species_t strucuture.
 *
 * @note This function returns the units of the Species_t expressed 
 * as a UnitDefinition_t. The units may be those explicitly declared 
 * or those derived from the default units of the Model_t containing
 * this Species_t and it's Compartment_t.
 *
 * Note that the functionality that facilitates unit analysis depends 
 * on the model as a whole.  Thus, in cases where the object has not 
 * been added to a model or the model itself is incomplete,
 * unit analysis is not possible and this method will return NULL.
 *
 */
LIBSBML_EXTERN
UnitDefinition_t * 
Species_getDerivedUnitDefinition(Species_t *s)
{
  return s->getDerivedUnitDefinition();
}


/**
 * Predicate returning @c true or @c false depending on whether
 * all the required attributes for this Species object
 * have been set.
 *
 * @param s the Species_t structure to check.
 *
 * @note The required attributes for a Species object are:
 * @li id (name L1)
 * @li compartment
 * @li initialAmount (L1 only)
 * @li hasOnlySubstanceUnits (L3 on)
 * @li boundaryCondition (L3 on)
 * @li constant (L3 on)
 *
 * @return a true if all the required
 * attributes for this object have been defined, false otherwise.
 */
LIBSBML_EXTERN
int
Species_hasRequiredAttributes(Species_t *s)
{
  return static_cast<int>(s->hasRequiredAttributes());
}


/**
 * @return item in this ListOfSpecies with the given id or NULL if no such
 * item exists.
 */
LIBSBML_EXTERN
Species_t *
ListOfSpecies_getById (ListOf_t *lo, const char *sid)
{
  return (sid != NULL) ? 
    static_cast <ListOfSpecies *> (lo)->get(sid) : NULL;
}


/**
 * Removes item in this ListOf items with the given id or NULL if no such
 * item exists.  The caller owns the returned item and is responsible for
 * deleting it.
 */
LIBSBML_EXTERN
Species_t *
ListOfSpecies_removeById (ListOf_t *lo, const char *sid)
{
  return (sid != NULL) ? 
    static_cast <ListOfSpecies *> (lo)->remove(sid) : NULL;
}


/** @endcond */

LIBSBML_CPP_NAMESPACE_END
