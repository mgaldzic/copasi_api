/**
 * @file    Compartment.cpp
 * @brief   Implementations of Compartment and ListOfCompartments.
 * @author  Ben Bornstein
 *
 * $Id: Compartment.cpp 11635 2010-08-03 04:05:22Z mhucka $
 * $HeadURL: https://sbml.svn.sourceforge.net/svnroot/sbml/trunk/libsbml/src/sbml/Compartment.cpp $
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
#include <sbml/SBMLDocument.h>
#include <sbml/SBMLError.h>
#include <sbml/Model.h>
#include <sbml/Compartment.h>

/** @cond doxygen-ignored */

using namespace std;

/** @endcond */

LIBSBML_CPP_NAMESPACE_BEGIN

Compartment::Compartment (unsigned int level, unsigned int version) :
   SBase             ( level, version )
 , mId               ( ""       )
 , mName             ( ""       )
 , mSpatialDimensions( 3        )
 , mSpatialDimensionsDouble( 3        )
 , mSize             ( 1.0      )
 , mConstant         ( true     )
 , mIsSetSize        ( false    )
 , mIsSetSpatialDimensions ( false    )
 , mIsSetConstant          ( false    )
{
  if (!hasValidLevelVersionNamespaceCombination())
    throw SBMLConstructorException();

  // if level 3 values have no defaults
  if (level == 3)
  {
    mSize = numeric_limits<double>::quiet_NaN();
    mSpatialDimensionsDouble = numeric_limits<double>::quiet_NaN();
  }
}

Compartment::Compartment(SBMLNamespaces * sbmlns) :
   SBase             ( sbmlns   )
 , mId               ( ""       )
 , mName             ( ""       )
 , mSpatialDimensions( 3        )
 , mSpatialDimensionsDouble( 3        )
 , mSize             ( 1.0      )
 , mConstant         ( true     )
 , mIsSetSize        ( false    )
 , mIsSetSpatialDimensions ( false    )
 , mIsSetConstant          ( false    )
{
  if (!hasValidLevelVersionNamespaceCombination())
    throw SBMLConstructorException();
 
  // if level 3 values have no defaults
  if (sbmlns->getLevel() == 3)
  {
    mSize = numeric_limits<double>::quiet_NaN();
    mSpatialDimensionsDouble = numeric_limits<double>::quiet_NaN();
  }
}

/** @cond doxygen-libsbml-internal */

/* constructor for validators */
Compartment::Compartment() :
  SBase()
{
}

/** @endcond */
                          
/*
 * Destroys this Compartment.
 */
Compartment::~Compartment ()
{
}


/*
 * Copy constructor. Creates a copy of this compartment.
 */
Compartment::Compartment(const Compartment& orig) :
   SBase             ( orig                    )
 , mId               ( orig.mId                )  
 , mName             ( orig.mName              )
 , mCompartmentType  ( orig.mCompartmentType   )
 , mSpatialDimensions( orig.mSpatialDimensions )
 , mSpatialDimensionsDouble( orig.mSpatialDimensionsDouble )
 , mSize             ( orig.mSize              )
 , mUnits            ( orig.mUnits             )
 , mOutside          ( orig.mOutside           )
 , mConstant         ( orig.mConstant          )
 , mIsSetSize        ( orig.mIsSetSize         )
 , mIsSetSpatialDimensions ( orig.mIsSetSpatialDimensions    )
 , mIsSetConstant          ( orig.mIsSetConstant    )
{
}


/*
 * Assignment operator
 */
Compartment& Compartment::operator=(const Compartment& rhs)
{
  if(&rhs!=this)
  {
    this->SBase::operator =(rhs);
    mSpatialDimensions= rhs.mSpatialDimensions  ;
    mSpatialDimensionsDouble= rhs.mSpatialDimensions  ;
    mSize             = rhs.mSize      ;
    mConstant         = rhs.mConstant     ;
    mIsSetSize        = rhs.mIsSetSize    ;
    mCompartmentType  = rhs.mCompartmentType;
    mUnits            = rhs.mUnits ;
    mOutside          = rhs.mOutside ;
    mId               = rhs.mId;
    mName             = rhs.mName;
    mIsSetSpatialDimensions = rhs.mIsSetSpatialDimensions;
    mIsSetConstant          = rhs.mIsSetConstant;
  }

  return *this;
}



/*
 * Accepts the given SBMLVisitor.
 *
 * @return the result of calling <code>v.visit()</code>, which indicates
 * whether or not the Visitor would like to visit the Model's next
 * Compartment (if available).
 */
bool
Compartment::accept (SBMLVisitor& v) const
{
  return v.visit(*this);
}


/*
 * @return a (deep) copy of this Compartment.
 */
Compartment*
Compartment::clone () const
{
  return new Compartment(*this);
}


/*
 * Initializes the fields of this Compartment to their defaults:
 *
 *   - volume            = 1.0          (L1 only)
 *   - spatialDimensions = 3            (L2 only)
 *   - constant          = 1    (true)  (L2 only)
 */
LIBSBML_EXTERN
void
Compartment::initDefaults ()
{
  //// level 3 has no defaults
  //if (getLevel() < 3)
//  {
    mSize      = 1.0;    // Actually, setting L1 volume not
    mIsSetSize = false;  // L2 size.

    unsigned int dims = 3;
    setSpatialDimensions(dims);
    setConstant(1);
  //}
}


/*
 * @return the id of this SBML object.
 */
const string&
Compartment::getId () const
{
  return mId;
}


/*
 * @return the name of this SBML object.
 */
const string&
Compartment::getName () const
{
  return (getLevel() == 1) ? mId : mName;
}


/*
 * @return the compartmentType of this Compartment.
 */
const string&
Compartment::getCompartmentType () const
{
  return mCompartmentType;
}


/*
 * @return the spatialDimensions of this Compartment.
 */
unsigned int
Compartment::getSpatialDimensions () const
{
  if (getLevel() < 3)
  {
    return mSpatialDimensions;
  }
  else
  {
    if (isSetSpatialDimensions())
    {
      if (ceil(mSpatialDimensionsDouble) == 
          floor(mSpatialDimensionsDouble))
      {
        return static_cast<unsigned int>(mSpatialDimensionsDouble);
      }
      else
      {
        return numeric_limits<unsigned int>::quiet_NaN();
      }
    }
    else
    {
      return static_cast<unsigned int>(mSpatialDimensionsDouble);
    }
  }
}


/*
 * @return the spatialDimensions of this Compartment.
 */
double
Compartment::getSpatialDimensionsAsDouble () const
{
  if (getLevel() > 2)
    return mSpatialDimensionsDouble;
  else
    return static_cast<double>(mSpatialDimensions);
}


/*
 * @return the size (volume in L1) of this Compartment.
 */
double
Compartment::getSize () const
{
  return mSize;
}


/*
 * @return the volume (size in L2) of this Compartment.
 */
double
Compartment::getVolume () const
{
  return getSize();
}


/*
 * @return the units of this Compartment.
 */
const string&
Compartment::getUnits () const
{
  return mUnits;
}


/*
 * @return the outside of this Compartment.
 */
const string&
Compartment::getOutside () const
{
  return mOutside;
}


/*
 * @return true if this Compartment is constant, false otherwise.
 */
bool
Compartment::getConstant () const
{
  return mConstant;
}


/*
 * @return true if the id of this SBML object has been set, false
 * otherwise.
 */
bool
Compartment::isSetId () const
{
  return (mId.empty() == false);
}


/*
 * @return true if the name of this SBML object has been set, false
 * otherwise.
 */
bool
Compartment::isSetName () const
{
  return (getLevel() == 1) ? (mId.empty() == false) : 
                            (mName.empty() == false);
}


/*
 * @return true if the compartmentType of this Compartment has been set,
 * false otherwise. 
 */
bool
Compartment::isSetCompartmentType () const
{
  return (mCompartmentType.empty() == false);
}


/*
 * @return true if the size (volume in L1) of this Compartment has been
 * set, false otherwise.
 */
bool
Compartment::isSetSize () const
{
  return mIsSetSize;
}


/*
 * @return true if the volume (size in L2) of this Compartment has been
 * set, false otherwise.
 *
 * In SBML L1, a Compartment volume has a default value (1.0) and therefore
 * <b>should always be set</b>.  In L2, volume (size) is optional with no
 * default value and as such may or may not be set.
 */
bool
Compartment::isSetVolume () const
{
  return (getLevel() == 1) ? true : isSetSize();
}


/*
 * @return true if the units of this Compartment has been set, false
 * otherwise.
 */
bool
Compartment::isSetUnits () const
{
  return (mUnits.empty() == false);
}


/*
 * @return true if the outside of this Compartment has been set, false
 * otherwise.
 */
bool
Compartment::isSetOutside () const
{
  return (mOutside.empty() == false);
}


/*
 * @return true if the spatialDimenions of this Compartment has been set, false
 * otherwise.
 */
bool
Compartment::isSetSpatialDimensions () const
{
  return mIsSetSpatialDimensions;
}


/*
 * @return true if the constant of this Compartment has been set, false
 * otherwise.
 */
bool
Compartment::isSetConstant () const
{
  return mIsSetConstant;
}


/*
 * Sets the id of this SBML object to a copy of sid.
 */
int
Compartment::setId (const std::string& sid)
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
Compartment::setName (const std::string& name)
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
 * Sets the compartmentType field of this Compartment to a copy of sid.
 */
int
Compartment::setCompartmentType (const std::string& sid)
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
    mCompartmentType = sid;
    return LIBSBML_OPERATION_SUCCESS;
  }
}


/*
 * Sets the spatialDimensions of this Compartment to value.
 *
 * If value is not one of [0, 1, 2, 3] the function will have no effect
 * (i.e. spatialDimensions will not be set).
 */
int
Compartment::setSpatialDimensions (unsigned int value)
{
  return setSpatialDimensions((double) value);
}


/*
 * Sets the spatialDimensions of this Compartment to value.
 */
int
Compartment::setSpatialDimensions (double value)
{
  bool representsInteger = true;
  if (floor(value) != value)
    representsInteger = false;

  switch (getLevel())
  {
  case 1:
    /* level 1 spatialDimensions was not an attribute */
    mSpatialDimensions = 3;
    return LIBSBML_UNEXPECTED_ATTRIBUTE;
    break;
  case 2:
    if (!representsInteger || value < 0 || value > 3)
    {
      return LIBSBML_INVALID_ATTRIBUTE_VALUE;
    }
    else
    {
      mSpatialDimensions = (int) value;
      mSpatialDimensionsDouble = value;
      mIsSetSpatialDimensions  = true;
      return LIBSBML_OPERATION_SUCCESS;
    }
    break;
  case 3:
  default:
      mSpatialDimensions = (int) value;
      mSpatialDimensionsDouble = value;
      mIsSetSpatialDimensions  = true;
      return LIBSBML_OPERATION_SUCCESS;
    break;
  }
}


/*
 * Sets the size (volume in L1) of this Compartment to value.
 */
int
Compartment::setSize (double value)
{
  /* since the setSize function has been used as an
   * alias for setVolume we cant require it to only
   * be used on a L2 model
   */
/*  if ( getLevel() < 2 )
  {
    return LIBSBML_UNEXPECTED_ATTRIBUTE;
  }
*/
  mSize      = value;
  mIsSetSize = true;
  return LIBSBML_OPERATION_SUCCESS;
}


/*
 * Sets the volume (size in L2) of this Compartment to value.
 */
int
Compartment::setVolume (double value)
{
  /* since the setVolume function has been used as an
   * alias for setSize we cant require it to only
   * be used on a L1 model
   */
/*  if ( getLevel() != 1 )
  {
    return LIBSBML_UNEXPECTED_ATTRIBUTE;
  }
*/
  return setSize(value);
}


/*
 * Sets the units of this Compartment to a copy of sid.
 */
int
Compartment::setUnits (const std::string& sid)
{
  if (!(SyntaxChecker::isValidUnitSId(sid)))
  {
    return LIBSBML_INVALID_ATTRIBUTE_VALUE;
  }
  else
  {
    mUnits = sid;
    return LIBSBML_OPERATION_SUCCESS;
  }
}


/*
 * Sets the outside of this Compartment to a copy of sid.
 */
int
Compartment::setOutside (const std::string& sid)
{
  if (!(SyntaxChecker::isValidSBMLSId(sid)))
  {
    return LIBSBML_INVALID_ATTRIBUTE_VALUE;
  }
  else
  {
    mOutside = sid;
    return LIBSBML_OPERATION_SUCCESS;
  }
}


/*
 * Sets the constant field of this Compartment to value.
 */
int
Compartment::setConstant (bool value)
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
 * Unsets the name of this SBML object.
 */
int
Compartment::unsetName ()
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
 * Unsets the compartmentType of this Compartment.
 */
int
Compartment::unsetCompartmentType ()
{
  if ( (getLevel() < 2)
    || (getLevel() == 2 && getVersion() == 1))
  {
    mCompartmentType.erase();
    return LIBSBML_UNEXPECTED_ATTRIBUTE;
  }

  mCompartmentType.erase();

  if (mCompartmentType.empty()) 
  {
    return LIBSBML_OPERATION_SUCCESS;
  }
  else
  {
    return LIBSBML_OPERATION_FAILED;
  }
}


/*
 * Unsets the size (volume in L1) of this Compartment.
 */
int
Compartment::unsetSize ()
{
  if (getLevel() == 1) 
  {
    mSize = 1.0;
  }
  else
  {
    mSize = numeric_limits<double>::quiet_NaN();
  }

  mIsSetSize = false;
  
  if (!isSetSize())
  {
    return LIBSBML_OPERATION_SUCCESS;
  }
  else
  {
    return LIBSBML_OPERATION_FAILED;
  }
}


/*
 * Unsets the volume (size in L2) of this Compartment.
 *
 * In SBML L1, a Compartment volume has a default value (1.0) and therefore
 * <b>should always be set</b>.  In L2, volume is optional with no default
 * value and as such may or may not be set.
 */
int
Compartment::unsetVolume ()
{
  return unsetSize();
}


/*
 * Unsets the units of this Compartment.
 */
int
Compartment::unsetUnits ()
{
  mUnits.erase();

  if (mUnits.empty()) 
  {
    return LIBSBML_OPERATION_SUCCESS;
  }
  else
  {
    return LIBSBML_OPERATION_FAILED;
  }
}


/*
 * Unsets the outside of this Compartment.
 */
int
Compartment::unsetOutside ()
{
  mOutside.erase();

  if (mOutside.empty()) 
  {
    return LIBSBML_OPERATION_SUCCESS;
  }
  else
  {
    return LIBSBML_OPERATION_FAILED;
  }
}

/*
 * Unsets the spatialDimensions of this Compartment.
 */
int
Compartment::unsetSpatialDimensions ()
{
  if (getLevel() < 3) 
  {
    mSpatialDimensions = 3;
    return LIBSBML_UNEXPECTED_ATTRIBUTE;
  }
  else
  {
    mSpatialDimensionsDouble = numeric_limits<double>::quiet_NaN();
  }

  mIsSetSpatialDimensions = false;
  
  if (!isSetSpatialDimensions())
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
  * Compartment.
  */
UnitDefinition *
Compartment::getDerivedUnitDefinition()
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
Compartment::getDerivedUnitDefinition() const
{
  return const_cast <Compartment *> (this)->getDerivedUnitDefinition();
}


/*
 * @return the SBMLTypeCode_t of this SBML object or SBML_UNKNOWN
 * (default).
 *
 * @see getElementName()
 */
SBMLTypeCode_t
Compartment::getTypeCode () const
{
  return SBML_COMPARTMENT;
}


/*
 * @return the name of this element ie "compartment".
 */
const string&
Compartment::getElementName () const
{
  static const string name = "compartment";
  return name;
}


bool 
Compartment::hasRequiredAttributes() const
{
  bool allPresent = true;

  /* required attributes for compartment: id (name in L1) 
   * constant (L3 -> )
   */

  if (!isSetId())
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
Compartment::readAttributes (const XMLAttributes& attributes)
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
  //const unsigned int version = getVersion();

  //std::vector<std::string> expectedAttributes;
  //expectedAttributes.clear();
  //expectedAttributes.push_back("name");
  //expectedAttributes.push_back("units");

  //if (level == 1)
  //{
  //  expectedAttributes.push_back("volume");
  //  expectedAttributes.push_back("outside");
  //}
  //else
  //{
  //  expectedAttributes.push_back("metaid");
  //  expectedAttributes.push_back("id");
  //  expectedAttributes.push_back("size");
  //  expectedAttributes.push_back("spatialDimensions");
  //  expectedAttributes.push_back("constant");

  //  if (!(level == 2 && version < 3))
  //  {
  //    expectedAttributes.push_back("sboTerm");
  //  }

  //  if (level == 2)
  //  {
  //    expectedAttributes.push_back("outside");
  //    
  //    if (version > 1)
  //    {
  //      expectedAttributes.push_back("compartmentType");
  //    }
  //  }
  //}


  //// check that all attributes are expected
  //for (int i = 0; i < attributes.getLength(); i++)
  //{
  //  std::vector<std::string>::const_iterator end = expectedAttributes.end();
  //  std::vector<std::string>::const_iterator begin = expectedAttributes.begin();
  //  std::string name = attributes.getName(i);
  //  if (std::find(begin, end, name) == end)
  //  {
  //    logUnknownAttribute(name, level, version, "<compartment>");
  //  }
  //}

  ////
  //// name: SName   { use="required" }  (L1v1, L1v2)
  ////   id: SId     { use="required" }  (L2v1 ->)
  ////
  //const string id = (level == 1) ? "name" : "id";
  //bool assigned = attributes.readInto(id, mId, getErrorLog(), true);
  //if (assigned && mId.size() == 0)
  //{
  //  logEmptyString(id, level, version, "<compartment>");
  //}
  //if (!SyntaxChecker::isValidSBMLSId(mId)) logError(InvalidIdSyntax);

  ////
  //// volume  { use="optional" default="1" }  (L1v1, L1v2)
  //// size    { use="optional" }              (L2v1 ->)
  ////
  //const string size = (level == 1) ? "volume" : "size";
  //mIsSetSize = attributes.readInto(size, mSize, getErrorLog(), false);

  ////
  //// units  { use="optional" }  (L1v1 ->)
  ////
  //assigned = attributes.readInto("units", mUnits, getErrorLog(), false);
  //if (assigned && mUnits.size() == 0)
  //{
  //  logEmptyString("units", level, version, "<compartment>");
  //}
  //if (!SyntaxChecker::isValidUnitSId(mUnits))
  //{
  //  logError(InvalidUnitIdSyntax);
  //}

  ////
  //// outside  { use="optional" }  (L1v1 -> L2v4)
  ////
  //if (level < 3)
  //{
  //  attributes.readInto("outside", mOutside, getErrorLog(), false);
  //}

  //if (level > 1)
  //{
  //  //
  //  // name: string  { use="optional" }  (L2v1 ->)
  //  //
  //  attributes.readInto("name", mName, getErrorLog(), false);
  // 
  //  //
  //  // spatialDimensions { maxInclusive="3" minInclusive="0" use="optional"
  //  //                     default="3" }  (L2v1 ->)
  //  // spatialDimensions { use="optional"}  (L3v1 ->)
  //  //
  //  if (level < 3)
  //  {
  //    attributes.readInto("spatialDimensions", mSpatialDimensions, 
  //                                                    getErrorLog(), false);
  //    if (mSpatialDimensions < 0 || mSpatialDimensions > 3)
  //    {
  //      std::string message = "The spatialDimensions attribute on ";
  //      message += "a <compartment> may only have values 0, 1, 2 or 3.";
  //      getErrorLog()->logError(NotSchemaConformant, level, version,
  //                                                            message);
  //    }
  //    else
  //    {
  //      // keep record as double
  //      mSpatialDimensionsDouble = (double)(mSpatialDimensions);
  //      mIsSetSpatialDimensions = true;
  //    }
  //  }
  //  else
  //  {
  //    mIsSetSpatialDimensions = attributes.readInto("spatialDimensions", 
  //                        mSpatialDimensionsDouble, getErrorLog(), false);
  //  }
  //  
  //  //
  //  // constant  { use="optional" default="true" }  (L2v1 ->)
  //  // constant  { use="required" }  (L3v1 ->)
  //  //
  //  if (level < 3)
  //  {
  //    attributes.readInto("constant", mConstant, getErrorLog(), false);
  //  }
  //  else
  //  {
  //    mIsSetConstant = attributes.readInto("constant", mConstant, 
  //                                          getErrorLog(), true);
  //  }

  //  //
  //  // compartmentType: SId  { use="optional" }  (L2v2 -> L2v4)
  //  //
  //  if ( level == 2 && version != 1)
  //  {
  //    attributes.readInto("compartmentType", mCompartmentType, 
  //                                       getErrorLog(), false);
  //  }

  //  //
  //  // sboTerm: SBOTerm { use="optional" }  (L2v3 ->)
  //  //
  //  if (!(level == 2 && version < 3)) 
  //  {
  //    mSBOTerm = SBO::readTerm(attributes, this->getErrorLog(), level, version);
  //  }
  //}
}
/** @endcond */


/** @cond doxygen-libsbml-internal */
/*
 * Subclasses should override this method to read values from the given
 * XMLAttributes set into their specific fields.  Be sure to call your
 * parents implementation of this method as well.
 */
void
Compartment::readL1Attributes (const XMLAttributes& attributes)
{
  const unsigned int level = 1;
  const unsigned int version = getVersion();

  std::vector<std::string> expectedAttributes;
  expectedAttributes.clear();
  expectedAttributes.push_back("name");
  expectedAttributes.push_back("units");
  expectedAttributes.push_back("volume");
  expectedAttributes.push_back("outside");

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
        logUnknownAttribute(name, level, version, "<compartment>");
      }
    }
  }

  //
  // name: SName   { use="required" }  (L1v1, L1v2)
  //
  bool assigned = attributes.readInto("name", mId, getErrorLog(), true);
  if (assigned && mId.size() == 0)
  {
    logEmptyString("name", level, version, "<compartment>");
  }
  if (!SyntaxChecker::isValidSBMLSId(mId)) logError(InvalidIdSyntax);

  //
  // volume  { use="optional" default="1" }  (L1v1, L1v2)
  //
  mIsSetSize = attributes.readInto("volume", mSize, getErrorLog(), false);

  //
  // units  { use="optional" }  (L1v1 ->)
  //
  assigned = attributes.readInto("units", mUnits, getErrorLog(), false);
  if (assigned && mUnits.size() == 0)
  {
    logEmptyString("units", level, version, "<compartment>");
  }
  if (!SyntaxChecker::isValidUnitSId(mUnits))
  {
    logError(InvalidUnitIdSyntax);
  }

  //
  // outside  { use="optional" }  (L1v1 -> L2v4)
  //
  attributes.readInto("outside", mOutside, getErrorLog(), false);
}
/** @endcond */


/** @cond doxygen-libsbml-internal */
/*
 * Subclasses should override this method to read values from the given
 * XMLAttributes set into their specific fields.  Be sure to call your
 * parents implementation of this method as well.
 */
void
Compartment::readL2Attributes (const XMLAttributes& attributes)
{
  const unsigned int level = 2;
  const unsigned int version = getVersion();

  std::vector<std::string> expectedAttributes;
  expectedAttributes.clear();
  expectedAttributes.push_back("name");
  expectedAttributes.push_back("units");
  expectedAttributes.push_back("metaid");
  expectedAttributes.push_back("id");
  expectedAttributes.push_back("size");
  expectedAttributes.push_back("spatialDimensions");
  expectedAttributes.push_back("constant");
  expectedAttributes.push_back("outside");

  if (version > 1)
  {
    expectedAttributes.push_back("compartmentType");
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
        logUnknownAttribute(name, level, version, "<compartment>");
      }
    }
  }

  //
  //   id: SId     { use="required" }  (L2v1 ->)
  //
  bool assigned = attributes.readInto("id", mId, getErrorLog(), true);
  if (assigned && mId.size() == 0)
  {
    logEmptyString("id", level, version, "<compartment>");
  }
  if (!SyntaxChecker::isValidSBMLSId(mId)) logError(InvalidIdSyntax);

  //
  // size    { use="optional" }              (L2v1 ->)
  //
  mIsSetSize = attributes.readInto("size", mSize, getErrorLog(), false);

  //
  // units  { use="optional" }  (L1v1 ->)
  //
  assigned = attributes.readInto("units", mUnits, getErrorLog(), false);
  if (assigned && mUnits.size() == 0)
  {
    logEmptyString("units", level, version, "<compartment>");
  }
  if (!SyntaxChecker::isValidUnitSId(mUnits))
  {
    logError(InvalidUnitIdSyntax);
  }

  //
  // outside  { use="optional" }  (L1v1 -> L2v4)
  //
  attributes.readInto("outside", mOutside, getErrorLog(), false);

  //
  // name: string  { use="optional" }  (L2v1 ->)
  //
  attributes.readInto("name", mName, getErrorLog(), false);
  
  //
  // spatialDimensions { maxInclusive="3" minInclusive="0" use="optional"
  //                     default="3" }  (L2v1 ->)
  attributes.readInto("spatialDimensions", mSpatialDimensions, 
                                                  getErrorLog(), false);
  if (mSpatialDimensions < 0 || mSpatialDimensions > 3)
  {
    std::string message = "The spatialDimensions attribute on ";
    message += "a <compartment> may only have values 0, 1, 2 or 3.";
    getErrorLog()->logError(NotSchemaConformant, level, version,
                                                          message);
  }
  else
  {
    // keep record as double
    mSpatialDimensionsDouble = (double)(mSpatialDimensions);
    mIsSetSpatialDimensions = true;
  }

  //
  // constant  { use="optional" default="true" }  (L2v1 ->)
  //
  attributes.readInto("constant", mConstant, getErrorLog(), false);

  //
  // compartmentType: SId  { use="optional" }  (L2v2 -> L2v4)
  //
  if (version != 1)
  {
    attributes.readInto("compartmentType", mCompartmentType, 
                                        getErrorLog(), false);
  }

  //
  // sboTerm: SBOTerm { use="optional" }  (L2v3 ->)
  //
  if (version > 2) 
  {
    mSBOTerm = SBO::readTerm(attributes, this->getErrorLog(), level, version);
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
Compartment::readL3Attributes (const XMLAttributes& attributes)
{
  const unsigned int level = 3;
  const unsigned int version = getVersion();

  std::vector<std::string> expectedAttributes;
  expectedAttributes.clear();
  expectedAttributes.push_back("name");
  expectedAttributes.push_back("units");
  expectedAttributes.push_back("metaid");
  expectedAttributes.push_back("id");
  expectedAttributes.push_back("size");
  expectedAttributes.push_back("spatialDimensions");
  expectedAttributes.push_back("constant");
  expectedAttributes.push_back("sboTerm");


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
        logUnknownAttribute(name, level, version, "<compartment>");
      }
    }
  }

  //
  //   id: SId     { use="required" }  (L2v1 ->)
  //
  bool assigned = attributes.readInto("id", mId, getErrorLog());
  if (!assigned)
  {
    getErrorLog()->logError(AllowedAttributesOnCompartment, level, version);
  }
  if (assigned && mId.size() == 0)
  {
    logEmptyString("id", level, version, "<compartment>");
  }
  if (!SyntaxChecker::isValidSBMLSId(mId)) logError(InvalidIdSyntax);

  //
  // size    { use="optional" }              (L2v1 ->)
  //
  mIsSetSize = attributes.readInto("size", mSize, getErrorLog(), false);

  //
  // units  { use="optional" }  (L1v1 ->)
  //
  assigned = attributes.readInto("units", mUnits, getErrorLog(), false);
  if (assigned && mUnits.size() == 0)
  {
    logEmptyString("units", level, version, "<compartment>");
  }
  if (!SyntaxChecker::isValidUnitSId(mUnits))
  {
    logError(InvalidUnitIdSyntax);
  }


  //
  // name: string  { use="optional" }  (L2v1 ->)
  //
  attributes.readInto("name", mName, getErrorLog(), false);
   
  //
  // spatialDimensions { use="optional"}  (L3v1 ->)
  //
  mIsSetSpatialDimensions = attributes.readInto("spatialDimensions", 
                        mSpatialDimensionsDouble, getErrorLog(), false);
    
  //
  // constant  { use="required" }  (L3v1 ->)
  //
  mIsSetConstant = attributes.readInto("constant", mConstant, 
                                          getErrorLog());
  if (!mIsSetConstant)
  {
    logError(AllowedAttributesOnCompartment, level, version);
  }

  mSBOTerm = SBO::readTerm(attributes, this->getErrorLog(), level, version);
}
/** @endcond */


/** @cond doxygen-libsbml-internal */
/*
 * Subclasses should override this method to write their XML attributes
 * to the XMLOutputStream.  Be sure to call your parents implementation
 * of this method as well.
 */
void
Compartment::writeAttributes (XMLOutputStream& stream) const
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
    // compartmentType: SId  { use="optional" }  (L2v2 -> L2v4)
    //
    if (level == 2 && version > 1)
    {
      stream.writeAttribute("compartmentType", mCompartmentType);
    }

    //
    // spatialDimensions { maxInclusive="3" minInclusive="0" use="optional"
    //                     default="3" }  (L2v1->L2v4)
    // spatialDimensions { use="optional"}  (L3v1 ->)
    //
    if (level == 2)
    {
      unsigned int sd = mSpatialDimensions;
      if (sd >= 0 && sd <= 2)
      {
        stream.writeAttribute("spatialDimensions", sd);
      }
    }
    else
    {
      if (isSetSpatialDimensions())
      {
        stream.writeAttribute("spatialDimensions", mSpatialDimensionsDouble);
      }
    }
  }

  //
  // volume  { use="optional" default="1" }  (L1v1, L1v2)
  // size    { use="optional" }              (L2v1->)
  //
  if (mIsSetSize)
  {
    const string size = (level == 1) ? "volume" : "size";
    stream.writeAttribute(size, mSize);
  }

  //
  // units  { use="optional" }  (L1v1, L1v2, L2v1->)
  //
  stream.writeAttribute("units", mUnits);

  //
  // outside  { use="optional" }  (L1v1-> L2v4)
  //
  if (level < 3)
  {
    stream.writeAttribute("outside", mOutside);
  }

  if (level > 1)
  {
    //
    // constant  { use="optional" default="true" }  (L2v1->)
    // constant  { use="required" }  (L3v1 ->)
    //
    if (level == 2)
    {
      if (mConstant != true)
      {
        stream.writeAttribute("constant", mConstant);
      }
    }
    else
    {
      stream.writeAttribute("constant", mConstant);
    }
    //
    // sboTerm: SBOTerm { use="optional" }  (L2v3 ->)
    //
    if (!(level == 2 && version < 3)) 
    {
      SBO::writeTerm(stream, mSBOTerm);
    }
  }
}
/** @endcond */


/*
 * @return a (deep) copy of this ListOfCompartments.
 */
ListOfCompartments*
ListOfCompartments::clone () const
{
  return new ListOfCompartments(*this);
}


/*
 * @return the SBMLTypeCode_t of SBML objects contained in this ListOf or
 * SBML_UNKNOWN (default).
 */
SBMLTypeCode_t
ListOfCompartments::getItemTypeCode () const
{
  return SBML_COMPARTMENT;
}


/*
 * @return the name of this element ie "listOfCompartments".
 */
const string&
ListOfCompartments::getElementName () const
{
  static const string name = "listOfCompartments";
  return name;
}

/* return nth item in list */
Compartment *
ListOfCompartments::get(unsigned int n)
{
  return static_cast<Compartment*>(ListOf::get(n));
}


/* return nth item in list */
const Compartment *
ListOfCompartments::get(unsigned int n) const
{
  return static_cast<const Compartment*>(ListOf::get(n));
}


/**
 * Used by ListOf::get() to lookup an SBase based by its id.
 */
struct IdEqComp : public unary_function<SBase*, bool>
{
  const string& id;

  IdEqComp (const string& id) : id(id) { }
  bool operator() (SBase* sb) 
       { return static_cast <Compartment *> (sb)->getId() == id; }
};


/* return item by id */
Compartment*
ListOfCompartments::get (const std::string& sid)
{
  return const_cast<Compartment*>( 
    static_cast<const ListOfCompartments&>(*this).get(sid) );
}


/* return item by id */
const Compartment*
ListOfCompartments::get (const std::string& sid) const
{
  vector<SBase*>::const_iterator result;

  result = find_if( mItems.begin(), mItems.end(), IdEqComp(sid) );
  return (result == mItems.end()) ? 0 : static_cast <Compartment*> (*result);
}


/* Removes the nth item from this list */
Compartment*
ListOfCompartments::remove (unsigned int n)
{
  return static_cast<Compartment*>(ListOf::remove(n));
}


/* Removes item in this list by id */
Compartment*
ListOfCompartments::remove (const std::string& sid)
{
  SBase* item = 0;
  vector<SBase*>::iterator result;

  result = find_if( mItems.begin(), mItems.end(), IdEqComp(sid) );

  if (result != mItems.end())
  {
    item = *result;
    mItems.erase(result);
  }

  return static_cast <Compartment*> (item);
}


/** @cond doxygen-libsbml-internal */
/*
 * @return the ordinal position of the element with respect to its siblings
 * or -1 (default) to indicate the position is not significant.
 */
int
ListOfCompartments::getElementPosition () const
{
  return 5;
}
/** @endcond */


/** @cond doxygen-libsbml-internal */
/*
 * @return the SBML object corresponding to next XMLToken in the
 * XMLInputStream or NULL if the token was not recognized.
 */
SBase*
ListOfCompartments::createObject (XMLInputStream& stream)
{
  const string& name   = stream.peek().getName();
  SBase*        object = 0;


  if (name == "compartment")
  {
    try
    {
      object = new Compartment(getSBMLNamespaces());
    }
    catch (SBMLConstructorException*)
    {
      object = new Compartment(SBMLDocument::getDefaultLevel(),
        SBMLDocument::getDefaultVersion());
    }
    catch ( ... )
    {
      object = new Compartment(SBMLDocument::getDefaultLevel(),
        SBMLDocument::getDefaultVersion());
    }

    if (object) mItems.push_back(object);
  }

  return object;
}
/** @endcond */


/** @cond doxygen-c-only */


/**
 * Creates a new Compartment_t structure using the given SBML @p level
 * and @p version values.
 *
 * @param level an unsigned int, the SBML Level to assign to this
 * Compartment
 *
 * @param version an unsigned int, the SBML Version to assign to this
 * Compartment
 *
 * @return a pointer to the newly created Compartment_t structure.
 *
 * @note Once a Compartment has been added to an SBMLDocument, the @p
 * level and @p version for the document @em override those used to create
 * the Compartment.  Despite this, the ability to supply the values at
 * creation time is an important aid to creating valid SBML.  Knowledge of
 * the intended SBML Level and Version  determine whether it is valid to
 * assign a particular value to an attribute, or whether it is valid to add
 * an object to an existing SBMLDocument.
 */
LIBSBML_EXTERN
Compartment_t *
Compartment_create (unsigned int level, unsigned int version)
{
  try
  {
    Compartment* obj = new Compartment(level,version);
    return obj;
  }
  catch (SBMLConstructorException)
  {
    return NULL;
  }
}


/**
 * Creates a new Compartment_t structure using the given
 * SBMLNamespaces_t structure.
 *
 * @param sbmlns SBMLNamespaces, a pointer to an SBMLNamespaces structure
 * to assign to this Compartment
 *
 * @return a pointer to the newly created Compartment_t structure.
 *
 * @note Once a Compartment has been added to an SBMLDocument, the
 * @p sbmlns namespaces for the document @em override those used to create
 * the Compartment.  Despite this, the ability to supply the values at creation 
 * time is an important aid to creating valid SBML.  Knowledge of the intended 
 * SBML Level and Version determine whether it is valid to assign a particular 
 * value to an attribute, or whether it is valid to add an object to an existing
 * SBMLDocument.
 */
LIBSBML_EXTERN
Compartment_t *
Compartment_createWithNS (SBMLNamespaces_t* sbmlns)
{
  try
  {
    Compartment* obj = new Compartment(sbmlns);
    return obj;
  }
  catch (SBMLConstructorException)
  {
    return NULL;
  }
}


/**
 * Frees the given Compartment_t structure.
 *
 * @param c the Compartment_t structure to be freed.
 */
LIBSBML_EXTERN
void
Compartment_free (Compartment_t *c)
{
  delete c;
}


/**
 * Creates a deep copy of the given Compartment_t structure
 * 
 * @param p the Compartment_t structure to be copied
 * 
 * @return a (deep) copy of the given Compartment_t structure.
 */
LIBSBML_EXTERN
Compartment_t *
Compartment_clone (const Compartment_t* c)
{
  return static_cast<Compartment*>( c->clone() );
}


/**
 * Initializes the attributes of this Compartment_t structure to their defaults.
 *
 * The exact results depends on the %SBML Level and Version in use.  The
 * cases are currently the following:
 * 
 * @li (SBML Level 1 only) sets attribute "volume" to @c 1.0
 * @li (SBML Level 2 only) sets attribute "spatialDimensions" to @c 3
 * @li (SBML Level 2 only) sets attribute "constant" to @c 1 (true)
 *
 * @param p the Compartment_t structure to initialize
 */
LIBSBML_EXTERN
void
Compartment_initDefaults (Compartment_t *c)
{
  c->initDefaults();
}


/**
 * Returns a list of XMLNamespaces_t associated with this Compartment_t
 * structure.
 *
 * @param c the Compartment_t structure
 * 
 * @return pointer to the XMLNamespaces_t structure associated with 
 * this SBML object
 */
LIBSBML_EXTERN
const XMLNamespaces_t *
Compartment_getNamespaces(Compartment_t *c)
{
  return c->getNamespaces();
}

/**
 * Takes a Compartment_t structure and returns its identifier.
 *
 * @param c the Compartment_t structure whose identifier is sought
 * 
 * @return the identifier of the Compartment_t structure @p c, as a pointer
 * to a string.
 */
LIBSBML_EXTERN
const char *
Compartment_getId (const Compartment_t *c)
{
  return c->isSetId() ? c->getId().c_str() : NULL;
}


/**
 * Takes a Compartment_t structure and returns its name.
 *
 * @param c the Compartment_t whose name is sought.
 *
 * @return the name of the Compartment_t structure @p c, as a pointer to a
 * string.
 */
LIBSBML_EXTERN
const char *
Compartment_getName (const Compartment_t *c)
{
  return c->isSetName() ? c->getName().c_str() : NULL;
}


/**
 * Get the compartment type of this Compartment, as indicated by the
 * Compartment_t structure's "compartmentType" attribute.
 *
 * @param c the Compartment_t structure
 * 
 * @return the value of the "compartmentType" attribute of the
 * Compartment_t structure @p c as a string.
 */
LIBSBML_EXTERN
const char *
Compartment_getCompartmentType (const Compartment_t *c)
{
  return c->isSetCompartmentType() ? c->getCompartmentType().c_str() : NULL;
}


/**
 * Get the number of spatial dimensions of this Compartment_t structure.
 *
 * @param c the Compartment_t structure
 * 
 * @return the value of the "spatialDimensions" attribute of the
 * Compartment_t structure @p c as an unsigned integer
 */
LIBSBML_EXTERN
unsigned int
Compartment_getSpatialDimensions (const Compartment_t *c)
{
  return c->getSpatialDimensions();
}


/**
 * Get the number of spatial dimensions of this Compartment_t structure.
 *
 * @param c the Compartment_t structure
 * 
 * @return the value of the "spatialDimensions" attribute of the
 * Compartment_t structure @p c as a double
 */
LIBSBML_EXTERN
double
Compartment_getSpatialDimensionsAsDouble (const Compartment_t *c)
{
  return c->getSpatialDimensionsAsDouble();
}


/**
 * Get the size of this Compartment.
 *
 * This method is identical to Compartment_getVolume().  In SBML Level 1,
 * compartments are always three-dimensional constructs and only have
 * volumes, whereas in SBML Level 2, compartments may be other than
 * three-dimensional and therefore the "volume" attribute is named "size"
 * in Level 2.  LibSBML provides both Compartment_getSize() and
 * Compartment_getVolume() for easier compatibility between SBML Levels.
 *
 * @param c the Compartment_t structure
 *
 * @return the value of the "size" attribute ("volume" in Level 1) of
 * the Compartment_t structure @p c as a float-point number.
 *
 * @see Compartment_isSetSize()
 */
LIBSBML_EXTERN
double
Compartment_getSize (const Compartment_t *c)
{
  return c->getSize();
}


/**
 * (For SBML Level 1) Get the volume of this Compartment
 * 
 * This method is identical to Compartment_getSize().  In SBML Level 1,
 * compartments are always three-dimensional constructs and only have
 * volumes, whereas in SBML Level 2, compartments may be other than
 * three-dimensional and therefore the "volume" attribute is named "size"
 * in Level 2.  LibSBML provides both Compartment_getSize() and
 * Compartment_getVolume() for easier compatibility between SBML Levels.
 *
 * @param c the Compartment_t structure
 *
 * @return the value of the "volume" attribute ("size" in Level 2) of
 * the Compartment_t structure @p c, as a floating-point number.
 *
 * @see Compartment_isSetVolume()
 */
LIBSBML_EXTERN
double
Compartment_getVolume (const Compartment_t *c)
{
  return c->getVolume();
}


/**
 * Get the units of this compartment's size or volume.
 *
 * @param c the Compartment_t structure
 * 
 * @return the value of the "units" attribute of the Compartment_t
 * structure @p c.
 */
LIBSBML_EXTERN
const char *
Compartment_getUnits (const Compartment_t *c)
{
  return c->isSetUnits() ? c->getUnits().c_str() : NULL;
}


/**
 * Get the identifier, if any, of the compartment that is designated
 * as being outside of this one.
 *
 * @param c the Compartment_t structure
 * 
 * @return the value of the "outside" attribute of the Compartment_t
 * structure @p c.
 */
LIBSBML_EXTERN
const char *
Compartment_getOutside (const Compartment_t *c)
{
  return c->isSetOutside() ? c->getOutside().c_str() : NULL;
}


/**
 * Get the value of the "constant" attribute of this Compartment.
 *
 * @param c the Compartment_t structure
 *
 * @return @c true if the Compartment_t structure's size is flagged as
 * being constant, @c false otherwise.
 */
LIBSBML_EXTERN
int
Compartment_getConstant (const Compartment_t *c)
{
  return static_cast<int>( c->getConstant() );
}


/**
 * Predicate indicating whether the identifier of the given Compartment_t
 * structure has been set.
 * 
 * @param c the Compartment_t structure
 * 
 * @return true (non-zero) if the "id" attribute of the Compartment_t
 * structure @p c has been set, false (0) otherwise.
 */
LIBSBML_EXTERN
int
Compartment_isSetId (const Compartment_t *c)
{
  return static_cast<int>( c->isSetId() );
}


/**
 * Predicate indicating whether the name of the given Compartment_t
 * structure has been set.
 * 
 * @param c the Compartment_t structure
 * 
 * @return true (non-zero) if the "name" attribute of the Compartment_t
 * structure @p c has been set, false (0) otherwise.
 */
LIBSBML_EXTERN
int
Compartment_isSetName (const Compartment_t *c)
{
  return static_cast<int>( c->isSetName() );
}


/**
 * Predicate returning @c true or @c false depending on whether the given
 * Compartment_t structures's "compartmentType" attribute has been set.
 *
 * @param c the Compartment_t structure
 * 
 * @return @c true if the "compartmentType" attribute of the Compartment_t
 * structure @p c has been set, @c false otherwise.
 */
LIBSBML_EXTERN
int
Compartment_isSetCompartmentType (const Compartment_t *c)
{
  return static_cast<int>( c->isSetCompartmentType() );
}


/**
 * Predicate returning @c true or @c false depending on whether the given
 * Compartment_t structure's "size" attribute has been set.
 *
 * This method is similar but not identical to Compartment_isSetVolume().
 * The latter should be used in the context of SBML Level 1 models instead
 * of Compartment_isSetSize() because Compartment_isSetVolume() performs
 * extra processing to take into account the difference in default values
 * between SBML Levels 1 and 2.
 *
 * @param c the Compartment_t structure
 * 
 * @return @c true if the "size" attribute ("volume" in Level) of the
 * Compartment_t structure @p c has been set, @c false otherwise.
 *
 * @see Compartment_isSetVolume(), Compartment_setSize()
 */
LIBSBML_EXTERN
int
Compartment_isSetSize (const Compartment_t *c)
{
  return static_cast<int>( c->isSetSize() );
}


/**
 * (For SBML Level 1) Predicate returning @c true or @c false depending on
 * whether the given Compartment_t structures's "volume" attribute has been
 * set.
 * 
 * This method is similar but not identical to Compartment_isSetSize().
 * The latter should not be used in the context of SBML Level 1 models
 * because this method (Compartment_isSetVolume()) performs extra
 * processing to take into account the difference in default values between
 * SBML Levels 1 and 2.
 *
 * @param c the Compartment_t structure
 * 
 * @return @c true if the "volume" attribute ("size" in L2) of the given
 * Compartment_t structure @p c has been set, @c false otherwise.
 *
 * @see Compartment_isSetSize(), Compartment_setVolume()
 *
 * @note In SBML Level 1, a compartment's volume has a default value (@c
 * 1.0) and therefore this method will always return @c true.  In Level
 * 2, a compartment's size (the equivalent of SBML Level 1's "volume") is
 * optional and has no default value, and therefore may or may not be
 * set.
 */
LIBSBML_EXTERN
int
Compartment_isSetVolume (const Compartment_t *c)
{
  return static_cast<int>( c->isSetVolume() );
}


/**
 * Predicate returning @c true or @c false depending on whether the given
 * Compartment_t structures's "units" attribute has been set.
 *
 * @param c the Compartment_t structure
 * 
 * @return @c true if the "units" attribute of the Compartment_t structure
 * @p c has been set, @c false otherwise.
 */
LIBSBML_EXTERN
int
Compartment_isSetUnits (const Compartment_t *c)
{
  return static_cast<int>( c->isSetUnits() );
}


/**
 * Predicate returning @c true or @c false depending on whether the given
 * Compartment_t structure's "outside" attribute has been set.
 *
 * @param c the Compartment_t structure
 * 
 * @return @c true if the "outside" attribute of the Compartment_t
 * structure @p c has been set, @c false otherwise.
 */
LIBSBML_EXTERN
int
Compartment_isSetOutside (const Compartment_t *c)
{
  return static_cast<int>( c->isSetOutside() );
}


/**
 * Predicate returning @c true or @c false depending on whether the given
 * Compartment_t structure's "spatialDimensions" attribute has been set.
 *
 * @param c the Compartment_t structure
 * 
 * @return @c true if the "spatialDimensions" attribute of the Compartment_t
 * structure @p c has been set, @c false otherwise.
 */
LIBSBML_EXTERN
int
Compartment_isSetSpatialDimensions (const Compartment_t *c)
{
  return static_cast<int>( c->isSetSpatialDimensions() );
}


/**
 * Predicate returning @c true or @c false depending on whether the given
 * Compartment_t structure's "constant" attribute has been set.
 *
 * @param c the Compartment_t structure
 * 
 * @return @c true if the "constant" attribute of the Compartment_t
 * structure @p c has been set, @c false otherwise.
 */
LIBSBML_EXTERN
int
Compartment_isSetConstant (const Compartment_t *c)
{
  return static_cast<int>( c->isSetConstant() );
}


/**
 * Sets the identifier of the given Compartment_t structure.
 *
 * This function copies the string given in @p sid.  If the string is
 * NULL, this function performs unsetId() instead.
 *
 * @param c the Compartment_t structure.
 * @param sid the identifier to which the structures "id" attribute should
 * be set.
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
Compartment_setId (Compartment_t *c, const char *sid)
{
  return (sid == NULL) ? c->setId("") : c->setId(sid);
}


/**
 * Sets the name of the given Compartment_t structure.
 *
 * This function copies the string given in @p string.  If the string is
 * NULL, this function performs unsetName() instead.
 *
 * @param c the Compartment_t structure
 *
 * @param string the identifier to which the structures "id" attribute
 * should be set.
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
Compartment_setName (Compartment_t *c, const char *name)
{
  return (name == NULL) ? c->unsetName() : c->setName(name);
}


/**
 * Sets the "compartmentType" attribute of the given Compartment_t
 * structure.
 *
 * This function copies the string given in @p string.  If the string is
 * NULL, this function performs unsetName() instead.
 *
 * @param c the Compartment_t structure
 * @param sid, the identifier of a CompartmentType object defined
 * elsewhere in the enclosing Model_t structure.
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
 * unsetting the "compartmentType" attribute.
 */
LIBSBML_EXTERN
int
Compartment_setCompartmentType (Compartment_t *c, const char *sid)
{
  return (sid == NULL) ? 
             c->unsetCompartmentType() : c->setCompartmentType(sid);
}


/**
 * Sets the "spatialDimensions" attribute of the given Compartment_t
 * structure.
 *
 * If @p value is not one of @c 0, @c 1, @c 2, or @c 3, this method will
 * have no effect (i.e., the "spatialDimensions" attribute will not be
 * set).
 * 
 *
 * @param c the Compartment_t structure
 * @param value an unsigned integer indicating the number of dimensions
 * of the given compartment.
 *
 * @return integer value indicating success/failure of the
 * function.  @if clike The value is drawn from the
 * enumeration #OperationReturnValues_t. @endif The possible values
 * returned by this function are:
 * @li LIBSBML_OPERATION_SUCCESS
 * @li LIBSBML_INVALID_ATTRIBUTE_VALUE
 * @li LIBSBML_UNEXPECTED_ATTRIBUTE
 */
LIBSBML_EXTERN
int
Compartment_setSpatialDimensions (Compartment_t *c, unsigned int value)
{
  return c->setSpatialDimensions(value);
}


/**
 * Sets the "spatialDimensions" attribute of the given Compartment_t
 * structure.
 *
 * @param c the Compartment_t structure
 * @param value a double indicating the number of dimensions
 * of the given compartment.
 *
 * @return integer value indicating success/failure of the
 * function.  @if clike The value is drawn from the
 * enumeration #OperationReturnValues_t. @endif The possible values
 * returned by this function are:
 * @li LIBSBML_OPERATION_SUCCESS
 * @li LIBSBML_INVALID_ATTRIBUTE_VALUE
 * @li LIBSBML_UNEXPECTED_ATTRIBUTE
 */
LIBSBML_EXTERN
int
Compartment_setSpatialDimensionsAsDouble (Compartment_t *c, double value)
{
  return c->setSpatialDimensions(value);
}


/**
 * Sets the "size" attribute (or "volume" in SBML Level 1) of the given
 * Compartment_t structure.
 *
 * This method is identical to Compartment_setVolume() and is provided for
 * compatibility between SBML Level 1 and Level 2.
 *
 * @param c the Compartment_t structure
 * @param value a @c double representing the size of the given
 * Compartment_t structure in whatever units are in effect
 *
 * @return integer value indicating success/failure of the
 * function.  @if clike The value is drawn from the
 * enumeration #OperationReturnValues_t. @endif The possible values
 * returned by this function are:
 * @li LIBSBML_OPERATION_SUCCESS
 */
LIBSBML_EXTERN
int
Compartment_setSize (Compartment_t *c, double value)
{
  return c->setSize(value);
}


/**
 * Sets the "volume" attribute (or "size" in SBML Level 2) of the givenq
 * Compartment_t structure.
 *
 * This method is identical to setVolume() and is provided for
 * compatibility between SBML Level 1 and Level 2.
 *
 * @param c the Compartment_t structure
 * 
 * @param value a @c double representing the volume of the given
 * Compartment_t structure in whatever units are in effect
 *
 * @return integer value indicating success/failure of the
 * function.  @if clike The value is drawn from the
 * enumeration #OperationReturnValues_t. @endif The possible values
 * returned by this function are:
 * @li LIBSBML_OPERATION_SUCCESS
 */
LIBSBML_EXTERN
int
Compartment_setVolume (Compartment_t *c, double value)
{
  return c->setVolume(value);
}


/**
 * Sets the "units" attribute of the given Compartment_t structure.
 *
 * @param c the Compartment_t structure
 * 
 * @param sid the identifier of the defined units to use.  The string will
 * be copied.
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
Compartment_setUnits (Compartment_t *c, const char *sid)
{
  return (sid == NULL) ? c->unsetUnits() : c->setUnits(sid);
}


/**
 * Sets the "outside" attribute of the given Compartment_t structure.
 *
 * @param c the Compartment_t structure
 * 
 * @param sid the identifier of a compartment that encloses this one.  The
 * string will be copied.
 *
 * @return integer value indicating success/failure of the
 * function.  @if clike The value is drawn from the
 * enumeration #OperationReturnValues_t. @endif The possible values
 * returned by this function are:
 * @li LIBSBML_OPERATION_SUCCESS
 * @li LIBSBML_INVALID_ATTRIBUTE_VALUE
 *
 * @note Using this function with an id of NULL is equivalent to
 * unsetting the "outside" attribute.
 */
LIBSBML_EXTERN
int
Compartment_setOutside (Compartment_t *c, const char *sid)
{
  return (sid == NULL) ? c->unsetOutside() : c->setOutside(sid);
}


/**
 * Sets the value of the "constant" attribute of the given Compartment_t
 * structure.
 *
 * @param c the Compartment_t structure
 * @param value an integer indicating whether the size/volume of the
 * compartment @p c should be considered constant (nonzero) or variable (zero).
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
Compartment_setConstant (Compartment_t *c, int value)
{
  return c->setConstant( static_cast<bool>(value) );
}


/**
 * Unsets the name of the given Compartment_t structure.
 *
 * @param c the Compartment_t structure
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
Compartment_unsetName (Compartment_t *c)
{
  return c->unsetName();
}


/**
 * Unsets the value of the "compartmentType" attribute of the given
 * Compartment_t structure.
 *
 * @param c the Compartment_t structure
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
Compartment_unsetCompartmentType (Compartment_t *c)
{
  return c->unsetCompartmentType();
}


/**
 * Unsets the value of the "size" attribute of the given Compartment_t
 * structure. 
 *
 * @param c the Compartment_t structure
 *
 * @return integer value indicating success/failure of the
 * function.  @if clike The value is drawn from the
 * enumeration #OperationReturnValues_t. @endif The possible values
 * returned by this function are:
 * @li LIBSBML_OPERATION_SUCCESS
 */
LIBSBML_EXTERN
int
Compartment_unsetSize (Compartment_t *c)
{
  return c->unsetSize();
}


/**
 * (For SBML Level 1) Unsets the value of the "volume" attribute of the 
 * given Compartment_t structure.
 *
 * In SBML Level 1, a Compartment_t structure's "volume" attribute has a
 * default value (1.0) and therefore <em>should always be set</em>.  In
 * Level 2, "size" is optional with no default value and as such may or may
 * not be set.
 *
 * @param c the Compartment_t structure
 *
 * @return integer value indicating success/failure of the
 * function.  @if clike The value is drawn from the
 * enumeration #OperationReturnValues_t. @endif The possible values
 * returned by this function are:
 * @li LIBSBML_OPERATION_SUCCESS
 */
LIBSBML_EXTERN
int
Compartment_unsetVolume (Compartment_t *c)
{
  return c->unsetVolume();
}


/**
 * Unsets the value of the "units" attribute of the given Compartment_t
 * structure.
 *
 * @param c the Compartment_t structure
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
Compartment_unsetUnits (Compartment_t *c)
{
  return c->unsetUnits();
}


/**
 * Unsets the value of the "outside" attribute of the given Compartment_t
 * structure.
 *
 * @param c the Compartment_t structure
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
Compartment_unsetOutside (Compartment_t *c)
{
  return c->unsetOutside();
}


/**
 * Unsets the value of the "spatialDimensions" attribute of the given Compartment_t
 * structure.
 *
 * @param c the Compartment_t structure
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
Compartment_unsetSpatialDimensions (Compartment_t *c)
{
  return c->unsetSpatialDimensions();
}


/**
 * Constructs and returns a UnitDefinition_t structure that expresses 
 * the units of this Compartment_t structure.
 *
 * @param c the Compartment_t structure whose units are to be returned.
 *
 * @return a UnitDefinition_t structure that expresses the units 
 * of this Compartment_t strucuture.
 *
 * @note This function returns the units of the Compartment_t expressed 
 * as a UnitDefinition_t. The units may be those explicitly declared 
 * or those derived from the default units of the Model_t containing
 * this Compartment_t.
 *
 * Note that the functionality that facilitates unit analysis depends 
 * on the model as a whole.  Thus, in cases where the object has not 
 * been added to a model or the model itself is incomplete,
 * unit analysis is not possible and this method will return NULL.
 *
 */
LIBSBML_EXTERN
UnitDefinition_t * 
Compartment_getDerivedUnitDefinition(Compartment_t *c)
{
  return c->getDerivedUnitDefinition();
}


/**
  * Predicate returning @c true or @c false depending on whether
  * all the required attributes for this Compartment object
  * have been set.
  *
 * @param c the Compartment_t structure to check.
 *
  * @note The required attributes for a Compartment object are:
  * @li id (name in L1)
  * @li constant (in L3 only)
  *
  * @return a true if all the required
  * attributes for this object have been defined, false otherwise.
  */
LIBSBML_EXTERN
int
Compartment_hasRequiredAttributes(Compartment_t *c)
{
  return static_cast<int>(c->hasRequiredAttributes());
}


/**
 * @return item in this ListOfCompartment with the given id or NULL if no such
 * item exists.
 */
LIBSBML_EXTERN
Compartment_t *
ListOfCompartments_getById (ListOf_t *lo, const char *sid)
{
  return (sid != NULL) ? 
    static_cast <ListOfCompartments *> (lo)->get(sid) : NULL;
}


/**
 * Removes item in this ListOf items with the given id or NULL if no such
 * item exists.  The caller owns the returned item and is responsible for
 * deleting it.
 */
LIBSBML_EXTERN
Compartment_t *
ListOfCompartments_removeById (ListOf_t *lo, const char *sid)
{
  return (sid != NULL) ? 
    static_cast <ListOfCompartments *> (lo)->remove(sid) : NULL;
}

/** @endcond */

LIBSBML_CPP_NAMESPACE_END
