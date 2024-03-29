/**
 * @file    StoichiometryMath.cpp
 * @brief   Implementation of StoichiometryMath.
 * @author  Sarah Keating
 *
 * $Id: StoichiometryMath.cpp 11635 2010-08-03 04:05:22Z mhucka $
 * $HeadURL: https://sbml.svn.sourceforge.net/svnroot/sbml/trunk/libsbml/src/sbml/StoichiometryMath.cpp $
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

#include <sbml/xml/XMLNode.h>
#include <sbml/xml/XMLAttributes.h>
#include <sbml/xml/XMLInputStream.h>
#include <sbml/xml/XMLOutputStream.h>

#include <sbml/math/FormulaFormatter.h>
#include <sbml/math/FormulaParser.h>
#include <sbml/math/MathML.h>
#include <sbml/math/ASTNode.h>

#include <sbml/SBO.h>
#include <sbml/SBMLVisitor.h>
#include <sbml/SBMLError.h>
#include <sbml/SBMLDocument.h>
#include <sbml/Model.h>
#include <sbml/Parameter.h>
#include <sbml/StoichiometryMath.h>

/** @cond doxygen-ignored */

using namespace std;

/** @endcond */

LIBSBML_CPP_NAMESPACE_BEGIN

StoichiometryMath::StoichiometryMath (unsigned int level, unsigned int version) :
   SBase ( level, version )
 , mMath      ( 0              )
 , mInternalId ( "" )
{
  if (!hasValidLevelVersionNamespaceCombination())
    throw SBMLConstructorException();
}


StoichiometryMath::StoichiometryMath (SBMLNamespaces * sbmlns) :
   SBase ( sbmlns )
 , mMath      ( 0              )
 , mInternalId ( "" )
{
  if (!hasValidLevelVersionNamespaceCombination())
    throw SBMLConstructorException();
}


/** @cond doxygen-libsbml-internal */

/* constructor for validators */
StoichiometryMath::StoichiometryMath() :
  SBase()
{
}

/** @endcond */
                          

/*
 * Destroys this StoichiometryMath.
 */
StoichiometryMath::~StoichiometryMath ()
{
  delete mMath;
}


/*
 * Copy constructor. Creates a copy of this StoichiometryMath.
 */
StoichiometryMath::StoichiometryMath (const StoichiometryMath& rhs) :
   SBase          ( rhs )
 , mMath          ( 0    )
 , mInternalId    ( rhs.mInternalId )
{
  if (rhs.mMath) 
  {
    mMath = rhs.mMath->deepCopy();
    mMath->setParentSBMLObject(this);
  }
}


/*
 * Assignment operator.
 */
StoichiometryMath& StoichiometryMath::operator=(const StoichiometryMath& rhs)
{
  if(&rhs!=this)
  {
    this->SBase::operator =(rhs);
    this->mInternalId = rhs.mInternalId;

    delete mMath;
    if (rhs.mMath) 
    {
      mMath = rhs.mMath->deepCopy();
      mMath->setParentSBMLObject(this);
    }
    else
    {
      mMath = 0;
    }
  }

  return *this;
}


/*
 * Accepts the given SBMLVisitor.
 */
bool
StoichiometryMath::accept (SBMLVisitor& v) const
{
  return v.visit(*this);
}


/*
 * @return a (deep) copy of this StoichiometryMath.
 */
StoichiometryMath*
StoichiometryMath::clone () const
{
  return new StoichiometryMath(*this);
}


/*
 * @return the math of this StoichiometryMath.
 */
const ASTNode*
StoichiometryMath::getMath () const
{
  return mMath;
}


/*
 * @return true if the math (or equivalently the formula) of this
 * StoichiometryMath has been set, false otherwise.
 */
bool
StoichiometryMath::isSetMath () const
{
  return (mMath != 0);
}



/*
 * Sets the math of this StoichiometryMath to a copy of the given ASTNode.
 */
int
StoichiometryMath::setMath (const ASTNode* math)
{
  if (mMath == math) 
  {
    return LIBSBML_OPERATION_SUCCESS;
  }
  else if (math == NULL)
  {
    delete mMath;
    mMath = 0;
    return LIBSBML_OPERATION_SUCCESS;
  }
  else if (!(math->isWellFormedASTNode()))
  {
    return LIBSBML_INVALID_OBJECT;
  }
  else
  {
    delete mMath;
    mMath = (math != 0) ? math->deepCopy() : 0;
    if (mMath) mMath->setParentSBMLObject(this);
    return LIBSBML_OPERATION_SUCCESS;
  }
}

/** @cond doxygen-libsbml-internal */

/*
 * Sets the parent SBMLDocument of this SBML object.
 */
void
StoichiometryMath::setSBMLDocument (SBMLDocument* d)
{
  mSBML = d;
}


/**
  * Sets the parent SBML object of this SBML object.
  *
  * @param sb the SBML object to use
  */
void 
StoichiometryMath::setParentSBMLObject (SBase* sb)
{
  mParentSBMLObject = sb;
}
/** @endcond */


/*
 * @return the SBMLTypeCode_t of this SBML object or SBML_UNKNOWN
 * (default).
 *
 * @see getElementName()
 */
SBMLTypeCode_t
StoichiometryMath::getTypeCode () const
{
  return SBML_STOICHIOMETRY_MATH;
}


/*
 * @return the name of this element ie "stoichiometryMath".
 */
const string&
StoichiometryMath::getElementName () const
{
  static const string name = "stoichiometryMath";
  return name;
}


/** @cond doxygen-libsbml-internal */
/*
 * @return the ordinal position of the element with respect to its siblings
 * or -1 (default) to indicate the position is not significant.
 */
int
StoichiometryMath::getElementPosition () const
{
  return 0;
}
/** @endcond */


bool 
StoichiometryMath::hasRequiredElements() const
{
  bool allPresent = true;

  /* required attributes for stoichiometryMath: math */

  if (!isSetMath())
    allPresent = false;

  return allPresent;
}


/** @cond doxygen-libsbml-internal */
/*
 * Subclasses should override this method to read (and store) XHTML,
 * MathML, etc. directly from the XMLInputStream.
 *
 * @return true if the subclass read from the stream, false otherwise.
 */
bool
StoichiometryMath::readOtherXML (XMLInputStream& stream)
{
  bool          read = false;
  const string& name = stream.peek().getName();

  if (name == "math")
  {
    // if this is level 1 there shouldnt be any math!!!
    if (getLevel() == 1) 
    {
      logError(NotSchemaConformant, getLevel(), getVersion(),
	       "SBML Level 1 does not support MathML.");
      delete mMath;
      return false;
    }

    /* check for MathML namespace 
     * this may be explicitly declared here
     * or implicitly declared on the whole document
     */
    const XMLToken elem = stream.peek();
    unsigned int match = 0;
    int n;
    if (elem.getNamespaces().getLength() != 0)
    {
      for (n = 0; n < elem.getNamespaces().getLength(); n++)
      {
        if (!strcmp(elem.getNamespaces().getURI(n).c_str(), "http://www.w3.org/1998/Math/MathML"))
        {
          match = 1;
          break;
        }
      }
    }
    if (match == 0)
    {
      if( mSBML->getNamespaces() != NULL)
      /* check for implicit declaration */
      {
        for (n = 0; n < mSBML->getNamespaces()->getLength(); n++)
        {
          if (!strcmp(mSBML->getNamespaces()->getURI(n).c_str(), 
                                                     "http://www.w3.org/1998/Math/MathML"))
          {
            match = 1;
            break;
          }
        }
      }
    }
    if (match == 0)
    {
      logError(InvalidMathElement);
    }
    delete mMath;
    mMath = readMathML(stream);
    if (mMath) mMath->setParentSBMLObject(this);
    read  = true;
  }

  return read;
}
/** @endcond */


/** @cond doxygen-libsbml-internal */
/*
 * Subclasses should override this method to read values from the given
 * XMLAttributes set into their specific fields.  Be sure to call your
 * parents implementation of this method as well.
 */
void
StoichiometryMath::readAttributes (const XMLAttributes& attributes)
{
  SBase::readAttributes(attributes);
  const unsigned int level = getLevel();
  const unsigned int version = getVersion();

  if (level < 2 )
  {
    logError(NotSchemaConformant, getLevel(), getVersion(),
	      "StoichiometryMath is not a valid component for this level/version.");
    return;
  }
  std::vector<std::string> expectedAttributes;
  expectedAttributes.clear();
  expectedAttributes.push_back("metaid");

  if (!(level == 2 && version < 3))
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
        logUnknownAttribute(name, level, version, "<stoichiometryMath>");
      }
    }
  }

  //
  // sboTerm: SBOTerm { use="optional" }  (L2v3->)
  //
  if (!(level == 2 && version < 3))
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
StoichiometryMath::writeAttributes (XMLOutputStream& stream) const
{
  SBase::writeAttributes(stream);
  const unsigned int level = getLevel();
  const unsigned int version = getVersion();

  /* invalid level/version */
  if (level < 2)
  {
    return;
  }
  //
  // sboTerm: SBOTerm { use="optional" }  (L2v3->)
  //
  if (!(level == 2 && version < 3))
    SBO::writeTerm(stream, mSBOTerm);
}
/** @endcond */


/*
  * Calculates and returns a UnitDefinition that expresses the units
  * returned by the math expression of this StoichiometryMath.
  */
UnitDefinition * 
StoichiometryMath::getDerivedUnitDefinition()
{
  if (!isSetMath())
    return NULL;
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
    
    if (m->getFormulaUnitsData(getInternalId(), getTypeCode()))
    {
      return m->getFormulaUnitsData(getInternalId(), getTypeCode())
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
StoichiometryMath::getDerivedUnitDefinition() const
{
  return const_cast <StoichiometryMath *> (this)->getDerivedUnitDefinition();
}


/*
 * Predicate returning @c true or @c false depending on whether 
 * the math expression of this StoichiometryMath contains
 * parameters/numbers with undeclared units that cannot be ignored.
 */
bool 
StoichiometryMath::containsUndeclaredUnits()
{
  if (!isSetMath())
    return false;
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
    
    if (m->getFormulaUnitsData(getInternalId(), getTypeCode()))
    {
      return m->getFormulaUnitsData(getInternalId(), getTypeCode())
      ->getContainsUndeclaredUnits();
    }
    else
    {
      return false;
    }  
  }
  else
  {
    return false;
  }
}


bool 
StoichiometryMath::containsUndeclaredUnits() const
{
  return const_cast<StoichiometryMath *> (this)->containsUndeclaredUnits();
}


/** @cond doxygen-libsbml-internal */
/*
 * Subclasses should override this method to write out their contained
 * SBML objects as XML elements.  Be sure to call your parents
 * implementation of this method as well.
 */
void
StoichiometryMath::writeElements (XMLOutputStream& stream) const
{
  SBase::writeElements(stream);

  if ( getLevel() == 2 && isSetMath() ) writeMathML(getMath(), stream);
}
/** @endcond */



/**
 * Creates a new StoichiometryMath_t structure using the given SBML @p level
 * and @p version values.
 *
 * @param level an unsigned int, the SBML Level to assign to this
 * StoichiometryMath
 *
 * @param version an unsigned int, the SBML Version to assign to this
 * StoichiometryMath
 *
 * @return a pointer to the newly created StoichiometryMath_t structure.
 *
 * @note Once a StoichiometryMath has been added to an SBMLDocument, the @p
 * level and @p version for the document @em override those used to create
 * the StoichiometryMath.  Despite this, the ability to supply the values at
 * creation time is an important aid to creating valid SBML.  Knowledge of
 * the intended SBML Level and Version  determine whether it is valid to
 * assign a particular value to an attribute, or whether it is valid to add
 * an object to an existing SBMLDocument.
 */
LIBSBML_EXTERN
StoichiometryMath_t *
StoichiometryMath_create (unsigned int level, unsigned int version)
{
  try
  {
    StoichiometryMath* obj = new StoichiometryMath(level,version);
    return obj;
  }
  catch (SBMLConstructorException)
  {
    return NULL;
  }
}


/**
 * Creates a new StoichiometryMath_t structure using the given
 * SBMLNamespaces_t structure.
 *
 * @param sbmlns SBMLNamespaces, a pointer to an SBMLNamespaces structure
 * to assign to this StoichiometryMath
 *
 * @return a pointer to the newly created StoichiometryMath_t structure.
 *
 * @note Once a StoichiometryMath has been added to an SBMLDocument, the
 * @p sbmlns namespaces for the document @em override those used to create
 * the StoichiometryMath.  Despite this, the ability to supply the values at 
 * creation time is an important aid to creating valid SBML.  Knowledge of the
 * intended SBML Level and Version determine whether it is valid to assign a 
 * particular value to an attribute, or whether it is valid to add an object 
 * to an existing SBMLDocument.
 */
LIBSBML_EXTERN
StoichiometryMath_t *
StoichiometryMath_createWithNS (SBMLNamespaces_t* sbmlns)
{
  try
  {
    StoichiometryMath* obj = new StoichiometryMath(sbmlns);
    return obj;
  }
  catch (SBMLConstructorException)
  {
    return NULL;
  }
}


/**
 * Frees the given StoichiometryMath.
 */
LIBSBML_EXTERN
void
StoichiometryMath_free (StoichiometryMath_t *stoichMath)
{
  delete stoichMath;
}


/**
 * @return a (deep) copy of this StoichiometryMath.
 */
LIBSBML_EXTERN
StoichiometryMath_t *
StoichiometryMath_clone (const StoichiometryMath_t *stoichMath)
{
  return stoichMath->clone();
}


/**
 * Returns a list of XMLNamespaces_t associated with this StoichiometryMath_t
 * structure.
 *
 * @param sm the StoichiometryMath_t structure
 * 
 * @return pointer to the XMLNamespaces_t structure associated with 
 * this SBML object
 */
LIBSBML_EXTERN
const XMLNamespaces_t *
StoichiometryMath_getNamespaces(StoichiometryMath_t *sm)
{
  return sm->getNamespaces();
}


/**
 * @return the stoichMath of this StoichiometryMath.
 */
LIBSBML_EXTERN
const ASTNode_t *
StoichiometryMath_getMath (const StoichiometryMath_t *stoichMath)
{
  return stoichMath->getMath();
}


/**
 * @return true (non-zero) if the stoichMath (or equivalently the formula) of
 * this StoichiometryMath has been set, false (0) otherwise.
 */
LIBSBML_EXTERN
int
StoichiometryMath_isSetMath (const StoichiometryMath_t *stoichMath)
{
  return static_cast<int>( stoichMath->isSetMath() );
}


/**
 * Sets the math of this StoichiometryMath to a copy of the given ASTNode.
 *
 * @return integer value indicating success/failure of the
 * function.  @if clike The value is drawn from the
 * enumeration #OperationReturnValues_t. @endif The possible values
 * returned by this function are:
 * @li LIBSBML_OPERATION_SUCCESS
 * @li LIBSBML_INVALID_OBJECT
 */
LIBSBML_EXTERN
int
StoichiometryMath_setMath (StoichiometryMath_t *stoichMath, const ASTNode_t *math)
{
  return stoichMath->setMath(math);
}


/**
 * Calculates and returns a UnitDefinition that expresses the units
 * returned by the math expression in this StoichiometryMath.
 *
 * @param stoichMath the StoichiometryMath structure to check
 *
 * @return A UnitDefinition that expresses the units of the math 
 *
 *
 * The units are calculated based on the mathematical expression in the
 * StoichiometryMath and the model quantities referenced by
 * <code>&lt;ci&gt;</code> elements used within that expression.  The
 * getDerivedUnitDefinition() method returns the calculated units.
 * 
 * @warning Note that it is possible the "math" expression in the
 * StoichiometryMath contains pure numbers or parameters with undeclared
 * units.  In those cases, it is not possible to calculate the units of
 * the overall expression without making assumptions.  LibSBML does not
 * make assumptions about the units, and getDerivedUnitDefinition() only
 * returns the units as far as it is able to determine them.  For
 * example, in an expression <em>X + Y</em>, if <em>X</em> has
 * unambiguously-defined units and <em>Y</em> does not, it will return
 * the units of <em>X</em>.  <strong>It is important that callers also
 * invoke the method</strong> containsUndeclaredUnits() <strong>to
 * determine whether this situation holds</strong>.  Callers may wish to
 * take suitable actions in those scenarios.
 *
 * @see containsUndeclaredUnits()
 */
LIBSBML_EXTERN
UnitDefinition_t * 
StoichiometryMath_getDerivedUnitDefinition(StoichiometryMath_t *stoichMath)
{
  return static_cast<StoichiometryMath*>(stoichMath)->getDerivedUnitDefinition();
}


/**
 * Predicate returning @c true or @c false depending on whether 
 * the math expression of this StoichiometryMath contains
 * parameters/numbers with undeclared units.
 *
 * @param stoichMath the StoichiometryMath structure to check
 * 
 * @return @c true if the math expression of this StoichiometryMath
 * includes parameters/numbers 
 * with undeclared units, @c false otherwise.
 *
 * @note A return value of @c true indicates that the UnitDefinition
 * returned by getDerivedUnitDefinition() may not accurately represent
 * the units of the expression.
 *
 * @see StoichiometryMath_getDerivedUnitDefinition()
 */
LIBSBML_EXTERN
int 
StoichiometryMath_containsUndeclaredUnits(StoichiometryMath_t *stoichMath)
{
  return static_cast<int>(static_cast<StoichiometryMath*>(stoichMath)
                                         ->containsUndeclaredUnits());
}


LIBSBML_CPP_NAMESPACE_END
