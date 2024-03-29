/**
 * @file    EventAssignment.cpp
 * @brief   Implementation of EventAssignment.
 * @author  Ben Bornstein
 *
 * $Id: EventAssignment.cpp 11635 2010-08-03 04:05:22Z mhucka $
 * $HeadURL: https://sbml.svn.sourceforge.net/svnroot/sbml/trunk/libsbml/src/sbml/EventAssignment.cpp $
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

#include <sbml/math/FormulaParser.h>
#include <sbml/math/MathML.h>
#include <sbml/math/ASTNode.h>

#include <sbml/SBO.h>
#include <sbml/SBMLVisitor.h>
#include <sbml/SBMLDocument.h>
#include <sbml/SBMLError.h>
#include <sbml/Model.h>
#include <sbml/EventAssignment.h>

/** @cond doxygen-ignored */

using namespace std;

/** @endcond */

LIBSBML_CPP_NAMESPACE_BEGIN

EventAssignment::EventAssignment (unsigned int level, unsigned int version) :
   SBase ( level, version )
 , mVariable ( "" )
 , mMath     ( 0  )
{
  if (!hasValidLevelVersionNamespaceCombination())
    throw SBMLConstructorException();
}


EventAssignment::EventAssignment (SBMLNamespaces * sbmlns) :
   SBase ( sbmlns )
 , mVariable ( "" )
 , mMath     ( 0  )
{
  if (!hasValidLevelVersionNamespaceCombination())
    throw SBMLConstructorException();
}


/** @cond doxygen-libsbml-internal */

/* constructor for validators */
EventAssignment::EventAssignment() :
  SBase()
{
}

/** @endcond */
                          

/*
 * Destroys this EventAssignment.
 */
EventAssignment::~EventAssignment ()
{
  delete mMath;
}


/*
 * Copy constructor. Creates a copy of this EventAssignment.
 */
EventAssignment::EventAssignment (const EventAssignment& orig) :
   SBase   ( orig )
  , mVariable (orig.mVariable)
  , mMath   ( 0    )
{
  if (orig.mMath) 
  {
    mMath = orig.mMath->deepCopy();
    mMath->setParentSBMLObject(this);
  }
}


/*
 * Assignment operator
 */
EventAssignment& EventAssignment::operator=(const EventAssignment& rhs)
{
  if(&rhs!=this)
  {
    this->SBase::operator =(rhs);
    this->mVariable = rhs.mVariable;

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
 *
 * @return the result of calling <code>v.visit()</code>, which indicates
 * whether or not the Visitor would like to visit the Event's next
 * EventAssignment (if available).
 */
bool
EventAssignment::accept (SBMLVisitor& v) const
{
  return v.visit(*this);
}


/*
 * @return a (deep) copy of this EventAssignment.
 */
EventAssignment*
EventAssignment::clone () const
{
  return new EventAssignment(*this);
}


/*
 * @return the variable of this EventAssignment.
 */
const string&
EventAssignment::getVariable () const
{
  return mVariable;
}


/*
 * @return the math of this EventAssignment.
 */
const ASTNode*
EventAssignment::getMath () const
{
  return mMath;
}


/*
 * @return true if the variable of this EventAssignment has been set, false
 * otherwise.
 */
bool
EventAssignment::isSetVariable () const
{
  return (mVariable.empty() == false);
}


/*
 * @return true if the math of this EventAssignment has been set, false
 * otherwise.
 */
bool
EventAssignment::isSetMath () const
{
  return (mMath != 0);
}


/*
 * Sets the variable of this EventAssignment to a copy of sid.
 */
int
EventAssignment::setVariable (const std::string& sid)
{
  if (!(SyntaxChecker::isValidSBMLSId(sid)))
  {
    return LIBSBML_INVALID_ATTRIBUTE_VALUE;
  }
  else
  {
    mVariable = sid;
    return LIBSBML_OPERATION_SUCCESS;
  }
}


/*
 * Sets the math of this EventAssignment to a copy of the given ASTNode.
 */
int
EventAssignment::setMath (const ASTNode* math)
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


/*
  * Calculates and returns a UnitDefinition that expresses the units
  * returned by the math expression of this EventAssignment.
  */
UnitDefinition * 
EventAssignment::getDerivedUnitDefinition()
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
    
    std::string id = getId() + getAncestorOfType(SBML_EVENT)->getId();
    if (m->getFormulaUnitsData(id, getTypeCode()))
    {
      return m->getFormulaUnitsData(id, getTypeCode())
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
EventAssignment::getDerivedUnitDefinition() const
{
  return const_cast <EventAssignment *> (this)->getDerivedUnitDefinition();
}


/*
 * Predicate returning @c true or @c false depending on whether 
 * the math expression of this EventAssignment contains
 * parameters/numbers with undeclared units that cannot be ignored.
 */
bool 
EventAssignment::containsUndeclaredUnits()
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

    std::string id = getId() + getAncestorOfType(SBML_EVENT)->getId();
    if (m->getFormulaUnitsData(id, getTypeCode()))
    {
      return m->getFormulaUnitsData(id, getTypeCode())
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
EventAssignment::containsUndeclaredUnits() const
{
  return const_cast<EventAssignment *> (this)->containsUndeclaredUnits();
}


/*
 * @return the SBMLTypeCode_t of this SBML object or SBML_UNKNOWN
 * (default).
 *
 * @see getElementName()
 */
SBMLTypeCode_t
EventAssignment::getTypeCode () const
{
  return SBML_EVENT_ASSIGNMENT;
}


/*
 * @return the name of this element ie "eventAssignment".
 */
const string&
EventAssignment::getElementName () const
{
  static const string name = "eventAssignment";
  return name;
}


bool 
EventAssignment::hasRequiredAttributes() const
{
  bool allPresent = true;

  /* required attributes for eventAssignment: variable */

  if (!isSetVariable())
    allPresent = false;

  return allPresent;
}


bool 
EventAssignment::hasRequiredElements() const
{
  bool allPresent = true;

  /* required attributes for eventAssignment: math */

  if (!isSetMath())
    allPresent = false;

  return allPresent;
}


/** @cond doxygen-libsbml-internal */
/*
 * Subclasses should override this method to write out their contained
 * SBML objects as XML elements.  Be sure to call your parents
 * implementation of this method as well.
 */
void
EventAssignment::writeElements (XMLOutputStream& stream) const
{
  SBase::writeElements(stream);

  if (mMath) writeMathML(mMath, stream);
}
/** @endcond */


/** @cond doxygen-libsbml-internal */
/*
 * Subclasses should override this method to read (and store) XHTML,
 * MathML, etc. directly from the XMLInputStream.
 *
 * @return true if the subclass read from the stream, false otherwise.
 */
bool
EventAssignment::readOtherXML (XMLInputStream& stream)
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

    if (mMath)
    {
      if (getLevel() < 3) 
      {
        logError(NotSchemaConformant, getLevel(), getVersion(),
	        "Only one <math> element is permitted inside a "
	        "particular containing element.");
      }
      else
      {
        logError(OneMathPerEventAssignment, getLevel(), getVersion());
      }
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
EventAssignment::readAttributes (const XMLAttributes& attributes)
{
  SBase::readAttributes(attributes);

  const unsigned int level   = getLevel  ();

  switch (level)
  {
  case 1:
    logError(NotSchemaConformant, getLevel(), getVersion(),
	      "EventAssignment is not a valid component for this level/version.");
    return;
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
EventAssignment::readL2Attributes (const XMLAttributes& attributes)
{
  const unsigned int level = getLevel();
  const unsigned int version = getVersion();

  std::vector<std::string> expectedAttributes;
  expectedAttributes.clear();
  expectedAttributes.push_back("metaid");
  expectedAttributes.push_back("variable");

  if (version > 1)
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
        logUnknownAttribute(name, level, version, "<eventAssignment>");
      }
    }
  }

  //
  // variable: SId  { use="required" }  (L2v1 ->)
  //
  bool assigned = attributes.readInto("variable", mVariable, getErrorLog(), true);
  if (assigned && mVariable.size() == 0)
  {
    logEmptyString("variable", level, version, "<eventAssignment>");
  }
  if (!SyntaxChecker::isValidSBMLSId(mVariable)) logError(InvalidIdSyntax);


  //
  // sboTerm: SBOTerm { use="optional" }  (L2v2 ->)
  //
  if (version > 1)
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
EventAssignment::readL3Attributes (const XMLAttributes& attributes)
{
  const unsigned int level = getLevel();
  const unsigned int version = getVersion();

  std::vector<std::string> expectedAttributes;
  expectedAttributes.clear();
  expectedAttributes.push_back("metaid");
  expectedAttributes.push_back("variable");
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
        logUnknownAttribute(name, level, version, "<eventAssignment>");
      }
    }
  }

  //
  // variable: SId  { use="required" }  (L2v1 ->)
  //
  bool assigned = attributes.readInto("variable", mVariable);
  if (!assigned)
  {
    getErrorLog()->logError(AllowedAttributesOnEventAssignment, level, version);
  }
  if (assigned && mVariable.size() == 0)
  {
    logEmptyString("variable", level, version, "<eventAssignment>");
  }
  if (!SyntaxChecker::isValidSBMLSId(mVariable)) logError(InvalidIdSyntax);


  //
  // sboTerm: SBOTerm { use="optional" }  (L2v2 ->)
  //
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
EventAssignment::writeAttributes (XMLOutputStream& stream) const
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
  // variable: SId  { use="required" }  (L2v1 ->)
  //
  stream.writeAttribute("variable", mVariable);


  //
  // sboTerm: SBOTerm { use="optional" }  (L2v2 ->)
  //
  if (!(level == 2 && version == 1)) 
    SBO::writeTerm(stream, mSBOTerm);
}
/** @endcond */


/*
 * @return a (deep) copy of this ListOfEventAssignments.
 */
ListOfEventAssignments*
ListOfEventAssignments::clone () const
{
  return new ListOfEventAssignments(*this);
}


/*
 * @return the SBMLTypeCode_t of SBML objects contained in this ListOf or
 * SBML_UNKNOWN (default).
 */
SBMLTypeCode_t
ListOfEventAssignments::getItemTypeCode () const
{
  return SBML_EVENT_ASSIGNMENT;
}


/*
 * @return the name of this element ie "listOfEventAssignments".
 */
const string&
ListOfEventAssignments::getElementName () const
{
  static const string name = "listOfEventAssignments";
  return name;
}


/* return nth item in list */
EventAssignment *
ListOfEventAssignments::get(unsigned int n)
{
  return static_cast<EventAssignment*>(ListOf::get(n));
}


/* return nth item in list */
const EventAssignment *
ListOfEventAssignments::get(unsigned int n) const
{
  return static_cast<const EventAssignment*>(ListOf::get(n));
}


/**
 * Used by ListOf::get() to lookup an SBase based by its id.
 */
struct IdEqEA : public unary_function<SBase*, bool>
{
  const string& id;

  IdEqEA (const string& id) : id(id) { }
  bool operator() (SBase* sb) 
       { return static_cast <EventAssignment *> (sb)->getId() == id; }
};


/* return item by id */
EventAssignment*
ListOfEventAssignments::get (const std::string& sid)
{
  return const_cast<EventAssignment*>( 
    static_cast<const ListOfEventAssignments&>(*this).get(sid) );
}


/* return item by id */
const EventAssignment*
ListOfEventAssignments::get (const std::string& sid) const
{
  vector<SBase*>::const_iterator result;

  result = find_if( mItems.begin(), mItems.end(), IdEqEA(sid) );
  return (result == mItems.end()) ? 0 : 
                     static_cast <EventAssignment*> (*result);
}


/* Removes the nth item from this list */
EventAssignment*
ListOfEventAssignments::remove (unsigned int n)
{
   return static_cast<EventAssignment*>(ListOf::remove(n));
}


/* Removes item in this list by id */
EventAssignment*
ListOfEventAssignments::remove (const std::string& sid)
{
  SBase* item = 0;
  vector<SBase*>::iterator result;

  result = find_if( mItems.begin(), mItems.end(), IdEqEA(sid) );

  if (result != mItems.end())
  {
    item = *result;
    mItems.erase(result);
  }

  return static_cast <EventAssignment*> (item);
}


/** @cond doxygen-libsbml-internal */
/*
 * @return the ordinal position of the element with respect to its
 * siblings or -1 (default) to indicate the position is not significant.
 */
int
ListOfEventAssignments::getElementPosition () const
{
  return 3;
}
/** @endcond */


/** @cond doxygen-libsbml-internal */
/*
 * @return the SBML object corresponding to next XMLToken in the
 * XMLInputStream or NULL if the token was not recognized.
 */
SBase*
ListOfEventAssignments::createObject (XMLInputStream& stream)
{
  const string& name   = stream.peek().getName();
  SBase*        object = 0;


  if (name == "eventAssignment")
  {
    try
    {
      object = new EventAssignment(getSBMLNamespaces());
    }
    catch (SBMLConstructorException*)
    {
      object = new EventAssignment(SBMLDocument::getDefaultLevel(),
        SBMLDocument::getDefaultVersion());
    }
    catch ( ... )
    {
      object = new EventAssignment(SBMLDocument::getDefaultLevel(),
        SBMLDocument::getDefaultVersion());
    }
    
    if (object) mItems.push_back(object);
  }

  return object;
}
/** @endcond */



/** @cond doxygen-c-only */


/**
 * Creates a new EventAssignment_t structure using the given SBML @p level
 * and @p version values.
 *
 * @param level an unsigned int, the SBML Level to assign to this
 * EventAssignment
 *
 * @param version an unsigned int, the SBML Version to assign to this
 * EventAssignment
 *
 * @return a pointer to the newly created EventAssignment_t structure.
 *
 * @note Once a EventAssignment has been added to an SBMLDocument, the @p
 * level and @p version for the document @em override those used to create
 * the EventAssignment.  Despite this, the ability to supply the values at
 * creation time is an important aid to creating valid SBML.  Knowledge of
 * the intended SBML Level and Version  determine whether it is valid to
 * assign a particular value to an attribute, or whether it is valid to add
 * an object to an existing SBMLDocument.
 */
LIBSBML_EXTERN
EventAssignment_t *
EventAssignment_create (unsigned int level, unsigned int version)
{
  try
  {
    EventAssignment* obj = new EventAssignment(level,version);
    return obj;
  }
  catch (SBMLConstructorException)
  {
    return NULL;
  }
}


/**
 * Creates a new EventAssignment_t structure using the given
 * SBMLNamespaces_t structure.
 *
 * @param sbmlns SBMLNamespaces, a pointer to an SBMLNamespaces structure
 * to assign to this EventAssignment
 *
 * @return a pointer to the newly created EventAssignment_t structure.
 *
 * @note Once a EventAssignment has been added to an SBMLDocument, the
 * @p sbmlns namespaces for the document @em override those used to create
 * the EventAssignment.  Despite this, the ability to supply the values at 
 * creation time is an important aid to creating valid SBML.  Knowledge of the 
 * intended SBML Level and Version determine whether it is valid to assign a 
 * particular value to an attribute, or whether it is valid to add an object to 
 * an existing SBMLDocument.
 */
LIBSBML_EXTERN
EventAssignment_t *
EventAssignment_createWithNS (SBMLNamespaces_t* sbmlns)
{
  try
  {
    EventAssignment* obj = new EventAssignment(sbmlns);
    return obj;
  }
  catch (SBMLConstructorException)
  {
    return NULL;
  }
}


/**
 * Frees the given EventAssignment_t structure.
 *
 * @param ea the EventAssignment_t to be freed.
 */
LIBSBML_EXTERN
void
EventAssignment_free (EventAssignment_t *ea)
{
  delete ea;
}


/**
 * Creates a (deep) copy of the given EventAssignment_t structure.
 *
 * @param ea the EventAssignment_t to be copied
 * 
 * @return a (deep) copy of @p ea.
 */
LIBSBML_EXTERN
EventAssignment_t *
EventAssignment_clone (const EventAssignment_t *ea)
{
  return static_cast<EventAssignment*>( ea->clone() );
}


/**
 * Returns a list of XMLNamespaces_t associated with this EventAssignment_t
 * structure.
 *
 * @param ea the EventAssignment_t structure
 * 
 * @return pointer to the XMLNamespaces_t structure associated with 
 * this SBML object
 */
LIBSBML_EXTERN
const XMLNamespaces_t *
EventAssignment_getNamespaces(EventAssignment_t *ea)
{
  return ea->getNamespaces();
}


/**
 * Gets the value of the "variable" attribute of this EventAssignment_t
 * structure.
 *
 * @param ea the EventAssignment_t structure to query.
 *
 * @return the identifier stored in the "variable" attribute of @p ea.
 */
LIBSBML_EXTERN
const char *
EventAssignment_getVariable (const EventAssignment_t *ea)
{
  return ea->isSetVariable() ? ea->getVariable().c_str() : NULL;
}


/**
 * Gets the mathematical formula stored in the given EventAssignment_t
 * structure.
 *
 * @param ea the EventAssignment_t structure to query.
 *
 * @return the ASTNode tree stored in @p ea.
 */
LIBSBML_EXTERN
const ASTNode_t *
EventAssignment_getMath (const EventAssignment_t *ea)
{
  return ea->getMath();
}


/**
 * Predicate for testing whether the attribute "variable" of the
 * given EventAssignment_t structure has been set.
 *
 * @param ea the EventAssignment_t structure to query.
 * 
 * @return nonzero (for true) if the "variable" attribute of @p ea
 * has been set, zero (0) otherwise.
 */
LIBSBML_EXTERN
int
EventAssignment_isSetVariable (const EventAssignment_t *ea)
{
  return static_cast<int>( ea->isSetVariable() );
}


/**
 * Predicate for testing whether the attribute "variable" of the
 * given EventAssignment_t structure has been set.
 *
 * @param ea the EventAssignment_t structure to query.
 * 
 * @return nonzero (for true) if the "variable" attribute of @p ea
 * has been set, zero (0) otherwise.
 */
LIBSBML_EXTERN
int
EventAssignment_isSetMath (const EventAssignment_t *ea)
{
  return static_cast<int>( ea->isSetMath() );
}


/**
 * Sets the attribute "variable" of the given EventAssignment_t structure
 * to a copy of the given identifier string.
 *
 * @param ea the EventAssignment_t to set.
 * @param sid the identifier of a Compartment, Species or (global)
 * Parameter defined in this model.
 *
 * @return integer value indicating success/failure of the
 * function.  @if clike The value is drawn from the
 * enumeration #OperationReturnValues_t. @endif The possible values
 * returned by this function are:
 * @li LIBSBML_OPERATION_SUCCESS
 * @li LIBSBML_INVALID_ATTRIBUTE_VALUE
 *
 * @note Using this function with an id of NULL is equivalent to
 * unsetting the "variable" attribute.
 */
LIBSBML_EXTERN
int
EventAssignment_setVariable (EventAssignment_t *ea, const char *sid)
{
  return ea->setVariable(sid ? sid : "");
}


/**
 * Sets the "math" subelement content of the given EventAssignment_t
 * structure to the given ASTNode.
 *
 * The given @p math ASTNode is copied.
 *
 * @param ea the EventAssignment_t to set.
 * @param math the ASTNode to copy into @p ea
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
EventAssignment_setMath (EventAssignment_t *ea, const ASTNode_t *math)
{
  return ea->setMath(math);
}

/**
  * Calculates and returns a UnitDefinition_t that expresses the units
  * returned by the math expression of this EventAssignment_t.
  *
  * @return a UnitDefinition_t that expresses the units of the math 
  * expression of this EventAssignment_t.
  *
  * Note that the functionality that facilitates unit analysis depends 
  * on the model as a whole.  Thus, in cases where the object has not 
  * been added to a model or the model itself is incomplete,
  * unit analysis is not possible and this method will return NULL.
  *
  * @note The units are calculated by applying the mathematics 
  * from the expression to the units of the <ci> elements used 
  * within the expression. Where there are parameters/numbers
  * with undeclared units the UnitDefinition_t returned by this
  * function may not accurately represent the units of the expression.
  *
  * @see EventAssignment_containsUndeclaredUnits()
  */
LIBSBML_EXTERN
UnitDefinition_t * 
EventAssignment_getDerivedUnitDefinition(EventAssignment_t *ea)
{
  return ea->getDerivedUnitDefinition();
}


/**
  * Predicate returning @c true or @c false depending on whether 
  * the math expression of this EventAssignment_t contains
  * parameters/numbers with undeclared units.
  * 
  * @return @c true if the math expression of this EventAssignment_t
  * includes parameters/numbers 
  * with undeclared units, @c false otherwise.
  *
  * @note a return value of @c true indicates that the UnitDefinition_t
  * returned by the getDerivedUnitDefinition function may not 
  * accurately represent the units of the expression.
  *
  * @see EventAssignment_getDerivedUnitDefinition()
  */
LIBSBML_EXTERN
int 
EventAssignment_containsUndeclaredUnits(EventAssignment_t *ea)
{
  return static_cast<int>(ea->containsUndeclaredUnits());
}


/**
 * @return item in this ListOfEventAssignment with the given id or NULL if no such
 * item exists.
 */
LIBSBML_EXTERN
EventAssignment_t *
ListOfEventAssignments_getById (ListOf_t *lo, const char *sid)
{
  return (sid != NULL) ? 
    static_cast <ListOfEventAssignments *> (lo)->get(sid) : NULL;
}


/**
 * Removes item in this ListOf items with the given id or NULL if no such
 * item exists.  The caller owns the returned item and is responsible for
 * deleting it.
 */
LIBSBML_EXTERN
EventAssignment_t *
ListOfEventAssignments_removeById (ListOf_t *lo, const char *sid)
{
  return (sid != NULL) ? 
    static_cast <ListOfEventAssignments *> (lo)->remove(sid) : NULL;
}

/** @endcond */

LIBSBML_CPP_NAMESPACE_END

