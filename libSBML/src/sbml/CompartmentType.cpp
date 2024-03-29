/**
 * @file    CompartmentType.cpp
 * @brief   Implementation of CompartmentType and ListOfCompartmentTypes.
 * @author  Ben Bornstein
 *
 * $Id: CompartmentType.cpp 11635 2010-08-03 04:05:22Z mhucka $
 * $HeadURL: https://sbml.svn.sourceforge.net/svnroot/sbml/trunk/libsbml/src/sbml/CompartmentType.cpp $
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

#include <sbml/SBO.h>
#include <sbml/SBMLVisitor.h>
#include <sbml/SBMLDocument.h>
#include <sbml/SBMLError.h>
#include <sbml/Model.h>
#include <sbml/CompartmentType.h>

/** @cond doxygen-ignored */

using namespace std;

/** @endcond */

LIBSBML_CPP_NAMESPACE_BEGIN

CompartmentType::CompartmentType (unsigned int level, unsigned int version) :
   SBase ( level, version )
 , mId   ( "" )
 , mName ( "" )
{
  if (!hasValidLevelVersionNamespaceCombination())
    throw SBMLConstructorException();
}


CompartmentType::CompartmentType (SBMLNamespaces* sbmlns) :
   SBase ( sbmlns )
 , mId   ( "" )
 , mName ( "" )
{
  if (!hasValidLevelVersionNamespaceCombination())
    throw SBMLConstructorException();
}

                          
/** @cond doxygen-libsbml-internal */

/* constructor for validators */
CompartmentType::CompartmentType() :
  SBase()
{
}

/** @endcond */
                          
/*
 * Destroys this CompartmentType.
 */
CompartmentType::~CompartmentType ()
{
}


/*
 * Copy constructor. Creates a copy of this CompartmentType.
 */
CompartmentType::CompartmentType(const CompartmentType& orig) :
   SBase             ( orig                    )
 , mId               ( orig.mId                )  
 , mName             ( orig.mName              )
{
}


/*
 * Assignment operator
 */
CompartmentType& CompartmentType::operator=(const CompartmentType& rhs)
{
  if(&rhs!=this)
  {
    this->SBase::operator =(rhs);
    mId = rhs.mId;
    mName = rhs.mName;
  }

  return *this;
}


/*
 * Accepts the given SBMLVisitor.
 *
 * @return the result of calling <code>v.visit()</code>, which indicates
 * whether or not the Visitor would like to visit the Model's next
 * CompartmentType (if available).
 */
bool
CompartmentType::accept (SBMLVisitor& v) const
{
  return v.visit(*this);
}


/*
 * @return a (deep) copy of this CompartmentType.
 */
CompartmentType*
CompartmentType::clone () const
{
  return new CompartmentType(*this);
}


/*
 * @return the id of this SBML object.
 */
const string&
CompartmentType::getId () const
{
  return mId;
}


/*
 * @return the name of this SBML object.
 */
const string&
CompartmentType::getName () const
{
  return (getLevel() == 1) ? mId : mName;
}


/*
 * @return true if the id of this SBML object has been set, false
 * otherwise.
 */
bool
CompartmentType::isSetId () const
{
  return (mId.empty() == false);
}


/*
 * @return true if the name of this SBML object has been set, false
 * otherwise.
 */
bool
CompartmentType::isSetName () const
{
  return (getLevel() == 1) ? (mId.empty() == false) : 
                            (mName.empty() == false);
}


/*
 * Sets the id of this SBML object to a copy of sid.
 */
int
CompartmentType::setId (const std::string& sid)
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
CompartmentType::setName (const std::string& name)
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
 * Unsets the name of this SBML object.
 */
int
CompartmentType::unsetName ()
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
 * @return the SBMLTypeCode_t of this SBML object or SBML_UNKNOWN
 * (default).
 *
 * @see getElementName()
 */
SBMLTypeCode_t
CompartmentType::getTypeCode () const
{
  return SBML_COMPARTMENT_TYPE;
}


/*
 * @return the name of this element ie "compartmentType".
 */
const string&
CompartmentType::getElementName () const
{
  static const string name = "compartmentType";
  return name;
}


bool 
CompartmentType::hasRequiredAttributes() const
{
  bool allPresent = true;

  /* required attributes for compartmentType: id */

  if (!isSetId())
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
CompartmentType::readAttributes (const XMLAttributes& attributes)
{
  SBase::readAttributes(attributes);

  const unsigned int level = getLevel();
  const unsigned int version = getVersion();

  if (level < 2 || (level == 2 && version == 1))
  {
    logError(NotSchemaConformant, getLevel(), getVersion(),
	      "CompartmentType is not a valid component for this level/version.");
    return;
  }

  std::vector<std::string> expectedAttributes;
  expectedAttributes.clear();

  expectedAttributes.push_back("name");
  expectedAttributes.push_back("id");
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
        logUnknownAttribute(name, level, version, "<compartmentType>");
      }
    }
  }

  //
  // id: SId  { use="required" }  (L2v2 ->)
  //
  bool assigned = attributes.readInto("id", mId, getErrorLog(), true);
  if (assigned && mId.size() == 0)
  {
    logEmptyString("id", level, version, "<compartmentType>");
  }
  if (!SyntaxChecker::isValidSBMLSId(mId)) logError(InvalidIdSyntax);

  //
  // name: string  { use="optional" }  (L2v2 ->)
  //
  attributes.readInto("name", mName);
  //
  // sboTerm: SBOTerm { use="optional" }  (L2v3 ->)
  //
  if (!(level == 2 && version < 3))
  {
    mSBOTerm = SBO::readTerm(attributes, this->getErrorLog(), level, version);
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
CompartmentType::writeAttributes (XMLOutputStream& stream) const
{
  SBase::writeAttributes(stream);

  const unsigned int level = getLevel();
  const unsigned int version = getVersion();

  /* invalid level/version */
  if (level < 2 || (level == 2 && version == 1))
  {
    return;
  }

  //
  // id: SId  { use="required" }  (L2v2 ->)
  //
  stream.writeAttribute("id", mId);

  //
  // name: string  { use="optional" }  (L2v2 ->)
  //
  stream.writeAttribute("name", mName);
  //
  // sboTerm: SBOTerm { use="optional" }  (L2v3 ->)
  //
  if (!(level == 2 && version < 3)) 
  {
    SBO::writeTerm(stream, mSBOTerm);
  }
}
/** @endcond */



/*
 * @return a (deep) copy of this ListOfCompartmentTypes.
 */
ListOfCompartmentTypes*
ListOfCompartmentTypes::clone () const
{
  return new ListOfCompartmentTypes(*this);
}


/*
 * @return the SBMLTypeCode_t of SBML objects contained in this ListOf or
 * SBML_UNKNOWN (default).
 */
SBMLTypeCode_t
ListOfCompartmentTypes::getItemTypeCode () const
{
  return SBML_COMPARTMENT_TYPE;
}


/*
 * @return the name of this element ie "listOfCompartmentTypes".
 */
const string&
ListOfCompartmentTypes::getElementName () const
{
  static const string name = "listOfCompartmentTypes";
  return name;
}


/* return nth item in list */
CompartmentType *
ListOfCompartmentTypes::get(unsigned int n)
{
  return static_cast<CompartmentType*>(ListOf::get(n));
}


/* return nth item in list */
const CompartmentType *
ListOfCompartmentTypes::get(unsigned int n) const
{
  return static_cast<const CompartmentType*>(ListOf::get(n));
}


/**
 * Used by ListOf::get() to lookup an SBase based by its id.
 */
struct IdEqCT : public unary_function<SBase*, bool>
{
  const string& id;

  IdEqCT (const string& id) : id(id) { }
  bool operator() (SBase* sb) 
       { return static_cast <CompartmentType *> (sb)->getId() == id; }
};


/* return item by id */
CompartmentType*
ListOfCompartmentTypes::get (const std::string& sid)
{
  return const_cast<CompartmentType*>( 
    static_cast<const ListOfCompartmentTypes&>(*this).get(sid) );
}


/* return item by id */
const CompartmentType*
ListOfCompartmentTypes::get (const std::string& sid) const
{
  vector<SBase*>::const_iterator result;

  result = find_if( mItems.begin(), mItems.end(), IdEqCT(sid) );
  return (result == mItems.end()) ? 0 : 
                     static_cast <CompartmentType*> (*result);
}


/* Removes the nth item from this list */
CompartmentType*
ListOfCompartmentTypes::remove (unsigned int n)
{
   return static_cast<CompartmentType*>(ListOf::remove(n));
}


/* Removes item in this list by id */
CompartmentType*
ListOfCompartmentTypes::remove (const std::string& sid)
{
  SBase* item = 0;
  vector<SBase*>::iterator result;

  result = find_if( mItems.begin(), mItems.end(), IdEqCT(sid) );

  if (result != mItems.end())
  {
    item = *result;
    mItems.erase(result);
  }

  return static_cast <CompartmentType*> (item);
}


/** @cond doxygen-libsbml-internal */
/*
 * @return the ordinal position of the element with respect to its siblings
 * or -1 (default) to indicate the position is not significant.
 */
int
ListOfCompartmentTypes::getElementPosition () const
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
ListOfCompartmentTypes::createObject (XMLInputStream& stream)
{
  const string& name   = stream.peek().getName();
  SBase*        object = 0;


  if (name == "compartmentType")
  {
    try
    {
      object = new CompartmentType(getSBMLNamespaces());
    }
    catch (SBMLConstructorException*)
    {
      object = new CompartmentType(SBMLDocument::getDefaultLevel(),
        SBMLDocument::getDefaultVersion());
    }
    catch ( ... )
    {
      object = new CompartmentType(SBMLDocument::getDefaultLevel(),
        SBMLDocument::getDefaultVersion());
    }
    
    if (object) mItems.push_back(object);
  }

  return object;
}
/** @endcond */


/** @cond doxygen-c-only */


/**
 * Creates a new CompartmentType_t structure using the given SBML @p level
 * and @p version values.
 *
 * @param level an unsigned int, the SBML Level to assign to this
 * CompartmentType
 *
 * @param version an unsigned int, the SBML Version to assign to this
 * CompartmentType
 *
 * @return a pointer to the newly created CompartmentType_t structure.
 *
 * @note Once a CompartmentType has been added to an SBMLDocument, the @p
 * level and @p version for the document @em override those used to create
 * the CompartmentType.  Despite this, the ability to supply the values at
 * creation time is an important aid to creating valid SBML.  Knowledge of
 * the intended SBML Level and Version  determine whether it is valid to
 * assign a particular value to an attribute, or whether it is valid to add
 * an object to an existing SBMLDocument.
 */
LIBSBML_EXTERN
CompartmentType_t *
CompartmentType_create (unsigned int level, unsigned int version)
{
  try
  {
    CompartmentType* obj = new CompartmentType(level,version);
    return obj;
  }
  catch (SBMLConstructorException)
  {
    return NULL;
  }
}


/**
 * Creates a new CompartmentType_t structure using the given
 * SBMLNamespaces_t structure.
 *
 * @param sbmlns SBMLNamespaces, a pointer to an SBMLNamespaces structure
 * to assign to this CompartmentType
 *
 * @return a pointer to the newly created CompartmentType_t structure.
 *
 * @note Once a CompartmentType has been added to an SBMLDocument, the
 * @p sbmlns namespaces for the document @em override those used to create
 * the CompartmentType.  Despite this, the ability to supply the values at 
 * creation time is an important aid to creating valid SBML.  Knowledge of the 
 * intended SBML Level and Version determine whether it is valid to assign a 
 * particular value to an attribute, or whether it is valid to add an object 
 * to an existing SBMLDocument.
 */
LIBSBML_EXTERN
CompartmentType_t *
CompartmentType_createWithNS (SBMLNamespaces_t* sbmlns)
{
  try
  {
    CompartmentType* obj = new CompartmentType(sbmlns);
    return obj;
  }
  catch (SBMLConstructorException)
  {
    return NULL;
  }
}


/**
 * Frees the given CompartmentType_t structure.
 *
 * @param ct the CompartmentType_t structure to be freed.
 */
LIBSBML_EXTERN
void
CompartmentType_free (CompartmentType_t *ct)
{
  delete ct;
}


/**
 * Creates a deep copy of the given CompartmentType_t structure
 * 
 * @param ct the CompartmentType_t structure to be copied
 * 
 * @return a (deep) copy of this CompartmentType_t structure.
 */
LIBSBML_EXTERN
CompartmentType_t *
CompartmentType_clone (const CompartmentType_t *ct)
{
  return static_cast<CompartmentType*>( ct->clone() );
}


/**
 * Returns a list of XMLNamespaces_t associated with this CompartmentType_t
 * structure.
 *
 * @param ct the CompartmentType_t structure
 * 
 * @return pointer to the XMLNamespaces_t structure associated with 
 * this SBML object
 */
LIBSBML_EXTERN
const XMLNamespaces_t *
CompartmentType_getNamespaces(CompartmentType_t *ct)
{
  return ct->getNamespaces();
}


/**
 * Takes a CompartmentType_t structure and returns its identifier.
 *
 * @param ct the CompartmentType_t structure whose identifier is sought
 * 
 * @return the identifier of this CompartmentType_t, as a pointer to a string.
 */
LIBSBML_EXTERN
const char *
CompartmentType_getId (const CompartmentType_t *ct)
{
  return ct->isSetId() ? ct->getId().c_str() : NULL;
}


/**
 * Takes a CompartmentType_t structure and returns its name.
 *
 * @param ct the CompartmentType_t whose name is sought.
 *
 * @return the name of this CompartmentType_t, as a pointer to a string.
 */
LIBSBML_EXTERN
const char *
CompartmentType_getName (const CompartmentType_t *ct)
{
  return ct->isSetName() ? ct->getName().c_str() : NULL;
}


/**
 * Predicate returning @c true or @c false depending on whether the given
 * CompartmentType_t structure's identifier has been set.
 *
 * @param ct the CompartmentType_t structure to query
 * 
 * @return @c non-zero (true) if the "id" field of the given
 * CompartmentType has been set, zero (false) otherwise.
 */
LIBSBML_EXTERN
int
CompartmentType_isSetId (const CompartmentType_t *ct)
{
  return static_cast<int>( ct->isSetId() );
}


/**
 * Predicate returning @c true or @c false depending on whether the given
 * CompartmentType_t structure's name has been set.
 *
 * @param ct the CompartmentType_t structure to query
 * 
 * @return @c non-zero (true) if the "name" field of the given
 * CompartmentType has been set, zero (false) otherwise.
 */
LIBSBML_EXTERN
int
CompartmentType_isSetName (const CompartmentType_t *ct)
{
  return static_cast<int>( ct->isSetName() );
}


/**
 * Assigns the identifier of a CompartmentType_t structure.
 *
 * This makes a copy of the string passed as the argument @p sid.
 *
 * @param ct the CompartmentType_t structure to set.
 * @param sid the string to use as the identifier.
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
CompartmentType_setId (CompartmentType_t *ct, const char *sid)
{
  return (sid == NULL) ? ct->setId("") : ct->setId(sid);
}


/**
 * Assign the name of a CompartmentType_t structure.
 *
 * This makes a copy of the string passed as the argument @p name.
 *
 * @param ct the CompartmentType_t structure to set.
 * @param name the string to use as the name.
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
CompartmentType_setName (CompartmentType_t *ct, const char *name)
{
  return (name == NULL) ? ct->unsetName() : ct->setName(name);
}


/**
 * Unsets the name of a CompartmentType.
 * 
 * @param ct the CompartmentType_t structure whose name is to be unset.
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
CompartmentType_unsetName (CompartmentType_t *ct)
{
  return ct->unsetName();
}


/**
 * @return item in this ListOfCompartmentType with the given id or NULL if no such
 * item exists.
 */
LIBSBML_EXTERN
CompartmentType_t *
ListOfCompartmentTypes_getById (ListOf_t *lo, const char *sid)
{
  return (sid != NULL) ? 
    static_cast <ListOfCompartmentTypes *> (lo)->get(sid) : NULL;
}


/**
 * Removes item in this ListOf items with the given id or NULL if no such
 * item exists.  The caller owns the returned item and is responsible for
 * deleting it.
 */
LIBSBML_EXTERN
CompartmentType_t *
ListOfCompartmentTypes_removeById (ListOf_t *lo, const char *sid)
{
  return (sid != NULL) ? 
    static_cast <ListOfCompartmentTypes *> (lo)->remove(sid) : NULL;
}

/** @endcond */

LIBSBML_CPP_NAMESPACE_END
