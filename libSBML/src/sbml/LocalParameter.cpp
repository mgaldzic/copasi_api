/**
 * @file    LocalParameter.cpp
 * @brief   Implementations of LocalLocalParameter and ListOfLocalLocalParameters.
 * @author  Sarah Keating
 *
 * $Id: LocalParameter.cpp 11635 2010-08-03 04:05:22Z mhucka $
 * $HeadURL: https://sbml.svn.sourceforge.net/svnroot/sbml/trunk/libsbml/src/sbml/LocalParameter.cpp $
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
 * in the file named "LICENSE.txt" included with this software distribution and
 * also available online as http://sbml.org/software/libsbml/license.html
 *----------------------------------------------------------------------- -->*/

#include <limits>

#include <sbml/xml/XMLNode.h>
#include <sbml/xml/XMLAttributes.h>
#include <sbml/xml/XMLInputStream.h>
#include <sbml/xml/XMLOutputStream.h>

#include <sbml/SBO.h>
#include <sbml/SBMLVisitor.h>
#include <sbml/SBMLDocument.h>
#include <sbml/SBMLError.h>
#include <sbml/Model.h>
#include <sbml/KineticLaw.h>
#include <sbml/LocalParameter.h>

/** @cond doxygen-ignored */

using namespace std;

/** @endcond */

LIBSBML_CPP_NAMESPACE_BEGIN

LocalParameter::LocalParameter (unsigned int level, unsigned int version) :
   Parameter ( level, version )
{
  if (!hasValidLevelVersionNamespaceCombination())
    throw SBMLConstructorException();

  // if level 3 values have no defaults
  if (level == 3)
  {
    mValue = numeric_limits<double>::quiet_NaN();
  }
}


LocalParameter::LocalParameter (SBMLNamespaces * sbmlns) :
   Parameter      ( sbmlns   )
{
  if (!hasValidLevelVersionNamespaceCombination())
    throw SBMLConstructorException();

  // if level 3 values have no defaults
  if (sbmlns->getLevel() == 3)
  {
    mValue = numeric_limits<double>::quiet_NaN();
  }
}


/** @cond doxygen-libsbml-internal */

/* constructor for validators */
LocalParameter::LocalParameter() :
  Parameter()
{
}

/** @endcond */
                          
/*
 * Destroys this LocalParameter.
 */
LocalParameter::~LocalParameter ()
{
}


/*
 * Copy constructor. Creates a copy of this LocalParameter.
 */
LocalParameter::LocalParameter(const LocalParameter& orig) :
    Parameter      ( orig             )
{
}


/*
 * Copy constructor. Creates a copy of this LocalParameter.
 */
LocalParameter::LocalParameter(const Parameter& orig) :
    Parameter      ( orig             )
{
}
/*
 * Assignment operator.
 */
LocalParameter& LocalParameter::operator=(const LocalParameter& rhs)
{
  if(&rhs!=this)
  {
    this->Parameter::operator =(rhs);
  }

  return *this;
}


/*
 * Accepts the given SBMLVisitor.
 *
 * @return the result of calling <code>v.visit()</code>, which indicates
 * whether or not the Visitor would like to visit the parent Model's or
 * KineticLaw's next LocalParameter (if available).
 */
bool
LocalParameter::accept (SBMLVisitor& v) const
{
  return v.visit(*this);
}


/*
 * @return a (deep) copy of this LocalParameter.
 */
LocalParameter*
LocalParameter::clone () const
{
  return new LocalParameter(*this);
}


/*
  * Constructs and returns a UnitDefinition that expresses the units of this 
  * LocalParameter.
  */
UnitDefinition *
LocalParameter::getDerivedUnitDefinition()
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

    UnitDefinition *ud = NULL;
    const char * units = getUnits().c_str();
    if (!strcmp(units, ""))
    {
      ud   = new UnitDefinition(getSBMLNamespaces());
      return ud;
    }
    else
    {
      if (UnitKind_isValidUnitKindString(units, 
                                getLevel(), getVersion()))
      {
        Unit * unit = new Unit(getSBMLNamespaces());
        unit->setKind(UnitKind_forName(units));
        unit->initDefaults();
        ud   = new UnitDefinition(getSBMLNamespaces());
        
        ud->addUnit(unit);

        delete unit;
      }
      else
      {
        /* must be a unit definition */
        ud = static_cast <Model *> (getAncestorOfType(SBML_MODEL))->getUnitDefinition(units);
      }
      return ud;
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
LocalParameter::getDerivedUnitDefinition() const
{
  return const_cast <LocalParameter *> (this)->getDerivedUnitDefinition();
}


/*
 * @return the SBMLTypeCode_t of this SBML object or SBML_UNKNOWN
 * (default).
 *
 * @see getElementName()
 */
SBMLTypeCode_t
LocalParameter::getTypeCode () const
{
  return SBML_LOCAL_PARAMETER;
}


/*
 * @return the name of this element ie "localparameter".
 */
const string&
LocalParameter::getElementName () const
{
  static const string name = "localParameter";
  return name;
}


bool 
LocalParameter::hasRequiredAttributes() const
{
  bool allPresent = true;

  /* required attributes for parameter: id (name in L1)
   * and value (in L1V1 only)*/

  if (!isSetId())
    allPresent = false;

  if (getLevel() == 1
    && getVersion() == 1
    && !isSetValue())
    allPresent = false;

  return allPresent;
}


/** @cond doxygen-libsbml-internal */
/*
 * Subclasses should override this method to read values from the given
 * XMLAttributes set into their specific fields.  Be sure to call your
 * parents implementation of this method as well.
 *
 * @param attributes the XMLAttributes object to use
 */
void
LocalParameter::readAttributes (const XMLAttributes& attributes)
{
  Parameter::readAttributes(attributes);
}
/** @endcond */


/** @cond doxygen-libsbml-internal */
/*
 * Subclasses should override this method to write their XML attributes
 * to the XMLOutputStream.  Be sure to call your parents implementation
 * of this method as well.
 *
 * @param stream the XMLOutputStream to use
 */
void
LocalParameter::writeAttributes (XMLOutputStream& stream) const
{
  Parameter::writeAttributes(stream);
}
/** @endcond */


/*
 * @return a (deep) copy of this ListOfLocalParameters.
 */
ListOfLocalParameters*
ListOfLocalParameters::clone () const
{
  return new ListOfLocalParameters(*this);
}


/*
 * @return the SBMLTypeCode_t of SBML objects contained in this ListOf or
 * SBML_UNKNOWN (default).
 */
SBMLTypeCode_t
ListOfLocalParameters::getItemTypeCode () const
{
  return SBML_LOCAL_PARAMETER;
}


/*
 * @return the name of this element ie "listOfLocalParameters".
 */
const string&
ListOfLocalParameters::getElementName () const
{
  static const string name = "listOfLocalParameters";
  return name;
}


/* return nth item in list */
LocalParameter *
ListOfLocalParameters::get(unsigned int n)
{
  return static_cast<LocalParameter*>(ListOf::get(n));
}


/* return nth item in list */
const LocalParameter *
ListOfLocalParameters::get(unsigned int n) const
{
  return static_cast<const LocalParameter*>(ListOf::get(n));
}


/**
 * Used by ListOf::get() to lookup an SBase based by its id.
 */
struct IdEqP : public unary_function<SBase*, bool>
{
  const string& id;

  IdEqP (const string& id) : id(id) { }
  bool operator() (SBase* sb) 
       { return static_cast <LocalParameter *> (sb)->getId() == id; }
};


/* return item by id */
LocalParameter*
ListOfLocalParameters::get (const std::string& sid)
{
  return const_cast<LocalParameter*>( 
    static_cast<const ListOfLocalParameters&>(*this).get(sid) );
}


/* return item by id */
const LocalParameter*
ListOfLocalParameters::get (const std::string& sid) const
{
  vector<SBase*>::const_iterator result;

  result = find_if( mItems.begin(), mItems.end(), IdEqP(sid) );
  return (result == mItems.end()) ? 0 : static_cast <LocalParameter*> (*result);
}


/* Removes the nth item from this list */
LocalParameter*
ListOfLocalParameters::remove (unsigned int n)
{
   return static_cast<LocalParameter*>(ListOf::remove(n));
}


/* Removes item in this list by id */
LocalParameter*
ListOfLocalParameters::remove (const std::string& sid)
{
  SBase* item = 0;
  vector<SBase*>::iterator result;

  result = find_if( mItems.begin(), mItems.end(), IdEqP(sid) );

  if (result != mItems.end())
  {
    item = *result;
    mItems.erase(result);
  }

  return static_cast <LocalParameter*> (item);
}


/** @cond doxygen-libsbml-internal */
/*
 * @return the ordinal position of the element with respect to its siblings
 * or -1 (default) to indicate the position is not significant.
 */
int
ListOfLocalParameters::getElementPosition () const
{
  return 7;
}
/** @endcond */


/** @cond doxygen-libsbml-internal */
/*
 * @return the SBML object corresponding to next XMLToken in the
 * XMLInputStream or NULL if the token was not recognized.
 *
 * @param stream the XMLInputStream to use
 */
SBase*
ListOfLocalParameters::createObject (XMLInputStream& stream)
{
  const string& name   = stream.peek().getName();
  SBase*        object = 0;


  if (name == "localParameter")
  {
    try
    {
      object = new LocalParameter(getSBMLNamespaces());
    }
    catch (SBMLConstructorException*)
    {
      object = new LocalParameter(SBMLDocument::getDefaultLevel(),
        SBMLDocument::getDefaultVersion());
    }
    catch ( ... )
    {
      object = new LocalParameter(SBMLDocument::getDefaultLevel(),
        SBMLDocument::getDefaultVersion());
    }
    
    if (object) mItems.push_back(object);
  }

  return object;
}
/** @endcond */


/** @cond doxygen-c-only */


/**
 * Creates a new LocalParameter_t structure using the given SBML @p level
 * and @p version values.
 *
 * @param level an unsigned int, the SBML Level to assign to this
 * LocalParameter
 *
 * @param version an unsigned int, the SBML Version to assign to this
 * LocalParameter
 *
 * @return a pointer to the newly created LocalParameter_t structure.
 *
 * @note Once a LocalParameter has been added to an SBMLDocument, the @p
 * level and @p version for the document @em override those used to create
 * the LocalParameter.  Despite this, the ability to supply the values at
 * creation time is an important aid to creating valid SBML.  Knowledge of
 * the intended SBML Level and Version  determine whether it is valid to
 * assign a particular value to an attribute, or whether it is valid to add
 * an object to an existing SBMLDocument.
 */
LIBSBML_EXTERN
LocalParameter_t *
LocalParameter_create (unsigned int level, unsigned int version)
{
  try
  {
    LocalParameter* obj = new LocalParameter(level,version);
    return obj;
  }
  catch (SBMLConstructorException)
  {
    return NULL;
  }
}


/**
 * Creates a new LocalParameter_t structure using the given
 * SBMLNamespaces_t structure.
 *
 * @param sbmlns SBMLNamespaces, a pointer to an SBMLNamespaces structure
 * to assign to this LocalParameter
 *
 * @return a pointer to the newly created LocalParameter_t structure.
 *
 * @note Once a LocalParameter has been added to an SBMLDocument, the
 * @p sbmlns namespaces for the document @em override those used to create
 * the LocalParameter.  Despite this, the ability to supply the values at creation time
 * is an important aid to creating valid SBML.  Knowledge of the intended SBML
 * Level and Version determine whether it is valid to assign a particular value
 * to an attribute, or whether it is valid to add an object to an existing
 * SBMLDocument.
 */
LIBSBML_EXTERN
LocalParameter_t *
LocalParameter_createWithNS (SBMLNamespaces_t* sbmlns)
{
  try
  {
    LocalParameter* obj = new LocalParameter(sbmlns);
    return obj;
  }
  catch (SBMLConstructorException)
  {
    return NULL;
  }
}


/**
 * Frees the given LocalParameter_t structure.
 *
 * @param p the LocalParameter_t structure to be freed.
 */
LIBSBML_EXTERN
void
LocalParameter_free (LocalParameter_t *p)
{
  delete p;
}


/**
 * Creates a deep copy of the given LocalParameter_t structure
 * 
 * @param p the LocalParameter_t structure to be copied
 * 
 * @return a (deep) copy of the given LocalParameter_t structure.
 */
LIBSBML_EXTERN
LocalParameter_t *
LocalParameter_clone (const LocalParameter_t *p)
{
  return static_cast<LocalParameter_t*>( p->clone() );
}


/**
 * Returns a list of XMLNamespaces_t associated with this LocalParameter_t
 * structure.
 *
 * @param p the LocalParameter_t structure
 * 
 * @return pointer to the XMLNamespaces_t structure associated with 
 * this SBML object
 */
LIBSBML_EXTERN
const XMLNamespaces_t *
LocalParameter_getNamespaces(LocalParameter_t *p)
{
  return p->getNamespaces();
}

/**
 * Takes a LocalParameter_t structure and returns its identifier.
 *
 * @param p the LocalParameter_t structure whose identifier is sought
 * 
 * @return the identifier of this LocalParameter_t, as a pointer to a string.
 */
LIBSBML_EXTERN
const char *
LocalParameter_getId (const LocalParameter_t *p)
{
  return p->isSetId() ? p->getId().c_str() : NULL;
}


/**
 * Takes a LocalParameter_t structure and returns its name.
 *
 * @param p the LocalParameter_t whose name is sought.
 *
 * @return the name of this LocalParameter_t, as a pointer to a string.
 */
LIBSBML_EXTERN
const char *
LocalParameter_getName (const LocalParameter_t *p)
{
  return p->isSetName() ? p->getName().c_str() : NULL;
}


/**
 * Takes a LocalParameter_t structure and returns its value.
 *
 * @param p the LocalParameter_t whose value is sought.
 *
 * @return the value assigned to this LocalParameter_t structure, as a @c double.
 */
LIBSBML_EXTERN
double
LocalParameter_getValue (const LocalParameter_t *p)
{
  return p->getValue();
}


/**
 * Takes a LocalParameter_t structure and returns its units.
 *
 * @param p the LocalParameter_t whose units are sought.
 *
 * @return the units assigned to this LocalParameter_t structure, as a pointer
 * to a string.  
 */
LIBSBML_EXTERN
const char *
LocalParameter_getUnits (const LocalParameter_t *p)
{
  return p->isSetUnits() ? p->getUnits().c_str() : NULL;
}


/**
 * Predicate returning @c true or @c false depending on whether the given
 * LocalParameter_t structure's identifier has been set.
 *
 * @param p the LocalParameter_t structure to query
 * 
 * @return @c non-zero (true) if the "id" attribute of the given
 * LocalParameter_t structure has been set, zero (false) otherwise.
 */
LIBSBML_EXTERN
int
LocalParameter_isSetId (const LocalParameter_t *p)
{
  return static_cast<int>( p->isSetId() );
}


/**
 * Predicate returning @c true or @c false depending on whether the given
 * LocalParameter_t structure's name has been set.
 *
 * @param p the LocalParameter_t structure to query
 * 
 * @return @c non-zero (true) if the "name" attribute of the given
 * LocalParameter_t structure has been set, zero (false) otherwise.
 */
LIBSBML_EXTERN
int
LocalParameter_isSetName (const LocalParameter_t *p)
{
  return static_cast<int>( p->isSetName() );
}


/**
 * Predicate returning @c true or @c false depending on whether the given
 * LocalParameter_t structure's value has been set.
 * 
 * @param p the LocalParameter_t structure to query
 * 
 * @return @c non-zero (true) if the "value" attribute of the given
 * LocalParameter_t structure has been set, zero (false) otherwise.
 *
 * @note In SBML Level 1 Version 1, a LocalParameter value is required and
 * therefore <em>should always be set</em>.  In Level 1 Version 2 and
 * later, the value is optional, and as such, may or may not be set.
 */
LIBSBML_EXTERN
int
LocalParameter_isSetValue (const LocalParameter_t *p)
{
  return static_cast<int>( p->isSetValue() );
}


/**
 * Predicate returning @c true or @c false depending on whether the given
 * LocalParameter_t structure's units have been set.
 *
 * @param p the LocalParameter_t structure to query
 * 
 * @return @c non-zero (true) if the "units" attribute of the given
 * LocalParameter_t structure has been set, zero (false) otherwise.
 */
LIBSBML_EXTERN
int
LocalParameter_isSetUnits (const LocalParameter_t *p)
{
  return static_cast<int>( p->isSetUnits() );
}


/**
 * Assigns the identifier of a LocalParameter_t structure.
 *
 * This makes a copy of the string passed in the param @p sid.
 *
 * @param p the LocalParameter_t structure to set.
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
LocalParameter_setId (LocalParameter_t *p, const char *sid)
{
  return (sid == NULL) ? p->setId("") : p->setId(sid);
}


/**
 * Assign the name of a LocalParameter_t structure.
 *
 * This makes a copy of the string passed in as the argument @p name.
 *
 * @param p the LocalParameter_t structure to set.
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
LocalParameter_setName (LocalParameter_t *p, const char *name)
{
  return (name == NULL) ? p->unsetName() : p->setName(name);
}


/**
 * Assign the value of a LocalParameter_t structure.
 *
 * @param p the LocalParameter_t structure to set.
 * @param value the @c double value to use.
 *
 * @return integer value indicating success/failure of the
 * function.  @if clike The value is drawn from the
 * enumeration #OperationReturnValues_t. @endif The possible values
 * returned by this function are:
 * @li LIBSBML_OPERATION_SUCCESS
 */
LIBSBML_EXTERN
int
LocalParameter_setValue (LocalParameter_t *p, double value)
{
  return p->setValue(value);
}


/**
 * Assign the units of a LocalParameter_t structure.
 *
 * This makes a copy of the string passed in as the argument @p units.
 *
 * @param p the LocalParameter_t structure to set.
 * @param units the string to use as the identifier of the units to assign.
 *
 * @return integer value indicating success/failure of the
 * function.  @if clike The value is drawn from the
 * enumeration #OperationReturnValues_t. @endif The possible values
 * returned by this function are:
 * @li LIBSBML_OPERATION_SUCCESS
 * @li LIBSBML_INVALID_ATTRIBUTE_VALUE
 *
 * @note Using this function with units set to NULL is equivalent to
 * unsetting the "units" attribute.
 */
LIBSBML_EXTERN
int
LocalParameter_setUnits (LocalParameter_t *p, const char *units)
{
  return (units == NULL) ? p->unsetUnits() : p->setUnits(units);
}


/**
 * Unsets the name of this LocalParameter_t structure.
 * 
 * @param p the LocalParameter_t structure whose name is to be unset.
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
LocalParameter_unsetName (LocalParameter_t *p)
{
  return p->unsetName();
}


/**
 * Unsets the value of this LocalParameter_t structure.
 *
 * In SBML Level 1 Version 1, a parameter is required to have a value and
 * therefore this attribute <em>should always be set</em>.  In Level 1
 * Version 2 and beyond, a value is optional, and as such, may or may not be
 * set.
 *
 * @param p the LocalParameter_t structure whose value is to be unset.
 *
 * @return integer value indicating success/failure of the
 * function.  @if clike The value is drawn from the
 * enumeration #OperationReturnValues_t. @endif The possible values
 * returned by this function are:
 * @li LIBSBML_OPERATION_SUCCESS
 */
LIBSBML_EXTERN
int
LocalParameter_unsetValue (LocalParameter_t *p)
{
  return p->unsetValue();
}


/**
 * Unsets the units of this LocalParameter_t structure.
 * 
 * @param p the LocalParameter_t structure whose units are to be unset.
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
LocalParameter_unsetUnits (LocalParameter_t *p)
{
  return p->unsetUnits();
}


/**
 * Constructs and returns a UnitDefinition_t structure that expresses 
 * the units of this LocalParameter_t structure.
 *
 * @param p the LocalParameter_t structure whose units are to be returned.
 *
 * @return a UnitDefinition_t structure that expresses the units 
 * of this LocalParameter_t strucuture.
 *
 * @note This function returns the units of the LocalParameter_t expressed 
 * as a UnitDefinition_t. The units may be those explicitly declared. 
 * In the case where no units have been declared, NULL is returned.
 */
LIBSBML_EXTERN
UnitDefinition_t * 
LocalParameter_getDerivedUnitDefinition(LocalParameter_t *p)
{
  return p->getDerivedUnitDefinition();
}


/**
 * @return item in this ListOfLocalParameter with the given id or NULL if no such
 * item exists.
 */
LIBSBML_EXTERN
LocalParameter_t *
ListOfLocalParameters_getById (ListOf_t *lo, const char *sid)
{
  return (sid != NULL) ? 
    static_cast <ListOfLocalParameters *> (lo)->get(sid) : NULL;
}


/**
 * Removes item in this ListOf items with the given id or NULL if no such
 * item exists.  The caller owns the returned item and is responsible for
 * deleting it.
 */
LIBSBML_EXTERN
LocalParameter_t *
ListOfLocalParameters_removeById (ListOf_t *lo, const char *sid)
{
  return (sid != NULL) ? 
    static_cast <ListOfLocalParameters *> (lo)->remove(sid) : NULL;
}

/** @endcond */
LIBSBML_CPP_NAMESPACE_END
