/**
 * @file    XMLAttributes.cpp
 * @brief   XMLAttributes are a list of name/value pairs for XML elements
 * @author  Ben Bornstein
 *
 * $Id: XMLAttributes.cpp 11633 2010-08-03 03:53:27Z mhucka $
 * $HeadURL: https://sbml.svn.sourceforge.net/svnroot/sbml/trunk/libsbml/src/xml/XMLAttributes.cpp $
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

#include <cerrno>
#include <clocale>
#include <cstdlib>
#include <limits>
#include <sstream>

#include <sbml/xml/XMLErrorLog.h>
#include <sbml/xml/XMLAttributes.h>
/** @cond doxygen-libsbml-internal */
#include <sbml/xml/XMLOutputStream.h>
#include <sbml/util/util.h>
/** @endcond */

/** @cond doxygen-ignored */

using namespace std;

/** @endcond */

LIBSBML_CPP_NAMESPACE_BEGIN

/*
 * @return s with whitespace removed from the beginning and end.
 */
static const std::string
trim (const std::string& s)
{
  static const std::string whitespace(" \t\r\n");

  std::string::size_type begin = s.find_first_not_of(whitespace);
  std::string::size_type end   = s.find_last_not_of (whitespace);

  return (begin == std::string::npos) ? "" : s.substr(begin, end - begin + 1);
}


/*
 * Creates a new empty XMLAttributes set.
 */
XMLAttributes::XMLAttributes () : mLog( 0 )
{
}


/*
 * Destroys this XMLAttributes set.
 */
XMLAttributes::~XMLAttributes ()
{
}

/*
 * Copy constructor; creates a copy of this XMLAttributes set.
 */
XMLAttributes::XMLAttributes(const XMLAttributes& orig)
{
  this->mNames.assign( orig.mNames.begin(), orig.mNames.end() ); 
  this->mValues.assign( orig.mValues.begin(), orig.mValues.end() ); 
  this->mElementName = orig.mElementName;
  this->mLog = orig.mLog;
}


/*
 * Assignment operator for XMLAttributes.
 */
XMLAttributes& 
XMLAttributes::operator=(const XMLAttributes& orig)
{
  if(&orig!=this)
  {
    this->mNames.assign( orig.mNames.begin(), orig.mNames.end() ); 
    this->mValues.assign( orig.mValues.begin(), orig.mValues.end() ); 
    this->mElementName = orig.mElementName;
    this->mLog = orig.mLog;
  }

  return *this;
}

/*
 * Creates and returns a deep copy of this XMLAttributes set.
 * 
 * @return a (deep) copy of this XMLAttributes set.
 */
XMLAttributes* 
XMLAttributes::clone () const
{
  return new XMLAttributes(*this);
}



/*
 * Adds an attribute (a name/value pair) to this XMLAttributes set.  
 * If name with the same namespace URI already exists in this attribute set, 
 * its value will be replaced.
 */
int
XMLAttributes::add (const std::string& name,
		    const std::string& value,
		    const std::string& namespaceURI,
		    const std::string& prefix)
{
  int index = getIndex(name, namespaceURI);

  // since in the old version of the method the XMLTriple was initialized
  // with empty strings for the prefix and the uri, I assume that only
  // attributes that are not from the default namespace should have a set
  // prefix and uri.

  if (index == -1)
  {
    mNames .push_back( XMLTriple(name, namespaceURI, prefix) );
    mValues.push_back( value );
  }
  else
  {
    mValues[index] = value;
    mNames[index]  = XMLTriple(name, namespaceURI, prefix);
  }
  return LIBSBML_OPERATION_SUCCESS;
}


/*
 * Adds an attribute with the given triple/value pair to this XMLAttributes set.
 * If name with the same namespaceURI already exists in this attribute set,
 * its value will be replaced.
 */
int 
XMLAttributes::add ( const XMLTriple& triple, const std::string& value)
{
  return add(triple.getName(), value, triple.getURI(), triple.getPrefix());
}


/** @cond doxygen-libsbml-internal */
/*
 * Adds an attribute with the given name/value pair to this XMLAttributes set.  
 * This is really the add function but an attribute with same name wont 
 * be overwritten - this is for annotations
 */
int
XMLAttributes::addResource (const std::string& name, const std::string& value)
{
  mNames .push_back( XMLTriple(name, "", "") );
  mValues.push_back( value );
  return LIBSBML_OPERATION_SUCCESS;
}
/** @endcond */


/*
 * Removes an attribute with the given index from this XMLAttributes set.  
 * This is for annotations
 */
int
XMLAttributes::removeResource (int n)
{
  if (n < 0 || n >= getLength()) 
  {
    return LIBSBML_INDEX_EXCEEDS_SIZE;
  }

  vector<XMLTriple>::iterator   names_iter  = mNames.begin()  + n;
  vector<std::string>::iterator values_iter = mValues.begin() + n;

  mNames.erase(names_iter);
  mValues.erase(values_iter);

  return LIBSBML_OPERATION_SUCCESS;
}


/*
 * Removes an attribute with the given index from this XMLAttributes set.  
 * This is for annotations
 */
int
XMLAttributes::remove (int n)
{
  return removeResource(n);
}


/*
 * Removes an attribute with the given name and namespace URI from this 
 * XMLAttributes set.
 */
int 
XMLAttributes::remove (const std::string& name, const std::string& uri)
{
  return remove(getIndex(name,uri));
}


/*
 * Removes an attribute with the given triple from this XMLAttributes set.
 */
int 
XMLAttributes::remove (const XMLTriple& triple)
{
  return remove(getIndex(triple));
}


/*
 * Clears (deletes) all attributes in this XMLAttributes object.
 */
int 
XMLAttributes::clear()
{
  mNames.clear();
  mValues.clear();
  return LIBSBML_OPERATION_SUCCESS;
}


/*
 * Lookup the index of an attribute with the given name.
 *
 * @return the index of an attribute with the given name, or -1 if not present.
 */
int
XMLAttributes::getIndex (const std::string& name) const
{
  for (int index = 0; index < getLength(); ++index)
  {
    if (getName(index) == name) return index;
  }
  
  return -1;
}


/*
 * Lookup the index of an attribute with the given name and namespace URI
 *
 * @return the index of an attribute with the given name and namespace URI, 
 * or -1 if not present.
 */
int
XMLAttributes::getIndex (const std::string& name, const std::string& uri) const
{
  for (int index = 0; index < getLength(); ++index)
  {
    if ( (getName(index) == name) && (getURI(index) == uri) ) return index;
  }
  
  return -1;
}


/*
 * Lookup the index of an attribute by XMLTriple.
 *
 * @return the index of an attribute with the given XMLTriple, or -1 if not present.
 */
int 
XMLAttributes::getIndex (const XMLTriple& triple) const
{

  for (int index = 0; index < getLength(); ++index)
  {
    if (mNames[index] == triple) return index;
  }
  
  return -1;
}


/*
 * @return the number of attributes in this list.
 */
int
XMLAttributes::getLength () const
{
  return mNames.size();
}


/*
 * @return the name of an attribute in this list (by position).  If index
 * is out of range, an empty string will be returned.  Use hasAttribute(index)
 * to test for attribute existence.
 */
std::string
XMLAttributes::getName (int index) const
{
  return (index < 0 || index >= getLength()) ? "" : mNames[index].getName();
}


/*
 * @return the namespace prefix of an attribute in this list (by
 * position).  If index is out of range, an empty string will be
 * returned.  Use hasAttribute(index) to test for attribute existence.
 */
std::string
XMLAttributes::getPrefix (int index) const
{
  return (index < 0 || index >= getLength()) ? "" : mNames[index].getPrefix();
}


/*
 * @return the prefixed name of an attribute in this list (by
 * position).  If index is out of range, an empty string will be
 * returned.  Use hasAttribute(index) to test for 
 * attribute existence.
 */
std::string
XMLAttributes::getPrefixedName (int index) const
{
  return (index < 0 || index >= getLength()) ? "" : mNames[index].getPrefixedName();
}


/*
 * @return the namespace URI of an attribute in this list (by position).
 * If index is out of range, an empty string will be returned.  Use
 * hasAttribute(index) to test for attribute existence.
 */
std::string
XMLAttributes::getURI (int index) const
{
  return (index < 0 || index >= getLength()) ? "" : mNames[index].getURI();
}


/*
 * @return the value of an attribute in the list (by position).  If index
 * is out of range, an empty string will be returned.  Use hasAttribute(index)
 * to test for attribute existence.
 */
std::string
XMLAttributes::getValue (int index) const
{
  return (index < 0 || index >= getLength()) ? "" : mValues[index];
}


/*
 * Lookup an attribute's value by name.
 *
 * @return The attribute value as a string.  If an attribute with the
 * given name does not exist, an empty string will be returned.  Use
 * hasAttribute(name) to test for attribute existence.
 */
std::string
XMLAttributes::getValue (const std::string name) const
{
  return getValue( getIndex(name) );
}


/*
 * Lookup an attribute's value with the name and namespace URI.
 *
 * @return The attribute value as a string.  If an attribute with the
 * given name does not exist, an empty string will be returned.  Use
 * hasAttribute(name,uri) to test for attribute existence.
 */
std::string
XMLAttributes::getValue (const std::string name, const std::string uri) const
{
  return getValue( getIndex(name,uri) );
}


/*
 * Return an attribute's value by XMLTriple.
 *
 * @return The attribute value as a string.
 * If an attribute with the
 * given XMLTriple does not exist, an empty string will be returned.
 * Use hasAttribute(triple) to test for attribute existence.
 */
std::string 
XMLAttributes::getValue (const XMLTriple& triple) const
{
  return getValue( getIndex(triple) );
}


/*
 * Predicate returning @c true or @c false depending on whether
 * an attribute with the given name exists in this XMLAttributes.
 */
bool 
XMLAttributes::hasAttribute (int index) const 
{ 
   return ( (index >= 0) && (index < getLength()) );
}


/*
 * Predicate returning @c true or @c false depending on whether
 * an attribute with the given name and namespace URI exists in this XMLAttributes.
 *
 * @param name a string, the name of the attribute 
 * @param uri  a string, the namespace URI of the attribute.
 *
 * @return @c true if an attribute with the given name exists in this
 * XMLAttributes, @c false otherwise.
 *
 */
bool 
XMLAttributes::hasAttribute (const std::string name, const std::string uri) const 
{ 
  return ( getIndex(name,uri) != -1 ); 
}


/*
 * Predicate returning @c true or @c false depending on whether
 * an attribute with the given XML triple exists in this XMLAttributes.
 *
 * @param triple an XMLTriple, the XML triple of the attribute 
 *
 * @return @c true if an attribute with the given XML triple exists in this
 * XMLAttributes, @c false otherwise.
 *
 */
bool 
XMLAttributes::hasAttribute (const XMLTriple& triple) const 
{ 
  return ( getIndex(triple) != -1 ); 
}


/*
 * @return true if this XMLAttributes set is empty, false otherwise.
 */
bool
XMLAttributes::isEmpty () const
{
  return (getLength() == 0);
}


/** @cond doxygen-libsbml-internal */
/*
 * Reads the value for the attribute with the index into value.  If attribute 
 * was not found or value could not be interpreted as a boolean, value is not 
 * modified.
 *
 * According to the W3C XML Schema, valid boolean values are: "true",
 * "false", "1", and "0" (case-insensitive).  For more information, see:
 * http://www.w3.org/TR/xmlschema-2/#boolean
 *
 * If an XMLErrorLog is passed in datatype format errors are logged.  If
 * required is true, missing attributes are also logged.
 *
 * @returns true if the attribute was read into value, false otherwise.
 */
bool
XMLAttributes::readInto (  int          index
			 , const std::string& name
                         , bool&        value
                         , XMLErrorLog* log
                         , bool         required ) const
{
  bool assigned = false;
  bool missing  = true;

  if ( index != -1 )
  {
    const string& trimmed = trim( getValue(index) );
    if ( !trimmed.empty() )
    {
      missing = false;

      if (trimmed == "0" || trimmed == "false")
      {
        value    = false;
        assigned = true;
      }
      else if (trimmed == "1" || trimmed == "true")
      {
        value    = true;
        assigned = true;
      }
    }
  }

  if ( !log ) log = mLog;

  if ( log && !assigned )
  {
    if ( !missing ) attributeTypeError(name, Boolean, log);
    else if ( required ) attributeRequiredError (name, log);
  }

  return assigned;
}
/** @endcond */


/*
 * Reads the value for the attribute name into value.  If the given local
 * name was not found or value could not be interpreted as a boolean, 
 * value is not modified.
 *
 * According to the W3C XML Schema, valid boolean values are: "true",
 * "false", "1", and "0" (case-insensitive).  For more information, see:
 * http://www.w3.org/TR/xmlschema-2/#boolean
 *
 * If an XMLErrorLog is passed in datatype format errors are logged.  If
 * required is true, missing attributes are also logged.
 *
 * @returns true if the attribute was read into value, false otherwise.
 */
bool
XMLAttributes::readInto (  const std::string&   name
                         , bool&                value
                         , XMLErrorLog*         log
                         , bool                 required ) const
{
   return readInto(getIndex(name), name, value, log, required);
}


/*
 * Reads the value for the attribute XMLTriple into value.  If XMLTriple was not
 * found or value could not be interpreted as a boolean, value is not modified.
 *
 * According to the W3C XML Schema, valid boolean values are: "true",
 * "false", "1", and "0" (case-insensitive).  For more information, see:
 * http://www.w3.org/TR/xmlschema-2/#boolean
 *
 * If an XMLErrorLog is passed in datatype format errors are logged.  If
 * required is true, missing attributes are also logged.
 *
 * @returns true if the attribute was read into value, false otherwise.
 */
bool
XMLAttributes::readInto (  const XMLTriple& triple
                         , bool&            value
                         , XMLErrorLog*     log
                         , bool             required ) const
{
   return readInto(getIndex(triple), triple.getPrefixedName(), value, log, required);
}


/** @cond doxygen-libsbml-internal */
/*
 * Reads the value for the attribute with the given index into value.  
 * If name was not found or value could be interpreted as a double, value 
 * is not modified.
 *
 * According to the W3C XML Schema, valid doubles are the same as valid
 * doubles for C and the special values "INF", "-INF", and "NaN"
 * (case-sensitive).  For more information, see:
 * http://www.w3.org/TR/xmlschema-2/#double
 *
 * If an XMLErrorLog is passed in datatype format errors are logged.  If
 * required is true, missing attributes are also logged.
 *
 * @returns @c true if the attribute was read into value, @c false otherwise.
 *
 */
bool
XMLAttributes::readInto (  int          index
			 , const std::string& name
                         , double&      value
                         , XMLErrorLog* log
                         , bool         required ) const
{
  bool assigned = false;
  bool missing  = true;

  if ( index != -1 )
  {
    const std::string& trimmed = trim( getValue(index) );
    if ( !trimmed.empty() )
    {
      if (trimmed == "-INF")
      {
        value    = - numeric_limits<double>::infinity();
        assigned = true;
      }
      else if (trimmed == "INF")
      {
        value    = numeric_limits<double>::infinity();
        assigned = true;
      }
      else if (trimmed == "NaN")
      {
        value    = numeric_limits<double>::quiet_NaN();
        assigned = true;
      }
      else
      {
        // Ensure C locale
        char*  ptr    =  setlocale(LC_ALL, NULL);
        std::string locale = (ptr) ? ptr : "";
        setlocale(LC_ALL, "C");

        errno               = 0;
        char*        endptr = 0;
        const char*  nptr   = trimmed.c_str();
        double       result = strtod(nptr, &endptr);
        unsigned int length = endptr - nptr;

        // Restore previous locale
        setlocale(LC_ALL, locale.empty() ? 0 : locale.c_str());

        if ((length == trimmed.size()) && (errno != ERANGE))
        {
          value    = result;
          assigned = true;
        }
      }
    }
  }

  if ( !log ) log = mLog;

  if ( log && !assigned )
  {
    if ( !missing ) attributeTypeError(name, Double, log);
    else if ( required ) attributeRequiredError (name, log);
  }

  return assigned;
}
/** @endcond */


/*
 * Reads the value for the attribute with the given XMLTriple into value.  
 * If the triple was not found or value could be interpreted as a double, 
 *value is not modified.
 *
 * According to the W3C XML Schema, valid doubles are the same as valid
 * doubles for C and the special values "INF", "-INF", and "NaN"
 * (case-sensitive).  For more information, see:
 * http://www.w3.org/TR/xmlschema-2/#double
 *
 * If an XMLErrorLog is passed in datatype format errors are logged.  If
 * required is true, missing attributes are also logged.
 *
 * @returns @c true if the attribute was read into value, @c false otherwise.
 *
 */
bool
XMLAttributes::readInto (  const XMLTriple& triple
                         , double&          value
                         , XMLErrorLog*     log
                         , bool             required ) const
{
  return readInto(getIndex(triple), triple.getPrefixedName(), value, log, required);
}


/*
 * Reads the value for the attribute name into value.  If name was not
 * found or value could not be interpreted as a double, value is not
 * modified.
 *
 * According to the W3C XML Schema, valid doubles are the same as valid
 * doubles for C and the special values "INF", "-INF", and "NaN"
 * (case-sensitive).  For more information, see:
 * http://www.w3.org/TR/xmlschema-2/#double
 *
 * If an XMLErrorLog is passed in datatype format errors are logged.  If
 * required is true, missing attributes are also logged.
 *
 * @returns true if the attribute was read into value, false otherwise.
 */
bool
XMLAttributes::readInto (  const std::string&   name
                         , double&          value
                         , XMLErrorLog*     log
                         , bool             required ) const
{
  return readInto(getIndex(name), name, value, log, required);
}


/** @cond doxygen-libsbml-internal */
/*
 * Reads the value for the attribute with the given index into value.  
 * If the attribute was not found or value could be interpreted as a long, 
 * value is not modified.
 *
 * According to the W3C XML Schema valid integers include zero, *all*
 * positive and *all* negative whole numbers.  For practical purposes, we
 * limit values to what can be stored in a long.  For more information,
 * see: http://www.w3.org/TR/xmlschema-2/#integer
 *
 * If an XMLErrorLog is passed in datatype format errors are logged.  If
 * required is true, missing attributes are also logged.
 *
 * @returns @c true if the attribute was read into value, @c false otherwise.
 *
 */
bool
XMLAttributes::readInto (  int          index
			 , const std::string& name
                         , long&        value
                         , XMLErrorLog* log
                         , bool         required ) const
{
  bool assigned = false;
  bool missing  = true;

  if ( index != -1 )
  {
    const std::string& trimmed = trim( getValue(index) );
    if ( !trimmed.empty() )
    {
      missing = false;

      errno               = 0;
      char*        endptr = 0;
      const char*  nptr   = trimmed.c_str();
      long         result = strtol(nptr, &endptr, 10);
      unsigned int length = endptr - nptr;

      if ((length == trimmed.size()) && (errno != ERANGE))
      {
        value    = result;
        assigned = true;
      }
    }
  }

  if ( !log ) log = mLog;

  if ( log && !assigned )
  {
    if ( !missing ) attributeTypeError(name, Integer, log);
    else if ( required ) attributeRequiredError (name, log);
  }

  return assigned;
}
/** @endcond */


/*
 * Reads the value for the attribute XMLTriple into value.  
 * If the XMLTriple was not found or value could be interpreted as a long, 
 * value is not modified.
 *
 * According to the W3C XML Schema valid integers include zero, *all*
 * positive and *all* negative whole numbers.  For practical purposes, we
 * limit values to what can be stored in a long.  For more information,
 * see: http://www.w3.org/TR/xmlschema-2/#integer
 *
 * If an XMLErrorLog is passed in datatype format errors are logged.  If
 * required is true, missing attributes are also logged.
 *
 * @returns @c true if the attribute was read into value, @c false otherwise.
 *
 */
bool
XMLAttributes::readInto (  const XMLTriple& triple
                         , long&            value
                         , XMLErrorLog*     log
                         , bool             required ) const
{
  return readInto(getIndex(triple), triple.getPrefixedName(), value, log, required);
}

/*
 * Reads the value for the attribute name into value.  If name was not
 * found or value could not be interpreted as an long, value is not modified.
 *
 * According to the W3C XML Schema valid integers include zero, *all*
 * positive and *all* negative whole numbers.  For practical purposes, we
 * limit values to what can be stored in a long.  For more information,
 * see: http://www.w3.org/TR/xmlschema-2/#integer
 *
 * If an XMLErrorLog is passed in datatype format errors are logged.  If
 * required is true, missing attributes are also logged.
 *
 * @returns true if the attribute was read into value, false otherwise.
 */
bool
XMLAttributes::readInto (  const std::string& name
                         , long&              value
                         , XMLErrorLog*       log
                         , bool               required ) const
{
  return readInto(getIndex(name), name, value, log, required);
}


/** @cond doxygen-libsbml-internal */
/*
 * Reads the value for the attribute with the given index into value.  
 * If the attribute was not found or value could be interpreted as an integer, 
 * value is not modified.
 *
 * According to the W3C XML Schema valid integers include zero, *all*
 * positive and *all* negative whole numbers.  For practical purposes, we
 * limit values to what can be stored in a int.  For more information,
 * see: http://www.w3.org/TR/xmlschema-2/#integer
 *
 * If an XMLErrorLog is passed in datatype format errors are logged.  If
 * required is true, missing attributes are also logged.
 *
 * @returns @c true if the attribute was read into value, @c false otherwise.
 *
 */
bool
XMLAttributes::readInto (  int          index
			 , const std::string& name
                         , int&         value
                         , XMLErrorLog* log
                         , bool         required ) const
{
  long  temp;
  bool  assigned = readInto(index, name, temp, log, required);

  if (assigned) value = temp;
  return assigned;
}
/** @endcond */


/*
 * Reads the value for the attribute with the given XMLTriple into value.  
 * If the XMLTriple was not found or value could be interpreted as an int, 
 * value is not modified.
 *
 * According to the W3C XML Schema valid integers include zero, *all*
 * positive and *all* negative whole numbers.  For practical purposes, we
 * limit values to what can be stored in a int.  For more information,
 * see: http://www.w3.org/TR/xmlschema-2/#integer
 *
 * If an XMLErrorLog is passed in datatype format errors are logged.  If
 * required is true, missing attributes are also logged.
 *
 * @returns @c true if the attribute was read into value, @c false otherwise.
 *
 */
bool
XMLAttributes::readInto (  const XMLTriple& triple
                         , int&             value
                         , XMLErrorLog*     log
                         , bool             required ) const
{
   return readInto(getIndex(triple), triple.getPrefixedName(), value, log, required);    
}


/*
 * Reads the value for the attribute name into value.  If name was not
 * found or value could not be interpreted as an int, value is not modified.
 *
 * According to the W3C XML Schema valid integers include zero, *all*
 * positive and *all* negative whole numbers.  For practical purposes, we
 * limit values to what can be stored in a int.  For more information,
 * see: http://www.w3.org/TR/xmlschema-2/#integer
 *
 * If an XMLErrorLog is passed in datatype format errors are logged.  If
 * required is true, missing attributes are also logged.
 *
 * @returns true if the attribute was read into value, false otherwise.
 */
bool
XMLAttributes::readInto (  const std::string&  name
                         , int&                value
                         , XMLErrorLog*        log
                         , bool                required ) const
{
  return readInto(getIndex(name), name, value, log, required);
}


/** @cond doxygen-libsbml-internal */
/*
 * Reads the value for the attribute with the given index into value.  
 * If the attribute was not found or value could be interpreted as an 
 * unsigned int, value is not modified.
 *
 * According to the W3C XML Schema valid integers include zero, *all*
 * positive and *all* negative whole numbers.  For practical purposes, we
 * limit values to what can be stored in a unsigned int.  For more
 * information, see: http://www.w3.org/TR/xmlschema-2/#integer
 *
 * If an XMLErrorLog is passed in datatype format errors are logged.  If
 * required is true, missing attributes are also logged.
 *
 * @returns @c true if the attribute was read into value, @c false otherwise.
 *
 */
bool
XMLAttributes::readInto (  int           index
			 , const std::string& name
                         , unsigned int& value
                         , XMLErrorLog*  log
                         , bool          required ) const
{
  long  temp;
  bool  assigned = readInto(index, name, temp, log, required);

  if (assigned && temp >= 0) value = temp;
  else assigned = false;

  return assigned;
}
/** @endcond */


/*
 * Reads the value for the attribute with the given XMLTriple into value.  
 * If the XMLTriple was not found or value could be interpreted as an unsigned int, 
 * value is not modified.
 *
 * According to the W3C XML Schema valid integers include zero, *all*
 * positive and *all* negative whole numbers.  For practical purposes, we
 * limit values to what can be stored in a unsigned int.  For more
 * information, see: http://www.w3.org/TR/xmlschema-2/#integer
 *
 * If an XMLErrorLog is passed in datatype format errors are logged.  If
 * required is true, missing attributes are also logged.
 *
 * @returns @c true if the attribute was read into value, @c false otherwise.
 *
 */
bool
XMLAttributes::readInto (  const XMLTriple& triple
                         , unsigned int&    value
                         , XMLErrorLog*     log
                         , bool             required ) const
{
  return readInto(getIndex(triple), triple.getPrefixedName(), value, log, required);
}


/*
 * Reads the value for the attribute name into value.  If name was not
 * found or value could be interpreted as an unsigned int, value is not
 * modified.
 *
 * According to the W3C XML Schema valid integers include zero, *all*
 * positive and *all* negative whole numbers.  For practical purposes, we
 * limit values to what can be stored in a unsigned int.  For more
 * information, see: http://www.w3.org/TR/xmlschema-2/#integer
 *
 * If an XMLErrorLog is passed in datatype format errors are logged.  If
 * required is true, missing attributes are also logged.
 *
 * @returns true if the attribute was read into value, false otherwise.
 */
bool
XMLAttributes::readInto (  const std::string&  name
                         , unsigned int&       value
                         , XMLErrorLog*        log
                         , bool                required ) const
{
  return readInto(getIndex(name), name, value, log, required);
}


/** @cond doxygen-libsbml-internal */
/*
 * Reads the value for the attribute with the given index into value.  
 * If the attribute was not found, value is not modified.
 *
 * If an XMLErrorLog is passed in and required is true, missing
 * attributes are logged.
 *
 * @returns @c true if the attribute was read into value, @c false otherwise.
 *
 */
bool
XMLAttributes::readInto (  int          index
			 , const std::string& name
			 , std::string& value
                         , XMLErrorLog* log
                         , bool         required ) const
{
  bool assigned = false;

  if ( index != -1 )
  {
    value    = getValue(index);
    assigned = true;
  }

  if ( !log ) log = mLog;

  if ( log && !assigned && required )
  {
    attributeRequiredError(name, log);
  }

  return assigned;
}
/** @endcond */


/*
 * Reads the value for the attribute with the given XMLTriple into value.  
 * If the XMLTriple was not found, value is not modified.
 *
 * If an XMLErrorLog is passed in and required is true, missing
 * attributes are logged.
 *
 * @returns @c true if the attribute was read into value, @c false otherwise.
 *
 */
bool
XMLAttributes::readInto (  const XMLTriple& triple
			 , std::string&     value
                         , XMLErrorLog*     log
                         , bool             required ) const
{
  return readInto(getIndex(triple), triple.getPrefixedName(), value, log, required);
}


/*
 * Reads the value for the attribute name into value.  If name was not
 * found, value is not modified.
 *
 * If an XMLErrorLog is passed in and required is true, missing
 * attributes are logged.
 *
 * @returns true if the attribute was read into value, false otherwise.
 */
bool
XMLAttributes::readInto (  const std::string& name
			 , std::string&       value
                         , XMLErrorLog*       log
                         , bool               required ) const
{
  return readInto(getIndex(name), name, value, log, required);
}


/** @cond doxygen-libsbml-internal */
/*
 * Writes this XMLAttributes set to stream.
 */
void
XMLAttributes::write (XMLOutputStream& stream) const
{
  for (int n = 0; n < getLength(); ++n)
  {
    if ( getPrefix(n).empty() )
    {
      stream.writeAttribute( getName(n), getValue(n) );
    }
    else
    {
      stream.writeAttribute( mNames[n], getValue(n) );
    }
  }
}
/** @endcond */


/** @cond doxygen-libsbml-internal */

/*
 * Logs an attribute format error.
 *
 * @param  name  Name of the attribute
 * @param  type  The datatype of the attribute value.
 */
void
XMLAttributes::attributeTypeError (  const std::string& name
				   , DataType           type
				   , XMLErrorLog*       log ) const
{
  ostringstream message;

  if ( !log ) log = mLog;
  if ( !log ) return;

  message << "The ";
  if ( !mElementName.empty() ) message << mElementName << ' ';
  message << name;

  switch ( type )
  {
    case XMLAttributes::Boolean:
      message <<
        " attribute must have a value of either \"true\" or \"false\""
        " (all lowercase).  The numbers \"1\" (true) and \"0\" (false) are"
        " also allowed, but not preferred.  For more information, see:"
        " http://www.w3.org/TR/xmlschema-2/#boolean.";
      break;

    case XMLAttributes::Double:
      message <<
        " attribute must be a double (decimal number).  To represent"
        " infinity use \"INF\", negative infinity use \"-INF\", and"
        " not-a-number use \"NaN\".  For more information, see:"
        " http://www.w3.org/TR/xmlschema-2/#double.";
      break;

    case XMLAttributes::Integer:
      message <<
        " attribute must be an integer (whole number).  For more"
        " information, see: http://www.w3.org/TR/xmlschema-2/#integer.";
      break;
  }

  log->add( XMLError(XMLAttributeTypeMismatch, message.str()) );
}


/*
 * Logs an error indicating a required attribute was missing.
 *
 * @param  name  Name of the attribute
 */
void
XMLAttributes::attributeRequiredError (const std::string&  name,
				       XMLErrorLog*        log) const
{
  ostringstream message;

  if ( !log ) log = mLog;
  if ( !log ) return;

  message << "The ";
  if ( !mElementName.empty() ) message << mElementName << ' ';
  message << "attribute '" << name << "' is required.";

  log->add( XMLError(MissingXMLRequiredAttribute, message.str()) );
}

/** @endcond */


/** @cond doxygen-libsbml-internal */
/*
 * Sets the XMLErrorLog this parser will use to log errors.
 */
int
XMLAttributes::setErrorLog (XMLErrorLog* log)
{
  if (mLog == log)
  {
    return LIBSBML_OPERATION_SUCCESS;
  }
  else if (log == NULL)
  {
    delete mLog;
    mLog = 0;
    return LIBSBML_OPERATION_SUCCESS;
  }
  else
  {
    mLog = log;
    return LIBSBML_OPERATION_SUCCESS;
  }
}
/** @endcond */


/** @cond doxygen-libsbml-internal */
/*
 * Inserts this XMLAttributes set into stream.
 */
LIBLAX_EXTERN
XMLOutputStream&
operator<< (XMLOutputStream& stream, const XMLAttributes& attributes)
{
  attributes.write(stream);
  return stream;
}
/** @endcond */


/** @cond doxygen-c-only */


/**
 * Creates a new empty XMLAttributes_t set.
 */
LIBLAX_EXTERN
XMLAttributes_t *
XMLAttributes_create (void)
{
  return new(nothrow) XMLAttributes;
}


/**
 * Frees the given XMLAttributes_t structure.
 *
 * @param xa the XMLAttributes_t structure to be freed.
 **/
LIBLAX_EXTERN
void
XMLAttributes_free (XMLAttributes_t *xa)
{
  delete static_cast<XMLAttributes*>(xa);
}


/**
 * Creates a deep copy of the given XMLAttributes_t structure.
 * 
 * @param att the XMLAttributes_t structure to be copied
 * 
 * @return a (deep) copy of the given XMLAttributes_t structure.
 */
LIBLAX_EXTERN
XMLAttributes_t *
XMLAttributes_clone (const XMLAttributes_t* att)
{
  return static_cast<XMLAttributes*>( att->clone() );
}


/**
 * Adds a name/value pair to this XMLAttributes_t structure.
 *
 * @param xa the XMLAttributes_t structure 
 * @param name a string, the local name of the attribute.
 * @param value a string, the value of the attribute.
 *
 * @return integer value indicating success/failure of the
 * function.  @if clike The value is drawn from the
 * enumeration #OperationReturnValues_t. @endif The possible values
 * returned by this function are:
 * @li LIBSBML_OPERATION_SUCCESS
 *
 * @note if local name already exists in this attribute set, its value 
 * will be replaced.
 **/
LIBLAX_EXTERN
int
XMLAttributes_add (XMLAttributes_t *xa,
                   const char *name,
                   const char *value)
{
  return xa->add(name, value);
}


/**
 * Adds a name/value pair to this XMLAttributes_t structure with a
 * prefix and URI defining a namespace.
 *
 * @param xa the XMLAttributes_t structure.
 * @param name a string, the local name of the attribute.
 * @param value a string, the value of the attribute.
 * @param namespaceURI a string, the namespace URI of the attribute.
 * @param prefix a string, the prefix of the namespace
 *
 * @return integer value indicating success/failure of the
 * function.  @if clike The value is drawn from the
 * enumeration #OperationReturnValues_t. @endif The possible values
 * returned by this function are:
 * @li LIBSBML_OPERATION_SUCCESS
 *
 * @note if local name with the same namespace URI already exists in this 
 * attribute set, its value will be replaced.
 **/
LIBLAX_EXTERN
int
XMLAttributes_addWithNamespace (XMLAttributes_t *xa,
                                const char *name,
                                const char *value,
                                const char* uri,
                                const char* prefix)
{
  return xa->add(name, value, uri, prefix);
}


/**
  * Adds an attribute with the given XMLtriple/value pair to this XMLAttributes_t structure.  
  *
  * @param xa the XMLAttributes_t structure.
  * @param triple an XMLTriple_t, the triple of the attribute.
  * @param value a string, the value of the attribute.
 *
 * @return integer value indicating success/failure of the
 * function.  @if clike The value is drawn from the
 * enumeration #OperationReturnValues_t. @endif The possible values
 * returned by this function are:
 * @li LIBSBML_OPERATION_SUCCESS
 **/
LIBLAX_EXTERN
int
XMLAttributes_addWithTriple ( XMLAttributes_t *xa 
			              , const XMLTriple_t* triple
			              , const char *value)
{
  return xa->add(*triple, value);
}


/**
  * Adds an attribute with the given name/value pair to this XMLAttributes_t structure.  
  *
  * This method is similar to the add method but an attribute with same name wont 
  * be overwritten. This facilitates the addition of multiple resource attributes 
  * to a annotations.
  *
  * @param xa the XMLAttributes_t structure.
  * @param name a string, the name of the attribute.
  * @param value a string, the value of the attribute.
LIBLAX_EXTERN
void
XMLAttributes_addResource (XMLAttributes_t *xa, 
			   const char *name, 
			   const char *value)
{
  xa->addResource(name, value);
}
 **/


/**
 * Removes an attribute (a name/value pair) from this XMLAttributes set.  
 *
 * @param xa the XMLAttributes_t structure.
 * @param n an integer the index of the resource to be deleted
 *
 * @return integer value indicating success/failure of the
 * function.  @if clike The value is drawn from the
 * enumeration #OperationReturnValues_t. @endif The possible values
 * returned by this function are:
 * @li LIBSBML_OPERATION_SUCCESS
 * @li LIBSBML_INDEX_EXCEEDS_SIZE
 */
LIBLAX_EXTERN
int
XMLAttributes_removeResource (XMLAttributes_t *xa, int n)
{
  return xa->removeResource(n);
}


/**
 * Removes an attribute (a name/value pair) from this XMLAttributes set.  
 *
 * @param xa the XMLAttributes_t structure.
 * @param n an integer the index of the resource to be deleted
 *
 * @return integer value indicating success/failure of the
 * function.  @if clike The value is drawn from the
 * enumeration #OperationReturnValues_t. @endif The possible values
 * returned by this function are:
 * @li LIBSBML_OPERATION_SUCCESS
 * @li LIBSBML_INDEX_EXCEEDS_SIZE
 */
LIBLAX_EXTERN
int
XMLAttributes_remove (XMLAttributes_t *xa, int n)
{
  return xa->remove(n);
}


/**
 * Removes an attribute with the given local name from this XMLAttributes set.  
 *
 * @param xa the XMLAttributes_t structure.
 * @param name a string, the local name of the attribute.
 *
 * @return integer value indicating success/failure of the
 * function.  @if clike The value is drawn from the
 * enumeration #OperationReturnValues_t. @endif The possible values
 * returned by this function are:
 * @li LIBSBML_OPERATION_SUCCESS
 * @li LIBSBML_INDEX_EXCEEDS_SIZE
 *
 * @note A prefix and namespace URI bound to the local name are set to empty 
 * in this function.
 * XMLAttributes_removeByNS(name,uri) or XMLAttributes_removeByTriple(triple) 
 * should be used to remove an attribute with the given local name and namespace.
 */
LIBLAX_EXTERN
int
XMLAttributes_removeByName (XMLAttributes_t *xa, const char* name)
{
  return xa->remove(name);
}


/**
 * Removes an attribute with the given name and namespace URI from this
 * XMLAttributes set.
 *
 * @param xa the XMLAttributes_t structure.
 * @param name a string, the local name of the attribute for which the index is required.
 * @param uri a string, the namespace URI of the attribute.
 *
 * @return integer value indicating success/failure of the
 * function.  @if clike The value is drawn from the
 * enumeration #OperationReturnValues_t. @endif The possible values
 * returned by this function are:
 * @li LIBSBML_OPERATION_SUCCESS
 * @li LIBSBML_INDEX_EXCEEDS_SIZE
 */
LIBLAX_EXTERN
int
XMLAttributes_removeByNS (XMLAttributes_t *xa, const char* name, const char* uri)
{
  return xa->remove(name, uri);
}


/**
 * Removes an attribute with the given triple from this XMLAttributes set.
 *
 * @param xa the XMLAttributes_t structure.
 * @param triple an XMLTriple, the XML triple of the attribute for which
 *        the index is required.
 *
 * @return integer value indicating success/failure of the
 * function.  @if clike The value is drawn from the
 * enumeration #OperationReturnValues_t. @endif The possible values
 * returned by this function are:
 * @li LIBSBML_OPERATION_SUCCESS
 * @li LIBSBML_INDEX_EXCEEDS_SIZE
 */
LIBLAX_EXTERN
int
XMLAttributes_removeByTriple (XMLAttributes_t *xa, const XMLTriple_t* triple)
{
  return xa->remove(*triple);
}


/**
 * Clears (deletes) all attributes in this XMLAttributes object.
 *
 * @param xa the XMLAttributes_t structure.
 *
 * @return integer value indicating success/failure of the
 * function.  @if clike The value is drawn from the
 * enumeration #OperationReturnValues_t. @endif The possible values
 * returned by this function are:
 * @li LIBSBML_OPERATION_SUCCESS
 */
LIBLAX_EXTERN
int 
XMLAttributes_clear(XMLAttributes_t *xa)
{
  return xa->clear();
}


/**
 * Return the index of an attribute with the given name.
 *
 * @param xa the XMLAttributes_t structure.
 * @param name a string, the local name of the attribute for which the index is required.
 *
 * @return the index of an attribute with the given local name, or -1 if not present.
 */
LIBLAX_EXTERN
int
XMLAttributes_getIndex (const XMLAttributes_t *xa, const char *name)
{
  return xa->getIndex(name);
}


/**
 * Return the index of an attribute with the given name and namespace URI.
 *
 * @param xa the XMLAttributes_t structure.
 * @param name a string, the local name of the attribute for which the index is required.
 * @param uri a string, the namespace URI of the attribute.
 *
 * @return the index of an attribute with the given local name and namespace URI, 
 * or -1 if not present.
 */
LIBLAX_EXTERN
int
XMLAttributes_getIndexByNS (const XMLAttributes_t *xa, const char *name, const char *uri)
{
  return xa->getIndex(name,uri);
}


/**
 * Return the index of an attribute with the given XML triple.
 *
 * @param xa the XMLAttributes_t structure.
 * @param triple an XMLTriple, the XML triple of the attribute for which
 *        the index is required.
 *
 * @return the index of an attribute with the given XMLTriple, or -1 if not present.
 */
LIBLAX_EXTERN
int 
XMLAttributes_getIndexByTriple (const XMLAttributes_t *xa, const XMLTriple_t* triple)
{
  return xa->getIndex(*triple);
}


/**
 * Return the number of attributes in the set.
 *
 * @param xa the XMLAttributes_t structure.
 *
 * @return the number of attributes in this XMLAttributes_t structure.
 **/
LIBLAX_EXTERN
int
XMLAttributes_getLength (const XMLAttributes_t *xa)
{
  return xa->getLength();
}


/**
 * Return the local name of an attribute in this XMLAttributes_t structure (by position).
 *
 * @param xa the XMLAttributes_t structure.
 * @param index an integer, the position of the attribute whose name is 
 * required.
 *
 * @return the local name of an attribute in this list (by position).  
 *         NULL will be returned if the name is empty.
 *
 * @note If index
 * is out of range, an empty string will be returned.  
 * Use XMLNamespaces_hasAttribute(...) > 0 to test for attribute existence.
 * to test for attribute existence.
 * Returned const char* should be freed with safe_free() by the caller.
 **/
LIBLAX_EXTERN
const char *
XMLAttributes_getName (const XMLAttributes_t *xa, int index)
{
  /**
   * I did this because MSVC and gcc handle .c_str() in different ways
   * meaning that with MSVC the actual string goes out of scope before 
   * the char * is returned and thus the char * is garbage once returned
   */
  if (xa->getName(index).empty())
    return NULL;
  else
    return safe_strdup(xa->getName(index).c_str());
}


/**
 * Return the value of an attribute in this XMLAttributes_t structure (by position).
 *
 * @param xa the XMLAttributes_t structure.
 * @param index an integer, the position of the attribute whose value is 
 * required.
 *
 * @return the value of an attribute in the list (by position).  
 *         NULL will be returned if the prefix is empty.
 *
 * @note If index
 * is out of range, an empty string will be returned.  
 * Use XMLNamespaces_hasAttribute(...) > 0 to test for attribute existence.
 * Returned const char* should be freed with safe_free() by the caller.
 **/
LIBLAX_EXTERN
const char *
XMLAttributes_getPrefix (const XMLAttributes_t *xa, int index)
{
  /**
   * I did this because MSVC and gcc handle .c_str() in different ways
   * meaning that with MSVC the actual string goes out of scope before 
   * the char * is returned and thus the char * is garbage once returned
   */
  if (xa->getPrefix(index).empty())
    return NULL;
  else
    return safe_strdup(xa->getPrefix(index).c_str());
}


/**
 * Return the namespace URI of an attribute in this XMLAttributes_t structure (by position).
 *
 * @param xa the XMLAttributes_t structure.
 * @param index an integer, the position of the attribute whose namespace URI is 
 * required.
 *
 * @return the namespace URI of an attribute in this list (by position).
 *         NULL will be returned if the URI is empty.
 *
 * @note If index is out of range, an empty string will be returned.  Use
 * XMLNamespaces_hasAttribute(...) > 0 to test for attribute existence.
 * Returned const char* should be freed with safe_free() by the caller.
 **/
LIBLAX_EXTERN
const char *
XMLAttributes_getURI (const XMLAttributes_t *xa, int index)
{
  /**
   * I did this because MSVC and gcc handle .c_str() in different ways
   * meaning that with MSVC the actual string goes out of scope before 
   * the char * is returned and thus the char * is garbage once returned
   */
  if (xa->getURI(index).empty())
    return NULL;
  else
    return safe_strdup(xa->getURI(index).c_str());
}


/**
 * Return the value of an attribute in this XMLAttributes_t structure (by position).
 *
 * @param xa the XMLAttributes_t structure.
 * @param index an integer, the position of the attribute whose value is 
 * required.
 *
 * @return the value of an attribute in the list (by position).  
 *         NULL will be returned if the value is empty.
 *
 * @note If index
 * is out of range, NULL will be returned.  
 * Use XMLAttributes_hasAttribute(...) > 0 to test for attribute existence.
 * Returned const char* should be freed with safe_free() by the caller.
 **/
LIBLAX_EXTERN
const char *
XMLAttributes_getValue (const XMLAttributes_t *xa, int index)
{
  /**
   * I did this because MSVC and gcc handle .c_str() in different ways
   * meaning that with MSVC the actual string goes out of scope before 
   * the char * is returned and thus the char * is garbage once returned
   */
  if (xa->getValue(index).empty())
    return NULL;
  else
    return safe_strdup(xa->getValue(index).c_str());
}


/**
 * Return an attribute's value by name.
 *
 * @param xa the XMLAttributes_t structure.
 * @param name a string, the local name of the attribute whose value is required.
 *
 * @return The attribute value as a string.  
 *         NULL will be returned if the value is empty.
 *
 * @note If an attribute with the
 * given local name does not exist, NULL will be returned.  Use
 * XMLAttributes_hasAttributeWithName(...) > 0 to test for attribute existence.
 * A namespace bound to the local name is not checked by this function.
 * Thus, if there are multiple attributes with the given local name and 
 * different namespaces, the value of an attribute with the smallest index 
 * among those attributes will be returned.
 * XMLAttributes_getValueByNS(...) or XMLAttributes_getValueByTriple(...) 
 * should be used to get a value of an attribute with the given local name 
 * and namespace.
 * Returned const char* should be freed with safe_free() by the caller.
 **/
LIBLAX_EXTERN
const char *
XMLAttributes_getValueByName (const XMLAttributes_t *xa, const char *name)
{
  /**
   * I did this because MSVC and gcc handle .c_str() in different ways
   * meaning that with MSVC the actual string goes out of scope before 
   * the char * is returned and thus the char * is garbage once returned
   */
  if (xa->getValue(name).empty())
    return NULL;
  else
    return safe_strdup(xa->getValue(name).c_str());
}


/**
 * Return a value of an attribute with the given local name and namespace URI.
 *
 * @param xa the XMLAttributes_t structure.
 * @param name a string, the local name of the attribute whose value is required.
 * @param uri  a string, the namespace URI of the attribute.
 *
 * @return The attribute value as a string.  
 * NULL will be returned if the value is empty.
 *
 * @note If an attribute with the 
 * given local name and namespace URI does not exist, an empty string will be 
 * returned.  
 * Use XMLAttributes_hasAttributeWithNS(...) to test for attribute existence.
 * Returned const char* should be freed with safe_free() by the caller.
 */
LIBLAX_EXTERN
const char *
XMLAttributes_getValueByNS (const XMLAttributes_t *xa, const char *name, const char* uri)
{
  /**
   * I did this because MSVC and gcc handle .c_str() in different ways
   * meaning that with MSVC the actual string goes out of scope before 
   * the char * is returned and thus the char * is garbage once returned
   */
  if (xa->getValue(name, uri).empty())
    return NULL;
  else
    return safe_strdup(xa->getValue(name, uri).c_str());
}


/**
 * Return an attribute's value by XMLTriple.
 *
 * @param xa the XMLAttributes_t structure.
 * @param triple an XMLTriple_t, the XML triple of the attribute 
 * whose value is required.
 *
 * @return The attribute value as a string.
 *         NULL will be returned if the value is empty.
 *
 * @note If an attribute with the
 * given XMLTriple does not exist, NULL will be returned.
 * Use XMLAttributes_hasAttributeWithTriple(..) > 0 to test for attribute existence.
 * Returned const char* should be freed with safe_free() by the caller.
 */
LIBLAX_EXTERN
const char *
XMLAttributes_getValueByTriple (const XMLAttributes_t *xa, const XMLTriple_t* triple)
{
  std::string val = xa->getValue(*triple);
  if (val.empty()) return NULL;

  return safe_strdup(val.c_str());
}


/**
 * Predicate returning @c true or @c false depending on whether
 * an attribute with the given index exists in this XMLAttributes_t
 * structure.
 *
 * @param xa the XMLAttributes_t structure.
 * @param index an integer, the position of the attribute.
 *
 * @return @c non-zero (true) if an attribute with the given index exists 
 * in this XMLAttributes_t structure, @c zero (false) otherwise.
 */
LIBLAX_EXTERN
int 
XMLAttributes_hasAttribute (const XMLAttributes_t *xa, int index)
{
  return static_cast<int>( xa->hasAttribute(index) );
}


/**
 * Predicate returning @c true or @c false depending on whether
 * an attribute with the given local name exists in this XMLAttributes_t
 * structure.
 *
 * @param xa the XMLAttributes_t structure.
 * @param name a string, the local name of the attribute.
 *
 * @return @c non-zero (true) if an attribute with the given local name 
 * exists in this XMLAttributes_t structure, @c zero (false) otherwise.
 */
LIBLAX_EXTERN
int 
XMLAttributes_hasAttributeWithName (const XMLAttributes_t *xa, const char* name)
{
  return static_cast<int>( xa->hasAttribute(name) );
}


/**
 * Predicate returning @c true or @c false depending on whether
 * an attribute with the given local name and namespace URI exists in this 
 * XMLAttributes_t structure.
 *
 * @param xa the XMLAttributes_t structure.
 * @param name a string, the local name of the attribute.
 * @param uri  a string, the namespace URI of the attribute.
 *
 * @return @c non-zero (true) if an attribute with the given local name 
 * and namespace URI exists in this XMLAttributes_t structure, @c zero (false) 
 * otherwise.
 */
LIBLAX_EXTERN
int 
XMLAttributes_hasAttributeWithNS (const XMLAttributes_t *xa, 
                                  const char* name, const char* uri)
{
  return static_cast<int>( xa->hasAttribute(name, uri) );
}


/**
 * Predicate returning @c true or @c false depending on whether
 * an attribute with the given XMLtriple_t exists in this XMLAttributes_t
 * structure.
 *
 * @param xa the XMLAttributes_t structure.
 * @param triple an XMLTriple_t, the XML triple of the attribute.
 *
 * @return @c non-zero (true) if an attribute with the given XMLTriple_t
 * exists in this XMLAttributes_t structure, @c zero (false) otherwise.
 */
LIBLAX_EXTERN
int 
XMLAttributes_hasAttributeWithTriple (const XMLAttributes_t *xa, const XMLTriple_t* triple)
{
  return static_cast<int>( xa->hasAttribute(*triple) );
}  


/**
 * Predicate returning @c true or @c false depending on whether 
 * this XMLAttributes_t structure is empty.
 *
 * @param xa the XMLAttributes_t structure.
 * 
 * @return @c non-zero (true) if this XMLAttributes_t structure is empty, 
 * @c zero (false) otherwise.
 **/
LIBLAX_EXTERN
int
XMLAttributes_isEmpty (const XMLAttributes_t *xa)
{
  return static_cast<int>( xa->isEmpty() );
}


/**
 * Reads the value for the attribute name into value.  If the given local
 * name was not found or value could be interpreted as a boolean, value is 
 * not modified.
 *
 * According to the W3C XML Schema, valid boolean values are: "true",
 * "false", "1", and "0" (case-insensitive).  For more information, see:
 * http://www.w3.org/TR/xmlschema-2/#boolean
 *
 * If an XMLErrorLog is passed in datatype format errors are logged.  If
 * required is true, missing attributes are also logged.
 *
 * @param xa the XMLAttributes_t structure.
 * @param name a string, the local name of the attribute.
 * @param value a boolean, the value of the attribute.
 * @param log an XMLErrorLog_t, the error log.
 * @param required a boolean, indicating whether the attribute is required.
 *
 * @returns @c non-zero (true) if the attribute was read into value, 
 * @c zero (false) otherwise.
 *
 * @note A namespace bound to the given local name is not checked by this 
 * function. readIntoBooleanByTriple(...) should be used to read a value for 
 * an attribute name with a prefix and namespace.
 *
 **/
LIBLAX_EXTERN
int
XMLAttributes_readIntoBoolean (XMLAttributes_t *xa,
			       const char *name,
			       int *value,
			       XMLErrorLog_t *log,
			       int required)
{
  bool temp;
  bool result = xa->readInto(name, temp, log, required);
  if (result)
  {
    *value = static_cast<int>(temp);
  }
  return static_cast<int>(result);
}


/**
 * Reads the value for the attribute with the given XMLTriple_t into value.  
 * If the XMLTriple_t was not found or value could be interpreted as a boolean, 
 * value is not modified.
 *
 * According to the W3C XML Schema, valid boolean values are: "true",
 * "false", "1", and "0" (case-insensitive).  For more information, see:
 * http://www.w3.org/TR/xmlschema-2/#boolean
 *
 * If an XMLErrorLog is passed in datatype format errors are logged.  If
 * required is true, missing attributes are also logged.
 *
 *
 * @param xa the XMLAttributes_t structure.
 * @param triple an XMLTriple_t, the XML triple of the attribute.
 * @param value a boolean, the value of the attribute.
 * @param log an XMLErrorLog_t, the error log.
 * @param required a boolean, indicating whether the attribute is required.
 *
 * @returns @c non-zero (true) if the attribute was read into value, 
 * @c zero (false) otherwise.
 **/
LIBLAX_EXTERN
int
XMLAttributes_readIntoBooleanByTriple (XMLAttributes_t *xa,
                               const XMLTriple_t* triple,
                               int *value,
                               XMLErrorLog_t *log,
                               int required)
{
  bool temp;
  bool result = xa->readInto(*triple, temp, log, required);
  if (result)
  {
    *value = static_cast<int>(temp);
  }
  return static_cast<int>(result);
}


/**
 * Reads the value for the attribute name into value.  If the given local
 * name was not found or value could be interpreted as a double, value is 
 * not modified.
 *
 * According to the W3C XML Schema, valid doubles are the same as valid
 * doubles for C and the special values "INF", "-INF", and "NaN"
 * (case-sensitive).  For more information, see:
 * http://www.w3.org/TR/xmlschema-2/#double
 *
 * If an XMLErrorLog is passed in datatype format errors are logged.  If
 * required is true, missing attributes are also logged.
 *
 *
 * @param xa the XMLAttributes_t structure.
 * @param name a string, the local name of the attribute.
 * @param value a boolean, the value of the attribute.
 * @param log an XMLErrorLog_t, the error log.
 * @param required a boolean, indicating whether the attribute is required.
 *
 * @returns @c non-zero (true) if the attribute was read into value, 
 * @c zero (false) otherwise.
 *
 * @note A namespace bound to the given local name is not checked by this 
 * function. readIntoDoubleByTriple(...) should be used to read a value for 
 * an attribute name with a prefix and namespace.
 **/
LIBLAX_EXTERN
int
XMLAttributes_readIntoDouble (XMLAttributes_t *xa,
			      const char *name,
			      double *value,
			      XMLErrorLog_t *log,
			      int required)
{
  return static_cast<int>( xa->readInto(name, *(value), log, required) );
}


/**
 * Reads the value for the attribute with the given XMLTriple_t into value.  
 * If the XMLTriple_t was not found or value could be interpreted as a double,
 * value is not modified.
 *
 * According to the W3C XML Schema, valid doubles are the same as valid
 * doubles for C and the special values "INF", "-INF", and "NaN"
 * (case-sensitive).  For more information, see:
 * http://www.w3.org/TR/xmlschema-2/#double
 *
 * If an XMLErrorLog is passed in datatype format errors are logged.  If
 * required is true, missing attributes are also logged.
 *
 *
 * @param xa the XMLAttributes_t structure.
 * @param triple an XMLTriple_t, the XML triple of the attribute.
 * @param value a boolean, the value of the attribute.
 * @param log an XMLErrorLog_t, the error log.
 * @param required a boolean, indicating whether the attribute is required.
 *
 * @returns @c non-zero (true) if the attribute was read into value, 
 * @c zero (false) otherwise.
 **/
LIBLAX_EXTERN
int
XMLAttributes_readIntoDoubleByTriple (XMLAttributes_t *xa,
                               const XMLTriple_t* triple,
                               double *value,
                               XMLErrorLog_t *log,
                               int required)
{
  return static_cast<int>( xa->readInto(*triple, *(value), log, required) );
}


/**
 * Reads the value for the attribute name into value.  If the given local
 * name was not found or value could be interpreted as a long, value is not
 * modified.
 *
 * According to the W3C XML Schema valid integers include zero, *all*
 * positive and *all* negative whole numbers.  For practical purposes, we
 * limit values to what can be stored in a long.  For more information,
 * see: http://www.w3.org/TR/xmlschema-2/#integer
 *
 * If an XMLErrorLog is passed in datatype format errors are logged.  If
 * required is true, missing attributes are also logged.
 *
 *
 * @param xa the XMLAttributes_t structure.
 * @param name a string, the local name of the attribute.
 * @param value a boolean, the value of the attribute.
 * @param log an XMLErrorLog_t, the error log.
 * @param required a boolean, indicating whether the attribute is required.
 *
 * @returns @c non-zero (true) if the attribute was read into value, 
 * @c zero (false) otherwise.
 *
 * @note A namespace bound to the given local name is not checked by this 
 * function. readIntoLongByTriple(...) should be used to read a value for 
 * an attribute name with a prefix and namespace.
 **/
LIBLAX_EXTERN
int
XMLAttributes_readIntoLong (XMLAttributes_t *xa,
			    const char *name,
			    long *value,
			    XMLErrorLog_t *log,
			    int required)
{
  return static_cast<int>( xa->readInto(name, *(value), log, required) );
}


/**
 * Reads the value for the attribute with the given XMLTriple_t into value.  
 * If the XMLTriple_t was not found or value could be interpreted as a long,
 * value is not modified.
 *
 * According to the W3C XML Schema valid integers include zero, *all*
 * positive and *all* negative whole numbers.  For practical purposes, we
 * limit values to what can be stored in a long.  For more information,
 * see: http://www.w3.org/TR/xmlschema-2/#integer
 *
 * If an XMLErrorLog is passed in datatype format errors are logged.  If
 * required is true, missing attributes are also logged.
 *
 *
 * @param xa the XMLAttributes_t structure.
 * @param triple an XMLTriple_t, the XML triple of the attribute.
 * @param value a boolean, the value of the attribute.
 * @param log an XMLErrorLog_t, the error log.
 * @param required a boolean, indicating whether the attribute is required.
 *
 * @returns @c non-zero (true) if the attribute was read into value, 
 * @c zero (false) otherwise.
 **/
LIBLAX_EXTERN
int
XMLAttributes_readIntoLongByTriple (XMLAttributes_t *xa,
                               const XMLTriple_t* triple,
                               long *value,
                               XMLErrorLog_t *log,
                               int required)
{
  return static_cast<int>( xa->readInto(*triple, *(value), log, required) );
}


/**
 * Reads the value for the attribute name into value.  If the given local
 * name was not found or value could be interpreted as an integer, value 
 * is not modified.
 *
 * According to the W3C XML Schema valid integers include zero, *all*
 * positive and *all* negative whole numbers.  For practical purposes, we
 * limit values to what can be stored in a int.  For more information,
 * see: http://www.w3.org/TR/xmlschema-2/#integer
 *
 * If an XMLErrorLog is passed in datatype format errors are logged.  If
 * required is true, missing attributes are also logged.
 *
 *
 * @param xa the XMLAttributes_t structure.
 * @param name a string, the local name of the attribute.
 * @param value a boolean, the value of the attribute.
 * @param log an XMLErrorLog_t, the error log.
 * @param required a boolean, indicating whether the attribute is required.
 *
 * @returns @c non-zero (true) if the attribute was read into value, 
 * @c zero (false) otherwise.
 *
 * @note A namespace bound to the given local name is not checked by this 
 * function. readIntoIntByTriple(...) should be used to read a value for 
 * an attribute name with a prefix and namespace.
 **/
LIBLAX_EXTERN
int
XMLAttributes_readIntoInt (XMLAttributes_t *xa,
			   const char *name,
			   int *value,
			   XMLErrorLog_t *log,
			   int required)
{
  return static_cast<int>( xa->readInto(name, *(value), log, required) );
}

/**
 * Reads the value for the attribute with the given XMLTriple_t into value.  
 * If the XMLTriple_t was not found or value could be interpreted as an integer,
 * value is not modified.
 *
 * According to the W3C XML Schema valid integers include zero, *all*
 * positive and *all* negative whole numbers.  For practical purposes, we
 * limit values to what can be stored in a int.  For more information,
 * see: http://www.w3.org/TR/xmlschema-2/#integer
 *
 * If an XMLErrorLog is passed in datatype format errors are logged.  If
 * required is true, missing attributes are also logged.
 *
 *
 * @param xa the XMLAttributes_t structure.
 * @param triple an XMLTriple_t, the XML triple of the attribute.
 * @param value a boolean, the value of the attribute.
 * @param log an XMLErrorLog_t, the error log.
 * @param required a boolean, indicating whether the attribute is required.
 *
 * @returns @c non-zero (true) if the attribute was read into value, 
 * @c zero (false) otherwise.
 **/
LIBLAX_EXTERN
int
XMLAttributes_readIntoIntByTriple (XMLAttributes_t *xa,
                               const XMLTriple_t* triple,
                               int *value,
                               XMLErrorLog_t *log,
                               int required)
{
  return static_cast<int>( xa->readInto(*triple, *(value), log, required) );
}


/**
 * Reads the value for the attribute name into value.  If the given local
 * name was not found or value could be interpreted as an unsigned int, 
 * value is not modified.
 *
 * According to the W3C XML Schema valid integers include zero, *all*
 * positive and *all* negative whole numbers.  For practical purposes, we
 * limit values to what can be stored in a unsigned int.  For more
 * information, see: http://www.w3.org/TR/xmlschema-2/#integer
 *
 * If an XMLErrorLog is passed in datatype format errors are logged.  If
 * required is true, missing attributes are also logged.
 *
 *
 * @param xa the XMLAttributes_t structure.
 * @param name a string, the local name of the attribute.
 * @param value an unsigned int, the value of the attribute.
 * @param log an XMLErrorLog_t, the error log.
 * @param required a boolean, indicating whether the attribute is required.
 *
 * @returns @c non-zero (true) if the attribute was read into value, 
 * @c zero (false) otherwise.
 *
 * @note A namespace bound to the given local name is not checked by this 
 * function. readIntoUnsignedIntByTriple(...) should be used to read a value for 
 * an attribute name with a prefix and namespace.
 **/
LIBLAX_EXTERN
int
XMLAttributes_readIntoUnsignedInt (XMLAttributes_t *xa,
				   const char *name,
				   unsigned int *value,
				   XMLErrorLog_t *log,
				   int required)
{
  return static_cast<int>( xa->readInto(name, *(value), log, required) );
}


/**
 * Reads the value for the attribute with the given XMLTriple_t into value.  
 * If the XMLTriple_t was not found or value could be interpreted as an unsigned 
 * integer, value is not modified.
 *
 * According to the W3C XML Schema valid integers include zero, *all*
 * positive and *all* negative whole numbers.  For practical purposes, we
 * limit values to what can be stored in a unsigned int.  For more
 * information, see: http://www.w3.org/TR/xmlschema-2/#integer
 *
 * If an XMLErrorLog is passed in datatype format errors are logged.  If
 * required is true, missing attributes are also logged.
 *
 *
 * @param xa the XMLAttributes_t structure.
 * @param triple an XMLTriple_t, the XML triple of the attribute.
 * @param value an unsigned int, the value of the attribute.
 * @param log an XMLErrorLog_t, the error log.
 * @param required a boolean, indicating whether the attribute is required.
 *
 * @returns @c non-zero (true) if the attribute was read into value, 
 * @c zero (false) otherwise.
 **/
LIBLAX_EXTERN
int
XMLAttributes_readIntoUnsignedIntByTriple (XMLAttributes_t *xa,
                               const XMLTriple_t* triple,
                               unsigned int *value,
                               XMLErrorLog_t *log,
                               int required)
{
  return static_cast<int>( xa->readInto(*triple, *(value), log, required) );
}



/**
 * Reads the value for the attribute name into value.  If the given local
 * name was not found, value is not modified.
 *
 * If an XMLErrorLog is passed in datatype format errors are logged.  If
 * required is true, missing attributes are also logged.
 *
 *
 * @param xa the XMLAttributes_t structure.
 * @param name a string, the local name of the attribute.
 * @param value a string, the value of the attribute.
 * @param log an XMLErrorLog_t, the error log.
 * @param required a boolean, indicating whether the attribute is required.
 *
 * @returns @c non-zero (true) if the attribute was read into value, 
 * @c zero (false) otherwise.
 *
 * @note A namespace bound to the given local name is not checked by this 
 * function. readIntoStringByTriple(...) should be used to read a value for 
 * an attribute name with a prefix and namespace.
 **/
LIBLAX_EXTERN
int
XMLAttributes_readIntoString (XMLAttributes_t *xa,
			      const char *name,
			      char **value,
			      XMLErrorLog_t *log,
			      int required)
{
  std::string temp;
  int result = static_cast<int>( xa->readInto(name, temp, log, required) );
  if(result)
  {
    *value = safe_strdup(temp.c_str());
  }
  return result;
}


/**
 * Reads the value for the attribute with the given XMLTriple_t into value.  
 * If the XMLTriple_t was not found, value is not modified.
 *
 * If an XMLErrorLog is passed in datatype format errors are logged.  If
 * required is true, missing attributes are also logged.
 *
 *
 * @param xa the XMLAttributes_t structure.
 * @param triple an XMLTriple_t, the XML triple of the attribute.
 * @param value a string, the value of the attribute.
 * @param log an XMLErrorLog_t, the error log.
 * @param required a boolean, indicating whether the attribute is required.
 *
 * @returns @c non-zero (true) if the attribute was read into value, 
 * @c zero (false) otherwise.
 **/
LIBLAX_EXTERN
int
XMLAttributes_readIntoStringByTriple (XMLAttributes_t *xa,
                               const XMLTriple_t* triple,
                               char **value,
                               XMLErrorLog_t *log,
                               int required)
{
  std::string temp;
  int result = static_cast<int>( xa->readInto(*triple, temp, log, required) );
  if(result)
  {
    *value = safe_strdup(temp.c_str());
  }
  return result;
}


LIBSBML_CPP_NAMESPACE_END


/** @endcond */
