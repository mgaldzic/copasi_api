/**
 * @file    XMLError.cpp
 * @brief   Represents errors (and messages) encountered during an XML parse
 * @author  Ben Bornstein
 *
 * $Id: XMLError.cpp 11633 2010-08-03 03:53:27Z mhucka $
 * $HeadURL: https://sbml.svn.sourceforge.net/svnroot/sbml/trunk/libsbml/src/xml/XMLError.cpp $
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

#include <iostream>
#include <iomanip>
#include <sstream>

#include <sbml/xml/XMLError.h>

/** @cond doxygen-ignored */

using namespace std;

/** @endcond */

LIBSBML_CPP_NAMESPACE_BEGIN

/** @cond doxygen-libsbml-internal **/

static const xmlErrorTableEntry errorTable[] =
{
  // 0
  { XMLUnknownError, LIBSBML_CAT_INTERNAL, LIBSBML_SEV_FATAL,
    "Unknown error",
    "Unrecognized error encountered internally." },


  // System diagnostics:

  // 0001
  { XMLOutOfMemory, LIBSBML_CAT_SYSTEM, LIBSBML_SEV_FATAL,
    "Out of memory",
    "Out of memory." },

  // 0002
  { XMLFileUnreadable, LIBSBML_CAT_SYSTEM, LIBSBML_SEV_ERROR,
    "File unreadable",
    "File unreadable." },

  // 0003
  { XMLFileUnwritable, LIBSBML_CAT_SYSTEM, LIBSBML_SEV_ERROR,
    "File unwritable",
    "File unwritable." },

  // 0004
  { XMLFileOperationError, LIBSBML_CAT_SYSTEM, LIBSBML_SEV_ERROR,
    "File operation error",
    "Error encountered while attempting file operation." },

  // 0005
  { XMLNetworkAccessError, LIBSBML_CAT_SYSTEM, LIBSBML_SEV_ERROR,
    "Network access error",
    "Network access error." },


  // Internal diagnostics:

  // 0101
  { InternalXMLParserError, LIBSBML_CAT_INTERNAL, LIBSBML_SEV_FATAL,
    "Internal XML parser error",
    "Internal XML parser state error." },

  // 0102
  { UnrecognizedXMLParserCode, LIBSBML_CAT_INTERNAL, LIBSBML_SEV_FATAL,
    "Unrecognized XML parser code",
    "XML parser returned an unrecognized error code." },

  // 0102
  { XMLTranscoderError, LIBSBML_CAT_INTERNAL, LIBSBML_SEV_FATAL,
    "Transcoder error",
    "Character transcoder error." },


  // Diagnostics about XML content:

  // 1001
  { MissingXMLDecl, LIBSBML_CAT_XML, LIBSBML_SEV_ERROR,
    "Missing XML declaration",
    "Missing XML declaration at beginning of XML input." },

  // 1002
  { MissingXMLEncoding, LIBSBML_CAT_XML, LIBSBML_SEV_ERROR,
    "Missing XML encoding attribute",
    "Missing encoding attribute in XML declaration." },

  // 1003
  { BadXMLDecl, LIBSBML_CAT_XML, LIBSBML_SEV_ERROR,
    "Bad XML declaration",
    "Invalid or unrecognized XML declaration or XML encoding." },

  // 1004
  { BadXMLDOCTYPE, LIBSBML_CAT_XML, LIBSBML_SEV_ERROR,
    "Bad XML DOCTYPE",
    "Invalid, malformed or unrecognized XML DOCTYPE declaration." },

  // 1005
  { InvalidCharInXML, LIBSBML_CAT_XML, LIBSBML_SEV_ERROR,
    "Invalid character",
    "Invalid character in XML content." },

  // 1006
  { BadlyFormedXML, LIBSBML_CAT_XML, LIBSBML_SEV_ERROR,
    "Badly formed XML",
    "XML content is not well-formed." },

  // 1007
  { UnclosedXMLToken, LIBSBML_CAT_XML, LIBSBML_SEV_ERROR,
    "Unclosed token",
    "Unclosed XML token." },

  // 1008
  { InvalidXMLConstruct, LIBSBML_CAT_XML, LIBSBML_SEV_ERROR,
    "Invalid XML construct",
    "XML construct is invalid or not permitted." },

  // 1009
  { XMLTagMismatch, LIBSBML_CAT_XML, LIBSBML_SEV_ERROR,
    "XML tag mismatch",
    "Element tag mismatch or missing tag." },

  // 1010
  { DuplicateXMLAttribute, LIBSBML_CAT_XML, LIBSBML_SEV_ERROR,
    "Duplicate attribute",
    "Duplicate XML attribute." },

  // 1011
  { UndefinedXMLEntity, LIBSBML_CAT_XML, LIBSBML_SEV_ERROR,
    "Undefined XML entity",
    "Undefined XML entity." },

  // 1012
  { BadProcessingInstruction, LIBSBML_CAT_XML, LIBSBML_SEV_ERROR,
    "Bad XML processing instruction",
    "Invalid, malformed or unrecognized XML processing instruction." },

  // 1013
  { BadXMLPrefix, LIBSBML_CAT_XML, LIBSBML_SEV_ERROR,
    "Bad XML prefix",
    "Invalid or undefined XML namespace prefix." },

  // 1014
  { BadXMLPrefixValue, LIBSBML_CAT_XML, LIBSBML_SEV_ERROR,
    "Bad XML prefix value",
    "Invalid XML namespace prefix value." },

  // 1015
  { MissingXMLRequiredAttribute, LIBSBML_CAT_XML, LIBSBML_SEV_ERROR,
    "Missing required attribute",
    "Missing a required XML attribute." },

  // 1016
  { XMLAttributeTypeMismatch, LIBSBML_CAT_XML, LIBSBML_SEV_ERROR,
    "Attribute type mismatch",
    "Data type mismatch for the value of an attribute." },

  // 1017
  { XMLBadUTF8Content, LIBSBML_CAT_XML, LIBSBML_SEV_ERROR,
    "Bad UTF8 content",
    "Invalid UTF8 content." },

  // 1018
  { MissingXMLAttributeValue, LIBSBML_CAT_XML, LIBSBML_SEV_ERROR,
    "Missing attribute value",
    "Missing or improperly formed attribute value." },

  // 1019
  { BadXMLAttributeValue, LIBSBML_CAT_XML, LIBSBML_SEV_ERROR,
    "Bad attribute value",
    "Invalid or unrecognizable attribute value." },

  // 1020
  { BadXMLAttribute, LIBSBML_CAT_XML, LIBSBML_SEV_ERROR,
    "Bad XML attribute",
    "Invalid, unrecognized or malformed attribute." },

  // 1021
  { UnrecognizedXMLElement, LIBSBML_CAT_XML, LIBSBML_SEV_ERROR,
    "Unrecognized XML element",
    "Element either not recognized or not permitted." },

  // 1022
  { BadXMLComment, LIBSBML_CAT_XML, LIBSBML_SEV_ERROR,
    "Bad XML comment",
    "Badly formed XML comment." },

  // 1023
  { BadXMLDeclLocation, LIBSBML_CAT_XML, LIBSBML_SEV_ERROR,
    "Bad XML declaration location",
    "XML declaration not permitted in this location." },

  // 1024
  { XMLUnexpectedEOF, LIBSBML_CAT_XML, LIBSBML_SEV_ERROR,
    "Unexpected EOF",
    "Reached end of input unexpectedly." },

  // 1025
  { BadXMLIDValue, LIBSBML_CAT_XML, LIBSBML_SEV_ERROR,
    "Bad XML ID value",
    "Value is invalid for XML ID, or has already been used." },

  // 1026
  { BadXMLIDRef, LIBSBML_CAT_XML, LIBSBML_SEV_ERROR,
    "Bad XML IDREF",
    "XML ID value was never declared." },

  // 1027
  { UninterpretableXMLContent, LIBSBML_CAT_XML, LIBSBML_SEV_ERROR,
    "Uninterpretable XML content",
    "Unable to interpret content." },

  // 1028
  { BadXMLDocumentStructure, LIBSBML_CAT_XML, LIBSBML_SEV_ERROR,
    "Bad XML document structure",
    "Bad XML document structure." },

  // 1029
  { InvalidAfterXMLContent, LIBSBML_CAT_XML, LIBSBML_SEV_ERROR,
    "Invalid content after XML content",
    "Encountered invalid content after expected content." },

  // 1031
  { XMLExpectedQuotedString, LIBSBML_CAT_XML, LIBSBML_SEV_ERROR,
    "Expected quoted string",
    "Expected to find a quoted string." },

  // 1032
  { XMLEmptyValueNotPermitted, LIBSBML_CAT_XML, LIBSBML_SEV_ERROR,
    "Empty value not permitted",
    "An empty value is not permitted in this context." },

  // 1033
  { XMLBadNumber, LIBSBML_CAT_XML, LIBSBML_SEV_ERROR,
    "Bad number",
    "Invalid or unrecognized number." },

  // 1034
  { XMLBadColon, LIBSBML_CAT_XML, LIBSBML_SEV_ERROR,
    "Colon character not permitted",
    "Colon characters are invalid in this context." },

  // 1035
  { MissingXMLElements, LIBSBML_CAT_XML, LIBSBML_SEV_ERROR,
    "Missing XML elements",
    "One or more expected elements are missing." },

  // 1036
  { XMLContentEmpty, LIBSBML_CAT_XML, LIBSBML_SEV_ERROR,
    "Empty XML content",
    "Main XML content is empty." },

};


/*
 * Table of strings corresponding to the values from XMLErrorSeverity_t.
 * Be careful that the enum is used to index into this table.
 * This approach is not very object oriented and someday we should do better.
 */
const char * xmlSeverityStringTable[] =
{
  "Informational"                       // LIBSBML_SEV_INFO
, "Warning"                             // LIBSBML_SEV_WARNING
, "Error"                               // LIBSBML_SEV_ERROR
, "Fatal"                               // LIBSBML_SEV_FATAL
};


/*
 * Table of strings corresponding to the values from XMLErrorCategory_t.
 * Be careful that the enum is used to index into this table.
 * This approach is not very object oriented and someday we should do better.
 */
const char * xmlCategoryStringTable[] =
{
  "Internal"                            // LIBSBML_CAT_INTERNAL
, "Operating system"                    // LIBSBML_CAT_SYSTEM
, "XML content"                         // LIBSBML_CAT_XML
};

/** @endcond **/


/*
 * XMLErrorLog will check if line & column = 0 and attempt to fill in
 * the line and column by consulting the parser.  This constructor
 * purposefully doesn't do that.
 */
XMLError::XMLError (  const int errorId
                    , const std::string& details
                    , const unsigned int line
                    , const unsigned int column
                    , const unsigned int severity
                    , const unsigned int category ) :
    mErrorId( errorId )
  , mLine   ( line    )
  , mColumn ( column  )
{
  // Check if the given id is one we have in our table of error codes.  If
  // it is, fill in the fields of the error object with the appropriate
  // content.  If it's not in the table, take the content as-is.

  if ( errorId >= 0 && errorId < XMLErrorCodesUpperBound )
  {
    unsigned int tableSize = sizeof(errorTable)/sizeof(errorTable[0]);    

    for ( unsigned int i = 0; i < tableSize; i++ )
    {
      if ( errorTable[i].code == errorId )
      {
        mMessage      = errorTable[i].message;
	mShortMessage = errorTable[i].shortMessage;

        if ( !details.empty() )
        {
          mMessage.append(" ");
          mMessage.append(details);
        }

        mSeverity = errorTable[i].severity;
        mCategory = errorTable[i].category;
        
        mSeverityString = stringForSeverity(mSeverity);
        mCategoryString = stringForCategory(mCategory);

        return;
      }
    }

    // The id is in the range of error numbers that are supposed to be in
    // the XML layer, but it's NOT in our table.  This is an internal error.
    // Unfortunately, we don't have an error log or anywhere to report it
    // except the measure of last resort: the standard error output.
    
    cerr << "Internal error: unknown error code '" << errorId
         << "' encountered while processing error" << endl;
  }

  // It's not an error code in the XML layer, so assume the caller has
  // filled in all the relevant additional data.  (If they didn't, the
  // following ends up merely assigning the defaults.)

  mMessage  = details;

  // The following is just a default that seems more sensible than setting
  // an empty string for the short message.

  mShortMessage = mMessage;

  // Other initializations.

  mSeverity       = severity;
  mSeverityString = stringForSeverity(severity);

  mCategory       = category;
  mCategoryString = stringForCategory(category);
}


/*
 * Copy Constructor
 */
XMLError::XMLError(const XMLError& orig) :
    mErrorId(orig.mErrorId)
  , mMessage(orig.mMessage)
  , mShortMessage(orig.mShortMessage)
  , mSeverity(orig.mSeverity)
  , mCategory(orig.mCategory)
  , mLine(orig.mLine)
  , mColumn(orig.mColumn)
  , mSeverityString(orig.mSeverityString)
  , mCategoryString(orig.mCategoryString)
{
}


/*
 * Destroys this XMLError.
 */
XMLError::~XMLError ()
{
}


/** @cond doxygen-libsbml-internal **/
/*
 * Clone function.
 */
XMLError* 
XMLError::clone() const
{
  return new XMLError(*this);
}
/** @endcond **/


/*
 * @return the id of this XMLError.
 */
const unsigned int
XMLError::getErrorId () const
{
  return mErrorId;
}


/*
 * @return the message text of this XMLError.
 */
const string&
XMLError::getMessage () const
{
  return mMessage;
}


/*
 * @return the short message text of this XMLError.
 */
const string&
XMLError::getShortMessage () const
{
  return mShortMessage;
}


/*
 * @return the line number where this XMLError ocurred.
 */
unsigned int
XMLError::getLine () const
{
  return mLine;
}


/*
 * @return the column number where this XMLError occurred.
 */
unsigned int
XMLError::getColumn () const
{
  return mColumn;
}


/*
 * @return the severity of this XMLError.  XMLError severity levels
 * correspond to those defined in the XML specification (with the addition
 * of Info for informational messages).
 */
unsigned int
XMLError::getSeverity () const
{
  return mSeverity;
}


/*
 * @return a string explaining the severity code for this XMLError.
 * XMLError severity levels correspond to those defined in the XML
 * specification (with the addition of Info for informational messages).
 */
const std::string& 
XMLError::getSeverityAsString() const
{
  return mSeverityString;
}


/*
 * @return the category of this XMLError.  A category is a string, similiar
 * in spirit to an XML namespace, which can be used to partition errors
 * into distinct groups.  Among other things, this can be used to prevent
 * id conflicts by uniquely identifying an XMLError by both id and
 * category.
 */
unsigned int
XMLError::getCategory () const
{
  return mCategory;
}


/*
 * @return a string explaining the category code for this XMLError.
 */
const std::string& 
XMLError::getCategoryAsString() const
{
  return mCategoryString;
}


/*
 * @return true if this XMLError is for informational purposes only,
 * false otherwise.
 */
bool
XMLError::isInfo () const
{
  return (mSeverity == LIBSBML_SEV_INFO);
}


/*
 * @return true if this XMLError is a warning, false otherwise.
 */
bool
XMLError::isWarning () const
{
  return (mSeverity == LIBSBML_SEV_WARNING);
}


/*
 * @return true if this XMLError is an error, false otherwise.
 */
bool
XMLError::isError () const
{
  return (mSeverity == LIBSBML_SEV_ERROR);
}


/*
 * @return true if this XMLError is a fatal error, false otherwise.
 */
bool
XMLError::isFatal () const
{
  return (mSeverity == LIBSBML_SEV_FATAL);
}


/*
 * Predicate returning @c true or @c false depending on whether 
 * this XMLError resulted from the operating system.
 *
 * @return @c true or @c false
 */
bool
XMLError::isSystem () const
{
  return (mCategory == LIBSBML_CAT_SYSTEM);
}


/*
 * Predicate returning @c true or @c false depending on whether this
 * XMLError resulted from a problem at the raw XML level (e.g., an XML
 * syntax error).
 *
 * @return @c true or @c false
 */
bool
XMLError::isXML () const
{
  return (mCategory == LIBSBML_CAT_XML);
}


/*
 * Predicate returning @c true or @c false depending on whether this
 * XMLError resulted from an internal program error.
 *
 * @return @c true or @c false
 */
bool
XMLError::isInternal () const
{
  return (mCategory == LIBSBML_CAT_INTERNAL);
}


/*
 * Sets the line number where this XMLError occurred.
 */
int
XMLError::setLine (unsigned int line)
{
  mLine = line;
  return LIBSBML_OPERATION_SUCCESS;
}


/*
 * Sets the column number where this XMLError occurred.
 */
int
XMLError::setColumn (unsigned int column)
{
  mColumn = column;
  return LIBSBML_OPERATION_SUCCESS;
}


/*
 * Given an XMLError::Code, return a copy of the error text.
 * 
 * @return the message text 
 */
const string
XMLError::getStandardMessage (const int code)
{
  string msg;

  if ( code >= 0 && code < XMLErrorCodesUpperBound )
  {
    unsigned int tableSize = sizeof(errorTable)/sizeof(errorTable[0]);    

    for ( unsigned int i = 0; i < tableSize; i++ )
      if ( errorTable[i].code == code )
        msg.append(errorTable[i].message);
  }
  
  return msg;
}


/** @cond doxygen-libsbml-internal **/
/*
 * @return the severity as a string for the given @n code.
 */
std::string XMLError::stringForSeverity(unsigned int code) const
{
  if ( code >= LIBSBML_SEV_INFO && code <= LIBSBML_SEV_FATAL)
    return xmlSeverityStringTable[code];
  else 
    return "";
}
/** @endcond **/


/** @cond doxygen-libsbml-internal **/
/*
 * @return the category as a string for the given @n code.
 */
std::string XMLError::stringForCategory(unsigned int code) const
{
  if ( code >= LIBSBML_CAT_INTERNAL && code <= LIBSBML_CAT_XML )
    return xmlCategoryStringTable[code];
  else 
    return "";
}
/** @endcond **/


/** @cond doxygen-libsbml-internal **/
/*
 * Prints this XMLError to stream @p s in the following format (and followed
 * by a newline):
 *
 *   line: (error_id [severity]) message
 *
 * Implementation note:
 * The following code was originally written before the implementation of
 * xmlSeverityStringTable[] and getSeverityAsString().  It previously used
 * hard-code strings, and for LIBSBML_SEV_INFO, it printed "Advisory"
 * instead of "Information".  However, an inspection of the libSBML code
 * showed that nothing ever used LIBSBML_SEV_INFO!  Therefore, changing the 
 * severity string used here to be the same as what getSeverityAsString()
 * returns should not break any caller code.
 */
void
XMLError::print(std::ostream& s) const
{
  s << "line " << getLine() << ": ("
    << setfill('0') << setw(5) << getErrorId()
    << " [" << getSeverityAsString() << "]) "
    << getMessage() << endl;
}
/** @endcond **/


/*
 * Outputs the given XMLError (or the derived class (e.g. SBMLError) ) to stream
 * by invoking the print function which is implemented as a virtual function in
 * the class.
 *
 */
ostream& operator<< (ostream& s, const XMLError& error)
{
  error.print(s);

  return s;
}


/** @cond doxygen-c-only */


/**
 * Creates a new XMLError to report that something occurred.
 */
LIBLAX_EXTERN
XMLError_t*
XMLError_create (void)
{
  return new(nothrow) XMLError;
}


/**
 * Creates a new XMLError with the identification number
 * and detailed message set.
 *
 * If the identifier is < 10000, it must be one of the predefined XML layer
 * error codes.
 *
 * @param errorId an unsigned int, the identification number of the error.
 * @param message a string, the error message.
 *
 */
LIBLAX_EXTERN
XMLError_t*
XMLError_createWithIdAndMessage (unsigned int errorId, const char * message)
{
  return new(nothrow) XMLError(errorId, message);
}

/*
 * Creates a new XMLError to report that something occurred at the given
 * line and column.  Each XMLError also has an identification number, a
 * category, and a severity level associated with it.
 *
 * @param errorId an unsigned int, the identification number of the error.
 * @param message a string, the error message.
 * @param severity XMLErrorSeverity_t, severity of the error.
 * @param category a string, the category to which the error belongs.
 * @param line an unsigned int, the line number at which the error occurs.
 * @param column an unsigned int, the column number at which the error occurs.
 *
 */
//LIBLAX_EXTERN
//XMLError_t*
//XMLError_createWithAll (unsigned int errorId, const char * message, XMLError_Severity severity,
//                        const char * category, unsigned int line, unsigned int column)
//{
//  XMLError::Severity s;
//  switch (severity)
//  {
//  case Info:
//    s = XMLError::Severity::Info;
//    break;
//  case Warning:
//    s = XMLError::Severity::Warning;
//    break;
//  case Error:
//    s = XMLError::Severity::Error;
//    break;
//  case Fatal:
//    s = XMLError::Severity::Fatal;
//    break;
//  default:
//    s = XMLError::Severity::Error;
//    break;
//  }
//  return new(nothrow) XMLError(errorId, message, s, category, line, column);
//}


/**
 * Frees the given XMLError_t structure.
 *
 * @param error the XMLError_t structure to be freed.
 **/
LIBLAX_EXTERN
void
XMLError_free(XMLError_t* error)
{
  delete static_cast<XMLError*>(error);
}

/**
 * Returns the id of this XMLError.
 *
 * @param error the XMLError_t from which to return the id.
 *
 * @return the id of this XMLError.
 */
LIBLAX_EXTERN
unsigned int
XMLError_getErrorId (const XMLError_t *error)
{
  return error->getErrorId();
}


/**
 * Returns the message text of this XMLError.
 *
 * @param error the XMLError_t from which to return the message.
 *
 * @return the message text of this XMLError.
 */
LIBLAX_EXTERN
const char *
XMLError_getMessage (const XMLError_t *error)
{
  return error->getMessage().empty() ? 0 : error->getMessage().c_str();
}


/**
 * Returns the short message text of this XMLError.
 *
 * @param error the XMLError_t from which to return the short message.
 *
 * @return the message text of this XMLError.
 */
LIBLAX_EXTERN
const char *
XMLError_getShortMessage (const XMLError_t *error)
{
  return error->getShortMessage().empty() ? 0 : error->getShortMessage().c_str();
}


/**
 * Return the line number where this XMLError occurred.
 *
 * @param error the XMLError_t from which to return the line number.
 *
 * @return the line number where this XMLError occurred.
 */
LIBLAX_EXTERN
unsigned int
XMLError_getLine (const XMLError_t *error)
{
  return error->getLine();
}


/**
 * Return the column number where this XMLError occurred.
 *
 * @param error the XMLError_t from which to return the column number.
 *
 * @return the column number where this XMLError occurred.
 */
LIBLAX_EXTERN
unsigned int
XMLError_getColumn (const XMLError_t *error)
{
  return error->getColumn();
}


/**
 * Return the severity of this XMLError.  The possible values (for the XML
 * layer) are those from the enumeration XMLErrorSeverity_t.
 *
 * @param error the XMLError_t from which to return the severity.
 *
 * @return the severity of this XMLError.
 */
LIBLAX_EXTERN
unsigned int
XMLError_getSeverity (const XMLError_t *error)
{
  return error->getSeverity();
}


/**
 * Return the severity of this XMLError as a string.erity_t.
 *
 * @param error the XMLError_t from which to return the severity.
 *
 * @return string representing the severity of this XMLError.
 */
LIBLAX_EXTERN
const char *
XMLError_getSeverityAsString (const XMLError_t *error)
{
  return error->getSeverityAsString().empty() ? 0 : 
                            error->getSeverityAsString().c_str();
}


/**
 * Return the category of this XMLError.  The possible values (for the XML
 * layers) are those from the enumeration XMLErrorCategory_t.
 *
 * @param error the XMLError_t from which to return the category.
 *
 * @return the category of this XMLError.
 */
LIBLAX_EXTERN
unsigned int
XMLError_getCategory (const XMLError_t *error)
{
  return error->getCategory();
}


/**
 * Return the category of this XMLError as a string.
 *
 * @param error the XMLError_t from which to return the category.
 *
 * @return string representing the category of this XMLError.
 */
LIBLAX_EXTERN
const char *
XMLError_getCategoryAsString (const XMLError_t *error)
{
  return error->getCategoryAsString().empty() ? 0 : 
                             error->getCategoryAsString().c_str();
}


/**
 * Predicate returning @c true or @c false depending on whether 
 * this XMLError_t structure is for information only.
 *
 * @param error the XMLError_t.
 *
 * @return @c non-zero (true) if this XMLError is for informational purposes
 * only, @c zero (false) otherwise.
 */
LIBLAX_EXTERN
int
XMLError_isInfo (const XMLError_t *error)
{
  return static_cast<int>( error->isInfo() );
}


/**
 * Predicate returning @c true or @c false depending on whether 
 * this XMLError_t structure is a warning.
 *
 * @param error the XMLError_t.
 *
 * @return @c non-zero (true) if this XMLError is a warning, @c zero (false) otherwise.
 */
LIBLAX_EXTERN
int
XMLError_isWarning (const XMLError_t *error)
{
  return static_cast<int>( error->isWarning() );
}


/**
 * Predicate returning @c true or @c false depending on whether 
 * this XMLError_t structure is an error.
 *
 * @param error the XMLError_t.
 *
 * @return @c non-zero (true) if this XMLError is an error, @c zero (false) otherwise.
 */
LIBLAX_EXTERN
int
XMLError_isError (const XMLError_t *error)
{
  return static_cast<int>( error->isError() );
}


/**
 * Predicate returning @c true or @c false depending on whether 
 * this XMLError_t structure is a fatal error.
 *
 * @param error the XMLError_t.
 *
 * @return @c non-zero (true) if this XMLError is a fatal error, @c zero (false) otherwise.
 */
LIBLAX_EXTERN
int
XMLError_isFatal (const XMLError_t *error)
{
  return static_cast<int>( error->isFatal() );
}


/**
 * Outputs this XMLError to stream in the following format (and
 * followed by a newline):
 *
 *   line: (id) message
 *
 * @param error, the XMLError_t structure to write.
 * @param stream, the stream to write to.
 */
LIBLAX_EXTERN
void
XMLError_print (const XMLError_t *error, FILE *stream)
{
  ostringstream os;
  os << *(static_cast<const XMLError*>(error));

  fprintf(stream, "%s", os.str().c_str());

}

LIBSBML_CPP_NAMESPACE_END

/** @endcond */
