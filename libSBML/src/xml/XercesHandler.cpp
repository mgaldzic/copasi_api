/**
 * @file    XercesHandler.cpp
 * @brief   Redirect Xerces-C++ SAX2 events to an XMLHandler
 * @author  Ben Bornstein
 *
 * $Id: XercesHandler.cpp 11633 2010-08-03 03:53:27Z mhucka $
 * $HeadURL: https://sbml.svn.sourceforge.net/svnroot/sbml/trunk/libsbml/src/xml/XercesHandler.cpp $
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

#include <xercesc/sax/Locator.hpp>
#include <xercesc/internal/ReaderMgr.hpp>

#include <sbml/xml/XMLHandler.h>
#include <sbml/xml/XMLToken.h>

#include <sbml/xml/XercesAttributes.h>
#include <sbml/xml/XercesNamespaces.h>
#include <sbml/xml/XercesTranscode.h>
#include <sbml/xml/XercesHandler.h>

/** @cond doxygen-ignored */

using namespace std;
using namespace xercesc;

/** @endcond */

LIBSBML_CPP_NAMESPACE_BEGIN

/** @cond doxygen-libsbml-internal */

/**
 * @return the prefix portion of the XML qualified name, or an empty
 * string if no prefix exists.
 */
const string
getPrefix (const string& qname)
{
  string::size_type pos = qname.find(':', 0);
  return (pos != string::npos) ? qname.substr(0, pos) : "";
}


/**
 * Creates a new XercesHandler.  Xerces-C++ SAX2 events will be
 * redirected to the given XMLHandler.
 */
XercesHandler::XercesHandler (XMLHandler& handler) :
   mHandler( handler )
 , mLocator( 0       )
{
}


/**
 * Destroys this XercesHandler.
 */
XercesHandler::~XercesHandler ()
{
}


/**
 * Receive notification of the beginning of the document.
 */
void
XercesHandler::startDocument ()
{
  mHandler.startDocument();
}


/**
 * Receive notification of the start of an element.
 *
 * @param  uri        The URI of the associated namespace for this element
 * @param  localname  The local part of the element name
 * @param  qname      The qualified name of this element
 * @param  attrs      The specified or defaulted attributes
 */
void
XercesHandler::startElement (  const XMLCh* const  uri
                             , const XMLCh* const  localname
                             , const XMLCh* const  qname
                             , const Attributes&   attrs )
{
  const string nsuri  = XercesTranscode( uri       );
  const string name   = XercesTranscode( localname );
  const string prefix = getPrefix( XercesTranscode(qname) );

  const XMLTriple         triple    ( name, nsuri, prefix );
  const XercesAttributes  attributes( attrs, name );
  const XercesNamespaces  namespaces( attrs );
  const XMLToken          element   ( triple, attributes, namespaces,
                                      getLine(), getColumn() );

  mHandler.startElement(element);
}


/**
 * Receive notification of the end of the document.
 */
void
XercesHandler::endDocument ()
{
  mHandler.endDocument();
}


/**
 * Receive notification of the end of an element.
 *
 * @param  uri        The URI of the associated namespace for this element
 * @param  localname  The local part of the element name
 * @param  qname      The qualified name of this element
 */
void
XercesHandler::endElement (  const XMLCh* const  uri
                           , const XMLCh* const  localname
                           , const XMLCh* const  qname )
{
  const string nsuri  = XercesTranscode( uri       );
  const string name   = XercesTranscode( localname );
  const string prefix = getPrefix( XercesTranscode(qname) );

  const XMLTriple  triple ( name, nsuri, prefix );
  const XMLToken   element( triple, getLine(), getColumn() );

  mHandler.endElement(element);
}


/**
 * Receive notification of character data inside an element.
 *
 * @param  chars   The characters
 * @param  length  The number of characters to use from the character array
 */
void
XercesHandler::characters (  const XMLCh* const  chars
                           , const XercesSize_t  length )
{
  const string   transcoded = XercesTranscode(chars);
  const XMLToken data(transcoded);

  mHandler.characters(transcoded);
}


/**
 * @return the column number of the current XML event.
 */
unsigned int
XercesHandler::getColumn () const
{
  unsigned int column = 0;


  if (mLocator
  &&  static_cast<const xercesc::ReaderMgr*>(mLocator)->getCurrentReader()
  &&  mLocator->getColumnNumber() > 0)
  {
    column = static_cast<unsigned int>( mLocator->getColumnNumber() );
  }

  return column;
}


/**
 * @return the line number of the current XML event.
 */
unsigned int
XercesHandler::getLine () const
{
  unsigned int line = 0;


  if (mLocator 
  &&  static_cast<const xercesc::ReaderMgr*>(mLocator)->getCurrentReader() 
  &&  mLocator->getLineNumber() > 0)
  {
    line = static_cast<unsigned int>( mLocator->getLineNumber() );
  }

  return line;
}


/**
 * Receive a Locator object for document events.
 */
void
XercesHandler::setDocumentLocator (const xercesc::Locator* const locator)
{
  mLocator = locator;
}

/** @endcond */

LIBSBML_CPP_NAMESPACE_END
