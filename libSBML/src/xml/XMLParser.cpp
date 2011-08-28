/**
 * @file    XMLParser.cpp
 * @brief   XMLParser interface and factory
 * @author  Ben Bornstein
 * @author  Michael Hucka
 *
 * $Id: XMLParser.cpp 11633 2010-08-03 03:53:27Z mhucka $
 * $HeadURL: https://sbml.svn.sourceforge.net/svnroot/sbml/trunk/libsbml/src/xml/XMLParser.cpp $
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


#ifdef USE_EXPAT
#include <sbml/xml/ExpatParser.h>
#endif

#ifdef USE_LIBXML
#include <sbml/xml/LibXMLParser.h>
#endif

#ifdef USE_XERCES
#include <sbml/xml/XercesParser.h>
#endif

#include <sbml/xml/XMLErrorLog.h>
#include <sbml/xml/XMLParser.h>

/** @cond doxygen-ignored */

using namespace std;

/** @endcond */

LIBSBML_CPP_NAMESPACE_BEGIN

/*
 * Creates a new XMLParser.  The parser will notify the given XMLHandler
 * of parse events and errors.
 */
XMLParser::XMLParser () : mErrorLog(0)
{
}


/*
 * Destroys this XMLParser.
 */
XMLParser::~XMLParser ()
{
}


/*
 * Creates a new XMLParser.  The parser will notify the given XMLHandler
 * of parse events and errors.
 *
 * The library parameter indicates the underlying XML library to use if
 * the XML compatibility layer has been linked against multiple XML
 * libraries.  It may be one of: "expat" (default), "libxml", or
 * "xerces".
 *
 * If the XML compatibility layer has been linked against only a single
 * XML library, the library parameter is ignored.
 */
XMLParser*
XMLParser::create (XMLHandler& handler, const string library)
{
#ifdef USE_EXPAT
  if (library.empty() || library == "expat")  return new ExpatParser(handler);
#endif

#ifdef USE_LIBXML
  if (library.empty() || library == "libxml") return new LibXMLParser(handler);
#endif

#ifdef USE_XERCES
  if (library.empty() || library == "xerces") return new XercesParser(handler);
#endif

  return 0;
}


/*
 * @return an XMLErrorLog which can be used to log XML parse errors and
 * other validation errors (and messages).
 */
XMLErrorLog*
XMLParser::getErrorLog ()
{
  return mErrorLog;
}


/*
 * Sets the XMLErrorLog this parser will use to log errors.
 */
int
XMLParser::setErrorLog (XMLErrorLog* log)
{
  mErrorLog = log;
  if (mErrorLog) 
  {
    return mErrorLog->setParser(this);
  }
  else
  {
    return LIBSBML_OPERATION_FAILED;
  }
}

LIBSBML_CPP_NAMESPACE_END

