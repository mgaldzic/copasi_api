/**
 * @file    XMLErrorLog.cpp
 * @brief   Stores errors (and messages) encountered while processing XML.
 * @author  Ben Bornstein
 *
 * $Id: XMLErrorLog.cpp 11633 2010-08-03 03:53:27Z mhucka $
 * $HeadURL: https://sbml.svn.sourceforge.net/svnroot/sbml/trunk/libsbml/src/xml/XMLErrorLog.cpp $
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

#include <algorithm>
#include <functional>
#include <sstream>

#include <sbml/xml/XMLError.h>
#include <sbml/xml/XMLParser.h>

#include <sbml/xml/XMLErrorLog.h>

/** @cond doxygen-ignored */

using namespace std;

/** @endcond */

LIBSBML_CPP_NAMESPACE_BEGIN

/** @cond doxygen-libsbml-internal */
/*
 * Creates a new empty XMLErrorLog.
 */
XMLErrorLog::XMLErrorLog ():mParser(NULL)
{
}
/** @endcond */


/** @cond doxygen-libsbml-internal */

/**
 * Used by the Destructor to delete each item in mErrors.
 */
struct Delete : public unary_function<XMLError*, void>
{
  void operator() (XMLError* error) { delete error; }
};


/*
 * Destroys this XMLErrorLog.
 */
XMLErrorLog::~XMLErrorLog ()
{
  for_each( mErrors.begin(), mErrors.end(), Delete() );
}
/** @endcond */


/** @cond doxygen-libsbml-internal */
/*
 * Logs the given XMLError.
 */
void
XMLErrorLog::add (const XMLError& error)
{
  XMLError* cerror;

  try
  {
    cerror = error.clone();
  }
  catch (...)
  {
    // Currently do nothing.
    // An error status would be returned in the 4.x.
    return;
  }

  mErrors.push_back(cerror);

  if (cerror->getLine() == 0 && cerror->getColumn() == 0)
  {
    unsigned int line, column;
    if (mParser != NULL)
    {
      try
      {
        line = mParser->getLine();
        column = mParser->getColumn();
      }
      catch (...)
      {
        line = 1;
        column = 1;
      }
    }
    else
    {
      line = 1;
      column = 1;
    }

    cerror->setLine(line);
    cerror->setColumn(column);
  }
}
/** @endcond */


/** @cond doxygen-libsbml-internal */
/*
 * Logs (copies) the XMLErrors in the given XMLError list to this
 * XMLErrorLog.
 */
void
XMLErrorLog::add (const std::list<XMLError>& errors)
{
  list<XMLError>::const_iterator end = errors.end();
  list<XMLError>::const_iterator iter;

  for (iter = errors.begin(); iter != end; ++iter) add( *iter );
}
/** @endcond */


/*
 * @return the nth XMLError in this log.
 */
const XMLError*
XMLErrorLog::getError (unsigned int n) const
{
  return (n < mErrors.size()) ? mErrors[n] : 0;
}


/*
 * @return the number of errors that have been logged.
 */
unsigned int
XMLErrorLog::getNumErrors () const
{
  return mErrors.size();
}


/*
 * Removes all errors from this log.
 */
void 
XMLErrorLog::clearLog()
{
  for_each( mErrors.begin(), mErrors.end(), Delete() );
  mErrors.clear();
}

/** @cond doxygen-libsbml-internal */
/*
 * Sets the XMLParser for this XMLErrorLog.
 *
 * The XMLParser will be used to obtain the current line and column
 * number as XMLErrors are logged (if they have a line and column number
 * of zero).
 */
int
XMLErrorLog::setParser (const XMLParser* p)
{
  mParser = p;

  if (mParser)
    return LIBSBML_OPERATION_SUCCESS;
  else
    return LIBSBML_OPERATION_FAILED;
}
/** @endcond */


/** @cond doxygen-c-only */


/**
 * Creates a new empty XMLErrorLog_t structure and returns it.
 *
 * @return the new XMLErrorLog_t structure.
 **/
LIBLAX_EXTERN
XMLErrorLog_t *
XMLErrorLog_create (void)
{
  return new(nothrow) XMLErrorLog;
}


/**
 * Frees the given XMLError_t structure.
 *
 * @param log XMLErrorLog_t, the error log to be freed.
 */
LIBLAX_EXTERN
void
XMLErrorLog_free (XMLErrorLog_t *log)
{
    delete static_cast<XMLErrorLog*>(log);
}


/**
 * Logs the given XMLError_t structure.
 *
 * @param log XMLErrorLog_t, the error log to be added to.
 * @param error XMLError_t, the error to be logged.
 */
LIBLAX_EXTERN
void
XMLErrorLog_add (XMLErrorLog_t *log, const XMLError_t *error)
{
  log->add(*error);
}


/**
 * Returns the nth XMLError_t in this log.
 *
 * @param log XMLErrorLog_t, the error log to be queried.
 * @param n unsigned int number of the error to retrieve.
 *
 * @return the nth XMLError_t in this log.
 */
LIBLAX_EXTERN
const XMLError_t *
XMLErrorLog_getError (const XMLErrorLog_t *log, unsigned int n)
{
  return log->getError(n);
}


/**
 * Returns the number of errors that have been logged.
 *
 * @param log XMLErrorLog_t, the error log to be queried.
 *
 * @return the number of errors that have been logged.
 */
LIBLAX_EXTERN
unsigned int
XMLErrorLog_getNumErrors (const XMLErrorLog_t *log)
{
  return log->getNumErrors();
}

/**
 * Removes all errors from this log.
 *
 * @param log XMLErrorLog_t, the error log to be cleared.
 */
LIBLAX_EXTERN
void
XMLErrorLog_clearLog (XMLErrorLog_t *log)
{
  log->clearLog();
}

LIBSBML_CPP_NAMESPACE_END

/** @endcond */
