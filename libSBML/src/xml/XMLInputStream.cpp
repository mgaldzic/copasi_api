/**
 * @file    XMLInputStream.cpp
 * @brief   XMLInputStream
 * @author  Ben Bornstein
 *
 * $Id: XMLInputStream.cpp 11633 2010-08-03 03:53:27Z mhucka $
 * $HeadURL: https://sbml.svn.sourceforge.net/svnroot/sbml/trunk/libsbml/src/xml/XMLInputStream.cpp $
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

#include <sbml/xml/XMLErrorLog.h>
#include <sbml/xml/XMLParser.h>

#include <sbml/xml/XMLInputStream.h>

/** @cond doxygen-ignored */

using namespace std;

/** @endcond */

LIBSBML_CPP_NAMESPACE_BEGIN

/** @cond doxygen-libsbml-internal */

class XMLParser;


/*
 * Creates a new XMLInputStream.
 */
XMLInputStream::XMLInputStream (  const char*   content
                                , bool          isFile
                                , const std::string  library 
		                , XMLErrorLog*  errorLog ) :


   mIsError ( false )
 , mParser  ( XMLParser::create( mTokenizer, library) )
 , mSBMLns  ( 0 )
{
  if ( !isGood() ) return;
  if ( errorLog ) setErrorLog(errorLog);
  // if this fails we should probably flag the stream as error
  if (!mParser->parseFirst(content, isFile))
    mIsError = true; 
}


/*
 * Destroys this XMLInputStream.
 */
XMLInputStream::~XMLInputStream ()
{
  if ( mParser )
  {
     /**
      *  set NULL to 'XMLErrorLog::mParser' (corresponding XMLErrorLog* 
      *  object was passed to the above constructer as 'errorLog') because 
      *  the corresponding 'mParser' is deleted here and can't be accessed 
      *  anymore.    
	   */
    XMLErrorLog* errorLog = mParser->getErrorLog();
    if ( errorLog ) errorLog->setParser(NULL);
  }
  delete mParser;
  delete mSBMLns;
}


/*
 * @return the encoding of the XML stream.
 */
const string&
XMLInputStream::getEncoding ()
{
  return mTokenizer.getEncoding();
}


/*
 * @return the version of the XML stream.
 */
const string&
XMLInputStream::getVersion ()
{
  return mTokenizer.getVersion();
}


/*
 * @return an XMLErrorLog which can be used to log XML parse errors and
 * other validation errors (and messages).
 */
XMLErrorLog*
XMLInputStream::getErrorLog ()
{
  return mParser->getErrorLog();
}


/*
 * @return true if end of file (stream) has been reached, false otherwise.
 */
bool
XMLInputStream::isEOF () const
{
  return mTokenizer.isEOF();
}


/*
 * @return true if a fatal error occurred while reading from this stream.
 */
bool
XMLInputStream::isError () const
{
  return (mIsError || mParser == 0);
}


/*
 * @return true if the stream is in a good state (i.e. isEOF() and
 * isError() are both false), false otherwise.
 */
bool
XMLInputStream::isGood () const
{
  return (isError() == false && isEOF() == false);
}


/*
 * Consumes the next XMLToken and return it.
 *
 * @return the next XMLToken or EOF (XMLToken.isEOF() == true).
 */
XMLToken
XMLInputStream::next ()
{
  queueToken();
  return mTokenizer.hasNext() ? mTokenizer.next() : XMLToken();
}


/*
 * Returns the next XMLToken without consuming it.  A subsequent call to
 * either peek() or next() will return the same token.
 *
 * @return the next XMLToken or EOF (XMLToken.isEOF() == true).
 */
const XMLToken&
XMLInputStream::peek ()
{
  queueToken();
  return mTokenizer.hasNext() ? mTokenizer.peek() : mEOF;
}


/*
 * Runs mParser until mTokenizer is ready to deliver at least one XMLToken
 * or a fatal error occurs.
 */
void
XMLInputStream::queueToken ()
{
  if ( !isGood() ) return;

  bool success = true;

  while ( success && mTokenizer.hasNext() == false )
  {
    success = mParser->parseNext();
  }

  if (success == false && isEOF() == false)
  {
    mIsError = true;
  }
}


/*
 * Sets the XMLErrorLog this stream will use to log errors.
 */
int
XMLInputStream::setErrorLog (XMLErrorLog* log)
{
  return mParser->setErrorLog(log);
}


/*
 * Consume zero or more XMLTokens up to and including the corresponding
 * end XML element or EOF.
 */
void
XMLInputStream::skipPastEnd (const XMLToken& element)
{
  if ( element.isEnd() ) return;

  while ( isGood() && !peek().isEndFor(element) ) next();
  next();
}


/*
 * Consume zero or more XMLTokens up to but not including the next XML
 * element or EOF.
 */
void
XMLInputStream::skipText ()
{
  while ( isGood() && peek().isText() ) next();
}


/*
 * Prints a string representation of the underlying token stream, for
 * debugging purposes.
 */
string
XMLInputStream::toString ()
{
  return mTokenizer.toString();
}

SBMLNamespaces *
XMLInputStream::getSBMLNamespaces()
{
  return mSBMLns;
}

void
XMLInputStream::setSBMLNamespaces(SBMLNamespaces * sbmlns)
{
  if (sbmlns)
    mSBMLns = sbmlns->clone();
  else
    mSBMLns = 0;
}
/**
 * 
 **/
LIBLAX_EXTERN
XMLInputStream_t *
XMLInputStream_create (const char* content, int isFile, const char *library)
{
  return new(nothrow) XMLInputStream(content, isFile, library);
}


/**
 * 
 **/
LIBLAX_EXTERN
void
XMLInputStream_free (XMLInputStream_t *stream)
{
  delete static_cast<XMLInputStream*>(stream);
}  


/**
 * 
 **/
LIBLAX_EXTERN
const char *
XMLInputStream_getEncoding (XMLInputStream_t *stream)
{
  return stream->getEncoding().empty() ? NULL : stream->getEncoding().c_str();
}


/**
 * 
 **/
LIBLAX_EXTERN
XMLErrorLog_t *
XMLInputStream_getErrorLog (XMLInputStream_t *stream)
{
  return stream->getErrorLog();
}


/**
 * 
 **/
LIBLAX_EXTERN
int
XMLInputStream_isEOF (XMLInputStream_t *stream)
{
  return static_cast<int>(stream->isEOF());
}


/**
 * 
 **/
LIBLAX_EXTERN
int
XMLInputStream_isError (XMLInputStream_t *stream)
{
  return static_cast<int>(stream->isError());
}


/**
 * 
 **/
LIBLAX_EXTERN
int
XMLInputStream_isGood (XMLInputStream_t *stream)
{
  return static_cast<int>(stream->isGood());
}


/**
 * 
 **/
LIBLAX_EXTERN
XMLToken_t *
XMLInputStream_next (XMLInputStream_t *stream)
{
  return new (nothrow) XMLToken(stream->next());
}


/**
 * 
 **/
LIBLAX_EXTERN
const XMLToken_t *
XMLInputStream_peek (XMLInputStream_t *stream)
{
  return &(stream->peek());
}


/**
 * 
 **/
LIBLAX_EXTERN
void
XMLInputStream_skipPastEnd (XMLInputStream_t *stream,
			    const XMLToken_t *element)
{
  stream->skipPastEnd(*element);
}


/**
 * 
 **/
LIBLAX_EXTERN
void
XMLInputStream_skipText (XMLInputStream_t *stream)
{
  stream->skipText();
}


/**
 * 
 **/
LIBLAX_EXTERN
int
XMLInputStream_setErrorLog (XMLInputStream_t *stream, XMLErrorLog_t *log)
{
  return stream->setErrorLog(log);
}

/** @endcond */

LIBSBML_CPP_NAMESPACE_END

