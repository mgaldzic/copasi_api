/**
 * @file    XercesTranscode.h
 * @brief   Transcodes a Xerces-C++ XMLCh* string to the local code page.
 * @author  Ben Bornstein
 *
 * $Id: XercesTranscode.h 11633 2010-08-03 03:53:27Z mhucka $
 * $HeadURL: https://sbml.svn.sourceforge.net/svnroot/sbml/trunk/libsbml/src/xml/XercesTranscode.h $
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

#ifndef XercesTranscode_h
#define XercesTranscode_h

#ifdef __cplusplus

#include <string>
#include <xercesc/util/XMLString.hpp>
#include <sbml/common/libsbml-namespace.h>

LIBSBML_CPP_NAMESPACE_BEGIN

#if XERCES_VERSION_MAJOR <= 2
typedef unsigned int XercesSize_t;
typedef XMLSSize_t   XercesFileLoc;
#else
typedef XMLSize_t    XercesSize_t;
typedef XMLFileLoc   XercesFileLoc;
#endif

/** @cond doxygen-libsbml-internal */

/**
 * Transcodes a Xerces-C++ XMLCh* string to the UTF8 string.  This
 * class offers implicit conversion to a C++ string and destroys the
 * transcoded buffer when the XercesTranscode object goes out of scope.
 */
class XercesTranscode
{
public:

  XercesTranscode (const XMLCh* s) :
    mBuffer( transcodeToUTF8(s) ) { }

  ~XercesTranscode     () { delete [] mBuffer; }
  operator std::string () { return std::string(mBuffer); }


private:

  char* mBuffer;

  XercesTranscode  ();
  XercesTranscode  (const XercesTranscode&);
  XercesTranscode& operator= (const XercesTranscode&);

 /**
  * convert the given internal XMLCh* string to the UTF-8 char* string.
  */
  char* transcodeToUTF8(const XMLCh* src_str);

};

/** @endcond */

LIBSBML_CPP_NAMESPACE_END

#endif  /* __cplusplus */
#endif  /* XercesTranscode_h */
