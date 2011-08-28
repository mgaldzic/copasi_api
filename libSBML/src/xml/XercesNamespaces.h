/**
 * @file    XercesNamespaces.h
 * @brief   Extracts XML namespace declarations from Xerces-C++ attributes.
 * @author  Ben Bornstein
 *
 * $Id: XercesNamespaces.h 11633 2010-08-03 03:53:27Z mhucka $
 * $HeadURL: https://sbml.svn.sourceforge.net/svnroot/sbml/trunk/libsbml/src/xml/XercesNamespaces.h $
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

#ifndef XercesNamespaces_h
#define XercesNamespaces_h

#ifdef __cplusplus

#include <string>

#include <xercesc/sax2/Attributes.hpp>
#include <sbml/xml/XMLNamespaces.h>

LIBSBML_CPP_NAMESPACE_BEGIN

/** @cond doxygen-libsbml-internal */

class XercesNamespaces : public XMLNamespaces
{
public:

  /**
   * Creates a new list of XML namespaces declarations from a "raw" Xerces-C++
   * Attributes set.
   */
  XercesNamespaces (const xercesc::Attributes& attrs);


  /**
   * Destroys this list of XML namespace declarations.
   */
  virtual ~XercesNamespaces ();
};

/** @endcond */

LIBSBML_CPP_NAMESPACE_END

#endif  /* __cplusplus */
#endif  /* XercesNamespaces_h */
