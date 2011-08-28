/**
 * @file    LibXMLAttributes.h
 * @brief   Creates new XMLAttributes from "raw" LibXML attributes.
 * @author  Ben Bornstein
 *
 * $Id: LibXMLAttributes.h 11633 2010-08-03 03:53:27Z mhucka $
 * $HeadURL: https://sbml.svn.sourceforge.net/svnroot/sbml/trunk/libsbml/src/xml/LibXMLAttributes.h $
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

#ifndef LibXMLAttributes_h
#define LibXMLAttributes_h

#ifdef __cplusplus

#include <string>

#include <libxml/parser.h>
#include <sbml/xml/XMLAttributes.h>

LIBSBML_CPP_NAMESPACE_BEGIN

/** @cond doxygen-libsbml-internal */

class LibXMLAttributes : public XMLAttributes
{
public:

  /**
   * Creates a new XMLAttributes set from the given "raw" LibXML attributes.
   */
  LibXMLAttributes (  const xmlChar** attributes
		    , const xmlChar*  elementName
		    , const unsigned int& size);


  /**
   * Destroys this LibXMLAttributes set.
   */
  virtual ~LibXMLAttributes ();
};

/** @endcond */

LIBSBML_CPP_NAMESPACE_END

#endif  /* __cplusplus */
#endif  /* LibXMLAttributes_h */
