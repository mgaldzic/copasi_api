/**
 * @file    ExpatAttributes.h
 * @brief   Creates new XMLAttributes from "raw" Expat attributes.
 * @author  Ben Bornstein
 *
 * $Id: ExpatAttributes.h 11701 2010-08-17 01:10:32Z mhucka $
 * $HeadURL: https://sbml.svn.sourceforge.net/svnroot/sbml/trunk/libsbml/src/xml/ExpatAttributes.h $
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

#ifndef ExpatAttributes_h
#define ExpatAttributes_h

#ifdef __cplusplus

#include <string>
#include <expat.h>

#include <sbml/xml/XMLAttributes.h>

LIBSBML_CPP_NAMESPACE_BEGIN

/** @cond doxygen-libsbml-internal */

class ExpatAttributes : public XMLAttributes
{
public:

  /**
   * Creates a new XMLAttributes set from the given "raw" Expat attributes.
   * The Expat attribute names are assumed to be in namespace triplet form
   * separated by sepchar.
   *
   * @if notcpp @docnote @htmlinclude warn-default-args-in-docs.html @endif
   */
  ExpatAttributes (const XML_Char** attrs,
		   const XML_Char* elementName,
		   const XML_Char sepchar = ' ');


  /**
   * Destroys this ExpatAttributes set.
   */
  virtual ~ExpatAttributes ();
};

/** @endcond */

LIBSBML_CPP_NAMESPACE_END

#endif  /* __cplusplus */
#endif  /* ExpatAttributes_h */
