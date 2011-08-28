/**
 * @file    ExpatAttributes.cpp
 * @brief   Creates new XMLAttributes from "raw" Expat attributes.
 * @author  Ben Bornstein
 *
 * $Id: ExpatAttributes.cpp 11633 2010-08-03 03:53:27Z mhucka $
 * $HeadURL: https://sbml.svn.sourceforge.net/svnroot/sbml/trunk/libsbml/src/xml/ExpatAttributes.cpp $
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

#include <sbml/xml/ExpatAttributes.h>

/** @cond doxygen-ignored */

using namespace std;

LIBSBML_CPP_NAMESPACE_BEGIN

/** @endcond */


/** @cond doxygen-libsbml-internal */

/**
 * Creates a new XMLAttributes set from the given "raw" Expat attributes.
 * The Expat attribute names are assumed to be in namespace triplet form
 * separated by sepchar.
 */
ExpatAttributes::ExpatAttributes (const XML_Char** attrs,
				  const XML_Char* elementName,
				  const XML_Char sep)
{
  unsigned int size = 0;
  while (attrs[2 * size]) ++size;

  mNames .reserve(size);
  mValues.reserve(size);

  for (unsigned int n = 0; n < size; ++n)
  {
    mNames .push_back( XMLTriple( attrs[2 * n], sep ) );
    mValues.push_back( string   ( attrs[2 * n + 1]  ) );
  }

  mElementName = elementName;
}


/**
 * Destroys this Attribute set.
 */
ExpatAttributes::~ExpatAttributes ()
{
}

LIBSBML_CPP_NAMESPACE_END

/** @endcond */
