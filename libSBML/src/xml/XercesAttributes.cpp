/**
 * @file    XercesAttributes.cpp
 * @brief   Creates new XMLAttributes from "raw" Xerces-C++ attributes.
 * @author  Ben Bornstein
 *
 * $Id: XercesAttributes.cpp 11633 2010-08-03 03:53:27Z mhucka $
 * $HeadURL: https://sbml.svn.sourceforge.net/svnroot/sbml/trunk/libsbml/src/xml/XercesAttributes.cpp $
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

#include <sbml/xml/XercesTranscode.h>
#include <sbml/xml/XercesAttributes.h>

/** @cond doxygen-ignored */

using namespace std;
using namespace xercesc;

LIBSBML_CPP_NAMESPACE_BEGIN

/** @endcond */


/** @cond doxygen-libsbml-internal */

/**
 * Creates a new XMLAttributes set that wraps the given "raw" Xerces-C++
 * Attributes set.
 */
XercesAttributes::XercesAttributes (const Attributes& attrs,
				    const string elementName)
{
  unsigned int size = attrs.getLength();

  mNames .reserve(size);
  mValues.reserve(size);

  for (unsigned int n = 0; n < size; ++n)
  {
    const string            uri   = XercesTranscode( attrs.getURI      (n) );
    const string            name  = XercesTranscode( attrs.getLocalName(n) );
    const string            qname = XercesTranscode( attrs.getQName    (n) );
    const string            value = XercesTranscode( attrs.getValue    (n) );
    const string::size_type pos   = qname.find(":", 0);

    const string prefix = (pos != string::npos) ? qname.substr(0, pos) : "";

    //
    // Skip XML namespace declarations
    //
    if (prefix != "xmlns" && name != "xmlns")
    {
      mNames .push_back( XMLTriple(name, uri, prefix) );
      mValues.push_back( value );
    }
  }

  mElementName = elementName;
}


/**
 * Destroys this XercesAttributes set.
 */
XercesAttributes::~XercesAttributes ()
{
}

LIBSBML_CPP_NAMESPACE_END

/** @endcond */
