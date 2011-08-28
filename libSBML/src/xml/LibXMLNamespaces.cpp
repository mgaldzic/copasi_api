/**
 * @file    LibXMLNamespaces.cpp
 * @brief   Extracts XML namespace declarations from LibXML prefix/URI pairs.
 * @author  Ben Bornstein
 *
 * $Id: LibXMLNamespaces.cpp 11633 2010-08-03 03:53:27Z mhucka $
 * $HeadURL: https://sbml.svn.sourceforge.net/svnroot/sbml/trunk/libsbml/src/xml/LibXMLNamespaces.cpp $
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

#include <sbml/xml/LibXMLTranscode.h>
#include <sbml/xml/LibXMLNamespaces.h>

/** @cond doxygen-ignored */

using namespace std;

/** @endcond */

LIBSBML_CPP_NAMESPACE_BEGIN

/** @cond doxygen-libsbml-internal */

/**
 * Creates a new list of XML namespaces declarations from a "raw" LibXML
 * prefix/URI pairs.
 */
LibXMLNamespaces::LibXMLNamespaces (  const xmlChar**     namespaces
                                    , const unsigned int& size )
{
  mNamespaces.reserve(size);

  for (unsigned int n = 0; n < size; ++n)
  {
    const string prefix = LibXMLTranscode( namespaces[2 * n]     );
    const string uri    = LibXMLTranscode( namespaces[2 * n + 1], true );

    add(uri, prefix);
  }
}


/**
 * Destroys this list of XML namespace declarations.
 */
LibXMLNamespaces::~LibXMLNamespaces ()
{
}

/** @endcond */

LIBSBML_CPP_NAMESPACE_END
