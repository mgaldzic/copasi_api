/**
 * @cond doxygen-libsbml-internal
 *
 * @file    L1CompatibilityValidator.cpp
 * @brief   Checks whether an SBML model can be converted from L2 to L1
 * @author  Sarah Keating
 *
 * $Id: L1CompatibilityValidator.cpp 11634 2010-08-03 03:57:18Z mhucka $
 * $HeadURL: https://sbml.svn.sourceforge.net/svnroot/sbml/trunk/libsbml/src/validator/L1CompatibilityValidator.cpp $
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

#ifndef doxygen_ignore

#include <sbml/SBMLTypes.h>
#include <sbml/validator/L1CompatibilityValidator.h>


/*
 * Compile L1CompatibilityConstraints
 */
#include "constraints/L1CompatibilityConstraints.cpp"

LIBSBML_CPP_NAMESPACE_BEGIN

/*
 * Initializes this Validator with a set of Constraints.
 */
void
L1CompatibilityValidator::init ()
{
#define  AddingConstraintsToValidator 1
#include "constraints/L1CompatibilityConstraints.cpp"
}

LIBSBML_CPP_NAMESPACE_END

#endif

/** @endcond */
