/**
 * @cond doxygen-libsbml-internal
 *
 * @file    ModelingPracticeConstraints.cpp
 * @brief   ModelingPractice check constraints.  
 * @author  Sarah Keating
 *
 * $Id: ModelingPracticeConstraints.cpp 11634 2010-08-03 03:57:18Z mhucka $
 * $HeadURL: https://sbml.svn.sourceforge.net/svnroot/sbml/trunk/libsbml/src/validator/constraints/ModelingPracticeConstraints.cpp $
 */
/* Copyright 2005 California Institute of Technology and Japan Science and
 * Technology Corporation.
 *
 * This library is free software; you can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation.  A copy of the license agreement is
 * provided in the file named "LICENSE.txt" included with this software
 * distribution.  It is also available online at
 * http://sbml.org/software/libsbml/license.html
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this library; if not, write to the Free Software Foundation,
 * Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA.
 */



#ifndef AddingConstraintsToValidator
#include <sbml/SBMLTypes.h>
#include <sbml/validator/VConstraint.h>
#endif

#include <sbml/validator/ConstraintMacros.h>

#include "LocalParameterShadowsIdInModel.h"
/** @cond doxygen-ignored */

using namespace std;

/** @endcond */

// Compartment validation

START_CONSTRAINT (80501, Compartment, c)
{
  pre( c.getLevel() > 1);
  pre( c.getSpatialDimensions() != 0 );
  
  //msg =
  //  "It is recommended that the size of a compartment is set.";

  inv( c.isSetSize() == true );
}
END_CONSTRAINT

// Parameters
EXTERN_CONSTRAINT( 81121, LocalParameterShadowsIdInModel             )


START_CONSTRAINT (80701, Parameter, p)
{
  inv(p.isSetUnits() == true);
}
END_CONSTRAINT



/** @endcond */
