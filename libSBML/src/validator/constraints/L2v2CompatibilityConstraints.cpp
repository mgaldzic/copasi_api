/**
 * @cond doxygen-libsbml-internal
 *
 * @file    L2v2CompatibilityConstraints.cpp
 * @brief   L1 compatibility for conversion from L2
 * @author  Sarah Keating
 *
 * $Id: L2v2CompatibilityConstraints.cpp 11634 2010-08-03 03:57:18Z mhucka $
 * $HeadURL: https://sbml.svn.sourceforge.net/svnroot/sbml/trunk/libsbml/src/validator/constraints/L2v2CompatibilityConstraints.cpp $
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
#include <sbml/SBase.h>
#include <sbml/validator/VConstraint.h>
#include <math.h>
#include "DuplicateTopLevelAnnotation.h"
#endif


#include <sbml/validator/ConstraintMacros.h>

/** @cond doxygen-ignored */

using namespace std;

/** @endcond */


START_CONSTRAINT (93001, UnitDefinition, ud)
{
  //msg =
  //  "In SBML Level 2 Version 2, an 'sboTerm' attribute is only permitted on "
  //  "the following elements: <model>, <functionDefinition>, <parameter>, "
  //  "<initialAssignment>, <rule>, <constraint>, <reaction>, "
  //  "<speciesReference>, <kineticLaw>, <event> and <eventAssignment>.";

  inv( ud.getSBOTerm() == -1 );
}
END_CONSTRAINT


START_CONSTRAINT (93001, Unit, u)
{
  //msg =
  //  "In SBML Level 2 Version 2, an 'sboTerm' attribute is only permitted on "
  //  "the following elements: <model>, <functionDefinition>, <parameter>, "
  //  "<initialAssignment>, <rule>, <constraint>, <reaction>, "
  //  "<speciesReference>, <kineticLaw>, <event> and <eventAssignment>.";

  inv( u.getSBOTerm() == -1 );
}
END_CONSTRAINT


START_CONSTRAINT (93001, CompartmentType, ct)
{
  //msg =
  //  "In SBML Level 2 Version 2, an 'sboTerm' attribute is only permitted on "
  //  "the following elements: <model>, <functionDefinition>, <parameter>, "
  //  "<initialAssignment>, <rule>, <constraint>, <reaction>, "
  //  "<speciesReference>, <kineticLaw>, <event> and <eventAssignment>.";

  inv( ct.getSBOTerm() == -1 );
}
END_CONSTRAINT


START_CONSTRAINT (93001, SpeciesType, st)
{
  //msg =
  //  "In SBML Level 2 Version 2, an 'sboTerm' attribute is only permitted on "
  //  "the following elements: <model>, <functionDefinition>, <parameter>, "
  //  "<initialAssignment>, <rule>, <constraint>, <reaction>, "
  //  "<speciesReference>, <kineticLaw>, <event> and <eventAssignment>.";

  inv( st.getSBOTerm() == -1 );
}
END_CONSTRAINT


START_CONSTRAINT (93001, Compartment, c)
{
  //msg =
  //  "In SBML Level 2 Version 2, an 'sboTerm' attribute is only permitted on "
  //  "the following elements: <model>, <functionDefinition>, <parameter>, "
  //  "<initialAssignment>, <rule>, <constraint>, <reaction>, "
  //  "<speciesReference>, <kineticLaw>, <event> and <eventAssignment>.";

  inv( c.getSBOTerm() == -1 );
}
END_CONSTRAINT


START_CONSTRAINT (93001, Species, s)
{
  //msg =
  //  "In SBML Level 2 Version 2, an 'sboTerm' attribute is only permitted on "
  //  "the following elements: <model>, <functionDefinition>, <parameter>, "
  //  "<initialAssignment>, <rule>, <constraint>, <reaction>, "
  //  "<speciesReference>, <kineticLaw>, <event> and <eventAssignment>.";

  inv( s.getSBOTerm() == -1 );
}
END_CONSTRAINT


START_CONSTRAINT (93001, Trigger, t)
{
  //msg =
  //  "In SBML Level 2 Version 2, an 'sboTerm' attribute is only permitted on "
  //  "the following elements: <model>, <functionDefinition>, <parameter>, "
  //  "<initialAssignment>, <rule>, <constraint>, <reaction>, "
  //  "<speciesReference>, <kineticLaw>, <event> and <eventAssignment>.";

  inv( t.getSBOTerm() == -1 );
}
END_CONSTRAINT


START_CONSTRAINT (93001, Delay, d)
{
  //msg =
  //  "In SBML Level 2 Version 2, an 'sboTerm' attribute is only permitted on "
  //  "the following elements: <model>, <functionDefinition>, <parameter>, "
  //  "<initialAssignment>, <rule>, <constraint>, <reaction>, "
  //  "<speciesReference>, <kineticLaw>, <event> and <eventAssignment>.";

  inv( d.getSBOTerm() == -1 );
}
END_CONSTRAINT


START_CONSTRAINT (93002, Unit, u)
{
  //msg =
  //  "The 'offset' attribute on <unit> previously available in SBML Level 2 "
  //  "Version 1, has been removed as of SBML Level 2 Version 2. (References: "
  //  "L2V2 Section 4.4.)";

  inv( u.getOffset() == 0.0 );
}
END_CONSTRAINT

START_CONSTRAINT (93003, KineticLaw, kl)
{
  //msg =
  //  "The 'timeUnits' attribute on <kineticLaw>, previously available in SBML "
  //  "Level 1 and Level 2 Version 1, has been removed as of SBML Level 2 "
  //  "Version 2. In SBML Level 2 Version 2, the time units of a reaction rate "
  //  "expression are those of the global 'time' units of the model. "
  //  "(References: L2V2 Section 4.13.5.)";

  inv( kl.isSetTimeUnits() == false );
}
END_CONSTRAINT

START_CONSTRAINT (93004, KineticLaw, kl)
{
  //msg =
  //  "The 'substanceUnits' attribute on <kineticLaw>, previously available in "
  //  "SBML Level 1 and Level 2 Version 1, has been removed as of SBML Level 2 "
  //  "Version 2. In SBML Level 2 Version 2, the substance units of a reaction "
  //  "rate expression are those of the global 'substance' units of the model. "
  //  "(References: L2V2 Section 4.13.5.)";

  inv( kl.isSetSubstanceUnits() == false );
}
END_CONSTRAINT


START_CONSTRAINT (93005, Event, e)
{
  inv( e.getUseValuesFromTriggerTime() == true);
}
END_CONSTRAINT

START_CONSTRAINT (93006, Model, m1)
{
  // if the model was L2V4 or above the model sbo term will not
  // be valid in l2v2
  pre( m1.getLevel() >1 );
  if (m1.getLevel() == 2)
  {
    pre( m1.getVersion() > 3);
  }

  inv( !m1.isSetSBOTerm());
}
END_CONSTRAINT

EXTERN_CONSTRAINT(93009, DuplicateTopLevelAnnotation)


/** @endcond */
