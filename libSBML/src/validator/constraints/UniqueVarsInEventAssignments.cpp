/**
 * @cond doxygen-libsbml-internal
 *
 * @file    UniqueVarsInEventAssignments.cpp
 * @brief   Ensures variables within a set of EventAssignments are unique
 * @author  Ben Bornstein
 *
 * $Id: UniqueVarsInEventAssignments.cpp 11634 2010-08-03 03:57:18Z mhucka $
 * $HeadURL: https://sbml.svn.sourceforge.net/svnroot/sbml/trunk/libsbml/src/validator/constraints/UniqueVarsInEventAssignments.cpp $
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


#include <sbml/Model.h>
#include <sbml/Event.h>
#include <sbml/EventAssignment.h>

#include "UniqueVarsInEventAssignments.h"

/** @cond doxygen-ignored */

using namespace std;

/** @endcond */

LIBSBML_CPP_NAMESPACE_BEGIN

static const char* PREAMBLE =
    "In each <event>, the value of the 'variable' field within every "
    "<eventAssignment> definition must be unique across the set of all "
    "<eventAssignment>s within that <event>. (References: L2V1 erratum 17; "
    "L2V2 Section 4.14.)";


/**
 * Creates a new Constraint with the given constraint id.
 */
UniqueVarsInEventAssignments::UniqueVarsInEventAssignments ( unsigned int id,
                                                             Validator& v ) :
  UniqueIdBase(id, v)
{
}


/**
 * Destroys this Constraint.
 */
UniqueVarsInEventAssignments::~UniqueVarsInEventAssignments ()
{
}


/**
 * @return the fieldname ("variable") to use when logging constraint
 * violations.
 */
const char*
UniqueVarsInEventAssignments::getFieldname ()
{
  return "variable";
}


/**
 * @return the preamble to use when logging constraint violations.
 */
const char*
UniqueVarsInEventAssignments::getPreamble ()
{
  return PREAMBLE;
}


/**
 * Checks that all ids on KineticLawParameters are unique.
 */
void
UniqueVarsInEventAssignments::doCheck (const Model& m)
{
  for (unsigned int n = 0; n < m.getNumEvents(); ++n)
  {
    const Event* e = m.getEvent(n);

    for (unsigned int ea = 0; ea < e->getNumEventAssignments(); ++ea)
    {
      checkId( *e->getEventAssignment(ea) );
    }

    mIdObjectMap.clear();
  }
}

LIBSBML_CPP_NAMESPACE_END

/** @endcond */
