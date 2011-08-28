/**
 * @cond doxygen-libsbml-internal
 *
 * @file    KineticLawUnitsCheck.cpp
 * @brief   Ensures math units are consistent.
 * @author  Sarah Keating
 *
 * $Id: KineticLawUnitsCheck.cpp 11634 2010-08-03 03:57:18Z mhucka $
 * $HeadURL: https://sbml.svn.sourceforge.net/svnroot/sbml/trunk/libsbml/src/validator/constraints/KineticLawUnitsCheck.cpp $
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
#include <sbml/Compartment.h>
#include <sbml/Species.h>
#include <sbml/Parameter.h>
#include <sbml/UnitDefinition.h>
#include <sbml/Event.h>
#include <sbml/Reaction.h>
#include <sbml/EventAssignment.h>
#include <sbml/SpeciesReference.h>
#include <sbml/Rule.h>
#include <sbml/math/FormulaFormatter.h>

#include <sbml/units/UnitFormulaFormatter.h>

#include "KineticLawUnitsCheck.h"
#include "MathMLBase.h"

/** @cond doxygen-ignored */

using namespace std;

/** @endcond */

LIBSBML_CPP_NAMESPACE_BEGIN

/**
 * Creates a new Constraint with the given id.
 */
KineticLawUnitsCheck::KineticLawUnitsCheck (unsigned int id, Validator& v) :
  TConstraint<Model>(id, v)
{
}


/**
 * Destroys this Constraint.
 */
KineticLawUnitsCheck::~KineticLawUnitsCheck ()
{
}

/**
 * @return the fieldname to use logging constraint violations.  If not
 * overridden, "id" is returned.
 */
const char*
KineticLawUnitsCheck::getFieldname ()
{
  return "math";
}

/**
 * @return the preamble to use when logging constraint violations.  The
 * preamble will be prepended to each log message.  If not overriden,
 * returns an empty string.
 */
const char*
KineticLawUnitsCheck::getPreamble ()
{
  return "";
}




/**
  * Checks that the units of the result of the assignment rule
  * are consistent with variable being assigned
  *
  * @return true if units are consistent, false otherwise.
  */
void
KineticLawUnitsCheck::check_ (const Model& m, const Model& object)
{
  unsigned int n, p;
  IdList matched;
  IdList unmatched;
  const UnitDefinition *ud1, *ud2;

  if (m.getLevel() < 3)
    return;

  if (m.getNumReactions() < 2)
    return;

  /* log first KL with units declared*/ 
  for (n = 0; n < m.getNumReactions(); n++)
  {
    if (m.getReaction(n)->isSetKineticLaw())
    {
      if (m.getReaction(n)->getKineticLaw()->isSetMath())
      {
        if (!(m.getReaction(n)->getKineticLaw()->containsUndeclaredUnits()))
        {
          ud1 = m.getReaction(n)->getKineticLaw()->getDerivedUnitDefinition();
          matched.append(m.getReaction(n)->getId());
          break;
        }
      }
    }
  }

  /* loop thru remaining kl - if they are fully declared check that they match
   * and add to matched or unmatch as appropriate
   */
  for (p = n+1; p < m.getNumReactions(); p++)
  {
    if (m.getReaction(p)->isSetKineticLaw())
    {
      if (m.getReaction(p)->getKineticLaw()->isSetMath())
      {
        if (!(m.getReaction(p)->getKineticLaw()->containsUndeclaredUnits()))
        {
          ud2 = m.getReaction(p)->getKineticLaw()->getDerivedUnitDefinition();
          if (UnitDefinition::areEquivalent(ud1, ud2))
          {
            matched.append(m.getReaction(p)->getId());
          }
          else
          {
            unmatched.append(m.getReaction(p)->getId());
          }
        }
      }
    }
  }

  /* see if we have any unmatched */
  for (n = 0; n < unmatched.size(); n++)
  {

    logKLConflict(*(m.getReaction(unmatched.at(n))->getKineticLaw()->getMath()),
                  *(static_cast<const SBase *>(m.getReaction(unmatched.at(n)))));
  }

}
/**
 * @return the typename of the given SBase object.
 */
const char*
KineticLawUnitsCheck::getTypename (const SBase& object)
{
  return SBMLTypeCode_toString( object.getTypeCode() );
}


/**
 * Logs a message that the given id (and its corresponding object) have
 * failed to satisfy this constraint.
 */
void
KineticLawUnitsCheck::logKLConflict (const ASTNode& node, const SBase& object)
{
  logFailure(object, getMessage(node, object));
}

const std::string
KineticLawUnitsCheck::getMessage (const ASTNode& node, const SBase& object)
{
  ostringstream msg;

  //msg << getPreamble();
  char * formula = SBML_formulaToString(&node);
  msg << "The formula '" << formula;
  msg << "' in the KineticLaw element of the Reaction with id " << object.getId();
  msg << " produces units that are inconsistent with units of earlier KineticLaw";
  msg << " elements.";
  safe_free(formula);

  return msg.str();
}

LIBSBML_CPP_NAMESPACE_END

/** @endcond */
