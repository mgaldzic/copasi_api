/**
 * @cond doxygen-libsbml-internal
 *
 * @file    CiElementMathCheck.cpp
 * @brief   checks <ci> element is the id of a component
 * @author  Sarah Keating
 *
 * $Id: CiElementMathCheck.cpp 11634 2010-08-03 03:57:18Z mhucka $
 * $HeadURL: https://sbml.svn.sourceforge.net/svnroot/sbml/trunk/libsbml/src/validator/constraints/CiElementMathCheck.cpp $
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
#include <sbml/SBMLTypeCodes.h>

#include <sbml/units/UnitFormulaFormatter.h>

#include "CiElementMathCheck.h"

/** @cond doxygen-ignored */

using namespace std;

/** @endcond */

LIBSBML_CPP_NAMESPACE_BEGIN

static const char* PREAMBLE =
    "Outside of a <functionDefinition>, if a 'ci' element is not the first "
    "element within a MathML 'apply', then the 'ci''s value can only be "
    "chosen from the set of identifiers of <species>, <compartment>, "
    "<parameter> or <reaction> objects defined in the SBML model. "
    "(References: L2V2 Section 3.5.3.)";

/**
 * Creates a new Constraint with the given id.
 */
CiElementMathCheck::CiElementMathCheck (unsigned int id, Validator& v) : MathMLBase(id, v)
{
}


/**
 * Destroys this Constraint.
 */
CiElementMathCheck::~CiElementMathCheck ()
{
}


/**
 * @return the preamble to use when logging constraint violations.
 */
const char*
CiElementMathCheck::getPreamble ()
{
  return PREAMBLE;
}


/**
  * Checks the MathML of the ASTnode 
  * is appropriate for the function being performed
  *
  * If an inconsistency is found, an error message is logged.
  */
void
CiElementMathCheck::checkMath (const Model& m, const ASTNode& node, const SBase & sb)
{
  ASTNodeType_t type = node.getType();
    
  /* if the node is a <ci> element it will have type AST_NAME
   * check that this name is an appropriate component of the model */
  switch (type) 
  {
    case AST_NAME:

      checkCiElement(m, node, sb);
      break;

    default:

      checkChildren(m, node, sb);
      break;

  }
}

  
/**
  * Checks any <ci> elements in the MathML of the ASTnode 
  * contain the id of an appropriate component of the model
  *
  * If an inconsistency is found, an error message is logged.
  */
void 
CiElementMathCheck::checkCiElement (const Model& m, 
                                        const ASTNode& node, 
                                        const SBase & sb)
{
  std::string name = node.getName();
  const KineticLaw * kl;

  /* leave out this check if the ci element is a local parameter in a kineticLaw
   * caught by LocalParameterMathCheck instead
   */
  if (!mLocalParameters.contains(name))
  {
    bool allowReactionId = true;
    bool allowSpeciesRef = false;

    if ( (m.getLevel() == 2) && (m.getVersion() == 1) )
      allowReactionId = false;

    if (m.getLevel() > 2)
      allowSpeciesRef = true;

    if (!m.getCompartment(name) &&
        !m.getSpecies(name)     &&
        !m.getParameter(name)   &&
        (!allowReactionId || !m.getReaction(name)) &&
        (!allowSpeciesRef || !m.getSpeciesReference(name)) )
    {
      /* check whether we are in a kinetic law since there
      * may be local parameters
      */

      if (sb.getTypeCode() == SBML_KINETIC_LAW)
      {
        kl = m.getReaction(mKLCount)->getKineticLaw();

        if (!kl->getParameter(name))
        {
          logMathConflict(node, sb);
        }
      }
      else
      {
          logMathConflict(node, sb);
      }
    }
  }
    
}


/**
 * @return the error message to use when logging constraint violations.
 * This method is called by logFailure.
 *
 * Returns a message that the given id and its corresponding object are
 * in  conflict with an object previously defined.
 */
const string
CiElementMathCheck::getMessage (const ASTNode& node, const SBase& object)
{

  ostringstream msg;

  //msg << getPreamble();
  char * formula = SBML_formulaToString(&node);
  msg << "\nThe formula '" << formula;
  msg << "' in the " << getFieldname() << " element of the " << getTypename(object);
  if (object.getLevel() == 2 && object.getVersion() == 1)
    msg << " uses '" << node.getName() << "' that is not the id of a species/compartment/parameter.";
  else
    msg << " uses '" << node.getName() << "' that is not the id of a species/compartment/parameter/reaction.";
  safe_free(formula);

  return msg.str();
}

LIBSBML_CPP_NAMESPACE_END

/** @endcond */
