/**
 * @cond doxygen-libsbml-internal
 *
 * @file    PieceBooleanMathCheck.cpp
 * @brief   Ensures piecewise piece element returns boolean.
 * @author  Sarah Keating
 *
 * $Id: PieceBooleanMathCheck.cpp 11634 2010-08-03 03:57:18Z mhucka $
 * $HeadURL: https://sbml.svn.sourceforge.net/svnroot/sbml/trunk/libsbml/src/validator/constraints/PieceBooleanMathCheck.cpp $
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

#include "PieceBooleanMathCheck.h"

/** @cond doxygen-ignored */

using namespace std;

/** @endcond */

LIBSBML_CPP_NAMESPACE_BEGIN

static const char* PREAMBLE =
    "The second argument of a MathML 'piece' operator must have a boolean "
    "value. (References: L2V2 Section 3.5.8.)";


/**
 * Creates a new Constraint with the given id.
 */
PieceBooleanMathCheck::PieceBooleanMathCheck (unsigned int id, Validator& v) : MathMLBase(id, v)
{
}


/**
 * Destroys this Constraint.
 */
PieceBooleanMathCheck::~PieceBooleanMathCheck ()
{
}


/**
 * @return the preamble to use when logging constraint violations.
 */
const char*
PieceBooleanMathCheck::getPreamble ()
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
PieceBooleanMathCheck::checkMath (const Model& m, const ASTNode& node, const SBase & sb)
{
  ASTNodeType_t type = node.getType();

  switch (type) 
  {
    case AST_FUNCTION_PIECEWISE:

      checkPiece(m, node, sb);
      break;


    case AST_FUNCTION:

      checkFunction(m, node, sb);
      break;

    default:

      checkChildren(m, node, sb);
      break;

  }
}

  
/**
  * Checks that the second argument of a piecewise returns a boolean
  *
  * If not, an error message is logged.
  */
void 
PieceBooleanMathCheck::checkPiece (const Model& m, const ASTNode& node, 
                                        const SBase & sb)
{
  unsigned int numChildren = node.getNumChildren();
  unsigned int numPieces = numChildren;

  if ((numChildren % 2) != 0) numPieces--;

  for (unsigned int n = 1; n < numPieces; n += 2)
  {
    if (!node.getChild(n)->isBoolean())
    {
      logMathConflict(node, sb);
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
PieceBooleanMathCheck::getMessage (const ASTNode& node, const SBase& object)
{

  ostringstream msg;

  //msg << getPreamble();

  char * formula = SBML_formulaToString(&node);
  msg << "The formula '" << formula;
  msg << "' in the " << getFieldname() << " element of the " << getTypename(object);
  msg << " uses an piecewise function that does not return a boolean.";
  safe_free(formula);

  return msg.str();
}

LIBSBML_CPP_NAMESPACE_END

/** @endcond */
