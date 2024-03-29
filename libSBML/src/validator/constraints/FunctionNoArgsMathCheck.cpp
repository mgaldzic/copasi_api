/**
 * @cond doxygen-libsbml-internal
 *
 * @file    FunctionNoArgsMathCheck.cpp
 * @brief   Ensures correct number of arguments to a function definition.
 * @author  Sarah Keating
 *
 * $Id: FunctionNoArgsMathCheck.cpp 11634 2010-08-03 03:57:18Z mhucka $
 * $HeadURL: https://sbml.svn.sourceforge.net/svnroot/sbml/trunk/libsbml/src/validator/constraints/FunctionNoArgsMathCheck.cpp $
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

#include "FunctionNoArgsMathCheck.h"

/** @cond doxygen-ignored */

using namespace std;

/** @endcond */

LIBSBML_CPP_NAMESPACE_BEGIN

static const char* PREAMBLE =
    "The number of arguments used in a call to a function defined by a "
    "<functionDefinition> must equal the number of arguments accepted by "
    "that function, or in other words, the number of <bvar> elements "
    "inside the <lambda> element of the function definition.  "
    "(References: L2V4 Section 4.3.4.)";


/**
 * Creates a new Constraint with the given id.
 */
FunctionNoArgsMathCheck::FunctionNoArgsMathCheck (unsigned int id, Validator& v) : MathMLBase(id, v)
{
}


/**
 * Destroys this Constraint.
 */
FunctionNoArgsMathCheck::~FunctionNoArgsMathCheck ()
{
}


/**
 * @return the preamble to use when logging constraint violations.
 */
const char*
FunctionNoArgsMathCheck::getPreamble ()
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
FunctionNoArgsMathCheck::checkMath (const Model& m, const ASTNode& node, const SBase & sb)
{
  ASTNodeType_t type = node.getType();

  switch (type) 
  {
    case AST_FUNCTION:

      checkNumArgs(m, node, sb);
      break;

    default:

      checkChildren(m, node, sb);
      break;

  }
}

  
/**
  * Checks that the functionDefinition referred to by a <ci> element 
  * has the appropriate number of arguments.
  *
  * If not, an error message is logged.
  */
void 
FunctionNoArgsMathCheck::checkNumArgs (const Model& m, const ASTNode& node, 
                                                const SBase & sb)
{
  /* this rule was only introduced level 2 version 4 */
  if (m.getLevel() > 2 || (m.getLevel() == 2 && m.getVersion() > 3))
  {
    if (m.getFunctionDefinition(node.getName()))
    {
      /* functiondefinition math */
      const ASTNode * fdMath = m.getFunctionDefinition(node.getName())->getMath();
      if (fdMath != NULL)
      {
      /* We have a definition for this function.  Does the defined number
	        of arguments equal the number used here? */

        if (node.getNumChildren() + 1 != fdMath->getNumChildren())
	      {
          logMathConflict(node, sb);
	      }
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
FunctionNoArgsMathCheck::getMessage (const ASTNode& node, const SBase& object)
{

  ostringstream msg;

  //msg << getPreamble();
  char * formula = SBML_formulaToString(&node);
  msg << "\nThe formula '" << formula;
  msg << "' in the " << getFieldname() << " element of the ";
  msg << getTypename(object);
  msg << " uses the function '" << node.getName() << "' which requires ";
  msg << "a different number of arguments than the number supplied.";
  safe_free(formula);

  return msg.str();
}

LIBSBML_CPP_NAMESPACE_END

/** @endcond */
