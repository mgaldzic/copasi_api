/**
 * @cond doxygen-libsbml-internal
 *
 * @file    StoichiometryMathVars.cpp
 * @brief   Ensures FunctionDefinitions contain no undefined variables.
 * @author  Ben Bornstein
 *
 * $Id: StoichiometryMathVars.cpp 11634 2010-08-03 03:57:18Z mhucka $
 * $HeadURL: https://sbml.svn.sourceforge.net/svnroot/sbml/trunk/libsbml/src/validator/constraints/StoichiometryMathVars.cpp $
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
#include <sbml/Reaction.h>
#include <sbml/util/List.h>
#include <sbml/math/ASTNode.h>

#include "StoichiometryMathVars.h"

/** @cond doxygen-ignored */

using namespace std;

/** @endcond */

LIBSBML_CPP_NAMESPACE_BEGIN

/**
 * Creates a new Constraint with the given id.
 */
StoichiometryMathVars::StoichiometryMathVars (unsigned int id, Validator& v) :
  TConstraint<Reaction>(id, v)
{
}


/**
 * Destroys this Constraint.
 */
StoichiometryMathVars::~StoichiometryMathVars ()
{
}


/**
 * Checks that all variables referenced in FunctionDefinition bodies are
 * bound variables (function arguments).
 */
void
StoichiometryMathVars::check_ (const Model& m, const Reaction& r)
{
  if (r.getLevel() == 1) return;

  unsigned int n, ns;
  
  for (n = 0; n < r.getNumReactants(); n++)
  {
    mSpecies.append(r.getReactant(n)->getSpecies());
  }
  for (n = 0; n < r.getNumProducts(); n++)
  {
    mSpecies.append(r.getProduct(n)->getSpecies());
  }
  for (n = 0; n < r.getNumModifiers(); n++)
  {
    mSpecies.append(r.getModifier(n)->getSpecies());
  }

  for (n = 0; n < r.getNumReactants(); n++)
  {
    const SpeciesReference *sr = r.getReactant(n);

    if (sr->isSetStoichiometryMath())
    {
      List* variables = sr->getStoichiometryMath()->getMath()
        ->getListOfNodes( ASTNode_isName );
      for (ns = 0; ns < variables->getSize(); ns++)
      {
        ASTNode* node = static_cast<ASTNode*>( variables->get(ns) );
        string   name = node->getName() ? node->getName() : "";
       
        if (m.getSpecies(name) && !mSpecies.contains(name))
          logUndefined(r, name);
      }

      delete variables;
      
    }
  }
  for (n = 0; n < r.getNumProducts(); n++)
  {
    const SpeciesReference *sr = r.getProduct(n);

    if (sr->isSetStoichiometryMath())
    {
      List* variables = sr->getStoichiometryMath()->getMath()
        ->getListOfNodes( ASTNode_isName );
      for (ns = 0; ns < variables->getSize(); ns++)
      {
        ASTNode* node = static_cast<ASTNode*>( variables->get(ns) );
        string   name = node->getName() ? node->getName() : "";
       
        if (m.getSpecies(name) && !mSpecies.contains(name))
          logUndefined(r, name);
      }

      delete variables;

    }
  }
}


/**
 * Logs a message about an undefined variable in the given
 * FunctionDefinition.
 */
void
StoichiometryMathVars::logUndefined ( const Reaction& r,
                                       const string& varname )
{
  msg =
    //"All species referenced in the <stoichiometryMath> formula of a given "
    //"reaction must first be declared using <speciesReference> or "
    //"<modifierSpeciesReference>. More formally, if a <species> identifier "
    //"appears in a <ci> element of a <reaction>'s <stoichiometryMath> "
    //"formula, that same identifier must also appear in at least one "
    //"<speciesReference> or <modifierSpeciesReference> in the <reaction> "
    //"definition. (References: L2V2 Sections 4.13.2 and 4.13.4; L2V3 Sections "
    //"4.13.2 and 4.13.4.)"
    "The species '";

  msg += varname;
  msg += "' is not listed as a product, reactant, or modifier of reaction '";
  msg += r.getId();
  msg += "'.";
  
  logFailure(r);
}

LIBSBML_CPP_NAMESPACE_END

/** @endcond */
