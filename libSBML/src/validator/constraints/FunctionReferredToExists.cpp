/**
 * @cond doxygen-libsbml-internal
 *
 * @file    FunctionReferredToExists.cpp
 * @brief   Ensures unique variables assigned by rules and events
 * @author  Sarah Keating
 *
 * $Id: FunctionReferredToExists.cpp 11634 2010-08-03 03:57:18Z mhucka $
 * $HeadURL: https://sbml.svn.sourceforge.net/svnroot/sbml/trunk/libsbml/src/validator/constraints/FunctionReferredToExists.cpp $
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
#include <sbml/Rule.h>
#include <sbml/Event.h>
#include <sbml/EventAssignment.h>

#include "FunctionReferredToExists.h"
#include "IdList.h"

/** @cond doxygen-ignored */

using namespace std;

/** @endcond */

LIBSBML_CPP_NAMESPACE_BEGIN


/**
 * Creates a new Constraint with the given constraint id.
 */
FunctionReferredToExists::FunctionReferredToExists (unsigned int id, Validator& v) :
  TConstraint<Model>(id, v)
{
}


/**
 * Destroys this Constraint.
 */
FunctionReferredToExists::~FunctionReferredToExists ()
{
}


/**
 * Checks that all ids on the following Model objects are unique:
 * event assignments and assignment rules.
 */
void
FunctionReferredToExists::check_ (const Model& m, const Model& object)
{
  // does not apply in l2v4 and beyond
  if (m.getLevel() == 2 && m.getVersion() < 4)
  {
    unsigned int n;

    for (n = 0; n < m.getNumFunctionDefinitions(); ++n)
    {
      mFunctions.append(m.getFunctionDefinition(n)->getId());

      checkCiElements(m.getFunctionDefinition(n));
    }
  }
}

/**
  * Checks that <ci> element after an apply is already listed as a FunctionDefinition.
  */
void FunctionReferredToExists::checkCiElements(const FunctionDefinition * fd)
{
  const ASTNode* node = fd->getBody();

  checkCiIsFunction(fd, node);

  //if (node != NULL && node->getType() == AST_FUNCTION)
  //{
  //  if (!mFunctions.contains(node->getName()))
  //  {
  //    logUndefined(*fd, node->getName());
  //  }
  //}

}

/**
  * Checks that <ci> element after an apply is already listed as a FunctionDefinition.
  */
void FunctionReferredToExists::checkCiIsFunction(const FunctionDefinition * fd,
                                                 const ASTNode * node)
{
  if (!node) return;
  if (node != NULL && node->getType() == AST_FUNCTION)
  {
    if (!mFunctions.contains(node->getName()))
    {
      logUndefined(*fd, node->getName());
    }
  }

  for (unsigned int i = 0; i < node->getNumChildren(); i++)
  {
    checkCiIsFunction(fd, node->getChild(i));
  }
}

/**
  * Logs a message about an undefined <ci> element in the given
  * FunctionDefinition.
  */
void
FunctionReferredToExists::logUndefined ( const FunctionDefinition& fd,
                                       const string& varname )
{
  //msg =
  //  "Inside the 'lambda' of a <functionDefinition>, if a 'ci' element is the "
  //  "first element within a MathML 'apply', then the 'ci''s value can only "
  //  "be chosen from the set of identifiers of other SBML "
  //  "<functionDefinition>s defined prior to that point in the SBML model. In "
  //  "other words, forward references to user-defined functions are not "
  //  "permitted. (References: L2V2 Section 4.3.2.)";

  msg = "'";
  msg += varname;
  msg += "' is not listed as the id of an existing FunctionDefinition.";

  
  logFailure(fd);
}

LIBSBML_CPP_NAMESPACE_END

/** @endcond */
