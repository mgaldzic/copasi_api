/**
 * @cond doxygen-libsbml-internal
 *
 * @file    AssignmentRuleOrdering.cpp
 * @brief   Checks rule ordering for l2v1 and l1
 * @author  Sarah Keating
 *
 * $Id: AssignmentRuleOrdering.cpp 11634 2010-08-03 03:57:18Z mhucka $
 * $HeadURL: https://sbml.svn.sourceforge.net/svnroot/sbml/trunk/libsbml/src/validator/constraints/AssignmentRuleOrdering.cpp $
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

#include <cstring>

#include <sbml/Model.h>
#include <sbml/Rule.h>
#include <sbml/Reaction.h>
#include <sbml/InitialAssignment.h>
#include <sbml/util/List.h>
#include <sbml/util/memory.h>

#include "AssignmentRuleOrdering.h"
#include "IdList.h"

/** @cond doxygen-ignored */

using namespace std;

/** @endcond */

LIBSBML_CPP_NAMESPACE_BEGIN

/**
 * Creates a new Constraint with the given constraint id.
 */
AssignmentRuleOrdering::AssignmentRuleOrdering (unsigned int id, Validator& v) :
  TConstraint<Model>(id, v)
{
}


/**
 * Destroys this Constraint.
 */
AssignmentRuleOrdering::~AssignmentRuleOrdering ()
{
}


/**
 * Checks that all ids on the following Model objects are unique:
 * event assignments and assignment rules.
 */
void
AssignmentRuleOrdering::check_ (const Model& m, const Model& object)
{
  //// this rule ony applies in l2v1 and l1
  //if (!(object.getLevel() == 1 
  //  || (object.getLevel() == 2 && object.getVersion() == 1)))
  //  return;

  unsigned int n;

  mVariableList.clear();

  // create a list of all assignment rule variables 
  // in the order they appear
  for (n = 0; n < m.getNumRules(); ++n)
  { 
    if (m.getRule(n)->isAssignment())
    {
      mVariableList.append(m.getRule(n)->getId());
    }
  }
 
  for (n = 0; n < m.getNumRules(); ++n)
  { 
    if (m.getRule(n)->isAssignment())
    {
      if (m.getRule(n)->isSetMath())
      {
        checkRuleForVariable(m, *m.getRule(n));
        checkRuleForLaterVariables(m, *m.getRule(n), n);
      }
    }
  }
}
 
void 
AssignmentRuleOrdering::checkRuleForVariable(const Model& m, const Rule& object)
{
  /* list the <ci> elements */
  List* variables = object.getMath()->getListOfNodes( ASTNode_isName );
  std::string variable = object.getVariable();

  if (variables)
  {
    for (unsigned int i = 0; i < variables->getSize(); i++)
    {
      ASTNode* node = static_cast<ASTNode*>( variables->get(i) );
      const char *   name = node->getName() ? node->getName() : "";
      if (!(strcmp(variable.c_str(), name)))
        logRuleRefersToSelf(*(object.getMath()), object);
    }
    // return value of ASTNode::getListOfNodes() needs to be
    // deleted by caller.
    delete variables;
  }

}


void 
AssignmentRuleOrdering::checkRuleForLaterVariables(const Model& m, 
                                                   const Rule& object,
                                                   unsigned int n)
{
  /* list the <ci> elements of this rule*/
  List* variables = object.getMath()->getListOfNodes( ASTNode_isName );

  if (variables)
  {
    unsigned int index;
    for (unsigned int i = 0; i < variables->getSize(); i++)
    {
      ASTNode* node = static_cast<ASTNode*>( variables->get(i) );
      const char *   name = node->getName() ? node->getName() : "";
  
      if (mVariableList.contains(name))
      {
        // this <ci> is a variable
        // check that it occurs later
        index = 0; 
        while(index < mVariableList.size())
        {
          if (!strcmp(name, mVariableList.at(index).c_str()))
            break;
          index++;
        }
        if (index > n)
          logForwardReference(*(object.getMath()), object, name);
      }
    }
    // return value of ASTNode::getListOfNodes() needs to be
    // deleted by caller.
    delete variables;
  }
}


void
AssignmentRuleOrdering::logRuleRefersToSelf (const ASTNode & node,
                                             const SBase& object)
{
  char * formula = SBML_formulaToString(&node);
  msg =
    "The AssignmentRule with variable '";
  msg += object.getId();
  msg += "' refers to that variable within the math formula '";
  msg += formula;
  msg += "'.";
  safe_free(formula);
  
  logFailure(object);

}

void
AssignmentRuleOrdering::logForwardReference (const ASTNode & node,
                                             const SBase& object,
                                             std::string name)
{
  char * formula = SBML_formulaToString(&node);
  msg =
    "The AssignmentRule with variable '";
  msg += object.getId();
  msg += "' refers to the variable '";
  msg += name;
  msg += "' within the math formula '";
  msg += formula;
  msg += "'. '";
  msg += name;
  msg += "' is the subject of a later assignment rule.";
  safe_free(formula);
  
  logFailure(object);

}

LIBSBML_CPP_NAMESPACE_END

/** @endcond */
