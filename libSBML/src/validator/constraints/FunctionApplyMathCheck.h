/**
 * @cond doxygen-libsbml-internal
 *
 * @file    FunctionApplyMathCheck.h
 * @brief   Ensures <ci> after apply refers to function definition.
 * @author  Sarah Keating
 *
 * $Id: FunctionApplyMathCheck.h 11634 2010-08-03 03:57:18Z mhucka $
 * $HeadURL: https://sbml.svn.sourceforge.net/svnroot/sbml/trunk/libsbml/src/validator/constraints/FunctionApplyMathCheck.h $
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


#ifndef FunctionApplyMathCheck_h
#define FunctionApplyMathCheck_h


#ifdef __cplusplus


#include <string>
#include <sstream>
#include <math.h>

#include <sbml/validator/VConstraint.h>

#include "MathMLBase.h"

LIBSBML_CPP_NAMESPACE_BEGIN

class ASTNode;


class FunctionApplyMathCheck: public MathMLBase
{
public:

  /**
   * Creates a new Constraint with the given id.
   */
  FunctionApplyMathCheck (unsigned int id, Validator& v);

  /**
   * Destroys this Constraint.
   */
  virtual ~FunctionApplyMathCheck ();


protected:

  /**
   * Checks the MathML of the ASTnode 
   * is appropriate for the function being performed
   *
   * If an inconsistency is found, an error message is logged.
   */
  virtual void checkMath (const Model& m, const ASTNode& node, const SBase & sb);
  
  /**
   * @return the preamble to use when logging constraint violations.  The
   * preamble will be prepended to each log message.  If not overriden,
   * returns an empty string.
   */
  virtual const char* getPreamble ();

  /**
   * Checks that the functionDefinition referred to by a <ci> element exists
   *
   * If <ci> does not refer to functionDefinition id, an error message is logged.
   */
  void checkExists (const Model& m, const ASTNode& node, const SBase & sb);

  /**
   * @return the error message to use when logging constraint violations.
   * This method is called by logFailure.
   *
   * If at all possible please use getPreamble() and getFieldname() when
   * constructing error messages.  This will help to make your constraint
   * easily customizable.
   */
  virtual const std::string
  getMessage (const ASTNode& node, const SBase& object);

};

LIBSBML_CPP_NAMESPACE_END

#endif  /* __cplusplus */
#endif  /* FunctionApplyMathCheck_h */

/** @endcond */
