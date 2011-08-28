/**
 * @cond doxygen-libsbml-internal
 *
 * @file    NumberArgsMathCheck.h
 * @brief   Ensures number of arguments to functions are appropriate.
 * @author  Sarah Keating
 *
 * $Id: NumberArgsMathCheck.h 11634 2010-08-03 03:57:18Z mhucka $
 * $HeadURL: https://sbml.svn.sourceforge.net/svnroot/sbml/trunk/libsbml/src/validator/constraints/NumberArgsMathCheck.h $
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


#ifndef NumberArgsMathCheck_h
#define NumberArgsMathCheck_h


#ifdef __cplusplus


#include <string>
#include <sstream>
#include <math.h>

#include <sbml/validator/VConstraint.h>

#include "MathMLBase.h"

LIBSBML_CPP_NAMESPACE_BEGIN

class ASTNode;


class NumberArgsMathCheck: public MathMLBase
{
public:

  /**
   * Creates a new Constraint with the given id.
   */
  NumberArgsMathCheck (unsigned int id, Validator& v);

  /**
   * Destroys this Constraint.
   */
  virtual ~NumberArgsMathCheck ();


protected:

  /**
   * Checks the MathML of the ASTnode 
   * is appropriate for the function being performed
   *
   * If an inconsistency is found, an error message is logged.
   */
  virtual void checkMath (const Model& m, const ASTNode& node, const SBase & sb);
 
  /**
   * Checks that the function has only one argument
   */
  void checkUnary(const Model& m, const ASTNode& node, const SBase & sb);

  /**
   * Checks that the function has exactly two arguments
   */
  void checkBinary(const Model& m, const ASTNode& node, const SBase & sb);

  /**
   * Checks that the function at least two arguments
   */
  void checkNary(const Model& m, const ASTNode& node, const SBase & sb);

  /**
   * Checks that the functions have either one or two arguments
   */
  void checkSpecialCases(const Model& m, const ASTNode& node, const SBase & sb);

  /**
   * @return the preamble to use when logging constraint violations.  The
   * preamble will be prepended to each log message.  If not overriden,
   * returns an empty string.
   */
  virtual const char* getPreamble ();

  /**
   * Checks that the arguments to logical operators are all boolean
   *
   * If not, an error message is logged.
   */
  void checkMathFromLogical (const Model& m, const ASTNode& node, const SBase & sb);

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
#endif  /* NumberArgsMathCheck_h */

/** @endcond */
