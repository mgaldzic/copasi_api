/**
 * @cond doxygen-libsbml-internal
 *
 * @file    NumericArgsMathCheck.h
 * @brief   Ensures arguments to numeric operators/functions are consistent.
 * @author  Sarah Keating
 *
 * $Id: NumericArgsMathCheck.h 11634 2010-08-03 03:57:18Z mhucka $
 * $HeadURL: https://sbml.svn.sourceforge.net/svnroot/sbml/trunk/libsbml/src/validator/constraints/NumericArgsMathCheck.h $
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


#ifndef NumericArgsMathCheck_h
#define NumericArgsMathCheck_h


#ifdef __cplusplus


#include <string>
#include <sstream>
#include <math.h>

#include <sbml/validator/VConstraint.h>

#include "MathMLBase.h"

LIBSBML_CPP_NAMESPACE_BEGIN

class ASTNode;


class NumericArgsMathCheck: public MathMLBase
{
public:

  /**
   * Creates a new Constraint with the given id.
   */
  NumericArgsMathCheck (unsigned int id, Validator& v);

  /**
   * Destroys this Constraint.
   */
  virtual ~NumericArgsMathCheck ();


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
   * Checks that the arguments of numeric functions are consistent
   *
   * If not, an error message is logged.
   */
  void checkNumericArgs (const Model& m, const ASTNode& node, const SBase & sb);

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
#endif  /* NumericArgsMathCheck_h */

/** @endcond */
