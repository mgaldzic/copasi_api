/**
 * @cond doxygen-libsbml-internal
 *
 * @file    LocalParameterShadowsIdInModel.h
 * @brief   Checks whether local parameters are shadowing another id
 * @author  Sarah Keating
 *
 * $Id: LocalParameterShadowsIdInModel.h 11634 2010-08-03 03:57:18Z mhucka $
 * $HeadURL: https://sbml.svn.sourceforge.net/svnroot/sbml/trunk/libsbml/src/validator/constraints/LocalParameterShadowsIdInModel.h $
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


#ifndef LocalParameterShadowsIdInModel_h
#define LocalParameterShadowsIdInModel_h


#ifdef __cplusplus

#include <string>

#include <sbml/validator/VConstraint.h>

#include "IdList.h"

LIBSBML_CPP_NAMESPACE_BEGIN

class LocalParameterShadowsIdInModel: public TConstraint<Model>
{
public:

  /**
   * Creates a new Constraint with the given constraint id.
   */
  LocalParameterShadowsIdInModel (unsigned int id, Validator& v);

  /**
   * Destroys this Constraint.
   */
  virtual ~LocalParameterShadowsIdInModel ();


protected:

  /**
   * Checks that any species with boundary condition false
   * is not set by reaction and rules
   */
  virtual void check_ (const Model& m, const Model& object);

  /**
   * Logs a message about local parameter being shadowed
   */
  void logConflict (const Parameter& p, const SBase& object);

  IdList mAll;

};

LIBSBML_CPP_NAMESPACE_END

#endif  /* __cplusplus */
#endif  /* LocalParameterShadowsIdInModel_h */

/** @endcond */
