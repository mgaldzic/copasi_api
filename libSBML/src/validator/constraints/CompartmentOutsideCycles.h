/**
 * @cond doxygen-libsbml-internal
 *
 * @file    CompartmentOutsideCycles.h
 * @brief   Ensures no cycles exist via a Compartment's 'outside' attribute.
 * @author  Ben Bornstein
 *
 * $Id: CompartmentOutsideCycles.h 11634 2010-08-03 03:57:18Z mhucka $
 * $HeadURL: https://sbml.svn.sourceforge.net/svnroot/sbml/trunk/libsbml/src/validator/constraints/CompartmentOutsideCycles.h $
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


#ifndef CompartmentOutsideCycles_h
#define CompartmentOutsideCycles_h


#ifdef __cplusplus



#include <string>
#include <vector>

#include <algorithm>
#include <functional>

#include <sbml/validator/VConstraint.h>
#include "IdList.h"

LIBSBML_CPP_NAMESPACE_BEGIN

class Model;
class Compartment;
class Validator;


class CompartmentOutsideCycles: public TConstraint<Model>
{
public:

  /**
   * Creates a new Constraint with the given id.
   */
  CompartmentOutsideCycles (unsigned int id, Validator& v);

  /**
   * Destroys this Constraint.
   */
  virtual ~CompartmentOutsideCycles ();


protected:

  /**
   * Checks that no Compartments in Model have a cycle via their 'outside'
   * attribute.
   *
   * Sets mHolds to true if no cycles are found, false otherwise.
   */
  virtual void check_ (const Model& m, const Model& object);

  /**
   * Checks for a cycle by following Compartment c's 'outside' attribute.
   * If a cycle is found, it is added to the list of found cycles, mCycles.
   */
  void checkForCycle (const Model& m, const Compartment* c);

  /**
   * @return true if Compartment c is contained in one of the already found
   * cycles, false otherwise.
   */
  bool isInCycle (const Compartment* c);

  /**
   * Logs a message about a cycle found starting at Compartment c.
   */
  void logCycle (const Compartment* c, const IdList& cycle);


  std::vector<IdList> mCycles;
};

LIBSBML_CPP_NAMESPACE_END

#endif  /* __cplusplus */
#endif  /* CompartmentOutsideCycles_h */

/** @endcond */
