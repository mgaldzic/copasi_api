/**
 * @cond doxygen-libsbml-internal
 *
 * @file    SBMLTransforms.h
 * @brief   Transform functions
 * @author  Sarah Keating
 *
 * $Id: SBMLTransforms.h 11672 2010-08-06 08:05:50Z mhucka $
 * $HeadURL: https://sbml.svn.sourceforge.net/svnroot/sbml/trunk/libsbml/src/sbml/SBMLTransforms.h $
 *
 *<!---------------------------------------------------------------------------
 * This file is part of libSBML.  Please visit http://sbml.org for more
 * information about SBML, and the latest version of libSBML.
 *
 * Copyright 2005-2010 California Institute of Technology.
 * Copyright 2002-2005 California Institute of Technology and
 *                     Japan Science and Technology Corporation.
 * 
 * This library is free software; you can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation.  A copy of the license agreement is provided
 * in the file named "LICENSE.txt" included with this software distribution
 * and also available online as http://sbml.org/software/libsbml/license.html
 *----------------------------------------------------------------------- -->
 *
 * @class SBMLTransforms
 * @brief Methods for transform elements of SBML
 *
 */

#ifndef SBMLTransforms_h
#define SBMLTransforms_h


#include <sbml/common/extern.h>
#include <sbml/SBase.h>
#include <sbml/util/util.h>
#include <sbml/math/ASTNode.h>
#include <sbml/FunctionDefinition.h>
#include <sbml/InitialAssignment.h>
#include <sbml/SBMLTypes.h>
#include <sbml/SpeciesReference.h>

#ifdef __cplusplus


#include <string>

LIBSBML_CPP_NAMESPACE_BEGIN

class IdList;

class LIBSBML_EXTERN SBMLTransforms
{
public:

#ifndef SWIG
  typedef std::pair<double, bool>   ValueSet;
  typedef std::map<const std::string, ValueSet> IdValueMap;
  typedef IdValueMap::iterator                  IdValueIter;
#endif

  /** @cond doxygen-libsbml-internal */
  
  /**
   * Expands the math represented by the ASTNode to implement the functionality
   * of the FunctionDefinition, if it occurs within the original
   * math.
   *
   * For example, an ASTNode represents the math expression: f(s, p) where
   * f is the id of a FunctionDefinition representing f(x, y) = x * y.
   * The outcome of the function is that the ASTNode now represents
   * the math expression: s * p
   *
   * @param math ASTNode representing the math to be transformed
   *
   * @param fd the FunctionDefinition to be expanded
   */
  static void replaceFD(ASTNode * math, const FunctionDefinition * fd);

  /** @endcond */

  /** @cond doxygen-libsbml-internal */

  /**
   * Expands the math represented by the ASTNode to implement the functionality
   * of all the FunctionDefinitions in the list, if they occur within the 
   * original math.
   *
   * For example, an ASTNode represents the math expression: f(s, g(p, q)) where
   * f is the id of a FunctionDefinition representing f(x, y) = x * y
   * and g is the id of a FunctionDefinition representing f(x, y) = x/y
   * The outcome of the function is that the ASTNode now represents
   * the math expression: s * p/q
   *
   * @param math ASTNode representing the math to be transformed
   *
   * @param lofd the ListOfFunctionDefinitions to be expanded
   */
  static void replaceFD(ASTNode * math, const ListOfFunctionDefinitions * lofd);

  /** @endcond */


  /** @cond doxygen-libsbml-internal */
  
  static bool expandInitialAssignments(Model * m);

  /** @endcond */

 
  /** @cond doxygen-libsbml-internal */
  
  static double evaluateASTNode(const ASTNode * node);

  /** @endcond */

  /** @cond doxygen-libsbml-internal */
  
  static IdList mapComponentValues(const Model * m);

  /** @endcond */

protected:  
  /** @cond doxygen-libsbml-internal */

  static void replaceBvars(ASTNode * math, const FunctionDefinition * fd);

  /** @endcond */

  /** @cond doxygen-libsbml-internal */
  
  static bool checkFunctionNodeForIds(ASTNode * node, IdList& ids);

  /** @endcond */
 
  /** @cond doxygen-libsbml-internal */
  
  static bool nodeContainsId(const ASTNode * node, IdList& ids);

  /** @endcond */
 
  /** @cond doxygen-libsbml-internal */
  
  static bool nodeContainsNameNotInList(const ASTNode * node, IdList& ids);

  /** @endcond */

  /** @cond doxygen-libsbml-internal */
  
  static bool expandInitialAssignment(Parameter * p, 
                                          const InitialAssignment *ia);
  /** @endcond */

  /** @cond doxygen-libsbml-internal */
  
  static bool expandInitialAssignment(Compartment * c, 
                                          const InitialAssignment *ia);
  /** @endcond */
  
  /** @cond doxygen-libsbml-internal */
  
  static bool expandInitialAssignment(SpeciesReference * sr, 
                                          const InitialAssignment *ia);
  /** @endcond */

  
  /** @cond doxygen-libsbml-internal */
  
  static bool expandInitialAssignment(Species * s, 
                                          const InitialAssignment *ia);
  /** @endcond */



  static IdValueMap mValues;

};

LIBSBML_CPP_NAMESPACE_END

#endif  /* __cplusplus */

#endif  /* SBMLTransforms_h */


/** @endcond */
