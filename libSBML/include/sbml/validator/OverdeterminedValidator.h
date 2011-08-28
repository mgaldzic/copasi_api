/**
 * @cond doxygen-libsbml-internal
 *
 * @file    OverdeterminedValidator.h
 * @brief   Performs consistency checks on an SBML model
 * @author  Sarah Keating
 *
 * $Id: OverdeterminedValidator.h 11634 2010-08-03 03:57:18Z mhucka $
 * $HeadURL: https://sbml.svn.sourceforge.net/svnroot/sbml/trunk/libsbml/src/validator/OverdeterminedValidator.h $
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
 * in the file named "LICENSE.txt" included with this software distribution and
 * also available online as http://sbml.org/software/libsbml/license.html
 *----------------------------------------------------------------------- -->*/

#ifndef OverdeterminedValidator_h
#define OverdeterminedValidator_h


#ifdef __cplusplus


#include <sbml/validator/Validator.h>
#include <sbml/SBMLError.h>

LIBSBML_CPP_NAMESPACE_BEGIN

class OverdeterminedValidator: public Validator
{
public:

  OverdeterminedValidator () :
    Validator( LIBSBML_CAT_OVERDETERMINED_MODEL ) { }

  virtual ~OverdeterminedValidator () { }

  /**
   * Initializes this Validator with a set of Constraints.
   */
  virtual void init ();
};

LIBSBML_CPP_NAMESPACE_END

#endif  /* __cplusplus */
#endif  /* OverdeterminedValidator_h */


/** @endcond */
