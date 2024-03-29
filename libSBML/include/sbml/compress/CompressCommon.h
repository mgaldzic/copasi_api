/**
 *@cond doxygen-libsbml-internal
 **
 *
 * @file    CompressCommon.h
 * @brief   common classes/functions for compression/decompression I/O
 * @author  Akiya Jouraku
 *
 * $Id: CompressCommon.h 11633 2010-08-03 03:53:27Z mhucka $
 * $HeadURL: https://sbml.svn.sourceforge.net/svnroot/sbml/trunk/libsbml/src/compress/CompressCommon.h $
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
 *----------------------------------------------------------------------- -->*/

#ifndef CompressCommon_h
#define CompressCommon_h

#include <exception>
#include <sbml/common/extern.h>


LIBSBML_CPP_NAMESPACE_BEGIN

/**
 *
 *  This exception will be thrown if a function which depends on
 *  some library invoked and underlying libSBML is not linked with 
 *  the library.
 *
 */
class LIBSBML_EXTERN NotLinked : public std::exception
{
public:
   NotLinked() throw() { }
   virtual ~NotLinked() throw() {}
};


/**
 *
 * This exception will be thrown if a function which depends on 
 * zlib invoked and underlying libSBML is not linked with zlib.
 *
 */
class LIBSBML_EXTERN ZlibNotLinked : public NotLinked
{
public:
   ZlibNotLinked() throw() { }
   virtual ~ZlibNotLinked() throw() {}
};


/**
 *
 *  This exception will be thrown if a function which depends on 
 *  bzip2 library invoked and underlying libSBML is not linked with
 *  bzip2.
 *
 */
class LIBSBML_EXTERN Bzip2NotLinked : public NotLinked
{
public:
   Bzip2NotLinked() throw() { }
   virtual ~Bzip2NotLinked() throw() {}
};


/**
 * Predicate returning @c true or @c false depending on whether
 * underlying libSBML is linked with zlib.
 *
 * @return @c true if libSBML is linked with zlib, @c false otherwise.
 */
LIBSBML_EXTERN
bool hasZlib(); 


/**
 * Predicate returning @c true or @c false depending on whether
 * underlying libSBML is linked with bzip2.
 *
 * @return @c true if libSBML is linked with bzip2, @c false otherwise.
 */
LIBSBML_EXTERN
bool hasBzip2();

LIBSBML_CPP_NAMESPACE_END

#endif //CompressCommon_h

/** @endcond */
