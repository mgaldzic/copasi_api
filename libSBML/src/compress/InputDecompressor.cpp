/**
 *@cond doxygen-libsbml-internal
 **
 *
 * @file    InputDecompressor.cpp
 * @brief   utility class for input decompression
 * @author  Akiya Jouraku
 *
 * $Id: InputDecompressor.cpp 11633 2010-08-03 03:53:27Z mhucka $
 * $HeadURL: https://sbml.svn.sourceforge.net/svnroot/sbml/trunk/libsbml/src/compress/InputDecompressor.cpp $
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

#include <istream>
#include <fstream>
#include <sstream>
#include <iostream>
#include <new>
#include <algorithm>
#include <cstring>

#include <sbml/compress/InputDecompressor.h>

#ifdef USE_ZLIB
#include <sbml/compress/zfstream.h>
#include <sbml/compress/zipfstream.h>
#endif //USE_ZLIB

#ifdef USE_BZ2
#include <sbml/compress/bzfstream.h>
#endif //USE_BZ2

using namespace std;

LIBSBML_CPP_NAMESPACE_BEGIN

/**
 * Opens the given gzip file as a gzifstream (subclass of std::ifstream class) object
 * for read access and returned the stream object.
 *
 * @return a istream* object bound to the given gzip file or NULL if the initialization
 * for the object failed.
 */
std::istream* 
InputDecompressor::openGzipIStream (const std::string& filename)
{
#ifdef USE_ZLIB
  return new(std::nothrow) gzifstream(filename.c_str(), ios_base::in | ios_base::binary);
#else
  throw ZlibNotLinked();
  return NULL; // never reached
#endif
}


/**
 * Opens the given bzip2 file as a bzifstream (subclass of std::ifstream class) object
 * for read access and returned the stream object.
 *
 * @return a istream* object bound to the given bzip2 file or NULL if the initialization
 * for the object failed.
 */
std::istream* 
InputDecompressor::openBzip2IStream (const std::string& filename)
{
#ifdef USE_BZ2
  return new(std::nothrow) bzifstream(filename.c_str(), ios_base::in | ios_base::binary);
#else
  throw Bzip2NotLinked();
  return NULL; // never reached
#endif
}


/**
 * Opens the given zip file as a zipifstream (subclass of std::ifstream class) object
 * for read access and returned the stream object.
 *
 * @return a istream* object bound to the given zip file or NULL if the initialization
 * for the object failed.
 *
 * @note The first file in the given zip archive file will be opened if the zip archive
 * contains two or more files.
 */
std::istream* 
InputDecompressor::openZipIStream (const std::string& filename)
{
#ifdef USE_ZLIB
  return new(std::nothrow) zipifstream(filename.c_str(), ios_base::in | ios_base::binary);
#else
  throw ZlibNotLinked();
  return NULL; // never reached
#endif
}


/**
 * Opens the given gzip file and returned the string in the file.
 *
 * @return a string, the string in the given file, or empty string if 
 * failed to open the file.
 */
char* 
InputDecompressor::getStringFromGzip (const std::string& filename) 
{
#ifdef USE_ZLIB
  std::ostringstream oss;
  gzifstream in(filename.c_str(), ios_base::in | ios_base::binary);
  istreambuf_iterator<char> in_itr(in);
  ostreambuf_iterator<char> out_itr(oss);

  std::copy(in_itr, istreambuf_iterator<char>(), out_itr);

  return strdup(oss.str().c_str());
#else
  throw ZlibNotLinked();
  return NULL; // never reached
#endif
}


/**
 * Opens the given bzip2 file and returned the string in the file.
 *
 * @return a string, the string in the given file, or empty string if failed to open
 * the file.
 */
char* 
InputDecompressor::getStringFromBzip2 (const std::string& filename) 
{
#ifdef USE_BZ2
  std::ostringstream oss;
  bzifstream in(filename.c_str(), ios_base::in | ios_base::binary);
  istreambuf_iterator<char> in_itr(in);
  ostreambuf_iterator<char> out_itr(oss);

  std::copy(in_itr, istreambuf_iterator<char>(), out_itr);

  return strdup(oss.str().c_str());
#else
  throw Bzip2NotLinked();
  return NULL; // never reached
#endif
}


/**
 * Opens the given zip file and returned the string in the file.
 *
 * @return a string, the string in the given file, or empty string if failed to open
 * the file.
 *
 * @note The first file in the given zip archive file will be opened if the zip archive
 * contains two or more files.
 */
char* 
InputDecompressor::getStringFromZip (const std::string& filename) 
{
#ifdef USE_ZLIB
  std::ostringstream oss;
  zipifstream in(filename.c_str(), ios_base::in | ios_base::binary);
  istreambuf_iterator<char> in_itr(in);
  ostreambuf_iterator<char> out_itr(oss);

  std::copy(in_itr, istreambuf_iterator<char>(), out_itr);

  return strdup(oss.str().c_str());
#else
  throw ZlibNotLinked();
  return NULL; // never reached
#endif
}

LIBSBML_CPP_NAMESPACE_END

/** @endcond */
