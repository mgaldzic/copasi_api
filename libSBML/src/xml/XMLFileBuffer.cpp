/**
 * @file    XMLFileBuffer.cpp
 * @brief   XMLFileBuffer implements the XMLBuffer interface
 * @author  Ben Bornstein
 * @author  Akiya Jouraku (replaced cstdio based code with std::istream based code
 * to support compressed files)
 *
 * $Id: XMLFileBuffer.cpp 11633 2010-08-03 03:53:27Z mhucka $
 * $HeadURL: https://sbml.svn.sourceforge.net/svnroot/sbml/trunk/libsbml/src/xml/XMLFileBuffer.cpp $
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

#include <cstdio>
#include<iostream>
#include<fstream>

#include <sbml/xml/XMLFileBuffer.h>
#include <sbml/compress/CompressCommon.h>
#include <sbml/compress/InputDecompressor.h>

/** @cond doxygen-ignored */

using namespace std;

/** @endcond */

LIBSBML_CPP_NAMESPACE_BEGIN

/** @cond doxygen-libsbml-internal */

/*
 * Creates a XMLBuffer based on the given file.  The file will be opened
 * for reading.
 */
XMLFileBuffer::XMLFileBuffer (const string& filename) :
   mFilename( filename )
{
  mStream = NULL;

  try
  {
    // open an uncompressed XML file
    if ( string::npos != filename.find(".xml", filename.length() -  4) )
    {
      mStream = new(std::nothrow) std::ifstream(filename.c_str());
    }
    // open a gzip file
    else if ( string::npos != filename.find(".gz", filename.length() -  3) )
    {
      mStream = InputDecompressor::openGzipIStream(filename);
    }
    // open a bz2 file
    else if ( string::npos != filename.find(".bz2", filename.length() - 4) )
    {
      mStream = InputDecompressor::openBzip2IStream(filename);
    }
    // open a zip file
    else if ( string::npos != filename.find(".zip", filename.length() - 4) )
    {
      mStream = InputDecompressor::openZipIStream(filename);
    }
    else
    {
      // open an uncompressed file
      mStream = new(std::nothrow) std::ifstream(filename.c_str());
    }
  }
  catch ( ZlibNotLinked& zlib)
  {
    // liBSBML is not linked with zlib.
    throw;
  }
  catch ( Bzip2NotLinked& bz2)
  {
    // liBSBML is not linked with bzip2.
    throw;
  }

  if(mStream != NULL)
  {
    // invoke peek() to set a badbit when the given compressed file is unreadable
    mStream->peek();
  }

}


/*
 * Destroys this XMLFileBuffer and closes the underlying file.
 */
XMLFileBuffer::~XMLFileBuffer ()
{
  if(mStream) delete mStream;
}


/*
 * Copies at most nbytes from this XMLFileBuffer to the memory pointed to
 * by destination.
 *
 * @return the number of bytes actually copied (may be 0).
 */
unsigned int
XMLFileBuffer::copyTo (void* destination, unsigned int bytes) 
{
  if (mStream)
  {
    mStream->read( static_cast<char*>(destination), bytes);
    return mStream->gcount();
  }
  else
  {
    return 0;
  }
}


/*
 * @return true if there was an error reading from the underlying buffer,
 * false otherwise.
 */
bool
XMLFileBuffer::error ()
{
  if (mStream) return (!mStream->eof() && mStream->fail());
  else return true;
}

/** @endcond */

LIBSBML_CPP_NAMESPACE_END

