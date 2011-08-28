/**
 * \file   TestFile.cpp
 * \brief  Enscapsulates an XML file in the test-data/ directory
 * \author Ben Bornstein
 * 
 * $Id: TestFile.cpp 11634 2010-08-03 03:57:18Z mhucka $
 * $HeadURL: https://sbml.svn.sourceforge.net/svnroot/sbml/trunk/libsbml/src/validator/test/TestFile.cpp $
 */
/* Copyright 2005 California Institute of Technology and
 * Japan Science and Technology Corporation.
 *
 * This library is free software; you can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation; either version 2.1 of the License, or
 * any later version.
 *
 * This library is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY, WITHOUT EVEN THE IMPLIED WARRANTY OF
 * MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE.  The software and
 * documentation provided hereunder is on an "as is" basis, and the
 * California Institute of Technology and Japan Science and Technology
 * Corporation have no obligations to provide maintenance, support,
 * updates, enhancements or modifications.  In no event shall the
 * California Institute of Technology or the Japan Science and Technology
 * Corporation be liable to any party for direct, indirect, special,
 * incidental or consequential damages, including lost profits, arising
 * out of the use of this software and its documentation, even if the
 * California Institute of Technology and/or Japan Science and Technology
 * Corporation have been advised of the possibility of such damage.  See
 * the GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this library; if not, write to the Free Software Foundation,
 * Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA.
 *
 * The original code contained here was initially developed by:
 *
 *     Ben Bornstein
 *     The Systems Biology Markup Language Development Group
 *     ERATO Kitano Symbiotic Systems Project
 *     Control and Dynamical Systems, MC 107-81
 *     California Institute of Technology
 *     Pasadena, CA, 91125, USA
 *
 *     http://www.sbml.org
 *     mailto:sbml-team@caltech.edu
 *
 * Contributor(s):
 */


#include <cstdlib>

#if defined(WIN32) && !defined(CYGWIN)
#  include "tps/dirent.h"
#else
#  include <sys/types.h>
#  include <dirent.h>
#endif  /* WIN32 */

#include <iostream>
#include "TestFile.h"

/** @cond doxygen-ignored */

using namespace std;

/** @endcond */


/**
 * TestFiles (e.g. in the test-data/ directory) have the following naming
 * convention:
 *
 *   ccccc-pass-00-nn.xml, or
 *   ccccc-fail-ff-nn.xml
 *
 *   ccccc-fail-ff-nn-xxxxx.xml
 *
 * Where:
 *
 *   ccccc  is the four digit constraint id the file is designed to test
 *
 *   pass   indicates the file must pass validation without error
 *
 *   fail   indicates the file must fail validation with extactly ff errors
 *          all with constraint id ccccc.
 *
 *   nn     is the sequence id (to allow multiple test files per constraint).
 *
 *   xxxxx  is the number of an additional constraint that this test fails
 *
 *
 * Offsets within mFilename:
 *
 *           1        1
 * 012345678901234567890123456
 * ccccc-pass-00-nn.xml
 * ccccc-fail-ff-nn.xml
 * ccccc-fail-ff-nn-xxxxx.xml
 */


std::string
TestFile::getFullname () const
{
  return mDirectory + "/" + mFilename;
}


unsigned int
TestFile::getConstraintId () const
{
  return atol( mFilename.substr(0, 5).c_str() );
}


unsigned int
TestFile::getSequenceId () const
{
  return atol( mFilename.substr(14, 2).c_str() );
}


unsigned int
TestFile::getNumFailures () const
{
  unsigned num = atol( mFilename.substr(11, 2).c_str() );
  if (mFilename.length() > 20)
    return num+1;
  else
    return num;
}

unsigned int
TestFile::getAdditionalFailId () const
{
  if (mFilename.length() > 20)
    return atol( mFilename.substr(17, 5).c_str() );
  else
    return 0;

}


/**
 * @return true if filename adheres to the TestFile naming convention,
 * false otherwise.
 */
bool
TestFile::isValid (const string& filename)
{
  return ((filename.length() == 20 && filename.substr(16, 4) == ".xml")
    || (filename.length() == 26 && filename.substr(22, 4) == ".xml"));
}


/**
 * @return the set of TestFiles in the given directory.
 *
 * You may optionally limit to the TestFiles returned to only those with
 * ConstraintIds in the range [begin, end] (if begin == end == 0, all
 * TestFiles in the given directory will be returned).
 */
set<TestFile>
TestFile::getFilesIn ( const string& directory,
                       unsigned int  begin,
                       unsigned int  end,
                       unsigned int  library)
{
  DIR*           dir;
  struct dirent* entry;
  set<TestFile>  result;

  dir = opendir( directory.c_str() );

  if (dir == NULL)
  {
    cerr << "Could not obtain a list of files in directory "
         << "[" << directory << "]." << endl;
    return result;
  }

  for (entry = readdir(dir); entry != NULL; entry = readdir(dir))
  {
    string filename(entry->d_name);

    if ( TestFile::isValid(filename) )
    {
      TestFile     file(directory, filename);
      unsigned int id = file.getConstraintId();

      //// leave out the following tests dependent on parser
      //if (library == 0)
      //{
      // // xerces bug in reading multibyte chars
      //  if (id == 10309) continue;

        // libxml bug for 2.6.16 on a Mac
#ifdef BUGGY_APPLE_LIBXML
      unsigned int num = file.getSequenceId();

      if (id == 1013 && num <= 10) continue;
      if (id == 1014 && num <= 4) continue;
#endif

      // leave out model constraints moved to units
      if (id == 20702 || id == 20616
        || id == 20511 || id == 20512 || id == 20513)
        continue;
      //}

      if ((begin == 0 && end == 0) || (id >= begin && id <= end))
      {
        result.insert( TestFile(directory, filename) );
      }
    }
  }

  closedir(dir);

  return result;
}
