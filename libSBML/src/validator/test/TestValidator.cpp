/**
 * \file   TestValidator.cpp
 * \brief  Validator unit tests
 * \author Ben Bornstein
 * 
 * $Id: TestValidator.cpp 11634 2010-08-03 03:57:18Z mhucka $
 * $HeadURL: https://sbml.svn.sourceforge.net/svnroot/sbml/trunk/libsbml/src/validator/test/TestValidator.cpp $
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


#include <iostream>
#include <sstream>

#include <vector>

#include <algorithm>
#include <functional>
#include <iterator>

#include <sbml/SBMLDocument.h>
#include <sbml/SBMLReader.h>
#include <sbml/validator/Validator.h>

#include "TestFile.h"

#include "TestValidator.h"

/** @cond doxygen-ignored */

using namespace std;
LIBSBML_CPP_NAMESPACE_USE

/** @endcond */


TestValidator::TestValidator (Validator& v) : mValidator(v), mVerboseAll(false)
{
  readVerboseFromEnvironment();
}


TestValidator::~TestValidator ()
{
}


/**
 * Function Object: Return true if the given SBMLError has the given
 * id, false otherwise.
 */
struct HasId : public unary_function<SBMLError, bool>
{
  unsigned int id;

  HasId (unsigned int id) : id(id) { }
  bool operator() (const SBMLError& msg) { return msg.getErrorId() == id; }
};


/**
 * Function Object: Takes a SBMLError and returns its integer id.
 */
struct ToId : public unary_function<SBMLError, unsigned int>
{
  unsigned int operator() (const SBMLError& msg) { return msg.getErrorId(); }
};


/**
 * @return true if the Validator behaved as expected when validating
 * TestFile, false otherwise.
 */
bool
TestValidator::test (const TestFile& file)
{
  bool error = false;

  unsigned int id       = file.getConstraintId();

  /* change numbers for specific units tests that report same number */
  if (id == 99502 || id == 99503 || id == 99504)
    id = 10501;

  /* change numbers for specific units tests that report same number */
  if (id == 90502 || id == 90503 || id == 90504)
    id = 90501;

  unsigned int expected = file.getNumFailures();
  unsigned int others   = file.getAdditionalFailId();

  unsigned int actual   = mValidator.validate( file.getFullname() );

  list<SBMLError>::const_iterator begin = mValidator.getFailures().begin();
  list<SBMLError>::const_iterator end   = mValidator.getFailures().end();



  if (expected != actual)
  {
    error = true;

    cout << endl;
    cout << "Error: " << file.getFilename() << endl;
    cout << "  - Failures:  Expected: "  << expected << "  Actual: " << actual;
    cout << endl << endl;
  }


  unsigned int same = count_if(begin, end, HasId(id));
  vector<unsigned int> ids;

  if (expected != same && actual != same)
  {
    // need to consider case where the test case has
    // an additional fail
    if (expected - same != 1)
    {
      error = true;
    }
    else
    {
      transform(begin, end, back_inserter(ids), ToId());

      unsigned int match = 0;
      for (unsigned int i = 0; i < ids.size(); i++)
      {
        if (ids.at(i) == others)
        {
          match = 1;
        }
      }

      if (match == 0)
      {
        error = true;
      }
    }
  }

  if (error)
  {      
    cout << endl;
    cout << "Error: " << file.getFilename() << endl;
    cout << "  - Constraints:  Expected: " << id << "  Actual: ";
    cout << endl;
    copy(ids.begin(), ids.end(), ostream_iterator<unsigned int>(cout, " "));
    cout << endl;
  }

  if ( error || isVerbose(id) )
  {
    copy(begin, end, ostream_iterator<SBMLError>(cout, "\n"));
  }

  mValidator.clearFailures();

  return error == false;
}


/**
 * @return true if the all Validator failures for the given Constraint id
 * test cases should be printed, false otherwise.
 */
bool
TestValidator::isVerbose (unsigned int id)
{
  vector<unsigned int>::iterator begin = mVerboseConstraintIds.begin();
  vector<unsigned int>::iterator end   = mVerboseConstraintIds.end();


  return mVerboseAll || find(begin, end, id) != end;
}


/**
 * Reads the environment variable 'LIBSBML_TEST_VALIDATOR_VERBOSE' and sets
 * the internal verbose state accordingly.
 *
 * LIBSBML_TEST_VALIDATOR_VERBOSE may either be set to 'all' (case
 * sensitive) or a space separated list of constraint ids.
 */
void
TestValidator::readVerboseFromEnvironment ()
{
  const char* s = getenv("LIBSBML_TEST_VALIDATOR_VERBOSE");


  if (!s)
  {
    return;
  }
  else if (string(s) == "all")
  {
    mVerboseAll = true;
  }
  else
  {
    unsigned int  id;
    istringstream is(s);

    while ( !is.eof() )
    {
      is >> id;
      mVerboseConstraintIds.push_back(id);
    }
  }
}
