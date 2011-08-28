/**
 * @file    TestXMLError.c
 * @brief   XMLError unit tests, C version
 * @author  Sarah Keating
 *
 * $Id: TestXMLErrorC.c 10866 2010-01-29 19:52:27Z mhucka $
 * $HeadURL: https://sbml.svn.sourceforge.net/svnroot/sbml/trunk/libsbml/src/xml/test/TestXMLErrorC.c $
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

#include <sbml/common/common.h>
#include <sbml/xml/XMLError.h>

#include <check.h>


START_TEST (test_XMLError_create_C)
{
  XMLError_t *error = XMLError_create();
  fail_unless(error != NULL);
  fail_unless(XMLError_isInfo(error) == 0);
  fail_unless(XMLError_isWarning(error) == 0);
  fail_unless(XMLError_isError(error) == 0);
  fail_unless(XMLError_isFatal(error) == 1);
  XMLError_free(error);

  error = XMLError_createWithIdAndMessage(12345, "My message");
  fail_unless( strcmp(XMLError_getMessage(error), "My message") == 0 );
  fail_unless( XMLError_getErrorId(error) == 12345 );
  XMLError_free(error);
}
END_TEST


START_TEST (test_XMLError_variablesAsStrings)
{
  XMLError_t *error = XMLError_createWithIdAndMessage(1003, "");
  
  fail_unless( XMLError_getErrorId(error)  == 1003 );
  fail_unless( XMLError_getSeverity(error) == LIBSBML_SEV_ERROR );
  fail_unless( !strcmp(XMLError_getSeverityAsString(error), "Error") );
  fail_unless( XMLError_getCategory(error) == LIBSBML_CAT_XML );
  fail_unless( !strcmp(XMLError_getCategoryAsString(error), "XML content"));

  XMLError_free(error);
}
END_TEST


Suite *
create_suite_XMLError_C (void)
{
  Suite *suite = suite_create("XMLErrorC");
  TCase *tcase = tcase_create("XMLErrorC");

  tcase_add_test( tcase, test_XMLError_create_C  );
  tcase_add_test( tcase, test_XMLError_variablesAsStrings  );
  suite_add_tcase(suite, tcase);

  return suite;
}

