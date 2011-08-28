/**
 * \file    TestCVTerms.cpp
 * \brief   CVTerms unit tests
 * \author  Sarah Keating
 *
 * $Id: TestCVTerms.c 10125 2009-08-28 12:14:18Z sarahkeating $
 * $HeadURL: https://sbml.svn.sourceforge.net/svnroot/sbml/trunk/libsbml/src/annotation/test/TestCVTerms.c $
 */
/* Copyright 2007 California Institute of Technology.
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

#include <sbml/common/common.h>
#include <sbml/annotation/CVTerm.h>
#include <sbml/xml/XMLNode.h>
#include <sbml/xml/XMLTriple.h>

#include <check.h>


START_TEST (test_CVTerm_create)
{
  CVTerm_t *term = CVTerm_createWithQualifierType(MODEL_QUALIFIER);

  fail_unless(term != NULL);
  fail_unless(CVTerm_getQualifierType(term) == MODEL_QUALIFIER);
  CVTerm_free(term);

}
END_TEST

START_TEST (test_CVTerm_set_get)
{
  CVTerm_t *term = CVTerm_createWithQualifierType(MODEL_QUALIFIER);

  fail_unless(term != NULL);
  fail_unless(CVTerm_getQualifierType(term) == MODEL_QUALIFIER);
  
  CVTerm_setModelQualifierType(term, BQM_IS);

  fail_unless(term != NULL);
  fail_unless(CVTerm_getQualifierType(term) == MODEL_QUALIFIER);
  fail_unless(CVTerm_getModelQualifierType(term) == BQM_IS);

  CVTerm_setQualifierType(term, BIOLOGICAL_QUALIFIER);
  CVTerm_setBiologicalQualifierType( term, BQB_IS);

  fail_unless(CVTerm_getQualifierType(term) == BIOLOGICAL_QUALIFIER);
  fail_unless(CVTerm_getBiologicalQualifierType(term) == BQB_IS);

  
  
  CVTerm_free(term);
}
END_TEST

START_TEST (test_CVTerm_createFromNode)
{
  XMLAttributes_t * xa;
  XMLTriple_t * qual_triple = XMLTriple_createWith ("is", "", "bqbiol");
  XMLTriple_t * bag_triple = XMLTriple_create ();
  XMLTriple_t * li_triple = XMLTriple_create();
  XMLAttributes_t * att = XMLAttributes_create ();
  XMLAttributes_add(att, "", "This is my resource");
  XMLAttributes_t *att1 = XMLAttributes_create();

  XMLToken_t * li_token = XMLToken_createWithTripleAttr(li_triple, att);
  XMLToken_t * bag_token = XMLToken_createWithTripleAttr(bag_triple, att1);
  XMLToken_t * qual_token = XMLToken_createWithTripleAttr(qual_triple, att1);

  XMLNode_t * li = XMLNode_createFromToken(li_token);
  XMLNode_t * bag = XMLNode_createFromToken(bag_token);
  XMLNode_t * node = XMLNode_createFromToken(qual_token);

  XMLNode_addChild(bag, li);
  XMLNode_addChild(node, bag);


  CVTerm_t *term = CVTerm_createFromNode(node);

  fail_unless(term != NULL);
  fail_unless(CVTerm_getQualifierType(term) == BIOLOGICAL_QUALIFIER);
  fail_unless(CVTerm_getBiologicalQualifierType(term) == BQB_IS);

  xa = CVTerm_getResources(term);

  fail_unless(XMLAttributes_getLength(xa) == 1);
  fail_unless(!strcmp(XMLAttributes_getName(xa, 0), "rdf:resource"));
  fail_unless(!strcmp(XMLAttributes_getValue(xa, 0), "This is my resource"));

  XMLTriple_free(qual_triple);
  XMLTriple_free(bag_triple);
  XMLTriple_free(li_triple);
  XMLToken_free(li_token);
  XMLToken_free(bag_token);
  XMLToken_free(qual_token);
  XMLAttributes_free(att);
  XMLAttributes_free(att1);
  CVTerm_free(term);
  XMLNode_free(node);
  XMLNode_free(bag);
  XMLNode_free(li);
  

}
END_TEST

START_TEST (test_CVTerm_addResource)
{
  CVTerm_t *term = CVTerm_createWithQualifierType(MODEL_QUALIFIER);
  const char * resource = "GO6666";
  XMLAttributes_t *xa;

  fail_unless(term != NULL);
  fail_unless(CVTerm_getQualifierType(term) == MODEL_QUALIFIER);
  
  CVTerm_addResource(term, resource);

  xa = CVTerm_getResources(term);

  fail_unless(XMLAttributes_getLength(xa) == 1);
  fail_unless(!strcmp(XMLAttributes_getName(xa, 0), "rdf:resource"));
  fail_unless(!strcmp(XMLAttributes_getValue(xa, 0), "GO6666"));


  
  CVTerm_free(term);

}
END_TEST


START_TEST (test_CVTerm_getResources)
{
  CVTerm_t *term = CVTerm_createWithQualifierType(MODEL_QUALIFIER);
  const char * resource = "GO6666";
  const char * resource1 = "OtherURI";
  unsigned int number;
  
  CVTerm_addResource(term, resource);
  CVTerm_addResource(term, resource1);

  number = CVTerm_getNumResources(term);

  fail_unless(number == 2);
  fail_unless(!strcmp(CVTerm_getResourceURI(term, 0), "GO6666"));
  fail_unless(!strcmp(CVTerm_getResourceURI(term, 1), "OtherURI"));


  
  CVTerm_free(term);

}
END_TEST


Suite *
create_suite_CVTerms (void)
{
  Suite *suite = suite_create("CVTerms");
  TCase *tcase = tcase_create("CVTerms");

  tcase_add_test( tcase, test_CVTerm_create  );
  tcase_add_test( tcase, test_CVTerm_set_get  );
  tcase_add_test( tcase, test_CVTerm_addResource  );
  tcase_add_test( tcase, test_CVTerm_createFromNode  );
  tcase_add_test( tcase, test_CVTerm_getResources  );
  suite_add_tcase(suite, tcase);

  return suite;
}


END_C_DECLS
