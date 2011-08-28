/**
 * \file    TestRDFAnnotation.cpp
 * \brief   fomula units data unit tests
 * \author  Ben Bornstein
 *
 * $Id: TestRDFAnnotation2.cpp 11376 2010-07-03 05:41:44Z mhucka $
 * $HeadURL: https://sbml.svn.sourceforge.net/svnroot/sbml/trunk/libsbml/src/annotation/test/TestRDFAnnotation2.cpp $
 */
/* Copyright 2002 California Institute of Technology and Japan Science and
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


#include <sbml/common/common.h>
#include <sbml/common/extern.h>

#include <sbml/SBMLReader.h>
#include <sbml/SBMLTypes.h>

#include <sbml/SBMLDocument.h>
#include <sbml/Model.h>
#include <sbml/SBMLTypeCodes.h>

#include <sbml/annotation/RDFAnnotation.h>
#include <sbml/annotation/ModelHistory.h>

#include <check.h>

LIBSBML_CPP_NAMESPACE_USE

Model *m2;
SBMLDocument* d2;

extern char *TestDataDirectory;

/* 
 * tests the results from rdf annotations
 */
CK_CPPSTART

void
RDFAnnotation2_setup (void)
{
  char *filename = safe_strcat(TestDataDirectory, "annotation2.xml");

  // The following will return a pointer to a new SBMLDocument.
  d2 = readSBML(filename);
  m2 = d2->getModel();
}


void
RDFAnnotation2_teardown (void)
{
  delete d2;
}

static bool
equals (const char* expected, const char* actual)
{
  if ( !strcmp(expected, actual) ) return true;

  printf( "\nStrings are not equal:\n"  );
  printf( "Expected:\n[%s]\n", expected );
  printf( "Actual:\n[%s]\n"  , actual   );

  return false;
}



START_TEST (test_RDFAnnotation2_getModelHistory)
{
  ModelHistory * history = m2->getModelHistory();

  fail_unless(history != NULL);

  ModelCreator * mc = (ModelCreator * )(history->getListCreators()->get(0));

  fail_unless(!strcmp(ModelCreator_getFamilyName(mc), "Hucka"));
  fail_unless(!strcmp(ModelCreator_getGivenName(mc), "Mike"));
  fail_unless(!strcmp(ModelCreator_getEmail(mc), "mhucka@caltech.edu"));
  fail_unless(!strcmp(ModelCreator_getOrganisation(mc), "BNMC"));

  ModelCreator * mc1 = (ModelCreator * )(history->getListCreators()->get(1));

  fail_unless(!strcmp(ModelCreator_getFamilyName(mc1), "Keating"));
  fail_unless(!strcmp(ModelCreator_getGivenName(mc1), "Sarah"));
  fail_unless(!strcmp(ModelCreator_getEmail(mc1), "skeating@caltech.edu"));
  fail_unless(!strcmp(ModelCreator_getOrganisation(mc1), "UH"));

  Date * date = history->getCreatedDate();
  fail_unless(Date_getYear(date) == 2005);
  fail_unless(Date_getMonth(date) == 2);
  fail_unless(Date_getDay(date) == 2);
  fail_unless(Date_getHour(date) == 14);
  fail_unless(Date_getMinute(date) == 56);
  fail_unless(Date_getSecond(date) == 11);
  fail_unless(Date_getSignOffset(date) == 0);
  fail_unless(Date_getHoursOffset(date) == 0);
  fail_unless(Date_getMinutesOffset(date) == 0);
  fail_unless(!strcmp(Date_getDateAsString(date), "2005-02-02T14:56:11Z"));

  date = history->getModifiedDate();
  fail_unless(Date_getYear(date) == 2006);
  fail_unless(Date_getMonth(date) == 5);
  fail_unless(Date_getDay(date) == 30);
  fail_unless(Date_getHour(date) == 10);
  fail_unless(Date_getMinute(date) == 46);
  fail_unless(Date_getSecond(date) == 2);
  fail_unless(Date_getSignOffset(date) == 0);
  fail_unless(Date_getHoursOffset(date) == 0);
  fail_unless(Date_getMinutesOffset(date) == 0);
  fail_unless(!strcmp(Date_getDateAsString(date), "2006-05-30T10:46:02Z"));

  date = history->getModifiedDate(1);
  fail_unless(Date_getYear(date) == 2007);
  fail_unless(Date_getMonth(date) == 1);
  fail_unless(Date_getDay(date) == 16);
  fail_unless(Date_getHour(date) == 15);
  fail_unless(Date_getMinute(date) == 31);
  fail_unless(Date_getSecond(date) == 52);
  fail_unless(Date_getSignOffset(date) == 0);
  fail_unless(Date_getHoursOffset(date) == 0);
  fail_unless(Date_getMinutesOffset(date) == 0);
  fail_unless(!strcmp(Date_getDateAsString(date), "2007-01-16T15:31:52Z"));
}
END_TEST


START_TEST (test_RDFAnnotation2_modelWithHistoryAndCVTerms)
{
  ModelHistory * h = new ModelHistory();

  ModelCreator *c = new ModelCreator();
  c->setFamilyName("Keating");
  c->setGivenName("Sarah");

  h->addCreator(c);

  Date *d = new Date(2008, 11, 17, 18, 37, 0, 0, 0, 0);
  h->setCreatedDate(d);
  h->setModifiedDate(d);

  m2->unsetModelHistory();

  m2->setModelHistory(h);

  CVTerm *cv = new CVTerm();
  cv->setQualifierType(BIOLOGICAL_QUALIFIER);
  cv->setBiologicalQualifierType(BQB_IS_VERSION_OF);
  cv->addResource("http://www.geneontology.org/#GO:0005892");

  m2->addCVTerm(cv);
  XMLNode *ann = RDFAnnotationParser::parseModelHistory(m2);

  const char * expected =
    "<annotation>\n"
		"  <rdf:RDF xmlns:rdf=\"http://www.w3.org/1999/02/22-rdf-syntax-ns#\" xmlns:dc=\"http://purl.org/dc/elements/1.1/\" xmlns:dcterms=\"http://purl.org/dc/terms/\" xmlns:vCard=\"http://www.w3.org/2001/vcard-rdf/3.0#\" xmlns:bqbiol=\"http://biomodels.net/biology-qualifiers/\" xmlns:bqmodel=\"http://biomodels.net/model-qualifiers/\">\n"
		"    <rdf:Description rdf:about=\"#_000001\">\n"
		"      <dc:creator rdf:parseType=\"Resource\">\n"
		"        <rdf:Bag>\n"
		"          <rdf:li rdf:parseType=\"Resource\">\n"
		"            <vCard:N rdf:parseType=\"Resource\">\n"
		"              <vCard:Family>Keating</vCard:Family>\n"
		"              <vCard:Given>Sarah</vCard:Given>\n"
		"            </vCard:N>\n"
		"          </rdf:li>\n"
		"        </rdf:Bag>\n"
		"      </dc:creator>\n"
    "      <dcterms:created rdf:parseType=\"Resource\">\n"
    "        <dcterms:W3CDTF>2008-11-17T18:37:00Z</dcterms:W3CDTF>\n"
    "      </dcterms:created>\n"
    "      <dcterms:modified rdf:parseType=\"Resource\">\n"
    "        <dcterms:W3CDTF>2008-11-17T18:37:00Z</dcterms:W3CDTF>\n"
    "      </dcterms:modified>\n"
		"      <bqbiol:isVersionOf>\n"
		"        <rdf:Bag>\n"
		"          <rdf:li rdf:resource=\"http://www.geneontology.org/#GO:0005892\"/>\n"
		"        </rdf:Bag>\n"
		"      </bqbiol:isVersionOf>\n"
		"    </rdf:Description>\n"
		"  </rdf:RDF>\n"
    "</annotation>";

  if (ann != NULL) 
  {
    fail_unless( equals(expected, ann->toXMLString().c_str()) );
  }
  else
  {
    fail("parseModelHistory failed");
  }
}
END_TEST


START_TEST (test_RDFAnnotation2_modelWithHistoryAndMultipleModifiedDates)
{
  ModelHistory * h = new ModelHistory();

  ModelCreator *c = new ModelCreator();
  c->setFamilyName("Keating");
  c->setGivenName("Sarah");

  h->addCreator(c);

  Date * d = new Date(2005, 2, 2, 14, 56, 11);
  h->setCreatedDate(d);
  h->addModifiedDate(d);
  h->addModifiedDate(d);
  m2->unsetModelHistory();

  m2->setModelHistory(h);

  XMLNode *ann = RDFAnnotationParser::parseModelHistory(m2);

  const char * expected =
    "<annotation>\n"
		"  <rdf:RDF xmlns:rdf=\"http://www.w3.org/1999/02/22-rdf-syntax-ns#\" xmlns:dc=\"http://purl.org/dc/elements/1.1/\" xmlns:dcterms=\"http://purl.org/dc/terms/\" xmlns:vCard=\"http://www.w3.org/2001/vcard-rdf/3.0#\" xmlns:bqbiol=\"http://biomodels.net/biology-qualifiers/\" xmlns:bqmodel=\"http://biomodels.net/model-qualifiers/\">\n"
		"    <rdf:Description rdf:about=\"#_000001\">\n"
		"      <dc:creator rdf:parseType=\"Resource\">\n"
		"        <rdf:Bag>\n"
		"          <rdf:li rdf:parseType=\"Resource\">\n"
		"            <vCard:N rdf:parseType=\"Resource\">\n"
		"              <vCard:Family>Keating</vCard:Family>\n"
		"              <vCard:Given>Sarah</vCard:Given>\n"
		"            </vCard:N>\n"
		"          </rdf:li>\n"
		"        </rdf:Bag>\n"
		"      </dc:creator>\n"
		"      <dcterms:created rdf:parseType=\"Resource\">\n"
		"        <dcterms:W3CDTF>2005-02-02T14:56:11Z</dcterms:W3CDTF>\n"
		"      </dcterms:created>\n"
		"      <dcterms:modified rdf:parseType=\"Resource\">\n"
		"        <dcterms:W3CDTF>2005-02-02T14:56:11Z</dcterms:W3CDTF>\n"
		"      </dcterms:modified>\n"
		"      <dcterms:modified rdf:parseType=\"Resource\">\n"
		"        <dcterms:W3CDTF>2005-02-02T14:56:11Z</dcterms:W3CDTF>\n"
		"      </dcterms:modified>\n"
		"    </rdf:Description>\n"
		"  </rdf:RDF>\n"
    "</annotation>";

  fail_unless( equals(expected, ann->toXMLString().c_str()) );
}
END_TEST

START_TEST (test_RDFAnnotation2_modelWithHistoryWithCharacterReference)

{
  ModelHistory * h = new ModelHistory();

  ModelCreator *c = new ModelCreator();
  c->setFamilyName("Dr&#228;ger");
  c->setGivenName("Andreas");

  h->addCreator(c);
  Date * d = new Date(2005, 2, 2, 14, 56, 11);
  h->setCreatedDate(d);
  h->addModifiedDate(d);

  m2->unsetModelHistory();

  m2->setModelHistory(h);

  XMLNode *ann = RDFAnnotationParser::parseModelHistory(m2);

  const char * expected =
    "<annotation>\n"
		"  <rdf:RDF xmlns:rdf=\"http://www.w3.org/1999/02/22-rdf-syntax-ns#\" xmlns:dc=\"http://purl.org/dc/elements/1.1/\" xmlns:dcterms=\"http://purl.org/dc/terms/\" xmlns:vCard=\"http://www.w3.org/2001/vcard-rdf/3.0#\" xmlns:bqbiol=\"http://biomodels.net/biology-qualifiers/\" xmlns:bqmodel=\"http://biomodels.net/model-qualifiers/\">\n"
		"    <rdf:Description rdf:about=\"#_000001\">\n"
		"      <dc:creator rdf:parseType=\"Resource\">\n"
		"        <rdf:Bag>\n"
		"          <rdf:li rdf:parseType=\"Resource\">\n"
		"            <vCard:N rdf:parseType=\"Resource\">\n"
		"              <vCard:Family>Dr&#228;ger</vCard:Family>\n"
		"              <vCard:Given>Andreas</vCard:Given>\n"
		"            </vCard:N>\n"
		"          </rdf:li>\n"
		"        </rdf:Bag>\n"
		"      </dc:creator>\n"
		"      <dcterms:created rdf:parseType=\"Resource\">\n"
		"        <dcterms:W3CDTF>2005-02-02T14:56:11Z</dcterms:W3CDTF>\n"
		"      </dcterms:created>\n"
		"      <dcterms:modified rdf:parseType=\"Resource\">\n"
		"        <dcterms:W3CDTF>2005-02-02T14:56:11Z</dcterms:W3CDTF>\n"
		"      </dcterms:modified>\n"
		"    </rdf:Description>\n"
		"  </rdf:RDF>\n"
    "</annotation>";

  fail_unless( equals(expected, ann->toXMLString().c_str()) );
}
END_TEST

Suite *
create_suite_RDFAnnotation2 (void)
{
  Suite *suite = suite_create("RDFAnnotation2");
  TCase *tcase = tcase_create("RDFAnnotation2");

  tcase_add_checked_fixture(tcase,
                            RDFAnnotation2_setup,
                            RDFAnnotation2_teardown);

  tcase_add_test(tcase, test_RDFAnnotation2_getModelHistory );
  tcase_add_test(tcase, test_RDFAnnotation2_modelWithHistoryAndCVTerms );
  tcase_add_test(tcase, test_RDFAnnotation2_modelWithHistoryAndMultipleModifiedDates );
  tcase_add_test(tcase, test_RDFAnnotation2_modelWithHistoryWithCharacterReference);
  suite_add_tcase(suite, tcase);

  return suite;
}


CK_CPPEND
