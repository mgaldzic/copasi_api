/**
 * \file    TestSBase.cpp
 * \brief   SBase unit tests
 * \author  Ben Bornstein
 *
 * $Id: TestSBase.cpp 11572 2010-07-26 07:02:55Z mhucka $
 * $HeadURL: https://sbml.svn.sourceforge.net/svnroot/sbml/trunk/libsbml/src/sbml/test/TestSBase.cpp $
 */
/* Copyright 2002 California Institute of Technology and
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
 *     http://www.cds.caltech.edu/erato
 *     mailto:sbml-team@caltech.edu
 *
 * Contributor(s):
 */


#include <sbml/common/common.h>
#include <sbml/common/extern.h>

#include <sbml/SBase.h>
#include <sbml/Model.h>
#include <sbml/annotation/CVTerm.h>

#include <check.h>

LIBSBML_CPP_NAMESPACE_USE

static SBase *S;


BEGIN_C_DECLS


void
SBaseTest_setup (void)
{
  S = new(std::nothrow) Model(2, 4);

  if (S == NULL)
  {
    fail("'new(std::nothrow) SBase;' returned a NULL pointer.");
  }

}


void
SBaseTest_teardown (void)
{
  delete S;
}


START_TEST (test_SBase_setMetaId)
{
  const char *metaid = "x12345";


  SBase_setMetaId(S, metaid);

  fail_unless( !strcmp(SBase_getMetaId(S), metaid), NULL );
  fail_unless( SBase_isSetMetaId(S)      , NULL );

  if (SBase_getMetaId(S) == metaid)
  {
    fail("SBase_setMetaId(...) did not make a copy of string.");
  }

  /* Reflexive case (pathological) */
  SBase_setMetaId(S, SBase_getMetaId(S));
  fail_unless( !strcmp(SBase_getMetaId(S), metaid), NULL );

  SBase_setMetaId(S, NULL);
  fail_unless( !SBase_isSetMetaId(S), NULL );

  if (SBase_getMetaId(S) != NULL)
  {
    fail("SBase_setMetaId(S, NULL) did not clear string.");
  }
}
END_TEST


START_TEST (test_SBase_setNotes)
{
  SBase_t *c = new(std::nothrow) Model(1, 2);
  XMLToken_t *token;
  XMLNode_t *node;

  token = XMLToken_createWithText("This is a test note");
  node = XMLNode_createFromToken(token);


  SBase_setNotes(c, node);

  fail_unless(SBase_isSetNotes(c) == 1);

  if (SBase_getNotes(c) == node)
  {
    fail("SBase_setNotes(...) did not make a copy of node.");
  }
  XMLNode_t *t1 = SBase_getNotes(c);

  fail_unless(XMLNode_getNumChildren(t1) == 1);
  fail_unless(!strcmp(XMLNode_getCharacters(XMLNode_getChild(t1,0)), "This is a test note"));


  /* Reflexive case (pathological)  */
  SBase_setNotes(c, SBase_getNotes(c));
  t1 = SBase_getNotes(c);
  fail_unless(XMLNode_getNumChildren(t1) == 1);
  const char * chars = XMLNode_getCharacters(XMLNode_getChild(t1,0));
  fail_unless(!strcmp(chars, "This is a test note"));

  SBase_setNotes(c, NULL);
  fail_unless(SBase_isSetNotes(c) == 0 );

  if (SBase_getNotes(c) != NULL)
  {
    fail("SBase_setNotes(c, NULL) did not clear string.");
  }

  SBase_setNotes(c, node);

  fail_unless(SBase_isSetNotes(c) == 1);

  /* test notes with character reference */

  token = XMLToken_createWithText("(CR) &#0168; &#x00a8; &#x00A8; (NOT CR) &#; &#x; &#00a8; &#0168 &#x00a8");
  node  = XMLNode_createFromToken(token);

  SBase_setNotes(c, node);
  t1 = SBase_getNotes(c);

  fail_unless(XMLNode_getNumChildren(t1) == 1);

  const char * s = XMLNode_toXMLString(XMLNode_getChild(t1,0));
  const char * expected = "(CR) &#0168; &#x00a8; &#x00A8; (NOT CR) &amp;#; &amp;#x; &amp;#00a8; &amp;#0168 &amp;#x00a8";

  fail_unless(!strcmp(s,expected));

  /* test notes with predefined entity */

  token = XMLToken_createWithText("& ' > < \" &amp; &apos; &gt; &lt; &quot;");
  node  = XMLNode_createFromToken(token);

  SBase_setNotes(c, node);
  t1 = SBase_getNotes(c);

  fail_unless(XMLNode_getNumChildren(t1) == 1);

  const char * s2 = XMLNode_toXMLString(XMLNode_getChild(t1,0));
  const char * expected2 = "&amp; &apos; &gt; &lt; &quot; &amp; &apos; &gt; &lt; &quot;";

  fail_unless(!strcmp(s2,expected2));

  XMLToken_free(token);
  XMLNode_free(node);
}
END_TEST


START_TEST (test_SBase_setAnnotation)
{
  XMLToken_t *token;
  XMLNode_t *node;

  token = XMLToken_createWithText("This is a test note");
  node = XMLNode_createFromToken(token);

  SBase_setAnnotation(S, node);

  fail_unless(SBase_isSetAnnotation(S) == 1);

  XMLNode_t *t1 = SBase_getAnnotation(S);

  fail_unless(XMLNode_getNumChildren(t1) == 1);
  fail_unless(!strcmp(XMLNode_getCharacters(XMLNode_getChild(t1,0)), "This is a test note"));

  if (SBase_getAnnotation(S) == node)
  {
    fail("SBase_setAnnotation(...) did not make a copy of node.");
  }

  /* Reflexive case (pathological) */
  SBase_setAnnotation(S, SBase_getAnnotation(S));
  fail_unless(!strcmp(XMLNode_getCharacters(XMLNode_getChild(SBase_getAnnotation(S),0)), "This is a test note"));

  SBase_setAnnotation(S, NULL);
  fail_unless(SBase_isSetAnnotation(S) == 0 );

  if (SBase_getAnnotation(S) != NULL)
  {
    fail("SBase_setAnnotation(S, NULL) did not clear string.");
  }

  SBase_setAnnotation(S, node);

  fail_unless(SBase_isSetAnnotation(S) == 1);

  SBase_unsetAnnotation(S);

  fail_unless(SBase_isSetAnnotation(S) == 0);

  /* test annotations with character reference */

  token = XMLToken_createWithText("(CR) &#0168; &#x00a8; &#x00A8; (NOT CR) &#; &#x; &#00a8; &#0168 &#x00a8");
  node  = XMLNode_createFromToken(token);

  SBase_setAnnotation(S, node);
  t1 = SBase_getAnnotation(S);

  fail_unless(XMLNode_getNumChildren(t1) == 1);

  const char * s = XMLNode_toXMLString(XMLNode_getChild(t1,0));
  const char * expected = "(CR) &#0168; &#x00a8; &#x00A8; (NOT CR) &amp;#; &amp;#x; &amp;#00a8; &amp;#0168 &amp;#x00a8";

  fail_unless(!strcmp(s,expected));

  /* test notes with predefined entity */

  token = XMLToken_createWithText("& ' > < \" &amp; &apos; &gt; &lt; &quot;");
  node  = XMLNode_createFromToken(token);

  SBase_setAnnotation(S, node);
  t1 = SBase_getAnnotation(S);

  fail_unless(XMLNode_getNumChildren(t1) == 1);

  const char * s2 = XMLNode_toXMLString(XMLNode_getChild(t1,0));
  const char * expected2 = "&amp; &apos; &gt; &lt; &quot; &amp; &apos; &gt; &lt; &quot;";

  fail_unless(!strcmp(s2,expected2));

  XMLToken_free(token);
  XMLNode_free(node);

}
END_TEST


START_TEST (test_SBase_unsetAnnotationWithCVTerms)
{
  CVTerm_t   *cv;

  const char *annt =
        "<annotation>\n"
        "  <test:test xmlns:test=\"http://test.org/test\">this is a test node</test:test>\n"
        "</annotation>";

  const char *annt_with_cvterm =
        "<annotation>\n"
        "  <test:test xmlns:test=\"http://test.org/test\">this is a test node</test:test>\n"
        "  <rdf:RDF xmlns:rdf=\"http://www.w3.org/1999/02/22-rdf-syntax-ns#\" "
				"xmlns:dc=\"http://purl.org/dc/elements/1.1/\" "
        "xmlns:dcterms=\"http://purl.org/dc/terms/\" "
				"xmlns:vCard=\"http://www.w3.org/2001/vcard-rdf/3.0#\" "
        "xmlns:bqbiol=\"http://biomodels.net/biology-qualifiers/\" "
				"xmlns:bqmodel=\"http://biomodels.net/model-qualifiers/\">\n"
        "    <rdf:Description rdf:about=\"#_000001\">\n"
        "      <bqbiol:is>\n"
        "        <rdf:Bag>\n"
        "          <rdf:li rdf:resource=\"http://www.geneontology.org/#GO:0005895\"/>\n"
        "        </rdf:Bag>\n"
        "      </bqbiol:is>\n"
        "    </rdf:Description>\n"
        "  </rdf:RDF>\n"
        "</annotation>";

  SBase_setAnnotationString(S, (char*)annt);
  fail_unless(SBase_isSetAnnotation(S) == 1);
  fail_unless(!strcmp(SBase_getAnnotationString(S), annt));

  SBase_unsetAnnotation(S);
  fail_unless(SBase_isSetAnnotation(S) == 0);
  fail_unless(SBase_getAnnotation(S) == NULL);

  SBase_setAnnotationString(S, (char*)annt);
  SBase_setMetaId(S, "_000001");
  cv = CVTerm_createWithQualifierType(BIOLOGICAL_QUALIFIER);
  CVTerm_setBiologicalQualifierType(cv, BQB_IS);
  CVTerm_addResource(cv, "http://www.geneontology.org/#GO:0005895");
  SBase_addCVTerm(S, cv);

  fail_unless(SBase_isSetAnnotation(S) == 1);
  fail_unless(!strcmp(SBase_getAnnotationString(S), annt_with_cvterm));

  SBase_unsetAnnotation(S);
  fail_unless(SBase_isSetAnnotation(S) == 0);
  fail_unless(SBase_getAnnotation(S) == NULL);

  CVTerm_free(cv);
}
END_TEST

START_TEST (test_SBase_unsetAnnotationWithModelHistory)
{
  ModelHistory_t *h = ModelHistory_create();
  ModelCreator_t *c = ModelCreator_create();
  Date_t *dc;
  Date_t *dm;

  const char *annt =
        "<annotation>\n"
        "  <test:test xmlns:test=\"http://test.org/test\">this is a test node</test:test>\n"
        "</annotation>";

  const char *annt_with_modelhistory =
       "<annotation>\n"
       "  <test:test xmlns:test=\"http://test.org/test\">this is a test node</test:test>\n"
       "  <rdf:RDF xmlns:rdf=\"http://www.w3.org/1999/02/22-rdf-syntax-ns#\" "
			 "xmlns:dc=\"http://purl.org/dc/elements/1.1/\" "
			 "xmlns:dcterms=\"http://purl.org/dc/terms/\" "
			 "xmlns:vCard=\"http://www.w3.org/2001/vcard-rdf/3.0#\" "
			 "xmlns:bqbiol=\"http://biomodels.net/biology-qualifiers/\" "
			 "xmlns:bqmodel=\"http://biomodels.net/model-qualifiers/\">\n"
       "    <rdf:Description rdf:about=\"#_000001\">\n"
       "      <dc:creator>\n"
       "        <rdf:Bag>\n"
       "          <rdf:li rdf:parseType=\"Resource\">\n"
       "            <vCard:N rdf:parseType=\"Resource\">\n"
       "              <vCard:Family>Keating</vCard:Family>\n"
       "              <vCard:Given>Sarah</vCard:Given>\n"
       "            </vCard:N>\n"
       "            <vCard:EMAIL>sbml-team@caltech.edu</vCard:EMAIL>\n"
       "          </rdf:li>\n"
       "        </rdf:Bag>\n"
       "      </dc:creator>\n"
       "      <dcterms:created rdf:parseType=\"Resource\">\n"
       "        <dcterms:W3CDTF>2005-12-29T12:15:45+02:00</dcterms:W3CDTF>\n"
       "      </dcterms:created>\n"
       "      <dcterms:modified rdf:parseType=\"Resource\">\n"
       "        <dcterms:W3CDTF>2005-12-30T12:15:45+02:00</dcterms:W3CDTF>\n"
       "      </dcterms:modified>\n"
       "    </rdf:Description>\n"
       "  </rdf:RDF>\n"
       "</annotation>";

  SBase_setAnnotationString(S, (char*)annt);
  fail_unless(SBase_isSetAnnotation(S) == 1);
  fail_unless(!strcmp(SBase_getAnnotationString(S), annt));

  SBase_unsetAnnotation(S);
  fail_unless(SBase_isSetAnnotation(S) == 0);
  fail_unless(SBase_getAnnotation(S) == NULL);

  SBase_setAnnotationString(S, (char*)annt);
  SBase_setMetaId(S, "_000001");

  ModelCreator_setFamilyName(c, (char *) "Keating");
  ModelCreator_setGivenName(c, (char *) "Sarah");
  ModelCreator_setEmail(c, (char *) "sbml-team@caltech.edu");

  ModelHistory_addCreator(h, c);
  dc =  Date_createFromValues(2005, 12, 29, 12, 15, 45, 1, 2, 0);
  ModelHistory_setCreatedDate(h, dc);
  dm =  Date_createFromValues(2005, 12, 30, 12, 15, 45, 1, 2, 0);
  ModelHistory_setModifiedDate(h, dm);
  Model_setModelHistory((Model_t*)S, h);

  fail_unless(SBase_isSetAnnotation(S) == 1);
  fail_unless(!strcmp(SBase_getAnnotationString(S), annt_with_modelhistory));

  SBase_unsetAnnotation(S);
  fail_unless(SBase_isSetAnnotation(S) == 0);
  fail_unless(SBase_getAnnotation(S) == NULL);

  ModelCreator_free(c);
  ModelHistory_free(h);
}
END_TEST


START_TEST (test_SBase_setNotesString)
{
  SBase_t *c = new(std::nothrow) Model(1, 2);
  char * notes = "This is a test note";
  char * taggednotes = "<notes>This is a test note</notes>";

  SBase_setNotesString(c, notes);

  fail_unless(SBase_isSetNotes(c) == 1);

  if (strcmp(SBase_getNotesString(c), taggednotes))
  {
    fail("SBase_setNotesString(...) did not make a copy of node.");
  }
  XMLNode_t *t1 = SBase_getNotes(c);
  fail_unless(XMLNode_getNumChildren(t1) == 1);

  const XMLNode_t *t2 = XMLNode_getChild(t1,0);
  fail_unless(!strcmp(XMLNode_getCharacters(t2), "This is a test note"));


  /* Reflexive case (pathological)  */
  SBase_setNotesString(c, SBase_getNotesString(c));
  t1 = SBase_getNotes(c);
  fail_unless(XMLNode_getNumChildren(t1) == 1);
  const char * chars = SBase_getNotesString(c);
  fail_unless(!strcmp(chars, taggednotes));

  SBase_setNotesString(c, (char *)"");
  fail_unless(SBase_isSetNotes(c) == 0 );

  if (SBase_getNotesString(c) != NULL)
  {
    fail("SBase_getNotesString(c, "") did not clear string.");
  }

  SBase_setNotesString(c, taggednotes);

  fail_unless(SBase_isSetNotes(c) == 1);

  if (strcmp(SBase_getNotesString(c), taggednotes))
  {
    fail("SBase_setNotesString(...) did not make a copy of node.");
  }
  t1 = SBase_getNotes(c);
  fail_unless(XMLNode_getNumChildren(t1) == 1);

  t2 = XMLNode_getChild(t1,0);
  fail_unless(!strcmp(XMLNode_getCharacters(t2), "This is a test note"));

}
END_TEST


START_TEST (test_SBase_setAnnotationString)
{
  char * annotation = "This is a test note";
  char * taggedannotation = "<annotation>This is a test note</annotation>";

  SBase_setAnnotationString(S, annotation);

  fail_unless(SBase_isSetAnnotation(S) == 1);

  if (strcmp(SBase_getAnnotationString(S), taggedannotation))
  {
    fail("SBase_setAnnotationString(...) did not make a copy of node.");
  }
  XMLNode_t *t1 = SBase_getAnnotation(S);
  fail_unless(XMLNode_getNumChildren(t1) == 1);

  fail_unless(!strcmp(XMLNode_getCharacters(XMLNode_getChild(t1,0)), "This is a test note"));


  /* Reflexive case (pathological)  */
  SBase_setAnnotationString(S, SBase_getAnnotationString(S));
  t1 = SBase_getAnnotation(S);
  fail_unless(XMLNode_getNumChildren(t1) == 1);
  const char * chars = SBase_getAnnotationString(S);
  fail_unless(!strcmp(chars, taggedannotation));

  SBase_setAnnotationString(S, "");
  fail_unless(SBase_isSetAnnotation(S) == 0 );

  if (SBase_getAnnotationString(S) != NULL)
  {
    fail("SBase_getAnnotationString(S, "") did not clear string.");
  }

  SBase_setAnnotationString(S, taggedannotation);

  fail_unless(SBase_isSetAnnotation(S) == 1);

  if (strcmp(SBase_getAnnotationString(S), taggedannotation))
  {
    fail("SBase_setAnnotationString(...) did not make a copy of node.");
  }
  t1 = SBase_getAnnotation(S);
  fail_unless(XMLNode_getNumChildren(t1) == 1);

  const XMLNode_t *t2 = XMLNode_getChild(t1,0);
  fail_unless(!strcmp(XMLNode_getCharacters(t2), "This is a test note"));
}
END_TEST


START_TEST (test_SBase_appendNotes)
{ // add a p tag to a p tag
  XMLToken_t *token;
  XMLNode_t *node;
  XMLToken_t *token1;
  XMLNode_t *node1;
  XMLNode_t * node2;
  XMLTriple_t *triple = XMLTriple_createWith("p", "", "");
  XMLAttributes_t * att = XMLAttributes_create ();
  XMLNamespaces_t *ns = XMLNamespaces_create();
  XMLNamespaces_add(ns, "http://www.w3.org/1999/xhtml", "");
  XMLToken_t *token4 = XMLToken_createWithText("This is my text");
  XMLNode_t *node4 = XMLNode_createFromToken(token4);
  XMLToken_t *token5 = XMLToken_createWithText("This is additional text");
  XMLNode_t *node5 = XMLNode_createFromToken(token5);

  token = XMLToken_createWithTripleAttrNS(triple, att, ns);
  node = XMLNode_createFromToken(token);
  XMLNode_addChild(node, node4);

  SBase_setNotes(S, node);

  fail_unless(SBase_isSetNotes(S) == 1);

  token1 = XMLToken_createWithTripleAttrNS(triple, att, ns);
  node1 = XMLNode_createFromToken(token1);
  XMLNode_addChild(node1, node5);
  
  SBase_appendNotes(S, node1);

  fail_unless(SBase_isSetNotes(S) == 1);

  node2 = SBase_getNotes(S);

  fail_unless(XMLNode_getNumChildren(node2) == 2);
  fail_unless(!strcmp(XMLNode_getName(XMLNode_getChild(node2, 0)), "p"));
  fail_unless(XMLNode_getNumChildren(XMLNode_getChild(node2, 0)) == 1);
  fail_unless(!strcmp(XMLNode_getName(XMLNode_getChild(node2, 1)), "p"));
  fail_unless(XMLNode_getNumChildren(XMLNode_getChild(node2, 1)) == 1);

  const char * chars1 = XMLNode_getCharacters(XMLNode_getChild(
    XMLNode_getChild(node2, 0), 0));
  const char * chars2 = XMLNode_getCharacters(XMLNode_getChild(
    XMLNode_getChild(node2, 1), 0));

  fail_unless(!strcmp(chars1, "This is my text"));
  fail_unless(!strcmp(chars2, "This is additional text"));

  XMLNode_free(node);
  XMLNode_free(node1);
}
END_TEST


START_TEST (test_SBase_appendNotes1)
{
  // add a html tag to an html tag
  XMLAttributes_t * att = XMLAttributes_create ();
  XMLNamespaces_t *ns = XMLNamespaces_create();
  XMLNamespaces_add(ns, "http://www.w3.org/1999/xhtml", "");
  XMLTriple_t *html_triple = XMLTriple_createWith("html", "", "");
  XMLTriple_t *head_triple = XMLTriple_createWith("head", "", "");
  XMLTriple_t *title_triple = XMLTriple_createWith("title", "", "");
  XMLTriple_t *body_triple = XMLTriple_createWith("body", "", "");
  XMLTriple_t *p_triple = XMLTriple_createWith("p", "", "");
  XMLToken_t *html_token = XMLToken_createWithTripleAttrNS(html_triple, att, ns);
  XMLToken_t *head_token = XMLToken_createWithTripleAttr(head_triple, att);
  XMLToken_t *title_token = XMLToken_createWithTripleAttr(title_triple, att);
  XMLToken_t *body_token = XMLToken_createWithTripleAttr(body_triple, att);
  XMLToken_t *p_token = XMLToken_createWithTripleAttr(p_triple, att);
  XMLToken_t *text_token = XMLToken_createWithText("This is my text");
  XMLNode_t *html_node = XMLNode_createFromToken(html_token);
  XMLNode_t *head_node = XMLNode_createFromToken(head_token);
  XMLNode_t *title_node = XMLNode_createFromToken(title_token);
  XMLNode_t *body_node = XMLNode_createFromToken(body_token);
  XMLNode_t *p_node = XMLNode_createFromToken(p_token);
  XMLNode_t *text_node = XMLNode_createFromToken(text_token);

  XMLToken_t *text_token1 = XMLToken_createWithText("This is more text");
  XMLNode_t *html_node1 = XMLNode_createFromToken(html_token);
  XMLNode_t *head_node1 = XMLNode_createFromToken(head_token);
  XMLNode_t *title_node1 = XMLNode_createFromToken(title_token);
  XMLNode_t *body_node1 = XMLNode_createFromToken(body_token);
  XMLNode_t *p_node1 = XMLNode_createFromToken(p_token);
  XMLNode_t *text_node1 = XMLNode_createFromToken(text_token1);

  XMLNode_t * notes;
  const XMLNode_t *child, *child1;

  XMLNode_addChild(p_node, text_node);
  XMLNode_addChild(body_node, p_node);
  XMLNode_addChild(head_node, title_node);
  XMLNode_addChild(html_node, head_node);
  XMLNode_addChild(html_node, body_node);

  XMLNode_addChild(p_node1, text_node1);
  XMLNode_addChild(body_node1, p_node1);
  XMLNode_addChild(head_node1, title_node1);
  XMLNode_addChild(html_node1, head_node1);
  XMLNode_addChild(html_node1, body_node1);

  SBase_setNotes(S, html_node);
  SBase_appendNotes(S, html_node1);

  notes = SBase_getNotes(S);

  fail_unless(!strcmp(XMLNode_getName(notes), "notes"));
  fail_unless(XMLNode_getNumChildren(notes) == 1);

  child = XMLNode_getChild(notes, 0);

  fail_unless(!strcmp(XMLNode_getName(child), "html"));
  fail_unless(XMLNode_getNumChildren(child) == 2);

  child = XMLNode_getChild(child, 1);

  fail_unless(!strcmp(XMLNode_getName(child), "body"));
  fail_unless(XMLNode_getNumChildren(child) == 2);

  child1 = XMLNode_getChild(child, 0);

  fail_unless(!strcmp(XMLNode_getName(child1), "p"));
  fail_unless(XMLNode_getNumChildren(child1) == 1);

  child1 = XMLNode_getChild(child1, 0);

  fail_unless(!strcmp(XMLNode_getCharacters(child1), "This is my text"));
  fail_unless(XMLNode_getNumChildren(child1) == 0);

  child1 = XMLNode_getChild(child, 1);

  fail_unless(!strcmp(XMLNode_getName(child1), "p"));
  fail_unless(XMLNode_getNumChildren(child1) == 1);

  child1 = XMLNode_getChild(child1, 0);

  fail_unless(!strcmp(XMLNode_getCharacters(child1), "This is more text"));
  fail_unless(XMLNode_getNumChildren(child1) == 0);

  XMLAttributes_free(att);
  XMLNamespaces_free(ns);
  XMLTriple_free(html_triple);
  XMLTriple_free(head_triple);
  XMLTriple_free(body_triple);
  XMLTriple_free(p_triple);
  XMLToken_free(html_token);
  XMLToken_free(head_token);
  XMLToken_free(body_token);
  XMLToken_free(p_token);
  XMLToken_free(text_token);
  XMLToken_free(text_token1);
  XMLNode_free(html_node);
  XMLNode_free(head_node);
  XMLNode_free(body_node);
  XMLNode_free(p_node);
  XMLNode_free(text_node);
  XMLNode_free(html_node1);
  XMLNode_free(head_node1);
  XMLNode_free(body_node1);
  XMLNode_free(p_node1);
  XMLNode_free(text_node1);
}
END_TEST


START_TEST (test_SBase_appendNotes2)
{// add a body tag to an html tag
  XMLAttributes_t * att = XMLAttributes_create ();
  XMLNamespaces_t *ns = XMLNamespaces_create();
  XMLNamespaces_add(ns, "http://www.w3.org/1999/xhtml", "");
  XMLTriple_t *html_triple = XMLTriple_createWith("html", "", "");
  XMLTriple_t *head_triple = XMLTriple_createWith("head", "", "");
  XMLTriple_t *title_triple = XMLTriple_createWith("title", "", "");
  XMLTriple_t *body_triple = XMLTriple_createWith("body", "", "");
  XMLTriple_t *p_triple = XMLTriple_createWith("p", "", "");
  XMLToken_t *html_token = XMLToken_createWithTripleAttrNS(html_triple, att, ns);
  XMLToken_t *head_token = XMLToken_createWithTripleAttr(head_triple, att);
  XMLToken_t *title_token = XMLToken_createWithTripleAttr(title_triple, att);
  XMLToken_t *body_token = XMLToken_createWithTripleAttr(body_triple, att);
  XMLToken_t *p_token = XMLToken_createWithTripleAttr(p_triple, att);
  XMLToken_t *text_token = XMLToken_createWithText("This is my text");
  XMLNode_t *html_node = XMLNode_createFromToken(html_token);
  XMLNode_t *head_node = XMLNode_createFromToken(head_token);
  XMLNode_t *title_node = XMLNode_createFromToken(title_token);
  XMLNode_t *body_node = XMLNode_createFromToken(body_token);
  XMLNode_t *p_node = XMLNode_createFromToken(p_token);
  XMLNode_t *text_node = XMLNode_createFromToken(text_token);

  XMLToken_t *body_token1 = XMLToken_createWithTripleAttrNS(body_triple, att, ns);
  XMLToken_t *text_token1 = XMLToken_createWithText("This is more text");
  XMLNode_t *body_node1 = XMLNode_createFromToken(body_token1);
  XMLNode_t *p_node1 = XMLNode_createFromToken(p_token);
  XMLNode_t *text_node1 = XMLNode_createFromToken(text_token1);

  XMLNode_t * notes;
  const XMLNode_t *child, *child1;

  XMLNode_addChild(p_node, text_node);
  XMLNode_addChild(body_node, p_node);
  XMLNode_addChild(head_node, title_node);
  XMLNode_addChild(html_node, head_node);
  XMLNode_addChild(html_node, body_node);

  XMLNode_addChild(p_node1, text_node1);
  XMLNode_addChild(body_node1, p_node1);

  SBase_setNotes(S, html_node);
  SBase_appendNotes(S, body_node1);

  notes = SBase_getNotes(S);

  fail_unless(!strcmp(XMLNode_getName(notes), "notes"));
  fail_unless(XMLNode_getNumChildren(notes) == 1);

  child = XMLNode_getChild(notes, 0);

  fail_unless(!strcmp(XMLNode_getName(child), "html"));
  fail_unless(XMLNode_getNumChildren(child) == 2);

  child = XMLNode_getChild(child, 1);

  fail_unless(!strcmp(XMLNode_getName(child), "body"));
  fail_unless(XMLNode_getNumChildren(child) == 2);

  child1 = XMLNode_getChild(child, 0);

  fail_unless(!strcmp(XMLNode_getName(child1), "p"));
  fail_unless(XMLNode_getNumChildren(child1) == 1);

  child1 = XMLNode_getChild(child1, 0);

  fail_unless(!strcmp(XMLNode_getCharacters(child1), "This is my text"));
  fail_unless(XMLNode_getNumChildren(child1) == 0);

  child1 = XMLNode_getChild(child, 1);

  fail_unless(!strcmp(XMLNode_getName(child1), "p"));
  fail_unless(XMLNode_getNumChildren(child1) == 1);

  child1 = XMLNode_getChild(child1, 0);

  fail_unless(!strcmp(XMLNode_getCharacters(child1), "This is more text"));
  fail_unless(XMLNode_getNumChildren(child1) == 0);

  XMLAttributes_free(att);
  XMLNamespaces_free(ns);
  XMLTriple_free(html_triple);
  XMLTriple_free(head_triple);
  XMLTriple_free(body_triple);
  XMLTriple_free(p_triple);
  XMLToken_free(html_token);
  XMLToken_free(head_token);
  XMLToken_free(body_token);
  XMLToken_free(p_token);
  XMLToken_free(text_token);
  XMLToken_free(text_token1);
  XMLToken_free(body_token1);
  XMLNode_free(html_node);
  XMLNode_free(head_node);
  XMLNode_free(body_node);
  XMLNode_free(p_node);
  XMLNode_free(text_node);
  XMLNode_free(body_node1);
  XMLNode_free(p_node1);
  XMLNode_free(text_node1);
}
END_TEST


START_TEST (test_SBase_appendNotes3)
{
  // add a p tag to an html tag
  XMLAttributes_t * att = XMLAttributes_create ();
  XMLNamespaces_t *ns = XMLNamespaces_create();
  XMLNamespaces_add(ns, "http://www.w3.org/1999/xhtml", "");
  XMLTriple_t *html_triple = XMLTriple_createWith("html", "", "");
  XMLTriple_t *head_triple = XMLTriple_createWith("head", "", "");
  XMLTriple_t *title_triple = XMLTriple_createWith("title", "", "");
  XMLTriple_t *body_triple = XMLTriple_createWith("body", "", "");
  XMLTriple_t *p_triple = XMLTriple_createWith("p", "", "");
  XMLToken_t *html_token = XMLToken_createWithTripleAttrNS(html_triple, att, ns);
  XMLToken_t *head_token = XMLToken_createWithTripleAttr(head_triple, att);
  XMLToken_t *title_token = XMLToken_createWithTripleAttr(title_triple, att);
  XMLToken_t *body_token = XMLToken_createWithTripleAttr(body_triple, att);
  XMLToken_t *p_token = XMLToken_createWithTripleAttr(p_triple, att);
  XMLToken_t *text_token = XMLToken_createWithText("This is my text");
  XMLNode_t *html_node = XMLNode_createFromToken(html_token);
  XMLNode_t *head_node = XMLNode_createFromToken(head_token);
  XMLNode_t *title_node = XMLNode_createFromToken(title_token);
  XMLNode_t *body_node = XMLNode_createFromToken(body_token);
  XMLNode_t *p_node = XMLNode_createFromToken(p_token);
  XMLNode_t *text_node = XMLNode_createFromToken(text_token);

  XMLToken_t *p_token1 = XMLToken_createWithTripleAttrNS(p_triple, att, ns);
  XMLToken_t *text_token1 = XMLToken_createWithText("This is more text");
  XMLNode_t *p_node1 = XMLNode_createFromToken(p_token1);
  XMLNode_t *text_node1 = XMLNode_createFromToken(text_token1);

  XMLNode_t * notes;
  const XMLNode_t *child, *child1;

  XMLNode_addChild(p_node, text_node);
  XMLNode_addChild(body_node, p_node);
  XMLNode_addChild(head_node, title_node);
  XMLNode_addChild(html_node, head_node);
  XMLNode_addChild(html_node, body_node);

  XMLNode_addChild(p_node1, text_node1);

  SBase_setNotes(S, html_node);
  SBase_appendNotes(S, p_node1);

  notes = SBase_getNotes(S);

  fail_unless(!strcmp(XMLNode_getName(notes), "notes"));
  fail_unless(XMLNode_getNumChildren(notes) == 1);

  child = XMLNode_getChild(notes, 0);

  fail_unless(!strcmp(XMLNode_getName(child), "html"));
  fail_unless(XMLNode_getNumChildren(child) == 2);

  child = XMLNode_getChild(child, 1);

  fail_unless(!strcmp(XMLNode_getName(child), "body"));
  fail_unless(XMLNode_getNumChildren(child) == 2);

  child1 = XMLNode_getChild(child, 0);

  fail_unless(!strcmp(XMLNode_getName(child1), "p"));
  fail_unless(XMLNode_getNumChildren(child1) == 1);

  child1 = XMLNode_getChild(child1, 0);

  fail_unless(!strcmp(XMLNode_getCharacters(child1), "This is my text"));
  fail_unless(XMLNode_getNumChildren(child1) == 0);

  child1 = XMLNode_getChild(child, 1);

  fail_unless(!strcmp(XMLNode_getName(child1), "p"));
  fail_unless(XMLNode_getNumChildren(child1) == 1);

  child1 = XMLNode_getChild(child1, 0);

  fail_unless(!strcmp(XMLNode_getCharacters(child1), "This is more text"));
  fail_unless(XMLNode_getNumChildren(child1) == 0);

  XMLAttributes_free(att);
  XMLNamespaces_free(ns);
  XMLTriple_free(html_triple);
  XMLTriple_free(head_triple);
  XMLTriple_free(body_triple);
  XMLTriple_free(p_triple);
  XMLToken_free(html_token);
  XMLToken_free(head_token);
  XMLToken_free(body_token);
  XMLToken_free(p_token);
  XMLToken_free(text_token);
  XMLToken_free(text_token1);
  XMLToken_free(p_token1);
  XMLNode_free(html_node);
  XMLNode_free(head_node);
  XMLNode_free(body_node);
  XMLNode_free(p_node);
  XMLNode_free(text_node);
  XMLNode_free(p_node1);
  XMLNode_free(text_node1);
}
END_TEST


START_TEST (test_SBase_appendNotes4)
{
  // add a html tag to a body tag
  XMLAttributes_t * att = XMLAttributes_create ();
  XMLNamespaces_t *ns = XMLNamespaces_create();
  XMLNamespaces_add(ns, "http://www.w3.org/1999/xhtml", "");
  XMLTriple_t *html_triple = XMLTriple_createWith("html", "", "");
  XMLTriple_t *head_triple = XMLTriple_createWith("head", "", "");
  XMLTriple_t *title_triple = XMLTriple_createWith("title", "", "");
  XMLTriple_t *body_triple = XMLTriple_createWith("body", "", "");
  XMLTriple_t *p_triple = XMLTriple_createWith("p", "", "");
  XMLToken_t *html_token = XMLToken_createWithTripleAttrNS(html_triple, att, ns);
  XMLToken_t *head_token = XMLToken_createWithTripleAttr(head_triple, att);
  XMLToken_t *title_token = XMLToken_createWithTripleAttr(title_triple, att);
  XMLToken_t *body_token = XMLToken_createWithTripleAttr(body_triple, att);
  XMLToken_t *p_token = XMLToken_createWithTripleAttr(p_triple, att);
  XMLToken_t *body_token1 = XMLToken_createWithTripleAttrNS(body_triple, att, ns);
  XMLToken_t *text_token = XMLToken_createWithText("This is my text");
  XMLNode_t *body_node = XMLNode_createFromToken(body_token1);
  XMLNode_t *p_node = XMLNode_createFromToken(p_token);
  XMLNode_t *text_node = XMLNode_createFromToken(text_token);

  XMLToken_t *text_token1 = XMLToken_createWithText("This is more text");
  XMLNode_t *html_node1 = XMLNode_createFromToken(html_token);
  XMLNode_t *head_node1 = XMLNode_createFromToken(head_token);
  XMLNode_t *title_node1 = XMLNode_createFromToken(title_token);
  XMLNode_t *body_node1 = XMLNode_createFromToken(body_token);
  XMLNode_t *p_node1 = XMLNode_createFromToken(p_token);
  XMLNode_t *text_node1 = XMLNode_createFromToken(text_token1);

  XMLNode_t * notes;
  const XMLNode_t *child, *child1;

  XMLNode_addChild(p_node, text_node);
  XMLNode_addChild(body_node, p_node);

  XMLNode_addChild(p_node1, text_node1);
  XMLNode_addChild(body_node1, p_node1);
  XMLNode_addChild(head_node1, title_node1);
  XMLNode_addChild(html_node1, head_node1);
  XMLNode_addChild(html_node1, body_node1);

  SBase_setNotes(S, body_node);
  SBase_appendNotes(S, html_node1);

  notes = SBase_getNotes(S);

  fail_unless(!strcmp(XMLNode_getName(notes), "notes"));
  fail_unless(XMLNode_getNumChildren(notes) == 1);

  child = XMLNode_getChild(notes, 0);

  fail_unless(!strcmp(XMLNode_getName(child), "html"));
  fail_unless(XMLNode_getNumChildren(child) == 2);

  child = XMLNode_getChild(child, 1);

  fail_unless(!strcmp(XMLNode_getName(child), "body"));
  fail_unless(XMLNode_getNumChildren(child) == 2);

  child1 = XMLNode_getChild(child, 0);

  fail_unless(!strcmp(XMLNode_getName(child1), "p"));
  fail_unless(XMLNode_getNumChildren(child1) == 1);

  child1 = XMLNode_getChild(child1, 0);

  fail_unless(!strcmp(XMLNode_getCharacters(child1), "This is my text"));
  fail_unless(XMLNode_getNumChildren(child1) == 0);

  child1 = XMLNode_getChild(child, 1);

  fail_unless(!strcmp(XMLNode_getName(child1), "p"));
  fail_unless(XMLNode_getNumChildren(child1) == 1);

  child1 = XMLNode_getChild(child1, 0);

  fail_unless(!strcmp(XMLNode_getCharacters(child1), "This is more text"));
  fail_unless(XMLNode_getNumChildren(child1) == 0);

  XMLAttributes_free(att);
  XMLNamespaces_free(ns);
  XMLTriple_free(html_triple);
  XMLTriple_free(head_triple);
  XMLTriple_free(body_triple);
  XMLTriple_free(p_triple);
  XMLToken_free(body_token);
  XMLToken_free(p_token);
  XMLToken_free(text_token);
  XMLToken_free(text_token1);
  XMLToken_free(body_token1);
  XMLNode_free(body_node);
  XMLNode_free(p_node);
  XMLNode_free(text_node);
  XMLNode_free(html_node1);
  XMLNode_free(head_node1);
  XMLNode_free(body_node1);
  XMLNode_free(p_node1);
  XMLNode_free(text_node1);
}
END_TEST


START_TEST (test_SBase_appendNotes5)
{
  // add a html tag to a p tag
  XMLAttributes_t * att = XMLAttributes_create ();
  XMLNamespaces_t *ns = XMLNamespaces_create();
  XMLNamespaces_add(ns, "http://www.w3.org/1999/xhtml", "");
  XMLTriple_t *html_triple = XMLTriple_createWith("html", "", "");
  XMLTriple_t *head_triple = XMLTriple_createWith("head", "", "");
  XMLTriple_t *title_triple = XMLTriple_createWith("title", "", "");
  XMLTriple_t *body_triple = XMLTriple_createWith("body", "", "");
  XMLTriple_t *p_triple = XMLTriple_createWith("p", "", "");
  XMLToken_t *html_token = XMLToken_createWithTripleAttrNS(html_triple, att, ns);
  XMLToken_t *head_token = XMLToken_createWithTripleAttr(head_triple, att);
  XMLToken_t *title_token = XMLToken_createWithTripleAttr(title_triple, att);
  XMLToken_t *body_token = XMLToken_createWithTripleAttr(body_triple, att);
  XMLToken_t *p_token = XMLToken_createWithTripleAttr(p_triple, att);
  XMLToken_t *p_token1 = XMLToken_createWithTripleAttrNS(p_triple, att, ns);
  XMLToken_t *text_token = XMLToken_createWithText("This is my text");
  XMLNode_t *p_node = XMLNode_createFromToken(p_token1);
  XMLNode_t *text_node = XMLNode_createFromToken(text_token);

  XMLToken_t *text_token1 = XMLToken_createWithText("This is more text");
  XMLNode_t *html_node1 = XMLNode_createFromToken(html_token);
  XMLNode_t *head_node1 = XMLNode_createFromToken(head_token);
  XMLNode_t *title_node1 = XMLNode_createFromToken(title_token);
  XMLNode_t *body_node1 = XMLNode_createFromToken(body_token);
  XMLNode_t *p_node1 = XMLNode_createFromToken(p_token);
  XMLNode_t *text_node1 = XMLNode_createFromToken(text_token1);

  XMLNode_t * notes;
  const XMLNode_t *child, *child1;

  XMLNode_addChild(p_node, text_node);

  XMLNode_addChild(p_node1, text_node1);
  XMLNode_addChild(body_node1, p_node1);
  XMLNode_addChild(head_node1, title_node1);
  XMLNode_addChild(html_node1, head_node1);
  XMLNode_addChild(html_node1, body_node1);

  SBase_setNotes(S, p_node);
  SBase_appendNotes(S, html_node1);

  notes = SBase_getNotes(S);

  fail_unless(!strcmp(XMLNode_getName(notes), "notes"));
  fail_unless(XMLNode_getNumChildren(notes) == 1);

  child = XMLNode_getChild(notes, 0);

  fail_unless(!strcmp(XMLNode_getName(child), "html"));
  fail_unless(XMLNode_getNumChildren(child) == 2);

  child = XMLNode_getChild(child, 1);

  fail_unless(!strcmp(XMLNode_getName(child), "body"));
  fail_unless(XMLNode_getNumChildren(child) == 2);

  child1 = XMLNode_getChild(child, 0);

  fail_unless(!strcmp(XMLNode_getName(child1), "p"));
  fail_unless(XMLNode_getNumChildren(child1) == 1);

  child1 = XMLNode_getChild(child1, 0);

  fail_unless(!strcmp(XMLNode_getCharacters(child1), "This is my text"));
  fail_unless(XMLNode_getNumChildren(child1) == 0);

  child1 = XMLNode_getChild(child, 1);

  fail_unless(!strcmp(XMLNode_getName(child1), "p"));
  fail_unless(XMLNode_getNumChildren(child1) == 1);

  child1 = XMLNode_getChild(child1, 0);

  fail_unless(!strcmp(XMLNode_getCharacters(child1), "This is more text"));
  fail_unless(XMLNode_getNumChildren(child1) == 0);

  XMLAttributes_free(att);
  XMLNamespaces_free(ns);
  XMLTriple_free(html_triple);
  XMLTriple_free(head_triple);
  XMLTriple_free(body_triple);
  XMLTriple_free(p_triple);
  XMLToken_free(body_token);
  XMLToken_free(p_token);
  XMLToken_free(p_token1);
  XMLToken_free(text_token);
  XMLToken_free(text_token1);
  XMLNode_free(p_node);
  XMLNode_free(text_node);
  XMLNode_free(html_node1);
  XMLNode_free(head_node1);
  XMLNode_free(body_node1);
  XMLNode_free(p_node1);
  XMLNode_free(text_node1);
}
END_TEST


START_TEST (test_SBase_appendNotes6)
{// add a body tag to an body tag
  XMLAttributes_t * att = XMLAttributes_create ();
  XMLNamespaces_t *ns = XMLNamespaces_create();
  XMLNamespaces_add(ns, "http://www.w3.org/1999/xhtml", "");
  XMLTriple_t *body_triple = XMLTriple_createWith("body", "", "");
  XMLTriple_t *p_triple = XMLTriple_createWith("p", "", "");
  XMLToken_t *body_token = XMLToken_createWithTripleAttrNS(body_triple, att, ns);
  XMLToken_t *p_token = XMLToken_createWithTripleAttr(p_triple, att);
  XMLToken_t *text_token = XMLToken_createWithText("This is my text");
  XMLNode_t *body_node = XMLNode_createFromToken(body_token);
  XMLNode_t *p_node = XMLNode_createFromToken(p_token);
  XMLNode_t *text_node = XMLNode_createFromToken(text_token);

  XMLToken_t *text_token1 = XMLToken_createWithText("This is more text");
  XMLNode_t *body_node1 = XMLNode_createFromToken(body_token);
  XMLNode_t *p_node1 = XMLNode_createFromToken(p_token);
  XMLNode_t *text_node1 = XMLNode_createFromToken(text_token1);

  XMLNode_t * notes;
  const XMLNode_t *child, *child1;

  XMLNode_addChild(p_node, text_node);
  XMLNode_addChild(body_node, p_node);

  XMLNode_addChild(p_node1, text_node1);
  XMLNode_addChild(body_node1, p_node1);

  SBase_setNotes(S, body_node);
  SBase_appendNotes(S, body_node1);

  notes = SBase_getNotes(S);

  fail_unless(!strcmp(XMLNode_getName(notes), "notes"));
  fail_unless(XMLNode_getNumChildren(notes) == 1);

  child = XMLNode_getChild(notes, 0);

  fail_unless(!strcmp(XMLNode_getName(child), "body"));
  fail_unless(XMLNode_getNumChildren(child) == 2);

  child1 = XMLNode_getChild(child, 0);

  fail_unless(!strcmp(XMLNode_getName(child1), "p"));
  fail_unless(XMLNode_getNumChildren(child1) == 1);

  child1 = XMLNode_getChild(child1, 0);

  fail_unless(!strcmp(XMLNode_getCharacters(child1), "This is my text"));
  fail_unless(XMLNode_getNumChildren(child1) == 0);

  child1 = XMLNode_getChild(child, 1);

  fail_unless(!strcmp(XMLNode_getName(child1), "p"));
  fail_unless(XMLNode_getNumChildren(child1) == 1);

  child1 = XMLNode_getChild(child1, 0);

  fail_unless(!strcmp(XMLNode_getCharacters(child1), "This is more text"));
  fail_unless(XMLNode_getNumChildren(child1) == 0);

  XMLAttributes_free(att);
  XMLNamespaces_free(ns);
  XMLTriple_free(body_triple);
  XMLTriple_free(p_triple);
  XMLToken_free(body_token);
  XMLToken_free(p_token);
  XMLToken_free(text_token);
  XMLToken_free(text_token1);
  XMLNode_free(body_node);
  XMLNode_free(p_node);
  XMLNode_free(text_node);
  XMLNode_free(body_node1);
  XMLNode_free(p_node1);
  XMLNode_free(text_node1);
}
END_TEST


START_TEST (test_SBase_appendNotes7)
{// add a body tag to an p tag
  XMLAttributes_t * att = XMLAttributes_create ();
  XMLNamespaces_t *ns = XMLNamespaces_create();
  XMLNamespaces_add(ns, "http://www.w3.org/1999/xhtml", "");
  XMLTriple_t *body_triple = XMLTriple_createWith("body", "", "");
  XMLTriple_t *p_triple = XMLTriple_createWith("p", "", "");
  XMLToken_t *body_token = XMLToken_createWithTripleAttrNS(body_triple, att, ns);
  XMLToken_t *p_token1 = XMLToken_createWithTripleAttrNS(p_triple, att, ns);
  XMLToken_t *text_token = XMLToken_createWithText("This is my text");
  XMLToken_t *p_token = XMLToken_createWithTripleAttr(p_triple, att);
  XMLNode_t *p_node = XMLNode_createFromToken(p_token1);
  XMLNode_t *text_node = XMLNode_createFromToken(text_token);

  XMLToken_t *text_token1 = XMLToken_createWithText("This is more text");
  XMLNode_t *body_node1 = XMLNode_createFromToken(body_token);
  XMLNode_t *p_node1 = XMLNode_createFromToken(p_token);
  XMLNode_t *text_node1 = XMLNode_createFromToken(text_token1);

  XMLNode_t * notes;
  const XMLNode_t *child, *child1;

  XMLNode_addChild(p_node, text_node);

  XMLNode_addChild(p_node1, text_node1);
  XMLNode_addChild(body_node1, p_node1);

  SBase_setNotes(S, p_node);
  SBase_appendNotes(S, body_node1);

  notes = SBase_getNotes(S);

  fail_unless(!strcmp(XMLNode_getName(notes), "notes"));
  fail_unless(XMLNode_getNumChildren(notes) == 1);

  child = XMLNode_getChild(notes, 0);

  fail_unless(!strcmp(XMLNode_getName(child), "body"));
  fail_unless(XMLNode_getNumChildren(child) == 2);

  child1 = XMLNode_getChild(child, 0);

  fail_unless(!strcmp(XMLNode_getName(child1), "p"));
  fail_unless(XMLNode_getNumChildren(child1) == 1);

  child1 = XMLNode_getChild(child1, 0);

  fail_unless(!strcmp(XMLNode_getCharacters(child1), "This is my text"));
  fail_unless(XMLNode_getNumChildren(child1) == 0);

  child1 = XMLNode_getChild(child, 1);

  fail_unless(!strcmp(XMLNode_getName(child1), "p"));
  fail_unless(XMLNode_getNumChildren(child1) == 1);

  child1 = XMLNode_getChild(child1, 0);

  fail_unless(!strcmp(XMLNode_getCharacters(child1), "This is more text"));
  fail_unless(XMLNode_getNumChildren(child1) == 0);

  XMLAttributes_free(att);
  XMLNamespaces_free(ns);
  XMLTriple_free(body_triple);
  XMLTriple_free(p_triple);
  XMLToken_free(body_token);
  XMLToken_free(p_token);
  XMLToken_free(p_token1);
  XMLToken_free(text_token);
  XMLToken_free(text_token1);
  XMLNode_free(p_node);
  XMLNode_free(text_node);
  XMLNode_free(body_node1);
  XMLNode_free(p_node1);
  XMLNode_free(text_node1);
}
END_TEST


START_TEST (test_SBase_appendNotes8)
{
  // add a p tag to an body tag
  XMLAttributes_t * att = XMLAttributes_create ();
  XMLNamespaces_t *ns = XMLNamespaces_create();
  XMLNamespaces_add(ns, "http://www.w3.org/1999/xhtml", "");
  XMLTriple_t *body_triple = XMLTriple_createWith("body", "", "");
  XMLTriple_t *p_triple = XMLTriple_createWith("p", "", "");
  XMLToken_t *body_token = XMLToken_createWithTripleAttrNS(body_triple, att, ns);
  XMLToken_t *p_token = XMLToken_createWithTripleAttr(p_triple, att);
  XMLToken_t *text_token = XMLToken_createWithText("This is my text");
  XMLNode_t *body_node = XMLNode_createFromToken(body_token);
  XMLNode_t *p_node = XMLNode_createFromToken(p_token);
  XMLNode_t *text_node = XMLNode_createFromToken(text_token);

  XMLToken_t *p_token1 = XMLToken_createWithTripleAttrNS(p_triple, att, ns);
  XMLToken_t *text_token1 = XMLToken_createWithText("This is more text");
  XMLNode_t *p_node1 = XMLNode_createFromToken(p_token1);
  XMLNode_t *text_node1 = XMLNode_createFromToken(text_token1);

  XMLNode_t * notes;
  const XMLNode_t *child, *child1;

  XMLNode_addChild(p_node, text_node);
  XMLNode_addChild(body_node, p_node);

  XMLNode_addChild(p_node1, text_node1);

  SBase_setNotes(S, body_node);
  SBase_appendNotes(S, p_node1);

  notes = SBase_getNotes(S);

  fail_unless(!strcmp(XMLNode_getName(notes), "notes"));
  fail_unless(XMLNode_getNumChildren(notes) == 1);

  child = XMLNode_getChild(notes, 0);

  fail_unless(!strcmp(XMLNode_getName(child), "body"));
  fail_unless(XMLNode_getNumChildren(child) == 2);

  child1 = XMLNode_getChild(child, 0);

  fail_unless(!strcmp(XMLNode_getName(child1), "p"));
  fail_unless(XMLNode_getNumChildren(child1) == 1);

  child1 = XMLNode_getChild(child1, 0);

  fail_unless(!strcmp(XMLNode_getCharacters(child1), "This is my text"));
  fail_unless(XMLNode_getNumChildren(child1) == 0);

  child1 = XMLNode_getChild(child, 1);

  fail_unless(!strcmp(XMLNode_getName(child1), "p"));
  fail_unless(XMLNode_getNumChildren(child1) == 1);

  child1 = XMLNode_getChild(child1, 0);

  fail_unless(!strcmp(XMLNode_getCharacters(child1), "This is more text"));
  fail_unless(XMLNode_getNumChildren(child1) == 0);

  XMLAttributes_free(att);
  XMLNamespaces_free(ns);
  XMLTriple_free(body_triple);
  XMLTriple_free(p_triple);
  XMLToken_free(body_token);
  XMLToken_free(p_token);
  XMLToken_free(text_token);
  XMLToken_free(text_token1);
  XMLToken_free(p_token1);
  XMLNode_free(body_node);
  XMLNode_free(p_node);
  XMLNode_free(text_node);
  XMLNode_free(p_node1);
  XMLNode_free(text_node1);
}
END_TEST


START_TEST (test_SBase_appendNotesString)
{
  // add p to p 
  char * notes = "<p xmlns=\"http://www.w3.org/1999/xhtml\">This is a test note </p>";
  char * taggednewnotes = "<notes>\n"
                       "  <p xmlns=\"http://www.w3.org/1999/xhtml\">This is a test note </p>\n"
                       "  <p xmlns=\"http://www.w3.org/1999/xhtml\">This is more test notes </p>\n"
                       "</notes>";
  char * taggednewnotes2 = "<notes>\n"
                       "  <p xmlns=\"http://www.w3.org/1999/xhtml\">This is a test note </p>\n"
                       "  <p xmlns=\"http://www.w3.org/1999/xhtml\">This is more test notes 1</p>\n"
                       "  <p xmlns=\"http://www.w3.org/1999/xhtml\">This is more test notes 2</p>\n"
                       "</notes>";
  char * newnotes = "<p xmlns=\"http://www.w3.org/1999/xhtml\">This is more test notes </p>";
  char * newnotes2 = "<p xmlns=\"http://www.w3.org/1999/xhtml\">This is more test notes 1</p>"
                     "<p xmlns=\"http://www.w3.org/1999/xhtml\">This is more test notes 2</p>";
  char * newnotes3= "<notes>\n"
                    "  <p xmlns=\"http://www.w3.org/1999/xhtml\">This is more test notes </p>\n"
                    "</notes>";
  char * newnotes4 = "<notes>\n"
                     "  <p xmlns=\"http://www.w3.org/1999/xhtml\">This is more test notes 1</p>\n"
                     "  <p xmlns=\"http://www.w3.org/1999/xhtml\">This is more test notes 2</p>\n"
                     "</notes>";

  SBase_setNotesString(S, notes);

  fail_unless(SBase_isSetNotes(S) == 1);

  //
  // add one p tag
  //
  SBase_appendNotesString(S, newnotes);

  const char * notes1 = SBase_getNotesString(S);

  fail_unless(SBase_isSetNotes(S) == 1);
  fail_unless(!strcmp(taggednewnotes, notes1));

  // 
  // add two p tags
  //
  SBase_setNotesString(S, notes);
  SBase_appendNotesString(S, newnotes2);

  const char * notes2 = SBase_getNotesString(S);

  fail_unless(SBase_isSetNotes(S) == 1);
  fail_unless(!strcmp(taggednewnotes2, notes2));

  //
  // add one p tag with notes tag
  //
  SBase_setNotesString(S, notes);
  SBase_appendNotesString(S, newnotes3);

  const char * notes3 = SBase_getNotesString(S);

  fail_unless(SBase_isSetNotes(S) == 1);
  fail_unless(!strcmp(taggednewnotes, notes3));

  //
  // add two p tags with notes tag
  //
  SBase_setNotesString(S, notes);
  SBase_appendNotesString(S, newnotes4);

  const char * notes4 = SBase_getNotesString(S);

  fail_unless(SBase_isSetNotes(S) == 1);
  fail_unless(!strcmp(taggednewnotes2, notes4));
}
END_TEST


START_TEST (test_SBase_appendNotesString1)
{ // add html to html
  char * notes = "<html xmlns=\"http://www.w3.org/1999/xhtml\">\n"
                 "  <head>\n"
                 "    <title/>\n"
                 "  </head>\n"
                 "  <body>\n"
                 "    <p>This is a test note </p>\n"
                 "  </body>\n"
                 "</html>";
  char * taggednewnotes = 
                 "<notes>\n"
                 "  <html xmlns=\"http://www.w3.org/1999/xhtml\">\n"
                 "    <head>\n"
                 "      <title/>\n"
                 "    </head>\n"
                 "    <body>\n"
                 "      <p>This is a test note </p>\n"
                 "      <p>This is more test notes </p>\n"
                 "    </body>\n"
                 "  </html>\n"
                 "</notes>";
  char * addnotes = "<html xmlns=\"http://www.w3.org/1999/xhtml\">\n"
                 "  <head>\n"
                 "    <title/>\n"
                 "  </head>\n"
                 "  <body>\n"
                 "    <p>This is more test notes </p>\n"
                 "  </body>\n"
                 "</html>";
  char * addnotes2 =
                 "<notes>\n"
                 "  <html xmlns=\"http://www.w3.org/1999/xhtml\">\n"
                 "    <head>\n"
                 "      <title/>\n"
                 "    </head>\n"
                 "    <body>\n"
                 "      <p>This is more test notes </p>\n"
                 "    </body>\n"
                 "  </html>\n"
                 "</notes>";

  // add html tag

  SBase_setNotesString(S, notes);
  SBase_appendNotesString(S, addnotes);

  const char *notes1 = SBase_getNotesString(S);

  fail_unless(SBase_isSetNotes(S) == 1);
  fail_unless(!strcmp(taggednewnotes, notes1));

  // add html tag with notes tag

  SBase_setNotesString(S, notes);
  SBase_appendNotesString(S, addnotes2);

  const char *notes2 = SBase_getNotesString(S);

  fail_unless(SBase_isSetNotes(S) == 1);
  fail_unless(!strcmp(taggednewnotes, notes2));

}
END_TEST


START_TEST (test_SBase_appendNotesString2)
{ // add body to html
  char * notes = "<html xmlns=\"http://www.w3.org/1999/xhtml\">\n"
                 "  <head>\n"
                 "    <title/>\n"
                 "  </head>\n"
                 "  <body>\n"
                 "    <p>This is a test note </p>\n"
                 "  </body>\n"
                 "</html>";
  char * taggednewnotes = 
                 "<notes>\n"
                 "  <html xmlns=\"http://www.w3.org/1999/xhtml\">\n"
                 "    <head>\n"
                 "      <title/>\n"
                 "    </head>\n"
                 "    <body>\n"
                 "      <p>This is a test note </p>\n"
                 "      <p>This is more test notes </p>\n"
                 "    </body>\n"
                 "  </html>\n"
                 "</notes>";
  char * addnotes = "<body xmlns=\"http://www.w3.org/1999/xhtml\">\n"
                    "  <p>This is more test notes </p>\n"
                    "</body>\n";
  char * addnotes2 =
                    "<notes>\n"
                    "  <body xmlns=\"http://www.w3.org/1999/xhtml\">\n"
                    "    <p>This is more test notes </p>\n"
                    "  </body>\n"
                    "</notes>";

  // add body tag 

  SBase_setNotesString(S, notes);
  SBase_appendNotesString(S, addnotes);

  const char *notes1 = SBase_getNotesString(S);

  fail_unless(SBase_isSetNotes(S) == 1);
  fail_unless(!strcmp(taggednewnotes, notes1));

  // add body tag with notes tag

  SBase_setNotesString(S, notes);
  SBase_appendNotesString(S, addnotes2);

  const char *notes2 = SBase_getNotesString(S);

  fail_unless(SBase_isSetNotes(S) == 1);
  fail_unless(!strcmp(taggednewnotes, notes2));
}
END_TEST


START_TEST (test_SBase_appendNotesString3)
{ // add p to html
  char * notes = "<html xmlns=\"http://www.w3.org/1999/xhtml\">\n"
                 "  <head>\n"
                 "    <title/>\n"
                 "  </head>\n"
                 "  <body>\n"
                 "    <p>This is a test note </p>\n"
                 "  </body>\n"
                 "</html>";
  char * taggednewnotes = 
                 "<notes>\n"
                 "  <html xmlns=\"http://www.w3.org/1999/xhtml\">\n"
                 "    <head>\n"
                 "      <title/>\n"
                 "    </head>\n"
                 "    <body>\n"
                 "      <p>This is a test note </p>\n"
                 "      <p xmlns=\"http://www.w3.org/1999/xhtml\">This is more test notes </p>\n"
                 "    </body>\n"
                 "  </html>\n"
                 "</notes>";
  char * taggednewnotes2 =
                 "<notes>\n"
                 "  <html xmlns=\"http://www.w3.org/1999/xhtml\">\n"
                 "    <head>\n"
                 "      <title/>\n"
                 "    </head>\n"
                 "    <body>\n"
                 "      <p>This is a test note </p>\n"
                 "      <p xmlns=\"http://www.w3.org/1999/xhtml\">This is more test notes 1</p>\n"
                 "      <p xmlns=\"http://www.w3.org/1999/xhtml\">This is more test notes 2</p>\n"
                 "    </body>\n"
                 "  </html>\n"
                 "</notes>";
  char * addnotes = "<p xmlns=\"http://www.w3.org/1999/xhtml\">This is more test notes </p>\n";
  char * addnotes2 = "<p xmlns=\"http://www.w3.org/1999/xhtml\">This is more test notes 1</p>\n"
                     "<p xmlns=\"http://www.w3.org/1999/xhtml\">This is more test notes 2</p>";
  char * addnotes3 = "<notes>\n"
                     "  <p xmlns=\"http://www.w3.org/1999/xhtml\">This is more test notes </p>\n"
                     "</notes>";
  char * addnotes4 = "<notes>\n"
                     "  <p xmlns=\"http://www.w3.org/1999/xhtml\">This is more test notes 1</p>\n"
                     "  <p xmlns=\"http://www.w3.org/1999/xhtml\">This is more test notes 2</p>\n"
                     "</notes>";

  // add one p tag

  SBase_setNotesString(S, notes);
  SBase_appendNotesString(S, addnotes);

  const char *notes1 = SBase_getNotesString(S);

  fail_unless(SBase_isSetNotes(S) == 1);
  fail_unless(!strcmp(taggednewnotes, notes1));

  // add two p tags

  SBase_setNotesString(S, notes);
  SBase_appendNotesString(S, addnotes2);

  const char *notes2 = SBase_getNotesString(S);

  fail_unless(SBase_isSetNotes(S) == 1);
  fail_unless(!strcmp(taggednewnotes2, notes2));

  // add one p tag with notes tag

  SBase_setNotesString(S, notes);
  SBase_appendNotesString(S, addnotes3);

  const char *notes3 = SBase_getNotesString(S);

  fail_unless(SBase_isSetNotes(S) == 1);
  fail_unless(!strcmp(taggednewnotes, notes3));

  // add two p tags with notes tags

  SBase_setNotesString(S, notes);
  SBase_appendNotesString(S, addnotes4);

  const char *notes4 = SBase_getNotesString(S);

  fail_unless(SBase_isSetNotes(S) == 1);
  fail_unless(!strcmp(taggednewnotes2, notes4));
}
END_TEST


START_TEST (test_SBase_appendNotesString4)
{ // add html to body
  char * notes = "<body xmlns=\"http://www.w3.org/1999/xhtml\">\n"
                 "  <p>This is a test note </p>\n"
                 "</body>";
  char * taggednewnotes = 
                 "<notes>\n"
                 "  <html xmlns=\"http://www.w3.org/1999/xhtml\">\n"
                 "    <head>\n"
                 "      <title/>\n"
                 "    </head>\n"
                 "    <body>\n"
                 "      <p>This is a test note </p>\n"
                 "      <p>This is more test notes </p>\n"
                 "    </body>\n"
                 "  </html>\n"
                 "</notes>";
  char * addnotes = "<html xmlns=\"http://www.w3.org/1999/xhtml\">\n"
                 "  <head>\n"
                 "    <title/>\n"
                 "  </head>\n"
                 "  <body>\n"
                 "    <p>This is more test notes </p>\n"
                 "  </body>\n"
                 "</html>";
  char * addnotes2 =
                 "<notes>\n"
                 "  <html xmlns=\"http://www.w3.org/1999/xhtml\">\n"
                 "    <head>\n"
                 "      <title/>\n"
                 "    </head>\n"
                 "    <body>\n"
                 "      <p>This is more test notes </p>\n"
                 "    </body>\n"
                 "  </html>\n"
                 "</notes>";

  // add html tag

  SBase_setNotesString(S, notes);
  SBase_appendNotesString(S, addnotes);

  const char *notes1 = SBase_getNotesString(S);

  fail_unless(SBase_isSetNotes(S) == 1);
  fail_unless(!strcmp(taggednewnotes, notes1));

  // add html tag with notes tag

  SBase_setNotesString(S, notes);
  SBase_appendNotesString(S, addnotes2);

  const char *notes2 = SBase_getNotesString(S);

  fail_unless(SBase_isSetNotes(S) == 1);
  fail_unless(!strcmp(taggednewnotes, notes2));
}
END_TEST


START_TEST (test_SBase_appendNotesString5)
{ // add html to p
  char * notes = "<p xmlns=\"http://www.w3.org/1999/xhtml\">This is a test note </p>";
  char * taggednewnotes = 
                 "<notes>\n"
                 "  <html xmlns=\"http://www.w3.org/1999/xhtml\">\n"
                 "    <head>\n"
                 "      <title/>\n"
                 "    </head>\n"
                 "    <body>\n"
                 "      <p xmlns=\"http://www.w3.org/1999/xhtml\">This is a test note </p>\n"
                 "      <p>This is more test notes </p>\n"
                 "    </body>\n"
                 "  </html>\n"
                 "</notes>";
  char * addnotes = "<html xmlns=\"http://www.w3.org/1999/xhtml\">\n"
                 "  <head>\n"
                 "    <title/>\n"
                 "  </head>\n"
                 "  <body>\n"
                 "    <p>This is more test notes </p>\n"
                 "  </body>\n"
                 "</html>";
  char * addnotes2 = 
                 "<notes>\n"
                 "  <html xmlns=\"http://www.w3.org/1999/xhtml\">\n"
                 "    <head>\n"
                 "      <title/>\n"
                 "    </head>\n"
                 "    <body>\n"
                 "      <p>This is more test notes </p>\n"
                 "    </body>\n"
                 "  </html>\n"
                 "</notes>";

  // add html tag

  SBase_setNotesString(S, notes);
  SBase_appendNotesString(S, addnotes);

  const char *notes1 = SBase_getNotesString(S);

  fail_unless(SBase_isSetNotes(S) == 1);
  fail_unless(!strcmp(taggednewnotes, notes1));

  // add html tag with notes tag

  SBase_setNotesString(S, notes);
  SBase_appendNotesString(S, addnotes2);

  const char *notes2 = SBase_getNotesString(S);

  fail_unless(SBase_isSetNotes(S) == 1);
  fail_unless(!strcmp(taggednewnotes, notes2));

}
END_TEST


START_TEST (test_SBase_appendNotesString6)
{ // add body to body
  char * notes = "<body xmlns=\"http://www.w3.org/1999/xhtml\">\n"
                 "  <p>This is a test note </p>\n"
                 "</body>";
  char * taggednewnotes = 
                 "<notes>\n"
                 "  <body xmlns=\"http://www.w3.org/1999/xhtml\">\n"
                 "    <p>This is a test note </p>\n"
                 "    <p>This is more test notes </p>\n"
                 "  </body>\n"
                 "</notes>";
  char * addnotes = "<body xmlns=\"http://www.w3.org/1999/xhtml\">\n"
                 "  <p>This is more test notes </p>\n"
                 "</body>";
  char * addnotes2 = 
                 "<notes>\n"
                 "  <body xmlns=\"http://www.w3.org/1999/xhtml\">\n"
                 "    <p>This is more test notes </p>\n"
                 "  </body>\n"
                 "</notes>";

  // add body tag

  SBase_setNotesString(S, notes);
  SBase_appendNotesString(S, addnotes);

  const char *notes1 = SBase_getNotesString(S);

  fail_unless(SBase_isSetNotes(S) == 1);
  fail_unless(!strcmp(taggednewnotes, notes1));

  // add body tag with notes tag

  SBase_setNotesString(S, notes);
  SBase_appendNotesString(S, addnotes2);

  const char *notes2 = SBase_getNotesString(S);

  fail_unless(SBase_isSetNotes(S) == 1);
  fail_unless(!strcmp(taggednewnotes, notes2));

}
END_TEST


START_TEST (test_SBase_appendNotesString7)
{ // add body to p
  char * notes = "<p xmlns=\"http://www.w3.org/1999/xhtml\">This is a test note </p>";
  char * taggednewnotes = 
                 "<notes>\n"
                 "  <body xmlns=\"http://www.w3.org/1999/xhtml\">\n"
                 "    <p xmlns=\"http://www.w3.org/1999/xhtml\">This is a test note </p>\n"
                 "    <p>This is more test notes </p>\n"
                 "  </body>\n"
                 "</notes>";
  char * addnotes = "<body xmlns=\"http://www.w3.org/1999/xhtml\">\n"
                 "  <p>This is more test notes </p>\n"
                 "</body>";
  char * addnotes2 = 
                 "<notes>\n"
                 "  <body xmlns=\"http://www.w3.org/1999/xhtml\">\n"
                 "    <p>This is more test notes </p>\n"
                 "  </body>\n"
                 "</notes>";

  // add body tag

  SBase_setNotesString(S, notes);
  SBase_appendNotesString(S, addnotes);

  const char *notes1 = SBase_getNotesString(S);

  fail_unless(SBase_isSetNotes(S) == 1);
  fail_unless(!strcmp(taggednewnotes, notes1));

  // add body tag with notes tag

  SBase_setNotesString(S, notes);
  SBase_appendNotesString(S, addnotes2);

  const char *notes2 = SBase_getNotesString(S);

  fail_unless(SBase_isSetNotes(S) == 1);
  fail_unless(!strcmp(taggednewnotes, notes2));

}
END_TEST


START_TEST (test_SBase_appendNotesString8)
{ // add p to body
  char * notes = "<body xmlns=\"http://www.w3.org/1999/xhtml\">\n"
                 "  <p>This is a test note </p>\n"
                 "</body>";
  char * taggednewnotes = 
                 "<notes>\n"
                 "  <body xmlns=\"http://www.w3.org/1999/xhtml\">\n"
                 "    <p>This is a test note </p>\n"
                 "    <p xmlns=\"http://www.w3.org/1999/xhtml\">This is more test notes </p>\n"
                 "  </body>\n"
                 "</notes>";
  char * taggednewnotes2 = 
                 "<notes>\n"
                 "  <body xmlns=\"http://www.w3.org/1999/xhtml\">\n"
                 "    <p>This is a test note </p>\n"
                 "    <p xmlns=\"http://www.w3.org/1999/xhtml\">This is more test notes 1</p>\n"
                 "    <p xmlns=\"http://www.w3.org/1999/xhtml\">This is more test notes 2</p>\n"
                 "  </body>\n"
                 "</notes>";
  char * addnotes = "<p xmlns=\"http://www.w3.org/1999/xhtml\">This is more test notes </p>";
  char * addnotes2 = "<p xmlns=\"http://www.w3.org/1999/xhtml\">This is more test notes 1</p>\n"
                     "<p xmlns=\"http://www.w3.org/1999/xhtml\">This is more test notes 2</p>";
  char * addnotes3 = 
                 "<notes>\n"
                 "  <p xmlns=\"http://www.w3.org/1999/xhtml\">This is more test notes </p>\n"
                 "</notes>";
  char * addnotes4 = 
                 "<notes>\n"
                 "  <p xmlns=\"http://www.w3.org/1999/xhtml\">This is more test notes 1</p>\n"
                 "  <p xmlns=\"http://www.w3.org/1999/xhtml\">This is more test notes 2</p>\n"
                 "</notes>";

  // add one p tag

  SBase_setNotesString(S, notes);
  SBase_appendNotesString(S, addnotes);

  const char *notes1 = SBase_getNotesString(S);

  fail_unless(SBase_isSetNotes(S) == 1);
  fail_unless(!strcmp(taggednewnotes, notes1));

  // add two p tags

  SBase_setNotesString(S, notes);
  SBase_appendNotesString(S, addnotes2);

  const char *notes2 = SBase_getNotesString(S);

  fail_unless(SBase_isSetNotes(S) == 1);
  fail_unless(!strcmp(taggednewnotes2, notes2));

  // add one p tag with notes tag

  SBase_setNotesString(S, notes);
  SBase_appendNotesString(S, addnotes3);

  const char *notes3 = SBase_getNotesString(S);

  fail_unless(SBase_isSetNotes(S) == 1);
  fail_unless(!strcmp(taggednewnotes, notes3));

  // add two p tags with notes tags

  SBase_setNotesString(S, notes);
  SBase_appendNotesString(S, addnotes4);

  const char *notes4 = SBase_getNotesString(S);

  fail_unless(SBase_isSetNotes(S) == 1);
  fail_unless(!strcmp(taggednewnotes2, notes4));

}
END_TEST


START_TEST(test_SBase_CVTerms)
{
  CVTerm_t * cv = CVTerm_createWithQualifierType(BIOLOGICAL_QUALIFIER);
  CVTerm_setBiologicalQualifierType(cv, BQB_IS);
  CVTerm_addResource(cv, "foo");
  
  fail_unless(SBase_getNumCVTerms(S) == 0);
  fail_unless(SBase_getCVTerms(S) == NULL);

  SBase_setMetaId(S, "_id");
  SBase_addCVTerm(S, cv);
  fail_unless(SBase_getNumCVTerms(S) == 1);
  fail_unless(SBase_getCVTerms(S) != NULL);

  fail_unless(SBase_getCVTerm(S, 0) != cv);

  CVTerm_free(cv);


}
END_TEST


START_TEST(test_SBase_addCVTerms)
{
  CVTerm_t * cv = CVTerm_createWithQualifierType(BIOLOGICAL_QUALIFIER);
  CVTerm_setBiologicalQualifierType(cv, BQB_ENCODES);
  CVTerm_addResource(cv, "foo");
  SBase_setMetaId(S, "sbase1");
  
  SBase_addCVTerm(S, cv);
  
  fail_unless(SBase_getNumCVTerms(S) == 1);
  fail_unless(SBase_getCVTerms(S) != NULL);

  XMLAttributes_t *res = CVTerm_getResources(SBase_getCVTerm(S, 0));
  fail_unless(!strcmp(res->getValue(0).c_str(), "foo"));

  CVTerm_t * cv1 = CVTerm_createWithQualifierType(BIOLOGICAL_QUALIFIER);
  CVTerm_setBiologicalQualifierType(cv1, BQB_IS);
  CVTerm_addResource(cv1, "bar");
  
  SBase_addCVTerm(S, cv1);
  
  fail_unless(SBase_getNumCVTerms(S) == 2);

  /* same qualifier so should just add to resources of existing term */
  CVTerm_t * cv2 = CVTerm_createWithQualifierType(BIOLOGICAL_QUALIFIER);
  CVTerm_setBiologicalQualifierType(cv2, BQB_IS);
  CVTerm_addResource(cv2, "bar1");
  
  SBase_addCVTerm(S, cv2);
  
  fail_unless(SBase_getNumCVTerms(S) == 2);
  
  res = CVTerm_getResources(SBase_getCVTerm(S, 1));

  fail_unless(XMLAttributes_getLength(res) == 2);
  fail_unless(!strcmp(res->getValue(0).c_str(), "bar"));
  fail_unless(!strcmp(res->getValue(1).c_str(), "bar1"));


  /* existing term shouldnt get added*/
  CVTerm_t * cv4 = CVTerm_createWithQualifierType(BIOLOGICAL_QUALIFIER);
  CVTerm_setBiologicalQualifierType(cv4, BQB_IS);
  CVTerm_addResource(cv4, "bar1");
  
  SBase_addCVTerm(S, cv4);
  
  fail_unless(SBase_getNumCVTerms(S) == 2);
  
  res = CVTerm_getResources(SBase_getCVTerm(S, 1));

  fail_unless(XMLAttributes_getLength(res) == 2);
  fail_unless(!strcmp(res->getValue(0).c_str(), "bar"));
  fail_unless(!strcmp(res->getValue(1).c_str(), "bar1"));
  
  /* existing term with different qualifier shouldnt get added*/
  CVTerm_t * cv5 = CVTerm_createWithQualifierType(BIOLOGICAL_QUALIFIER);
  CVTerm_setBiologicalQualifierType(cv5, BQB_HAS_PART);
  CVTerm_addResource(cv5, "bar1");
  
  SBase_addCVTerm(S, cv5);
  
  fail_unless(SBase_getNumCVTerms(S) == 2);
  
  res = CVTerm_getResources(SBase_getCVTerm(S, 1));

  fail_unless(XMLAttributes_getLength(res) == 2);
  fail_unless(!strcmp(res->getValue(0).c_str(), "bar"));
  fail_unless(!strcmp(res->getValue(1).c_str(), "bar1"));
 
  CVTerm_free(cv);
  CVTerm_free(cv2);
  CVTerm_free(cv1);
  CVTerm_free(cv4);


}
END_TEST


START_TEST(test_SBase_unsetCVTerms)
{
  CVTerm_t * cv = CVTerm_createWithQualifierType(BIOLOGICAL_QUALIFIER);
  CVTerm_setBiologicalQualifierType(cv, BQB_ENCODES);
  CVTerm_addResource(cv, "foo");
  
  SBase_setMetaId(S, "sbase1");
  SBase_addCVTerm(S, cv);
  CVTerm_t * cv1 = CVTerm_createWithQualifierType(BIOLOGICAL_QUALIFIER);
  CVTerm_setBiologicalQualifierType(cv1, BQB_IS);
  CVTerm_addResource(cv1, "bar");
  
  SBase_addCVTerm(S, cv1);
  CVTerm_t * cv2 = CVTerm_createWithQualifierType(BIOLOGICAL_QUALIFIER);
  CVTerm_setBiologicalQualifierType(cv2, BQB_IS);
  CVTerm_addResource(cv2, "bar1");
  
  SBase_addCVTerm(S, cv2);
  CVTerm_t * cv4 = CVTerm_createWithQualifierType(BIOLOGICAL_QUALIFIER);
  CVTerm_setBiologicalQualifierType(cv4, BQB_IS);
  CVTerm_addResource(cv4, "bar1");
  
  SBase_addCVTerm(S, cv4);
  
  fail_unless(SBase_getNumCVTerms(S) == 2);

  SBase_unsetCVTerms(S);

  fail_unless(SBase_getNumCVTerms(S) == 0);
  fail_unless(SBase_getCVTerms(S) == NULL);
  
  CVTerm_free(cv);
  CVTerm_free(cv2);
  CVTerm_free(cv1);
  CVTerm_free(cv4);
}
END_TEST


START_TEST(test_SBase_getQualifiersFromResources)
{
  CVTerm_t * cv = CVTerm_createWithQualifierType(BIOLOGICAL_QUALIFIER);
  CVTerm_setBiologicalQualifierType(cv, BQB_ENCODES);
  CVTerm_addResource(cv, "foo");
  
  SBase_setMetaId(S, "sbase1");
  SBase_addCVTerm(S, cv);

  fail_unless(SBase_getResourceBiologicalQualifier(S, "foo") == BQB_ENCODES);
  
  CVTerm_t * cv1 = CVTerm_createWithQualifierType(MODEL_QUALIFIER);
  CVTerm_setModelQualifierType(cv1, BQM_IS);
  CVTerm_addResource(cv1, "bar");
  
  SBase_addCVTerm(S, cv1);

  fail_unless(SBase_getResourceModelQualifier(S, "bar") == BQM_IS);
  
  CVTerm_free(cv);
  CVTerm_free(cv1);


}
END_TEST


Suite *
create_suite_SBase (void)
{
  Suite *suite = suite_create("SBase");
  TCase *tcase = tcase_create("SBase");


  tcase_add_checked_fixture(tcase, SBaseTest_setup, SBaseTest_teardown);

  tcase_add_test(tcase, test_SBase_setMetaId     );
 // tcase_add_test(tcase, test_SBase_setNotes      );
  tcase_add_test(tcase, test_SBase_setAnnotation );
  tcase_add_test(tcase, test_SBase_setNotesString);
  tcase_add_test(tcase, test_SBase_setAnnotationString);
  tcase_add_test(tcase, test_SBase_unsetAnnotationWithCVTerms );
  tcase_add_test(tcase, test_SBase_unsetAnnotationWithModelHistory );

  tcase_add_test(tcase, test_SBase_appendNotes );
  tcase_add_test(tcase, test_SBase_appendNotes1 );
  tcase_add_test(tcase, test_SBase_appendNotes2 );
  tcase_add_test(tcase, test_SBase_appendNotes3 );
  tcase_add_test(tcase, test_SBase_appendNotes4 );
  tcase_add_test(tcase, test_SBase_appendNotes5 );
  tcase_add_test(tcase, test_SBase_appendNotes6 );
  tcase_add_test(tcase, test_SBase_appendNotes7 );
  tcase_add_test(tcase, test_SBase_appendNotes8 );
  tcase_add_test(tcase, test_SBase_appendNotesString );
  tcase_add_test(tcase, test_SBase_appendNotesString1);
  tcase_add_test(tcase, test_SBase_appendNotesString2);
  tcase_add_test(tcase, test_SBase_appendNotesString3);
  tcase_add_test(tcase, test_SBase_appendNotesString4);
  tcase_add_test(tcase, test_SBase_appendNotesString5);
  tcase_add_test(tcase, test_SBase_appendNotesString6);
  tcase_add_test(tcase, test_SBase_appendNotesString7);
  tcase_add_test(tcase, test_SBase_appendNotesString8);

  tcase_add_test(tcase, test_SBase_CVTerms );
  tcase_add_test(tcase, test_SBase_addCVTerms );
  tcase_add_test(tcase, test_SBase_unsetCVTerms );
  tcase_add_test(tcase, test_SBase_getQualifiersFromResources );

  suite_add_tcase(suite, tcase);

  return suite;
}


END_C_DECLS
