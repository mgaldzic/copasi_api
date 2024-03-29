# @file    TestRDFAnnotation2.rb
# @brief   fomula units data unit tests
#
# @author  Akiya Jouraku (Ruby conversion)
# @author  Ben Bornstein 
#
# $Id$
# $HeadURL$
#
# ====== WARNING ===== WARNING ===== WARNING ===== WARNING ===== WARNING ======
#
# DO NOT EDIT THIS FILE.
#
# This file was generated automatically by converting the file located at
# src/annotation/test/TestRDFAnnotation2.cpp
# using the conversion program dev/utilities/translateTests/translateTests.pl.
# Any changes made here will be lost the next time the file is regenerated.
#
# -----------------------------------------------------------------------------
# This file is part of libSBML.  Please visit http://sbml.org for more
# information about SBML, and the latest version of libSBML.
#
# Copyright 2005-2010 California Institute of Technology.
# Copyright 2002-2005 California Institute of Technology and
#                     Japan Science and Technology Corporation.
# 
# This library is free software; you can redistribute it and/or modify it
# under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation.  A copy of the license agreement is provided
# in the file named "LICENSE.txt" included with this software distribution
# and also available online as http://sbml.org/software/libsbml/license.html
# -----------------------------------------------------------------------------
require 'test/unit'
require 'libSBML'

class TestRDFAnnotation2 < Test::Unit::TestCase

  def util_NaN
    z = 0.0
    return 0.0/z
  end

  def util_PosInf
    z = 0.0
    return 1.0/z
  end

  def util_NegInf
    z = 0.0
    return -1.0/z
  end

  def equals(*x)
    case x.size
    when 2
      e, s = x
      return e == s
    when 1
      e, = x
      return e == @@oss.str()
    end
  end

  def setup
    filename = "../../annotation/test/test-data/annotation2.xml"
    @@d2 = LibSBML::readSBML(filename)
    @@m2 = @@d2.getModel()
  end

  def teardown
    @@d2 = nil
  end

  def test_RDFAnnotation2_getModelHistory
    history = @@m2.getModelHistory()
    assert( history != nil )
    mc = history.getCreator(0)
    assert ((  "Hucka" == mc.getFamilyName() ))
    assert ((  "Mike" == mc.getGivenName() ))
    assert ((  "mhucka@caltech.edu" == mc.getEmail() ))
    assert ((  "BNMC" == mc.getOrganisation() ))
    mc1 = history.getCreator(1)
    assert ((  "Keating" == mc1.getFamilyName() ))
    assert ((  "Sarah" == mc1.getGivenName() ))
    assert ((  "skeating@caltech.edu" == mc1.getEmail() ))
    assert ((  "UH" == mc1.getOrganisation() ))
    date = history.getCreatedDate()
    assert( date.getYear() == 2005 )
    assert( date.getMonth() == 2 )
    assert( date.getDay() == 2 )
    assert( date.getHour() == 14 )
    assert( date.getMinute() == 56 )
    assert( date.getSecond() == 11 )
    assert( date.getSignOffset() == 0 )
    assert( date.getHoursOffset() == 0 )
    assert( date.getMinutesOffset() == 0 )
    assert ((  "2005-02-02T14:56:11Z" == date.getDateAsString() ))
    date = history.getModifiedDate()
    assert( date.getYear() == 2006 )
    assert( date.getMonth() == 5 )
    assert( date.getDay() == 30 )
    assert( date.getHour() == 10 )
    assert( date.getMinute() == 46 )
    assert( date.getSecond() == 2 )
    assert( date.getSignOffset() == 0 )
    assert( date.getHoursOffset() == 0 )
    assert( date.getMinutesOffset() == 0 )
    assert ((  "2006-05-30T10:46:02Z" == date.getDateAsString() ))
    date = history.getModifiedDate(1)
    assert( date.getYear() == 2007 )
    assert( date.getMonth() == 1 )
    assert( date.getDay() == 16 )
    assert( date.getHour() == 15 )
    assert( date.getMinute() == 31 )
    assert( date.getSecond() == 52 )
    assert( date.getSignOffset() == 0 )
    assert( date.getHoursOffset() == 0 )
    assert( date.getMinutesOffset() == 0 )
    assert ((  "2007-01-16T15:31:52Z" == date.getDateAsString() ))
  end

  def test_RDFAnnotation2_modelWithHistoryAndCVTerms
    h = LibSBML::ModelHistory.new()
    c = LibSBML::ModelCreator.new()
    c.setFamilyName("Keating")
    c.setGivenName("Sarah")
    h.addCreator(c)
    d = LibSBML::Date.new(2008,11,17,18,37,0,0,0,0)
    h.setCreatedDate(d)
    h.setModifiedDate(d)
    @@m2.unsetModelHistory()
    @@m2.setModelHistory(h)
    cv = LibSBML::CVTerm.new()
    cv.setQualifierType(LibSBML::BIOLOGICAL_QUALIFIER)
    cv.setBiologicalQualifierType(LibSBML::BQB_IS_VERSION_OF)
    cv.addResource("http://www.geneontology.org/#GO:0005892")
    @@m2.addCVTerm(cv)
    ann = LibSBML::RDFAnnotationParser.parseModelHistory(@@m2)
    expected = "<annotation>\n" + 
    "  <rdf:RDF xmlns:rdf=\"http://www.w3.org/1999/02/22-rdf-syntax-ns#\" xmlns:dc=\"http://purl.org/dc/elements/1.1/\" xmlns:dcterms=\"http://purl.org/dc/terms/\" xmlns:vCard=\"http://www.w3.org/2001/vcard-rdf/3.0#\" xmlns:bqbiol=\"http://biomodels.net/biology-qualifiers/\" xmlns:bqmodel=\"http://biomodels.net/model-qualifiers/\">\n" + 
    "    <rdf:Description rdf:about=\"#_000001\">\n" + 
    "      <dc:creator rdf:parseType=\"Resource\">\n" + 
    "        <rdf:Bag>\n" + 
    "          <rdf:li rdf:parseType=\"Resource\">\n" + 
    "            <vCard:N rdf:parseType=\"Resource\">\n" + 
    "              <vCard:Family>Keating</vCard:Family>\n" + 
    "              <vCard:Given>Sarah</vCard:Given>\n" + 
    "            </vCard:N>\n" + 
    "          </rdf:li>\n" + 
    "        </rdf:Bag>\n" + 
    "      </dc:creator>\n" + 
    "      <dcterms:created rdf:parseType=\"Resource\">\n" + 
    "        <dcterms:W3CDTF>2008-11-17T18:37:00Z</dcterms:W3CDTF>\n" + 
    "      </dcterms:created>\n" + 
    "      <dcterms:modified rdf:parseType=\"Resource\">\n" + 
    "        <dcterms:W3CDTF>2008-11-17T18:37:00Z</dcterms:W3CDTF>\n" + 
    "      </dcterms:modified>\n" + 
    "      <bqbiol:isVersionOf>\n" + 
    "        <rdf:Bag>\n" + 
    "          <rdf:li rdf:resource=\"http://www.geneontology.org/#GO:0005892\"/>\n" + 
    "        </rdf:Bag>\n" + 
    "      </bqbiol:isVersionOf>\n" + 
    "    </rdf:Description>\n" + 
    "  </rdf:RDF>\n" + 
    "</annotation>"
    if (ann != nil)
      assert_equal true, equals(expected,ann.toXMLString())
    end
    end
  end

  def test_RDFAnnotation2_modelWithHistoryAndMultipleModifiedDates
    h = LibSBML::ModelHistory.new()
    c = LibSBML::ModelCreator.new()
    c.setFamilyName("Keating")
    c.setGivenName("Sarah")
    h.addCreator(c)
    d = LibSBML::Date.new(2005,2,2,14,56,11)
    h.setCreatedDate(d)
    h.addModifiedDate(d)
    h.addModifiedDate(d)
    @@m2.unsetModelHistory()
    @@m2.setModelHistory(h)
    ann = LibSBML::RDFAnnotationParser.parseModelHistory(@@m2)
    expected = "<annotation>\n" + 
    "  <rdf:RDF xmlns:rdf=\"http://www.w3.org/1999/02/22-rdf-syntax-ns#\" xmlns:dc=\"http://purl.org/dc/elements/1.1/\" xmlns:dcterms=\"http://purl.org/dc/terms/\" xmlns:vCard=\"http://www.w3.org/2001/vcard-rdf/3.0#\" xmlns:bqbiol=\"http://biomodels.net/biology-qualifiers/\" xmlns:bqmodel=\"http://biomodels.net/model-qualifiers/\">\n" + 
    "    <rdf:Description rdf:about=\"#_000001\">\n" + 
    "      <dc:creator rdf:parseType=\"Resource\">\n" + 
    "        <rdf:Bag>\n" + 
    "          <rdf:li rdf:parseType=\"Resource\">\n" + 
    "            <vCard:N rdf:parseType=\"Resource\">\n" + 
    "              <vCard:Family>Keating</vCard:Family>\n" + 
    "              <vCard:Given>Sarah</vCard:Given>\n" + 
    "            </vCard:N>\n" + 
    "          </rdf:li>\n" + 
    "        </rdf:Bag>\n" + 
    "      </dc:creator>\n" + 
    "      <dcterms:created rdf:parseType=\"Resource\">\n" + 
    "        <dcterms:W3CDTF>2005-02-02T14:56:11Z</dcterms:W3CDTF>\n" + 
    "      </dcterms:created>\n" + 
    "      <dcterms:modified rdf:parseType=\"Resource\">\n" + 
    "        <dcterms:W3CDTF>2005-02-02T14:56:11Z</dcterms:W3CDTF>\n" + 
    "      </dcterms:modified>\n" + 
    "      <dcterms:modified rdf:parseType=\"Resource\">\n" + 
    "        <dcterms:W3CDTF>2005-02-02T14:56:11Z</dcterms:W3CDTF>\n" + 
    "      </dcterms:modified>\n" + 
    "    </rdf:Description>\n" + 
    "  </rdf:RDF>\n" + 
    "</annotation>"
    assert_equal true, equals(expected,ann.toXMLString())
  end

  def test_RDFAnnotation2_modelWithHistoryWithCharacterReference
    h = LibSBML::ModelHistory.new()
    c = LibSBML::ModelCreator.new()
    c.setFamilyName("Dr&#228;ger")
    c.setGivenName("Andreas")
    h.addCreator(c)
    d = LibSBML::Date.new(2005,2,2,14,56,11)
    h.setCreatedDate(d)
    h.addModifiedDate(d)
    @@m2.unsetModelHistory()
    @@m2.setModelHistory(h)
    ann = LibSBML::RDFAnnotationParser.parseModelHistory(@@m2)
    expected = "<annotation>\n" + 
    "  <rdf:RDF xmlns:rdf=\"http://www.w3.org/1999/02/22-rdf-syntax-ns#\" xmlns:dc=\"http://purl.org/dc/elements/1.1/\" xmlns:dcterms=\"http://purl.org/dc/terms/\" xmlns:vCard=\"http://www.w3.org/2001/vcard-rdf/3.0#\" xmlns:bqbiol=\"http://biomodels.net/biology-qualifiers/\" xmlns:bqmodel=\"http://biomodels.net/model-qualifiers/\">\n" + 
    "    <rdf:Description rdf:about=\"#_000001\">\n" + 
    "      <dc:creator rdf:parseType=\"Resource\">\n" + 
    "        <rdf:Bag>\n" + 
    "          <rdf:li rdf:parseType=\"Resource\">\n" + 
    "            <vCard:N rdf:parseType=\"Resource\">\n" + 
    "              <vCard:Family>Dr&#228;ger</vCard:Family>\n" + 
    "              <vCard:Given>Andreas</vCard:Given>\n" + 
    "            </vCard:N>\n" + 
    "          </rdf:li>\n" + 
    "        </rdf:Bag>\n" + 
    "      </dc:creator>\n" + 
    "      <dcterms:created rdf:parseType=\"Resource\">\n" + 
    "        <dcterms:W3CDTF>2005-02-02T14:56:11Z</dcterms:W3CDTF>\n" + 
    "      </dcterms:created>\n" + 
    "      <dcterms:modified rdf:parseType=\"Resource\">\n" + 
    "        <dcterms:W3CDTF>2005-02-02T14:56:11Z</dcterms:W3CDTF>\n" + 
    "      </dcterms:modified>\n" + 
    "    </rdf:Description>\n" + 
    "  </rdf:RDF>\n" + 
    "</annotation>"
    assert_equal true, equals(expected,ann.toXMLString())
  end

$end
