# @file    TestXMLToken.rb
# @brief   XMLToken unit tests
#
# @author  Akiya Jouraku (Ruby conversion)
# @author  Michael Hucka <mhucka@caltech.edu> 
#
# $Id$
# $HeadURL$
#
# ====== WARNING ===== WARNING ===== WARNING ===== WARNING ===== WARNING ======
#
# DO NOT EDIT THIS FILE.
#
# This file was generated automatically by converting the file located at
# src/xml/test/TestXMLToken.c
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

class TestXMLToken < Test::Unit::TestCase

  def test_XMLToken_attribute_add_remove
    triple = LibSBML::XMLTriple.new("test","","")
    attr = LibSBML::XMLAttributes.new()
    token = LibSBML::XMLToken.new(triple,attr)
    xt1 = LibSBML::XMLTriple.new("name1", "http://name1.org/", "p1")
    xt2 = LibSBML::XMLTriple.new("name2", "http://name2.org/", "p2")
    xt3 = LibSBML::XMLTriple.new("name3", "http://name3.org/", "p3")
    xt1a = LibSBML::XMLTriple.new("name1", "http://name1a.org/", "p1a")
    xt2a = LibSBML::XMLTriple.new("name2", "http://name2a.org/", "p2a")
    token.addAttr( "name1", "val1", "http://name1.org/", "p1")
    token.addAttr(xt2, "val2")
    assert( token.getAttributesLength() == 2 )
    assert( token.isAttributesEmpty() == false )
    assert( (  "name1" != token.getAttrName(0) ) == false )
    assert( (  "val1"  != token.getAttrValue(0) ) == false )
    assert( (  "http://name1.org/" != token.getAttrURI(0) ) == false )
    assert( (  "p1"    != token.getAttrPrefix(0) ) == false )
    assert( (  "name2" != token.getAttrName(1) ) == false )
    assert( (  "val2"  != token.getAttrValue(1) ) == false )
    assert( (  "http://name2.org/" != token.getAttrURI(1) ) == false )
    assert( (  "p2"    != token.getAttrPrefix(1) ) == false )
    assert( token.getAttrValue( "name1") == "" )
    assert( token.getAttrValue( "name2") == "" )
    assert( (  "val1"  != token.getAttrValue( "name1", "http://name1.org/") ) == false )
    assert( (  "val2"  != token.getAttrValue( "name2", "http://name2.org/") ) == false )
    assert( (  "val1"  != token.getAttrValue(xt1) ) == false )
    assert( (  "val2"  != token.getAttrValue(xt2) ) == false )
    assert( token.hasAttr(-1) == false )
    assert( token.hasAttr(2) == false )
    assert( token.hasAttr(0) == true )
    assert( token.hasAttr( "name1", "http://name1.org/") == true )
    assert( token.hasAttr( "name2", "http://name2.org/") == true )
    assert( token.hasAttr( "name3", "http://name3.org/") == false )
    assert( token.hasAttr(xt1) == true )
    assert( token.hasAttr(xt2) == true )
    assert( token.hasAttr(xt3) == false )
    token.addAttr( "noprefix", "val3")
    assert( token.getAttributesLength() == 3 )
    assert( token.isAttributesEmpty() == false )
    assert( (  "noprefix" != token.getAttrName(2) ) == false )
    assert( (  "val3"     != token.getAttrValue(2) ) == false )
    assert( token.getAttrURI(2) == "" )
    assert( token.getAttrPrefix(2) == "" )
    assert( (      "val3"  != token.getAttrValue( "noprefix") ) == false )
    assert( (  "val3"  != token.getAttrValue( "noprefix", "") ) == false )
    assert( token.hasAttr( "noprefix"    ) == true )
    assert( token.hasAttr( "noprefix", "") == true )
    token.addAttr(xt1, "mval1")
    token.addAttr( "name2", "mval2", "http://name2.org/", "p2")
    assert( token.getAttributesLength() == 3 )
    assert( token.isAttributesEmpty() == false )
    assert( (  "name1" != token.getAttrName(0) ) == false )
    assert( (  "mval1" != token.getAttrValue(0) ) == false )
    assert( (  "http://name1.org/" != token.getAttrURI(0) ) == false )
    assert( (  "p1"    != token.getAttrPrefix(0) ) == false )
    assert( (  "name2"    != token.getAttrName(1) ) == false )
    assert( (  "mval2"    != token.getAttrValue(1) ) == false )
    assert( (  "http://name2.org/" != token.getAttrURI(1) ) == false )
    assert( (  "p2"       != token.getAttrPrefix(1) ) == false )
    assert( token.hasAttr(xt1) == true )
    assert( token.hasAttr( "name1", "http://name1.org/") == true )
    token.addAttr( "noprefix", "mval3")
    assert( token.getAttributesLength() == 3 )
    assert( token.isAttributesEmpty() == false )
    assert( (  "noprefix" != token.getAttrName(2) ) == false )
    assert( (  "mval3"    != token.getAttrValue(2) ) == false )
    assert( token.getAttrURI(2) == "" )
    assert( token.getAttrPrefix(2) == "" )
    assert( token.hasAttr( "noprefix") == true )
    assert( token.hasAttr( "noprefix", "") == true )
    token.addAttr(xt1a, "val1a")
    token.addAttr(xt2a, "val2a")
    assert( token.getAttributesLength() == 5 )
    assert( (  "name1" != token.getAttrName(3) ) == false )
    assert( (  "val1a" != token.getAttrValue(3) ) == false )
    assert( (  "http://name1a.org/" != token.getAttrURI(3) ) == false )
    assert( (  "p1a" != token.getAttrPrefix(3) ) == false )
    assert( (  "name2" != token.getAttrName(4) ) == false )
    assert( (  "val2a" != token.getAttrValue(4) ) == false )
    assert( (  "http://name2a.org/" != token.getAttrURI(4) ) == false )
    assert( (  "p2a" != token.getAttrPrefix(4) ) == false )
    assert( (  "val1a"  != token.getAttrValue( "name1", "http://name1a.org/") ) == false )
    assert( (  "val2a"  != token.getAttrValue( "name2", "http://name2a.org/") ) == false )
    assert( (  "val1a"  != token.getAttrValue(xt1a) ) == false )
    assert( (  "val2a"  != token.getAttrValue(xt2a) ) == false )
    token.removeAttr(xt1a)
    token.removeAttr(xt2a)
    assert( token.getAttributesLength() == 3 )
    token.removeAttr( "name1", "http://name1.org/")
    assert( token.getAttributesLength() == 2 )
    assert( token.isAttributesEmpty() == false )
    assert( (  "name2" != token.getAttrName(0) ) == false )
    assert( (  "mval2" != token.getAttrValue(0) ) == false )
    assert( (  "http://name2.org/" != token.getAttrURI(0) ) == false )
    assert( (  "p2" != token.getAttrPrefix(0) ) == false )
    assert( (  "noprefix" != token.getAttrName(1) ) == false )
    assert( (  "mval3" != token.getAttrValue(1) ) == false )
    assert( token.getAttrURI(1) == "" )
    assert( token.getAttrPrefix(1) == "" )
    assert( token.hasAttr( "name1", "http://name1.org/") == false )
    token.removeAttr(xt2)
    assert( token.getAttributesLength() == 1 )
    assert( token.isAttributesEmpty() == false )
    assert( (  "noprefix" != token.getAttrName(0) ) == false )
    assert( (  "mval3" != token.getAttrValue(0) ) == false )
    assert( token.getAttrURI(0) == "" )
    assert( token.getAttrPrefix(0) == "" )
    assert( token.hasAttr(xt2) == false )
    assert( token.hasAttr( "name2", "http://name2.org/") == false )
    token.removeAttr( "noprefix")
    assert( token.getAttributesLength() == 0 )
    assert( token.isAttributesEmpty() == true )
    assert( token.hasAttr( "noprefix"    ) == false )
    assert( token.hasAttr( "noprefix", "") == false )
    token = nil
    xt1 = nil
    xt2 = nil
    xt3 = nil
    xt1a = nil
    xt2a = nil
    triple = nil
    attr = nil
  end

  def test_XMLToken_attribute_set_clear
    triple = LibSBML::XMLTriple.new("test","","")
    attr = LibSBML::XMLAttributes.new()
    token = LibSBML::XMLToken.new(triple,attr)
    nattr = LibSBML::XMLAttributes.new()
    xt1 = LibSBML::XMLTriple.new("name1", "http://name1.org/", "p1")
    xt2 = LibSBML::XMLTriple.new("name2", "http://name2.org/", "p2")
    xt3 = LibSBML::XMLTriple.new("name3", "http://name3.org/", "p3")
    xt4 = LibSBML::XMLTriple.new("name4", "http://name4.org/", "p4")
    xt5 = LibSBML::XMLTriple.new("name5", "http://name5.org/", "p5")
    nattr.add(xt1, "val1")
    nattr.add(xt2, "val2")
    nattr.add(xt3, "val3")
    nattr.add(xt4, "val4")
    nattr.add(xt5, "val5")
    token.setAttributes(nattr)
    assert( token.getAttributesLength() == 5 )
    assert( token.isAttributesEmpty() == false )
    assert( (  "name1" != token.getAttrName(0) ) == false )
    assert( (  "val1"  != token.getAttrValue(0) ) == false )
    assert( (  "http://name1.org/" != token.getAttrURI(0) ) == false )
    assert( (  "p1"    != token.getAttrPrefix(0) ) == false )
    assert( (  "name2" != token.getAttrName(1) ) == false )
    assert( (  "val2"  != token.getAttrValue(1) ) == false )
    assert( (  "http://name2.org/" != token.getAttrURI(1) ) == false )
    assert( (  "p2"    != token.getAttrPrefix(1) ) == false )
    assert( (  "name3" != token.getAttrName(2) ) == false )
    assert( (  "val3"  != token.getAttrValue(2) ) == false )
    assert( (  "http://name3.org/" != token.getAttrURI(2) ) == false )
    assert( (  "p3"    != token.getAttrPrefix(2) ) == false )
    assert( (  "name4" != token.getAttrName(3) ) == false )
    assert( (  "val4"  != token.getAttrValue(3) ) == false )
    assert( (  "http://name4.org/" != token.getAttrURI(3) ) == false )
    assert( (  "p4"    != token.getAttrPrefix(3) ) == false )
    assert( (  "name5" != token.getAttrName(4) ) == false )
    assert( (  "val5"  != token.getAttrValue(4) ) == false )
    assert( (  "http://name5.org/" != token.getAttrURI(4) ) == false )
    assert( (  "p5"    != token.getAttrPrefix(4) ) == false )
    ntriple = LibSBML::XMLTriple.new("test2","http://test2.org/","p2")
    token.setTriple(ntriple)
    assert( (    "test2" != token.getName() ) == false )
    assert( (     "http://test2.org/" != token.getURI() ) == false )
    assert( (  "p2" != token.getPrefix() ) == false )
    token.clearAttributes()
    assert( token.getAttributesLength() == 0 )
    assert( token.isAttributesEmpty() != false )
    nattr = nil
    triple = nil
    ntriple = nil
    attr = nil
    token = nil
    xt1 = nil
    xt2 = nil
    xt3 = nil
    xt4 = nil
    xt5 = nil
  end

  def test_XMLToken_chars
    token = LibSBML::XMLToken.new("This is text")
    assert( token.isElement() == false )
    assert( token.isEnd() == false )
    assert( token.isStart() == false )
    assert( token.isText() == true )
    assert( token.isEOF() == false )
    assert( (  "This is text" != token.getCharacters() ) == false )
    token = nil
  end

  def test_XMLToken_create
    token = LibSBML::XMLToken.new()
    assert( token != nil )
    token = nil
    triple = LibSBML::XMLTriple.new("attr", "uri", "prefix")
    token = LibSBML::XMLToken.new(triple)
    assert( token != nil )
    assert( (  "attr" != token.getName() ) == false )
    assert( (  "prefix" != token.getPrefix() ) == false )
    assert( (  "uri" != token.getURI() ) == false )
    token = nil
    attr = LibSBML::XMLAttributes.new()
    assert( attr != nil )
    attr.add( "attr2", "value")
    token = LibSBML::XMLToken.new(triple,attr)
    assert( token != nil )
    returnattr = token.getAttributes()
    assert( (  "attr2" != returnattr.getName(0) ) == false )
    token = nil
    triple = nil
    attr = nil
  end

  def test_XMLToken_fields
    triple = LibSBML::XMLTriple.new("attr", "uri", "prefix")
    token = LibSBML::XMLToken.new(triple)
    assert( token.isElement() == true )
    assert( token.isEnd() == true )
    assert( token.isStart() == false )
    assert( token.isText() == false )
    assert( token.isEOF() == false )
    assert( (  "attr" != token.getName() ) == false )
    assert( (  "uri" != token.getURI() ) == false )
    assert( (  "prefix" != token.getPrefix() ) == false )
    token = nil
    triple = nil
  end

  def test_XMLToken_namespace_add
    triple = LibSBML::XMLTriple.new("test","","")
    attr = LibSBML::XMLAttributes.new()
    token = LibSBML::XMLToken.new(triple,attr)
    assert( token.getNamespacesLength() == 0 )
    assert( token.isNamespacesEmpty() == true )
    token.addNamespace( "http://test1.org/", "test1")
    assert( token.getNamespacesLength() == 1 )
    assert( token.isNamespacesEmpty() == false )
    token.addNamespace( "http://test2.org/", "test2")
    assert( token.getNamespacesLength() == 2 )
    assert( token.isNamespacesEmpty() == false )
    token.addNamespace( "http://test1.org/", "test1a")
    assert( token.getNamespacesLength() == 3 )
    assert( token.isNamespacesEmpty() == false )
    token.addNamespace( "http://test1.org/", "test1a")
    assert( token.getNamespacesLength() == 3 )
    assert( token.isNamespacesEmpty() == false )
    assert( !( token.getNamespaceIndex( "http://test1.org/") == -1) )
    token = nil
    triple = nil
    attr = nil
  end

  def test_XMLToken_namespace_get
    triple = LibSBML::XMLTriple.new("test","","")
    attr = LibSBML::XMLAttributes.new()
    token = LibSBML::XMLToken.new(triple,attr)
    token.addNamespace( "http://test1.org/", "test1")
    token.addNamespace( "http://test2.org/", "test2")
    token.addNamespace( "http://test3.org/", "test3")
    token.addNamespace( "http://test4.org/", "test4")
    token.addNamespace( "http://test5.org/", "test5")
    token.addNamespace( "http://test6.org/", "test6")
    token.addNamespace( "http://test7.org/", "test7")
    token.addNamespace( "http://test8.org/", "test8")
    token.addNamespace( "http://test9.org/", "test9")
    assert( token.getNamespacesLength() == 9 )
    assert( token.getNamespaceIndex( "http://test1.org/") == 0 )
    assert( (  "test2" != token.getNamespacePrefix(1) ) == false )
    assert( ( 		      "test1" != token.getNamespacePrefix( "http://test1.org/") ) == false )
    assert( (  "http://test2.org/" != token.getNamespaceURI(1) ) == false )
    assert( ( 		      "http://test2.org/" != token.getNamespaceURI( "test2") ) == false )
    assert( token.getNamespaceIndex( "http://test1.org/") == 0 )
    assert( token.getNamespaceIndex( "http://test2.org/") == 1 )
    assert( token.getNamespaceIndex( "http://test5.org/") == 4 )
    assert( token.getNamespaceIndex( "http://test9.org/") == 8 )
    assert( token.getNamespaceIndex( "http://testX.org/") == -1 )
    assert( token.hasNamespaceURI( "http://test1.org/") != false )
    assert( token.hasNamespaceURI( "http://test2.org/") != false )
    assert( token.hasNamespaceURI( "http://test5.org/") != false )
    assert( token.hasNamespaceURI( "http://test9.org/") != false )
    assert( token.hasNamespaceURI( "http://testX.org/") == false )
    assert( token.getNamespaceIndexByPrefix( "test1") == 0 )
    assert( token.getNamespaceIndexByPrefix( "test5") == 4 )
    assert( token.getNamespaceIndexByPrefix( "test9") == 8 )
    assert( token.getNamespaceIndexByPrefix( "testX") == -1 )
    assert( token.hasNamespacePrefix( "test1") != false )
    assert( token.hasNamespacePrefix( "test5") != false )
    assert( token.hasNamespacePrefix( "test9") != false )
    assert( token.hasNamespacePrefix( "testX") == false )
    assert( token.hasNamespaceNS( "http://test1.org/", "test1") != false )
    assert( token.hasNamespaceNS( "http://test5.org/", "test5") != false )
    assert( token.hasNamespaceNS( "http://test9.org/", "test9") != false )
    assert( token.hasNamespaceNS( "http://testX.org/", "testX") == false )
    token = nil
    triple = nil
    attr = nil
  end

  def test_XMLToken_namespace_remove
    triple = LibSBML::XMLTriple.new("test","","")
    attr = LibSBML::XMLAttributes.new()
    token = LibSBML::XMLToken.new(triple,attr)
    token.addNamespace( "http://test1.org/", "test1")
    token.addNamespace( "http://test2.org/", "test2")
    token.addNamespace( "http://test3.org/", "test3")
    token.addNamespace( "http://test4.org/", "test4")
    token.addNamespace( "http://test5.org/", "test5")
    assert( token.getNamespacesLength() == 5 )
    token.removeNamespace(4)
    assert( token.getNamespacesLength() == 4 )
    token.removeNamespace(3)
    assert( token.getNamespacesLength() == 3 )
    token.removeNamespace(2)
    assert( token.getNamespacesLength() == 2 )
    token.removeNamespace(1)
    assert( token.getNamespacesLength() == 1 )
    token.removeNamespace(0)
    assert( token.getNamespacesLength() == 0 )
    token.addNamespace( "http://test1.org/", "test1")
    token.addNamespace( "http://test2.org/", "test2")
    token.addNamespace( "http://test3.org/", "test3")
    token.addNamespace( "http://test4.org/", "test4")
    token.addNamespace( "http://test5.org/", "test5")
    assert( token.getNamespacesLength() == 5 )
    token.removeNamespace(0)
    assert( token.getNamespacesLength() == 4 )
    token.removeNamespace(0)
    assert( token.getNamespacesLength() == 3 )
    token.removeNamespace(0)
    assert( token.getNamespacesLength() == 2 )
    token.removeNamespace(0)
    assert( token.getNamespacesLength() == 1 )
    token.removeNamespace(0)
    assert( token.getNamespacesLength() == 0 )
    token = nil
    triple = nil
    attr = nil
  end

  def test_XMLToken_namespace_remove_by_prefix
    triple = LibSBML::XMLTriple.new("test","","")
    attr = LibSBML::XMLAttributes.new()
    token = LibSBML::XMLToken.new(triple,attr)
    token.addNamespace( "http://test1.org/", "test1")
    token.addNamespace( "http://test2.org/", "test2")
    token.addNamespace( "http://test3.org/", "test3")
    token.addNamespace( "http://test4.org/", "test4")
    token.addNamespace( "http://test5.org/", "test5")
    assert( token.getNamespacesLength() == 5 )
    token.removeNamespace( "test1")
    assert( token.getNamespacesLength() == 4 )
    token.removeNamespace( "test2")
    assert( token.getNamespacesLength() == 3 )
    token.removeNamespace( "test3")
    assert( token.getNamespacesLength() == 2 )
    token.removeNamespace( "test4")
    assert( token.getNamespacesLength() == 1 )
    token.removeNamespace( "test5")
    assert( token.getNamespacesLength() == 0 )
    token.addNamespace( "http://test1.org/", "test1")
    token.addNamespace( "http://test2.org/", "test2")
    token.addNamespace( "http://test3.org/", "test3")
    token.addNamespace( "http://test4.org/", "test4")
    token.addNamespace( "http://test5.org/", "test5")
    assert( token.getNamespacesLength() == 5 )
    token.removeNamespace( "test5")
    assert( token.getNamespacesLength() == 4 )
    token.removeNamespace( "test4")
    assert( token.getNamespacesLength() == 3 )
    token.removeNamespace( "test3")
    assert( token.getNamespacesLength() == 2 )
    token.removeNamespace( "test2")
    assert( token.getNamespacesLength() == 1 )
    token.removeNamespace( "test1")
    assert( token.getNamespacesLength() == 0 )
    token.addNamespace( "http://test1.org/", "test1")
    token.addNamespace( "http://test2.org/", "test2")
    token.addNamespace( "http://test3.org/", "test3")
    token.addNamespace( "http://test4.org/", "test4")
    token.addNamespace( "http://test5.org/", "test5")
    assert( token.getNamespacesLength() == 5 )
    token.removeNamespace( "test3")
    assert( token.getNamespacesLength() == 4 )
    token.removeNamespace( "test1")
    assert( token.getNamespacesLength() == 3 )
    token.removeNamespace( "test4")
    assert( token.getNamespacesLength() == 2 )
    token.removeNamespace( "test5")
    assert( token.getNamespacesLength() == 1 )
    token.removeNamespace( "test2")
    assert( token.getNamespacesLength() == 0 )
    token = nil
    triple = nil
    attr = nil
  end

  def test_XMLToken_namespace_set_clear
    triple = LibSBML::XMLTriple.new("test","","")
    attr = LibSBML::XMLAttributes.new()
    token = LibSBML::XMLToken.new(triple,attr)
    ns = LibSBML::XMLNamespaces.new()
    assert( token.getNamespacesLength() == 0 )
    assert( token.isNamespacesEmpty() == true )
    ns.add( "http://test1.org/", "test1")
    ns.add( "http://test2.org/", "test2")
    ns.add( "http://test3.org/", "test3")
    ns.add( "http://test4.org/", "test4")
    ns.add( "http://test5.org/", "test5")
    token.setNamespaces(ns)
    assert( token.getNamespacesLength() == 5 )
    assert( token.isNamespacesEmpty() == false )
    assert( (  "test1" != token.getNamespacePrefix(0) ) == false )
    assert( (  "test2" != token.getNamespacePrefix(1) ) == false )
    assert( (  "test3" != token.getNamespacePrefix(2) ) == false )
    assert( (  "test4" != token.getNamespacePrefix(3) ) == false )
    assert( (  "test5" != token.getNamespacePrefix(4) ) == false )
    assert( (  "http://test1.org/" != token.getNamespaceURI(0) ) == false )
    assert( (  "http://test2.org/" != token.getNamespaceURI(1) ) == false )
    assert( (  "http://test3.org/" != token.getNamespaceURI(2) ) == false )
    assert( (  "http://test4.org/" != token.getNamespaceURI(3) ) == false )
    assert( (  "http://test5.org/" != token.getNamespaceURI(4) ) == false )
    token.clearNamespaces()
    assert( token.getNamespacesLength() == 0 )
    ns = nil
    token = nil
    triple = nil
    attr = nil
  end

end
