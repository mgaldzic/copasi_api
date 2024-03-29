/*
 * @file    TestXMLToken.java
 * @brief   XMLToken unit tests
 *
 * @author  Akiya Jouraku (Java conversion)
 * @author  Michael Hucka <mhucka@caltech.edu> 
 *
 * $Id$
 * $HeadURL$
 *
 * ====== WARNING ===== WARNING ===== WARNING ===== WARNING ===== WARNING ======
 *
 * DO NOT EDIT THIS FILE.
 *
 * This file was generated automatically by converting the file located at
 * src/xml/test/TestXMLToken.c
 * using the conversion program dev/utilities/translateTests/translateTests.pl.
 * Any changes made here will be lost the next time the file is regenerated.
 *
 * -----------------------------------------------------------------------------
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
 * -----------------------------------------------------------------------------
 */

package org.sbml.libsbml.test.xml;

import org.sbml.libsbml.*;

import java.io.File;
import java.lang.AssertionError;

public class TestXMLToken {

  static void assertTrue(boolean condition) throws AssertionError
  {
    if (condition == true)
    {
      return;
    }
    throw new AssertionError();
  }

  static void assertEquals(Object a, Object b) throws AssertionError
  {
    if ( (a == null) && (b == null) )
    {
      return;
    }
    else if ( (a == null) || (b == null) )
    {
      throw new AssertionError();
    }
    else if (a.equals(b))
    {
      return;
    }

    throw new AssertionError();
  }

  static void assertNotEquals(Object a, Object b) throws AssertionError
  {
    if ( (a == null) && (b == null) )
    {
      throw new AssertionError();
    }
    else if ( (a == null) || (b == null) )
    {
      return;
    }
    else if (a.equals(b))
    {
      throw new AssertionError();
    }
  }

  static void assertEquals(boolean a, boolean b) throws AssertionError
  {
    if ( a == b )
    {
      return;
    }
    throw new AssertionError();
  }

  static void assertNotEquals(boolean a, boolean b) throws AssertionError
  {
    if ( a != b )
    {
      return;
    }
    throw new AssertionError();
  }

  static void assertEquals(int a, int b) throws AssertionError
  {
    if ( a == b )
    {
      return;
    }
    throw new AssertionError();
  }

  static void assertNotEquals(int a, int b) throws AssertionError
  {
    if ( a != b )
    {
      return;
    }
    throw new AssertionError();
  }

  public void test_XMLToken_attribute_add_remove()
  {
    XMLTriple triple = new  XMLTriple("test","","");
    XMLAttributes attr = new  XMLAttributes();
    XMLToken token = new  XMLToken(triple,attr);
    XMLTriple xt1 = new  XMLTriple("name1", "http://name1.org/", "p1");
    XMLTriple xt2 = new  XMLTriple("name2", "http://name2.org/", "p2");
    XMLTriple xt3 = new  XMLTriple("name3", "http://name3.org/", "p3");
    XMLTriple xt1a = new  XMLTriple("name1", "http://name1a.org/", "p1a");
    XMLTriple xt2a = new  XMLTriple("name2", "http://name2a.org/", "p2a");
    token.addAttr( "name1", "val1", "http://name1.org/", "p1");
    token.addAttr(xt2, "val2");
    assertTrue( token.getAttributesLength() == 2 );
    assertTrue( token.isAttributesEmpty() == false );
    assertTrue( !token.getAttrName(0).equals( "name1") == false );
    assertTrue( !token.getAttrValue(0).equals( "val1" ) == false );
    assertTrue( !token.getAttrURI(0).equals( "http://name1.org/") == false );
    assertTrue( !token.getAttrPrefix(0).equals( "p1"   ) == false );
    assertTrue( !token.getAttrName(1).equals( "name2") == false );
    assertTrue( !token.getAttrValue(1).equals( "val2" ) == false );
    assertTrue( !token.getAttrURI(1).equals( "http://name2.org/") == false );
    assertTrue( !token.getAttrPrefix(1).equals( "p2"   ) == false );
    assertTrue( token.getAttrValue( "name1").equals("") == true );
    assertTrue( token.getAttrValue( "name2").equals("") == true );
    assertTrue( !token.getAttrValue( "name1", "http://name1.org/").equals( "val1" ) == false );
    assertTrue( !token.getAttrValue( "name2", "http://name2.org/").equals( "val2" ) == false );
    assertTrue( !token.getAttrValue(xt1).equals( "val1" ) == false );
    assertTrue( !token.getAttrValue(xt2).equals( "val2" ) == false );
    assertTrue( token.hasAttr(-1) == false );
    assertTrue( token.hasAttr(2) == false );
    assertTrue( token.hasAttr(0) == true );
    assertTrue( token.hasAttr( "name1", "http://name1.org/") == true );
    assertTrue( token.hasAttr( "name2", "http://name2.org/") == true );
    assertTrue( token.hasAttr( "name3", "http://name3.org/") == false );
    assertTrue( token.hasAttr(xt1) == true );
    assertTrue( token.hasAttr(xt2) == true );
    assertTrue( token.hasAttr(xt3) == false );
    token.addAttr( "noprefix", "val3");
    assertTrue( token.getAttributesLength() == 3 );
    assertTrue( token.isAttributesEmpty() == false );
    assertTrue( !token.getAttrName(2).equals( "noprefix") == false );
    assertTrue( !token.getAttrValue(2).equals( "val3"    ) == false );
    assertTrue( token.getAttrURI(2).equals("") == true );
    assertTrue( token.getAttrPrefix(2).equals("") == true );
    assertTrue( !token.getAttrValue( "noprefix").equals(     "val3" ) == false );
    assertTrue( !token.getAttrValue( "noprefix", "").equals( "val3" ) == false );
    assertTrue( token.hasAttr( "noprefix"    ) == true );
    assertTrue( token.hasAttr( "noprefix", "") == true );
    token.addAttr(xt1, "mval1");
    token.addAttr( "name2", "mval2", "http://name2.org/", "p2");
    assertTrue( token.getAttributesLength() == 3 );
    assertTrue( token.isAttributesEmpty() == false );
    assertTrue( !token.getAttrName(0).equals( "name1") == false );
    assertTrue( !token.getAttrValue(0).equals( "mval1") == false );
    assertTrue( !token.getAttrURI(0).equals( "http://name1.org/") == false );
    assertTrue( !token.getAttrPrefix(0).equals( "p1"   ) == false );
    assertTrue( !token.getAttrName(1).equals( "name2"   ) == false );
    assertTrue( !token.getAttrValue(1).equals( "mval2"   ) == false );
    assertTrue( !token.getAttrURI(1).equals( "http://name2.org/") == false );
    assertTrue( !token.getAttrPrefix(1).equals( "p2"      ) == false );
    assertTrue( token.hasAttr(xt1) == true );
    assertTrue( token.hasAttr( "name1", "http://name1.org/") == true );
    token.addAttr( "noprefix", "mval3");
    assertTrue( token.getAttributesLength() == 3 );
    assertTrue( token.isAttributesEmpty() == false );
    assertTrue( !token.getAttrName(2).equals( "noprefix") == false );
    assertTrue( !token.getAttrValue(2).equals( "mval3"   ) == false );
    assertTrue( token.getAttrURI(2).equals("") == true );
    assertTrue( token.getAttrPrefix(2).equals("") == true );
    assertTrue( token.hasAttr( "noprefix") == true );
    assertTrue( token.hasAttr( "noprefix", "") == true );
    token.addAttr(xt1a, "val1a");
    token.addAttr(xt2a, "val2a");
    assertTrue( token.getAttributesLength() == 5 );
    assertTrue( !token.getAttrName(3).equals( "name1") == false );
    assertTrue( !token.getAttrValue(3).equals( "val1a") == false );
    assertTrue( !token.getAttrURI(3).equals( "http://name1a.org/") == false );
    assertTrue( !token.getAttrPrefix(3).equals( "p1a") == false );
    assertTrue( !token.getAttrName(4).equals( "name2") == false );
    assertTrue( !token.getAttrValue(4).equals( "val2a") == false );
    assertTrue( !token.getAttrURI(4).equals( "http://name2a.org/") == false );
    assertTrue( !token.getAttrPrefix(4).equals( "p2a") == false );
    assertTrue( !token.getAttrValue( "name1", "http://name1a.org/").equals( "val1a" ) == false );
    assertTrue( !token.getAttrValue( "name2", "http://name2a.org/").equals( "val2a" ) == false );
    assertTrue( !token.getAttrValue(xt1a).equals( "val1a" ) == false );
    assertTrue( !token.getAttrValue(xt2a).equals( "val2a" ) == false );
    token.removeAttr(xt1a);
    token.removeAttr(xt2a);
    assertTrue( token.getAttributesLength() == 3 );
    token.removeAttr( "name1", "http://name1.org/");
    assertTrue( token.getAttributesLength() == 2 );
    assertTrue( token.isAttributesEmpty() == false );
    assertTrue( !token.getAttrName(0).equals( "name2") == false );
    assertTrue( !token.getAttrValue(0).equals( "mval2") == false );
    assertTrue( !token.getAttrURI(0).equals( "http://name2.org/") == false );
    assertTrue( !token.getAttrPrefix(0).equals( "p2") == false );
    assertTrue( !token.getAttrName(1).equals( "noprefix") == false );
    assertTrue( !token.getAttrValue(1).equals( "mval3") == false );
    assertTrue( token.getAttrURI(1).equals("") == true );
    assertTrue( token.getAttrPrefix(1).equals("") == true );
    assertTrue( token.hasAttr( "name1", "http://name1.org/") == false );
    token.removeAttr(xt2);
    assertTrue( token.getAttributesLength() == 1 );
    assertTrue( token.isAttributesEmpty() == false );
    assertTrue( !token.getAttrName(0).equals( "noprefix") == false );
    assertTrue( !token.getAttrValue(0).equals( "mval3") == false );
    assertTrue( token.getAttrURI(0).equals("") == true );
    assertTrue( token.getAttrPrefix(0).equals("") == true );
    assertTrue( token.hasAttr(xt2) == false );
    assertTrue( token.hasAttr( "name2", "http://name2.org/") == false );
    token.removeAttr( "noprefix");
    assertTrue( token.getAttributesLength() == 0 );
    assertTrue( token.isAttributesEmpty() == true );
    assertTrue( token.hasAttr( "noprefix"    ) == false );
    assertTrue( token.hasAttr( "noprefix", "") == false );
    token = null;
    xt1 = null;
    xt2 = null;
    xt3 = null;
    xt1a = null;
    xt2a = null;
    triple = null;
    attr = null;
  }

  public void test_XMLToken_attribute_set_clear()
  {
    XMLTriple triple = new  XMLTriple("test","","");
    XMLAttributes attr = new  XMLAttributes();
    XMLToken token = new  XMLToken(triple,attr);
    XMLAttributes nattr = new  XMLAttributes();
    XMLTriple xt1 = new  XMLTriple("name1", "http://name1.org/", "p1");
    XMLTriple xt2 = new  XMLTriple("name2", "http://name2.org/", "p2");
    XMLTriple xt3 = new  XMLTriple("name3", "http://name3.org/", "p3");
    XMLTriple xt4 = new  XMLTriple("name4", "http://name4.org/", "p4");
    XMLTriple xt5 = new  XMLTriple("name5", "http://name5.org/", "p5");
    nattr.add(xt1, "val1");
    nattr.add(xt2, "val2");
    nattr.add(xt3, "val3");
    nattr.add(xt4, "val4");
    nattr.add(xt5, "val5");
    token.setAttributes(nattr);
    assertTrue( token.getAttributesLength() == 5 );
    assertTrue( token.isAttributesEmpty() == false );
    assertTrue( !token.getAttrName(0).equals( "name1") == false );
    assertTrue( !token.getAttrValue(0).equals( "val1" ) == false );
    assertTrue( !token.getAttrURI(0).equals( "http://name1.org/") == false );
    assertTrue( !token.getAttrPrefix(0).equals( "p1"   ) == false );
    assertTrue( !token.getAttrName(1).equals( "name2") == false );
    assertTrue( !token.getAttrValue(1).equals( "val2" ) == false );
    assertTrue( !token.getAttrURI(1).equals( "http://name2.org/") == false );
    assertTrue( !token.getAttrPrefix(1).equals( "p2"   ) == false );
    assertTrue( !token.getAttrName(2).equals( "name3") == false );
    assertTrue( !token.getAttrValue(2).equals( "val3" ) == false );
    assertTrue( !token.getAttrURI(2).equals( "http://name3.org/") == false );
    assertTrue( !token.getAttrPrefix(2).equals( "p3"   ) == false );
    assertTrue( !token.getAttrName(3).equals( "name4") == false );
    assertTrue( !token.getAttrValue(3).equals( "val4" ) == false );
    assertTrue( !token.getAttrURI(3).equals( "http://name4.org/") == false );
    assertTrue( !token.getAttrPrefix(3).equals( "p4"   ) == false );
    assertTrue( !token.getAttrName(4).equals( "name5") == false );
    assertTrue( !token.getAttrValue(4).equals( "val5" ) == false );
    assertTrue( !token.getAttrURI(4).equals( "http://name5.org/") == false );
    assertTrue( !token.getAttrPrefix(4).equals( "p5"   ) == false );
    XMLTriple ntriple = new  XMLTriple("test2","http://test2.org/","p2");
    token.setTriple(ntriple);
    assertTrue( !token.getName().equals(   "test2") == false );
    assertTrue( !token.getURI().equals(    "http://test2.org/") == false );
    assertTrue( !token.getPrefix().equals( "p2") == false );
    token.clearAttributes();
    assertTrue( token.getAttributesLength() == 0 );
    assertTrue( token.isAttributesEmpty() != false );
    nattr = null;
    triple = null;
    ntriple = null;
    attr = null;
    token = null;
    xt1 = null;
    xt2 = null;
    xt3 = null;
    xt4 = null;
    xt5 = null;
  }

  public void test_XMLToken_chars()
  {
    XMLToken token;
    token = new  XMLToken("This is text");
    assertTrue( token.isElement() == false );
    assertTrue( token.isEnd() == false );
    assertTrue( token.isStart() == false );
    assertTrue( token.isText() == true );
    assertTrue( token.isEOF() == false );
    assertTrue( !token.getCharacters().equals( "This is text") == false );
    token = null;
  }

  public void test_XMLToken_create()
  {
    XMLToken token;
    XMLTriple triple;
    XMLAttributes attr;
    token = new  XMLToken();
    assertTrue( token != null );
    token = null;
    triple = new  XMLTriple("attr", "uri", "prefix");
    token = new  XMLToken(triple);
    assertTrue( token != null );
    assertTrue( !token.getName().equals( "attr") == false );
    assertTrue( !token.getPrefix().equals( "prefix") == false );
    assertTrue( !token.getURI().equals( "uri") == false );
    token = null;
    attr = new  XMLAttributes();
    assertTrue( attr != null );
    attr.add( "attr2", "value");
    token = new  XMLToken(triple,attr);
    assertTrue( token != null );
    XMLAttributes returnattr = token.getAttributes();
    assertTrue( !returnattr.getName(0).equals( "attr2") == false );
    token = null;
    triple = null;
    attr = null;
  }

  public void test_XMLToken_fields()
  {
    XMLTriple triple;
    XMLToken token;
    triple = new  XMLTriple("attr", "uri", "prefix");
    token = new  XMLToken(triple);
    assertTrue( token.isElement() == true );
    assertTrue( token.isEnd() == true );
    assertTrue( token.isStart() == false );
    assertTrue( token.isText() == false );
    assertTrue( token.isEOF() == false );
    assertTrue( !token.getName().equals( "attr") == false );
    assertTrue( !token.getURI().equals( "uri") == false );
    assertTrue( !token.getPrefix().equals( "prefix") == false );
    token = null;
    triple = null;
  }

  public void test_XMLToken_namespace_add()
  {
    XMLTriple triple = new  XMLTriple("test","","");
    XMLAttributes attr = new  XMLAttributes();
    XMLToken token = new  XMLToken(triple,attr);
    assertTrue( token.getNamespacesLength() == 0 );
    assertTrue( token.isNamespacesEmpty() == true );
    token.addNamespace( "http://test1.org/", "test1");
    assertTrue( token.getNamespacesLength() == 1 );
    assertTrue( token.isNamespacesEmpty() == false );
    token.addNamespace( "http://test2.org/", "test2");
    assertTrue( token.getNamespacesLength() == 2 );
    assertTrue( token.isNamespacesEmpty() == false );
    token.addNamespace( "http://test1.org/", "test1a");
    assertTrue( token.getNamespacesLength() == 3 );
    assertTrue( token.isNamespacesEmpty() == false );
    token.addNamespace( "http://test1.org/", "test1a");
    assertTrue( token.getNamespacesLength() == 3 );
    assertTrue( token.isNamespacesEmpty() == false );
    assertTrue( ! (token.getNamespaceIndex( "http://test1.org/") == -1) );
    token = null;
    triple = null;
    attr = null;
  }

  public void test_XMLToken_namespace_get()
  {
    XMLTriple triple = new  XMLTriple("test","","");
    XMLAttributes attr = new  XMLAttributes();
    XMLToken token = new  XMLToken(triple,attr);
    token.addNamespace( "http://test1.org/", "test1");
    token.addNamespace( "http://test2.org/", "test2");
    token.addNamespace( "http://test3.org/", "test3");
    token.addNamespace( "http://test4.org/", "test4");
    token.addNamespace( "http://test5.org/", "test5");
    token.addNamespace( "http://test6.org/", "test6");
    token.addNamespace( "http://test7.org/", "test7");
    token.addNamespace( "http://test8.org/", "test8");
    token.addNamespace( "http://test9.org/", "test9");
    assertTrue( token.getNamespacesLength() == 9 );
    assertTrue( token.getNamespaceIndex( "http://test1.org/") == 0 );
    assertTrue( !token.getNamespacePrefix(1).equals( "test2") == false );
    assertTrue( !token.getNamespacePrefix( "http://test1.org/").equals(		      "test1") == false );
    assertTrue( !token.getNamespaceURI(1).equals( "http://test2.org/") == false );
    assertTrue( !token.getNamespaceURI( "test2").equals(		      "http://test2.org/") == false );
    assertTrue( token.getNamespaceIndex( "http://test1.org/") == 0 );
    assertTrue( token.getNamespaceIndex( "http://test2.org/") == 1 );
    assertTrue( token.getNamespaceIndex( "http://test5.org/") == 4 );
    assertTrue( token.getNamespaceIndex( "http://test9.org/") == 8 );
    assertTrue( token.getNamespaceIndex( "http://testX.org/") == -1 );
    assertTrue( token.hasNamespaceURI( "http://test1.org/") != false );
    assertTrue( token.hasNamespaceURI( "http://test2.org/") != false );
    assertTrue( token.hasNamespaceURI( "http://test5.org/") != false );
    assertTrue( token.hasNamespaceURI( "http://test9.org/") != false );
    assertTrue( token.hasNamespaceURI( "http://testX.org/") == false );
    assertTrue( token.getNamespaceIndexByPrefix( "test1") == 0 );
    assertTrue( token.getNamespaceIndexByPrefix( "test5") == 4 );
    assertTrue( token.getNamespaceIndexByPrefix( "test9") == 8 );
    assertTrue( token.getNamespaceIndexByPrefix( "testX") == -1 );
    assertTrue( token.hasNamespacePrefix( "test1") != false );
    assertTrue( token.hasNamespacePrefix( "test5") != false );
    assertTrue( token.hasNamespacePrefix( "test9") != false );
    assertTrue( token.hasNamespacePrefix( "testX") == false );
    assertTrue( token.hasNamespaceNS( "http://test1.org/", "test1") != false );
    assertTrue( token.hasNamespaceNS( "http://test5.org/", "test5") != false );
    assertTrue( token.hasNamespaceNS( "http://test9.org/", "test9") != false );
    assertTrue( token.hasNamespaceNS( "http://testX.org/", "testX") == false );
    token = null;
    triple = null;
    attr = null;
  }

  public void test_XMLToken_namespace_remove()
  {
    XMLTriple triple = new  XMLTriple("test","","");
    XMLAttributes attr = new  XMLAttributes();
    XMLToken token = new  XMLToken(triple,attr);
    token.addNamespace( "http://test1.org/", "test1");
    token.addNamespace( "http://test2.org/", "test2");
    token.addNamespace( "http://test3.org/", "test3");
    token.addNamespace( "http://test4.org/", "test4");
    token.addNamespace( "http://test5.org/", "test5");
    assertTrue( token.getNamespacesLength() == 5 );
    token.removeNamespace(4);
    assertTrue( token.getNamespacesLength() == 4 );
    token.removeNamespace(3);
    assertTrue( token.getNamespacesLength() == 3 );
    token.removeNamespace(2);
    assertTrue( token.getNamespacesLength() == 2 );
    token.removeNamespace(1);
    assertTrue( token.getNamespacesLength() == 1 );
    token.removeNamespace(0);
    assertTrue( token.getNamespacesLength() == 0 );
    token.addNamespace( "http://test1.org/", "test1");
    token.addNamespace( "http://test2.org/", "test2");
    token.addNamespace( "http://test3.org/", "test3");
    token.addNamespace( "http://test4.org/", "test4");
    token.addNamespace( "http://test5.org/", "test5");
    assertTrue( token.getNamespacesLength() == 5 );
    token.removeNamespace(0);
    assertTrue( token.getNamespacesLength() == 4 );
    token.removeNamespace(0);
    assertTrue( token.getNamespacesLength() == 3 );
    token.removeNamespace(0);
    assertTrue( token.getNamespacesLength() == 2 );
    token.removeNamespace(0);
    assertTrue( token.getNamespacesLength() == 1 );
    token.removeNamespace(0);
    assertTrue( token.getNamespacesLength() == 0 );
    token = null;
    triple = null;
    attr = null;
  }

  public void test_XMLToken_namespace_remove_by_prefix()
  {
    XMLTriple triple = new  XMLTriple("test","","");
    XMLAttributes attr = new  XMLAttributes();
    XMLToken token = new  XMLToken(triple,attr);
    token.addNamespace( "http://test1.org/", "test1");
    token.addNamespace( "http://test2.org/", "test2");
    token.addNamespace( "http://test3.org/", "test3");
    token.addNamespace( "http://test4.org/", "test4");
    token.addNamespace( "http://test5.org/", "test5");
    assertTrue( token.getNamespacesLength() == 5 );
    token.removeNamespace( "test1");
    assertTrue( token.getNamespacesLength() == 4 );
    token.removeNamespace( "test2");
    assertTrue( token.getNamespacesLength() == 3 );
    token.removeNamespace( "test3");
    assertTrue( token.getNamespacesLength() == 2 );
    token.removeNamespace( "test4");
    assertTrue( token.getNamespacesLength() == 1 );
    token.removeNamespace( "test5");
    assertTrue( token.getNamespacesLength() == 0 );
    token.addNamespace( "http://test1.org/", "test1");
    token.addNamespace( "http://test2.org/", "test2");
    token.addNamespace( "http://test3.org/", "test3");
    token.addNamespace( "http://test4.org/", "test4");
    token.addNamespace( "http://test5.org/", "test5");
    assertTrue( token.getNamespacesLength() == 5 );
    token.removeNamespace( "test5");
    assertTrue( token.getNamespacesLength() == 4 );
    token.removeNamespace( "test4");
    assertTrue( token.getNamespacesLength() == 3 );
    token.removeNamespace( "test3");
    assertTrue( token.getNamespacesLength() == 2 );
    token.removeNamespace( "test2");
    assertTrue( token.getNamespacesLength() == 1 );
    token.removeNamespace( "test1");
    assertTrue( token.getNamespacesLength() == 0 );
    token.addNamespace( "http://test1.org/", "test1");
    token.addNamespace( "http://test2.org/", "test2");
    token.addNamespace( "http://test3.org/", "test3");
    token.addNamespace( "http://test4.org/", "test4");
    token.addNamespace( "http://test5.org/", "test5");
    assertTrue( token.getNamespacesLength() == 5 );
    token.removeNamespace( "test3");
    assertTrue( token.getNamespacesLength() == 4 );
    token.removeNamespace( "test1");
    assertTrue( token.getNamespacesLength() == 3 );
    token.removeNamespace( "test4");
    assertTrue( token.getNamespacesLength() == 2 );
    token.removeNamespace( "test5");
    assertTrue( token.getNamespacesLength() == 1 );
    token.removeNamespace( "test2");
    assertTrue( token.getNamespacesLength() == 0 );
    token = null;
    triple = null;
    attr = null;
  }

  public void test_XMLToken_namespace_set_clear()
  {
    XMLTriple triple = new  XMLTriple("test","","");
    XMLAttributes attr = new  XMLAttributes();
    XMLToken token = new  XMLToken(triple,attr);
    XMLNamespaces ns = new  XMLNamespaces();
    assertTrue( token.getNamespacesLength() == 0 );
    assertTrue( token.isNamespacesEmpty() == true );
    ns.add( "http://test1.org/", "test1");
    ns.add( "http://test2.org/", "test2");
    ns.add( "http://test3.org/", "test3");
    ns.add( "http://test4.org/", "test4");
    ns.add( "http://test5.org/", "test5");
    token.setNamespaces(ns);
    assertTrue( token.getNamespacesLength() == 5 );
    assertTrue( token.isNamespacesEmpty() == false );
    assertTrue( !token.getNamespacePrefix(0).equals( "test1") == false );
    assertTrue( !token.getNamespacePrefix(1).equals( "test2") == false );
    assertTrue( !token.getNamespacePrefix(2).equals( "test3") == false );
    assertTrue( !token.getNamespacePrefix(3).equals( "test4") == false );
    assertTrue( !token.getNamespacePrefix(4).equals( "test5") == false );
    assertTrue( !token.getNamespaceURI(0).equals( "http://test1.org/") == false );
    assertTrue( !token.getNamespaceURI(1).equals( "http://test2.org/") == false );
    assertTrue( !token.getNamespaceURI(2).equals( "http://test3.org/") == false );
    assertTrue( !token.getNamespaceURI(3).equals( "http://test4.org/") == false );
    assertTrue( !token.getNamespaceURI(4).equals( "http://test5.org/") == false );
    token.clearNamespaces();
    assertTrue( token.getNamespacesLength() == 0 );
    ns = null;
    token = null;
    triple = null;
    attr = null;
  }

  /**
   * Loads the SWIG-generated libSBML Java module when this class is
   * loaded, or reports a sensible diagnostic message about why it failed.
   */
  static
  {
    String varname;
    String shlibname;

    if (System.getProperty("mrj.version") != null)
    {
      varname = "DYLD_LIBRARY_PATH";    // We're on a Mac.
      shlibname = "libsbmlj.jnilib and/or libsbml.dylib";
    }
    else
    {
      varname = "LD_LIBRARY_PATH";      // We're not on a Mac.
      shlibname = "libsbmlj.so and/or libsbml.so";
    }

    try
    {
      System.loadLibrary("sbmlj");
      // For extra safety, check that the jar file is in the classpath.
      Class.forName("org.sbml.libsbml.libsbml");
    }
    catch (SecurityException e)
    {
      e.printStackTrace();
      System.err.println("Could not load the libSBML library files due to a"+
                         " security exception.\n");
      System.exit(1);
    }
    catch (UnsatisfiedLinkError e)
    {
      e.printStackTrace();
      System.err.println("Error: could not link with the libSBML library files."+
                         " It is likely\nyour " + varname +
                         " environment variable does not include the directories\n"+
                         "containing the " + shlibname + " library files.\n");
      System.exit(1);
    }
    catch (ClassNotFoundException e)
    {
      e.printStackTrace();
      System.err.println("Error: unable to load the file libsbmlj.jar."+
                         " It is likely\nyour -classpath option and CLASSPATH" +
                         " environment variable\n"+
                         "do not include the path to libsbmlj.jar.\n");
      System.exit(1);
    }
  }
}
