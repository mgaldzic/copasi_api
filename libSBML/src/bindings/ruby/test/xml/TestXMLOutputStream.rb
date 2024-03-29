# @file    TestXMLOutputStream.rb
# @brief   XMLOutputStream unit tests
#
# @author  Akiya Jouraku (Ruby conversion)
# @author  Sarah Keating 
#
# $Id$
# $HeadURL$
#
# ====== WARNING ===== WARNING ===== WARNING ===== WARNING ===== WARNING ======
#
# DO NOT EDIT THIS FILE.
#
# This file was generated automatically by converting the file located at
# src/xml/test/TestXMLOutputStream.c
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

class TestXMLOutputStream < Test::Unit::TestCase

  def test_XMLOutputStream_CharacterReference
    oss = LibSBML::Ostringstream.new
    stream = LibSBML::XMLOutputStream.new(oss,"",false)
    stream.startElement( "testcr")
    stream.writeAttribute( "chars",    "one"     )
    stream.writeAttribute( "amp",      "&"       )
    stream.writeAttribute( "deccr",    "&#0168;"  )
    stream.writeAttribute( "hexcr",    "&#x00a8;")
    stream.writeAttribute( "lhexcr",   "&#x00A8;")
    stream.writeAttribute( "nodeccr1", "&#01688"  )
    stream.writeAttribute( "nodeccr2", "&#;"     )
    stream.writeAttribute( "nodeccr3", "&#00a8;" )
    stream.writeAttribute( "nodeccr4", "&#00A8;" )
    stream.writeAttribute( "nohexcr1", "&#x;"    )
    stream.writeAttribute( "nohexcr2", "&#xABCD" )
    stream.endElement( "testcr")
    expected =  "<testcr chars=\"one\" amp=\"&amp;\" deccr=\"&#0168;\" hexcr=\"&#x00a8;\" lhexcr=\"&#x00A8;\" nodeccr1=\"&amp;#01688\" nodeccr2=\"&amp;#;\" nodeccr3=\"&amp;#00a8;\" nodeccr4=\"&amp;#00A8;\" nohexcr1=\"&amp;#x;\" nohexcr2=\"&amp;#xABCD\"/>";
    s = oss.str()
    assert (( expected == s ))
    stream = nil
  end

  def test_XMLOutputStream_Elements
    d = 2.4
    l = 123456789
    ui = 5
    i = -3
    oss = LibSBML::Ostringstream.new
    stream = LibSBML::XMLOutputStream.new(oss,"",false)
    stream.startElement( "fred")
    stream.writeAttribute( "chars", "two")
    stream.writeAttributeBool( "bool",true)
    stream.writeAttribute( "double",d)
    stream.writeAttribute( "long",l)
    stream.writeAttribute( "uint",ui)
    stream.writeAttribute( "int",i)
    stream.endElement( "fred")
    expected =  "<fred chars=\"two\" bool=\"true\" double=\"2.4\" long=\"123456789\" uint=\"5\" int=\"-3\"/>";
    s = oss.str()
    assert (( expected == s ))
    stream = nil
  end

  def test_XMLOutputStream_PredefinedEntity
    oss = LibSBML::Ostringstream.new
    stream = LibSBML::XMLOutputStream.new(oss,"",false)
    stream.startElement( "testpde")
    stream.writeAttribute( "amp",     "&"     )
    stream.writeAttribute( "apos",    "'"     )
    stream.writeAttribute( "gt",      ">"     )
    stream.writeAttribute( "lt",      "<"     )
    stream.writeAttribute( "quot",    "\""    )
    stream.writeAttribute( "pdeamp",  "&amp;" )
    stream.writeAttribute( "pdeapos", "&apos;")
    stream.writeAttribute( "pdegt",   "&gt;"  )
    stream.writeAttribute( "pdelt",   "&lt;"  )
    stream.writeAttribute( "pdequot", "&quot;")
    stream.endElement( "testpde")
    expected =  "<testpde amp=\"&amp;\" apos=\"&apos;\" gt=\"&gt;\" lt=\"&lt;\" quot=\"&quot;\" pdeamp=\"&amp;\" pdeapos=\"&apos;\" pdegt=\"&gt;\" pdelt=\"&lt;\" pdequot=\"&quot;\"/>";
    s = oss.str()
    assert (( expected == s ))
    stream = nil
  end

  def test_XMLOutputStream_createStdout
    stream = LibSBML::XMLOutputStream.new(LibSBML::cout,"UTF-8",false)
    assert( stream != nil )
    stream = nil
  end

  def test_XMLOutputStream_createStdoutWithProgramInfo
    stream = LibSBML::XMLOutputStream.new(LibSBML::cout,"UTF-8",false, "foo", "bar")
    assert( stream != nil )
    stream = nil
  end

  def test_XMLOutputStream_createString
    expected =  "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n";
    oss = LibSBML::Ostringstream.new
    stream = LibSBML::XMLOutputStream.new(oss,"UTF-8",true)
    assert( stream != nil )
    str = oss.str()
    assert (( expected == str ))
    stream = nil
  end

  def test_XMLOutputStream_createStringWithProgramInfo
    expected =  "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n";
    oss = LibSBML::Ostringstream.new
    stream = LibSBML::XMLOutputStream.new(oss,"UTF-8",true, "", "")
    assert( stream != nil )
    str = oss.str()
    assert (( expected == str ))
    stream = nil
  end

  def test_XMLOutputStream_startEnd
    oss = LibSBML::Ostringstream.new
    stream = LibSBML::XMLOutputStream.new(oss,"",false)
    assert( stream != nil )
    stream.startEndElement( "id")
    str = oss.str()
    assert ((  "<id/>" == str ))
    stream = nil
  end

end
