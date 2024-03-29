# @file    TestModelCreator_newSetters.rb
# @brief   ModelCreator unit tests
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
# src/annotation/test/TestModelCreator_newSetters.c
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

class TestModelCreator_newSetters < Test::Unit::TestCase

  def test_ModelCreator_setEmail
    mc = LibSBML::ModelCreator.new()
    assert( mc != nil )
    i = mc.setEmail( "Keating")
    assert( i == LibSBML::LIBSBML_OPERATION_SUCCESS )
    assert( mc.isSetEmail() == true )
    assert ((  "Keating" == mc.getEmail() ))
    i = mc.setEmail( "")
    assert( i == LibSBML::LIBSBML_OPERATION_SUCCESS )
    assert( mc.isSetEmail() == false )
    i = mc.setEmail( "Keating")
    assert( mc.isSetEmail() == true )
    i = mc.unsetEmail()
    assert( i == LibSBML::LIBSBML_OPERATION_SUCCESS )
    assert( mc.isSetEmail() == false )
    mc = nil
  end

  def test_ModelCreator_setFamilyName
    mc = LibSBML::ModelCreator.new()
    assert( mc != nil )
    i = mc.setFamilyName( "Keating")
    assert( i == LibSBML::LIBSBML_OPERATION_SUCCESS )
    assert( mc.isSetFamilyName() == true )
    assert ((  "Keating" == mc.getFamilyName() ))
    i = mc.setFamilyName( "")
    assert( i == LibSBML::LIBSBML_OPERATION_SUCCESS )
    assert( mc.isSetFamilyName() == false )
    i = mc.setFamilyName( "Keating")
    assert( mc.isSetFamilyName() == true )
    i = mc.unsetFamilyName()
    assert( i == LibSBML::LIBSBML_OPERATION_SUCCESS )
    assert( mc.isSetFamilyName() == false )
    mc = nil
  end

  def test_ModelCreator_setGivenName
    mc = LibSBML::ModelCreator.new()
    assert( mc != nil )
    i = mc.setGivenName( "Sarah")
    assert( i == LibSBML::LIBSBML_OPERATION_SUCCESS )
    assert( mc.isSetGivenName() == true )
    assert ((  "Sarah" == mc.getGivenName() ))
    i = mc.setGivenName( "")
    assert( i == LibSBML::LIBSBML_OPERATION_SUCCESS )
    assert( mc.isSetGivenName() == false )
    i = mc.setGivenName( "Sarah")
    assert( mc.isSetGivenName() == true )
    i = mc.unsetGivenName()
    assert( i == LibSBML::LIBSBML_OPERATION_SUCCESS )
    assert( mc.isSetGivenName() == false )
    mc = nil
  end

  def test_ModelCreator_setOrganization
    mc = LibSBML::ModelCreator.new()
    assert( mc != nil )
    i = mc.setOrganization( "Caltech")
    assert( i == LibSBML::LIBSBML_OPERATION_SUCCESS )
    assert( mc.isSetOrganization() == true )
    assert ((  "Caltech" == mc.getOrganization() ))
    i = mc.setOrganization( "")
    assert( i == LibSBML::LIBSBML_OPERATION_SUCCESS )
    assert( mc.isSetOrganization() == false )
    i = mc.setOrganization( "Caltech")
    assert( mc.isSetOrganization() == true )
    i = mc.unsetOrganization()
    assert( i == LibSBML::LIBSBML_OPERATION_SUCCESS )
    assert( mc.isSetOrganization() == false )
    mc = nil
  end

end
