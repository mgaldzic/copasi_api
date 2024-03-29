# @file    TestSBMLError.rb
# @brief   SBMLError unit tests, C++ version
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
# src/sbml/test/TestSBMLError.cpp
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

class TestSBMLError < Test::Unit::TestCase

  def test_SBMLError_create
    error = LibSBML::SBMLError.new()
    assert( error != nil )
    error = nil
    error = LibSBML::SBMLError.new(LibSBML::EmptyListInReaction)
    assert( error.getErrorId() == LibSBML::EmptyListInReaction )
    assert( error.getSeverity() == LibSBML::LIBSBML_SEV_ERROR )
    assert( error.getSeverityAsString() ==  "Error"  )
    assert( error.getCategory() == LibSBML::LIBSBML_CAT_SBML )
    assert( error.getCategoryAsString() ==  "General SBML conformance" )
    error = nil
    error = LibSBML::SBMLError.new(LibSBML::OverdeterminedSystem,2,1)
    assert( error.getErrorId() == LibSBML::OverdeterminedSystem )
    assert( error.getSeverity() == LibSBML::LIBSBML_SEV_WARNING )
    assert( error.getSeverityAsString() ==  "Warning"  )
    assert( error.getCategory() == LibSBML::LIBSBML_CAT_SBML )
    assert( error.getCategoryAsString() ==  "General SBML conformance" )
    error = nil
    error = LibSBML::SBMLError.new(LibSBML::OffsetNoLongerValid,2,2)
    assert( error.getErrorId() == LibSBML::OffsetNoLongerValid )
    assert( error.getSeverity() == LibSBML::LIBSBML_SEV_ERROR )
    assert( error.getSeverityAsString() ==  "Error"  )
    assert( error.getCategory() == LibSBML::LIBSBML_CAT_GENERAL_CONSISTENCY )
    assert( error.getCategoryAsString() ==  "SBML component consistency" )
    error = nil
    error = LibSBML::SBMLError.new(LibSBML::NoSBOTermsInL1,2,2)
    assert( error.getErrorId() == LibSBML::NoSBOTermsInL1 )
    assert( error.getSeverity() == LibSBML::LIBSBML_SEV_WARNING )
    assert( error.getSeverityAsString() ==  "Warning"  )
    assert( error.getCategory() == LibSBML::LIBSBML_CAT_SBML_L1_COMPAT )
    assert( error.getCategoryAsString() ==  "Translation to SBML L1V2" )
    error = nil
    error = LibSBML::SBMLError.new(LibSBML::DisallowedMathMLEncodingUse,2,2)
    assert( error.getErrorId() == LibSBML::DisallowedMathMLEncodingUse )
    assert( error.getSeverity() == LibSBML::LIBSBML_SEV_ERROR )
    assert( error.getSeverityAsString() ==  "Error"  )
    assert( error.getCategory() == LibSBML::LIBSBML_CAT_MATHML_CONSISTENCY )
    assert( error.getShortMessage() ==  "Disallowed use of MathML 'encoding' attribute" )
    error = nil
    error = LibSBML::SBMLError.new(LibSBML::DisallowedMathMLEncodingUse,1,2)
    assert( error.getErrorId() == LibSBML::DisallowedMathMLEncodingUse )
    assert( error.getSeverity() == LibSBML::LIBSBML_SEV_NOT_APPLICABLE )
    assert( error.getCategory() == LibSBML::LIBSBML_CAT_MATHML_CONSISTENCY )
    error = nil
    error = LibSBML::SBMLError.new(LibSBML::UnknownError,2,4)
    assert( error.getErrorId() == LibSBML::UnknownError )
    assert( error.getSeverity() == LibSBML::LIBSBML_SEV_FATAL )
    assert( error.getSeverityAsString() ==  "Fatal"  )
    assert( error.getCategory() == LibSBML::LIBSBML_CAT_INTERNAL )
    assert( error.getShortMessage() ==  "Unknown internal libSBML error" )
    error = nil
  end

end
