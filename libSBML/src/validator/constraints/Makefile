## src/validator/constraints/Makefile.  Generated from Makefile.in by configure.
##
## Filename    : Makefile.in
## Description : Makefile template for libsbml src/validator/constraints
## Author(s)   : SBML Team <sbml-team@caltech.edu>
## Organization: California Institute of Technology
## Created     : 2005-03-02
## Revision    : $Id: Makefile.in 11269 2010-05-27 11:27:49Z sarahkeating $
## $HeadURL: https://sbml.svn.sourceforge.net/svnroot/sbml/trunk/libsbml/src/validator/constraints/Makefile.in $
##
## Copyright 2005 California Institute of Technology and
## Japan Science and Technology Corporation.
##
## This library is free software; you can redistribute it and/or modify it
## under the terms of the GNU Lesser General Public License as published
## by the Free Software Foundation; either version 2.1 of the License, or
## any later version.
##
## This library is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY, WITHOUT EVEN THE IMPLIED WARRANTY OF
## MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE.  The software and
## documentation provided hereunder is on an "as is" basis, and the
## California Institute of Technology and Japan Science and Technology
## Corporation have no obligations to provide maintenance, support,
## updates, enhancements or modifications.  In no event shall the
## California Institute of Technology or the Japan Science and Technology
## Corporation be liable to any party for direct, indirect, special,
## incidental or consequential damages, including lost profits, arising
## out of the use of this software and its documentation, even if the
## California Institute of Technology and/or Japan Science and Technology
## Corporation have been advised of the possibility of such damage.  See
## the GNU Lesser General Public License for more details.
##
## You should have received a copy of the GNU Lesser General Public License
## along with this library; if not, write to the Free Software Foundation,
## Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA.
##
## The original code contained here was initially developed by:
##
##     Ben Bornstein
##     SBML Team
##     Control and Dynamical Systems, MC 107-81
##     California Institute of Technology
##     Pasadena, CA, 91125, USA
##
##     http://www.sbml.org
##     mailto:sbml-team@caltech.edu
##
## Contributor(s):
##   Michael Hucka <mhucka@caltech.edu> Wrote this Makefile.in.
##


# -----------------------------------------------------------------------------
# Configuration variables
# -----------------------------------------------------------------------------
# Some of the following are substituted automatically by `configure'.  If
# you are looking at "Makefile", do not edit these values; instead, run the
# configure script at the top level of the src tree.  It will recreate
# "Makefile".

include ../../../config/makefile-common-vars.mk

# `srcdir' points to the current directory, but should be set by configure.
# `subdir' must be set manually to the relative dir under `srcdir'.  Don't 
# set `subdir' to an absolute path, or some `make' actions will fail.

srcdir    = .
thisdir   = src/validator/constraints

# Variables `headers', `sources', `libraries', `extra_CPPFLAGS',
# `extra_CXXFLAGS', `extra_LDFLAGS' and `distfiles' are used by the default
# rules in `makefile-common-actions.mk' included at the end.

sources =                                     \
  ArgumentsUnitsCheck.cpp                     \
  ArgumentsUnitsCheckWarnings.cpp             \
  AssignmentCycles.cpp                        \
  AssignmentRuleOrdering.cpp                  \
  CiElementMathCheck.cpp                      \
  CompartmentOutsideCycles.cpp                \
  ConsistencyConstraints.cpp                  \
  DuplicateTopLevelAnnotation.cpp             \
  EqualityArgsMathCheck.cpp                   \
  ExponentUnitsCheck.cpp		          			  \
  FunctionApplyMathCheck.cpp                  \
  FunctionDefinitionVars.cpp                  \
  FunctionNoArgsMathCheck.cpp                 \
  FunctionReferredToExists.cpp                \
  IdBase.cpp                                  \
  IdentifierConsistencyConstraints.cpp        \
  IdList.cpp                                  \
  InternalConsistencyConstraints.cpp          \
  KineticLawVars.cpp                          \
  KineticLawUnitsCheck.cpp                    \
  L1CompatibilityConstraints.cpp              \
  L2v1CompatibilityConstraints.cpp            \
  L2v2CompatibilityConstraints.cpp            \
  L2v3CompatibilityConstraints.cpp            \
  L2v4CompatibilityConstraints.cpp            \
  L3v1CompatibilityConstraints.cpp            \
  LambdaMathCheck.cpp                         \
  LocalParameterMathCheck.cpp                 \
  LocalParameterShadowsIdInModel.cpp          \
  LogicalArgsMathCheck.cpp                    \
  MathMLBase.cpp                              \
  MathMLConsistencyConstraints.cpp            \
  ModelingPracticeConstraints.cpp             \
  NumberArgsMathCheck.cpp                     \
  NumericArgsMathCheck.cpp                    \
  NumericReturnMathCheck.cpp                  \
  OverDeterminedCheck.cpp                     \
  OverdeterminedConstraints.cpp				  \
  PieceBooleanMathCheck.cpp                   \
  PiecewiseValueMathCheck.cpp                 \
  PowerUnitsCheck.cpp						  \
  SBOConsistencyConstraints.cpp               \
  SpeciesReactionOrRule.cpp                   \
  StoichiometryMathVars.cpp                   \
  UniqueIdBase.cpp                            \
  UniqueIdsForUnitDefinitions.cpp             \
  UniqueIdsInKineticLaw.cpp                   \
  UniqueIdsInModel.cpp                        \
  UniqueMetaId.cpp                            \
  UniqueSpeciesTypesInCompartment.cpp         \
  UniqueSymbolsInInitialAssignments.cpp       \
  UniqueVarsInEventAssignments.cpp            \
  UniqueVarsInEventsAndRules.cpp              \
  UniqueVarsInInitialAssignmentsAndRules.cpp  \
  UniqueVarsInRules.cpp                       \
  UnitConsistencyConstraints.cpp              \
  ValidCnUnitsValue.cpp                       \
  UnitsBase.cpp

header_inst_prefix = validator/constraints

headers =               \
  IdBase.h              \
  IdList.h              \
  MathMLBase.h          \
  OverDeterminedCheck.h \
  UniqueIdBase.h        \
  UniqueMetaId.h        \
  UnitsBase.h

other_headers = \
  ArgumentsUnitsCheck.h                    \
  ArgumentsUnitsCheckWarnings.h            \
  AssignmentCycles.h                       \
  AssignmentRuleOrdering.h                 \
  CiElementMathCheck.h                     \
  CompartmentOutsideCycles.h               \
  DuplicateTopLevelAnnotation.h            \
  EqualityArgsMathCheck.h                  \
  ExponentUnitsCheck.h                     \
  FunctionApplyMathCheck.h                 \
  FunctionDefinitionVars.h                 \
  FunctionNoArgsMathCheck.h                \
  FunctionReferredToExists.h               \
  KineticLawUnitsCheck.h                   \
  KineticLawVars.h                         \
  LambdaMathCheck.h                        \
  LocalParameterMathCheck.h                \
  LocalParameterShadowsIdInModel.h         \
  LogicalArgsMathCheck.h                   \
  Makefile                                 \
  Makefile.in                              \
  NumberArgsMathCheck.h                    \
  NumericArgsMathCheck.h                   \
  NumericReturnMathCheck.h                 \
  PieceBooleanMathCheck.h                  \
  PiecewiseValueMathCheck.h                \
  PowerUnitsCheck.h                        \
  SpeciesReactionOrRule.h                  \
  StoichiometryMathVars.h                  \
  UniqueIdsForUnitDefinitions.h            \
  UniqueIdsInKineticLaw.h                  \
  UniqueIdsInModel.h                       \
  UniqueSpeciesTypesInCompartment.h        \
  UniqueSymbolsInInitialAssignments.h      \
  UniqueVarsInEventAssignments.h           \
  UniqueVarsInEventsAndRules.h             \
  UniqueVarsInInitialAssignmentsAndRules.h \
  UniqueVarsInRules.h                      \
  ValidCnUnitsValue.h

# `distfiles' determines the files and directories included in a distribution.
# `distfiles_exclude' is used to filter out specific items, so that even if
# they appear in `distfiles', they will not be copied in the distribution.

distfiles = $(sources) $(headers) $(other_headers) Makefile.in

ifeq "$(HOST_TYPE)" "cygwin"

# Under cygwin, at least in gcc 3.4.4, the standard gcc STL class files
# generate compiler warnings about use of uninitialized variables.  To
# avoid confusing libSBML users, ignore this particular warning.  This is
# not very safe because it can mask a real problem in OUR code, but it's
# not easy to figure out what else to do, and anyway, if there are problems
# we will probably catch them when testing on other platforms.

problemversions = $(shell expr match `gcc -dumpversion` "3.4.4")

CXXFLAGS += $(shell if (($(problemversions))); then echo -Wno-uninitialized; fi)

endif


# -----------------------------------------------------------------------------
# Primary build actions
# -----------------------------------------------------------------------------

# The default action is to remake everything.

all: Makefile default


# -----------------------------------------------------------------------------
# Checking.
# -----------------------------------------------------------------------------

check:


# -----------------------------------------------------------------------------
# Formatting documentation.
# -----------------------------------------------------------------------------

dvi:

pdf:

ps:

info:

html:


# -----------------------------------------------------------------------------
# Installation
# -----------------------------------------------------------------------------

# This library doesn't actually get installed.  The main libsbml library
# (in ../src) slurps in the contents when it's built.

install: all install-headers

installcheck: all installcheck-headers

uninstall: uninstall-headers


# -----------------------------------------------------------------------------
# Cleaning.
# -----------------------------------------------------------------------------

clean: clean-normal

distclean: distclean-normal

mostlyclean: mostlyclean-normal

maintainer-clean: maintainer-clean-normal


# -----------------------------------------------------------------------------
# Creating distribution (for libSBML maintainers only)
# -----------------------------------------------------------------------------

dist: all dist-normal

distcheck: all distcheck-normal


# -----------------------------------------------------------------------------
# Miscellaneous
# -----------------------------------------------------------------------------

include ../../../config/makefile-common-actions.mk


# -----------------------------------------------------------------------------
# End.
# -----------------------------------------------------------------------------
