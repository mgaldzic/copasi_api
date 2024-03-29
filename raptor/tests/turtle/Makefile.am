# -*- Mode: Makefile -*-
#
# Makefile.am - automake file for Raptor Turtle tests
#
# Copyright (C) 2003-2009, David Beckett http://www.dajobe.org/
# Copyright (C) 2003-2004, University of Bristol, UK http://www.bristol.ac.uk/
# 
# This package is Free Software and part of Redland http://librdf.org/
# 
# It is licensed under the following three licenses as alternatives:
#   1. GNU Lesser General Public License (LGPL) V2.1 or any newer version
#   2. GNU General Public License (GPL) V2 or any newer version
#   3. Apache License, V2.0 or any newer version
# 
# You may not use this file except in compliance with at least one of
# the above three licenses.
# 
# See LICENSE.html or LICENSE.txt at the top of this package for the
# complete terms and further detail along with the license texts for
# the licenses in COPYING.LIB, COPYING and LICENSE-2.0.txt respectively.
# 
# 

RDFXML_DIR = ../rdfxml
NTRIPLES_DIR = ../ntriples

TEST_FILES=test-00.ttl test-01.ttl test-02.ttl test-03.ttl \
test-04.ttl test-05.ttl test-06.ttl test-07.ttl test-08.ttl \
test-09.ttl test-10.ttl test-11.ttl test-12.ttl test-13.ttl \
test-14.ttl test-15.ttl test-16.ttl test-17.ttl test-18.ttl \
test-19.ttl test-20.ttl test-21.ttl test-22.ttl test-23.ttl \
test-24.ttl test-25.ttl test-26.ttl test-27.ttl \
test-29.ttl test-30.ttl test-33.ttl \
rdf-schema.ttl \
rdfs-namespace.ttl \
rdfq-results.ttl

TEST_BAD_FILES=bad-00.ttl bad-01.ttl bad-02.ttl bad-03.ttl \
bad-04.ttl bad-05.ttl bad-06.ttl bad-07.ttl bad-08.ttl bad-09.ttl \
bad-10.ttl bad-11.ttl bad-12.ttl bad-13.ttl bad-14.ttl \
bad-17.ttl bad-18.ttl bad-19.ttl bad-20.ttl bad-21.ttl bad-22.ttl

TEST_WARN_FILES=

TEST_OUT_FILES=test-00.out test-01.out test-02.out test-03.out \
test-04.out test-05.out test-06.out test-07.out test-08.out \
test-09.out test-10.out test-11.out test-12.out test-13.out \
test-14.out test-15.out test-16.out test-17.out test-18.out \
test-19.out test-20.out test-21.out test-22.out test-23.out \
test-24.out test-25.out test-26.out test-27.out test-28.out \
test-29.out test-30.out test-33.out \
rdf-schema.out \
rdfs-namespace.out \
rdfq-results.out

TEST_SERIALIZE_FILES= \
test-28.ttl test-31.ttl test-32.ttl test-34.ttl

TEST_SERIALIZE_OUT_FILES= \
test-28-out.ttl test-31-out.ttl test-32-out.ttl test-34-out.ttl

TURTLE_HACK_OUT_FILES= \
ex-38-turtle.out

TURTLE_SERIALIZE_RDF_FILES=ex-62.rdf


TEST_WARN_OUT_FILES=

TEST_MANIFEST_FILES=manifest.ttl manifest-bad.ttl

# Used to make N-triples output consistent
BASE_URI=http://www.w3.org/2001/sw/DataAccess/df1/tests/
# for rdf-schema.ttl
RDF_NS_URI=http://www.w3.org/1999/02/22-rdf-syntax-ns
# for rdfs-namespace.ttl
RDFS_NS_URI=http://www.w3.org/2000/01/rdf-schema

ALL_TEST_FILES= README.txt \
	$(TEST_FILES) \
	$(TEST_BAD_FILES) \
	$(TEST_WARN_FILES) \
	$(TEST_OUT_FILES) \
	$(TEST_WARN_OUT_FILES) \
	$(TEST_MANIFEST_FILES) \
	$(TEST_SERIALIZE_FILES) \
	$(TEST_SERIALIZE_OUT_FILES) \
	$(TURTLE_SERIALIZE_RDF_TEST_FILES) \
	$(TURTLE_SERIALIZE_RDF_FILES) \
	$(TURTLE_HACK_OUT_FILES)

EXTRA_DIST = $(ALL_TEST_FILES)

build-rapper:
	@(cd $(top_builddir)/utils ; $(MAKE) rapper$(EXEEXT))

build-rdfdiff:
	@(cd $(top_builddir)/utils ; $(MAKE) rdfdiff$(EXEEXT))

check-local: check-rdf check-bad-rdf check-turtle-serialize \
check-turtle-serialize-syntax check-turtle-parse-ntriples \
check-turtle-serialize-rdf

# check-warn-rdf

check-rdf: build-rapper $(TEST_FILES)
	@result=0; \
	$(ECHO) "Testing legal Turtle"; \
	for test in $(TEST_FILES); do \
	  name=`basename $$test .ttl` ; \
	  baseuri=$(BASE_URI)$$test; \
	  $(ECHO) $(ECHO_N) "Checking $$test $(ECHO_C)"; \
	  $(top_builddir)/utils/rapper -q -i turtle -o ntriples $(srcdir)/$$test $$baseuri > $$name.res 2> $$name.err; \
	  status=$$?; \
	  if test $$status != 0 -a $$status != 2 ; then \
	    $(ECHO) "FAILED returned status $$status"; result=1; \
	    $(ECHO) $(top_builddir)/utils/rapper -q -i turtle -o ntriples $(srcdir)/$$test $$baseuri '>' $$name.res; \
	    cat $$name.err; \
	  elif cmp $(srcdir)/$$name.out $$name.res >/dev/null 2>&1; then \
	    if test $$status = 2 ; then \
	      $(ECHO) "ok with warnings"; grep Warning $$name.err; \
	    else \
	      $(ECHO) "ok"; \
	    fi; \
	  else \
	    $(ECHO) "FAILED"; result=1; \
	    $(ECHO) $(top_builddir)/utils/rapper -q -i turtle -o ntriples $(srcdir)/$$test $$baseuri '>' $$name.res; \
	    cat $$name.err; \
	    diff $(srcdir)/$$name.out $$name.res; \
	  fi; \
	  rm -f $$name.res $$name.err; \
	done; \
	exit $$result

check-bad-rdf: build-rapper $(TEST_BAD_FILES)
	@set +e; result=0; \
	$(ECHO) "Testing bad Turtle fails"; \
	for test in $(TEST_BAD_FILES); do \
	  name=`basename $$test .ttl` ; \
	  baseuri=$(BASE_URI)$$test; \
	  $(ECHO) $(ECHO_N) "Checking $$test $(ECHO_C)"; \
	  $(top_builddir)/utils/rapper -q -i turtle -o ntriples file:$(srcdir)/$$test $$baseuri > $$name.res 2> $$name.err; \
	  status=$$?; \
	  if test $$status -eq 1 ; then \
	    $(ECHO) "ok"; \
	  elif test $$status -eq 2 ; then \
	    $(ECHO) "FAILED - parsing succeeded with a warning";  result=1; \
	    $(ECHO) $(top_builddir)/utils/rapper -q -i turtle -o ntriples file:$(srcdir)/$$test $$baseuri '>' $$name.res; \
	    cat $$name.res; grep Warning $$name.err; \
	  elif test $$status -eq 0 ; then \
	    $(ECHO) "FAILED - parsing succeeded but should have failed"; result=1; \
	    $(ECHO) $(top_builddir)/utils/rapper -q -i turtle -o ntriples file:$(srcdir)/$$test $$baseuri '>' $$name.res; \
	    cat $$name.res; \
	  else \
	    $(ECHO) "FAILED - parsing failed with unknown status $$status"; result=1; \
	    $(ECHO) $(top_builddir)/utils/rapper -q -i turtle -o ntriples file:$(srcdir)/$$test $$baseuri '>' $$name.res; \
	    cat $$name.err; \
	    cat $$name.res; \
	  fi; \
	  rm -f $$name.res $$name.err ; \
	done; \
	set -e; exit $$result

check-warn-rdf: build-rapper $(TEST_WARN_FILES)
	@set +e; result=0; \
	$(ECHO) "Testing Turtle with warnings"; \
	for test in $(TEST_WARN_FILES); do \
	  name=`basename $$test .ttl` ; \
	  baseuri=$(BASE_URI)$$test; \
	  $(ECHO) $(ECHO_N) "Checking $$test $(ECHO_C)"; \
	  $(top_builddir)/utils/rapper -q -i turtle -o ntriples file:$(srcdir)/$$test $$baseuri > $$name.res 2> $$name.err; \
	  status=$$?; \
	  if test $$status -eq 1 ; then \
	    $(ECHO) "FAILED - parsing failed when should have warned"; result=1; \
	    $(ECHO) $(top_builddir)/utils/rapper -q -i turtle -o ntriples file:$(srcdir)/$$test $$baseuri '>' $$name.res; \
	    cat $$name.res; grep Error $$name.err; \
	  elif test $$status -eq 2 ; then \
	    $(ECHO) "ok"; \
	  else \
	    $(ECHO) "FAILED - parsing failed with unknown status $$status"; result=1; \
	    $(ECHO) $(top_builddir)/utils/rapper -q -i turtle -o ntriples file:$(srcdir)/$$test $$baseuri '>' $$name.res; \
	    cat $$name.err; \
	    cat $$name.res; \
	  fi; \
	  rm -f $$name.res $$name.err ; \
	done; \
	set -e; exit $$result

check-turtle-serialize: build-rdfdiff build-rapper $(TEST_FILES)
	@set +e; result=0; \
	$(ECHO) "Testing turtle serialization with legal turtle"; \
	for test in $(TEST_FILES); do \
	  name=`basename $$test .ttl` ; \
	  if test $$name = rdf-schema; then \
	    baseuri=$(RDF_NS_URI); \
	  elif test $$name = rdfs-namespace; then \
	    baseuri=$(RDFS_NS_URI); \
	  else \
	    baseuri=$(BASE_URI)$$test; \
	  fi; \
	  $(ECHO) $(ECHO_N) "Checking $$test $(ECHO_C)"; \
	  $(top_builddir)/utils/rapper -q -i turtle -o turtle $(srcdir)/$$test $$baseuri > $$name-turtle.ttl 2> $$name.err; \
	  status1=$$?; \
	  $(top_builddir)/utils/rdfdiff -f turtle -u $$baseuri -t turtle $(srcdir)/$$test $$name-turtle.ttl > $$name.res 2> $$name.err; \
	  status2=$$?; \
	  if test $$status1 = 0 -a $$status2 = 0; then \
	    $(ECHO) "ok"; \
	  else \
	    $(ECHO) "FAILED"; result=1; \
	    $(ECHO) $(top_builddir)/utils/rapper -q -i turtle -o turtle $(srcdir)/$$test $$baseuri '>' $$name-turtle.ttl; \
	    $(ECHO) $(top_builddir)/utils/rdfdiff -f turtle -u $$baseuri -t turtle $(srcdir)/$$test $$name-turtle.ttl '>' $$name.res; \
	    cat $$name-turtle.ttl; cat $$name.err; \
	  fi; \
	  rm -f $$name-turtle.ttl $$name.res $$name.err; \
	done; \
	set -e; exit $$result

check-turtle-serialize-syntax: build-rapper $(TEST_SERIALIZE_FILES)
	@set +e; result=0; \
	$(ECHO) "Testing turtle exact serialization output"; \
	for test in $(TEST_SERIALIZE_FILES); do \
	  name=`basename $$test .ttl` ; \
	  if test $$name = rdf-schema; then \
	    baseuri=$(RDF_NS_URI); \
	  elif test $$name = rdfs-namespace; then \
	    baseuri=$(RDFS_NS_URI); \
	  elif test $$name = test-31 -o $$name = test-32; then \
	    baseuri=-; \
	  else \
	    baseuri=$(BASE_URI)$$test; \
	  fi; \
	  $(ECHO) $(ECHO_N) "Checking $$test $(ECHO_C)"; \
	  $(top_builddir)/utils/rapper -q -i turtle -o turtle $(srcdir)/$$test $$baseuri > $$name-turtle.ttl 2> $$name.err; \
	  cmp -s $(srcdir)/$$name-out.ttl $$name-turtle.ttl; \
	  status=$$?; \
	  if test $$status = 0; then \
	    $(ECHO) "ok"; \
	  else \
	    $(ECHO) "FAILED"; result=1; \
	    $(ECHO) $(top_builddir)/utils/rapper -q -i turtle -o turtle $(srcdir)/$$test $$baseuri '>' $$name-turtle.ttl; \
	    cat $$name.err; \
	    diff -a -u $(srcdir)/$$name-out.ttl $$name-turtle.ttl; \
	  fi; \
	  rm -f $$name-turtle.ttl $$name.err; \
	done; \
	set -e; exit $$result

check-turtle-parse-ntriples: build-rapper
	@set +e; result=0; \
	$(ECHO) "Testing Turtle parsing with N-Triples tests"; \
	NT_TEST_FILES=`unset MAKELEVEL MAKEFLAGS; cd $(NTRIPLES_DIR) && $(MAKE) print-nt-test-files | sed -e 's,^, ,' -e 's, ,$(NTRIPLES_DIR)/,g'`; \
	for test in $$NT_TEST_FILES; do \
	  name=`basename $$test .nt` ; \
	  baseuri=-; \
	  $(ECHO) $(ECHO_N) "Checking $$test $(ECHO_C)"; \
	  $(top_builddir)/utils/rapper -q -i turtle -o ntriples -n file:$(srcdir)/$$test $$baseuri > $$name.res 2> $$name.err; \
	  test_expected=$(srcdir)/$(NTRIPLES_DIR)/$$name.out; \
	  if cmp $$test_expected $$name.res >/dev/null 2>&1; then \
	    $(ECHO) "ok"; \
	  else \
	    $(ECHO) "FAILED";  result=1; \
	    $(ECHO) $(top_builddir)/utils/rapper -q -i turtle -o ntriples -n file:$(srcdir)/$$test $$baseuri '>' $$name.res; \
	    cat $$name.err; \
	    diff $$test_expected $$name.res; \
	  fi; \
	  rm -f $$name.res $$name.err; \
	done; \
	set -e; exit $$result

check-turtle-serialize-rdf: build-rdfdiff build-rapper
	@set +e; result=0; \
	RDF_TEST_FILES=`unset MAKELEVEL MAKEFLAGS; cd $(RDFXML_DIR) && $(MAKE) print-rdf-test-files | sed -e 's,^, ,' -e 's, ,$(RDFXML_DIR)/,g'`; \
	$(ECHO) "Testing turtle serialization with legal rdf/xml"; \
	for test in $$RDF_TEST_FILES $(TURTLE_SERIALIZE_RDF_FILES); do \
	  name=`basename $$test .rdf` ; \
	  baseuri=-; \
	  $(ECHO) $(ECHO_N) "Checking $$test $(ECHO_C)"; \
	  $(top_builddir)/utils/rapper -q -o turtle $(srcdir)/$$test $$baseuri > $$name-turtle.ttl 2> $$name.err; \
	  status1=$$?; \
	  if test $$test = "$(RDFXML_DIR)/ex-38.rdf"; then \
	    diff $(srcdir)/ex-38-turtle.out $$name-turtle.ttl > $$name.res 2> $$name.err; \
	    status2=$$?; \
	  else \
	    $(top_builddir)/utils/rdfdiff -t turtle $(srcdir)/$$test $$name-turtle.ttl > $$name.res 2> $$name.err; \
	    status2=$$?; \
	  fi; \
	  if test $$status1 = 0 -a $$status2 = 0; then \
	    $(ECHO) "ok"; \
	  else \
	    $(ECHO) "FAILED"; result=1; \
	    $(ECHO) $(top_builddir)/utils/rapper -q -o turtle $(srcdir)/$$test $$baseuri '>' $$name-turtle.ttl; \
	    $(ECHO) $(top_builddir)/utils/rdfdiff -t turtle $(srcdir)/$$test $$name-turtle.ttl '>' $$name.res; \
	    cat $$name.err; \
	  fi; \
	  rm -f $$name-turtle.ttl $$name.res $$name.err; \
	done; \
	set -e; exit $$result


if MAINTAINER_MODE
zip: tests.zip

tests.zip: $(ALL_TEST_FILES)
	rm -f $@
	zip $@ $(ALL_TEST_FILES)

endif
