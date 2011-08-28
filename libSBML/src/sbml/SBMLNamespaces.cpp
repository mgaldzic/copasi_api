/**
 * @file    SBMLNamespaces.cpp
 * @brief   SBMLNamespaces class to store level/version and namespace 
 * @author  Sarah Keating
 *
 * $Id: SBMLNamespaces.cpp 11635 2010-08-03 04:05:22Z mhucka $
 * $HeadURL: https://sbml.svn.sourceforge.net/svnroot/sbml/trunk/libsbml/src/sbml/SBMLNamespaces.cpp $
 *
 *<!---------------------------------------------------------------------------
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
 *----------------------------------------------------------------------- -->
 */

#include <sbml/SBMLNamespaces.h>

/** @cond doxygen-ignored */

using namespace std;

/** @endcond */

LIBSBML_CPP_NAMESPACE_BEGIN

SBMLNamespaces::SBMLNamespaces(unsigned int level, unsigned int version)
{
  /* default is l3v1 */
  mLevel = level;
  mVersion = version;
  mNamespaces = new XMLNamespaces();

  switch (level)
  {
  case 1:
    mNamespaces->add(SBML_XMLNS_L1);
    break;
  case 2:
    switch (version)
    {
    case 1:
      mNamespaces->add(SBML_XMLNS_L2V1);
      break;
    case 2:
      mNamespaces->add(SBML_XMLNS_L2V2);
      break;
    case 3:
      mNamespaces->add(SBML_XMLNS_L2V3);
      break;
    case 4:
    default:
      mNamespaces->add(SBML_XMLNS_L2V4);
      break;
    }
    break;
  case 3:
  default:
    switch (version)
    {
    case 1:
    default:
      mNamespaces->add(SBML_XMLNS_L3V1);
      break;
    }
    break;
  }
}

SBMLNamespaces::~SBMLNamespaces()
{
  if (mNamespaces)
    delete mNamespaces;
}


/*
 * Copy constructor; creates a copy of a SBMLNamespaces.
 */
SBMLNamespaces::SBMLNamespaces(const SBMLNamespaces& orig) :
          mLevel    (orig.mLevel)
      ,   mVersion  (orig.mVersion) 
{
  if(orig.mNamespaces)
    this->mNamespaces = 
           new XMLNamespaces(*const_cast<SBMLNamespaces&>(orig).mNamespaces);
  else
    this->mNamespaces = 0;
}



/*
 * Assignment operator for SBMLNamespaces.
 */
SBMLNamespaces&
SBMLNamespaces::operator=(const SBMLNamespaces& orig)
{
  if (&orig != this)
  {
    mLevel = orig.mLevel;
    mVersion = orig.mVersion;
    delete this->mNamespaces;
    if(orig.mNamespaces)
      this->mNamespaces = 
            new XMLNamespaces(*const_cast<SBMLNamespaces&>(orig).mNamespaces);
    else
      this->mNamespaces = 0;
  }

  return *this;
}



/*
 * Creates and returns a deep copy of this SBMLNamespaces.
 */
SBMLNamespaces *
SBMLNamespaces::clone () const
{
  return new SBMLNamespaces(*this);
}


std::string 
SBMLNamespaces::getSBMLNamespaceURI(unsigned int level,
                                 unsigned int version)
{
  switch (level)
  {
  case 1:
    return SBML_XMLNS_L1;
    break;
  case 3:
    switch(version)
    {
    case 1:
    default:
      return SBML_XMLNS_L3V1;
      break;
    }
    break;
  case 2:
  default:
    switch (version)
    {
    case 1:
      return SBML_XMLNS_L2V1;
      break;
    case 2:
      return SBML_XMLNS_L2V2;
      break;
    case 3:
      return SBML_XMLNS_L2V3;
      break;
    case 4:
    default:
      return SBML_XMLNS_L2V4;
      break;
    }
    break;
  }
}


unsigned int 
SBMLNamespaces::getLevel()
{
  return mLevel;
}


unsigned int 
SBMLNamespaces::getLevel() const
{
  return mLevel;
}


unsigned int 
SBMLNamespaces::getVersion()
{
  return mVersion;
}


unsigned int 
SBMLNamespaces::getVersion() const
{
  return mVersion;
}


XMLNamespaces * 
SBMLNamespaces::getNamespaces()
{
  return mNamespaces;
}


const XMLNamespaces * 
SBMLNamespaces::getNamespaces() const
{
  return mNamespaces;
}


void
SBMLNamespaces::addNamespaces(XMLNamespaces * xmlns)
{
  if (xmlns == NULL)
    return;

  /* check whether the namespace already exists
   * add if it does not
   */
  for (int i = 0; i < xmlns->getLength(); i++)
  {
    if (!(mNamespaces->hasNS(xmlns->getURI(i), xmlns->getPrefix(i))))
    {
      mNamespaces->add(xmlns->getURI(i), xmlns->getPrefix(i));
    }
  }
}

/** @cond doxygen-libsbml-internal */
void 
SBMLNamespaces::setLevel(unsigned int level)
{
  mLevel = level;
}


void 
SBMLNamespaces::setVersion(unsigned int version)
{
  mVersion = version;
}


void 
SBMLNamespaces::setNamespaces(XMLNamespaces * xmlns)
{
  delete mNamespaces;
  if (xmlns)
    mNamespaces = xmlns->clone();
  else
    mNamespaces = NULL;
}
/** @endcond */
/** @cond doxygen-c-only */

/**
 * Creates a new SBMLNamespaces_t structure corresponding to the given SBML
 * @p level and @p version.
 *
 * SBMLNamespaces objects are used in libSBML to communicate SBML Level
 * and Version data between constructors and other methods.  The
 * SBMLNamespaces object class tracks 3-tuples (triples) consisting of
 * SBML Level, Version, and the corresponding SBML XML namespace.  Most
 * constructors for SBML objects in libSBML take a SBMLNamespaces object
 * as an argument, thereby allowing the constructor to produce the proper
 * combination of attributes and other internal data structures for the
 * given SBML Level and Version.
 *
 * The plural name "SBMLNamespaces" is not a mistake, because in SBML
 * Level&nbsp;3, objects may have extensions added by Level&nbsp;3
 * packages used by a given model; however, until the introduction of
 * SBML Level&nbsp;3, the SBMLNamespaces object only records one SBML
 * Level/Version/namespace combination at a time.
 *
 * @param level the SBML level
 * @param version the SBML version
 *
 * @return SBMLNamespaces_t structure created
 * 
 * @docnote The native C++ implementation of this method defines a
 * default argument value.  In the documentation generated for different
 * libSBML language bindings, you may or may not see corresponding
 * arguments in the method declarations.  For example, in Java, a default
 * argument is handled by declaring two separate methods, with one of
 * them having the argument and the other one lacking the argument.
 * However, the libSBML documentation will be @em identical for both
 * methods.  Consequently, if you are reading this and do not see an
 * argument even though one is described, please look for descriptions of
 * other variants of this method near where this one appears in the
 * documentation.
 */

LIBSBML_EXTERN
SBMLNamespaces_t *
SBMLNamespaces_create(unsigned int level, unsigned int version)
{
  return new SBMLNamespaces(level, version);
}


/**
 * Get the SBML Level of this SBMLNamespaces_t structure.
 *
 * @param sbmlns the SBMLNamespaces_t structure to query
 *
 * @return the SBML Level of this SBMLNamespaces_t structure.
 */
LIBSBML_EXTERN
unsigned int
SBMLNamespaces_getLevel(SBMLNamespaces_t *sbmlns)
{
  return sbmlns->getLevel();
}


/**
 * Get the SBML Version of this SBMLNamespaces_t structure.
 *
 * @param sbmlns the SBMLNamespaces_t structure to query
 *
 * @return the SBML Version of this SBMLNamespaces_t structure.
 */
LIBSBML_EXTERN
unsigned int
SBMLNamespaces_getVersion(SBMLNamespaces_t *sbmlns)
{
  return sbmlns->getVersion();
}


/**
 * Get the SBML Version of this SBMLNamespaces_t structure.
 *
 * @param sbmlns the SBMLNamespaces_t structure to query
 *
 * @return the XMLNamespaces_t structure of this SBMLNamespaces_t structure.
 */
LIBSBML_EXTERN
XMLNamespaces_t *
SBMLNamespaces_getNamespaces(SBMLNamespaces_t *sbmlns)
{
  return sbmlns->getNamespaces();
}


/**
 * Returns a string representing the SBML XML namespace for the 
 * given @p level and @p version of SBML.
 *
 * @param level the SBML level
 * @param version the SBML version
 *
 * @return a string representing the SBML namespace that reflects the
 * SBML Level and Version specified.
 */
LIBSBML_EXTERN
const char *
SBMLNamespaces_getSBMLNamespaceURI(unsigned int level, unsigned int version)
{
  return SBMLNamespaces::getSBMLNamespaceURI(level, version).c_str();
}


/**
 * Add the XML namespaces list to the set of namespaces
 * within this SBMLNamespaces_t structure.
 * 
 * @param sbmlns the SBMLNamespaces_t structure to add to
 * @param xmlns the XML namespaces to be added.
 */
LIBSBML_EXTERN
void
SBMLNamespaces_addNamespaces(SBMLNamespaces_t *sbmlns,
                             XMLNamespaces_t * xmlns)
{
  sbmlns->addNamespaces(xmlns);
}

/** @endcond */
LIBSBML_CPP_NAMESPACE_END

