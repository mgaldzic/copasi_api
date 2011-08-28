// Begin CVS Header
//   $Source: /fs/turing/cvs/copasi_dev/copasi/MIRIAM/CModelMIRIAMInfo.cpp,v $
//   $Revision: 1.35.2.1 $
//   $Name: Build-33 $
//   $Author: shoops $
//   $Date: 2010/09/28 16:09:40 $
// End CVS Header

// Copyright (C) 2010 by Pedro Mendes, Virginia Tech Intellectual
// Properties, Inc., University of Heidelberg, and The University
// of Manchester.
// All rights reserved.

// Copyright (C) 2008 by Pedro Mendes, Virginia Tech Intellectual
// Properties, Inc., EML Research, gGmbH, University of Heidelberg,
// and The University of Manchester.
// All rights reserved.

// Copyright (C) 2001 - 2007 by Pedro Mendes, Virginia Tech Intellectual
// Properties, Inc. and EML Research, gGmbH.
// All rights reserved.

#include <iostream>
#include <fstream>
#include <string>

#include "copasi.h"

#include "CModelMIRIAMInfo.h"
#include "CRDFWriter.h"
#include "CRDFLiteral.h"
#include "CRDFParser.h"
#include "CConstants.h"
#include "CRDFObject.h"
#include "CRDFPredicate.h"
#include "CRDFGraph.h"

#include "model/CModelValue.h"
#include "model/CEvent.h"
#include "model/CReaction.h"
#include "function/CFunction.h"
#include "report/CKeyFactory.h"
#include "copasi/report/CCopasiRootContainer.h"

CMIRIAMInfo::CMIRIAMInfo() :
    CCopasiContainer("CMIRIAMInfoObject", NULL, "CMIRIAMInfo"),
    mKey(""),
    mCreators("Creators", this),
    mReferences("References", this),
    mModifications("Modifieds", this),
    mBiologicalDescriptions("BiologicalDescriptions", this),
    mCreatedObj(),
    mpRDFGraph(NULL),
    mTriplet(NULL, CRDFPredicate::about, NULL),
    mCreated()
{}

CMIRIAMInfo::~CMIRIAMInfo()
{pdelete(mpRDFGraph);}

CRDFGraph * CMIRIAMInfo::getRDFGraph()
{return mpRDFGraph;}

const CCopasiVector <CCreator> & CMIRIAMInfo::getCreators() const
{return mCreators;}

CCreator* CMIRIAMInfo::createCreator(const std::string & /* objectName */)
{
  const CRDFSubject & Subject = mpRDFGraph->getAboutNode()->getSubject();
  CRDFObject Object;
  Object.setType(CRDFObject::BLANK_NODE);
  std::string Id = mpRDFGraph->generatedNodeId();
  Object.setBlankNodeId(Id);

  CRDFTriplet Triplet =
    mpRDFGraph->addTriplet(Subject,
                           CRDFPredicate::getURI(CRDFPredicate::dcterms_creator),
                           Object);

  if (!Triplet)
    return NULL;

  CCreator * pCreator = new CCreator(Triplet);

  if (!mCreators.add(pCreator, true))
    {
      delete pCreator;
      return NULL;
    }

  return pCreator;
}

bool CMIRIAMInfo::removeCreator(int position)
{
  CCreator * pCreator = mCreators[position];

  if (!pCreator)
    return false;

  const CRDFTriplet & Triplet = pCreator->getTriplet();

  mpRDFGraph->removeTriplet(Triplet.pSubject,
                            CRDFPredicate::getURI(Triplet.Predicate),
                            Triplet.pObject);

  return mCreators.remove(pCreator);
}

void CMIRIAMInfo::loadCreators()
{
  mCreators.cleanup();

  CRDFPredicate::ePredicateType Predicates[] =
  {
    CRDFPredicate::dcterms_creator,
    CRDFPredicate::dc_creator,
    CRDFPredicate::end
  };

  CRDFPredicate::Path Path = mTriplet.pObject->getPath();
  std::set< CRDFTriplet > Triples;

  CRDFPredicate::ePredicateType * pPredicate = Predicates;
  std::set< CRDFTriplet >::iterator it;
  std::set< CRDFTriplet >::iterator end;

  for (; *pPredicate != CRDFPredicate::end; ++pPredicate)
    {
      Triples =
        mTriplet.pObject->getDescendantsWithPredicate(*pPredicate);
      it = Triples.begin();
      end = Triples.end();

      for (; it != end; ++it)
        mCreators.add(new CCreator(*it), true);
    }

  return;
}

const CCopasiVector <CReference> & CMIRIAMInfo::getReferences() const
{return mReferences;}

CReference* CMIRIAMInfo::createReference(const std::string & /* objectName */)
{
  const CRDFSubject & Subject = mpRDFGraph->getAboutNode()->getSubject();
  CRDFObject Object;
  Object.setType(CRDFObject::BLANK_NODE);
  std::string Id = mpRDFGraph->generatedNodeId();
  Object.setBlankNodeId(Id);

  CRDFTriplet Triplet =
    mpRDFGraph->addTriplet(Subject,
                           CRDFPredicate::getURI(CRDFPredicate::dcterms_bibliographicCitation),
                           Object);

  if (!Triplet)
    return NULL;

  CReference * pReference = new CReference(Triplet);

  if (!mReferences.add(pReference, true))
    {
      delete pReference;
      return NULL;
    }

  return pReference;
}

bool CMIRIAMInfo::removeReference(int position)
{
  CReference * pReference = mReferences[position];

  if (!pReference)
    return false;

  const CRDFTriplet & Triplet = pReference->getTriplet();

  mpRDFGraph->removeTriplet(Triplet.pSubject,
                            CRDFPredicate::getURI(Triplet.Predicate),
                            Triplet.pObject);

  return mReferences.remove(pReference);
}

void CMIRIAMInfo::loadReferences()
{
  mReferences.cleanup();

  CRDFPredicate::ePredicateType Predicates[] =
  {
    CRDFPredicate::dcterms_bibliographicCitation,
    CRDFPredicate::bqbiol_isDescribedBy,
    CRDFPredicate::bqmodel_isDescribedBy,
    CRDFPredicate::end
  };

  CRDFPredicate::Path Path = mTriplet.pObject->getPath();
  std::set< CRDFTriplet > Triples;

  CRDFPredicate::ePredicateType * pPredicate = Predicates;
  std::set< CRDFTriplet >::iterator it;
  std::set< CRDFTriplet >::iterator end;

  for (; *pPredicate != CRDFPredicate::end; ++pPredicate)
    {
      Triples = mTriplet.pObject->getDescendantsWithPredicate(*pPredicate);
      it = Triples.begin();
      end = Triples.end();

      for (; it != end; ++it)
        mReferences.add(new CReference(*it), true);
    }

  return;
}

const std::string CMIRIAMInfo::getCreatedDT() const
{
  if (!mCreated)
    return "";

  return mCreated.pObject->getFieldValue(CRDFPredicate::dcterms_W3CDTF);
}

void CMIRIAMInfo::setCreatedDT(const std::string& dt)
{
  std::string Date = dt;

  if (Date == "0000-00-00T00:00:00")
    Date = ""; // This causes deletion of the edge

  if (!mCreated)
    {
      const CRDFSubject & Subject = mTriplet.pObject->getSubject();
      CRDFObject Object;
      Object.setType(CRDFObject::BLANK_NODE);
      std::string Id = mpRDFGraph->generatedNodeId();
      Object.setBlankNodeId(Id);

      mCreated = mpRDFGraph->addTriplet(Subject,
                                        CRDFPredicate::getURI(CRDFPredicate::dcterms_created),
                                        Object);
      // Debugging
      assert(!mCreated == false);
    }

  mCreated.pObject->setFieldValue(Date, CRDFPredicate::dcterms_W3CDTF, mCreated.pObject->getPath());
}

const CCopasiVector <CModification> & CMIRIAMInfo::getModifications() const
{return mModifications;}

CModification * CMIRIAMInfo::createModification(const std::string& dateTime)
{
  const CRDFSubject & Subject = mpRDFGraph->getAboutNode()->getSubject();
  CRDFObject Object;
  Object.setType(CRDFObject::BLANK_NODE);
  std::string Id = mpRDFGraph->generatedNodeId();
  Object.setBlankNodeId(Id);

  CRDFTriplet Triplet =
    mpRDFGraph->addTriplet(Subject,
                           CRDFPredicate::getURI(CRDFPredicate::dcterms_modified),
                           Object);

  if (!Triplet)
    return NULL;

  CModification * pModification = new CModification(Triplet);

  if (dateTime.size())
    pModification->setDate(dateTime);

  if (!mModifications.add(pModification, true))
    {
      delete pModification;
      return NULL;
    }

  return pModification;
}

bool CMIRIAMInfo::removeModification(int position)
{
  CModification * pModified = mModifications[position];

  if (!pModified)
    return false;

  const CRDFTriplet & Triplet = pModified->getTriplet();

  mpRDFGraph->removeTriplet(Triplet.pSubject,
                            CRDFPredicate::getURI(Triplet.Predicate),
                            Triplet.pObject);

  return mModifications.remove(pModified);
}

void CMIRIAMInfo::loadModifications()
{
  mModifications.cleanup();

  std::set< CRDFTriplet > Triples =
    mTriplet.pObject->getDescendantsWithPredicate(CRDFPredicate::dcterms_modified);
  std::set< CRDFTriplet >::iterator it = Triples.begin();
  std::set< CRDFTriplet >::iterator end = Triples.end();

  for (; it != end; ++it)
    mModifications.add(new CModification(*it), true);

  return;
}

const CCopasiVector <CBiologicalDescription> & CMIRIAMInfo::getBiologicalDescriptions() const
{return mBiologicalDescriptions;}

CBiologicalDescription* CMIRIAMInfo::createBiologicalDescription()
{
  const CRDFSubject & Subject = mpRDFGraph->getAboutNode()->getSubject();
  CRDFObject Object;
  Object.setType(CRDFObject::RESOURCE);
  Object.setResource("", false);

  CRDFTriplet Triplet = mpRDFGraph->addTriplet(Subject, std::string("---"), Object);

  if (!Triplet)
    return NULL;

  CBiologicalDescription * pBiologicalDescription =
    new CBiologicalDescription(Triplet);

  if (!mBiologicalDescriptions.add(pBiologicalDescription, true))
    {
      delete pBiologicalDescription;
      return NULL;
    }

  return pBiologicalDescription;
}

bool CMIRIAMInfo::removeBiologicalDescription(int position)
{
  CBiologicalDescription * pBiologicalDescription =
    mBiologicalDescriptions[position];

  if (!pBiologicalDescription)
    return false;

  const CRDFTriplet & Triplet = pBiologicalDescription->getTriplet();

  mpRDFGraph->removeTriplet(Triplet.pSubject,
                            Triplet.Predicate,//CRDFPredicate::getURI(Triplet.Predicate),
                            Triplet.pObject);

  return mBiologicalDescriptions.remove(pBiologicalDescription);
}

void CMIRIAMInfo::loadBiologicalDescriptions()
{
  mBiologicalDescriptions.cleanup();

  CRDFPredicate::ePredicateType Predicates[] =
  {
    CRDFPredicate::copasi_encodes,
    CRDFPredicate::copasi_hasPart,
    CRDFPredicate::copasi_hasVersion,
    CRDFPredicate::copasi_is,
    CRDFPredicate::copasi_isEncodedBy,
    CRDFPredicate::copasi_isHomologTo,
    CRDFPredicate::copasi_isPartOf,
    CRDFPredicate::copasi_isVersionOf,
    CRDFPredicate::copasi_occursIn,
    CRDFPredicate::bqbiol_encodes,
    CRDFPredicate::bqbiol_hasPart,
    CRDFPredicate::bqbiol_hasVersion,
    CRDFPredicate::bqbiol_is,
    CRDFPredicate::bqbiol_isEncodedBy,
    CRDFPredicate::bqbiol_isHomologTo,
    CRDFPredicate::bqbiol_isPartOf,
    CRDFPredicate::bqbiol_isVersionOf,
    CRDFPredicate::bqbiol_occursIn,
    CRDFPredicate::bqmodel_is,
    CRDFPredicate::end
  };

  CRDFPredicate::Path Path = mTriplet.pObject->getPath();
  std::set< CRDFTriplet > Triples;

  CRDFPredicate::ePredicateType * pPredicate = Predicates;
  std::set< CRDFTriplet >::iterator it;
  std::set< CRDFTriplet >::iterator end;

  for (; *pPredicate != CRDFPredicate::end; ++pPredicate)
    {
      Triples = mTriplet.pObject->getDescendantsWithPredicate(*pPredicate);
      it = Triples.begin();
      end = Triples.end();

      for (; it != end; ++it)
        mBiologicalDescriptions.add(new CBiologicalDescription(*it), true);
    }
}

void CMIRIAMInfo::load(const std::string& key)
{
  pdelete(mpRDFGraph);

  mKey = key;
  CCopasiObject * pCopasiObject = dynamic_cast< CCopasiObject * >(CCopasiRootContainer::getKeyFactory()->get(mKey));

  if (pCopasiObject != NULL)
    {
      const std::string * pMiriamAnnotation = NULL;

      if (dynamic_cast< CModelEntity * >(pCopasiObject))
        pMiriamAnnotation = &static_cast< CModelEntity * >(pCopasiObject)->getMiriamAnnotation();
      else if (dynamic_cast< CEvent * >(pCopasiObject))
        pMiriamAnnotation = &static_cast< CEvent * >(pCopasiObject)->getMiriamAnnotation();
      else if (dynamic_cast< CReaction * >(pCopasiObject))
        pMiriamAnnotation = &static_cast< CReaction * >(pCopasiObject)->getMiriamAnnotation();
      else if (dynamic_cast< CFunction * >(pCopasiObject))
        pMiriamAnnotation = &static_cast< CFunction * >(pCopasiObject)->getMiriamAnnotation();

      if (pMiriamAnnotation && *pMiriamAnnotation != "")
        mpRDFGraph = CRDFParser::graphFromXml(*pMiriamAnnotation);
    }

  if (mpRDFGraph == NULL)
    mpRDFGraph = new CRDFGraph;

  // We make sure that we always have an about node.

  if (pCopasiObject != NULL)
    mTriplet.pObject = mpRDFGraph->createAboutNode(pCopasiObject->getKey());
  else
    mTriplet.pObject = mpRDFGraph->createAboutNode("");

  // Load the created date if set;
  CRDFPredicate::Path Path = mTriplet.pObject->getPath();
  std::set< CRDFTriplet > Triples =
    mTriplet.pObject->getDescendantsWithPredicate(CRDFPredicate::dcterms_created);

  if (Triples.size() > 0)
    mCreated = *Triples.begin();
  else
    mCreated = CRDFTriplet(); // This is an invalid triplet, i.e., !mCreated is true.

  loadCreators();
  loadReferences();
  loadModifications();
  loadBiologicalDescriptions();

  return;
}

bool CMIRIAMInfo::save()
{
  CCopasiObject * pCopasiObject = dynamic_cast< CCopasiObject * >(CCopasiRootContainer::getKeyFactory()->get(mKey));

  if (pCopasiObject && mpRDFGraph)
    {
      mpRDFGraph->clean();
      mpRDFGraph->updateNamespaces();

      std::string XML = CRDFWriter::xmlFromGraph(mpRDFGraph);

      CModelEntity * pEntity = NULL;
      CEvent * pEvent = NULL;
      CReaction * pReaction = NULL;
      CFunction * pFunction = NULL;

      if ((pEntity = dynamic_cast< CModelEntity * >(pCopasiObject)) != NULL)
        pEntity->setMiriamAnnotation(XML, pEntity->getKey(), pEntity->getKey());
      else if ((pEvent = dynamic_cast< CEvent * >(pCopasiObject)) != NULL)
        pEvent->setMiriamAnnotation(XML, pEvent->getKey(), pEvent->getKey());
      else if ((pReaction = dynamic_cast< CReaction * >(pCopasiObject)) != NULL)
        pReaction->setMiriamAnnotation(XML, pReaction->getKey(), pReaction->getKey());
      else if ((pFunction = dynamic_cast< CFunction * >(pCopasiObject)) != NULL)
        pFunction->setMiriamAnnotation(XML, pFunction->getKey(), pFunction->getKey());
      else
        return false;

      return true;
    }

  return false;
}

const std::string & CMIRIAMInfo::getKey() const
{return mKey;}
