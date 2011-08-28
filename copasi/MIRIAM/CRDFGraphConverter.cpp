// Begin CVS Header
//   $Source: /fs/turing/cvs/copasi_dev/copasi/MIRIAM/CRDFGraphConverter.cpp,v $
//   $Revision: 1.10 $
//   $Name: Build-33 $
//   $Author: shoops $
//   $Date: 2010/07/16 19:00:07 $
// End CVS Header

// Copyright (C) 2010 by Pedro Mendes, Virginia Tech Intellectual
// Properties, Inc., University of Heidelberg, and The University
// of Manchester.
// All rights reserved.

// Copyright (C) 2008 by Pedro Mendes, Virginia Tech Intellectual
// Properties, Inc., EML Research, gGmbH, University of Heidelberg,
// and The University of Manchester.
// All rights reserved.

#include "copasi.h"

#include "CRDFGraphConverter.h"
#include "CRDFParser.h"
#include "CRDFWriter.h"
#include "CRDFGraph.h"
#include "CRDFUtilities.h"

#include "utilities/CCopasiMessage.h"

// static
CRDFGraphConverter::sChange CRDFGraphConverter::SBML2CopasiChanges[] =
{
  {
    CRDFPredicate::bqbiol_encodes,
    {
      CRDFPredicate::about, CRDFPredicate::copasi_encodes, CRDFPredicate::end
    }
  },
  {
    CRDFPredicate::bqbiol_hasPart,
    {
      CRDFPredicate::about, CRDFPredicate::copasi_hasPart, CRDFPredicate::end
    }
  },
  {
    CRDFPredicate::bqbiol_hasVersion,
    {
      CRDFPredicate::about, CRDFPredicate::copasi_hasVersion, CRDFPredicate::end
    }
  },
  {
    CRDFPredicate::bqbiol_is,
    {
      CRDFPredicate::about, CRDFPredicate::copasi_is, CRDFPredicate::end
    }
  },
  {
    CRDFPredicate::bqbiol_isDescribedBy,
    {
      CRDFPredicate::about, CRDFPredicate::dcterms_bibliographicCitation, CRDFPredicate::copasi_isDescribedBy, CRDFPredicate::end
    }
  },
  {
    CRDFPredicate::bqbiol_isEncodedBy,
    {
      CRDFPredicate::about, CRDFPredicate::copasi_isEncodedBy, CRDFPredicate::end
    }
  },
  {
    CRDFPredicate::bqbiol_isHomologTo,
    {
      CRDFPredicate::about, CRDFPredicate::copasi_isHomologTo, CRDFPredicate::end
    }
  },
  {
    CRDFPredicate::bqbiol_isPartOf,
    {
      CRDFPredicate::about, CRDFPredicate::copasi_isPartOf, CRDFPredicate::end
    }
  },
  {
    CRDFPredicate::bqbiol_isVersionOf,
    {
      CRDFPredicate::about, CRDFPredicate::copasi_isVersionOf, CRDFPredicate::end
    }
  },
  {
    CRDFPredicate::bqbiol_occursIn,
    {
      CRDFPredicate::about, CRDFPredicate::copasi_occursIn, CRDFPredicate::end
    }
  },
  {
    CRDFPredicate::bqmodel_is,
    {
      CRDFPredicate::about, CRDFPredicate::copasi_is, CRDFPredicate::end
    }
  },
  {
    CRDFPredicate::bqmodel_isDescribedBy,
    {
      CRDFPredicate::about, CRDFPredicate::dcterms_bibliographicCitation, CRDFPredicate::copasi_isDescribedBy, CRDFPredicate::end
    }
  },
  {
    CRDFPredicate::dc_creator,
    {
      CRDFPredicate::about, CRDFPredicate::dcterms_creator, CRDFPredicate::end
    },
  },
  {
    CRDFPredicate::end,
    {
      CRDFPredicate::end
    }
  }
};

// static
bool CRDFGraphConverter::SBML2Copasi(std::string & XML)
{
  // Fix the the broken SBML RDF
  if (CRDFUtilities::fixSBMLRdf(XML))
    CCopasiMessage(CCopasiMessage::WARNING_FILTERED, MCSBML + 75);

  // Create the RDF graph
  CRDFGraph * pGraph = CRDFParser::graphFromXml(XML);

  if (pGraph == NULL)
    return false;

  // Convert the graph.
  bool success = convert(pGraph, SBML2CopasiChanges);

  // Serialize the graph
  pGraph->clean();
  pGraph->updateNamespaces();

  XML = CRDFWriter::xmlFromGraph(pGraph);

  if (pGraph != NULL)
    {
      delete pGraph;
      pGraph = NULL;
    }

  // It is possible that the graph is still corrupt due secondary errors.
  // Another parse and write should take care of this.
  unsigned C_INT32 Size = CCopasiMessage::size();

  pGraph = CRDFParser::graphFromXml(XML);

  if (pGraph == NULL)
    return false;

  XML = CRDFWriter::xmlFromGraph(pGraph);

  if (pGraph != NULL)
    {
      delete pGraph;
      pGraph = NULL;
    }

  // Remove error messages created by the second parse as these are secondary.
  while (CCopasiMessage::size() > Size)
    CCopasiMessage::getLastMessage();

  return success;
}

// static
bool CRDFGraphConverter::convert(CRDFGraph * pGraph, const CRDFGraphConverter::sChange * changes)
{
  bool success = true;

  std::set< CRDFTriplet> Triplets;
  std::set< CRDFTriplet>::iterator it;
  std::set< CRDFTriplet>::iterator end;

  const sChange * pChange = changes;

  // We iterate over each predicate which needs to be changed.
  for (; pChange->Source != CRDFPredicate::end; ++pChange)
    {
      std::set< CRDFTriplet > Failed;

      // Create the new path
      CRDFPredicate::Path NewPath;
      const CRDFPredicate::ePredicateType * pNewPath = pChange->Target;

      while (*pNewPath != CRDFPredicate::end)
        NewPath.push_back(*pNewPath++);

      // Each change may break the triplets, i.e. we need to refresh every time
      while ((Triplets = pGraph->getTriplets(pChange->Source, false)).size() > Failed.size())
        {
          // Skip all failed triplets
          for (it = Triplets.begin(), end = Triplets.end(); it != end; ++it)
            if (Failed.find(*it) == Failed.end())
              break;

          if (!convert(pGraph, *it, NewPath))
            {
              Failed.insert(*it);
              success = false;
            }
        }
    }

  return success;
}

// static
bool CRDFGraphConverter::convert(CRDFGraph * pGraph,
                                 const CRDFTriplet & triplet,
                                 const CRDFPredicate::Path & newPath)
{
  CRDFPredicate::Path CurrentPath = triplet.pObject->getPath();

  unsigned C_INT32 SubPathIndex = C_INVALID_INDEX;

  while (SubPathIndex == C_INVALID_INDEX && SubPathIndex > 0)
    {
      // We differ at least at the final predicate.
      CurrentPath.pop_back();
      SubPathIndex = CRDFPredicate::getSubPathIndex(newPath, CurrentPath);
    }

  if (SubPathIndex == 0)
    return false;

  CurrentPath = triplet.pObject->getPath();

  bool success = true;

  CRDFTriplet Triplet;

  if (CurrentPath.size() < newPath.size())
    {
      // We need to create the additional node pointed to by newPath[SubPathIndex]
      // Verify that the node to be created is a blank node.
      CRDFPredicate::AllowedLocationList List =
        CRDFPredicate::getAllowedLocationList(newPath[SubPathIndex]);
      CRDFPredicate::AllowedLocationList::const_iterator it = List.begin();
      CRDFPredicate::AllowedLocationList::const_iterator end = List.end();

      for (; it != end; ++it)
        {
          if (it->Type == CRDFObject::BLANK_NODE &&
              CRDFPredicate::getSubPathIndex(newPath, it->Location))
            break;
        }

      if (it == end)
        return false;

      // We are now sure that the new predicate points to a blank node.
      CRDFObject Object;
      Object.setType(CRDFObject::BLANK_NODE);
      Object.setBlankNodeId(pGraph->generatedNodeId());

      Triplet = pGraph->addTriplet(triplet.pSubject->getSubject(),
                                   CRDFPredicate::getURI(newPath[SubPathIndex]),
                                   Object);

      if (Triplet)
        Triplet = pGraph->moveTriplet(Triplet.pObject, triplet);
    }
  else if (CurrentPath.size() > newPath.size())
    {
      // TODO This needs to be implemented
      // Removing an edge will always result in loosing information!
    }
  else
    {
      // We just have to rename the predicate
      success &= triplet.pSubject->addEdge(newPath[SubPathIndex], triplet.pObject);
      triplet.pSubject->removeEdge(CurrentPath[SubPathIndex], triplet.pObject);

      if (success)
        {
          Triplet = triplet;
          Triplet.Predicate = newPath[SubPathIndex];
        }
    }

  if (!Triplet)
    return false;

  // Check whether we are finished
  if (SubPathIndex == newPath.size() - 1)
    return success;

  return convert(pGraph, Triplet, newPath);
}
