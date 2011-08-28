// Begin CVS Header
//   $Source: /fs/turing/cvs/copasi_dev/copasi/layout/CLayout.cpp,v $
//   $Revision: 1.17 $
//   $Name: Build-31 $
//   $Author: shoops $
//   $Date: 2009/10/27 16:52:20 $
// End CVS Header

// Copyright (C) 2008 by Pedro Mendes, Virginia Tech Intellectual
// Properties, Inc., EML Research, gGmbH, University of Heidelberg,
// and The University of Manchester.
// All rights reserved.

// Copyright (C) 2001 - 2007 by Pedro Mendes, Virginia Tech Intellectual
// Properties, Inc. and EML Research, gGmbH.
// All rights reserved.

#define USE_LAYOUT 1

#include "iostream"
#include "sbml/layout/Layout.h"

#include "copasi.h"

#include "CLayout.h"

#include "report/CKeyFactory.h"
#include "sbml/CSBMLExporter.h"
#include "copasi/report/CCopasiRootContainer.h"

CLayout::CLayout(const std::string & name,
                 const CCopasiContainer * pParent)
    : CLBase(),
    CCopasiContainer(name, pParent, "Layout"),
    mKey(CCopasiRootContainer::getKeyFactory()->add("Layout", this)),
    mDimensions(),
    mvCompartments("ListOfCompartmentGlyphs", this),
    mvMetabs("ListOfMetaboliteGlyphs", this),
    mvReactions("ListOfReactionGlyphs", this),
    mvLabels("ListOfTextGlyphs", this),
    mvGraphicalObjects("ListOfGraphicalObjects", this)
{}

CLayout::CLayout(const CLayout & src,
                 const CCopasiContainer * pParent)
    : CLBase(src),
    CCopasiContainer(src, pParent),
    mKey(CCopasiRootContainer::getKeyFactory()->add("Layout", this)),
    mDimensions(src.mDimensions),
    mvCompartments(src.mvCompartments, this),
    mvMetabs(src.mvMetabs, this),
    mvReactions(src.mvReactions, this),
    mvLabels(src.mvLabels, this),
    mvGraphicalObjects(src.mvGraphicalObjects, this)
{
  //TODO references from one glyph to another have to be reconstructed after
  //     copying. This applies to Labels and species reference glyphs
}

CLayout::CLayout(const Layout & sbml,
                 std::map<std::string, std::string> & layoutmap,
                 const CCopasiContainer * pParent)
    : CLBase(sbml),
    CCopasiContainer(sbml.getId(), pParent, "Layout"),
    mKey(CCopasiRootContainer::getKeyFactory()->add("Layout", this)),
    mDimensions(*sbml.getDimensions()),
    mvCompartments("ListOfCompartmentGlyphs", this),
    mvMetabs("ListOfMetaboliteGlyphs", this),
    mvReactions("ListOfReactionGlyphs", this),
    mvLabels("ListOfTextGlyphs", this),
    mvGraphicalObjects("ListOfGraphicalObjects", this)
{
  //add the copasi key to the map
  layoutmap[sbml.getId()] = mKey;
}

CLayout::~CLayout()
{
  CCopasiRootContainer::getKeyFactory()->remove(mKey);
}

void CLayout::addCompartmentGlyph(CLCompartmentGlyph * glyph)
{
  if (glyph)
    mvCompartments.add(glyph, true); //true means vector takes ownership
}

void CLayout::addMetaboliteGlyph(CLMetabGlyph * glyph)
{
  if (glyph)
    mvMetabs.add(glyph, true); //true means vector takes ownership
}

void CLayout::addReactionGlyph(CLReactionGlyph * glyph)
{
  if (glyph)
    mvReactions.add(glyph, true); //true means vector takes ownership
}

void CLayout::addTextGlyph(CLTextGlyph * glyph)
{
  if (glyph)
    mvLabels.add(glyph, true); //true means vector takes ownership
}

void CLayout::addGraphicalObject(CLGraphicalObject * glyph)
{
  if (glyph)
    mvGraphicalObjects.add(glyph, true); //true means vector takes ownership
}

std::ostream & operator<<(std::ostream &os, const CLayout & l)
{
  C_INT32 i, imax;

  os << "Layout  \"" << l.getObjectName() << "\" " << l.mDimensions << "\n\n";

  imax = l.mvCompartments.size();

  if (imax)
    {
      os << "List of compartment glyphs: \n\n";

      for (i = 0; i < imax; ++i)
        os << *l.mvCompartments[i];
    }

  imax = l.mvMetabs.size();

  if (imax)
    {
      os << "\nList of species glyphs: \n\n";

      for (i = 0; i < imax; ++i)
        os << *l.mvMetabs[i];
    }

  imax = l.mvReactions.size();

  if (imax)
    {
      os << "\nList of reaction glyphs: \n\n";

      for (i = 0; i < imax; ++i)
        os << *l.mvReactions[i];
    }

  imax = l.mvLabels.size();

  if (imax)
    {
      os << "\nList of labels: \n\n";

      for (i = 0; i < imax; ++i)
        os << *l.mvLabels[i];
    }

  imax = l.mvGraphicalObjects.size();

  if (imax)
    {
      os << "\nList of graphical objects: \n\n";

      for (i = 0; i < imax; ++i)
        os << *l.mvGraphicalObjects[i];
    }

  return os;
}

void CLayout::print(std::ostream * os) const
{*os << *this;}

void CLayout::exportToDotFile(std::ostream & os) const
{
  os << "digraph G {\n";

  //species glyphs
  unsigned C_INT32 i, imax = mvMetabs.size();

  for (i = 0; i < imax; ++i)
    {
      writeDotNode(os, mvMetabs[i]->getKey(), mvMetabs[i]->getModelObjectDisplayName());
    }

  //reaction glyphs
  imax = mvReactions.size();

  for (i = 0; i < imax; ++i)
    {
      writeDotNode(os, mvReactions[i]->getKey() + "_S", "", 1);
      writeDotNode(os, mvReactions[i]->getKey() + "_P", "", 1);
      writeDotEdge(os, mvReactions[i]->getKey() + "_S", mvReactions[i]->getKey() + "_P", 1);

      unsigned C_INT j, jmax = mvReactions[i]->getListOfMetabReferenceGlyphs().size();

      for (j = 0; j < jmax; ++j)
        {
          CLMetabReferenceGlyph* mrg = mvReactions[i]->getListOfMetabReferenceGlyphs()[j];

          if (mrg->getRole() == CLMetabReferenceGlyph::SUBSTRATE)
            writeDotEdge(os, mrg->getMetabGlyphKey(), mvReactions[i]->getKey() + "_S");
          else if (mrg->getRole() == CLMetabReferenceGlyph::PRODUCT)
            writeDotEdge(os, mvReactions[i]->getKey() + "_P", mrg->getMetabGlyphKey());
        }
    }

  os << "}" << std::endl;
}

void CLayout::writeDotNode(std::ostream & os, const std::string & id,
                           const std::string & label,
                           int t) const
{
  std::string tmp;

  if (t == 1)
    tmp = " shape=point ";

  os << id << " [" << tmp << " label=\"" << label << "\"] \n";
}

void CLayout::writeDotEdge(std::ostream & os, const std::string & id1,
                           const std::string & id2,
                           int t) const
{
  std::string tmp;

  if (t == 1)
    tmp = " [len=0.2] ";

  os << id1 << " -> " << id2 << tmp << "\n"; //[label=\"" << label << "\"] \n";
}

void CLayout::exportToSBML(Layout * layout, const std::map<CCopasiObject*, SBase*> & copasimodelmap,
                           std::map<std::string, const SBase*>& sbmlIDs) const
{
  if (!layout) return;

  //Name and ID
  std::string id = CSBMLExporter::createUniqueId(sbmlIDs, "layout_");
  layout->setId(id);
  sbmlIDs.insert(std::pair<const std::string, const SBase*>(id, layout));
  //we do not check if the layout is already present in the libsbml data
  //structures. This is no big deal since at the moment no software
  //relies on persistent IDs for layout elements.

  //Dimensions
  Dimensions tmpDim = mDimensions.getSBMLDimensions();
  layout->setDimensions(&tmpDim);

  //some of the following code is not used at the moment:  the COPASI model map
  //does not contain glyphs. Since this may change in the future I leave the code
  //below.

  // create a map from COPASI layout object to SBML objects. We do not put
  //the layout objects into the global map (copasimodelmap) but we need to have
  //access to all objects in the current layout since speciesReferenceGlyph and
  //textGlyph need to reference other graphical objects.
  std::map<const CLBase*, const SBase*> layoutmap;

  //Compartment glyphs
  unsigned C_INT32 i, imax = mvCompartments.size();

  for (i = 0; i < imax; ++i)
    {
      CLCompartmentGlyph * tmp = mvCompartments[i];

      //check if the compartment glyph exists in the libsbml data
      std::map<CCopasiObject*, SBase*>::const_iterator it;
      it = copasimodelmap.find(tmp);

      CompartmentGlyph * pCG;

      if (it == copasimodelmap.end()) //not found
        {
          pCG = new CompartmentGlyph;
          layout->getListOfCompartmentGlyphs()->appendAndOwn(pCG);
        }
      else
        {
          pCG = dynamic_cast<CompartmentGlyph*>(it->second);
        }

      layoutmap.insert(std::pair<const CLBase*, const SBase*>(tmp, pCG));
      tmp->exportToSBML(pCG, copasimodelmap, sbmlIDs);
    }

  //Species glyphs
  imax = mvMetabs.size();

  for (i = 0; i < imax; ++i)
    {
      CLMetabGlyph * tmp = mvMetabs[i];

      //check if the glyph exists in the libsbml data
      std::map<CCopasiObject*, SBase*>::const_iterator it;
      it = copasimodelmap.find(tmp);

      SpeciesGlyph * pG;

      if (it == copasimodelmap.end()) //not found
        {
          pG = new SpeciesGlyph;
          layout->getListOfSpeciesGlyphs()->appendAndOwn(pG);
        }
      else
        {
          pG = dynamic_cast<SpeciesGlyph*>(it->second);
        }

      layoutmap.insert(std::pair<const CLBase*, const SBase*>(tmp, pG));
      tmp->exportToSBML(pG, copasimodelmap, sbmlIDs);
    }

  //Reaction glyphs
  imax = mvReactions.size();

  for (i = 0; i < imax; ++i)
    {
      CLReactionGlyph * tmp = mvReactions[i];

      //check if the glyph exists in the libsbml data
      std::map<CCopasiObject*, SBase*>::const_iterator it;
      it = copasimodelmap.find(tmp);

      ReactionGlyph * pG;

      if (it == copasimodelmap.end()) //not found
        {
          pG = new ReactionGlyph;
          layout->getListOfReactionGlyphs()->appendAndOwn(pG);
        }
      else
        {
          pG = dynamic_cast<ReactionGlyph*>(it->second);
        }

      layoutmap.insert(std::pair<const CLBase*, const SBase*>(tmp, pG));
      //we need to pass the layoutmap here for 2 reasons:
      //1. the metabreferenceglyphs need to be added
      //2. the metabreferenceglyphs need to resolve the reference to the metabglyph
      tmp->exportToSBML(pG, copasimodelmap, sbmlIDs, layoutmap);
    }

  //Text glyphs
  imax = mvLabels.size();

  for (i = 0; i < imax; ++i)
    {
      CLTextGlyph * tmp = mvLabels[i];

      //check if the glyph exists in the libsbml data
      std::map<CCopasiObject*, SBase*>::const_iterator it;
      it = copasimodelmap.find(tmp);

      TextGlyph * pG;

      if (it == copasimodelmap.end()) //not found
        {
          pG = new TextGlyph;
          layout->getListOfTextGlyphs()->appendAndOwn(pG);
        }
      else
        {
          pG = dynamic_cast<TextGlyph*>(it->second);
        }

      layoutmap.insert(std::pair<const CLBase*, const SBase*>(tmp, pG));
      tmp->exportToSBML(pG, copasimodelmap, sbmlIDs);
    }

  //generic glyphs
  imax = mvGraphicalObjects.size();

  for (i = 0; i < imax; ++i)
    {
      CLGraphicalObject * tmp = mvGraphicalObjects[i];

      //check if the glyph exists in the libsbml data
      std::map<CCopasiObject*, SBase*>::const_iterator it;
      it = copasimodelmap.find(tmp);

      GraphicalObject * pG;

      if (it == copasimodelmap.end()) //not found
        {
          pG = new GraphicalObject;
          layout->getListOfAdditionalGraphicalObjects()->appendAndOwn(pG);
        }
      else
        {
          pG = dynamic_cast<GraphicalObject*>(it->second);
        }

      layoutmap.insert(std::pair<const CLBase*, const SBase*>(tmp, pG));
      tmp->exportToSBML(pG, copasimodelmap, sbmlIDs);
    }

  //now that we have all graphical objects in the layoutmap we can resolve the references
  //in the text glyphs
  imax = mvLabels.size();

  for (i = 0; i < imax; ++i)
    {
      const CLTextGlyph * tmp = mvLabels[i];

      //find the corresponding SBML object
      std::map<const CLBase*, const SBase*>::const_iterator it = layoutmap.find(tmp);

      if (it != layoutmap.end() && it->second && dynamic_cast<const TextGlyph*>(it->second))
        {
          tmp->exportReferenceToSBML(const_cast<TextGlyph*>(dynamic_cast<const TextGlyph*>(it->second)), layoutmap);
        }
    }
}
