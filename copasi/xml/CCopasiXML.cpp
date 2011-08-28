// Begin CVS Header
//   $Source: /fs/turing/cvs/copasi_dev/copasi/xml/CCopasiXML.cpp,v $
//   $Revision: 1.130 $
//   $Name: Build-33 $
//   $Author: shoops $
//   $Date: 2010/09/21 16:48:01 $
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

/*!
 \file CCopasiXML.cpp
 \brief Implementation file of class CCopasiXML.
 */

/**
 * CCopasiXML class.
 * This class implements a CCopasiXMLInterface to the COPASI XML specified in
 * http://www.copasi.org/schema/copasi.xsd
 *
 * Created for COPASI by Stefan Hoops 2003
 * Copyright Stefan Hoops
 */
#include <iostream>
#include <map>
#include <locale>

#include "copasi.h"

#include "CCopasiXML.h"
#include "CCopasiXMLParser.h"
#include "CCopasiXMLVersion.h"
#include "CFixLocalReactionParameters.h"

#include "CopasiDataModel/CCopasiDataModel.h"
#include "report/CCopasiRootContainer.h"
#include "utilities/CCopasiVector.h"
#include "utilities/CSlider.h"
#include "model/CModel.h"
#include "model/CState.h"
#include "function/CFunctionDB.h"
#include "function/CEvaluationTree.h"
#include "report/CReportDefinitionVector.h"
#include "utilities/CCopasiTask.h"
#include "utilities/CCopasiMethod.h"
#include "utilities/CCopasiProblem.h"
#include "utilities/CDirEntry.h"
#include "utilities/utility.h"
#include "plot/COutputDefinitionVector.h"
#include "plot/CPlotItem.h"
#include "layout/CListOfLayouts.h"

#ifdef USE_CRENDER_EXTENSION

#include <copasi/layout/CLLocalRenderInformation.h>
#include <copasi/layout/CLGlobalRenderInformation.h>
#include <copasi/layout/CLRenderInformationBase.h>
#include <copasi/layout/CLLocalStyle.h>
#include <copasi/layout/CLGlobalStyle.h>
#include <copasi/layout/CLColorDefinition.h>
#include <copasi/layout/CLGradientBase.h>
#include <copasi/layout/CLLinearGradient.h>
#include <copasi/layout/CLRadialGradient.h>
#include <copasi/layout/CLLineEnding.h>
#include <copasi/layout/CLRenderPoint.h>
#include <copasi/layout/CLRenderCubicBezier.h>
#include <copasi/layout/CLGroup.h>
#include <copasi/layout/CLTransformation2D.h>
#include <copasi/layout/CLImage.h>
#include <copasi/layout/CLGraphicalPrimitive1D.h>
#include <copasi/layout/CLText.h>
#include <copasi/layout/CLRenderCurve.h>
#include <copasi/layout/CLGraphicalPrimitive2D.h>
#include <copasi/layout/CLRectangle.h>
#include <copasi/layout/CLEllipse.h>
#include <copasi/layout/CLPolygon.h>
#include <copasi/layout/CLGradientStop.h>
#include <copasi/layout/CLLineEnding.h>

#endif /* USE_CRENDER_EXTENSION */

// class CCopasiTask;
// class CCopasiReport;

CCopasiXML::CCopasiXML():
    CCopasiXMLInterface(),
    mpModel(NULL),
    mpFunctionList(NULL),
    mpTaskList(NULL),
    mpReportList(NULL),
    mpPlotList(NULL),
    mpGUI(NULL),
    mpLayoutList(NULL)
{
  mVersion.setVersion(COPASI_XML_VERSION_MAJOR,
                      COPASI_XML_VERSION_MINOR,
                      COPASI_XML_VERSION_BUILD,
                      COPASI_XML_VERSION_COMMENT);
}

CCopasiXML::~CCopasiXML() {}

bool CCopasiXML::save(std::ostream & os,
                      const std::string & relativeTo)
{
  mFilename = relativeTo;

  os.imbue(std::locale::classic());
  os.precision(16);

  mpOstream = &os;
  bool success = true;

  *mpOstream << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>"
  << std::endl;

  *mpOstream << "<!-- generated with COPASI "
  << CVersion::VERSION.getVersion()
  << " (http://www.copasi.org) at "
  << UTCTimeStamp()
  << " UTC -->"
  << std::endl;

  *mpOstream << "<?oxygen RNGSchema=\"http://www.copasi.org/static/schema/CopasiML.rng\" type=\"xml\"?>" << std::endl;

  CXMLAttributeList Attributes;
  Attributes.add("xmlns", "http://www.copasi.org/static/schema");
  Attributes.add("versionMajor", mVersion.getVersionMajor());
  Attributes.add("versionMinor", mVersion.getVersionMinor());
  Attributes.add("versionDevel", mVersion.getVersionDevel());

  /*
  Attributes.add("xmlns:xsi", "http://www.w3.org/2001/XMLSchema-instance");
  Attributes.add("xsi:noNamespaceSchemaLocation",
                 "http://www.copasi.org/static/schema/CopasiML.rng");
  */

  startSaveElement("COPASI", Attributes);

  if (haveModel() && !haveFunctionList())
    {
      if (!buildFunctionList()) success = false;

      if (!saveFunctionList()) success = false;

      if (!freeFunctionList()) success = false;
    }
  else if (!saveFunctionList()) success = false;

  if (!saveModel()) success = false;

  if (!saveTaskList()) success = false;

  if (!saveReportList()) success = false;

  if (!savePlotList()) success = false;

  if (!saveGUI()) success = false;

  if (!saveLayoutList()) success = false;

  if (!saveSBMLReference()) success = false;

  endSaveElement("COPASI");

  return success;
}

bool CCopasiXML::load(std::istream & is,
                      const std::string & relativeTo)
{
  mFilename = relativeTo;

  is.imbue(std::locale::classic());
  is.precision(16);

  mpIstream = &is;
  bool success = true;
  bool done = false;

  CVersion FileVersion;
  CCopasiXMLParser Parser(FileVersion);

  Parser.setFunctionList(mpFunctionList);
  Parser.setGUI(mpGUI);
  Parser.setLayoutList(mpLayoutList);
  Parser.setDatamodel(this->mpDataModel);

#define BUFFER_SIZE 0xfffe
  char * pBuffer = new char[BUFFER_SIZE + 1];

  try
    {
      while (!done)
        {
          mpIstream->get(pBuffer, BUFFER_SIZE, 0);

          if (mpIstream->eof()) done = true;

          if (mpIstream->fail() && !done) fatalError();

          if (!Parser.parse(pBuffer, -1, done))
            {
              CCopasiMessage Message(CCopasiMessage::RAW, MCXML + 2,
                                     Parser.getCurrentLineNumber(),
                                     Parser.getCurrentColumnNumber(),
                                     Parser.getErrorString());
              done = true;
              success = false;
            }
        }
    }

  catch (...)
    {
      success = false;
    }

  delete [] pBuffer;
#undef BUFFER_SIZE

  mpModel = Parser.getModel();
  mpReportList = Parser.getReportList();
  mpTaskList = Parser.getTaskList();
  mpPlotList = Parser.getPlotList();
  mpLayoutList = Parser.getLayoutList();

  if (!success)
    {
      pdelete(mpModel);
      pdelete(mpReportList);
      pdelete(mpTaskList);
      pdelete(mpPlotList);
      pdelete(mpLayoutList);
    }

  if (FileVersion.getVersionDevel() > mVersion.getVersionDevel())
    CCopasiMessage(CCopasiMessage::WARNING, MCXML + 9,
                   mFilename.c_str(), FileVersion.getVersion().c_str());

  return success;
}

const CVersion & CCopasiXML::getVersion() const
{return mVersion;}

bool CCopasiXML::setModel(CModel * pModel)
{
  mpModel = pModel;
  return true;
}

CModel * CCopasiXML::getModel() const {return mpModel;}

bool CCopasiXML::haveModel() const {return mpModel != NULL;}

bool CCopasiXML::freeModel()
{
  pdelete(mpModel);
  return true;
}

bool CCopasiXML::setFunctionList(CCopasiVectorN< CEvaluationTree > *pFunctionList)
{
  mpFunctionList = pFunctionList;
  return true;
}

CCopasiVectorN< CEvaluationTree > * CCopasiXML::getFunctionList() const
{return mpFunctionList;}

bool CCopasiXML::haveFunctionList() const
{return mpFunctionList != NULL;}

bool CCopasiXML::freeFunctionList()
{
  pdelete(mpFunctionList);
  return true;
}

bool CCopasiXML::setTaskList(CCopasiVectorN< CCopasiTask > * pTaskList)
{
  mpTaskList = pTaskList;
  return true;
}

/**
 * Set the datamodel.
 * @param CCopasiDataModel* pDataModel
 * @return bool success
 */
bool CCopasiXML::setDatamodel(CCopasiDataModel* pDataModel)
{
  this->mpDataModel = pDataModel;
  return true;
}

CCopasiVectorN< CCopasiTask > * CCopasiXML::getTaskList() const
{return mpTaskList;}

bool CCopasiXML::haveTaskList() const
{return mpTaskList != NULL;}

bool CCopasiXML::freeTaskList()
{
  pdelete(mpTaskList);
  return true;
}

//************

bool CCopasiXML::setPlotList(COutputDefinitionVector * pPlotList)
{
  mpPlotList = pPlotList;
  return true;
}

COutputDefinitionVector * CCopasiXML::getPlotList() const
{return mpPlotList;}

bool CCopasiXML::havePlotList() const
{return mpPlotList != NULL;}

bool CCopasiXML::freePlotList()
{
  pdelete(mpPlotList);
  return true;
}

//************

bool CCopasiXML::setReportList(CReportDefinitionVector * pReportList)
{
  mpReportList = pReportList;
  return true;
}

CReportDefinitionVector * CCopasiXML::getReportList() const
{return mpReportList;}

bool CCopasiXML::haveReportList() const
{return mpReportList != NULL;}

bool CCopasiXML::freeReportList()
{
  pdelete(mpReportList);
  return true;
}

//************

bool CCopasiXML::setGUI(SCopasiXMLGUI * pGUI)
{
  mpGUI = pGUI;
  return true;
}

SCopasiXMLGUI * CCopasiXML::getGUI() const {return mpGUI;}

bool CCopasiXML::haveGUI() const {return mpGUI != NULL;}

bool CCopasiXML::freeGUI()
{
  pdelete(mpGUI);
  return true;
}

//************

bool CCopasiXML::setLayoutList(const CListOfLayouts & layoutList)
{
  mpLayoutList = const_cast<CListOfLayouts *>(&layoutList);
  return true;
}

CListOfLayouts * CCopasiXML::getLayoutList() const
{return mpLayoutList;}

bool CCopasiXML::haveLayoutList() const
{return mpLayoutList != NULL;}

bool CCopasiXML::freeLayoutList()
{
  pdelete(mpLayoutList);
  return true;
}

bool CCopasiXML::saveModel()
{
  bool success = true;

  if (!haveModel()) return success;

  CXMLAttributeList Attributes;
  Attributes.add("key", mpModel->getKey());
  Attributes.add("name", mpModel->getObjectName());
  Attributes.add("simulationType", CModelEntity::XMLStatus[mpModel->getStatus()]);
  Attributes.add("timeUnit", mpModel->getTimeUnitName());
  Attributes.add("volumeUnit", mpModel->getVolumeUnitName());
  Attributes.add("areaUnit", mpModel->getAreaUnitName());
  Attributes.add("lengthUnit", mpModel->getLengthUnitName());
  Attributes.add("quantityUnit", mpModel->getQuantityUnitName());
  Attributes.add("type", CModel::ModelTypeNames[mpModel->getModelType()]);
  Attributes.add("avogadroConstant", mpModel->getAvogadro());

  // This is now optional
  // Attributes.add("stateVariable",
  //                mpModel->getStateTemplate().getKey(mpModel->getKey()));

  startSaveElement("Model", Attributes);

  if (mpModel->getMiriamAnnotation() != "")
    {
      startSaveElement("MiriamAnnotation");
      *mpOstream << mpModel->getMiriamAnnotation() << std::endl;
      endSaveElement("MiriamAnnotation");
    }

  if (mpModel->getNotes() != "")
    {
      startSaveElement("Comment");
      saveXhtml(mpModel->getNotes());
      endSaveElement("Comment");
    }

  if (mpModel->getInitialExpression() != "")
    {
      startSaveElement("InitialExpression");
      saveData(mpModel->getInitialExpression());
      endSaveElement("InitialExpression");
    }

  unsigned C_INT32 i, imax;

  // Compartment
  if ((imax = mpModel->getCompartments().size()) > 0)
    {
      startSaveElement("ListOfCompartments");

      Attributes.erase();
      Attributes.add("key", "");
      Attributes.add("name", "");
      Attributes.add("simulationType", "");
      Attributes.add("dimensionality", "");

      unsigned C_INT32 i, imax = mpModel->getCompartments().size();

      for (i = 0; i < imax; i++)
        {
          CCompartment * pComp = mpModel->getCompartments()[i];

          Attributes.setValue(0, pComp->getKey());
          Attributes.setValue(1, pComp->getObjectName());
          CModelEntity::Status SimulationType = pComp->getStatus();
          Attributes.setValue(2, CModelEntity::XMLStatus[SimulationType]);
          Attributes.setValue(3, pComp->getDimensionality());

          startSaveElement("Compartment", Attributes);

          if (pComp->getMiriamAnnotation() != "")
            {
              startSaveElement("MiriamAnnotation");
              *mpOstream << pComp->getMiriamAnnotation() << std::endl;
              endSaveElement("MiriamAnnotation");
            }

          if (pComp->getNotes() != "")
            {
              startSaveElement("Comment");
              saveXhtml(pComp->getNotes());
              endSaveElement("Comment");
            }

          if (SimulationType != CModelEntity::FIXED &&
              pComp->getExpression() != "")
            {
              startSaveElement("Expression");
              saveData(pComp->getExpression());
              endSaveElement("Expression");
            }

          if (pComp->getInitialExpression() != "")
            {
              startSaveElement("InitialExpression");
              saveData(pComp->getInitialExpression());
              endSaveElement("InitialExpression");
            }

          endSaveElement("Compartment");

          if (pComp->getSBMLId() != "")
            mSBMLReference[pComp->getSBMLId()] = pComp->getKey();
        }

      endSaveElement("ListOfCompartments");
    }

  // Metabolites (aka. Species)
  if ((imax = mpModel->getMetabolites().size()) > 0)
    {
      startSaveElement("ListOfMetabolites");

      Attributes.erase();
      Attributes.add("key", "");
      Attributes.add("name", "");
      Attributes.add("simulationType", "");
      Attributes.add("compartment", "");

      for (i = 0; i < imax; i++)
        {
          CMetab * pMetab = mpModel->getMetabolites()[i];

          Attributes.setValue(0, pMetab->getKey());
          Attributes.setValue(1, pMetab->getObjectName());
          CModelEntity::Status SimulationType = pMetab->getStatus();
          Attributes.setValue(2, CModelEntity::XMLStatus[SimulationType]);
          Attributes.setValue(3, pMetab->getCompartment()->getKey());

          startSaveElement("Metabolite", Attributes);

          if (pMetab->getMiriamAnnotation() != "")
            {
              startSaveElement("MiriamAnnotation");
              *mpOstream << pMetab->getMiriamAnnotation() << std::endl;
              endSaveElement("MiriamAnnotation");
            }

          if (pMetab->getNotes() != "")
            {
              startSaveElement("Comment");
              saveXhtml(pMetab->getNotes());
              endSaveElement("Comment");
            }

          if (SimulationType != CModelEntity::FIXED &&
              SimulationType != CModelEntity::REACTIONS &&
              pMetab->getExpression() != "")
            {
              startSaveElement("Expression");
              saveData(pMetab->getExpression());
              endSaveElement("Expression");
            }

          if (pMetab->getInitialExpression() != "")
            {
              startSaveElement("InitialExpression");
              saveData(pMetab->getInitialExpression());
              endSaveElement("InitialExpression");
            }

          endSaveElement("Metabolite");

          if (pMetab->getSBMLId() != "")
            mSBMLReference[pMetab->getSBMLId()] = pMetab->getKey();
        }

      endSaveElement("ListOfMetabolites");
    }

  // Model Values (aka. Global Quantities)
  if ((imax = mpModel->getModelValues().size()) > 0)
    {
      startSaveElement("ListOfModelValues");

      Attributes.erase();
      Attributes.add("key", "");
      Attributes.add("name", "");
      Attributes.add("simulationType", "");

      for (i = 0; i < imax; i++)
        {
          CModelValue * pMV = mpModel->getModelValues()[i];

          Attributes.setValue(0, pMV->getKey());
          Attributes.setValue(1, pMV->getObjectName());
          CModelEntity::Status SimulationType = pMV->getStatus();
          Attributes.setValue(2, CModelEntity::XMLStatus[SimulationType]);

          startSaveElement("ModelValue", Attributes);

          if (pMV->getMiriamAnnotation() != "")
            {
              startSaveElement("MiriamAnnotation");
              *mpOstream << pMV->getMiriamAnnotation() << std::endl;
              endSaveElement("MiriamAnnotation");
            }

          if (pMV->getNotes() != "")
            {
              startSaveElement("Comment");
              saveXhtml(pMV->getNotes());
              endSaveElement("Comment");
            }

          if (SimulationType != CModelEntity::FIXED &&
              pMV->getExpression() != "")
            {
              startSaveElement("Expression");
              saveData(pMV->getExpression());
              endSaveElement("Expression");
            }

          if (pMV->getInitialExpression() != "")
            {
              startSaveElement("InitialExpression");
              saveData(pMV->getInitialExpression());
              endSaveElement("InitialExpression");
            }

          endSaveElement("ModelValue");

          if (pMV->getSBMLId() != "")
            mSBMLReference[pMV->getSBMLId()] = pMV->getKey();
        }

      endSaveElement("ListOfModelValues");
    }

  // Reactions
  if ((imax = mpModel->getReactions().size()) > 0)
    {
      startSaveElement("ListOfReactions");

      CXMLAttributeList Attr;
      const CCopasiVector< CChemEqElement > * pReactantList;
      unsigned C_INT32 j, jmax;

      std::vector< const CCopasiObject * > ObjectList;
      unsigned C_INT32 k, kmax;

      Attributes.erase();
      Attributes.add("key", "");
      Attributes.add("name", "");
      Attributes.add("reversible", "");

      for (i = 0; i < imax; i++)
        {
          CReaction * pReaction = mpModel->getReactions()[i];

          Attributes.setValue(0, pReaction->getKey());
          Attributes.setValue(1, pReaction->getObjectName());
          Attributes.setValue(2, pReaction->isReversible() ? "true" : "false");

          if (pReaction->getSBMLId() != "")
            mSBMLReference[pReaction->getSBMLId()] = pReaction->getKey();

          startSaveElement("Reaction", Attributes);

          if (pReaction->getMiriamAnnotation() != "")
            {
              startSaveElement("MiriamAnnotation");
              *mpOstream << pReaction->getMiriamAnnotation() << std::endl;
              endSaveElement("MiriamAnnotation");
            }

          if (pReaction->getNotes() != "")
            {
              startSaveElement("Comment");
              saveXhtml(pReaction->getNotes());
              endSaveElement("Comment");
            }

          Attr.erase();
          Attr.add("metabolite", "");
          Attr.add("stoichiometry", "");

          pReactantList = & pReaction->getChemEq().getSubstrates();

          if ((jmax = pReactantList->size()) > 0)
            {
              startSaveElement("ListOfSubstrates");

              for (j = 0; j < jmax; j++)
                {
                  Attr.setValue(0, (*pReactantList)[j]->getMetaboliteKey());
                  Attr.setValue(1, (*pReactantList)[j]->getMultiplicity());

                  saveElement("Substrate", Attr);
                }

              endSaveElement("ListOfSubstrates");
            }

          //startSaveElement("ListOfProducts"); // this seems to belong further down

          pReactantList = & pReaction->getChemEq().getProducts();

          if ((jmax = pReactantList->size()) > 0)
            {
              startSaveElement("ListOfProducts"); //this seems to belong here

              for (j = 0; j < jmax; j++)
                {
                  Attr.setValue(0, (*pReactantList)[j]->getMetaboliteKey());
                  Attr.setValue(1, (*pReactantList)[j]->getMultiplicity());

                  saveElement("Product", Attr);
                }

              endSaveElement("ListOfProducts");
            }

          pReactantList = & pReaction->getChemEq().getModifiers();

          if ((jmax = pReactantList->size()) > 0)
            {
              startSaveElement("ListOfModifiers");

              for (j = 0, jmax = pReactantList->size(); j < jmax; j++)
                {
                  Attr.setValue(0, (*pReactantList)[j]->getMetaboliteKey());
                  Attr.setValue(1, (*pReactantList)[j]->getMultiplicity());

                  saveElement("Modifier", Attr);
                }

              endSaveElement("ListOfModifiers");
            }

          const CCopasiParameterGroup * pParamList;

          pParamList = & pReaction->getParameters();

          if ((jmax = pParamList->size()) > 0)
            {
              startSaveElement("ListOfConstants");

              Attr.erase();
              Attr.add("key", "");
              Attr.add("name", "");
              Attr.add("value", "");

              for (j = 0; j < jmax; j++)
                {
                  Attr.setValue(0, pParamList->getKey(j));
                  Attr.setValue(1, pParamList->getName(j));
                  Attr.setValue(2, * pParamList->getValue(j).pDOUBLE);

                  saveElement("Constant", Attr);
                }

              endSaveElement("ListOfConstants");
            }

          if (pReaction->getFunction() !=
              dynamic_cast<CFunction *>(CCopasiRootContainer::getKeyFactory()->get("UndefinedFunction_0")))
            {
              Attr.erase();
              Attr.add("function", pReaction->getFunction()->getKey());
              startSaveElement("KineticLaw", Attr);

              if ((jmax = pReaction->getFunctionParameters().size()))
                {
                  startSaveElement("ListOfCallParameters");
                  const std::vector< std::vector<std::string> > & rMap =
                    pReaction->getParameterMappings();

                  for (j = 0; j < jmax; j++)
                    {
                      Attr.erase();
                      Attr.add("functionParameter",
                               pReaction->
                               getFunction()->getVariables()[j]->getKey());

                      startSaveElement("CallParameter", Attr);

                      Attr.erase();
                      Attr.add("reference", "");

                      for (k = 0, kmax = rMap[j].size(); k < kmax; k++)
                        {
                          Attr.setValue(0, rMap[j][k]);
                          saveElement("SourceParameter", Attr);
                        }

                      endSaveElement("CallParameter");
                    }

                  endSaveElement("ListOfCallParameters");
                }

              endSaveElement("KineticLaw");
            }

          endSaveElement("Reaction");
        }

      endSaveElement("ListOfReactions");
    }

  // Events (added 07.04.08)
  if ((imax = mpModel->getEvents().size()) > 0)
    {
      startSaveElement("ListOfEvents");

      Attributes.erase();
      Attributes.add("key", "");
      Attributes.add("name", "");
      Attributes.add("order", "");
      Attributes.add("delayAssignment", "");

      for (i = 0; i < imax; i++)
        {
          CEvent * pEvent = mpModel->getEvents()[i];

          Attributes.setValue(0, pEvent->getKey());
          Attributes.setValue(1, pEvent->getObjectName());
          Attributes.setValue(2, pEvent->getOrder());

          if (pEvent->getDelayExpression() != "")
            {
              Attributes.setValue(3, pEvent->getDelayAssignment() ? "true" : "false");
            }
          else
            {
              Attributes.skip(3);
            }

          startSaveElement("Event", Attributes);

          if (pEvent->getMiriamAnnotation() != "")
            {
              startSaveElement("MiriamAnnotation");
              *mpOstream << pEvent->getMiriamAnnotation() << std::endl;
              endSaveElement("MiriamAnnotation");
            }

          if (pEvent->getNotes() != "")
            {
              startSaveElement("Comment");
              saveXhtml(pEvent->getNotes());
              endSaveElement("Comment");
            }

          if (pEvent->getTriggerExpression() != "")
            {
              startSaveElement("TriggerExpression");
              saveData(pEvent->getTriggerExpression());
              endSaveElement("TriggerExpression");
            }

          if (pEvent->getDelayExpression() != "")
            {
              startSaveElement("DelayExpression");
              saveData(pEvent->getDelayExpression());
              endSaveElement("DelayExpression");
            }

          const CCopasiVectorN< CEventAssignment > & Assignments = pEvent->getAssignments();

          if (Assignments.size() > 0)
            {
              startSaveElement("ListOfAssignments");

              CXMLAttributeList Attr;
              Attr.add("targetKey", "");

              CCopasiVectorN< CEventAssignment >::const_iterator it = Assignments.begin();
              CCopasiVectorN< CEventAssignment >::const_iterator end = Assignments.end();

              for (; it != end; ++it)
                {
                  Attr.setValue(0, (*it)->getTargetKey());

                  startSaveElement("Assignment", Attr);

                  startSaveElement("Expression");
                  saveData((*it)->getExpression());
                  endSaveElement("Expression");

                  endSaveElement("Assignment");
                }

              endSaveElement("ListOfAssignments");
            }

          endSaveElement("Event");

          if (pEvent->getSBMLId() != "")
            mSBMLReference[pEvent->getSBMLId()] = pEvent->getKey();
        }

      endSaveElement("ListOfEvents");
    }

  startSaveElement("StateTemplate");

  Attributes.erase();
  // This is now optional.
  // Attributes.add("key", "");
  Attributes.add("objectReference", "");
  std::pair< std::string, std::string > Variable;

  CModelEntity *const* ppEntity = mpModel->getStateTemplate().getEntities();
  CModelEntity *const* ppEntityEnd = ppEntity + mpModel->getStateTemplate().size();

  for (; ppEntity != ppEntityEnd; ++ppEntity)
    {
      Attributes.setValue(0, (*ppEntity)->getKey());

      saveElement("StateTemplateVariable", Attributes);
    }

  endSaveElement("StateTemplate");

  Attributes.erase();
  Attributes.add("type", "initialState");
  startSaveElement("InitialState", Attributes);
  *mpOstream << mIndent;
  ppEntity = mpModel->getStateTemplate().getEntities();

  for (; ppEntity != ppEntityEnd; ++ppEntity)
    {
      *mpOstream << (DBL)(*ppEntity)->getInitialValue() << " ";
    }

  *mpOstream << std::endl;

  endSaveElement("InitialState");

  endSaveElement("Model");

  return success;
}

bool CCopasiXML::saveFunctionList()
{
  bool success = true;

  if (!haveFunctionList()) return success;

  unsigned C_INT32 i, imax = mpFunctionList->size();

  if (!imax) return success;

  CXMLAttributeList Attributes;
  CEvaluationTree * pEvaluationTree = NULL;
  CFunction * pFunction = NULL;

  startSaveElement("ListOfFunctions");

  for (i = 0; i < imax; i++)
    {
      pEvaluationTree = (*mpFunctionList)[i];
      pFunction = dynamic_cast<CFunction *>(pEvaluationTree);

      Attributes.erase();
      Attributes.add("key", pEvaluationTree->getKey());
      Attributes.add("name", pEvaluationTree->getObjectName());
      Attributes.add("type", CEvaluationTree::XMLType[pEvaluationTree->getType()]);

      if (pFunction)
        {
          switch (pFunction->isReversible())
            {
              case TriUnspecified:
                Attributes.add("reversible", "unspecified");
                break;

              case TriFalse:
                Attributes.add("reversible", "false");
                break;

              case TriTrue:
                Attributes.add("reversible", "true");
                break;
            }

          if (pFunction->getSBMLId() != "")
            mSBMLReference[pFunction->getSBMLId()] = pFunction->getKey();
        }

      startSaveElement("Function", Attributes);

      if (pEvaluationTree->getMiriamAnnotation() != "")
        {
          startSaveElement("MiriamAnnotation");
          *mpOstream << pEvaluationTree->getMiriamAnnotation() << std::endl;
          endSaveElement("MiriamAnnotation");
        }

      if (pEvaluationTree->getNotes() != "")
        {
          startSaveElement("Comment");
          saveXhtml(pEvaluationTree->getNotes());
          endSaveElement("Comment");
        }

      startSaveElement("Expression");
      saveData(pEvaluationTree->getInfix());
      endSaveElement("Expression");

      if (pFunction)
        {
          startSaveElement("ListOfParameterDescriptions");

          unsigned C_INT32 j, jmax = pFunction->getVariables().size();
          CFunctionParameter * pParameter;

          Attributes.erase();
          Attributes.add("key", "");
          Attributes.add("name", "");
          Attributes.add("order", "");
          Attributes.add("role", "");

          for (j = 0; j < jmax; j++)
            {
              pParameter = pFunction->getVariables()[j];
              Attributes.setValue(0, pParameter->getKey());
              Attributes.setValue(1, pParameter->getObjectName());
              Attributes.setValue(2, j);
              Attributes.setValue(3, CFunctionParameter::RoleNameXML[pParameter->getUsage()]);

              saveElement("ParameterDescription", Attributes);
            }

          endSaveElement("ListOfParameterDescriptions");
        }

      endSaveElement("Function");
    }

  endSaveElement("ListOfFunctions");

  return success;
}

bool CCopasiXML::savePlotList()
{
  //std::cerr << "Saving plot list. " << std::endl;
  bool success = true;

  if (!havePlotList())
    {
      //std::cerr << "No plot list defined." << std::endl;
      return success;
    }

  unsigned C_INT32 i, imax = mpPlotList->size();

  //std::cerr << "Saving " << imax << " plots." << std::endl;
  if (!imax) return success;

  CXMLAttributeList Attributes;

  startSaveElement("ListOfPlots");

  for (i = 0; i < imax; i++)
    {
      const CPlotSpecification* pPlot = (*mpPlotList)[i];

      Attributes.erase();
      Attributes.add("name", pPlot->getObjectName());
      Attributes.add("type", CPlotSpecification::XMLType[pPlot->getType()]);
      Attributes.add("active", pPlot->isActive());
      startSaveElement("PlotSpecification", Attributes);
      saveParameterGroup(* pPlot->CCopasiParameter::getValue().pGROUP);
      startSaveElement("ListOfPlotItems");
      unsigned C_INT32 j, jmax = pPlot->getItems().size();

      //std::cerr << "Saving " << jmax << "PlotItems." << std::endl;
      for (j = 0; j < jmax; j++)
        {
          const CPlotItem* pPlotItem = pPlot->getItems()[j];
          Attributes.erase();
          Attributes.add("name", pPlotItem->getObjectName());
          Attributes.add("type", CPlotItem::XMLType[pPlotItem->getType()]);
          startSaveElement("PlotItem", Attributes);
          saveParameterGroup(* pPlotItem->CCopasiParameter::getValue().pGROUP);
          startSaveElement("ListOfChannels");
          unsigned C_INT32 k, kmax = pPlotItem->getNumChannels();

          //std::cerr << "Saving " << kmax << " Channels." << std::endl;
          for (k = 0; k < kmax; k++)
            {
              const CPlotDataChannelSpec pDataChannelSpec = pPlotItem->getChannels()[k];
              Attributes.erase();
              Attributes.add("cn", pDataChannelSpec);

              if (!pDataChannelSpec.minAutoscale)
                {
                  Attributes.add("min", pDataChannelSpec.min);
                }

              if (!pDataChannelSpec.maxAutoscale)
                {
                  Attributes.add("max", pDataChannelSpec.max);
                }

              saveElement("ChannelSpec", Attributes);
            }

          endSaveElement("ListOfChannels");
          endSaveElement("PlotItem");
        }

      endSaveElement("ListOfPlotItems");
      endSaveElement("PlotSpecification");
    }

  endSaveElement("ListOfPlots");
  return success;
}

//Mrinmayee
bool CCopasiXML::saveTaskList()
{
  bool success = true;

  if (!haveTaskList()) return success;

  unsigned C_INT32 i, imax = mpTaskList->size();

  if (!imax) return success;

  CXMLAttributeList Attributes;
  CCopasiTask * pTask = NULL;

  startSaveElement("ListOfTasks");

  for (i = 0; i < imax; i++)
    {
      pTask = (*mpTaskList)[i];

      Attributes.erase();
      Attributes.add("key", pTask->getKey());
      Attributes.add("name", pTask->getObjectName());
      Attributes.add("type", CCopasiTask::XMLType[pTask->getType()]);
      Attributes.add("scheduled", pTask->isScheduled() ? "true" : "false");
      Attributes.add("updateModel", pTask->isUpdateModel() ? "true" : "false");

      startSaveElement("Task", Attributes);

      // Report Element
      CReport & tReport = pTask->getReport();

      if (tReport.getReportDefinition())
        {
          Attributes.erase();
          Attributes.add("reference", tReport.getReportDefinition()->getKey());

          std::string Target = tReport.getTarget();

          if (!CDirEntry::isRelativePath(Target) &&
              !CDirEntry::makePathRelative(Target, mFilename))
            Target = CDirEntry::fileName(Target);

          Attributes.add("target", Target);
          Attributes.add("append", tReport.append());
          saveElement("Report", Attributes);
        }

      //Problem Element
      CCopasiProblem *tProblem = pTask->getProblem();

      Attributes.erase();
      startSaveElement("Problem");
      saveParameterGroup(* tProblem->CCopasiParameter::getValue().pGROUP);
      endSaveElement("Problem");

      // Method Element
      CCopasiMethod *tMethod = pTask->getMethod();

      Attributes.erase();
      Attributes.add("name", tMethod->CCopasiParameter::getObjectName());
      Attributes.add("type", CCopasiMethod::XMLSubType[tMethod->getSubType()]);
      startSaveElement("Method", Attributes);
      saveParameterGroup(* tMethod->CCopasiParameter::getValue().pGROUP);
      endSaveElement("Method");

      endSaveElement("Task");
    }

  endSaveElement("ListOfTasks");

  return success;
}

//Mrinmayee
bool CCopasiXML::saveReportSection(const std::string & name,
                                   const std::vector <CRegisteredObjectName> & section)
{
  CXMLAttributeList Attributes;
  Attributes.add("NoName", "");

  unsigned C_INT32 j, jmax = section.size();

  if (jmax)
    {
      startSaveElement(name);

      for (j = 0; j < jmax; j++)
        {
          if (section[j].getObjectType() == "html")
            {
              //Write in Text
              startSaveElement("html");
              Attributes.set(0, "xmlns", "http://www.w3.org/1999/xhtml");
              startSaveElement("body", Attributes);
              saveData(section[j].getObjectName()); //TODO check
              endSaveElement("body");
              endSaveElement("html");
            }
          else
            {
              //Write in Object
              Attributes.set(0, "cn", section[j]);
              saveElement("Object", Attributes);
            }
        }

      endSaveElement(name);
    }

  return true;
}

bool CCopasiXML::saveReportList()
{
  bool success = true;

  if (!haveReportList()) return success;

  unsigned C_INT32 i, imax = mpReportList->size();

  if (!imax) return success;

  CXMLAttributeList Attributes;
  CReportDefinition * pReport = NULL;

  startSaveElement("ListOfReports");

  for (i = 0; i < imax; i++)
    {
      pReport = (*mpReportList)[i];

      Attributes.erase();
      Attributes.add("key", pReport->getKey());
      Attributes.add("name", pReport->getObjectName());
      Attributes.add("taskType", CCopasiTask::XMLType[pReport->getTaskType()]);
      Attributes.add("separator", pReport->getSeparator().getStaticString());
      Attributes.add("precision", pReport->getPrecision());

      startSaveElement("Report", Attributes);

      startSaveElement("Comment");
      saveXhtml(pReport->getComment());
      endSaveElement("Comment");

      if (pReport->isTable())
        {
          Attributes.erase();
          Attributes.add("printTitle", pReport->getTitle());
          startSaveElement("Table", Attributes);

          const std::vector <CRegisteredObjectName> & Table = * pReport->getTableAddr();
          unsigned C_INT32 j, jmax = Table.size();

          Attributes.erase();
          Attributes.add("cn", "");

          for (j = 0; j < jmax; j++)
            {
              Attributes.setValue(0, Table[j]);
              saveElement("Object", Attributes);
            }

          endSaveElement("Table");
        }
      else
        {
          saveReportSection("Header", * pReport->getHeaderAddr());
          saveReportSection("Body", * pReport->getBodyAddr());
          saveReportSection("Footer", * pReport->getFooterAddr());
        }

      endSaveElement("Report");
    }

  endSaveElement("ListOfReports");

  return success;
}

void CCopasiXML::savePosition(const CLPoint& p, const std::string & tag)
{
  CXMLAttributeList Attributes;
  Attributes.erase();
  Attributes.add("x", p.getX());
  Attributes.add("y", p.getY());
  saveElement(tag, Attributes);
}

void CCopasiXML::saveDimensions(const CLDimensions& d)
{
  CXMLAttributeList Attributes;
  Attributes.erase();
  Attributes.add("width", d.getWidth());
  Attributes.add("height", d.getHeight());
  saveElement("Dimensions", Attributes);
}

void CCopasiXML::saveBoundingBox(const CLBoundingBox& bb)
{
  startSaveElement("BoundingBox");
  savePosition(bb.getPosition());
  saveDimensions(bb.getDimensions());
  endSaveElement("BoundingBox");
}

void CCopasiXML::saveCurve(const CLCurve& c)
{
  CXMLAttributeList Attributes;
  startSaveElement("Curve");

  if (c.getNumCurveSegments() > 0)
    {
      startSaveElement("ListOfCurveSegments");
      unsigned C_INT32 i, imax = c.getNumCurveSegments();

      for (i = 0; i < imax; ++i)
        {
          const CLLineSegment & cs = c.getCurveSegments()[i];

          Attributes.erase();

          if (cs.isBezier())
            Attributes.add("xsi:type", "CubicBezier");
          else
            Attributes.add("xsi:type", "LineSegment");

          startSaveElement("CurveSegment", Attributes);

          savePosition(cs.getStart(), "Start");
          savePosition(cs.getEnd(), "End");

          if (cs.isBezier())
            {
              savePosition(cs.getBase1(), "BasePoint1");
              savePosition(cs.getBase2(), "BasePoint2");
            }

          endSaveElement("CurveSegment");
        }

      endSaveElement("ListOfCurveSegments");
    }

  endSaveElement("Curve");
}

bool CCopasiXML::saveLayoutList()
{
  bool success = true;

  if (!haveLayoutList()) return success;

  unsigned C_INT32 i, imax = mpLayoutList->size();

  if (!imax) return success;

  CXMLAttributeList Attributes;
  CLayout * pLayout = NULL;
  Attributes.add("xmlns:xsi", "http://www.w3.org/2001/XMLSchema-instance");
  startSaveElement("ListOfLayouts", Attributes);

  for (i = 0; i < imax; i++)
    {
      pLayout = (*mpLayoutList)[i];

      Attributes.erase();
      Attributes.add("key", pLayout->getKey());
      Attributes.add("name", pLayout->getObjectName());
      startSaveElement("Layout", Attributes);

      Attributes.erase();
      Attributes.add("width", pLayout->getDimensions().getWidth());
      Attributes.add("height", pLayout->getDimensions().getHeight());
      saveElement("Dimensions", Attributes);

      unsigned C_INT32 j, jmax;

      //compartment glyphs
      if (pLayout->getListOfCompartmentGlyphs().size() > 0)
        {
          startSaveElement("ListOfCompartmentGlyphs");

          jmax = pLayout->getListOfCompartmentGlyphs().size();

          for (j = 0; j < jmax; ++j)
            {
              CLCompartmentGlyph* cg = pLayout->getListOfCompartmentGlyphs()[j];
              Attributes.erase();
              Attributes.add("key", cg->getKey());
              Attributes.add("name", cg->getObjectName());
              Attributes.add("compartment", cg->getModelObjectKey());
#ifdef USE_CRENDER_EXTENSION

              if (cg->getObjectRole().find_first_not_of(" \t\r\n") != std::string::npos)
                {
                  Attributes.add("objectRole", cg->getObjectRole());
                }

#endif // USE_CRENDER_EXTENSION
              startSaveElement("CompartmentGlyph", Attributes);

              saveBoundingBox(cg->getBoundingBox());

              endSaveElement("CompartmentGlyph");
            }

          endSaveElement("ListOfCompartmentGlyphs");
        }

      //species glyphs
      if (pLayout->getListOfMetaboliteGlyphs().size() > 0)
        {
          startSaveElement("ListOfMetabGlyphs");

          jmax = pLayout->getListOfMetaboliteGlyphs().size();

          for (j = 0; j < jmax; ++j)
            {
              CLMetabGlyph* cg = pLayout->getListOfMetaboliteGlyphs()[j];
              Attributes.erase();
              Attributes.add("key", cg->getKey());
              Attributes.add("name", cg->getObjectName());
              Attributes.add("metabolite", cg->getModelObjectKey());
#ifdef USE_CRENDER_EXTENSION

              if (cg->getObjectRole().find_first_not_of(" \t\r\n") != std::string::npos)
                {
                  Attributes.add("objectRole", cg->getObjectRole());
                }

#endif // USE_CRENDER_EXTENSION
              startSaveElement("MetaboliteGlyph", Attributes);

              saveBoundingBox(cg->getBoundingBox());

              endSaveElement("MetaboliteGlyph");
            }

          endSaveElement("ListOfMetabGlyphs");
        }

      //reaction glyphs
      if (pLayout->getListOfReactionGlyphs().size() > 0)
        {
          startSaveElement("ListOfReactionGlyphs");

          jmax = pLayout->getListOfReactionGlyphs().size();

          for (j = 0; j < jmax; ++j)
            {
              CLReactionGlyph* cg = pLayout->getListOfReactionGlyphs()[j];
              Attributes.erase();
              Attributes.add("key", cg->getKey());
              Attributes.add("name", cg->getObjectName());
              Attributes.add("reaction", cg->getModelObjectKey());
#ifdef USE_CRENDER_EXTENSION

              if (cg->getObjectRole().find_first_not_of(" \t\r\n") != std::string::npos)
                {
                  Attributes.add("objectRole", cg->getObjectRole());
                }

#endif // USE_CRENDER_EXTENSION
              startSaveElement("ReactionGlyph", Attributes);

              if (cg->getCurve().getNumCurveSegments() == 0)
                saveBoundingBox(cg->getBoundingBox());
              else
                saveCurve(cg->getCurve());

              // metab reference glyphs
              startSaveElement("ListOfMetaboliteReferenceGlyphs");
              unsigned C_INT32 k, kmax = cg->getListOfMetabReferenceGlyphs().size();

              for (k = 0; k < kmax; ++k)
                {
                  CLMetabReferenceGlyph * mrg = cg->getListOfMetabReferenceGlyphs()[k];
                  Attributes.erase();
                  Attributes.add("key", mrg->getKey());
                  Attributes.add("name", mrg->getObjectName());
                  Attributes.add("metaboliteGlyph", mrg->getMetabGlyphKey());
                  //Attributes.add("metaboliteReference", mrg->getXXX());
                  Attributes.add("role", CLMetabReferenceGlyph::XMLRole[mrg->getRole()]);
#ifdef USE_CRENDER_EXTENSION

                  if (mrg->getObjectRole().find_first_not_of(" \t\r\n") != std::string::npos)
                    {
                      Attributes.add("objectRole", mrg->getObjectRole());
                    }

#endif // USE_CRENDER_EXTENSION
                  startSaveElement("MetaboliteReferenceGlyph", Attributes);

                  if (mrg->getCurve().getNumCurveSegments() == 0)
                    saveBoundingBox(mrg->getBoundingBox());
                  else
                    saveCurve(mrg->getCurve());

                  endSaveElement("MetaboliteReferenceGlyph");
                }

              endSaveElement("ListOfMetaboliteReferenceGlyphs");

              endSaveElement("ReactionGlyph");
            }

          endSaveElement("ListOfReactionGlyphs");
        }

      //text Glyphs
      if (pLayout->getListOfTextGlyphs().size() > 0)
        {
          startSaveElement("ListOfTextGlyphs");

          jmax = pLayout->getListOfTextGlyphs().size();

          for (j = 0; j < jmax; ++j)
            {
              CLTextGlyph* cg = pLayout->getListOfTextGlyphs()[j];
              // we only export the text glyph if it either has a text
              // or a valid originOfText
              std::string id = cg->getModelObjectKey();

              if (cg->isTextSet() || id.find_first_not_of(" \t\r\n") != std::string::npos)
                {
                  Attributes.erase();
                  Attributes.add("key", cg->getKey());
                  Attributes.add("name", cg->getObjectName());
                  Attributes.add("graphicalObject", cg->getGraphicalObjectKey());

                  if (cg->isTextSet())
                    Attributes.add("text", cg->getText());
                  else
                    {
                      Attributes.add("originOfText", id);
                    }

#ifdef USE_CRENDER_EXTENSION

                  if (cg->getObjectRole().find_first_not_of(" \t\r\n") != std::string::npos)
                    {
                      Attributes.add("objectRole", cg->getObjectRole());
                    }

#endif // USE_CRENDER_EXTENSION

                  startSaveElement("TextGlyph", Attributes);

                  saveBoundingBox(cg->getBoundingBox());

                  endSaveElement("TextGlyph");
                }
            }

          endSaveElement("ListOfTextGlyphs");
        }

      //additional graphical objects
      if (pLayout->getListOfGraphicalObjects().size() > 0)
        {
          startSaveElement("ListOfAdditionalGraphicalObjects");

          jmax = pLayout->getListOfGraphicalObjects().size();

          for (j = 0; j < jmax; ++j)
            {
              CLGraphicalObject* cg = pLayout->getListOfGraphicalObjects()[j];
              Attributes.erase();
              Attributes.add("key", cg->getKey());
              Attributes.add("name", cg->getObjectName());
#ifdef USE_CRENDER_EXTENSION

              if (cg->getObjectRole().find_first_not_of(" \t\r\n") != std::string::npos)
                {
                  Attributes.add("objectRole", cg->getObjectRole());
                }

#endif // USE_CRENDER_EXTENSION
              startSaveElement("AdditionalGraphicalObject", Attributes);

              saveBoundingBox(cg->getBoundingBox());

              endSaveElement("AdditionalGraphicalObject");
            }

          endSaveElement("ListOfAdditionalGraphicalObjects");
        }

#ifdef USE_CRENDER_EXTENSION

      // save the local render information
      if (pLayout->getListOfLocalRenderInformationObjects().size() > 0)
        {
          saveListOfLocalRenderInformation(pLayout->getListOfLocalRenderInformationObjects());
        }

#endif /* USE_CRENDER_EXTENSION */
      endSaveElement("Layout");
    }

#ifdef USE_CRENDER_EXTENSION

  // save the global render information list
  if (mpLayoutList->getListOfGlobalRenderInformationObjects().size() > 0)
    {
      saveListOfGlobalRenderInformation(mpLayoutList->getListOfGlobalRenderInformationObjects());
    }

#endif /* USE_CRENDER_EXTENSION */
  endSaveElement("ListOfLayouts");

  return success;
}

bool CCopasiXML::saveGUI()
{
  bool success = true;

  if (!haveGUI()) return success;

  startSaveElement("GUI");

  if (mpGUI->getSliderList() && mpGUI->getSliderList()->size())
    {
      startSaveElement("ListOfSliders");

      CSlider * pSlider;
      CXMLAttributeList Attributes;

      Attributes.add("key", "");
      Attributes.add("associatedEntityKey", "");
      Attributes.add("objectCN", "");
      Attributes.add("objectType", "");
      Attributes.add("objectValue", "");
      Attributes.add("minValue", "");
      Attributes.add("maxValue", "");
      Attributes.add("tickNumber", "");
      Attributes.add("tickFactor", "");
      Attributes.add("scaling", "");

      unsigned C_INT32 i, imax = mpGUI->getSliderList()->size();

      for (i = 0; i < imax; i++)
        {
          pSlider = (*mpGUI->getSliderList())[i];
          Attributes.setValue(0, pSlider->getKey());
          Attributes.setValue(1, pSlider->getAssociatedEntityKey());
          Attributes.setValue(2, pSlider->getSliderObjectCN());
          Attributes.setValue(3, CSlider::TypeName[pSlider->getSliderType()]);
          Attributes.setValue(4, pSlider->getSliderValue());
          Attributes.setValue(5, pSlider->getMinValue());
          Attributes.setValue(6, pSlider->getMaxValue());
          Attributes.setValue(7, pSlider->getTickNumber());
          Attributes.setValue(8, pSlider->getTickFactor());
          Attributes.setValue(9, pSlider->convertScaleToScaleName(pSlider->getScaling()));
          saveElement("Slider", Attributes);
        }

      endSaveElement("ListOfSliders");
    }

  endSaveElement("GUI");

  return success;
}

bool CCopasiXML::saveSBMLReference()
{
  assert(this->mpDataModel != NULL);

  if (!this->mpDataModel) return false;

  if (this->mpDataModel->getSBMLFileName() == "" ||
      mSBMLReference.size() == 0)
    return true;

  CXMLAttributeList Attributes;

  std::string SBMLFile = this->mpDataModel->getSBMLFileName();

  if (!CDirEntry::isRelativePath(SBMLFile) &&
      !CDirEntry::makePathRelative(SBMLFile, mFilename))
    SBMLFile = CDirEntry::fileName(SBMLFile);

  Attributes.add("file", SBMLFile);

  startSaveElement("SBMLReference", Attributes);
  Attributes.erase();
  Attributes.add("SBMLid", "");
  Attributes.add("COPASIkey", "");

  std::map<std::string, std::string>::const_iterator it = mSBMLReference.begin();
  std::map<std::string, std::string>::const_iterator end = mSBMLReference.end();

  for (; it != end; ++it)
    {
      Attributes.setValue(0, it->first);
      Attributes.setValue(1, it->second);

      saveElement("SBMLMap", Attributes);
    }

  endSaveElement("SBMLReference");

  return true;
}

bool CCopasiXML::buildFunctionList()
{
  bool success = true;

  CCopasiVectorN< CEvaluationTree > * pFunctionList
  = new CCopasiVectorN< CEvaluationTree >;

  *pFunctionList = CCopasiRootContainer::getFunctionList()->getUsedFunctions(this->mpDataModel->getModel());

  if (!setFunctionList(pFunctionList)) success = false;

  return success;
}

#ifdef USE_CRENDER_EXTENSION

/**
 * Saves the list of global render information objects.
 */
void CCopasiXML::saveListOfGlobalRenderInformation(const CCopasiVector<CLGlobalRenderInformation>& list)
{
  startSaveElement("ListOfGlobalRenderInformation");
  unsigned int i, iMax = list.size();

  for (i = 0; i < iMax; ++i)
    {
      saveGlobalRenderInformation(*list[i]);
    }

  endSaveElement("ListOfGlobalRenderInformation");
}

/**
 * Saves the list of local render information objects.
 */
void CCopasiXML::saveListOfLocalRenderInformation(const CCopasiVector<CLLocalRenderInformation>& list)
{
  startSaveElement("ListOfRenderInformation");
  unsigned int i, iMax = list.size();

  for (i = 0; i < iMax; ++i)
    {
      saveLocalRenderInformation(*list[i]);
    }

  endSaveElement("ListOfRenderInformation");
}

/**
 * Saves a single global render information object.
 */
void CCopasiXML::saveGlobalRenderInformation(const CLGlobalRenderInformation& renderInfo)
{
  // first we create the attributes that are common to all render information objects
  CXMLAttributeList attributes;
  saveRenderInformationAttributes(renderInfo, attributes);
  startSaveElement("RenderInformation", attributes);
  // now we save the definition that are the same for all render information objects
  saveRenderInformationDefinitionElements(renderInfo);
  // last we save the global styles
  unsigned int i, iMax = renderInfo.getNumStyles();

  if (iMax > 0)
    {
      startSaveElement("ListOfStyles");

      for (i = 0; i < iMax; ++i)
        {
          saveGlobalStyle(*(dynamic_cast<const CLGlobalStyle*>(renderInfo.getStyle(i))));
        }

      endSaveElement("ListOfStyles");
    }

  endSaveElement("RenderInformation");
}

/**
 * Saves a single local render information object.
 */
void CCopasiXML::saveLocalRenderInformation(const CLLocalRenderInformation& renderInfo)
{
  // first we create the attributes that are common to all render information objects
  CXMLAttributeList attributes;
  saveRenderInformationAttributes(renderInfo, attributes);
  startSaveElement("RenderInformation", attributes);
  // now we save the definition that are the same for all render information objects
  saveRenderInformationDefinitionElements(renderInfo);
  // last we save the global styles
  unsigned int i, iMax = renderInfo.getNumStyles();

  if (iMax > 0)
    {
      startSaveElement("ListOfStyles");

      for (i = 0; i < iMax; ++i)
        {
          saveLocalStyle(*(dynamic_cast<const CLLocalStyle*>(renderInfo.getStyle(i))));
        }

      endSaveElement("ListOfStyles");
    }

  endSaveElement("RenderInformation");
}

void CCopasiXML::saveRenderInformationAttributes(const CLRenderInformationBase& renderInfo, CXMLAttributeList& attributes)
{
  // save the key
  attributes.add("key", renderInfo.getKey());
  // save the name
  std::string s = renderInfo.getName();
  const char* ws = " \t\n\r";

  if (s.find_first_not_of(ws) != std::string::npos)
    {
      attributes.add("name", s);
    }

  // save the background color
  s = renderInfo.getBackgroundColor();

  if (s.find_first_not_of(ws) != std::string::npos)
    {
      // save the name
      attributes.add("backgroundColor", s);
    }

  // save the reference render information key
  s = renderInfo.getReferenceRenderInformationKey();

  if (s.find_first_not_of(ws) != std::string::npos)
    {
      // save the name
      attributes.add("referenceRenderInformation", s);
    }
}

/**
 * Saves color definitions , gradient definitions and line endings.
 */
void CCopasiXML::saveRenderInformationDefinitionElements(const CLRenderInformationBase& renderInfo)
{
  unsigned int i, iMax = renderInfo.getNumColorDefinitions();

  if (iMax > 0)
    {
      startSaveElement("ListOfColorDefinitions");

      for (i = 0; i < iMax; ++i)
        {
          saveColorDefinition(*(renderInfo.getColorDefinition(i)));
        }

      endSaveElement("ListOfColorDefinitions");
    }

  iMax = renderInfo.getNumGradientDefinitions();

  if (iMax > 0)
    {
      startSaveElement("ListOfGradientDefinitions");
      const CLGradientBase* pGradient = NULL;

      for (i = 0; i < iMax; ++i)
        {
          pGradient = renderInfo.getGradientDefinition(i);

          if (dynamic_cast<const CLRadialGradient*>(pGradient))
            {
              saveRadialGradient(*static_cast<const CLRadialGradient*>(pGradient));
            }
          else if (dynamic_cast<const CLLinearGradient*>(pGradient))
            {
              saveLinearGradient(*static_cast<const CLLinearGradient*>(pGradient));
            }
        }

      endSaveElement("ListOfGradientDefinitions");
    }

  iMax = renderInfo.getNumLineEndings();

  if (iMax > 0)
    {
      startSaveElement("ListOfLineEndings");

      for (i = 0; i < iMax; ++i)
        {
          saveLineEnding(*(renderInfo.getLineEnding(i)));
        }

      endSaveElement("ListOfLineEndings");
    }
}

/**
 * Save a single color definition element.
 */
void CCopasiXML::saveColorDefinition(const CLColorDefinition& color)
{
  CXMLAttributeList attributes;
  attributes.add("id", color.getId());
  attributes.add("value", color.createValueString());
  saveElement("ColorDefinition", attributes);
}

/**
 * Saves a single linear gradient definition.
 */
void CCopasiXML::saveLinearGradient(const CLLinearGradient& gradient)
{
  CXMLAttributeList attributes;
  // first we create the common attributes
  saveGradientAttributes(gradient, attributes);
  // now we add the attributes specific to the radial gradient
  attributes.add("x1", gradient.getXPoint1().toString());
  attributes.add("y1", gradient.getYPoint1().toString());

  if (gradient.getZPoint1() != CLRelAbsVector(0.0, 0.0))
    {
      attributes.add("z1", gradient.getZPoint1().toString());
    }

  attributes.add("x2", gradient.getXPoint2().toString());
  attributes.add("y2", gradient.getYPoint2().toString());

  if (gradient.getZPoint2() != CLRelAbsVector(0.0, 0.0))
    {
      attributes.add("z2", gradient.getZPoint2().toString());
    }

  startSaveElement("LinearGradient", attributes);
  // save the gradient stops
  saveGradientElements(gradient);
  endSaveElement("LinearGradient");
}

/**
 * Saves a single radial gradient definition.
 */
void CCopasiXML::saveRadialGradient(const CLRadialGradient& gradient)
{
  CXMLAttributeList attributes;
  // first we create the common attributes
  saveGradientAttributes(gradient, attributes);
  // now we add the attributes specific to the radial gradient
  attributes.add("cx", gradient.getCenterX().toString());
  attributes.add("cy", gradient.getCenterY().toString());
  attributes.add("cz", gradient.getCenterZ().toString());
  attributes.add("r", gradient.getRadius().toString());
  attributes.add("fx", gradient.getFocalPointX().toString());
  attributes.add("fy", gradient.getFocalPointY().toString());
  attributes.add("fz", gradient.getFocalPointZ().toString());
  startSaveElement("RadialGradient", attributes);
  // save the gradient stops
  saveGradientElements(gradient);
  endSaveElement("RadialGradient");
}

/**
 * Adds the attributes common to radial and linear gradient.
 */
void CCopasiXML::saveGradientAttributes(const CLGradientBase& gradient, CXMLAttributeList& attributes)
{
  attributes.add("id", gradient.getId());

  switch (gradient.getSpreadMethod())
    {
      case CLGradientBase::REPEAT:
        attributes.add("spreadMethod", "repeat");
        break;
      case CLGradientBase::REFLECT:
        attributes.add("spreadMethod", "reflect");
        break;
      case CLGradientBase::PAD:
      default:
        attributes.add("spreadMethod", "pad");
        break;
    }
}

/**
 * Saves the elements that are common to linear and radial gradients.
 */
void CCopasiXML::saveGradientElements(const CLGradientBase& gradient)
{
  unsigned int i, iMax = gradient.getNumGradientStops();

  if (iMax > 0)
    {
      for (i = 0; i < iMax; ++i)
        {
          saveGradientStop(*(gradient.getGradientStop(i)));
        }
    }
}

/**
 * Saves a single gradient stop element.
 */
void CCopasiXML::saveGradientStop(const CLGradientStop& stop)
{
  CXMLAttributeList attributes;
  attributes.add("offset", stop.getOffset().toString());
  attributes.add("stop-color", stop.getStopColor());
  saveElement("Stop", attributes);
}

void CCopasiXML::saveLineEnding(const CLLineEnding& lineEnding)
{
  // create the attributes
  CXMLAttributeList attributes;
  attributes.add("id", lineEnding.getId());
  attributes.add("enableRotationalMapping", lineEnding.getIsEnabledRotationalMapping() ? std::string("true") : std::string("false"));
  startSaveElement("LineEnding", attributes);
  // a line ending has a bounding box
  saveBoundingBox(*(lineEnding.getBoundingBox()));
  // and a group element
  assert(lineEnding.getGroup() != NULL);
  saveGroupElement(*lineEnding.getGroup());
  endSaveElement("LineEnding");
}

/**
 * Saves a single local style element.
 */
void CCopasiXML::saveLocalStyle(const CLLocalStyle& style)
{
  // first we create the attributes
  CXMLAttributeList attributes;
  saveStyleAttributes(style, attributes);

  // now we add the attributes that are specific to local styles
  if (style.getNumKeys() > 0)
    {
      attributes.add("keyList", CLStyle::createStringFromSet(style.getKeyList()));
    }

  startSaveElement("Style", attributes);
  saveStyleElements(style);
  endSaveElement("Style");
}

/**
 * Saves a single local style element.
 */
void CCopasiXML::saveGlobalStyle(const CLGlobalStyle& style)
{
  // first we create the attributes
  CXMLAttributeList attributes;
  saveStyleAttributes(style, attributes);
  startSaveElement("Style", attributes);
  saveStyleElements(style);
  endSaveElement("Style");
}

/**
 * Adds the attributes common to both style types.
 */
void CCopasiXML::saveStyleAttributes(const CLStyle& style, CXMLAttributeList& attributes)
{
  attributes.add("key", style.getKey());

  if (style.getNumRoles() > 0)
    {
      attributes.add("roleList", CLStyle::createStringFromSet(style.getRoleList()));
    }

  if (style.getNumTypes() > 0)
    {
      attributes.add("typeList", CLStyle::createStringFromSet(style.getTypeList()));
    }
}

/**
 * Saves the elements common to both style types.
 */
void CCopasiXML::saveStyleElements(const CLStyle& style)
{
  assert(style.getGroup() != NULL);
  saveGroupElement(*(style.getGroup()));
}

/**
 * Saves a group element.
 */
void CCopasiXML::saveGroupElement(const CLGroup& group)
{
  CXMLAttributeList attributes;
  save2DAttributes(group, attributes);
  saveTextAttributes<CLGroup>(group, attributes);
  saveArrowHeadAttributes<CLGroup>(group, attributes);
  startSaveElement("Group", attributes);
  unsigned int i, iMax = group.getNumElements();

  if (iMax > 0)
    {
      for (i = 0; i < iMax; ++i)
        {
          saveTransformation2D(*dynamic_cast<const CLTransformation2D*>(group.getElement(i)));
        }
    }

  endSaveElement("Group");
}

/**
 * Saves the attributes for a transformation.
 */
void CCopasiXML::saveTransformationAttributes(const CLTransformation2D& transformation, CXMLAttributeList& attributes)
{
  // transformation
  if (!transformation.isIdentityMatrix())
    {
      // check if it is a 2D or a 3D transformation
      if (transformation.is2DTransformation())
        {
          if (transformation.isSetMatrix())
            {
              attributes.add("transform", transformation.get2DTransformationString());
            }
        }
      else
        {
          if (transformation.isSetMatrix())
            {
              attributes.add("transform", transformation.get3DTransformationString());
            }
        }
    }
}

/**
 * Saves the attributes for a 1D element
 */
void CCopasiXML::save1DAttributes(const CLGraphicalPrimitive1D& primitive, CXMLAttributeList& attributes)
{
  // first we go and add the transformation attributes because each 1D element is automatically a transformation element.
  saveTransformationAttributes(primitive, attributes);

  // stroke
  if (primitive.isSetStroke())
    {
      attributes.add("stroke", primitive.getStroke());
    }

  // stroke size
  if (primitive.isSetStrokeWidth())
    {
      std::ostringstream os;
      os << primitive.getStrokeWidth();
      attributes.add("stroke-width", os.str());
    }

  // stroke dash array
  if (primitive.isSetDashArray())
    {
      std::ostringstream os;
      unsigned int i, iMax = primitive.getDashArray().size();
      os << primitive.getDashArray()[0];

      for (i = 1; i < iMax; ++i)
        {
          os << ", " << primitive.getDashArray()[i];
        }

      attributes.add("stroke-dasharray", os.str());
    }
}

/**
 * Saves the attributes for a 2D element
 */
void CCopasiXML::save2DAttributes(const CLGraphicalPrimitive2D& primitive, CXMLAttributeList& attributes)
{
  // first we go and add the 1D attributes because each 2D element is automatically a 1D element.
  save1DAttributes(primitive, attributes);

  // fill
  if (primitive.isSetFill())
    {
      attributes.add("fill", primitive.getFillColor());
    }

  // fill rule
  if (primitive.isSetFillRule())
    {
      switch (primitive.getFillRule())
        {
          case CLGraphicalPrimitive2D::EVENODD:
            attributes.add("fill-rule", "evenodd");
            break;
          case CLGraphicalPrimitive2D::NONZERO:
          default:
            attributes.add("fill-rule", "nonzero");
            break;
        }
    }
}

/**
 * Saves the attributes for a text element.
 * We make this a template so that we can use it for a group as well as a text element.
 */
template<typename TEXTELEMENT>
void CCopasiXML::saveTextAttributes(const TEXTELEMENT& text, CXMLAttributeList& attributes)
{
  // text size
  if (text.isSetFontSize())
    {
      attributes.add("font-size", text.getFontSize().toString());
    }

  // font family
  if (text.isSetFontFamily())
    {
      attributes.add("font-family", text.getFontFamily());
    }

  // font weight
  if (text.isSetFontWeight())
    {
      switch (text.getFontWeight())
        {
          case CLText::WEIGHT_BOLD:
            attributes.add("font-weight", "bold");
            break;
          default:
            break;
        }
    }

  // font style
  if (text.isSetFontStyle())
    {
      switch (text.getFontStyle())
        {
          case CLText::STYLE_ITALIC:
            attributes.add("font-style", "italic");
            break;
          default:
            break;
        }
    }

  // text anchor
  if (text.isSetTextAnchor())
    {
      switch (text.getTextAnchor())
        {
          case CLText::ANCHOR_MIDDLE:
            attributes.add("text-anchor", "middle");
            break;
          case CLText::ANCHOR_END:
            attributes.add("text-anchor", "end");
            break;
          case CLText::ANCHOR_START:
            attributes.add("text-anchor", "start");
            break;
          default:
            break;
        }
    }

  // vertical text anchor
  if (text.isSetVTextAnchor())
    {
      switch (text.getVTextAnchor())
        {
          case CLText::ANCHOR_MIDDLE:
            attributes.add("vtext-anchor", "middle");
            break;
          case CLText::ANCHOR_BOTTOM:
            attributes.add("vtext-anchor", "bottom");
            break;
          case CLText::ANCHOR_TOP:
            attributes.add("vtext-anchor", "top");
            break;
          default:
            break;
        }
    }
}

/**
 * Saves the startHead and endHead attribute as found in group and curves.
 * We write it as a template so that it can be used on curves and group elements.
 */
template<typename HEADELEMENT>
void CCopasiXML::saveArrowHeadAttributes(const HEADELEMENT& element, CXMLAttributeList& attributes)
{
  // start head
  if (element.isSetStartHead())
    {
      attributes.add("startHead", element.getStartHead());
    }

  // start head
  if (element.isSetEndHead())
    {
      attributes.add("endHead", element.getEndHead());
    }
}

/**
 * Saves a class that is subclasses from Transformation2D.
 * This covers images, curves, rectangles, ellipses, polygons, text elements and groups.
 */
void CCopasiXML::saveTransformation2D(const CLTransformation2D& transformation)
{
  if (dynamic_cast<const CLGraphicalPrimitive1D*>(&transformation))
    {
      if (dynamic_cast<const CLRenderCurve*>(&transformation))
        {
          saveRenderCurveElement(static_cast<const CLRenderCurve&>(transformation));
        }
      else if (dynamic_cast<const CLGraphicalPrimitive2D*>(&transformation))
        {
          if (dynamic_cast<const CLRectangle*>(&transformation))
            {
              saveRectangleElement(static_cast<const CLRectangle&>(transformation));
            }
          else if (dynamic_cast<const CLEllipse*>(&transformation))
            {
              saveEllipseElement(static_cast<const CLEllipse&>(transformation));
            }
          else if (dynamic_cast<const CLPolygon*>(&transformation))
            {
              savePolygonElement(static_cast<const CLPolygon&>(transformation));
            }
          else if (dynamic_cast<const CLGroup*>(&transformation))
            {
              saveGroupElement(static_cast<const CLGroup&>(transformation));
            }
          else
            {
              // we should never end up here.
              assert(false);
            }
        }
      else if (dynamic_cast<const CLText*>(&transformation))
        {
          saveRenderTextElement(static_cast<const CLText&>(transformation));
        }
      else
        {
          // we should never end up here.
          assert(false);
        }
    }
  else if (dynamic_cast<const CLImage*>(&transformation))
    {
      saveImageElement(static_cast<const CLImage&>(transformation));
    }
  else
    {
      // we should never end up here.
      assert(false);
    }
}

/**
 * saves a single image element.
 */
void CCopasiXML::saveImageElement(const CLImage& image)
{
  CXMLAttributeList attributes;
  // first we add the attributes for a transformation
  saveTransformationAttributes(image, attributes);
  // now we add the image specific attributes
  attributes.add("x", image.getX().toString());
  attributes.add("y", image.getY().toString());

  if (image.getZ() != CLRelAbsVector(0.0, 0.0))
    {
      attributes.add("z", image.getZ().toString());
    }

  attributes.add("width", image.getWidth().toString());
  attributes.add("height", image.getHeight().toString());
  attributes.add("href", image.getImageReference());
  saveElement("Image", attributes);
}

/**
 * saves a single rectangle element.
 */
void CCopasiXML::saveRectangleElement(const CLRectangle& rectangle)
{
  CXMLAttributeList attributes;
  // first we add the attributes for a 2D object
  save2DAttributes(rectangle, attributes);
  // now we add the rectangle specific attributes
  attributes.add("x", rectangle.getX().toString());
  attributes.add("y", rectangle.getY().toString());

  if (rectangle.getZ() != CLRelAbsVector(0.0, 0.0))
    {
      attributes.add("z", rectangle.getZ().toString());
    }

  attributes.add("width", rectangle.getWidth().toString());
  attributes.add("height", rectangle.getHeight().toString());

  if (rectangle.getRadiusX() != CLRelAbsVector(0.0, 0.0))
    {
      attributes.add("rx", rectangle.getRadiusX().toString());
    }

  if (rectangle.getRadiusY() != CLRelAbsVector(0.0, 0.0))
    {
      attributes.add("ry", rectangle.getRadiusY().toString());
    }

  saveElement("Rectangle", attributes);
}

/**
 * saves a single ellipse element.
 */
void CCopasiXML::saveEllipseElement(const CLEllipse& ellipse)
{
  CXMLAttributeList attributes;
  // first we add the attributes for a 2D object
  save2DAttributes(ellipse, attributes);
  // now we add the ellipse specific attributes
  attributes.add("cx", ellipse.getCX().toString());
  attributes.add("cy", ellipse.getCY().toString());

  if (ellipse.getCZ() != CLRelAbsVector(0.0, 0.0))
    {
      attributes.add("cz", ellipse.getCZ().toString());
    }

  if (ellipse.getRX() != CLRelAbsVector(0.0, 0.0))
    {
      attributes.add("rx", ellipse.getRX().toString());
    }

  if (ellipse.getRY() != CLRelAbsVector(0.0, 0.0))
    {
      attributes.add("ry", ellipse.getRY().toString());
    }

  saveElement("Ellipse", attributes);
}

/**
 * saves a single text element.
 */
void CCopasiXML::saveRenderTextElement(const CLText& text)
{
  CXMLAttributeList attributes;
  // first we add the attributes for a 1D object
  save1DAttributes(text, attributes);
  // save text attributes
  attributes.add("x", text.getX().toString());
  attributes.add("y", text.getY().toString());

  if (text.getZ() != CLRelAbsVector(0.0, 0.0))
    {
      attributes.add("z", text.getZ().toString());
    }

  // next we add the text specific attributes
  saveTextAttributes<CLText>(text, attributes);
  startSaveElement("Text", attributes);
  saveData(text.getText());
  endSaveElement("Text");
}

/**
 * saves a single image element.
 */
void CCopasiXML::savePolygonElement(const CLPolygon& polygon)
{
  CXMLAttributeList attributes;
  // first we add the attributes for a 2D object
  save2DAttributes(polygon, attributes);
  startSaveElement("Polygon", attributes);
  saveCurveElements(*polygon.getListOfElements());
  endSaveElement("Polygon");
}

/**
 * saves a single image element.
 */
void CCopasiXML::saveRenderCurveElement(const CLRenderCurve& curve)
{
  CXMLAttributeList attributes;
  // first we add the attributes for a 1D object
  save1DAttributes(curve, attributes);
  // next we add the arrow head attributes
  saveArrowHeadAttributes<CLRenderCurve>(curve, attributes);
  startSaveElement("Curve", attributes);
  saveCurveElements(*curve.getListOfCurveElements());
  endSaveElement("Curve");
}

/**
 * saves a vector of curve elements. This can be called from the polygon as well as the curve.
 */
void CCopasiXML::saveCurveElements(const std::vector<CLRenderPoint*>& curveElements)
{
  startSaveElement("ListOfElements");
  unsigned int i, iMax = curveElements.size();

  for (i = 0; i < iMax; ++i)
    {
      saveRenderPoint(*curveElements[i]);
    }

  endSaveElement("ListOfElements");
}

/**
 * saves a single render point element.
 */
void CCopasiXML::saveRenderPoint(const CLRenderPoint& point)
{
  CXMLAttributeList attributes;
  attributes.add("x", point.x().toString());
  attributes.add("y", point.y().toString());

  if (point.z() != CLRelAbsVector(0.0, 0.0))
    {
      attributes.add("z", point.z().toString());
    }

  const CLRenderCubicBezier* pCB = dynamic_cast<const CLRenderCubicBezier*>(&point);

  if (pCB != NULL)
    {
      attributes.add("basePoint1_x", pCB->basePoint1_X().toString());
      attributes.add("basePoint1_y", pCB->basePoint1_Y().toString());

      if (pCB->basePoint1_Z() != CLRelAbsVector(0.0, 0.0))
        {
          attributes.add("basePoint1_z", pCB->basePoint1_Z().toString());
        }

      attributes.add("basePoint2_x", pCB->basePoint2_X().toString());
      attributes.add("basePoint2_y", pCB->basePoint2_Y().toString());

      if (pCB->basePoint2_Z() != CLRelAbsVector(0.0, 0.0))
        {
          attributes.add("basePoint2_z", pCB->basePoint2_Z().toString());
        }
    }

  saveElement("Element", attributes);
}

#endif /* USE_CRENDER_EXTENSION */
