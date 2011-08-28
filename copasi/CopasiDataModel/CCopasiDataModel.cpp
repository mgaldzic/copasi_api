// Begin CVS Header
//   $Source: /fs/turing/cvs/copasi_dev/copasi/CopasiDataModel/CCopasiDataModel.cpp,v $
//   $Revision: 1.152.2.1 $
//   $Name: Build-33 $
//   $Author: shoops $
//   $Date: 2010/09/30 17:02:32 $
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

#define USE_LAYOUT 1

#include <sbml/SBMLDocument.h>

#include "copasi.h"
#include "copasiversion.h"
#include "CCopasiDataModel.h"

#include "copasi/report/CCopasiTimer.h"
#include "commandline/COptions.h"
#include "commandline/CConfigurationFile.h"
#include "function/CFunctionDB.h"
#include "model/CModel.h"
#include "optimization/COptTask.h"
#include "parameterFitting/CFitTask.h"
#include "plot/COutputDefinitionVector.h"
#include "report/CKeyFactory.h"
#include "report/CReportDefinitionVector.h"
#include "sbml/CSBMLExporter.h"
#include "sbml/SBMLImporter.h"
#include "sbml/SBMLIncompatibility.h"
#include "scan/CScanTask.h"
#include "elementaryFluxModes/CEFMTask.h"
#include "steadystate/CMCATask.h"
#include "steadystate/CMCAProblem.h"
#include "steadystate/CSteadyStateTask.h"
#include "trajectory/CTrajectoryTask.h"
#ifdef COPASI_TSS
# include "tss/CTSSTask.h"
#endif
#include "sensitivities/CSensTask.h"
#include "tssanalysis/CTSSATask.h"
#include "crosssection/CCrossSectionTask.h"
#include "lyap/CLyapTask.h"
#include "tss/CODEExporter.h"
#include "tss/CODEExporterC.h"
#include "tss/CODEExporterBM.h"
#include "tss/CODEExporterXPPAUT.h"
#include "moieties/CMoietiesTask.h"

#include "utilities/CCopasiException.h"
#include "utilities/CCopasiProblem.h"
#include "utilities/CCopasiTask.h"
#include "utilities/CCopasiVector.h"
#include "utilities/CDirEntry.h"
#include "xml/CCopasiXML.h"

#include "layout/CListOfLayouts.h"
#include "layout/CLayoutInitializer.h"
#include "report/CCopasiRootContainer.h"

#ifdef COPASI_NONLIN_DYN
#include "crosssection/CCrossSectionTask.h"
#endif

CDataModelRenameHandler::CDataModelRenameHandler(CCopasiDataModel* dm)
    : mpDataModel(dm)
{}

bool CDataModelRenameHandler::handle(const std::string & oldCN, const std::string & newCN) const
{
  const std::set<CRegisteredObjectName*> nameSet = CRegisteredObjectName::getSet();

  std::set<CRegisteredObjectName*>::const_iterator it, itEnd = nameSet.end();

  for (it = nameSet.begin(); it != itEnd; ++it)
    {
      // We need to make sure that we not change partial names
      if (((*it)->size() == oldCN.size() ||
           ((*it)->size() > oldCN.size() && (**it)[oldCN.size()] == ',')) &&
          oldCN.compare(0, oldCN.size(), **it, 0, oldCN.size()) == 0)
        {
          (**it).replace(0, oldCN.size(), newCN);
        }
    }

  return true;
}

//********************************************************************

CCopasiDataModel::CCopasiDataModel(const bool withGUI):
    CCopasiContainer("Root", NULL, "CN", CCopasiObject::DataModel),
    COutputHandler(),
    mpModel(NULL),
    mpTaskList(NULL),
    mpReportDefinitionList(NULL),
    mpPlotDefinitionList(NULL),
    mpListOfLayouts(NULL),
    mWithGUI(withGUI),
    mpGUI(NULL),
    mChanged(false),
    mAutoSaveNeeded(false),
    mRenameHandler(this),
    mpCurrentSBMLDocument(NULL),
    mSBMLFileName(""),
#ifdef USE_CRENDER_EXTENSION
    mReferenceDir(""),
#endif // USE_CRENDER_EXTENSION
    pOldMetabolites(new CCopasiVectorS < CMetabOld >)
{

  newModel(NULL, NULL);
  CCopasiObject::setRenameHandler(&mRenameHandler); //TODO where in the constructor should this be called?
  new CCopasiTimer(CCopasiTimer::WALL, this);
  new CCopasiTimer(CCopasiTimer::CPU, this);
}

CCopasiDataModel::CCopasiDataModel(const std::string & name,
                                   const CCopasiContainer * pParent,
                                   const std::string & type,
                                   bool withGUI):
    CCopasiContainer(name, pParent, type, CCopasiObject::DataModel),
    COutputHandler(),
    mpModel(NULL),
    mpTaskList(NULL),
    mpReportDefinitionList(NULL),
    mpPlotDefinitionList(NULL),
    mpListOfLayouts(NULL),
    mWithGUI(withGUI),
    mpGUI(NULL),
    mChanged(false),
    mAutoSaveNeeded(false),
    mRenameHandler(this),
    mpCurrentSBMLDocument(NULL),
    mSBMLFileName(""),
    pOldMetabolites(new CCopasiVectorS < CMetabOld >)
{
  newModel(NULL, NULL);
  CCopasiObject::setRenameHandler(&mRenameHandler); //TODO where in the contructor should this be called?
  new CCopasiTimer(CCopasiTimer::WALL, this);
  new CCopasiTimer(CCopasiTimer::CPU, this);
}

CCopasiDataModel::~CCopasiDataModel()
{
  CCopasiObject::setRenameHandler(NULL);
  pdelete(mpModel);
  pdelete(mpTaskList);
  pdelete(mpReportDefinitionList);
  pdelete(mpPlotDefinitionList);
  pdelete(mpListOfLayouts);
  pdelete(mpGUI);
  pdelete(mpCurrentSBMLDocument);
  pdelete(pOldMetabolites);
}

bool CCopasiDataModel::loadModel(const std::string & fileName, CProcessReport* pProcessReport)
{
  CCopasiMessage::clearDeque();

  std::string PWD;
  COptions::getValue("PWD", PWD);

  std::string FileName = fileName;

  if (CDirEntry::isRelativePath(FileName) &&
      !CDirEntry::makePathAbsolute(FileName, PWD))
    FileName = CDirEntry::fileName(FileName);

  std::ifstream File(utf8ToLocale(FileName).c_str());

  if (File.fail())
    {
      CCopasiMessage Message(CCopasiMessage::RAW,
                             "File error when opening '%s'.",
                             FileName.c_str());
      return false;
    }

  std::string Line;
  File >> Line;

  if (!Line.compare(0, 8, "Version="))
    {
      File.close();
      CReadConfig inbuf(FileName.c_str());

      if (inbuf.getVersion() >= "4")
        {
          CCopasiMessage(CCopasiMessage::ERROR,
                         "Can't handle Gepasi Files with Version>=4.");
          return false;
        }

      newModel(NULL, NULL);
      mpModel->load(inbuf);

      dynamic_cast<CSteadyStateTask *>((*mpTaskList)["Steady-State"])->load(inbuf);
      dynamic_cast<CTrajectoryTask *>((*mpTaskList)["Time-Course"])->load(inbuf);

      mSaveFileName = CDirEntry::dirName(FileName)
                      + CDirEntry::Separator
                      + CDirEntry::baseName(FileName);

      std::string Suffix = CDirEntry::suffix(FileName);

      if (strcasecmp(Suffix.c_str(), ".gps") != 0)
        mSaveFileName += Suffix;

      mSaveFileName += ".cps";
      mSaveFileName = CDirEntry::normalize(mSaveFileName);
#ifdef USE_CRENDER_EXTENSION
      // we have to store the reference directory
      mReferenceDir = CDirEntry::dirName(mSaveFileName);
#endif // USE_CRENDER_EXTENSION
      mSBMLFileName = "";

      pdelete(mpCurrentSBMLDocument);

      this->mCopasi2SBMLMap.clear();
    }
  else if (!Line.find("<?xml") != std::string::npos)
    {
      File.seekg(0, std::ios_base::beg);

      CCopasiXML XML;
      XML.setFunctionList(&CCopasiRootContainer::getFunctionList()->loadedFunctions());
      XML.setDatamodel(this);

      SCopasiXMLGUI *pGUI = NULL;
      std::string SBMLFileNameBkp = mSBMLFileName;

      if (mWithGUI)
        {
          pGUI = new SCopasiXMLGUI("GUI", this);
          XML.setGUI(pGUI);
        }

      // save the copasi2sbml map somewhere and clear it
      std::map<CCopasiObject*, SBase*> mapBackup(mCopasi2SBMLMap);
      mCopasi2SBMLMap.clear();

      try
        {
          if (!XML.load(File, FileName))
            {
              XML.freeModel();
              XML.freeTaskList();
              XML.freeReportList();
              XML.freePlotList();
              XML.freeGUI();
              XML.freeLayoutList();

              // restore the copasi2sbml map
              mCopasi2SBMLMap = mapBackup;
              mSBMLFileName = SBMLFileNameBkp;

              return false;
            }
        }
      catch (...)
        {
          XML.freeModel();
          XML.freeTaskList();
          XML.freeReportList();
          XML.freePlotList();
          XML.freeGUI();
          XML.freeLayoutList();

          // restore the copasi2sbml map
          mCopasi2SBMLMap = mapBackup;
          mSBMLFileName = SBMLFileNameBkp;
          // rethrow the exception so the program flow should still be
          // the same as before
          throw;
        }

      newModel(XML.getModel(), pProcessReport);

      pdelete(mpCurrentSBMLDocument);

      this->mCopasi2SBMLMap.clear();

      if (XML.getTaskList())
        {
          pdelete(mpTaskList);
          mpTaskList = XML.getTaskList();
          mpTaskList->setObjectName("TaskList");
          add(mpTaskList, true);
          addDefaultTasks();

          // We need to initialize all the task so that results are available

          // We suppress all errors and warnings
          unsigned C_INT32 Size = CCopasiMessage::size();

          CCopasiVectorN< CCopasiTask >::iterator it = mpTaskList->begin();
          CCopasiVectorN< CCopasiTask >::iterator end = mpTaskList->end();

          for (; it != end; ++it)
            {
              try
                {
                  (*it)->initialize(CCopasiTask::NO_OUTPUT, NULL, NULL);
                }

              catch (...) {}
            }

          // Remove error messages created by the task initialization as this may fail
          // due to incomplete task specification at this time.
          while (CCopasiMessage::size() > Size)
            CCopasiMessage::getLastMessage();
        }

      if (XML.getReportList())
        {
          pdelete(mpReportDefinitionList);
          mpReportDefinitionList = XML.getReportList();
          add(mpReportDefinitionList, true);
          addDefaultReports();
        }

      if (XML.getPlotList())
        {
          pdelete(mpPlotDefinitionList);
          mpPlotDefinitionList = XML.getPlotList();
          add(mpPlotDefinitionList, true);
        }

      //TODO: layouts
      if (XML.getLayoutList())
        {
          pdelete(mpListOfLayouts);
          mpListOfLayouts = XML.getLayoutList();
          add(mpListOfLayouts, true);
        }

      // for debugging create a template layout
      //mpListOfLayouts->add(CLayoutInitializer::createLayoutFromCModel(mpModel), true);

      if (mWithGUI)
        {
          pdelete(mpGUI);
          mpGUI = pGUI;
        }

      mSaveFileName = CDirEntry::normalize(FileName);
#ifdef USE_CRENDER_EXTENSION
      // we have to store the reference directory
      mReferenceDir = CDirEntry::dirName(mSaveFileName);
#endif // USE_CRENDER_EXTENSION
    }
  else
    {
      CCopasiMessage(CCopasiMessage::ERROR, MCXML + 13, FileName.c_str());
      return false;
    }

  if (mpModel)
    {
      mpModel->compileIfNecessary(pProcessReport);
      mpModel->updateInitialValues();
    }

  changed(false);

  return true;
}

bool CCopasiDataModel::saveModel(const std::string & fileName, CProcessReport* pProcessReport,
                                 bool overwriteFile,
                                 const bool & autoSave)
{
  CCopasiMessage::clearDeque();

  std::string FileName = (fileName != "") ? fileName : mSaveFileName;

  std::string PWD;
  COptions::getValue("PWD", PWD);

  if (CDirEntry::isRelativePath(FileName) &&
      !CDirEntry::makePathAbsolute(FileName, PWD))
    FileName = CDirEntry::fileName(FileName);

  if (CDirEntry::exist(FileName))
    {
      if (!overwriteFile)
        {
          CCopasiMessage(CCopasiMessage::ERROR,
                         MCDirEntry + 1,
                         FileName.c_str());
          return false;
        }

      if (!CDirEntry::isWritable(FileName))
        {
          CCopasiMessage(CCopasiMessage::ERROR,
                         MCDirEntry + 2,
                         FileName.c_str());
          return false;
        }
    }

  try
    {
      if (!mpModel->compileIfNecessary(pProcessReport))
        return false;
    }

  catch (...)
    {
      return false;
    }

  CCopasiXML XML;

  XML.setModel(mpModel);
  XML.setTaskList(mpTaskList);
  XML.setReportList(mpReportDefinitionList);
  XML.setPlotList(mpPlotDefinitionList);
  XML.setGUI(mpGUI);
  XML.setLayoutList(*mpListOfLayouts);
  XML.setDatamodel(this);
  bool success = true;

  if (!autoSave)
    {
      // We are first writing to a temporary file to prevent accidental
      // destruction of an existing file in case the save command fails.
      std::string TmpFileName;
      COptions::getValue("Tmp", TmpFileName);
      TmpFileName = CDirEntry::createTmpName(TmpFileName, ".cps");

      try
        {
          if (!XML.CCopasiXMLInterface::save(TmpFileName, FileName))
            {
              CDirEntry::remove(TmpFileName);
              success = false;
            }
        }

      catch (...)
        {
          CDirEntry::remove(TmpFileName);
          return false;
        }

      if (success && !CDirEntry::move(TmpFileName, FileName))
        success = false;
    }

  if (autoSave || !success)
    {
      try
        {
          if (!XML.CCopasiXMLInterface::save(FileName, FileName))
            return false;
        }

      catch (...)
        {
          return false;
        }
    }

  if (!autoSave)
    {
      changed(false);
      mSaveFileName = CDirEntry::normalize(FileName);
    }

  return true;
}

bool CCopasiDataModel::autoSave()
{
  if (!mAutoSaveNeeded) return true;

  std::string AutoSave;

  COptions::getValue("Tmp", AutoSave);

  if (AutoSave == "") return false;

  AutoSave += CDirEntry::Separator + "tmp_";

  if (mSaveFileName != "")
    AutoSave += CDirEntry::baseName(mSaveFileName);
  else
    AutoSave += "untitled";

  AutoSave += ".cps";

  try
    {
      if (!saveModel(AutoSave, NULL, true, true)) return false;
    }

  catch (...)
    {
      return false;
    }

  mAutoSaveNeeded = false;
  return true;
}

bool CCopasiDataModel::newModel(CModel * pModel, CProcessReport* pProcessReport, CListOfLayouts * pLol)
{
  //deal with the CModel
  pdelete(mpModel);

  if (pModel)
    mpModel = pModel;
  else
    {
      mpModel = new CModel(this);
      mSaveFileName = "";
      mSBMLFileName = "";
#ifdef USE_CRENDER_EXTENSION
      // we have to reset the reference directory
      mReferenceDir = "";
#endif // USE_CRENDER_EXTENSION

      pdelete(mpCurrentSBMLDocument);

      this->mCopasi2SBMLMap.clear();
    }

  //now do the same for the ListOfLayouts
  pdelete(mpListOfLayouts);

  if (pLol)
    mpListOfLayouts = pLol;
  else
    mpListOfLayouts = new CListOfLayouts("ListOflayouts", this);

  pdelete(mpTaskList);
  mpTaskList = new CCopasiVectorN< CCopasiTask >("TaskList", this);

  // We have at least one task of every type
  addDefaultTasks();

  // We need to initialize all the task so that results are available

  // We suppress all errors and warnings
  unsigned C_INT32 Size = CCopasiMessage::size();

  CCopasiVectorN< CCopasiTask >::iterator it = mpTaskList->begin();
  CCopasiVectorN< CCopasiTask >::iterator end = mpTaskList->end();

  for (; it != end; ++it)
    {
      try
        {
          (*it)->initialize(CCopasiTask::NO_OUTPUT, NULL, NULL);
        }

      catch (...) {}
    }

  // Remove error messages created by the task initialization as this may fail
  // due to incomplete task specification at this time.
  while (CCopasiMessage::size() > Size)
    CCopasiMessage::getLastMessage();

  pdelete(mpReportDefinitionList);
  mpReportDefinitionList = new CReportDefinitionVector("ReportDefinitions", this);
  addDefaultReports();

  pdelete(mpPlotDefinitionList);
  mpPlotDefinitionList = new COutputDefinitionVector("OutputDefinitions", this);

  if (mWithGUI)
    {
      pdelete(mpGUI);
      mpGUI = new SCopasiXMLGUI("GUI", this);
    }

  if (mpModel)
    {
      mpModel->compileIfNecessary(pProcessReport);
      mpModel->updateInitialValues();
    }

  changed(false);

  return true;
}

bool CCopasiDataModel::importSBMLFromString(const std::string& sbmlDocumentText, CProcessReport* pImportHandler)
{
  CCopasiMessage::clearDeque();

  SBMLImporter importer;
  // Right now we always import the COPASI MIRIAM annotation if it is there.
  // Later this will be settable by the user in the preferences dialog
  importer.setImportCOPASIMIRIAM(true);
  importer.setImportHandler(pImportHandler);
  //mCopasi2SBMLMap.clear();
  CModel* pModel = NULL;

  SBMLDocument * pSBMLDocument = NULL;
  std::map<CCopasiObject*, SBase*> Copasi2SBMLMap;

  CListOfLayouts * pLol = NULL; //

  try
    {
      pModel = importer.parseSBML(sbmlDocumentText, CCopasiRootContainer::getFunctionList(),
                                  pSBMLDocument, Copasi2SBMLMap, pLol, this);
    }
  catch (CCopasiException except)
    {
      importer.restoreFunctionDB();
      throw except;
    }

  if (pModel == NULL)
    {
      importer.restoreFunctionDB();
      return false;
    }

  pdelete(mpCurrentSBMLDocument);

  mpCurrentSBMLDocument = pSBMLDocument;
  mCopasi2SBMLMap = Copasi2SBMLMap;

  return newModel(pModel, pImportHandler, pLol);
}

bool CCopasiDataModel::importSBML(const std::string & fileName, CProcessReport* pImportHandler)
{
  CCopasiMessage::clearDeque();

  std::string PWD;
  COptions::getValue("PWD", PWD);

  std::string FileName = fileName;

  if (CDirEntry::isRelativePath(FileName) &&
      !CDirEntry::makePathAbsolute(FileName, PWD))
    FileName = CDirEntry::fileName(FileName);

  std::ifstream File(utf8ToLocale(FileName).c_str());

  SBMLImporter importer;
  // Right now we always import the COPASI MIRIAM annotation if it is there.
  // Later this will be settable by the user in the preferences dialog
  importer.setImportCOPASIMIRIAM(true);
  importer.setImportHandler(pImportHandler);
  CModel* pModel = NULL;

  SBMLDocument * pSBMLDocument = NULL;
  std::map<CCopasiObject*, SBase*> Copasi2SBMLMap;

  CListOfLayouts * pLol = NULL;

  try
    {
      pModel = importer.readSBML(FileName, CCopasiRootContainer::getFunctionList(),
                                 pSBMLDocument, Copasi2SBMLMap, pLol, this);
    }
  catch (CCopasiException except)
    {
      importer.restoreFunctionDB();
      throw except;
    }

  if (pModel == NULL)
    {
      importer.restoreFunctionDB();
      return false;
    }

  mSaveFileName = CDirEntry::dirName(FileName)
                  + CDirEntry::Separator
                  + CDirEntry::baseName(FileName);

  std::string Suffix = CDirEntry::suffix(FileName);

  if (strcasecmp(Suffix.c_str(), ".xml") != 0)
    mSaveFileName += Suffix;

  mSaveFileName += ".cps";
  mSaveFileName = CDirEntry::normalize(mSaveFileName);
#ifdef USE_CRENDER_EXTENSION
  // store the reference directory
  mReferenceDir = CDirEntry::dirName(mSaveFileName);
#endif // USE_CRENDER_EXTENSION
  mSBMLFileName = CDirEntry::normalize(FileName);

  pdelete(mpCurrentSBMLDocument);

  mpCurrentSBMLDocument = pSBMLDocument;
  mCopasi2SBMLMap = Copasi2SBMLMap;
  return newModel(pModel, pImportHandler, pLol);
}

std::string CCopasiDataModel::exportSBMLToString(CProcessReport* /*pExportHandler*/, int sbmlLevel, int sbmlVersion)
{
  CCopasiMessage::clearDeque();
  SBMLDocument* pOrigSBMLDocument = NULL;

#if LIBSBML_VERSION >= 40100

  // if we export an L2 model to L3 or vice versa, we have to throw away any prior information
  // about the current sbml document because libsbml does not support the conversion
  // so we need to make sure that all model elements are created from scratch from the corresponding COPASI elements
  if (this->mpCurrentSBMLDocument != NULL &&
      ((this->mpCurrentSBMLDocument->getLevel() > 2 && sbmlLevel < 3) ||
       (this->mpCurrentSBMLDocument->getLevel() < 3 && sbmlLevel > 2)
      )
     )
    {
      pOrigSBMLDocument = this->mpCurrentSBMLDocument;
      this->mpCurrentSBMLDocument = NULL;
    }

#endif // LIBSBML_VERSION


  CSBMLExporter exporter;
  // Per default export COPASIs MIRIAM annotation.
  // This should eventually be determined by a setting in the preferences
  // dialog.
  exporter.setExportCOPASIMIRIAM(true);
  std::string str = exporter.exportModelToString(*this, sbmlLevel, sbmlVersion);

  // only get the new model if it is not a Level 1 model
  // During export to Level 1 the function definitions have been deleted and therefore
  // all information assiociated with the function definitions will be gone if the user exports
  // to Level 2 after having exported to Level 1
  // This is actual vital to get around Bug 1086 as well.
  // Once I have a Level 1 model, all calls to setName on an
  // SBML object in that model also resets the id, which does not work with the current exporter
  if ((sbmlLevel != 1 || mpCurrentSBMLDocument == NULL) && pOrigSBMLDocument == NULL)
    {
      if (mpCurrentSBMLDocument != exporter.getSBMLDocument())
        {
          pdelete(mpCurrentSBMLDocument);
        }

      // disown the SBML Document from the exporter so we don't have to copy it
      exporter.disownSBMLDocument();
      mpCurrentSBMLDocument = exporter.getSBMLDocument();
      // we also need to get the new copasi2sbml map otherwise it contains invalid pointers
      // since the objects
      this->mCopasi2SBMLMap.clear();
      std::map<const CCopasiObject*, SBase*>::const_iterator it = exporter.getCOPASI2SBMLMap().begin();
      std::map<const CCopasiObject*, SBase*>::const_iterator endit = exporter.getCOPASI2SBMLMap().end();

      while (it != endit)
        {
          this->mCopasi2SBMLMap.insert(std::pair<CCopasiObject*, SBase*>(const_cast<CCopasiObject*>(it->first), it->second));
          ++it;
        }
    }
  // if we have saved the original SBML model somewhere
  // we have to reset it
  else if (pOrigSBMLDocument != NULL)
    {
      this->mpCurrentSBMLDocument = pOrigSBMLDocument;
    }


  return str;
}

bool CCopasiDataModel::exportSBML(const std::string & fileName, bool overwriteFile, int sbmlLevel, int sbmlVersion, bool /*exportIncomplete*/, bool exportCOPASIMIRIAM, CProcessReport* pExportHandler)
{
  CCopasiMessage::clearDeque();

  if (fileName == "") return false;

  std::string PWD;
  COptions::getValue("PWD", PWD);

  std::string FileName = fileName;

  if (CDirEntry::isRelativePath(FileName) &&
      !CDirEntry::makePathAbsolute(FileName, PWD))
    FileName = CDirEntry::fileName(FileName);

  if (CDirEntry::exist(FileName))
    {
      if (!overwriteFile)
        {
          CCopasiMessage(CCopasiMessage::ERROR,
                         MCDirEntry + 1,
                         FileName.c_str());
          return false;
        }

      if (!CDirEntry::isWritable(FileName))
        {
          CCopasiMessage(CCopasiMessage::ERROR,
                         MCDirEntry + 2,
                         FileName.c_str());
          return false;
        }
    }

  try
    {
      if (!mpModel->compileIfNecessary(pExportHandler))
        return false;
    }

  catch (...)
    {
      return false;
    }

  CSBMLExporter exporter;
  exporter.setExportCOPASIMIRIAM(exportCOPASIMIRIAM);
  SBMLDocument* pOrigSBMLDocument = NULL;

#if LIBSBML_VERSION >= 40100

  // if we export an L2 model to L3 or vice versa, we have to throw away any prior information
  // about the current sbml document because libsbml does not support the conversion
  // so we need to make sure that all model elements are created from scratch from the corresponding COPASI elements
  if (this->mpCurrentSBMLDocument != NULL &&
      ((this->mpCurrentSBMLDocument->getLevel() > 2 && sbmlLevel < 3) ||
       (this->mpCurrentSBMLDocument->getLevel() < 3 && sbmlLevel > 2)
      )
     )
    {
      pOrigSBMLDocument = this->mpCurrentSBMLDocument;
      this->mpCurrentSBMLDocument = NULL;
    }

#endif // LIBSBML_VERSION

  //exporter.setExportHandler(pExportHandler);
  if (!exporter.exportModel(*this, FileName, sbmlLevel, sbmlVersion, overwriteFile)) return false;

  // only get the new model if it is not a Level 1 model
  // During export to Level 1 the function definitions have been deleted and therefore
  // all information assiociated with the function definitions will be gone if the user exports
  // to Level 2 after having exported to Level 1
  // This is actual vital to get around Bug 1086 as well.
  // Once I have a Level 1 model, all calls to setName on an
  // SBML object in that model also resets the id, which does not work with the current exporter
  if ((sbmlLevel != 1 || mpCurrentSBMLDocument == NULL) && pOrigSBMLDocument == NULL)
    {

      if (mpCurrentSBMLDocument != exporter.getSBMLDocument())
        pdelete(mpCurrentSBMLDocument);

      // disown the SBML Document from the exporter so we don't have to copy it
      exporter.disownSBMLDocument();
      mpCurrentSBMLDocument = exporter.getSBMLDocument();
      // we also need to get the new copasi2sbml map otherwise it contains invalid pointers
      // since the objects
      this->mCopasi2SBMLMap.clear();
      std::map<const CCopasiObject*, SBase*>::const_iterator it = exporter.getCOPASI2SBMLMap().begin();
      std::map<const CCopasiObject*, SBase*>::const_iterator endit = exporter.getCOPASI2SBMLMap().end();

      while (it != endit)
        {
          this->mCopasi2SBMLMap.insert(std::pair<CCopasiObject*, SBase*>(const_cast<CCopasiObject*>(it->first), it->second));
          ++it;
        }
    }
  // if we have saved the original SBML model somewhere
  // we have to reset it
  else if (pOrigSBMLDocument != NULL)
    {
      this->mpCurrentSBMLDocument = pOrigSBMLDocument;
    }

  mSBMLFileName = FileName;
  return true;
}

bool CCopasiDataModel::exportMathModel(const std::string & fileName, CProcessReport* pProcessReport,
                                       const std::string & filter, bool overwriteFile)
{
  CCopasiMessage::clearDeque();

  if (fileName == "") return false;

  if (CDirEntry::exist(fileName))
    {
      if (!overwriteFile)
        {
          CCopasiMessage(CCopasiMessage::ERROR,
                         MCDirEntry + 1,
                         fileName.c_str());
          return false;
        }

      if (!CDirEntry::isWritable(fileName))
        {
          CCopasiMessage(CCopasiMessage::ERROR,
                         MCDirEntry + 2,
                         fileName.c_str());
          return false;
        }
    }

  try
    {
      if (!mpModel->compileIfNecessary(pProcessReport))
        return false;
    }

  catch (...)
    {
      return false;
    }

  CCopasiVector< CModelValue >::const_iterator it = mpModel->getModelValues().begin();
  CCopasiVector< CModelValue >::const_iterator end = mpModel->getModelValues().end();

  for (; it != end; ++it)
    if ((*it)->isUsed()) break;

  if (it != end)
    CCopasiMessage(CCopasiMessage::WARNING, MCODEExporter + 2);

  if (filter == "C Files (*.c)")
    {
      CODEExporterC exporter;

      return exporter.exportMathModel(this, fileName.c_str(), filter.c_str(), overwriteFile);
    }

  if (filter == "Berkeley Madonna Files (*.mmd)")
    {
      CODEExporterBM exporter;

      return exporter.exportMathModel(this, fileName.c_str(), filter.c_str(), overwriteFile);
    }

  if (filter == "XPPAUT (*.ode)")
    {
      CODEExporterXPPAUT exporter;

      return exporter.exportMathModel(this, fileName.c_str(), filter.c_str(), overwriteFile);
    }

  return false;
}

const CModel * CCopasiDataModel::getModel() const
{return mpModel;}

CModel * CCopasiDataModel::getModel()
{return mpModel;}

CCopasiVectorN< CCopasiTask > * CCopasiDataModel::getTaskList()
{return mpTaskList;}

const CCopasiVectorN< CCopasiTask > * CCopasiDataModel::getTaskList() const
{return mpTaskList;}

CCopasiTask * CCopasiDataModel::addTask(const CCopasiTask::Type & taskType)
{
  CCopasiTask * pTask = NULL;

  switch (taskType)
    {
      case CCopasiTask::steadyState:
        pTask = new CSteadyStateTask(mpTaskList);
        break;

      case CCopasiTask::timeCourse:
        pTask = new CTrajectoryTask(mpTaskList);
        break;

      case CCopasiTask::scan:
        pTask = new CScanTask(mpTaskList);
        break;

      case CCopasiTask::fluxMode:
        pTask = new CEFMTask(mpTaskList);
        break;

      case CCopasiTask::optimization:
        pTask = new COptTask(taskType, mpTaskList);
        break;

      case CCopasiTask::parameterFitting:
        pTask = new CFitTask(taskType, mpTaskList);
        break;

      case CCopasiTask::mca:
        pTask = new CMCATask(mpTaskList);
        static_cast< CMCAProblem * >(pTask->getProblem())->setSteadyStateRequested(true);
        break;

      case CCopasiTask::lyap:
        pTask = new CLyapTask(mpTaskList);
        break;

#ifdef COPASI_TSS
      case CCopasiTask::tss:
        pTask = new CTSSTask(mpTaskList);
        break;
#endif // COPASI_TSS

      case CCopasiTask::sens:
        pTask = new CSensTask(mpTaskList);
        break;

      case CCopasiTask::tssAnalysis:
        pTask = new CTSSATask(mpTaskList);
        break;

      case CCopasiTask::moieties:
        pTask = new CMoietiesTask(taskType, mpTaskList);
        break;

#ifdef COPASI_NONLIN_DYN
      case CCopasiTask::crosssection:
        pTask = new CCrossSectionTask(mpTaskList);
        break;
#endif

      default:
        return pTask;
    }

  pTask->getProblem()->setModel(mpModel);
  mpTaskList->add(pTask);

  return pTask;
}

bool CCopasiDataModel::addDefaultTasks()
{
  unsigned C_INT32 i;

  for (i = 0; CCopasiTask::TypeName[i] != ""; i++)
    if (mpTaskList->getIndex(CCopasiTask::TypeName[i]) == C_INVALID_INDEX)
      addTask((CCopasiTask::Type) i);

  return true;
}

bool CCopasiDataModel::appendDependentTasks(std::set< const CCopasiObject * > candidates,
    std::set< const CCopasiObject * > & dependentTasks) const
{
  size_t Size = dependentTasks.size();

  std::set< const CCopasiObject * >::const_iterator it = candidates.begin();
  std::set< const CCopasiObject * >::const_iterator end = candidates.end();

  CCopasiVectorN< CCopasiTask >::const_iterator itTask = mpTaskList->begin();
  CCopasiVectorN< CCopasiTask >::const_iterator endTask = mpTaskList->end();


  for (; it != end; ++it)
    {
      const CReportDefinition * pReportDefinition = dynamic_cast< const CReportDefinition * >(*it);

      if (pReportDefinition == NULL)
        continue;

      itTask = mpTaskList->begin();

      for (; itTask != endTask; ++itTask)
        {
          if ((*itTask)->getReport().getReportDefinition() == pReportDefinition)
            {
              dependentTasks.insert(*itTask);
            }
        }
    }

  return Size < dependentTasks.size();
}

CReportDefinition * CCopasiDataModel::addReport(const CCopasiTask::Type & taskType)
{
  CReportDefinition * pReport = NULL;

  switch (taskType)
    {
      case CCopasiTask::steadyState:
        pReport = new CReportDefinition(CCopasiTask::TypeName[taskType]);
        pReport->setTaskType(taskType);
        pReport->setComment("Automatically generated report.");
        pReport->setIsTable(false);
        pReport->setSeparator(CCopasiReportSeparator("\t"));
        pReport->getFooterAddr()->push_back(CCopasiObjectName("CN=Root,Vector=TaskList[Steady-State]"));
        break;

      case CCopasiTask::timeCourse:
        // No default report available.
        break;

      case CCopasiTask::scan:
        // No default report available.
        break;

      case CCopasiTask::fluxMode:
        pReport = new CReportDefinition(CCopasiTask::TypeName[taskType]);
        pReport->setTaskType(taskType);
        pReport->setComment("Automatically generated report.");
        pReport->setIsTable(false);
        pReport->setSeparator(CCopasiReportSeparator("\t"));
        pReport->getFooterAddr()->push_back(CCopasiObjectName("CN=Root,Vector=TaskList[Elementary Flux Modes],Object=Result"));
        break;

      case CCopasiTask::optimization:
        pReport = new CReportDefinition(CCopasiTask::TypeName[taskType]);
        pReport->setTaskType(taskType);
        pReport->setComment("Automatically generated report.");
        pReport->setIsTable(false);
        pReport->setTitle(false);
        pReport->setSeparator(CCopasiReportSeparator("\t"));

        // Header
        pReport->getHeaderAddr()->push_back(CCopasiObjectName("CN=Root,Vector=TaskList[Optimization],Object=Description"));
        pReport->getHeaderAddr()->push_back(CCopasiObjectName("String=\\[Function Evaluations\\]"));
        pReport->getHeaderAddr()->push_back(CCopasiObjectName("Separator=\t"));
        pReport->getHeaderAddr()->push_back(CCopasiObjectName("String=\\[Best Value\\]"));
        pReport->getHeaderAddr()->push_back(CCopasiObjectName("Separator=\t"));
        pReport->getHeaderAddr()->push_back(CCopasiObjectName("String=\\[Best Parameters\\]"));

        // Body
        pReport->getBodyAddr()->push_back(CCopasiObjectName("CN=Root,Vector=TaskList[Optimization],Problem=Optimization,Reference=Function Evaluations"));
        pReport->getBodyAddr()->push_back(CCopasiObjectName("Separator=\t"));
        pReport->getBodyAddr()->push_back(CCopasiObjectName("CN=Root,Vector=TaskList[Optimization],Problem=Optimization,Reference=Best Value"));
        pReport->getBodyAddr()->push_back(CCopasiObjectName("Separator=\t"));
        pReport->getBodyAddr()->push_back(CCopasiObjectName("CN=Root,Vector=TaskList[Optimization],Problem=Optimization,Reference=Best Parameters"));

        // Footer
        pReport->getFooterAddr()->push_back(CCopasiObjectName("String=\n"));
        pReport->getFooterAddr()->push_back(CCopasiObjectName("CN=Root,Vector=TaskList[Optimization],Object=Result"));
        break;

        //**************************************************************************
      case CCopasiTask::parameterFitting:
        pReport = new CReportDefinition(CCopasiTask::TypeName[taskType]);
        pReport->setTaskType(taskType);
        pReport->setComment("Automatically generated report.");
        pReport->setIsTable(false);
        pReport->setTitle(false);
        pReport->setSeparator(CCopasiReportSeparator("\t"));

        // Header
        pReport->getHeaderAddr()->push_back(CCopasiObjectName("CN=Root,Vector=TaskList[Parameter Estimation],Object=Description"));
        pReport->getHeaderAddr()->push_back(CCopasiObjectName("String=\\[Function Evaluations\\]"));
        pReport->getHeaderAddr()->push_back(CCopasiObjectName("Separator=\t"));
        pReport->getHeaderAddr()->push_back(CCopasiObjectName("String=\\[Best Value\\]"));
        pReport->getHeaderAddr()->push_back(CCopasiObjectName("Separator=\t"));
        pReport->getHeaderAddr()->push_back(CCopasiObjectName("String=\\[Best Parameters\\]"));

        // Body
        pReport->getBodyAddr()->push_back(CCopasiObjectName("CN=Root,Vector=TaskList[Parameter Estimation],Problem=Parameter Estimation,Reference=Function Evaluations"));
        pReport->getBodyAddr()->push_back(CCopasiObjectName("Separator=\t"));
        pReport->getBodyAddr()->push_back(CCopasiObjectName("CN=Root,Vector=TaskList[Parameter Estimation],Problem=Parameter Estimation,Reference=Best Value"));
        pReport->getBodyAddr()->push_back(CCopasiObjectName("Separator=\t"));
        pReport->getBodyAddr()->push_back(CCopasiObjectName("CN=Root,Vector=TaskList[Parameter Estimation],Problem=Parameter Estimation,Reference=Best Parameters"));

        // Footer
        pReport->getFooterAddr()->push_back(CCopasiObjectName("String=\n"));
        pReport->getFooterAddr()->push_back(CCopasiObjectName("CN=Root,Vector=TaskList[Parameter Estimation],Object=Result"));
        break;

        //**************************************************************************
      case CCopasiTask::lyap:
        pReport = new CReportDefinition(CCopasiTask::TypeName[taskType]);
        pReport->setTaskType(taskType);
        pReport->setComment("Automatically generated report.");
        pReport->setIsTable(false);
        pReport->setTitle(false);
        pReport->setSeparator(CCopasiReportSeparator("\t"));

        // Header
        pReport->getHeaderAddr()->push_back(CCopasiObjectName("CN=Root,Vector=TaskList[Lyapunov Exponents],Object=Description"));

        // Footer
        pReport->getFooterAddr()->push_back(CCopasiObjectName("String=\n"));
        pReport->getFooterAddr()->push_back(CCopasiObjectName("CN=Root,Vector=TaskList[Lyapunov Exponents],Object=Result"));
        break;

        //**************************************************************************
      case CCopasiTask::mca:
        pReport = new CReportDefinition(CCopasiTask::TypeName[taskType]);
        pReport->setTaskType(taskType);
        pReport->setComment("Automatically generated report.");
        pReport->setIsTable(false);
        pReport->setTitle(false);
        pReport->setSeparator(CCopasiReportSeparator("\t"));

        // Header
        pReport->getHeaderAddr()->push_back(CCopasiObjectName("CN=Root,Vector=TaskList[Metabolic Control Analysis],Object=Description"));

        // Footer
        pReport->getFooterAddr()->push_back(CCopasiObjectName("String=\n"));
        pReport->getFooterAddr()->push_back(CCopasiObjectName("CN=Root,Vector=TaskList[Metabolic Control Analysis],Object=Result"));
        break;

        //**************************************************************************
      case CCopasiTask::sens:
        pReport = new CReportDefinition(CCopasiTask::TypeName[taskType]);
        pReport->setTaskType(taskType);
        pReport->setComment("Automatically generated report.");
        pReport->setIsTable(false);
        pReport->setTitle(false);
        pReport->setSeparator(CCopasiReportSeparator("\t"));

        // Header
        pReport->getHeaderAddr()->push_back(CCopasiObjectName("CN=Root,Vector=TaskList[Sensitivities],Object=Description"));

        // Footer
        pReport->getFooterAddr()->push_back(CCopasiObjectName("String=\n"));
        pReport->getFooterAddr()->push_back(CCopasiObjectName("CN=Root,Vector=TaskList[Sensitivities],Object=Result"));
        break;

        //**************************************************************************
      case CCopasiTask::tssAnalysis:
        pReport = new CReportDefinition(CCopasiTask::TypeName[taskType]);
        pReport->setTaskType(taskType);
        pReport->setComment("Automatically generated report.");
        pReport->setIsTable(false);
        pReport->setTitle(false);
        pReport->setSeparator(CCopasiReportSeparator("\t"));

        // Header
        pReport->getHeaderAddr()->push_back(CCopasiObjectName("CN=Root,Vector=TaskList[Time Scale Separation Analysis],Object=Description"));

        // Footer
        pReport->getFooterAddr()->push_back(CCopasiObjectName("String=\n"));
        pReport->getFooterAddr()->push_back(CCopasiObjectName("CN=Root,Vector=TaskList[Time Scale Separation Analysis],Object=Result"));
        break;

      default:
        return pReport;
    }

  if (pReport) mpReportDefinitionList->add(pReport, true);

  return pReport;
}

bool CCopasiDataModel::addDefaultReports()
{
  unsigned C_INT32 i;

  for (i = 0; CCopasiTask::TypeName[i] != ""; i++)
    {
      //try to create the report if it doesn't exist
      if (mpReportDefinitionList->getIndex(CCopasiTask::TypeName[i]) == C_INVALID_INDEX)
        {
          addReport((CCopasiTask::Type) i);
        }

      //see if the report exists now
      CReportDefinition* pReportDef = NULL;

      if (mpReportDefinitionList->getIndex(CCopasiTask::TypeName[i]) != C_INVALID_INDEX)
        pReportDef = (*mpReportDefinitionList)[CCopasiTask::TypeName[i]];

      //see if the task exists
      CCopasiTask* pTask = NULL;

      if (mpTaskList->getIndex(CCopasiTask::TypeName[i]) != C_INVALID_INDEX)
        pTask = (*mpTaskList)[CCopasiTask::TypeName[i]];

      if (pTask && pReportDef) //task and report definition exist
        {
          //if there is no report definition set the default
          if (!pTask->getReport().getReportDefinition())
            {
              pTask->getReport().setReportDefinition(pReportDef);
            }

          //TODO: also set the default report if no target file is set
          //even if a report is already set?
        }
    }

  return true;
}

CReportDefinitionVector * CCopasiDataModel::getReportDefinitionList()
{return mpReportDefinitionList;}

COutputDefinitionVector * CCopasiDataModel::getPlotDefinitionList()
{return mpPlotDefinitionList;}

const CListOfLayouts * CCopasiDataModel::getListOfLayouts() const
{return mpListOfLayouts;}

CListOfLayouts * CCopasiDataModel::getListOfLayouts()
{return mpListOfLayouts;}

SCopasiXMLGUI * CCopasiDataModel::getGUI()
{return mpGUI;}

const std::string & CCopasiDataModel::getFileName() const
{return mSaveFileName;}

bool CCopasiDataModel::isChanged() const
{return mChanged;}

void CCopasiDataModel::changed(const bool & changed)
{
  mChanged = changed;
  mAutoSaveNeeded = changed;
}

SBMLDocument* CCopasiDataModel::getCurrentSBMLDocument()
{
  return this->mpCurrentSBMLDocument;
}

bool CCopasiDataModel::setSBMLFileName(const std::string & fileName)
{
  mSBMLFileName = CDirEntry::normalize(fileName);

  if (CDirEntry::isRelativePath(mSBMLFileName) &&
      !CDirEntry::makePathAbsolute(mSBMLFileName, mSaveFileName))
    mSBMLFileName = CDirEntry::fileName(mSBMLFileName);

  return true;
}

const std::string & CCopasiDataModel::getSBMLFileName() const
{return mSBMLFileName;}

std::map<CCopasiObject*, SBase*>& CCopasiDataModel::getCopasi2SBMLMap()
{
  return this->mCopasi2SBMLMap;
}

void CCopasiDataModel::removeSBMLIdFromFunctions()
{
  CFunctionDB* pFunDB = CCopasiRootContainer::getFunctionList();
  unsigned int i, iMax = pFunDB->loadedFunctions().size();

  for (i = 0; i < iMax; ++i)
    {
      pFunDB->loadedFunctions()[i]->setSBMLId("");
    }
}

bool CCopasiDataModel::removeLayout(const std::string & key)
{
  CLayout *pLayout =
    dynamic_cast< CLayout * >(CCopasiRootContainer::getKeyFactory()->get(key));

  if (!pLayout)
    return false;

  //Check if Layout with that name exists
  unsigned C_INT32 index =
    mpListOfLayouts->CCopasiVector< CLayout >::getIndex(pLayout);

  if (index == C_INVALID_INDEX)
    return false;

  mpListOfLayouts->CCopasiVector< CLayout >::remove(index);

  return true;
}

const CCopasiObject * CCopasiDataModel::ObjectFromName(const std::vector< CCopasiContainer * > & listOfContainer,
    const CCopasiObjectName & objName) const
{
  const CCopasiObject * pObject = NULL;
  const CCopasiContainer* pContainer;
  CCopasiObjectName ContainerName;
  unsigned C_INT32 containerIndex;
  std::string::size_type pos;

  //favor to search the list of container first
  for (containerIndex = 0;
       containerIndex < listOfContainer.size() && !pObject;
       containerIndex++)
    {
      pContainer = listOfContainer[containerIndex];
      ContainerName = pContainer->getCN();

      while (ContainerName.getRemainder() != "")
        {
          ContainerName = ContainerName.getRemainder();
        }

      if ((pos = objName.find(ContainerName)) == std::string::npos)
        continue;

      if (pos + ContainerName.length() == objName.length())
        pObject = pContainer;
      else
        pObject = pContainer->getObject(objName.substr(pos + ContainerName.length() + 1));
    }

  // If we have not found the object in the context we search the whole data model.
  if (pObject == NULL)
    {
      pObject = getObject(objName);
    }

  // if still not found search the function database in the root container
  if (!pObject)
    pObject = CCopasiRootContainer::getFunctionList()->getObject(objName);

  return pObject;
}

CCopasiObject * CCopasiDataModel::ObjectFromName(const std::vector< CCopasiContainer * > & listOfContainer,
    const CCopasiObjectName & objName)
{
  return const_cast<CCopasiObject *>(const_cast<const CCopasiDataModel*>(this)->ObjectFromName(listOfContainer, objName));
}

#ifdef USE_CRENDER_EXTENSION
const std::string& CCopasiDataModel::getReferenceDirectory() const
{
  return this->mReferenceDir;
}
#endif // USE_CRENDER_EXTENSION
