// Begin CVS Header
//   $Source: /fs/turing/cvs/copasi_dev/copasi/report/CReport.cpp,v $
//   $Revision: 1.62 $
//   $Name: Build-33 $
//   $Author: shoops $
//   $Date: 2010/07/16 19:02:49 $
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

#include "copasi.h"

#include "CReportDefinition.h"
#include "CReport.h"
#include "CCopasiContainer.h"
#include "CCopasiTimer.h"

#include "CopasiDataModel/CCopasiDataModel.h"
#include "report/CCopasiRootContainer.h"
#include "utilities/CCopasiMessage.h"
#include "utilities/CDirEntry.h"
#include "utilities/utility.h"

//////////////////////////////////////////////////
//
//class CReport
//
//////////////////////////////////////////////////
CReport::CReport():
    COutputInterface(),
    mpDataModel(NULL),
    mpOstream(NULL),
    mStreamOwner(false),
    mpReportDef(NULL),
    mTarget(""),
    mAppend(true),
    mFooterObjectList(),
    mBodyObjectList(),
    mHeaderObjectList(),
    mpHeader(NULL),
    mpBody(NULL),
    mpFooter(NULL),
    mState(Invalid)
{}

CReport::CReport(const CReport & src):
    COutputInterface(src),
    mpDataModel(src.mpDataModel),
    mpOstream(src.mpOstream),
    mStreamOwner(false),
    mpReportDef(src.mpReportDef),
    mTarget(src.mTarget),
    mAppend(src.mAppend),
    mFooterObjectList(src.mFooterObjectList),
    mBodyObjectList(src.mBodyObjectList),
    mHeaderObjectList(src.mHeaderObjectList),
    mpHeader(src.mpHeader),
    mpBody(src.mpBody),
    mpFooter(src.mpFooter),
    mState(Invalid)
{}

CReport::~CReport()
{cleanup();}

void CReport::cleanup()
{
  mHeaderObjectList.clear();
  mBodyObjectList.clear();
  mFooterObjectList.clear();

  finish();
}

CReportDefinition* CReport::getReportDefinition()
{return mpReportDef;}

void CReport::setReportDefinition(CReportDefinition* reportDef)
{mpReportDef = reportDef;}

const std::string& CReport::getTarget() const
{return mTarget;}

void CReport::setTarget(std::string target)
{mTarget = target;}

bool CReport::append() const
{return mAppend;}

void CReport::setAppend(bool append)
{mAppend = append;}

void CReport::output(const Activity & activity)
{
  switch (activity)
    {
      case COutputInterface::BEFORE:
        printHeader();
        break;

      case COutputInterface::DURING:
        printBody();
        break;

      case COutputInterface::AFTER:
        printFooter();
        break;
    }
}

void CReport::separate(const Activity & /* activity */)
{
  if (!mpOstream) return;

  (*mpOstream) << std::endl;
}

void CReport::finish()
{
  mState = FooterFooter;

  printFooter();

  if (mStreamOwner) pdelete(mpOstream);

  mpOstream = NULL;
  mStreamOwner = false;

  pdelete(mpHeader);
  pdelete(mpBody);
  pdelete(mpFooter);

  mState = Invalid;
}

void CReport::printHeader()
{
  if (!mpOstream) return;

  if (mpHeader)
    switch (mState)
      {
        case Compiled:
          mpHeader->printHeader();
          mState = HeaderHeader;
          return;

        case HeaderHeader:
          mpHeader->printBody();
          mState = HeaderBody;
          return;

        case HeaderBody:
          mpHeader->printBody();
          return;

        case HeaderFooter:
          mpHeader->printFooter();
          return;

        default:
          return;
      }

  if (mState == HeaderFooter) return;

  mState = HeaderFooter;

  std::vector< CCopasiObject * >::iterator it = mHeaderObjectList.begin();
  std::vector< CCopasiObject * >::iterator end = mHeaderObjectList.end();

  if (it == end) return;

  for (; it != end; ++it)(*it)->print(mpOstream);

  (*mpOstream) << std::endl;
}

void CReport::printBody()
{
  if (!mpOstream) return;

  // Close the header part
  if (mState < HeaderFooter)
    {
      mState = HeaderFooter;

      if (mpHeader) mpHeader->printFooter();
    }

  if (mpBody)
    switch (mState)
      {
        case HeaderFooter:
          mpBody->printHeader();
          mState = BodyHeader;
          return;

        case BodyHeader:
          mpBody->printBody();
          mState = BodyBody;
          return;

        case BodyBody:
          mpBody->printBody();
          return;

        case BodyFooter:
          mpBody->printFooter();
          return;

        default:
          return;
      }

  if (mState == BodyFooter) return;

  mState = BodyBody;

  std::vector< CCopasiObject * >::iterator it = mBodyObjectList.begin();
  std::vector< CCopasiObject * >::iterator end = mBodyObjectList.end();

  if (it == end) return;

  for (; it != end; ++it)(*it)->print(mpOstream);

  (*mpOstream) << std::endl;
}

void CReport::printFooter()
{
  if (!mpOstream) return;

  // Close the body part
  if (mState < BodyFooter)
    {
      mState = BodyFooter;

      if (mpBody) mpBody->printFooter();
    }

  if (mpFooter)
    switch (mState)
      {
        case BodyFooter:
          mpFooter->printHeader();
          mState = FooterHeader;
          return;

        case FooterHeader:
          mpFooter->printBody();
          mState = FooterBody;
          return;

        case FooterBody:
          mpFooter->printBody();
          return;

        case FooterFooter:
          mpFooter->printFooter();
          return;

        default:
          return;
      }

  if (mState != FooterFooter) return;

  std::vector< CCopasiObject * >::iterator it = mFooterObjectList.begin();
  std::vector< CCopasiObject * >::iterator end = mFooterObjectList.end();

  if (it == end) return;

  for (; it != end; ++it)(*it)->print(mpOstream);

  (*mpOstream) << std::endl;
}

// Compile the List of Report Objects;
// Support Parellel

bool CReport::compile(std::vector< CCopasiContainer * > listOfContainer,
                      const CCopasiDataModel* pDataModel)
{
  assert(mpDataModel == pDataModel);

  bool success = true;
  COutputInterface::mObjects.clear();

  // check if there is a Report Definition Defined
  if (!mpReportDef) return false;

  if (mpReportDef->isTable())
    if (!mpReportDef->preCompileTable(listOfContainer)) success = false;

  generateObjectsFromName(&listOfContainer, mHeaderObjectList, mpHeader,
                          mpReportDef->getHeaderAddr());

  if (mpHeader)
    success &= compileChildReport(mpHeader, listOfContainer);

  generateObjectsFromName(&listOfContainer, mBodyObjectList, mpBody,
                          mpReportDef->getBodyAddr());

  if (mpBody)
    success &= compileChildReport(mpBody, listOfContainer);

  generateObjectsFromName(&listOfContainer, mFooterObjectList, mpFooter,
                          mpReportDef->getFooterAddr());

  if (mpFooter)
    success &= compileChildReport(mpFooter, listOfContainer);

  mState = Compiled;

  return success;
}

std::ostream * CReport::open(const CCopasiDataModel * pDataModel,
                             std::ostream * pOstream)
{
  mpDataModel = pDataModel;
  assert(mpDataModel != NULL);

  // If an ostream is given and it is the currently assigned one
  // we do nothing.
  if (pOstream != NULL &&
      pOstream == mpOstream)
    return mpOstream;

  if (mStreamOwner)
    pdelete(mpOstream);

  mpOstream = pOstream;

  if (pOstream)
    {
      mStreamOwner = false;
    }
  else if (mTarget != "" && mpReportDef != NULL)
    {
      if (CDirEntry::isRelativePath(mTarget) &&
          !CDirEntry::makePathAbsolute(mTarget, mpDataModel->getFileName()))
        mTarget = CDirEntry::fileName(mTarget);

      mpOstream = new std::ofstream;
      mStreamOwner = true;

      if (mAppend)
        {
          ((std::ofstream *) mpOstream)->
          open(utf8ToLocale(mTarget).c_str(), std::ios_base::out | std::ios_base::app);
        }
      else
        {
          ((std::ofstream *) mpOstream)->
          open(utf8ToLocale(mTarget).c_str(), std::ios_base::out);
        }

      if (!((std::ofstream *) mpOstream)->is_open())
        {
          CCopasiMessage(CCopasiMessage::ERROR, MCDirEntry + 3, mTarget.c_str());
          pdelete(mpOstream);
        }

      if (mpOstream)
        {
          mpOstream->precision(mpReportDef->getPrecision());
        }
    }

  return mpOstream;
}

std::ostream * CReport::getStream() const {return mpOstream;}

// make to support parallel tasks
void CReport::generateObjectsFromName(const std::vector< CCopasiContainer * > * pListOfContainer,
                                      std::vector<CCopasiObject*> & objectList,
                                      CReport *& pReport,
                                      const std::vector<CRegisteredObjectName>* nameVector)
{
  objectList.clear();

  unsigned C_INT32 i;
  CCopasiObject* pSelected;
  CReportDefinition * pReportDefinition;

  for (i = 0; i < nameVector->size(); i++)
    {
      pSelected = const_cast<CCopasiObject *>(mpDataModel->ObjectFromName(*pListOfContainer, (*nameVector)[i]));

      if (!pSelected)
        {
          CCopasiMessage(CCopasiMessage::WARNING, MCCopasiTask + 6, (*nameVector)[i].c_str());
          continue;
        }

      if (!i && (pReportDefinition = dynamic_cast< CReportDefinition * >(pSelected)) != NULL)
        {
          pReport = new CReport();
          pReport->setReportDefinition(pReportDefinition);

          return;
        }

      COutputInterface::mObjects.insert(pSelected);
      objectList.push_back(pSelected);
    }
}

bool CReport::compileChildReport(CReport * pReport, std::vector< CCopasiContainer * > listOfContainer)
{
  pReport->open(mpDataModel, mpOstream);
  bool success = pReport->compile(listOfContainer, mpDataModel);

  const std::set< const CCopasiObject * > & Objects = pReport->COutputInterface::getObjects();
  std::set< const CCopasiObject * >::const_iterator it = Objects.begin();
  std::set< const CCopasiObject * >::const_iterator end = Objects.end();

  for (; it != end; ++it)
    COutputInterface::mObjects.insert(*it);

  return success;
}
