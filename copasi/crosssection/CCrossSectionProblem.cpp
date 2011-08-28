// Begin CVS Header
//   $Source: /fs/turing/cvs/copasi_dev/copasi/crosssection/CCrossSectionProblem.cpp,v $
//   $Revision: 1.3 $
//   $Name: Build-33 $
//   $Author: pwilly $
//   $Date: 2010/05/26 18:51:05 $
// End CVS Header

// Copyright (C) 2010 by Pedro Mendes, Virginia Tech Intellectual
// Properties, Inc., University of Heidelberg, and The University
// of Manchester.
// All rights reserved.

#include "CCrossSectionProblem.h"

CCrossSectionProblem::CCrossSectionProblem(const CCopasiContainer * pParent):
    CCopasiProblem(CCopasiTask::crosssection, pParent),
    mpFlagLimitCrossings(NULL),
    mpCrossingsLimit(NULL),
    mpFlagLimitTime(NULL),
    mpTimeLimit(NULL),
    mpOutputStartTime(NULL)
{
  addParameter("LimitCrossings", CCopasiParameter::BOOL, false);
  mpFlagLimitCrossings = getValue("LimitCrossings").pBOOL;

  addParameter("NumCrossingsLimit", CCopasiParameter::UINT, (unsigned C_INT32)0);
  mpCrossingsLimit = getValue("NumCrossingsLimit").pUINT;

  addParameter("LimitTime", CCopasiParameter::BOOL, true);
  mpFlagLimitTime = getValue("LimitTime").pBOOL;

  addParameter("TimeLimit", CCopasiParameter::DOUBLE, 100.0);
  mpTimeLimit = getValue("TimeLimit").pDOUBLE;

  addParameter("OutputStartTime", CCopasiParameter::DOUBLE, 0.0);
  mpOutputStartTime = getValue("OutputStartTime").pDOUBLE;

  mpTriggerExpression =
    assertParameter("TriggerExpression", CCopasiParameter::EXPRESSION, std::string(""))->getValue().pEXPRESSION;


  initObjects();
  CONSTRUCTOR_TRACE;
}

CCrossSectionProblem::CCrossSectionProblem(const CCrossSectionProblem & src,
    const CCopasiContainer * pParent):
    CCopasiProblem(src, pParent),
    mpFlagLimitCrossings(NULL),
    mpCrossingsLimit(NULL),
    mpFlagLimitTime(NULL),
    mpTimeLimit(NULL),
    mpOutputStartTime(NULL)
{
  mpFlagLimitCrossings = getValue("LimitCrossings").pBOOL;
  mpCrossingsLimit = getValue("NumCrossingsLimit").pUINT;
  mpFlagLimitTime = getValue("LimitTime").pBOOL;
  mpTimeLimit = getValue("TimeLimit").pDOUBLE;
  mpOutputStartTime = getValue("OutputStartTime").pDOUBLE;
  mpTriggerExpression =
    assertParameter("TriggerExpression", CCopasiParameter::EXPRESSION, std::string(""))->getValue().pEXPRESSION;

  initObjects();
  CONSTRUCTOR_TRACE;
}

void CCrossSectionProblem::initObjects()
{
  //here we should create things like object references to results of the calculation
}

/**
 *  Destructor.
 */
CCrossSectionProblem::~CCrossSectionProblem()
{DESTRUCTOR_TRACE;}

std::ostream &operator<<(std::ostream &os, const CCrossSectionProblem & o)
{
  os << "Cross Section Problem description: Not implemented yet." << std::endl;
  const CCopasiDataModel* pDataModel = o.getObjectDataModel();
  assert(pDataModel != NULL);

  return os;
}


void CCrossSectionProblem::print(std::ostream * ostream) const
{*ostream << *this;}

bool CCrossSectionProblem::getFlagLimitCrossings() const
{ return *mpFlagLimitCrossings; }

const unsigned C_INT32 & CCrossSectionProblem::getCrossingsLimit() const
{ return *mpCrossingsLimit; }

bool CCrossSectionProblem::getFlagLimitTime() const
{ return *mpFlagLimitTime; }

const C_FLOAT64 & CCrossSectionProblem::getTimeLimit() const
{ return *mpTimeLimit; }

const C_FLOAT64 & CCrossSectionProblem::getOutputStartTime() const
{ return *mpOutputStartTime; }

void CCrossSectionProblem::setFlagLimitCrossings(bool flagLimitCrossing)
{ *mpFlagLimitCrossings = flagLimitCrossing; }

void CCrossSectionProblem::setCrossingsLimit(const unsigned C_INT32 &crossingLimit)
{ *mpCrossingsLimit = crossingLimit; }

void CCrossSectionProblem::setFlagLimitTime(bool flagLimitTime)
{ *mpFlagLimitTime = flagLimitTime; }

void CCrossSectionProblem::setTimeLimit(const C_FLOAT64 &timeLimit)
{ *mpTimeLimit = timeLimit; }

void CCrossSectionProblem::setOutputStartTime(const C_FLOAT64 &outputStartTime)
{ *mpOutputStartTime = outputStartTime; }

