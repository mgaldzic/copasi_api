// Begin CVS Header
//   $Source: /fs/turing/cvs/copasi_dev/copasi/crosssection/CCrossSectionMethod.cpp,v $
//   $Revision: 1.2 $
//   $Name: Build-33 $
//   $Author: ssahle $
//   $Date: 2010/05/24 12:45:49 $
// End CVS Header

// Copyright (C) 2010 by Pedro Mendes, Virginia Tech Intellectual
// Properties, Inc., University of Heidelberg, and The University
// of Manchester.
// All rights reserved.

#include "CCrossSectionMethod.h"

#include "CCrossSectionProblem.h"

CCrossSectionMethod *
CCrossSectionMethod::createCrossSectionMethod(CCopasiMethod::SubType subType)
{
  CCrossSectionMethod * pMethod = NULL;

  switch (subType)
    {
      case unset:
      case crossSectionMethod:
        pMethod = new CCrossSectionMethod(subType);
        break;

      default:
        fatalError();
    }

  return pMethod;
}

/**
 *  Default constructor.
 */
CCrossSectionMethod::CCrossSectionMethod(CCopasiMethod::SubType subType,
    const CCopasiContainer * pParent):
    CCopasiMethod(CCopasiTask::crosssection, subType, pParent)
    //mpProblem(NULL)
{

  CONSTRUCTOR_TRACE;
}

/**
 *  Copy constructor.
 *  @param "const CCrossSectionMethod &" src
 */
CCrossSectionMethod::CCrossSectionMethod(const CCrossSectionMethod & src,
    const CCopasiContainer * pParent):
    CCopasiMethod(src, pParent)
    //mpProblem(src.mpProblem)
{CONSTRUCTOR_TRACE;}

/**
 *  Destructor.
 */
CCrossSectionMethod::~CCrossSectionMethod()
{DESTRUCTOR_TRACE;}


bool CCrossSectionMethod::process(CProcessReport * handler)
{
  return false;
}

bool CCrossSectionMethod::initialize(CCrossSectionProblem* problem)
{
  return false;
}

bool CCrossSectionMethod::isValidProblem(const CCopasiProblem * pProblem)
{
  if (!CCopasiMethod::isValidProblem(pProblem)) return false;

  const CCrossSectionProblem * pP = dynamic_cast<const CCrossSectionProblem *>(pProblem);

  if (!pP)
    {
      //not a correct Problem
      CCopasiMessage(CCopasiMessage::EXCEPTION, "Problem is not a Cross section problem.");
      return false;
    }

  return true;
}
