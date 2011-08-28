/* Begin CVS Header
   $Source: /fs/turing/cvs/copasi_dev/copasi/parameterFitting/CFitMethod.cpp,v $
   $Revision: 1.5 $
   $Name: Build-33 $
   $Author: shoops $
   $Date: 2006/04/27 01:30:29 $
   End CVS Header */

// Copyright � 2005 by Pedro Mendes, Virginia Tech Intellectual
// Properties, Inc. and EML Research, gGmbH.
// All rights reserved.

#include "copasi.h"

#include "CFitTask.h"
#include "CFitMethod.h"
#include "CFitProblem.h"

#include "utilities/CCopasiMessage.h"

const std::string CFitMethod::TypeName[] =
  {
    ""
  };

COptMethod * CFitMethod::createMethod(CCopasiMethod::SubType subType)
{
  COptMethod * pMethod = NULL;

  /*
  switch (subType)
    {
    case GeneticAlgorithmSR:
      pMethod = new COptMethodGASR();
     break;

    default:
  */
  pMethod = COptMethod::createMethod(subType);
  /*
      break;
    }
  */

  return pMethod;
}

// Default constructor
CFitMethod::CFitMethod():
    COptMethod(CCopasiTask::parameterFitting, CCopasiMethod::unset),
    mpFitProblem(NULL),
    mpFitTask(NULL)
{CONSTRUCTOR_TRACE;}

CFitMethod::CFitMethod(CCopasiMethod::SubType subType,
                       const CCopasiContainer * pParent):
    COptMethod(CCopasiTask::parameterFitting, subType, pParent),
    mpFitProblem(NULL),
    mpFitTask(NULL)
{CONSTRUCTOR_TRACE;}

CFitMethod::CFitMethod(const CFitMethod & src,
                       const CCopasiContainer * pParent):
    COptMethod(src, pParent),
    mpFitProblem(src.mpFitProblem),
    mpFitTask(src.mpFitTask)
{CONSTRUCTOR_TRACE;}

CFitMethod::~CFitMethod()
{}

bool CFitMethod::initialize()
{
  if (!COptMethod::initialize()) return false;

  mpFitTask = dynamic_cast<CFitTask *>(getObjectParent());
  if (!mpFitTask) return false;

  return true;
}

bool CFitMethod::isValidProblem(const CCopasiProblem * pProblem)
{
  if (!COptMethod::isValidProblem(pProblem)) return false;

  const CFitProblem * pTP = dynamic_cast<const CFitProblem *>(pProblem);
  if (!pTP)
    {
      CCopasiMessage(CCopasiMessage::EXCEPTION, "Problem is not a parameter estimation problem.");
      return false;
    }

  return true;
}
