// Begin CVS Header
//   $Source: /fs/turing/cvs/copasi_dev/copasi/UI/CQFittingWidget.cpp,v $
//   $Revision: 1.19 $
//   $Name: Build-33 $
//   $Author: aekamal $
//   $Date: 2010/05/10 16:12:14 $
// End CVS Header

// Copyright (C) 2010 by Pedro Mendes, Virginia Tech Intellectual
// Properties, Inc., University of Heidelberg, and The University
// of Manchester.
// All rights reserved.

// Copyright (C) 2008 by Pedro Mendes, Virginia Tech Intellectual
// Properties, Inc., EML Research, gGmbH, University of Heidelberg,
// and The University of Manchester.
// All rights reserved.

#include "CQFittingWidget.h"

#include <qlabel.h>
#include <qtoolbutton.h>
//Added by qt3to4:
//#include <Q3HBoxLayout>

#include "UI/CQTaskBtnWidget.h"
#include "UI/CQTaskHeaderWidget.h"
#include "UI/CQFittingItemWidget.h"
#include "UI/CProgressBar.h"
#include "UI/CQExperimentData.h"

#include "report/CKeyFactory.h"
#include "parameterFitting/CFitTask.h"
#include "parameterFitting/CFitItem.h"
#include "parameterFitting/CFitMethod.h"
#include "parameterFitting/CFitProblem.h"
#include "parameterFitting/CExperimentSet.h"
#include "parameterFitting/CExperiment.h"
#include "CopasiDataModel/CCopasiDataModel.h"
#include "utilities/CCopasiException.h"
#include "report/CCopasiRootContainer.h"

/*
 *  Constructs a CQFittingWidget which is a child of 'parent', with the
 *  name 'name'.'
 */
CQFittingWidget::CQFittingWidget(QWidget* parent, const char* name)
    : TaskWidget(parent, name)
{
  setupUi(this);

  init();
}

/*
 *  Destroys the object and frees any allocated resources
 */
CQFittingWidget::~CQFittingWidget()
{
  destroy();
  // no need to delete child widgets, Qt does it all for us
}

/*
 *  Sets the strings of the subwidgets using the current
 *  language.
 */
void CQFittingWidget::languageChange()
{
  retranslateUi(this);
}

bool CQFittingWidget::saveTask()
{
  CFitTask * pTask =
    dynamic_cast< CFitTask * >(mpTask);

  if (!pTask) return false;

  saveCommon();
  saveMethod();

  CFitProblem * pProblem =
    dynamic_cast<CFitProblem *>(mpTask->getProblem());

  if (!pProblem) return false;

  std::map<std::string, std::string> ExperimentMap;
  std::map<std::string, std::string> CrossValidationMap;

  // Save experiment set
  CExperimentSet * pExperimentSet =
    dynamic_cast<CExperimentSet *>(pProblem->getGroup("Experiment Set"));

  std::string Name;
  CExperiment * pExperiment;
  unsigned C_INT32 i, imax;

  std::set<CExperiment *> DealtWith;

  for (i = pExperimentSet->getExperimentCount() - 1; i != C_INVALID_INDEX; i--)
    {
      pExperiment =
        dynamic_cast<CExperiment *>(CCopasiRootContainer::getKeyFactory()->get(mExperimentKeyMap[pExperimentSet->getExperiment(i)->CCopasiParameter::getKey()]));

      if (pExperiment && pExperiment == mpExperimentSet->getExperiment(pExperiment->getObjectName()))
        {
          if (!(*pExperimentSet->getExperiment(i) == *pExperiment))
            {
              *pExperimentSet->getExperiment(i) = *pExperiment;
              mChanged = true;
            }

          DealtWith.insert(pExperiment);
        }
      else
        {
          mExperimentKeyMap.erase(pExperimentSet->getExperiment(i)->CCopasiParameter::getKey());
          pExperimentSet->removeExperiment(i);
          mChanged = true;
        }
    }

  for (i = 0, imax = mpExperimentSet->getExperimentCount(); i < imax; i++)
    {
      pExperiment = mpExperimentSet->getExperiment(i);

      if (DealtWith.count(pExperiment)) continue;

      pExperiment = pExperimentSet->addExperiment(*pExperiment);
      mExperimentKeyMap[pExperiment->CCopasiParameter::getKey()] =
        mpExperimentSet->getExperiment(i)->CCopasiParameter::getKey();
      mChanged = true;
    }

  // We need to invert the key map for saving!
  std::map<std::string, std::string>::iterator it = mExperimentKeyMap.begin();
  std::map<std::string, std::string>::iterator end = mExperimentKeyMap.end();

  for (; it != end; ++it)
    ExperimentMap[it->second] = it->first;


  mChanged |= mpParameters->save(&ExperimentMap, &CrossValidationMap);
  mChanged |= mpConstraints->save(&ExperimentMap, &CrossValidationMap);

  assert(CCopasiRootContainer::getDatamodelList()->size() > 0);

  if (mChanged)(*CCopasiRootContainer::getDatamodelList())[0]->changed();

  mChanged = false;
  return true;
}

bool CQFittingWidget::loadTask()
{
  CFitTask * pTask =
    dynamic_cast< CFitTask * >(mpTask);

  if (!pTask) return false;

  loadCommon();
  loadMethod();

  CFitProblem * pProblem =
    dynamic_cast<CFitProblem *>(mpTask->getProblem());

  if (!pProblem) return false;

  pdelete(mpExperimentSet)
  CExperimentSet * pExperimentSet =
    dynamic_cast<CExperimentSet *>(pProblem->getGroup("Experiment Set"));
  mpExperimentSet = new CExperimentSet(*pExperimentSet);

  mExperimentKeyMap.clear();
  unsigned C_INT32 i, imax = mpExperimentSet->getExperimentCount();

  for (i = 0; i < imax; i++)
    mExperimentKeyMap[pExperimentSet->getExperiment(i)->CCopasiParameter::getKey()] =
      mpExperimentSet->getExperiment(i)->CCopasiParameter::getKey();


  mpParameters->load(mpDataModel, pProblem->getGroup("OptimizationItemList"), &mExperimentKeyMap, &mCrossValidationKeyMap);
  mpParameters->setExperimentSet(const_cast<const CExperimentSet *&>(mpExperimentSet));

  mpConstraints->load(mpDataModel, pProblem->getGroup("OptimizationConstraintList"), &mExperimentKeyMap, &mCrossValidationKeyMap);
  mpConstraints->setExperimentSet(const_cast<const CExperimentSet *&>(mpExperimentSet));

  mChanged = false;

  return true;
}

CCopasiMethod * CQFittingWidget::createMethod(const CCopasiMethod::SubType & type)
{return CFitMethod::createMethod(type);}

bool CQFittingWidget::runTask()
{
  CFitTask * pTask =
    dynamic_cast< CFitTask * >(CCopasiRootContainer::getKeyFactory()->get(mKey));

  if (!pTask) return false;

  if (!commonBeforeRunTask()) return false;

  bool success = commonRunTask();

  return success;
}

void CQFittingWidget::slotExperimentData()
{
  CQExperimentData * pDialog = new CQExperimentData(this);
  pDialog->load(mpExperimentSet, mpDataModel);

  connect(pDialog, SIGNAL(experimentChanged()), mpParameters, SLOT(slotExperimentChanged()));
  connect(pDialog, SIGNAL(experimentChanged()), mpConstraints, SLOT(slotExperimentChanged()));

  pDialog->exec();

  pdelete(pDialog);
}

void CQFittingWidget::slotPageChange(QWidget * currentPage)
{
  if (mpTabWidget->tabLabel(currentPage).contains("Parameters", true))
    mpCurrentList = mpParameters;
  else
    mpCurrentList = mpConstraints;
}

void CQFittingWidget::init()
{
  mpHeaderWidget->setTaskName("Parameter Estimation");

  vboxLayout->insertWidget(0, mpHeaderWidget);
  vboxLayout->insertSpacing(1, 14);      // space between header and body
  vboxLayout->addWidget(mpBtnWidget);

  addMethodSelectionBox(CFitTask::ValidMethods, 0);
  addMethodParameterTable(1);

//  mpParameterPageLayout = new Q3HBoxLayout(mpParametersPage, 0, 6, "mpParameterPageLayout");
  mpParameterPageLayout = new QHBoxLayout(mpParametersPage);
  mpParameters = new CQFittingItemWidget(mpParametersPage);
  mpParameters->setItemType(CQFittingItemWidget::FIT_ITEM);
  mpParameterPageLayout->addWidget(mpParameters);
  connect(mpParameters, SIGNAL(numberChanged(int)), this, SLOT(slotParameterNumberChanged(int)));

//  mpConstraintPageLayout = new Q3HBoxLayout(mpConstraintsPage, 0, 6, "mpConstraintsPageLayout");
  mpConstraintPageLayout = new QHBoxLayout(mpConstraintsPage);
  mpConstraints = new CQFittingItemWidget(mpConstraintsPage);
  mpConstraints->setItemType(CQFittingItemWidget::FIT_CONSTRAINT);
  mpConstraintPageLayout->addWidget(mpConstraints);
  connect(mpConstraints, SIGNAL(numberChanged(int)), this, SLOT(slotConstraintNumberChanged(int)));

  mpCurrentList = mpParameters;
  mpExperimentSet = NULL;
  mpCrossValidationSet = NULL;

#ifndef COPASI_CROSSVALIDATION
  mpBtnCrossValidation->hide();
#endif
}

void CQFittingWidget::slotParameterNumberChanged(int number)
{
  QString TabLabel = "Parameters (" + QString::number(number) + ")";
  mpTabWidget->setTabLabel(mpParametersPage, TabLabel);
}

void CQFittingWidget::slotConstraintNumberChanged(int number)
{
  QString TabLabel = "Constraints (" + QString::number(number) + ")";
  mpTabWidget->setTabLabel(mpConstraintsPage, TabLabel);
}

void CQFittingWidget::destroy()
{
  pdelete(mpExperimentSet);

}

void CQFittingWidget::slotCrossValidationData()
{
}
