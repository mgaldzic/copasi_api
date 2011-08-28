// Begin CVS Header
//   $Source: /fs/turing/cvs/copasi_dev/copasi/UI/CQTSSAWidget.cpp,v $
//   $Revision: 1.16 $
//   $Name: Build-33 $
//   $Author: nsimus $
//   $Date: 2010/07/05 13:25:08 $
// End CVS Header

// Copyright (C) 2010 by Pedro Mendes, Virginia Tech Intellectual
// Properties, Inc., University of Heidelberg, and The University
// of Manchester.
// All rights reserved.

// Copyright (C) 2008 by Pedro Mendes, Virginia Tech Intellectual
// Properties, Inc., EML Research, gGmbH, University of Heidelberg,
// and The University of Manchester.
// All rights reserved.

#include "CQTSSAWidget.h"

#include "copasi.h"

//#include <q3table.h>
#include <qcombobox.h>
//#include <q3header.h>
#include <qtabwidget.h>

#include "CQTSSAResultSubWidget.h"
#include "CQTSSAResultWidget.h"
#include "CQTaskBtnWidget.h"
#include "CQTaskHeaderWidget.h"
#include "CProgressBar.h"
#include "CQValidator.h"
#include "CQMessageBox.h"
#include "qtUtilities.h"

#include "tssanalysis/CTSSATask.h"
#include "tssanalysis/CTSSAProblem.h"
#include "model/CModel.h"
#include "report/CKeyFactory.h"
#include "utilities/CCopasiException.h"
#include "tssanalysis/CCSPMethod.h"
#include "tssanalysis/CILDMMethod.h"
#include "tssanalysis/CILDMModifiedMethod.h"
#include "report/CCopasiRootContainer.h"

/*
 *  Constructs a CQTSSAWidget which is a child of 'parent', with the
 *  name 'name'.'
 */
CQTSSAWidget::CQTSSAWidget(QWidget* parent, const char* name)
    : TaskWidget(parent, name)
{
  setupUi(this);

  init();
}

/*
 *  Destroys the object and frees any allocated resources
 */
CQTSSAWidget::~CQTSSAWidget()
{
  destroy();
  // no need to delete child widgets, Qt does it all for us
}

/*
 *  Sets the strings of the subwidgets using the current
 *  language.
 */
void CQTSSAWidget::languageChange()
{
  retranslateUi(this);
}

CTSSAMethod* pTSSMethod;

CILDMMethod *pILDM_Method;
CILDMModifiedMethod *pILDMModiMethod;

CQTSSAResultSubWidget* pTSSResultSubWidget;
CTSSATask * pCTSSATask;

class mpTSSResultSubWidget;
class QTabWidget;

void CQTSSAWidget::init()
{
  mpTSSAProblem = NULL;

  mpHeaderWidget->setTaskName("Time Scale Separation Analysis");

  vboxLayout->insertWidget(0, mpHeaderWidget);  // header
  vboxLayout->insertSpacing(1, 14);       // space between header and body
  vboxLayout->addWidget(mpBtnWidget);     // 'footer'

  addMethodSelectionBox(CTSSATask::ValidMethods, 0);
  addMethodParameterTable(1);

  mpValidatorDuration = new CQValidatorDouble(mpEditDuration);
  mpEditDuration->setValidator(mpValidatorDuration);

  mpValidatorIntervalSize = new CQValidatorDouble(mpEditIntervalSize);
  mpValidatorIntervalSize->setRange(0, DBL_MAX);
  mpEditIntervalSize->setValidator(mpValidatorIntervalSize);
}

void CQTSSAWidget::destroy()
{
  pdelete(mpTSSAProblem);
}

void CQTSSAWidget::slotDuration()
{
  try
    {
      mpTSSAProblem->setDuration(mpEditDuration->text().toDouble());
    }
  catch (...)
    {
      CQMessageBox::information(this, QString("Information"),
                                FROM_UTF8(CCopasiMessage::getAllMessageText()),
                                QMessageBox::Ok, QMessageBox::Ok);
    }

  mpEditIntervalSize->setText(QString::number(mpTSSAProblem->getStepSize()));
  mpValidatorIntervalSize->revalidate();
  mpEditIntervals->setText(QString::number(mpTSSAProblem->getStepNumber()));

  checkTimeSeries();
}

void CQTSSAWidget::slotIntervalSize()
{
  try
    {
      mpTSSAProblem->setStepSize(mpEditIntervalSize->text().toDouble());
    }

  catch (...)
    {
      CQMessageBox::information(this, QString("Information"),
                                FROM_UTF8(CCopasiMessage::getAllMessageText()),
                                QMessageBox::Ok, QMessageBox::Ok);
    }

  mpEditIntervalSize->setText(QString::number(mpTSSAProblem->getStepSize()));
  mpValidatorIntervalSize->revalidate();
  mpEditIntervals->setText(QString::number(mpTSSAProblem->getStepNumber()));

  checkTimeSeries();
}

void CQTSSAWidget::slotIntervals()
{
  try
    {
      mpTSSAProblem->setStepNumber(mpEditIntervals->text().toULong());
    }
  catch (...)
    {
      CQMessageBox::information(this, QString("Information"),
                                FROM_UTF8(CCopasiMessage::getAllMessageText()),
                                QMessageBox::Ok, QMessageBox::Ok);
    }

  mpEditIntervalSize->setText(QString::number(mpTSSAProblem->getStepSize()));
  mpValidatorIntervalSize->revalidate();

  checkTimeSeries();
}

bool CQTSSAWidget::saveTask()
{
  CTSSATask * pTask =
    dynamic_cast< CTSSATask * >(mpTask);

  if (!pTask) return false;

  saveCommon();
  saveMethod();

  CTSSAProblem* tssaproblem =
    dynamic_cast<CTSSAProblem *>(pTask->getProblem());
  assert(tssaproblem);

  //numbers
  if (tssaproblem->getStepSize() != mpEditIntervalSize->text().toDouble())
    {
      tssaproblem->setStepSize(mpEditIntervalSize->text().toDouble());
      mChanged = true;
    }
  else if (tssaproblem->getStepNumber() != mpEditIntervals->text().toULong())
    {
      tssaproblem->setStepNumber(mpEditIntervals->text().toLong());
      mChanged = true;
    }

  if (tssaproblem->getDuration() != mpEditDuration->text().toDouble())
    {
      tssaproblem->setDuration(mpEditDuration->text().toDouble());
      mChanged = true;
    }

  if (tssaproblem->timeSeriesRequested() != mpCheckSave->isChecked())
    {
      tssaproblem->setTimeSeriesRequested(mpCheckSave->isChecked());
      mChanged = true;
    }

  mpValidatorDuration->saved();
  mpValidatorIntervalSize->saved();

  return true;
}

bool CQTSSAWidget::loadTask()
{
  CTSSATask * pTask =
    dynamic_cast< CTSSATask * >(mpTask);

  if (!pTask) return false;

  loadCommon();
  loadMethod();

  CTSSAProblem* tssaproblem =
    dynamic_cast<CTSSAProblem *>(pTask->getProblem());
  assert(tssaproblem);

  pdelete(mpTSSAProblem);
  mpTSSAProblem = new CTSSAProblem(*tssaproblem);

  //numbers
  mpEditIntervalSize->setText(QString::number(tssaproblem->getStepSize()));
  mpEditIntervals->setText(QString::number(tssaproblem->getStepNumber()));
  mpEditDuration->setText(QString::number(tssaproblem->getDuration()));

  //store time series checkbox
  mpCheckSave->setChecked(tssaproblem->timeSeriesRequested());
  checkTimeSeries();

  mpValidatorDuration->saved();
  mpValidatorIntervalSize->saved();

  return true;
}

CCopasiMethod * CQTSSAWidget::createMethod(const CCopasiMethod::SubType & type)
{
  return CTSSAMethod::createTSSAMethod(type);
}

bool CQTSSAWidget::runTask()
{
  assert(CCopasiRootContainer::getDatamodelList()->size() > 0);
  pCTSSATask =
    dynamic_cast<CTSSATask *>((*(*CCopasiRootContainer::getDatamodelList())[0]->getTaskList())["Time Scale Separation Analysis"]);

  if (!pCTSSATask) return false;

  pTSSMethod = dynamic_cast<CTSSAMethod*>(pCTSSATask->getMethod());

  if (!pTSSMethod)
    pTSSMethod->emptyVectors();

  checkTimeSeries();

  if (!commonBeforeRunTask()) return false;

  bool success = true;

  if (!commonRunTask()) success = false;

  return success;
}

bool CQTSSAWidget::taskFinishedEvent()
{
  bool success = true;
  // We need to load the result here as this is the only place where
  // we know that it is correct.
  CQTSSAResultWidget * pResult =
    dynamic_cast< CQTSSAResultWidget * >(mpListView->findWidgetFromId(271));

  if (pResult == NULL)
    {
      return false;
    }

  success &= pResult->loadFromBackend();

  pTSSResultSubWidget = pResult->getSubWidget();

  if (!pTSSResultSubWidget)
    return false;

  pTSSResultSubWidget->discardOldResults();



  if (success)
    {

      pTSSResultSubWidget->displayResult();

      mpListView->switchToOtherWidget(271, ""); //change to the results window
    }

  return success;
}

void CQTSSAWidget::checkTimeSeries()
{
  assert(CCopasiRootContainer::getDatamodelList()->size() > 0);

  if (mpEditIntervals->text().toLong() *(*CCopasiRootContainer::getDatamodelList())[0]->getModel()->getStateTemplate().getNumVariable() > TSSAMAX)
    {
      mpCheckSave->setChecked(false);
      mpCheckSave->setEnabled(false);
    }
  else
    {
      mpCheckSave->setEnabled(true);
    }
}
