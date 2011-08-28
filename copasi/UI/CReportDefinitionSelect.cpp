// Begin CVS Header
//   $Source: /fs/turing/cvs/copasi_dev/copasi/UI/CReportDefinitionSelect.cpp,v $
//   $Revision: 1.52.2.2 $
//   $Name: Build-33 $
//   $Author: shoops $
//   $Date: 2010/09/29 15:51:30 $
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

/********************************************************
Author: Liang Xu
Version : 1.xx  <first>
Description:
Date: 08/15
Comment : CReportDefinitionSelect to select the report definition for one task
Contact: Please contact lixu1@vt.edu.
 *********************************************************/

#include <qvariant.h>
#include <qpushbutton.h>
#include <q3frame.h>
#include <qlabel.h>
#include <qcombobox.h>
#include <qlineedit.h>
#include <qcheckbox.h>
#include <qlayout.h>
#include <qtooltip.h>
#include <q3whatsthis.h>
#include <q3filedialog.h>
//Added by qt3to4:
#include <Q3GridLayout>

#include "copasi.h"
#include "qtUtilities.h"
#include "CQMessageBox.h"
#include "CReportDefinitionSelect.h"
#include "listviews.h"
#include "DataModelGUI.h"
#include "CopasiDataModel/CCopasiDataModel.h"
#include "report/CCopasiRootContainer.h"
#include "utilities/CCopasiException.h"
#include "report/CReportDefinitionVector.h"
#include "report/CReport.h"

/*
 *  Constructs a CReportDefinitionSelect as a child of 'parent', with the
 *  name 'name' and widget flags set to 'f'.
 */

CReportDefinitionSelect::CReportDefinitionSelect(QWidget* parent, const char* name, Qt::WFlags fl)
    : QDialog(parent, name, fl),
    pListView((ListViews*)parent),
    mpReport(NULL),
    bShow(true)
{
  if (!name)
    setName("CReportDefinitionSelect");

  CReportDefinitionSelectLayout = new Q3GridLayout(this, 1, 1, 11, 6, "CReportDefinitionSelectLayout");

  frame5 = new Q3Frame(this, "frame5");
  frame5->setFrameShape(Q3Frame::Box);
  frame5->setFrameShadow(Q3Frame::Sunken);
  frame5Layout = new Q3GridLayout(frame5, 1, 1, 11, 6, "frame5Layout");

  reportLabel = new QLabel(frame5, "reportLabel");

  frame5Layout->addWidget(reportLabel, 0, 0);

  targetLabel = new QLabel(frame5, "targetLabel");

  frame5Layout->addWidget(targetLabel, 1, 0);

  appendChecked = new QCheckBox(frame5, "appendChecked");

  frame5Layout->addMultiCellWidget(appendChecked, 2, 2, 1, 2);

  reportDefinitionNameList = new QComboBox(false, frame5, "reportDefinitionNameList");

  frame5Layout->addWidget(reportDefinitionNameList, 0, 1);

  jumpButton = new QPushButton(frame5, "jumpButton");

  frame5Layout->addWidget(jumpButton, 0, 2);

  targetEdit = new QLineEdit(frame5, "targetEdit");
  targetEdit->setFrameShape(QLineEdit::LineEditPanel);
  targetEdit->setFrameShadow(QLineEdit::Sunken);

  frame5Layout->addWidget(targetEdit, 1, 1);

  browseButton = new QPushButton(frame5, "browseButton");

  frame5Layout->addWidget(browseButton, 1, 2);

  CReportDefinitionSelectLayout->addMultiCellWidget(frame5, 0, 0, 0, 1);

  confirmButton = new QPushButton(this, "confirmButton");

  CReportDefinitionSelectLayout->addWidget(confirmButton, 1, 0);

  cancelButton = new QPushButton(this, "cancelButton");

  CReportDefinitionSelectLayout->addWidget(cancelButton, 1, 1);
  languageChange();
  //clearWState(WState_Polished);

  // tab order
  setTabOrder(reportDefinitionNameList, confirmButton);
  setTabOrder(confirmButton, cancelButton);
  setTabOrder(cancelButton, jumpButton);
  setTabOrder(jumpButton, targetEdit);
  setTabOrder(targetEdit, browseButton);
  setTabOrder(browseButton, appendChecked);

  connect(cancelButton, SIGNAL(clicked()), this, SLOT(cancelClicked()));
  connect(confirmButton, SIGNAL(clicked()), this, SLOT(confirmClicked()));
  connect(jumpButton, SIGNAL(clicked()), this, SLOT(jumpToReportDefinitionEdit()));
  connect(browseButton, SIGNAL(clicked()), this, SLOT(jumpToFileBrowser()));
}

/*
 *  Destroys the object and frees any allocated resources
 */
CReportDefinitionSelect::~CReportDefinitionSelect()
{
  cleanup();
  // no need to delete child widgets, Qt does it all for us
}

/*
 *  Sets the strings of the subwidgets using the current
 *  language.
 */
void CReportDefinitionSelect::languageChange()
{
  setCaption(tr("CReportDefinitionSelect"));
  reportLabel->setText(tr("ReportDefinitions"));
  targetLabel->setText(tr("Target"));
  appendChecked->setText(tr("Append"));
  jumpButton->setText(tr("edit"));
  browseButton->setText(tr("browse"));
  confirmButton->setText(tr("Confirm"));
  cancelButton->setText(tr("Cancel"));
}

void CReportDefinitionSelect::loadReportDefinitionVector()
{
  assert(CCopasiRootContainer::getDatamodelList()->size() > 0);
  CReportDefinitionVector* pReportDefinitionVector = (*CCopasiRootContainer::getDatamodelList())[0]->getReportDefinitionList();
  unsigned C_INT32 i;

  for (i = 0; i < pReportDefinitionVector->size(); i++)
    reportDefinitionNameList->
    insertItem(FROM_UTF8((*(pReportDefinitionVector))[i]->getObjectName()));

  // if it is an empty list
  if (reportDefinitionNameList->count() == 0)
    {
      std::string name = "ReportDefinition_0";
      (*CCopasiRootContainer::getDatamodelList())[0]->getReportDefinitionList()->createReportDefinition(name, "");
      reportDefinitionNameList->insertItem(FROM_UTF8(name));
      reportDefinitionNameList->setCurrentItem(1);
      mpReport->setReportDefinition((*(*CCopasiRootContainer::getDatamodelList())[0]->getReportDefinitionList())[0]); //first one report definition
      mpReport->setAppend(appendChecked->isChecked());
      mpReport->setTarget(TO_UTF8(targetEdit->text()));
      pListView->getDataModel()->notify(ListViews::REPORT, ListViews::CHANGE, ""); //notify Table Definition to

      if (CQMessageBox::question(NULL, "No Report Definition Defined",
                                 "No report definition defined, COPASI has already created a new one for you.\n Do you want to switch to the GUI to edit it?",
                                 QMessageBox::Ok | QMessageBox::No, QMessageBox::Ok) == QMessageBox::Ok)
        jumpToReportDefinitionEdit();

      return;
    }

  if (!mpReport->getReportDefinition())
    {
      C_INT32 row;
      row = reportDefinitionNameList->currentItem();
      mpReport->setReportDefinition((*(pReportDefinitionVector))[row]);
      mpReport->setAppend(appendChecked->isChecked());
      mpReport->setTarget(TO_UTF8(targetEdit->text()));
      return;
    }
  else
    {
      C_INT32 i;

      // no use to compare the last one
      for (i = reportDefinitionNameList->count() - 1; i >= 1; i--)
        if (reportDefinitionNameList->text(i) == FROM_UTF8(mpReport->getReportDefinition()->getObjectName()))
          break;

      reportDefinitionNameList->setCurrentItem(i);
      appendChecked->setChecked(mpReport->append());
      targetEdit->setText(FROM_UTF8(mpReport->getTarget()));
    }
}

void CReportDefinitionSelect::cancelClicked()
{
  cleanup();
  QDialog::done(QDialog::Rejected);
  // delete this;
}

void CReportDefinitionSelect::confirmClicked()
{
  if (!mpReport)
    //exception made here
    return;

  assert(CCopasiRootContainer::getDatamodelList()->size() > 0);
  CReportDefinitionVector* pReportDefinitionVector = (*CCopasiRootContainer::getDatamodelList())[0]->getReportDefinitionList();
  C_INT32 row;
  row = reportDefinitionNameList->currentItem();
  mpReport->setReportDefinition((*(pReportDefinitionVector))[row]);
  mpReport->setAppend(appendChecked->isChecked());
  mpReport->setTarget(TO_UTF8(targetEdit->text()));
  cleanup();
  QDialog::done(QDialog::Accepted);
  //  delete this;
}

void CReportDefinitionSelect::cleanup()
{
  mpReport = NULL;
}

void CReportDefinitionSelect::jumpToReportDefinitionEdit()
{
  assert(CCopasiRootContainer::getDatamodelList()->size() > 0);
  CReportDefinitionVector* pReportDefinitionVector = (*CCopasiRootContainer::getDatamodelList())[0]->getReportDefinitionList();
  C_INT32 row;
  row = reportDefinitionNameList->currentItem();
  pListView->switchToOtherWidget(0, (*pReportDefinitionVector)[row]->getKey());
  confirmClicked(); // if shown then close
  bShow = false; // if not shown then close
}

void CReportDefinitionSelect::jumpToFileBrowser()
{
  QString reportFile =
    CopasiFileDialog::getSaveFileName(this, "Save File Dialog", "untitled.txt", "TEXT Files (*.txt)",
                                      "Choose to create a new a file");

  if (!reportFile.isNull())
    {
      targetEdit->setText(reportFile);
    }
}

void CReportDefinitionSelect::setReport(CReport* newReport)
{
  mpReport = newReport;
}

int CReportDefinitionSelect::exec()
{
  if (bShow)
    return QDialog::exec();
  else
    return QDialog::Accepted;
}
