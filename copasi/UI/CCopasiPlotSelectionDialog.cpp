// Begin CVS Header
//   $Source: /fs/turing/cvs/copasi_dev/copasi/UI/CCopasiPlotSelectionDialog.cpp,v $
//   $Revision: 1.12 $
//   $Name: Build-33 $
//   $Author: shoops $
//   $Date: 2009/04/21 16:20:31 $
// End CVS Header

// Copyright (C) 2008 by Pedro Mendes, Virginia Tech Intellectual
// Properties, Inc., EML Research, gGmbH, University of Heidelberg,
// and The University of Manchester.
// All rights reserved.

// Copyright (C) 2001 - 2007 by Pedro Mendes, Virginia Tech Intellectual
// Properties, Inc. and EML Research, gGmbH.
// All rights reserved.

#include <qlayout.h>
#include <qwidget.h>
#include <qpushbutton.h>
#include <qcheckbox.h>
#include <qlabel.h>
#include <qsplitter.h>
#include <q3hbox.h>
#include <q3vbox.h>
//Added by qt3to4:
#include <Q3HBoxLayout>
#include <Q3VBoxLayout>

#include "copasi.h"

#include "CCopasiPlotSelectionDialog.h"
#include "CCopasiSelectionWidget.h"
#include "qtUtilities.h"
#include "CQMessageBox.h"
#include "model/CModel.h"
#include "report/CCopasiObject.h"

CCopasiPlotSelectionDialog::CCopasiPlotSelectionDialog(QWidget* parent, const char* name, bool modal, Qt::WFlags f):
    QDialog(parent, name, modal, f)
    , mpOKButton(NULL)
    , mpCancelButton(NULL)
    , mpExpertCheckBox(NULL)
    , mpXAxisSelectionWidget(NULL)
    , mpYAxisSelectionWidget(NULL)
    , mpSplitter(NULL)
    , mpButtonBox(NULL)
    , mpMainLayout(NULL)
    , mpXAxisLabel(NULL)
    , mpYAxisLabel(NULL)
    , mpXAxisSelectionBox(NULL)
    , mpYAxisSelectionBox(NULL)
    , mpXAxisOutputVector(NULL)
    , mpYAxisOutputVector(NULL)
{
  mpMainLayout = new Q3VBoxLayout(this);

  mpSplitter = new QSplitter(this);
  mpSplitter->setOrientation(Qt::Horizontal);
  mpMainLayout->addWidget(mpSplitter);

  mpButtonBox = new Q3HBoxLayout(mpMainLayout);

  mpOKButton = new QPushButton(this);
  mpOKButton->setText("OK");
  mpOKButton->setDefault(true);
  mpButtonBox->add(mpOKButton);

  mpCancelButton = new QPushButton(this);
  mpCancelButton->setText("Cancel");
  mpButtonBox->add(mpCancelButton);

  mpExpertCheckBox = new QCheckBox(this);
  mpExpertCheckBox->setText("Expert Mode");
  mpExpertCheckBox->setChecked(false);
  mpButtonBox->add(mpExpertCheckBox);

  mpXAxisSelectionBox = new Q3VBox(mpSplitter);
  mpXAxisSelectionBox->layout()->setMargin(5);
  mpXAxisSelectionBox->layout()->setAutoAdd(false);

  mpYAxisSelectionBox = new Q3VBox(mpSplitter);
  mpYAxisSelectionBox->layout()->setMargin(5);
  mpYAxisSelectionBox->layout()->setAutoAdd(false);

  mpXAxisLabel = new QLabel("X-Axis:", mpXAxisSelectionBox);
  mpXAxisLabel->setSizePolicy(QSizePolicy::MinimumExpanding, QSizePolicy::Fixed);
  mpXAxisSelectionBox->layout()->add(mpXAxisLabel);

  mpYAxisLabel = new QLabel("Y-Axis:", mpYAxisSelectionBox);
  mpYAxisLabel->setSizePolicy(QSizePolicy::MinimumExpanding, QSizePolicy::Fixed);
  mpYAxisSelectionBox->layout()->add(mpYAxisLabel);

  mpXAxisSelectionWidget = new CCopasiSelectionWidget(mpXAxisSelectionBox);
  mpXAxisSelectionWidget->setSingleSelection(true);
  mpXAxisSelectionWidget->setOutputVector(mpXAxisOutputVector);
  mpXAxisSelectionBox->layout()->add(mpXAxisSelectionWidget);

  mpYAxisSelectionWidget = new CCopasiSelectionWidget(mpYAxisSelectionBox);
  mpYAxisSelectionWidget->setSingleSelection(false);
  mpYAxisSelectionWidget->setOutputVector(mpYAxisOutputVector);
  mpYAxisSelectionBox->layout()->add(mpYAxisSelectionWidget);

  connect(mpOKButton, SIGNAL(clicked()), this, SLOT(slotOKButtonClicked()));
  connect(mpCancelButton, SIGNAL(clicked()), this, SLOT(slotCancelButtonClicked()));
  connect(mpExpertCheckBox, SIGNAL(toggled(bool)), this, SLOT(slotExpertCheckBoxToggled(bool)));

  setTabOrder();
}

CCopasiPlotSelectionDialog::~CCopasiPlotSelectionDialog()
{}

void CCopasiPlotSelectionDialog::setTabOrder()
{
  QWidget::setTabOrder(this->mpOKButton, this->mpCancelButton);
  QWidget::setTabOrder(this->mpCancelButton, this->mpExpertCheckBox);
  QWidget::setTabOrder(this->mpExpertCheckBox, this->mpXAxisSelectionWidget);
  QWidget::setTabOrder(this->mpXAxisSelectionWidget, this->mpYAxisSelectionWidget);
  this->mpXAxisSelectionWidget->setFocus();
}

void CCopasiPlotSelectionDialog::slotOKButtonClicked()
{
  // fill the selection vectors
  this->mpXAxisSelectionWidget->commit();
  this->mpYAxisSelectionWidget->commit();
  std::string message = "";
  bool showWarning = false;

  if (this->mpXAxisOutputVector->empty())
    {
      message += "X Axis";
      showWarning = true;
    }

  if (this->mpYAxisOutputVector->empty())
    {
      if (showWarning)
        {
          message += " and ";
        }

      showWarning = true;
      message += "Y Axis";
    }

  if (showWarning)
    {
      message = "You did not select anything for the " + message + "!\nDo you want to proceed anyway?";

      if (CQMessageBox::question(this, "Empty Selection", FROM_UTF8(message), QMessageBox::Yes | QMessageBox::No, QMessageBox::Yes) == QMessageBox::No)
        {
          return;
        }
    }

  QDialog::done(QDialog::Accepted);
}

void CCopasiPlotSelectionDialog::slotCancelButtonClicked()
{
  QDialog::done(QDialog::Rejected);
}

void CCopasiPlotSelectionDialog::slotExpertCheckBoxToggled(bool checked)
{
  this->mpXAxisSelectionWidget->setExpertMode(checked);
  this->mpYAxisSelectionWidget->setExpertMode(checked);
}

void CCopasiPlotSelectionDialog::setOutputVectors(std::vector< const CCopasiObject * > * outputVector1,
    std::vector< const CCopasiObject * > * outputVector2)
{
  this->mpXAxisOutputVector = outputVector1;
  this->mpXAxisSelectionWidget->setOutputVector(this->mpXAxisOutputVector);
  this->mpYAxisOutputVector = outputVector2;
  this->mpYAxisSelectionWidget->setOutputVector(this->mpYAxisOutputVector);
}

void CCopasiPlotSelectionDialog::setModel(CModel* pModel,
    const CCopasiSimpleSelectionTree::ObjectClasses & classes)
{
  this->mpXAxisSelectionWidget->populateTree(pModel, classes);
  this->mpYAxisSelectionWidget->populateTree(pModel, classes);
}
