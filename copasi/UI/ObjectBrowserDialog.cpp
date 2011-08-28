// Begin CVS Header
//   $Source: /fs/turing/cvs/copasi_dev/copasi/UI/ObjectBrowserDialog.cpp,v $
//   $Revision: 1.11 $
//   $Name: Build-33 $
//   $Author: shoops $
//   $Date: 2008/12/18 19:57:54 $
// End CVS Header

// Copyright (C) 2008 by Pedro Mendes, Virginia Tech Intellectual
// Properties, Inc., EML Research, gGmbH, University of Heidelberg,
// and The University of Manchester.
// All rights reserved.

// Copyright (C) 2001 - 2007 by Pedro Mendes, Virginia Tech Intellectual
// Properties, Inc. and EML Research, gGmbH.
// All rights reserved.

#include "ObjectBrowserDialog.h"
#include "qpushbutton.h"
#include "qlayout.h"
#include "q3frame.h"
//Added by qt3to4:
#include <Q3GridLayout>
#include <QCloseEvent>

#include "copasi.h"

#include "copasiui3window.h"
#include "report/CCopasiObject.h"
#include "ObjectBrowserWidget.h"
//#include "qtUtilities.h"

/*
 *  Constructs a ObjectBrowserDialog which is a child of 'parent', with the
 *  name 'name' and widget flags set to 'f'.
 */
ObjectBrowserDialog::ObjectBrowserDialog(QWidget* parent, const char* name, bool modal, int state)
    : QDialog(parent, name, modal),
    clearButton(NULL),
    toggleViewButton(NULL),
    cancelButton(NULL),
    okButton(NULL),
    Line1(NULL),
    spacer(NULL),
    objectBrowserWidget(NULL),
    ObjectBrowserDialogLayout(NULL) //,
    //    mparent(NULL)
{
  this->setWindowFlags(windowFlags() | Qt::WDestructiveClose);
  ObjectBrowserDialogLayout = new Q3GridLayout(this, 3, 4, 11, 6, "ObjectBrowserDialogLayout");
  ObjectBrowserDialogLayout->setAutoAdd(false);

  //  objectBrowserWidget = new ObjectBrowserWidget(this, NULL, 0, 1);
  objectBrowserWidget = new ObjectBrowserWidget(this);
  ObjectBrowserDialogLayout->addMultiCellWidget(objectBrowserWidget, 0, 0, 0, 3);

  Line1 = new Q3Frame(this, "Line1");
  Line1->setFrameShape(Q3Frame::HLine);
  Line1->setFrameShadow(Q3Frame::Sunken);
  Line1->setFrameShape(Q3Frame::HLine);
  ObjectBrowserDialogLayout->addMultiCellWidget(Line1, 1, 1, 0, 3);

  //  spacer = new QSpacerItem(131, 31, QSizePolicy::Expanding, QSizePolicy::Minimum);
  //  ObjectBrowserDialogLayout->addItem(spacer, 2, 1);

  if (state == 0)
    {
      clearButton = new QPushButton(this, "clearButton");
      clearButton->setText("Clear");

      toggleViewButton = new QPushButton(this, "toggleViewButton");
      toggleViewButton->setText("Selected Items");

      cancelButton = new QPushButton(this, "cancelButton");
      cancelButton->setText("Cancel");

      okButton = new QPushButton(this, "okButton");
      okButton->setText("OK");
      okButton->setDefault(true);
      okButton->setAutoDefault(true);

      connect(clearButton, SIGNAL(clicked()), objectBrowserWidget, SLOT(clearClicked()));
      connect(toggleViewButton, SIGNAL(clicked()), this, SLOT(toggleViewClicked()));
      connect(cancelButton, SIGNAL(clicked()), this, SLOT(cancelClicked()));
      connect(okButton, SIGNAL(clicked()), this, SLOT(okClicked()));

      ObjectBrowserDialogLayout->addWidget(clearButton, 2, 0);
      ObjectBrowserDialogLayout->addWidget(toggleViewButton, 2, 1);
      ObjectBrowserDialogLayout->addWidget(cancelButton, 2, 2);
      ObjectBrowserDialogLayout->addWidget(okButton, 2, 3);

      setTabOrder(okButton, clearButton);
      setTabOrder(clearButton, toggleViewButton);
      setTabOrder(toggleViewButton, cancelButton);
      setTabOrder(cancelButton, objectBrowserWidget);
    }
  else
    {
      okButton = new QPushButton(this, "okButton");
      okButton->setText("Close");
      okButton->setDefault(true);
      okButton->setAutoDefault(true);
      ObjectBrowserDialogLayout->addWidget(okButton, 2, 3);

      connect(okButton, SIGNAL(clicked()), this, SLOT(okClicked()));

      setTabOrder(okButton, objectBrowserWidget);
    }

  if (!name)
    setName("ObjectBrowser");
  resize(420, 460);
  setCaption(trUtf8("Object Browser"));
}

ObjectBrowserDialog::~ObjectBrowserDialog()
{
  /*  CopasiUI3Window* uiWindowParent = dynamic_cast<CopasiUI3Window*>(this->parent());
    if (uiWindowParent)
      {
        uiWindowParent->slotObjectBrowserDialogWasClosed();
      }*/
  cleanup();
}

void ObjectBrowserDialog::closeEvent(QCloseEvent* e)
{
  QDialog::closeEvent(e);

  //the code to notify the parent window explicitly is now only called if the
  //dialog is closed by the user, not if it is deleted.
  CopasiUI3Window* uiWindowParent = dynamic_cast<CopasiUI3Window*>(this->parent());
  if (uiWindowParent)
    {
      uiWindowParent->slotObjectBrowserDialogWasClosed();
    }
}

void ObjectBrowserDialog::cleanup()
{
  return;
}

void ObjectBrowserDialog::setOutputVector(std::vector< const CCopasiObject * > * pObjectVector)
{
  if (objectBrowserWidget) objectBrowserWidget->setOutputVector(pObjectVector);
}

void ObjectBrowserDialog::cancelClicked()
{
  close();
  //  QDialog::done(QDialog::Accepted);
}

void ObjectBrowserDialog::okClicked()
{
  objectBrowserWidget->commitClicked();
  done(QDialog::Accepted);
  //  QDialog::done(QDialog::Accepted);
}

void ObjectBrowserDialog::toggleViewClicked()
{
  if (toggleViewButton->text() == "Selected Items")
    toggleViewButton->setText("Tree View");
  else
    toggleViewButton->setText("Selected Items");
  objectBrowserWidget->toggleViewClicked();
}

void ObjectBrowserDialog::selectObjects(std::vector< const CCopasiObject * > * pObjectVector)
{
  if (objectBrowserWidget) objectBrowserWidget->selectObjects(pObjectVector);
}
