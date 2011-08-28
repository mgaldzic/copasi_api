// Begin CVS Header
//   $Source: /fs/turing/cvs/copasi_dev/copasi/UI/CQProgressDialog.cpp,v $
//   $Revision: 1.8 $
//   $Name: Build-33 $
//   $Author: shoops $
//   $Date: 2010/09/02 14:30:56 $
// End CVS Header

// Copyright (C) 2010 by Pedro Mendes, Virginia Tech Intellectual
// Properties, Inc., University of Heidelberg, and The University
// of Manchester.
// All rights reserved.

// Copyright (C) 2008 by Pedro Mendes, Virginia Tech Intellectual
// Properties, Inc., EML Research, gGmbH, University of Heidelberg,
// and The University of Manchester.
// All rights reserved.

#include "CQProgressDialog.h"

#include <qvariant.h>
#include <qapplication.h>
#include <qtimer.h>

/*
 *  Constructs a CQProgressDialog as a child of 'parent', with the
 *  name 'name' and widget flags set to 'f'.
 *
 *  The dialog will by default be modeless, unless you set 'modal' to
 *  true to construct a modal dialog.
 */
CQProgressDialog::CQProgressDialog(QWidget* parent, const char* name, bool modal, Qt::WindowFlags fl)
    : QDialog(parent, name, modal, fl)
{
  setupUi(this);

  init();
}

/*
 *  Destroys the object and frees any allocated resources
 */
CQProgressDialog::~CQProgressDialog()
{
  // no need to delete child widgets, Qt does it all for us
}

/*
 *  Sets the strings of the subwidgets using the current
 *  language.
 */
void CQProgressDialog::languageChange()
{
  retranslateUi(this);
}

bool CQProgressDialog::insertProgressItem(CQProgressItem * pItem)
{
  if (mItemCount == 0) mpLine->show();

  vboxLayout->removeItem(mpSpacer);
  vboxLayout->addWidget(pItem);
  vboxLayout->addItem(mpSpacer);

  pItem->show();

  mItemCount++;
  return true;
}

bool CQProgressDialog::removeProgressItem(CQProgressItem * pItem)
{
  vboxLayout->remove(pItem);

  mItemCount--;

  if (mItemCount == 0) mpLine->hide();

  return true;
}

void CQProgressDialog::init()
{
  mpBtnPause->setIcon(CQIcons::getIcon(CQIcons::Pause));
  mpBtnContinue->setIcon(CQIcons::getIcon(CQIcons::Continue));
  mpBtnStop->setIcon(CQIcons::getIcon(CQIcons::Stop));

  mpLine->hide();
  //  mpBtnPause->hide();
  //  mpBtnContinue->hide();

  mpBtnContinue->setEnabled(false);

  mItemCount = 0;

  mPause = false;
  mProceed = true;

  mpTimer = new QTimer(this);
  mpTimer->start(1500, true);
  connect(mpTimer, SIGNAL(timeout()), this, SLOT(timerShow()));

  return;
}

void CQProgressDialog::btnContinuePressed()
{
  mPause = false;

  mpBtnPause->setEnabled(true);
  mpBtnContinue->setEnabled(false);
}

void CQProgressDialog::btnPausePressed()
{
  mPause = true;

  mpBtnPause->setEnabled(false);
  mpBtnContinue->setEnabled(true);
}

void CQProgressDialog::btnStopPressed()
{
  mPause = false;
  mProceed = false;

  mpBtnPause->setEnabled(false);
  mpBtnContinue->setEnabled(false);
  mpBtnStop->setEnabled(false);
}

void CQProgressDialog::timerShow()
{this->show();}
