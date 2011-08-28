// Begin CVS Header
//   $Source: /fs/turing/cvs/copasi_dev/copasi/UI/CQNotes.cpp,v $
//   $Revision: 1.12.2.2 $
//   $Name: Build-33 $
//   $Author: shoops $
//   $Date: 2010/09/29 15:59:09 $
// End CVS Header

// Copyright (C) 2010 by Pedro Mendes, Virginia Tech Intellectual
// Properties, Inc., University of Heidelberg, and The University
// of Manchester.
// All rights reserved.

/*
 * CQNotes.cpp
 *
 *  Created on: Aug 11, 2010
 *      Author: shoops
 */

#include <QWebFrame>
#include <QProcess>
#include <QXmlInputSource>
#include <QXmlSimpleReader>

#include "CQNotes.h"
#include "CQIcons.h"
#include "CQMessageBox.h"
#include "qtUtilities.h"

#include "model/CModelValue.h"
#include "model/CReaction.h"
#include "model/CEvent.h"
#include "function/CFunction.h"
#include "report/CKeyFactory.h"
#include "copasi/report/CCopasiRootContainer.h"
#include "commandline/CConfigurationFile.h"

CQValidatorXML::CQValidatorXML(QPlainTextEdit * parent, const char * name):
    CQValidator< QPlainTextEdit >(parent, &QPlainTextEdit::toPlainText, name)
{}

// virtual
QValidator::State CQValidatorXML::validate(QString & input, int & pos) const
{
  QXmlSimpleReader Validator;
  QXmlInputSource Input;

  // We like to allow free text and therefore wrap the text to create valid XML.
  Input.setData("<ValidateXML>" + input + "</ValidateXML>");

  if (Validator.parse(Input))
    return CQValidator< QPlainTextEdit >::validate(input, pos);

  setColor(Invalid);
  return Intermediate;
}

CQNotes::CQNotes(QWidget* parent, const char* name) :
    CopasiWidget(parent, name),
    mEditMode(false),
    mChanged(false),
    mpValidatorXML(NULL),
    mValidity(QValidator::Acceptable)

{
  setupUi(this);

  mpValidatorXML = new CQValidatorXML(mpEdit);

  mEditMode = false;
  mpEdit->hide();
  mpWebView->show();
  mpBtnToggleEdit->setIcon(CQIcons::getIcon(CQIcons::Edit));

  mpWebView->page()->setLinkDelegationPolicy(QWebPage::DelegateAllLinks);
}

CQNotes::~CQNotes()
{}

// virtual
bool CQNotes::update(ListViews::ObjectType /* objectType */, ListViews::Action action, const std::string & key)
{
  if (key == mKey)
    {
      switch (action)
        {
          case ListViews::CHANGE:
            load();
            break;

          case ListViews::DELETE:
            mpObject = NULL;
            mKey = "";
            break;

          default:
            break;
        }
    }

  return true;
}

// virtual
bool CQNotes::leave()
{
  mpBtnToggleEdit->setFocus();

  mpObject = CCopasiRootContainer::getKeyFactory()->get(mKey);

  if (mpObject != NULL)
    {
      save();
    }
  else
    {
      mKey = "";
      mpDataModel = NULL;
    }

  return true;
}

// virtual
bool CQNotes::enterProtected()
{
  load();

  return true;
}

void CQNotes::slotToggleMode()
{
  mEditMode = !mEditMode;

  if (mEditMode)
    {
      mpWebView->hide();
      mpEdit->show();
      mpBtnToggleEdit->setIcon(CQIcons::getIcon(CQIcons::View));
    }
  else
    {
      save();
      load();

      mpEdit->hide();
      mpWebView->show();
      mpBtnToggleEdit->setIcon(CQIcons::getIcon(CQIcons::Edit));
    }
}

void CQNotes::slotValidateXML()
{
  QString Input = mpEdit->toPlainText();
  int pos = 0;

  mValidity = mpValidatorXML->validate(Input, pos);
}

void CQNotes::load()
{
  if (mpObject != NULL)
    {
      const std::string * pNotes = NULL;

      if (dynamic_cast< CModelEntity * >(mpObject))
        pNotes = &static_cast< CModelEntity * >(mpObject)->getNotes();
      else if (dynamic_cast< CEvent * >(mpObject))
        pNotes = &static_cast< CEvent * >(mpObject)->getNotes();
      else if (dynamic_cast< CReaction * >(mpObject))
        pNotes = &static_cast< CReaction * >(mpObject)->getNotes();
      else if (dynamic_cast< CFunction * >(mpObject))
        pNotes = &static_cast< CFunction * >(mpObject)->getNotes();

      if (pNotes != NULL)
        {
          // The notes are UTF8 encoded however the html does not specify an encoding
          // thus Qt uses locale settings.
          mpWebView->setHtml(FROM_UTF8(*pNotes));
          mpEdit->setPlainText(FROM_UTF8(*pNotes));
          mpValidatorXML->saved();

          mValidity = QValidator::Acceptable;
        }
    }

  mChanged = false;

  return;
}

void CQNotes::save()
{
  if (mpObject != NULL &&
      mValidity == QValidator::Acceptable)
    {
      const std::string * pNotes = NULL;

      if (dynamic_cast< CModelEntity * >(mpObject))
        pNotes = &static_cast< CModelEntity * >(mpObject)->getNotes();
      else if (dynamic_cast< CEvent * >(mpObject))
        pNotes = &static_cast< CEvent * >(mpObject)->getNotes();
      else if (dynamic_cast< CReaction * >(mpObject))
        pNotes = &static_cast< CReaction * >(mpObject)->getNotes();
      else if (dynamic_cast< CFunction * >(mpObject))
        pNotes = &static_cast< CFunction * >(mpObject)->getNotes();

      if (pNotes &&
          mpEdit->toPlainText() != FROM_UTF8(*pNotes))
        {
          if (dynamic_cast< CModelEntity * >(mpObject))
            static_cast< CModelEntity * >(mpObject)->setNotes(TO_UTF8(mpEdit->toPlainText()));
          else if (dynamic_cast< CEvent * >(mpObject))
            static_cast< CEvent * >(mpObject)->setNotes(TO_UTF8(mpEdit->toPlainText()));
          else if (dynamic_cast< CReaction * >(mpObject))
            static_cast< CReaction * >(mpObject)->setNotes(TO_UTF8(mpEdit->toPlainText()));
          else if (dynamic_cast< CFunction * >(mpObject))
            static_cast< CFunction * >(mpObject)->setNotes(TO_UTF8(mpEdit->toPlainText()));

          mChanged = true;
        }
    }

  if (mChanged)
    {
      if (mpDataModel != NULL)
        {
          mpDataModel->changed();
        }

      protectedNotify(ListViews::MODEL, ListViews::CHANGE, mKey);
      mChanged = false;
    }

  return;
}

void CQNotes::slotOpenUrl(const QUrl & url)
{
  QString Commandline = FROM_UTF8(CCopasiRootContainer::getConfiguration()->getWebBrowser());

  if (Commandline == "")
    {
#ifdef Q_WS_MAC
      Commandline = "open %1";
#else
# ifdef Q_WS_WIN
      Commandline = "cmd /c start %1";
# else
      CQMessageBox::critical(this, "Unable to open link",
                             "COPASI requires you to specify an application for opening URLs for links to work.\n\nPlease go to the preferences and set an appropriate application in the format:\n  command [options] %1");

      return;
# endif  // Q_WS_WIN
#endif // Q_WS_MAC

      CCopasiRootContainer::getConfiguration()->setWebBrowser(TO_UTF8(Commandline));
    }

#ifdef Q_WS_WIN

  if (Commandline == "cmd /c start %1")
    {
      if (QProcess::execute(Commandline.arg(url.toString())) != 0)
        {
          CQMessageBox::critical(this, "Unable to open link",
                                 "COPASI requires you to specify an application for opening links. The currently provided command:\n  " +
                                 Commandline + "\nis not working properly.\n\nPlease go to the preferences and set an appropriate application in the format:\n  application [options] %1");
        }

      return;
    }

#endif // Q_WS_WIN

  if (!QProcess::startDetached(Commandline.arg(url.toString())))
    {
      CQMessageBox::critical(this, "Unable to open link",
                             "COPASI requires you to specify an application for opening links. The currently provided command:\n  " +
                             Commandline + "\nis not working properly.\n\nPlease go to the preferences and set an appropriate application in the format:\n  application [options] %1");
    }

  return;
}

