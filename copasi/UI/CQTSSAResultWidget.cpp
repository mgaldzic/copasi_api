// Begin CVS Header
//   $Source: /fs/turing/cvs/copasi_dev/copasi/UI/CQTSSAResultWidget.cpp,v $
//   $Revision: 1.10.2.1 $
//   $Name: Build-33 $
//   $Author: shoops $
//   $Date: 2010/09/29 11:46:54 $
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

//#include <qpushbutton.h>
#include <qlayout.h>
//Added by qt3to4:
#include <Q3GridLayout>
#include "copasi.h"

#include "CQTSSAResultWidget.h"
#include "CQTSSAResultSubWidget.h"
#include "CopasiDataModel/CCopasiDataModel.h"
#include "report/CCopasiRootContainer.h"
#include "utilities/CCopasiVector.h"
#include "tssanalysis/CTSSATask.h"
#include "CTimeSeriesTable.h"

//#include "report/CKeyFactory.h"
#include "qtUtilities.h"
#include "tssanalysis/CTSSAMethod.h"

#include "tssanalysis/CTSSATask.h"

CTSSATask* pTask;
CTSSAMethod* pTSSILDM;

/*
 *  Constructs a CQTSSAResultWidget which is a child of 'parent', with the
 *  name 'name' and widget flags set to 'f'.
 */
CQTSSAResultWidget::CQTSSAResultWidget(QWidget* parent, const char* name, Qt::WFlags fl)
    : CopasiWidget(parent, name, fl)
{
  if (!name)
    setName("CQTSSAResultWidget");

  setCaption(trUtf8("CQTSSAResultWidget"));

  mWidgetLayout = new Q3GridLayout(this, 1, 1, 0, -1, "Layout");

  mCentralWidget = new CQTSSAResultSubWidget(this, "CQTSSAResultSubWidget");
  mWidgetLayout->addWidget(mCentralWidget, 0, 0);

  /*commitChanges = new QPushButton(this, "commitChanges");
  commitChanges->setText(trUtf8("Commit"));
  Layout5->addWidget(commitChanges);

  cancelChanges = new QPushButton(this, "cancelChanges");
  cancelChanges->setText(trUtf8("Revert"));
  Layout5->addWidget(cancelChanges);*/

  // signals and slots connections
  //connect(commitChanges, SIGNAL(clicked()), this, SLOT(slotBtnOKClicked()));
  //connect(cancelChanges, SIGNAL(clicked()), this, SLOT(slotBtnCancelClicked()));
}

/*
 *  Destroys the object and frees any allocated resources
 */
CQTSSAResultWidget::~CQTSSAResultWidget()
{}

bool CQTSSAResultWidget::loadFromBackend()
{

#if 0
  mCentralWidget->displayOptimizationTab(false);
  assert(CCopasiRootContainer::getDatamodelList()->size() > 0);
  mCentralWidget->table()->setTimeSeries(dynamic_cast<CTSSATask *>((*(*CCopasiRootContainer::getDatamodelList())[0]->getTaskList())["Time Scale Separation Analysis"])->getTimeSeries());
#endif

  return true;
}

bool CQTSSAResultWidget::saveToBackend()
{
  return true;
}

/*void CQTSSAResultWidget::slotBtnCancelClicked()
{
  enter(objKey); // reload
}*/

/*void CQTSSAResultWidget::slotBtnOKClicked()
{
  saveToCompartment();
}*/

bool CQTSSAResultWidget::update(ListViews::ObjectType /* objectType */,
                                ListViews::Action /* action */,
                                const std::string & /* key */)
{

#if 0

  if (objectType == ListViews::MODEL &&
      action == ListViews::ADD)
    mCentralWidget->table()->setTimeSeries(CTimeSeries());

#endif

  return true;
}

bool CQTSSAResultWidget::leave()
{
  //return saveToCompartment();
  return true;
}

bool CQTSSAResultWidget::enterProtected()
{
  pTask =
    dynamic_cast<CTSSATask *>((*(*CCopasiRootContainer::getDatamodelList())[0]->getTaskList())["Time Scale Separation Analysis"]);
  pTSSILDM = dynamic_cast<CTSSAMethod*>(pTask->getMethod());

  if (!pTSSILDM->mCurrentStep)
    {
#if 0
      mCentralWidget->setStepSelectionDisabled(true);
#endif
      mCentralWidget->discardOldResults();
    }

#if 0
  else
    mCentralWidget->setStepSelectionDisabled(false);

#endif

  return loadFromBackend();
  /*objKey = key;
  CCompartment* comp = dynamic_cast< CCompartment * >(CCopasiRootContainer::getKeyFactory()->get(key));

  if (comp) return loadFromCompartment(comp);
  else return false;*/
}
