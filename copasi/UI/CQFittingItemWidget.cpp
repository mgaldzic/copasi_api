// Begin CVS Header
//   $Source: /fs/turing/cvs/copasi_dev/copasi/UI/CQFittingItemWidget.cpp,v $
//   $Revision: 1.33 $
//   $Name: Build-33 $
//   $Author: shoops $
//   $Date: 2010/07/16 19:05:18 $
// End CVS Header

// Copyright (C) 2010 by Pedro Mendes, Virginia Tech Intellectual
// Properties, Inc., University of Heidelberg, and The University
// of Manchester.
// All rights reserved.

// Copyright (C) 2008 by Pedro Mendes, Virginia Tech Intellectual
// Properties, Inc., EML Research, gGmbH, University of Heidelberg,
// and The University of Manchester.
// All rights reserved.

#include "CQFittingItemWidget.h"

#include <qapplication.h>

#include <QHeaderView>

#include "CCopasiSelectionDialog.h"
#include "CQValidator.h"
#include "CQExperimentSelection.h"
#include "CQStartValueReset.h"
#include "qtUtilities.h"

#include "CopasiDataModel/CCopasiDataModel.h"
#include "report/CKeyFactory.h"
#include "parameterFitting/CFitItem.h"
#include "parameterFitting/CFitProblem.h"
#include "parameterFitting/CExperiment.h"
#include "parameterFitting/CExperimentSet.h"
#include "utilities/utility.h"
#include "copasi/report/CCopasiRootContainer.h"

#include <QtDebug>

//#include "UI/icons/Copasi16-Alpha.xpm"

/*
 *  Constructs a CQFittingItemWidget as a child of 'parent', with the
 *  name 'name' and widget flags set to 'f'.
 */
CQFittingItemWidget::CQFittingItemWidget(QWidget* parent, const char* name, Qt::WindowFlags fl)
    : QWidget(parent, name, fl)
{
  setupUi(this);
  /*
    mpBtnObject->setIcon(QPixmap(Copasi16_Alpha_xpm));
    mpBtnLowerEdit->setIcon(QPixmap(Copasi16_Alpha_xpm));
    mpBtnUpperEdit->setIcon(QPixmap(Copasi16_Alpha_xpm));
  //  mpBtnObject->setIcon(QPixmap(Copasi16_Alpha_xpm));
  */
  // set icon

  const QIcon icon = qt_get_icon(image0_ID);
  mpBtnDel->setIcon(icon);

  const QIcon icon1 = qt_get_icon(image1_ID);
  mpBtnUp->setIcon(icon1);

  const QIcon icon2 = qt_get_icon(image2_ID);
  mpBtnUpperEdit->setIcon(icon2);
  mpBtnLowerEdit->setIcon(icon2);
  mpBtnObject->setIcon(icon2);

  const QIcon icon3 = qt_get_icon(image3_ID);
  mpBtnNew->setIcon(icon3);

  const QIcon icon4 = qt_get_icon(image4_ID);
  mpBtnDown->setIcon(icon4);

  const QIcon icon5 = qt_get_icon(image5_ID);
  mpBtnCopy->setIcon(icon5);

  const QIcon icon6 = qt_get_icon(image6_ID);
  mpBtnReset->setIcon(icon6);
  mpBtnExperiments->setIcon(icon6);
  mpBtnCrossValidations->setIcon(icon6);

  init();
}

/*
 *  Destroys the object and frees any allocated resources
 */
CQFittingItemWidget::~CQFittingItemWidget()
{
  destroy();
  // no need to delete child widgets, Qt does it all for us
}

/*
 *  Sets the strings of the subwidgets using the current
 *  language.
 */
void CQFittingItemWidget::languageChange()
{
  retranslateUi(this);
}

void CQFittingItemWidget::init()
{
  mppExperimentSet = NULL;
  mppCrossValidationSet = NULL;

  mpItemsCopy = new std::vector<COptItem *>;
  mCurrentRow = (unsigned int) - 1;

  setItemType(OPT_ITEM);

  mLowerInfChanged = false;
  mUpperInfChanged = false;

  mpTable->horizontalHeader()->hide();
  mpTable->horizontalHeader()->setResizeMode(QHeaderView::Stretch);
  mpTable->setAlternatingRowColors(true);
  mpTable->setSelectionMode(QAbstractItemView::ExtendedSelection);
  mpTable->setSelectionBehavior(QAbstractItemView::SelectRows);
  mpTable->verticalHeader()->setResizeMode(QHeaderView::ResizeToContents);

  int h, s, v;

  mSavedColor = paletteBackgroundColor();
  mSavedColor.getHsv(&h, &s, &v);

  if (s < 20) s = 20;

  mChangedColor.setHsv(240, s, v);

  mpObjectValidator = new CQValidatorNotEmpty(mpEditObject);
  mpEditObject->setValidator(mpObjectValidator);
  mpObjectValidator->revalidate();

  mpLowerValidator = new CQValidatorBound(mpEditLower, "-");
  mpEditLower->setValidator(mpLowerValidator);
  mpLowerValidator->revalidate();

  mpUpperValidator = new CQValidatorBound(mpEditUpper, "+");
  mpEditUpper->setValidator(mpUpperValidator);
  mpUpperValidator->revalidate();

#ifndef COPASI_CROSSVALIDATION
  mpBtnCrossValidations->hide();
  mpCheckCrossValidationsAll->hide();
  mpBoxCrossValidations->hide();
  mpLblCrossValidations->hide();
#endif // not COPASI_CROSSVALIDATION
}

void CQFittingItemWidget::destroy()
{
  std::vector< COptItem * >::iterator it = mpItemsCopy->begin();
  std::vector< COptItem * >::iterator end = mpItemsCopy->end();

  for (; it != end; ++it)
    pdelete(*it);

  pdelete(mpItemsCopy);
}

void CQFittingItemWidget::slotCheckLowerInf(bool checked)
{
  mpEditLower->setEnabled(!checked);

  mLowerInfChanged = !mLowerInfChanged;

  if (mLowerInfChanged)
    {
      mpCheckLowerInf->setPaletteBackgroundColor(mChangedColor);

      if (mpEditLower->isEnabled()) mpLowerValidator->revalidate();
    }
  else
    mpCheckLowerInf->setPaletteBackgroundColor(mSavedColor);

  std::string Number;

  if (checked)
    Number = "-inf";
  else if (isNumber(TO_UTF8(mpEditLower->text())))
    Number = TO_UTF8(mpEditLower->text());
  else if (mpLowerObject)
    Number = mpLowerObject->getCN();
  else return;

  std::set< unsigned int >::const_iterator it = mSelection.begin();
  std::set< unsigned int >::const_iterator end = mSelection.end();

  for (; it != end; ++it)
    {
      (*mpItemsCopy)[*it]->setLowerBound(Number);
      setTableText(*it, (*mpItemsCopy)[*it]);
    }

  mpTable->resizeColumnsToContents();
}

void CQFittingItemWidget::slotCheckUpperInf(bool checked)
{
  mpEditUpper->setEnabled(!checked);

  mUpperInfChanged = !mUpperInfChanged;

  if (mUpperInfChanged)
    {
      mpCheckUpperInf->setPaletteBackgroundColor(mChangedColor);

      if (mpEditUpper->isEnabled()) mpUpperValidator->revalidate();
    }
  else
    mpCheckUpperInf->setPaletteBackgroundColor(mSavedColor);

  std::string Number;

  if (checked)
    Number = "inf";
  else if (isNumber(TO_UTF8(mpEditUpper->text())))
    Number = TO_UTF8(mpEditUpper->text());
  else if (mpUpperObject)
    Number = mpUpperObject->getCN();
  else return;

  std::set< unsigned int >::const_iterator it = mSelection.begin();
  std::set< unsigned int >::const_iterator end = mSelection.end();

  for (; it != end; ++it)
    {
      (*mpItemsCopy)[*it]->setUpperBound(Number);
      setTableText(*it, (*mpItemsCopy)[*it]);
    }

  mpTable->resizeColumnsToContents();
}

void CQFittingItemWidget::slotLowerEdit()
{
  CCopasiSimpleSelectionTree::ObjectClasses Classes;

  switch (mItemType)
    {
      case OPT_ITEM:
      case FIT_ITEM:
        Classes =
          CCopasiSimpleSelectionTree::InitialTime |
          CCopasiSimpleSelectionTree::Parameters |
          CCopasiSimpleSelectionTree::ObservedConstants;
        break;

      case OPT_CONSTRAINT:
      case FIT_CONSTRAINT:
        Classes =
          CCopasiSimpleSelectionTree::InitialTime |
          CCopasiSimpleSelectionTree::Parameters |
          CCopasiSimpleSelectionTree::ObservedConstants |
          CCopasiSimpleSelectionTree::Time |
          CCopasiSimpleSelectionTree::Variables |
          CCopasiSimpleSelectionTree::ObservedValues;
        break;
    }

  const CCopasiObject * pObject =
    CCopasiSelectionDialog::getObjectSingle(this, Classes);

  if (pObject)
    {
      mpCheckLowerInf->setChecked(false);
      mpEditLower->setEnabled(true);

      if (!mSelection.size())
        slotNew();

      std::set< unsigned int >::const_iterator it = mSelection.begin();
      std::set< unsigned int >::const_iterator end = mSelection.end();

      mpLowerObject = pObject;
      CCopasiObjectName CN = mpLowerObject->getCN();

      for (; it != end; ++it)
        {
          (*mpItemsCopy)[*it]->setLowerBound(CN);
          setTableText(*it, (*mpItemsCopy)[*it]);
        }

      mpTable->resizeColumnsToContents();

      QString Value = FROM_UTF8(pObject->getObjectDisplayName());
      mpLowerValidator->force(Value);
      mpEditLower->setText(Value);
    }
}

void CQFittingItemWidget::slotUpperEdit()
{
  CCopasiSimpleSelectionTree::ObjectClasses Classes;

  switch (mItemType)
    {
      case OPT_ITEM:
      case FIT_ITEM:
        Classes =
          CCopasiSimpleSelectionTree::InitialTime |
          CCopasiSimpleSelectionTree::Parameters |
          CCopasiSimpleSelectionTree::ObservedConstants;
        break;

      case OPT_CONSTRAINT:
      case FIT_CONSTRAINT:
        Classes =
          CCopasiSimpleSelectionTree::InitialTime |
          CCopasiSimpleSelectionTree::Parameters |
          CCopasiSimpleSelectionTree::ObservedConstants |
          CCopasiSimpleSelectionTree::Time |
          CCopasiSimpleSelectionTree::Variables |
          CCopasiSimpleSelectionTree::ObservedValues;
        break;
    }

  const CCopasiObject * pObject =
    CCopasiSelectionDialog::getObjectSingle(this, Classes);

  if (pObject)
    {
      mpCheckUpperInf->setChecked(false);
      mpEditUpper->setEnabled(true);

      if (!mSelection.size())
        slotNew();

      std::set< unsigned int >::const_iterator it = mSelection.begin();
      std::set< unsigned int >::const_iterator end = mSelection.end();

      mpUpperObject = pObject;
      CCopasiObjectName CN = mpUpperObject->getCN();

      for (; it != end; ++it)
        {
          (*mpItemsCopy)[*it]->setUpperBound(CN);
          setTableText(*it, (*mpItemsCopy)[*it]);
        }

      mpTable->resizeColumnsToContents();

      QString Value = FROM_UTF8(pObject->getObjectDisplayName());
      mpUpperValidator->force(Value);
      mpEditUpper->setText(Value);
    }
}

void CQFittingItemWidget::slotParamEdit()
{
  std::vector< const CCopasiObject * > Selection;

  CCopasiSimpleSelectionTree::ObjectClasses Classes;

  switch (mItemType)
    {
      case OPT_ITEM:
      case FIT_ITEM:
        Classes =
          CCopasiSimpleSelectionTree::InitialTime |
          CCopasiSimpleSelectionTree::Parameters;
        break;

      case OPT_CONSTRAINT:
      case FIT_CONSTRAINT:
        Classes =
          CCopasiSimpleSelectionTree::Time |
          CCopasiSimpleSelectionTree::Variables |
          CCopasiSimpleSelectionTree::ObservedValues |
          CCopasiSimpleSelectionTree::ObservedConstants;
        break;
    }

  if (mSelection.size() > 1)
    {
      const CCopasiObject * pObject =
        CCopasiSelectionDialog::getObjectSingle(this, Classes);

      if (pObject)
        Selection.push_back(pObject);
    }
  else
    Selection =
      CCopasiSelectionDialog::getObjectVector(this, Classes);

  if (Selection.size() != 0)
    {
      // We need to loop through the selection.
      unsigned C_INT32 current = currentRow();

      // We must not trigger slotSelectionChanged!
      disconnect(mpTable, SIGNAL(itemSelectionChanged()), this, SLOT(slotSelectionChanged()));

      // If no item is selected create one.
      if (!mSelection.size())
        {
          COptItem * pItem;

          switch (mItemType)
            {
              case OPT_ITEM:
              case OPT_CONSTRAINT:
                pItem = new COptItem(mpDataModel);
                break;

              case FIT_ITEM:
                pItem = new CFitItem(mpDataModel);
                break;

              case FIT_CONSTRAINT:
                pItem = new CFitConstraint(mpDataModel);
                break;
            }

          mpItemsCopy->push_back(pItem);
          current = mpTable->rowCount();
          mSelection.insert(current);
          mpTable->insertRow(current);
          setTableText(current, pItem);
        }

      // Update the selected items
      std::set< unsigned int >::const_iterator it = mSelection.begin();
      std::set< unsigned int >::const_iterator end = mSelection.end();

      for (; it != end; ++it)
        {
          (*mpItemsCopy)[*it]->setObjectCN(Selection[0]->getCN());
          setTableText(*it, (*mpItemsCopy)[*it]);
        }

      mpTable->resizeColumnsToContents();

      saveSelection();

      // If one item is selected we may have no several objects
      // and must create the appropriate items
      unsigned C_INT32 i, imax = Selection.size();

      COptItem * pSrc = (*mpItemsCopy)[*mSelection.begin()];

      for (i = 1; i != imax; i++)
        {
          COptItem * pItem;

          switch (mItemType)
            {
              case OPT_ITEM:
              case OPT_CONSTRAINT:
                pItem = new COptItem(*pSrc);
                break;

              case FIT_ITEM:
                pItem = new CFitItem(*pSrc);
                break;

              case FIT_CONSTRAINT:
                pItem = new CFitConstraint(*pSrc);
                break;
            }

          pItem->setObjectCN(Selection[i]->getCN());

          // Add the new item to the list.
          mpItemsCopy->insert(mpItemsCopy->begin() + current + i, pItem);

          // Update the table
          mpTable->insertRow(current + i);
          setTableText(current + i, pItem);
        }

      mpTable->resizeColumnsToContents();

      connect(mpTable, SIGNAL(itemSelectionChanged()), this, SLOT(slotSelectionChanged()));

      // Update the selection;
      if (current != C_INVALID_INDEX)
        selectRow(current);
      else
        loadSelection();

      emit numberChanged(mpItemsCopy->size());
    }
}

void CQFittingItemWidget::slotExperiments()
{
  if (mItemType == FIT_ITEM || mItemType == FIT_CONSTRAINT)
    {
      CQExperimentSelection * pDialog = new CQExperimentSelection(this);
      pDialog->load(mpBoxExperiments, * mppExperimentSet);

      if (pDialog->exec() == QDialog::Accepted)
        {
          if (!mSelection.size())
            slotNew();

          std::set< unsigned int >::const_iterator it = mSelection.begin();
          std::set< unsigned int >::const_iterator end = mSelection.end();

          for (; it != end; ++it)
            {
              (*mpItemsCopy)[*it]->getGroup("Affected Experiments")->clear();
              unsigned C_INT32 i, imax = mpBoxExperiments->count();

              for (i = 0; i < imax && imax < (*mppExperimentSet)->getExperimentCount(); i++)
                static_cast<CFitItem *>((*mpItemsCopy)[*it])->addExperiment((*mppExperimentSet)->getExperiment(TO_UTF8(mpBoxExperiments->text(i)))->CCopasiParameter::getKey());

              setTableText(*it, (*mpItemsCopy)[*it]);
            }

          mpTable->resizeColumnsToContents();
        }

      loadSelection();

      delete pDialog;
    }
}

bool CQFittingItemWidget::load(CCopasiDataModel * pDataModel,
                               CCopasiParameterGroup * pItems,
                               const std::map<std::string, std::string> * pExperimentMap,
                               const std::map<std::string, std::string> * pCrossValidationMap)
{
  mpDataModel = pDataModel;

  mpItems = pItems;
  mSelection.clear();

  std::vector< COptItem * >::iterator it = mpItemsCopy->begin();
  std::vector< COptItem * >::iterator end = mpItemsCopy->end();

  for (; it != end; ++it)
    pdelete(*it);

  // We must not trigger slotSelectionChanged!
  disconnect(mpTable, SIGNAL(itemSelectionChanged()), this, SLOT(slotSelectionChanged()));
  mpTable->setRowCount(mpItems->size());
  connect(mpTable, SIGNAL(itemSelectionChanged()), this, SLOT(slotSelectionChanged()));
  mCurrentRow = (unsigned int) - 1;

  mpItemsCopy->resize(mpItems->size());

  it = mpItemsCopy->begin();
  end = mpItemsCopy->end();

  std::vector< COptItem * >::const_iterator src =
    static_cast<std::vector< COptItem * > *>(mpItems->CCopasiParameter::getValue().pVOID)->begin();

  unsigned C_INT32 i;

  for (i = 0; it != end; ++it, ++src, ++i)
    {
      switch (mItemType)
        {
          case OPT_ITEM:
          case OPT_CONSTRAINT:
            *it = new COptItem(**src);
            break;

          case FIT_ITEM:
            *it = new CFitItem(**src);
            break;

          case FIT_CONSTRAINT:
            *it = new CFitConstraint(**src);
            break;
        }

      if (mItemType == FIT_ITEM || mItemType == FIT_CONSTRAINT)
        {
          if (!pExperimentMap) return false;

          // Change the key to reflect the local copy *mppExperimentSet
          unsigned C_INT32 j, jmax = static_cast<CFitItem *>(*it)->getExperimentCount();

          for (j = 0; j < jmax; j++)
            {
              std::string & Key =
                *const_cast<std::string *>(&static_cast<CFitItem *>(*it)->getExperiment(j));
              Key = pExperimentMap->find(Key)->second;
            }

        }

      setTableText(i, *it);
    }

  if (mpItemsCopy->size())
    selectRow(0);
  else
    selectRow((unsigned int) - 1);

  mpTable->resizeColumnsToContents();

  emit numberChanged(mpItemsCopy->size());

  return true;
}

bool CQFittingItemWidget::save(const std::map<std::string, std::string> * pExperimentMap, const std::map<std::string, std::string> * pCrossValidationMap)
{
  // Make sure that the current items is saved.
  saveSelection();

  bool changed = false;

  std::vector< COptItem * >::const_iterator it = mpItemsCopy->begin();
  std::vector< COptItem * >::const_iterator end = mpItemsCopy->end();

  std::vector< COptItem * >::iterator target =
    static_cast<std::vector< COptItem * > *>(mpItems->CCopasiParameter::getValue().pVOID)->begin();
  std::vector< COptItem * >::iterator targetEnd =
    static_cast<std::vector< COptItem * > *>(mpItems->CCopasiParameter::getValue().pVOID)->end();

  unsigned C_INT32 i;
  unsigned C_INT32 imax = std::max<unsigned C_INT32>(mpItemsCopy->size(), mpItems->size());

  for (i = 0; it != end && target != targetEnd; ++it, ++target, ++i)
    {
      if ((*target)->getObjectCN() != (*it)->getObjectCN())
        {
          changed = true;

          if (!(*target)->setObjectCN((*it)->getObjectCN()))
            (*target)->setValue("ObjectCN", (*it)->getObjectCN());
        }

      if ((*target)->getLowerBound() != (*it)->getLowerBound())
        {
          changed = true;
          (*target)->setLowerBound((*it)->getLowerBound());
        }

      if ((*target)->getUpperBound() != (*it)->getUpperBound())
        {
          changed = true;
          (*target)->setUpperBound((*it)->getUpperBound());
        }

      if ((*target)->getStartValue() != (*it)->getStartValue())
        {
          changed = true;
          (*target)->setStartValue((*it)->getStartValue());
        }

      if (mItemType == FIT_ITEM || mItemType == FIT_CONSTRAINT)
        {
          if (pExperimentMap == NULL) return false;

          unsigned C_INT32 j, jmax =
            std::min(static_cast<CFitItem *>(*it)->getExperimentCount(),
                     static_cast<CFitItem *>(*target)->getExperimentCount());

          // Compare and assign common experiments
          for (j = 0; j < jmax; j++)
            {
              std::string & Target =
                *const_cast<std::string *>(&static_cast<CFitItem *>(*target)->getExperiment(j));
              const std::string & Source =
                pExperimentMap->find(static_cast<CFitItem *>(*it)->getExperiment(j))->second;

              if (Target != Source)
                {
                  changed = true;
                  Target = Source;
                }
            }

          // Remove exceeding experiments starting from the last
          for (jmax = static_cast<CFitItem *>(*target)->getExperimentCount() - 1;
               j <= jmax && jmax != C_INVALID_INDEX; jmax--)
            {
              changed = true;
              static_cast<CFitItem *>(*target)->removeExperiment(jmax);
            }

          // Add missing experiments
          for (jmax = static_cast<CFitItem *>(*it)->getExperimentCount(); j < jmax; j++)
            {
              changed = true;

              const std::string & Source =
                pExperimentMap->find(static_cast<CFitItem *>(*it)->getExperiment(j))->second;

              static_cast<CFitItem *>(*target)->addExperiment(Source);
            }

        }
    }

  // Remove exceeding items starting from the last.
  for (imax--; i <= imax && imax != C_INVALID_INDEX; imax--)
    {
      changed = true;
      mpItems->removeParameter(imax);
    }

  // Add missing items
  for (; it != end; ++it)
    {
      changed = true;

      COptItem * pItem = NULL;

      switch (mItemType)
        {
          case OPT_ITEM:
          case OPT_CONSTRAINT:
            pItem = new COptItem(**it);
            break;

          case FIT_ITEM:
            pItem = new CFitItem(**it);
            break;

          case FIT_CONSTRAINT:
            pItem = new CFitConstraint(**it);
            break;
        }

      if (mItemType == FIT_ITEM || mItemType == FIT_CONSTRAINT)
        {
          if (pExperimentMap == NULL) return false;

          unsigned C_INT32 j, jmax =
            static_cast<CFitItem *>(pItem)->getExperimentCount();

          // Remap the experiment keys
          for (j = 0; j < jmax; j++)
            {
              std::string & Target =
                *const_cast<std::string *>(&static_cast<CFitItem *>(pItem)->getExperiment(j));
              const std::string & Source =
                pExperimentMap->find(static_cast<CFitItem *>(*it)->getExperiment(j))->second;

              Target = Source;
            }
        }

      mpItems->addParameter(pItem);
    }

  return changed;
}

void CQFittingItemWidget::setItemType(const ItemType & type)
{
  mItemType = type;

  if (mItemType == FIT_ITEM || mItemType == FIT_CONSTRAINT)
    {
      mpLblExperiments->show();
      mpCheckAll->show();
      mpBoxExperiments->show();
      mpBtnExperiments->show();
      mpBtnPerExperiment->show();

    }
  else
    {
      mpLblExperiments->hide();
      mpCheckAll->hide();
      mpBoxExperiments->hide();
      mpBtnExperiments->hide();
      mpBtnPerExperiment->hide();

    }

  if (mItemType == OPT_CONSTRAINT || mItemType == FIT_CONSTRAINT)
    {
      mpLblStart->hide();
      mpEditStart->hide();
      mpBtnReset->hide();
    }
  else
    {
      mpLblStart->show();
      mpEditStart->show();
      mpBtnReset->show();
    }

  qApp->processEvents();
}

void CQFittingItemWidget::slotExperimentChanged()
{
  // This slot is triggered when an experiment is deleted or changed,
  // but before new experiments are created.

  if (mItemType == OPT_ITEM || mItemType == OPT_CONSTRAINT) return;

  // Remove all references to deleted experiments.
  std::vector< COptItem * >::iterator it = mpItemsCopy->begin();
  std::vector< COptItem * >::iterator end = mpItemsCopy->end();

  unsigned C_INT32 i, imax;
  unsigned C_INT32 Row;

  for (Row = 0; it != end; ++it, ++Row)
    {
      for (i = 0, imax = static_cast<CFitItem *>(*it)->getExperimentCount(); i < imax; ++i)
        if (!CCopasiRootContainer::getKeyFactory()->get(static_cast<CFitItem *>(*it)->getExperiment(i)))
          static_cast<CFitItem *>(*it)->removeExperiment(i);

      setTableText(Row, *it);
    }

  mpTable->resizeColumnsToContents();

  // Reload the current item.
  loadSelection();

  return;
}

void CQFittingItemWidget::setExperimentSet(const CExperimentSet * & pExperimentSet)
{mppExperimentSet = &pExperimentSet;}

void CQFittingItemWidget::slotDelete()
{
  unsigned C_INT32 row = mCurrentRow;

  COptItem * pItem;

  std::set< unsigned int >::const_iterator it = mSelection.begin();
  std::set< unsigned int >::const_iterator end = mSelection.end();

  // We must not trigger slotSelectionChanged!
  disconnect(mpTable, SIGNAL(itemSelectionChanged()), this, SLOT(slotSelectionChanged()));

  unsigned int Deleted = 0;

  for (; it != end; ++it, ++Deleted)
    {
      row = *it - Deleted;

      pItem = (*mpItemsCopy)[row];

      // Erase the item from the vector.
      mpItemsCopy->erase(mpItemsCopy->begin() + row);
      // Update the table
      mpTable->removeRow(row);
      pdelete(pItem);
    }

  mpTable->clearSelection();
  mSelection.clear();

  connect(mpTable, SIGNAL(itemSelectionChanged()), this, SLOT(slotSelectionChanged()));

  // Change the selection
  if (row >= mpItemsCopy->size())
    row = mpItemsCopy->size() - 1;

  selectRow(row);

  emit numberChanged(mpItemsCopy->size());
}

void CQFittingItemWidget::slotCopy()
{
  if (!mSelection.size())
    {
      // We have nothing to copy, therefore we create a new item.
      slotNew();
      return;
    }

  // Save the changes to  the current items
  saveSelection();

  // Create the new item.
  COptItem * pItem;
  COptItem * pSrc = (*mpItemsCopy)[*mSelection.begin()];

  switch (mItemType)
    {
      case OPT_ITEM:
      case OPT_CONSTRAINT:
        pItem = new COptItem(*pSrc);
        break;
      case FIT_ITEM:
        pItem = new CFitItem(*pSrc);
        break;
      case FIT_CONSTRAINT:
        pItem = new CFitConstraint(*pSrc);
        break;
    }

  unsigned C_INT32 row = mCurrentRow + 1;
  mpItemsCopy->insert(mpItemsCopy->begin() + row, pItem);

  // Update the table
  // We must not trigger slotSelectionChanged!
  disconnect(mpTable, SIGNAL(itemSelectionChanged()), this, SLOT(slotSelectionChanged()));
//  mpTable->insertRows(row);
  mpTable->insertRow(row);
  connect(mpTable, SIGNAL(itemSelectionChanged()), this, SLOT(slotSelectionChanged()));

  setTableText(row, pItem);
  mpTable->resizeColumnsToContents();

  // Update the selection
  selectRow(row);

  emit numberChanged(mpItemsCopy->size());
}

void CQFittingItemWidget::slotUp()
{
  if (!mSelection.size()) return; // Nothing to move.

  saveSelection();

  std::set< unsigned int > NewSelection;
  std::set< unsigned int >::const_iterator it = mSelection.begin();
  std::set< unsigned int >::const_iterator end = mSelection.end();
  COptItem * pItem;

  if (*it == 0) ++it; // The first row can not be moved up.

  for (; it != end; ++it)
    {
      NewSelection.insert(*it - 1);

      // Swap the items
      pItem = (*mpItemsCopy)[*it];
      (*mpItemsCopy)[*it] = (*mpItemsCopy)[*it - 1];
      (*mpItemsCopy)[*it - 1] = pItem;

      // Swap the table rows
      QTableWidgetItem *item1 = mpTable->takeItem(*it - 1, 0);
      QTableWidgetItem *item2 = mpTable->takeItem(*it, 0);
      mpTable->setItem(*it, 0, item1);
      mpTable->setItem(*it - 1, 0, item2);
    }

  setItemSelection(NewSelection);
  loadSelection();
}

void CQFittingItemWidget::slotDown()
{
  if (!mSelection.size()) return; // Nothing to move.

  saveSelection();

  std::set< unsigned int > NewSelection;
  std::set< unsigned int >::reverse_iterator it = mSelection.rbegin();
  std::set< unsigned int >::reverse_iterator end = mSelection.rend();
  COptItem * pItem;

  if (*it == (unsigned int)(mpTable->rowCount() - 1))
    ++it; // The last row can not be moved down.

  for (; it != end; ++it)
    {
      NewSelection.insert(*it + 1);

      // Swap the items
      pItem = (*mpItemsCopy)[*it];
      (*mpItemsCopy)[*it] = (*mpItemsCopy)[*it + 1];
      (*mpItemsCopy)[*it + 1] = pItem;

      // Swap the table rows
      QTableWidgetItem *item1 = mpTable->takeItem(*it + 1, 0);
      QTableWidgetItem *item2 = mpTable->takeItem(*it, 0);
      mpTable->setItem(*it, 0, item1);
      mpTable->setItem(*it + 1, 0, item2);
    }

  setItemSelection(NewSelection);
  loadSelection();
}

void CQFittingItemWidget::slotDuplicatePerExperiment()
{
  if (mItemType == OPT_ITEM || mItemType == OPT_CONSTRAINT) return;

  // Save the changes to  the current items
  saveSelection();

  if (!mSelection.size())
    slotNew();

  unsigned C_INT32 row = mCurrentRow + 1;
  unsigned C_INT32 i, imax =
    static_cast<CFitItem *>((*mpItemsCopy)[*mSelection.begin()])->getExperimentCount();

  CFitItem * pItem = NULL, * pTemplate = NULL;

  // We must not trigger slotSelectionChanged!
  disconnect(mpTable, SIGNAL(itemSelectionChanged()), this, SLOT(slotSelectionChanged()));

  if (imax)
    {
      // We have a list of experiments
      switch (mItemType)
        {
          case FIT_ITEM:
            pTemplate = new CFitItem(*(*mpItemsCopy)[*mSelection.begin()]);
            break;
          case FIT_CONSTRAINT:
            pTemplate = new CFitConstraint(*(*mpItemsCopy)[*mSelection.begin()]);
            break;
          default:
            break;
        }

      // Remove all experiments from the template
      for (i = imax - 1; i != C_INVALID_INDEX; i--)
        pTemplate->removeExperiment(i);

      // Create imax - 1 copies
      for (i = imax - 1; i != 0; i--)
        {
          switch (mItemType)
            {
              case FIT_ITEM:
                pItem = new CFitItem(*pTemplate);
                break;
              case FIT_CONSTRAINT:
                pItem = new CFitConstraint(*pTemplate);
                break;
              default:
                break;
            }

          pItem->addExperiment(static_cast<CFitItem *>((*mpItemsCopy)[*mSelection.begin()])->getExperiment(i));

          static_cast<CFitItem *>((*mpItemsCopy)[*mSelection.begin()])->removeExperiment(i);

          mpItemsCopy->insert(mpItemsCopy->begin() + row, pItem);

          // Update the table
          mpTable->insertRow(row);
          setTableText(row, pItem);
        }
    }
  else
    {
      // We have ALL
      imax = (*mppExperimentSet)->getExperimentCount();

      if (imax > 0)
        {
          // Create imax - 1 copies
          for (i = imax - 1; i != 0; i--)
            {
              pItem = new CFitItem(*(*mpItemsCopy)[*mSelection.begin()]);
              pItem->addExperiment((*mppExperimentSet)->getExperiment(i)->CCopasiParameter::getKey());

              mpItemsCopy->insert(mpItemsCopy->begin() + row, pItem);

              // Update the table
              mpTable->insertRow(row);
              setTableText(row, pItem);
            }
        }

      // Update the current item only when we have more than one experiment.
      if (imax > 1)
        static_cast<CFitItem *>((*mpItemsCopy)[*mSelection.begin()])->addExperiment((*mppExperimentSet)->getExperiment(0)->CCopasiParameter::getKey());
    }

  connect(mpTable, SIGNAL(itemSelectionChanged()), this, SLOT(slotSelectionChanged()));

  // Update the current item
  setTableText(row - 1, (*mpItemsCopy)[*mSelection.begin()]);

  // Update the selection
  selectRow(row - 1);

  mpTable->resizeColumnsToContents();

  emit numberChanged(mpItemsCopy->size());
}

void CQFittingItemWidget::slotNew()
{
  // Save the changes to  the current items
  saveSelection();

  // Create the new item.
  COptItem * pItem;

  switch (mItemType)
    {
      case OPT_ITEM:
      case OPT_CONSTRAINT:
        pItem = new COptItem(mpDataModel);
        break;

      case FIT_ITEM:
        pItem = new CFitItem(mpDataModel);
        break;

      case FIT_CONSTRAINT:
        pItem = new CFitConstraint(mpDataModel);
        break;
    }

  mpItemsCopy->push_back(pItem);

  unsigned C_INT32 row = mpTable->rowCount();

  // Update the table
  // We must not trigger slotSelectionChanged!
  disconnect(mpTable, SIGNAL(itemSelectionChanged()), this, SLOT(slotSelectionChanged()));
  // mpTable->insertRows(row);
  mpTable->insertRow(row);
  connect(mpTable, SIGNAL(itemSelectionChanged()), this, SLOT(slotSelectionChanged()));

  setTableText(row, pItem);
  mpTable->resizeColumnsToContents();

  // Update the selection
  selectRow(row);

  emit numberChanged(mpItemsCopy->size());
}

void CQFittingItemWidget::setTableText(const int & row, const COptItem * pItem)
{
  QString Item = "   ";
  const CCopasiObject *pObject;

  if (pItem->getLowerBound() == "-inf" ||
      isNumber(pItem->getLowerBound()))
    Item += FROM_UTF8(pItem->getLowerBound());
  else if ((pObject = mpDataModel->getObject(pItem->getLowerBound())))
    Item += FROM_UTF8(pObject->getObjectDisplayName());
  else
    Item += "Not found: " + FROM_UTF8(pItem->getLowerBound());

  // Insert less than character
  Item += FROM_UTF8(std::string(" \xe2\x89\xa4 "));

  pObject = mpDataModel->getObject(pItem->getObjectCN());

  if (pObject)
    Item += FROM_UTF8(pObject->getObjectDisplayName());
  else
    Item += "Not found: " + FROM_UTF8(pItem->getObjectCN());

  if (mItemType == FIT_ITEM || mItemType == FIT_CONSTRAINT)
    {
      QString Experiments =
        FROM_UTF8(static_cast<const CFitItem *>(pItem)->getExperiments());

      if (Experiments != "")
        Item += "; {" + Experiments + "}";
    }

  // Insert less than character
  Item += FROM_UTF8(std::string(" \xe2\x89\xa4 "));

  if (pItem->getUpperBound() == "inf" ||
      isNumber(pItem->getUpperBound()))
    Item += FROM_UTF8(pItem->getUpperBound());
  else if ((pObject = mpDataModel->getObject(pItem->getUpperBound())))
    Item += FROM_UTF8(pObject->getObjectDisplayName());
  else
    Item += "Not found: " + FROM_UTF8(pItem->getUpperBound());

  if (mItemType == OPT_ITEM || mItemType == FIT_ITEM)
    {
      Item += "; Start Value = " + QString::number(pItem->getStartValue());
    }

  QTableWidgetItem *itemValue = new QTableWidgetItem(Item);
  itemValue->setFlags(Qt::ItemIsSelectable | Qt::ItemIsEnabled);
  itemValue->setTextAlignment(Qt::AlignLeft | Qt::AlignVCenter);
  mpTable->setItem(row, 0, itemValue);
}

unsigned int CQFittingItemWidget::currentRow()
{
  if (mSelection.size() != 1)
    return (unsigned int) - 1;

  return *mSelection.begin();
}

void CQFittingItemWidget::slotSelectionChanged()
{
  saveSelection();

  unsigned int i, imax = mpTable->rowCount();

  // Update mSelection;
  for (i = 0; i != imax; i++)

    if (mpTable->item(i, 0)->isSelected())
      mSelection.insert(i);
    else
      mSelection.erase(i);

  mCurrentRow = currentRow();

  loadSelection();
}

void CQFittingItemWidget::loadSelection()
{
  bool Enabled = (mSelection.size() <= 1);

  mpBtnCopy->setEnabled(Enabled);
  mpBtnPerExperiment->setEnabled(Enabled);

  std::set< unsigned int >::const_iterator it = mSelection.begin();
  std::set< unsigned int >::const_iterator end = mSelection.end();

  if (it == end)
    {
      mpEditObject->setText("");

      mpEditLower->setText("");
      mpCheckLowerInf->setChecked(true);
      mpLowerObject = NULL;

      mpEditUpper->setText("");
      mpCheckUpperInf->setChecked(true);
      mpUpperObject = NULL;

      mpEditStart->setText("");

      mpCheckAll->setChecked(true);
      mpBoxExperiments->setEnabled(false);

      mpCheckCrossValidationsAll->setChecked(true);
      mpBoxCrossValidations->setEnabled(false);

      QString Empty("");
      int Pos = 0;
      mpObjectValidator->validate(Empty, Pos);
    }
  else
    {
      disconnect(mpCheckAll, SIGNAL(toggled(bool)), this, SLOT(slotCheckAllExperiments(bool)));
      disconnect(mpCheckCrossValidationsAll, SIGNAL(toggled(bool)), this, SLOT(slotCheckAllCrossValidations(bool)));
      QString Value;

      COptItem * pItem = (*mpItemsCopy)[*it];

      const CCopasiObject *pObject = mpDataModel->getObject(pItem->getObjectCN());

      if (pObject)
        {
          Value = FROM_UTF8(pObject->getObjectDisplayName());
          mpEditObject->setText(Value);
          mpObjectValidator->saved();
        }
      else
        {
          if (pItem->getObjectCN() != "")
            Value = "Not found: " + FROM_UTF8(pItem->getObjectCN());
          else
            Value = "";

          mpEditObject->setText(Value);

          int Pos = 0;
          mpObjectValidator->validate(Value, Pos);
        }

      mpLowerObject = NULL;

      if (pItem->getLowerBound() == "-inf" ||
          isNumber(pItem->getLowerBound()))
        Value = FROM_UTF8(pItem->getLowerBound());
      else if ((mpLowerObject = mpDataModel->getObject(pItem->getLowerBound())))
        Value = FROM_UTF8(mpLowerObject->getObjectDisplayName());
      else
        Value = "Not found: " + FROM_UTF8(pItem->getLowerBound());

      mpEditLower->setText(Value);
      mpCheckLowerInf->setChecked(Value == "-inf");

      mpUpperObject = NULL;

      if (pItem->getUpperBound() == "inf" ||
          isNumber(pItem->getUpperBound()))
        Value = FROM_UTF8(pItem->getUpperBound());
      else if ((mpUpperObject = mpDataModel->getObject(pItem->getUpperBound())))
        Value = FROM_UTF8(mpUpperObject->getObjectDisplayName());
      else
        Value = "Not found: " + FROM_UTF8(pItem->getUpperBound());

      mpEditUpper->setText(Value);
      mpCheckUpperInf->setChecked(Value == "inf");

      mpEditStart->setText(QString::number(pItem->getStartValue()));

      std::string Experiments;
      std::string CrossValidations;

      if (mItemType == FIT_ITEM || mItemType == FIT_CONSTRAINT)
        {
          Experiments = static_cast<CFitItem *>(pItem)->getExperiments();

          mpBoxExperiments->clear();

          unsigned C_INT32 i, imax = static_cast<CFitItem *>(pItem)->getExperimentCount();

          for (i = 0; i < imax; i++)
            {
              const CCopasiObject * pObject =
                CCopasiRootContainer::getKeyFactory()->get(static_cast<CFitItem *>(pItem)->getExperiment(i));

              if (pObject)
                mpBoxExperiments->insertItem(FROM_UTF8(pObject->getObjectName()));
            }

          mpCheckAll->setChecked(imax == 0);
          mpBoxExperiments->setEnabled(imax != 0);

        }

      for (++it; it != end; ++it)
        {
          pItem = (*mpItemsCopy)[*it];

          const CCopasiObject *pObject = mpDataModel->getObject(pItem->getObjectCN());

          if (pObject)
            Value = FROM_UTF8(pObject->getObjectDisplayName());
          else
            {
              if (pItem->getObjectCN() != "")
                Value = "Not found: " + FROM_UTF8(pItem->getObjectCN());
              else
                Value = "";
            }

          if (Value != mpEditObject->text())
            {
              mpEditObject->setText("");
              mpObjectValidator->saved();
            }

          if (pItem->getLowerBound() == "-inf" ||
              isNumber(pItem->getLowerBound()))
            Value = FROM_UTF8(pItem->getLowerBound());
          else if ((pObject = mpDataModel->getObject(pItem->getLowerBound())))
            Value = FROM_UTF8(pObject->getObjectDisplayName());
          else
            Value = "Not found: " + FROM_UTF8(pItem->getLowerBound());

          if (Value != mpEditLower->text())
            {
              mpEditLower->setText("");
              mpCheckLowerInf->setChecked(false);
            }

          if (pItem->getUpperBound() == "inf" ||
              isNumber(pItem->getUpperBound()))
            Value = FROM_UTF8(pItem->getUpperBound());
          else if ((pObject = mpDataModel->getObject(pItem->getUpperBound())))
            Value = FROM_UTF8(pObject->getObjectDisplayName());
          else
            Value = "Not found: " + FROM_UTF8(pItem->getUpperBound());

          if (Value != mpEditUpper->text())
            {
              mpEditUpper->setText("");
              mpCheckUpperInf->setChecked(false);
            }

          if (QString::number(pItem->getStartValue()) != mpEditStart->text())
            mpEditStart->setText("");

          if ((mItemType == FIT_ITEM || mItemType == FIT_CONSTRAINT) &&
              Experiments != static_cast<CFitItem *>(pItem)->getExperiments())
            {
              mpCheckAll->setChecked(false);
              mpBoxExperiments->setEnabled(false);
            }

        }

      connect(mpCheckAll, SIGNAL(toggled(bool)), this, SLOT(slotCheckAllExperiments(bool)));
      connect(mpCheckCrossValidationsAll, SIGNAL(toggled(bool)), this, SLOT(slotCheckAllCrossValidations(bool)));
    }

  mpLowerValidator->saved();
  mpUpperValidator->saved();
  mpCheckLowerInf->setPaletteBackgroundColor(mSavedColor);
  mLowerInfChanged = false;
  mpCheckUpperInf->setPaletteBackgroundColor(mSavedColor);
  mUpperInfChanged = false;
}

void CQFittingItemWidget::saveSelection()
{
  std::set< unsigned int >::const_iterator it = mSelection.begin();
  std::set< unsigned int >::const_iterator end = mSelection.end();

  COptItem * pItem;

  for (; it != end; ++it)
    {
      pItem = (*mpItemsCopy)[*it];

      if (mpCheckLowerInf->isChecked())
        pItem->setLowerBound(CCopasiObjectName("-inf"));
      else if (isNumber(TO_UTF8(mpEditLower->text())))
        pItem->setLowerBound(CCopasiObjectName(TO_UTF8(mpEditLower->text())));

      if (mpCheckUpperInf->isChecked())
        pItem->setUpperBound(CCopasiObjectName("inf"));
      else if (isNumber(TO_UTF8(mpEditUpper->text())))
        pItem->setUpperBound(CCopasiObjectName(TO_UTF8(mpEditUpper->text())));

      if (isNumber(TO_UTF8(mpEditStart->text())))
        pItem->setStartValue(mpEditStart->text().toDouble());

      mpCheckLowerInf->setPaletteBackgroundColor(mSavedColor);
      mLowerInfChanged = false;
      mpCheckUpperInf->setPaletteBackgroundColor(mSavedColor);
      mUpperInfChanged = false;
      mpObjectValidator->saved();
      mpLowerValidator->saved();
      mpUpperValidator->saved();

      setTableText(*it, pItem);
    }

  mpTable->resizeColumnsToContents();
}

void CQFittingItemWidget::selectRow(const unsigned int & row)
{
  // This sets the focus and creates a selection
  mpTable->setCurrentCell(row, 0);

  // We must not trigger slotSelectionChanged!
  disconnect(mpTable, SIGNAL(itemSelectionChanged()), this, SLOT(slotSelectionChanged()));
  mpTable->clearSelection();
  connect(mpTable, SIGNAL(itemSelectionChanged()), this, SLOT(slotSelectionChanged()));

  mpTable->selectRow(row);
}

void CQFittingItemWidget::setItemSelection(const std::set< unsigned int > & selection)
{
  mSelection = selection;

  std::set< unsigned int >::const_iterator it = selection.begin();
  std::set< unsigned int >::const_iterator end = selection.end();

  // We must not trigger slotSelectionChanged!
  disconnect(mpTable, SIGNAL(itemSelectionChanged()), this, SLOT(slotSelectionChanged()));

  // This sets the focus and creates a selection
  if (it != end)
    mpTable->setCurrentCell(*it, 0);
  else
    mpTable->setCurrentCell(-1, 0);

  // Update selection
  for (; it != end; ++it)
    mpTable->item(*it, 0)->setSelected(true);

  connect(mpTable, SIGNAL(itemSelectionChanged()), this, SLOT(slotSelectionChanged()));
}

void CQFittingItemWidget::slotLowerLostFocus()
{
  std::string Number = TO_UTF8(mpEditLower->text());

  if (!isNumber(Number) &&
      !(Number[0] == '-' &&
        Number[Number.length() - 1] == '%' &&
        isNumber(Number.substr(1, Number.length() - 2)))) return;

  mpLowerObject = NULL;

  bool first = true;
  std::string NewValue = "";

  std::set< unsigned int >::const_iterator it = mSelection.begin();
  std::set< unsigned int >::const_iterator end = mSelection.end();

  for (; it != end; ++it)
    {
      (*mpItemsCopy)[*it]->setLowerBound(Number);

      if (first)
        {
          NewValue = (*mpItemsCopy)[*it]->getLowerBound();
          first = false;
        }
      else if (NewValue != (*mpItemsCopy)[*it]->getLowerBound())
        NewValue = "";

      setTableText(*it, (*mpItemsCopy)[*it]);
    }

  mpTable->resizeColumnsToContents();

  mpEditLower->setText(FROM_UTF8(NewValue));
}

void CQFittingItemWidget::slotUpperLostFocus()
{
  std::string Number = TO_UTF8(mpEditUpper->text());

  if (!isNumber(Number) &&
      !(Number[0] == '+' &&
        Number[Number.length() - 1] == '%' &&
        isNumber(Number.substr(1, Number.length() - 2)))) return;

  mpUpperObject = NULL;

  std::set< unsigned int >::const_iterator it = mSelection.begin();
  std::set< unsigned int >::const_iterator end = mSelection.end();

  bool first = true;
  std::string NewValue = "";

  for (; it != end; ++it)
    {
      (*mpItemsCopy)[*it]->setUpperBound(Number);

      if (first)
        {
          NewValue = (*mpItemsCopy)[*it]->getUpperBound();
          first = false;
        }
      else if (NewValue != (*mpItemsCopy)[*it]->getUpperBound())
        NewValue = "";

      setTableText(*it, (*mpItemsCopy)[*it]);
    }

  mpTable->resizeColumnsToContents();

  mpEditUpper->setText(FROM_UTF8(NewValue));
}

void CQFittingItemWidget::slotReset()
{
  CQStartValueReset * pDialog = new CQStartValueReset(this);
  pDialog->exec();

  std::set< unsigned int >::const_iterator it = mSelection.begin();
  std::set< unsigned int >::const_iterator end = mSelection.end();

  switch (pDialog->result())
    {
      case QDialog::Rejected:
        break;

      case CQStartValueReset::MODEL:

        for (; it != end; ++it)
          {
            (*mpItemsCopy)[*it]->setStartValue(std::numeric_limits<C_FLOAT64>::quiet_NaN());
            mpEditStart->setText(QString::number((*mpItemsCopy)[*it]->getStartValue()));
            setTableText(*it, (*mpItemsCopy)[*it]);
          }

        break;

      case CQStartValueReset::RANDOM:

        for (; it != end; ++it)
          {
            (*mpItemsCopy)[*it]->setStartValue((*mpItemsCopy)[*it]->getRandomValue());
            mpEditStart->setText(QString::number((*mpItemsCopy)[*it]->getStartValue()));
            setTableText(*it, (*mpItemsCopy)[*it]);
          }

        break;

      case CQStartValueReset::SOLUTION:
      {
        COptProblem * pProblem = dynamic_cast< COptProblem * >(mpItems->getObjectParent());
        const CVector< C_FLOAT64 > & Solution = pProblem->getSolutionVariables();

        if (Solution.size() == mpItems->size())
          for (; it != end; ++it)
            {
              (*mpItemsCopy)[*it]->setStartValue(Solution[*it]);
              mpEditStart->setText(QString::number((*mpItemsCopy)[*it]->getStartValue()));
              setTableText(*it, (*mpItemsCopy)[*it]);
            }
      }
      break;
    }

  mpTable->resizeColumnsToContents();

  pdelete(pDialog);

  loadSelection();
}

void CQFittingItemWidget::slotStartLostFocus()
{
  if (!isNumber(TO_UTF8(mpEditStart->text()))) return;

  C_FLOAT64 Number = mpEditStart->text().toDouble();

  std::set< unsigned int >::const_iterator it = mSelection.begin();
  std::set< unsigned int >::const_iterator end = mSelection.end();

  for (; it != end; ++it)
    {
      (*mpItemsCopy)[*it]->setStartValue(Number);
      setTableText(*it, (*mpItemsCopy)[*it]);
    }

  mpTable->resizeColumnsToContents();
}

void CQFittingItemWidget::slotCrossValidations()
{
}

void CQFittingItemWidget::slotCrossValidationChanged()
{
  return;
}


void CQFittingItemWidget::slotCheckAllCrossValidations(bool checked)
{
  if (mItemType == OPT_ITEM || mItemType == OPT_CONSTRAINT) return;

  if (!checked && mpBoxCrossValidations->count() == 0)
    slotCrossValidations();
  else if (checked)
    {
      std::set< unsigned int >::const_iterator it = mSelection.begin();
      std::set< unsigned int >::const_iterator end = mSelection.end();

      for (; it != end; ++it)
        {
          (*mpItemsCopy)[*it]->getGroup("Affected Cross Validation Experiments")->clear();
          setTableText(*it, (*mpItemsCopy)[*it]);
        }
    }

  mpTable->resizeColumnsToContents();

  loadSelection();
}

void CQFittingItemWidget::slotCheckAllExperiments(bool checked)
{
  if (mItemType == OPT_ITEM || mItemType == OPT_CONSTRAINT) return;

  if (!checked && mpBoxExperiments->count() == 0)
    slotExperiments();
  else if (checked)
    {
      std::set< unsigned int >::const_iterator it = mSelection.begin();
      std::set< unsigned int >::const_iterator end = mSelection.end();

      for (; it != end; ++it)
        {
          (*mpItemsCopy)[*it]->getGroup("Affected Experiments")->clear();
          setTableText(*it, (*mpItemsCopy)[*it]);
        }
    }

  mpTable->resizeColumnsToContents();

  loadSelection();
}
