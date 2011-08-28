// Begin CVS Header
//   $Source: /fs/turing/cvs/copasi_dev/copasi/UI/CQEventsWidget.cpp,v $
//   $Revision: 1.24 $
//   $Name: Build-33 $
//   $Author: shoops $
//   $Date: 2010/09/08 17:35:02 $
// End CVS Header

// Copyright (C) 2010 by Pedro Mendes, Virginia Tech Intellectual
// Properties, Inc., University of Heidelberg, and The University
// of Manchester.
// All rights reserved.

// Copyright (C) 2008 by Pedro Mendes, Virginia Tech Intellectual
// Properties, Inc., EML Research, gGmbH, University of Heidelberg,
// and The University of Manchester.
// All rights reserved.

#include <QHeaderView>
#include <QClipboard>

#include "model/CModel.h"
#include "CopasiDataModel/CCopasiDataModel.h"
#include "report/CCopasiRootContainer.h"

#include "CQEventsWidget.h"
#include "qtUtilities.h"
#include "copasi.h"
#include "CQMessageBox.h"

/*
 *  Constructs a CQEventsWidget which is a child of 'parent', with the
 *  name 'name'.'
 */
CQEventsWidget::CQEventsWidget(QWidget* parent, const char* name)
    : CopasiWidget(parent, name)
{
  setupUi(this);

  //Create Source Data Model.
  mpEventDM = new CQEventDM(this);

  //Create the Proxy Model for sorting/filtering and set its properties.
  mpProxyModel = new CQSortFilterProxyModel();
  mpProxyModel->setDynamicSortFilter(true);
  mpProxyModel->setSortCaseSensitivity(Qt::CaseInsensitive);
  mpProxyModel->setFilterKeyColumn(-1);

  mpOrderDelegate = new CQSpinBoxDelegate(this);
  mpTblEvents->setItemDelegateForColumn(COL_ORDER_EVENTS, mpOrderDelegate);

  mpTblEvents->verticalHeader()->setResizeMode(QHeaderView::ResizeToContents);
  mpTblEvents->verticalHeader()->hide();
  mpTblEvents->sortByColumn(COL_ROW_NUMBER, Qt::AscendingOrder);

  // Connect the table widget
  connect(mpEventDM, SIGNAL(notifyGUI(ListViews::ObjectType, ListViews::Action, const std::string)),
          this, SLOT(protectedNotify(ListViews::ObjectType, ListViews::Action, const std::string)));
  connect(mpEventDM, SIGNAL(dataChanged(const QModelIndex&, const QModelIndex&)),
          this, SLOT(dataChanged(const QModelIndex&, const QModelIndex&)));
  connect(mpLEFilter, SIGNAL(textChanged(const QString &)),
          this, SLOT(slotFilterChanged()));
}

/*
 *  Destroys the object and frees any allocated resources
 */
CQEventsWidget::~CQEventsWidget()
{
  pdelete(mpOrderDelegate);
  pdelete(mpProxyModel);
  pdelete(mpEventDM);
  // no need to delete child widgets, Qt does it all for us
}

/*
 *  Sets the strings of the subwidgets using the current
 *  language.
 */
void CQEventsWidget::languageChange()
{
  retranslateUi(this);
}

void CQEventsWidget::slotBtnNewClicked()
{
  mpEventDM->insertRow();
  updateDeleteBtns();
}

void CQEventsWidget::slotBtnDeleteClicked()
{
  if (mpTblEvents->hasFocus())
    {deleteSelectedEvents();}

  updateDeleteBtns();
}

void CQEventsWidget::deleteSelectedEvents()
{
  const QItemSelectionModel * pSelectionModel = mpTblEvents->selectionModel();

  QModelIndexList mappedSelRows;
  size_t i, imax = mpEventDM->rowCount();

  for (i = 0; i < imax; i++)
    {
      if (pSelectionModel->isRowSelected(i, QModelIndex()))
        {
          mappedSelRows.append(mpProxyModel->mapToSource(mpProxyModel->index(i, 0)));
        }
    }

  if (mappedSelRows.empty())
    {return;}


  mpEventDM->removeRows(mappedSelRows);
}

void CQEventsWidget::slotBtnClearClicked()
{

  int ret = QMessageBox::question(this, tr("Confirm Delete"), "Delete all Events?",
                                  QMessageBox::Yes | QMessageBox::No, QMessageBox::No);

  if (ret == QMessageBox::Yes)
    {
      mpEventDM->clear();
    }

  updateDeleteBtns();
}

bool CQEventsWidget::update(ListViews::ObjectType C_UNUSED(objectType), ListViews::Action C_UNUSED(action), const std::string & C_UNUSED(key))
{
  if (!mIgnoreUpdates)
    {
      enterProtected();
    }

  return true;
}

bool CQEventsWidget::leave()
{
  return true;
}

bool CQEventsWidget::enterProtected()
{
  if (mpTblEvents->selectionModel() != NULL)
    {
      disconnect(mpTblEvents->selectionModel(), SIGNAL(selectionChanged(const QItemSelection&, const QItemSelection&)),
                 this, SLOT(slotSelectionChanged(const QItemSelection&, const QItemSelection&)));
    }

  mpProxyModel->setSourceModel(mpEventDM);
  //Set Model for the TableView
  mpTblEvents->setModel(NULL);
  mpTblEvents->setModel(mpProxyModel);
  connect(mpTblEvents->selectionModel(), SIGNAL(selectionChanged(const QItemSelection&, const QItemSelection&)),
          this, SLOT(slotSelectionChanged(const QItemSelection&, const QItemSelection&)));
  updateDeleteBtns();
  mpTblEvents->resizeColumnsToContents();

  return true;
}

void CQEventsWidget::updateDeleteBtns()
{
  bool selected = false;

  QModelIndexList selRows = mpTblEvents->selectionModel()->selectedRows();

  if (selRows.size() == 0)
    selected = false;
  else
    {
      if (selRows.size() == 1)
        {
          if (mpEventDM->isDefaultRow(mpProxyModel->mapToSource(selRows[0])))
            selected = false;
          else
            selected = true;
        }
      else
        selected = true;
    }

  mpBtnDelete->setEnabled(selected);

  if (mpProxyModel->rowCount() - 1)
    mpBtnClear->setEnabled(true);
  else
    mpBtnClear->setEnabled(false);
}

void CQEventsWidget::slotSelectionChanged(const QItemSelection& C_UNUSED(selected),
    const QItemSelection& C_UNUSED(deselected))
{
  updateDeleteBtns();
}

void CQEventsWidget::dataChanged(const QModelIndex& C_UNUSED(topLeft),
                                 const QModelIndex& C_UNUSED(bottomRight))
{
  mpTblEvents->resizeColumnsToContents();
  updateDeleteBtns();
}

void CQEventsWidget::slotDoubleClicked(const QModelIndex proxyIndex)
{
  QModelIndex index = mpProxyModel->mapToSource(proxyIndex);

  if (mpEventDM->isDefaultRow(index))
    return;

  assert(CCopasiRootContainer::getDatamodelList()->size() > 0);
  CCopasiDataModel* pDataModel = (*CCopasiRootContainer::getDatamodelList())[0];
  assert(pDataModel != NULL);
  CModel * pModel = pDataModel->getModel();

  if (pModel == NULL)
    return;

  std::string key = pModel->getEvents()[index.row()]->getKey();

  if (CCopasiRootContainer::getKeyFactory()->get(key))
    mpListView->switchToOtherWidget(0, key);
}

void CQEventsWidget::keyPressEvent(QKeyEvent* ev)
{
  if (ev->key() == Qt::Key_Delete)
    slotBtnDeleteClicked();
  else if (ev->key() == Qt::Key_C && ev->modifiers() & Qt::ControlModifier)
    {
      QModelIndexList selRows = mpTblEvents->selectionModel()->selectedRows(0);

      if (selRows.empty())
        {return;}

      QString str;
      QModelIndexList::const_iterator i;

      for (i = selRows.begin(); i != selRows.end(); ++i)
        {
          for (int x = 0; x < mpEventDM->columnCount(); ++x)
            {
              if (!mpTblEvents->isColumnHidden(x))
                {
                  if (!str.isEmpty())
                    str += "\t";

                  str += mpEventDM->index(mpProxyModel->mapToSource(*i).row(), x).data().toString();
                }
            }

          str += "\n";
        }

      QApplication::clipboard()->setText(str);
    }
}

void CQEventsWidget::slotFilterChanged()
{
  QRegExp regExp(mpLEFilter->text() + "|New Event", Qt::CaseInsensitive, QRegExp::RegExp);
  mpProxyModel->setFilterRegExp(regExp);
}
