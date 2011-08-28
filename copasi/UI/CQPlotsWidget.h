// Begin CVS Header
//   $Source: /fs/turing/cvs/copasi_dev/copasi/UI/CQPlotsWidget.h,v $
//   $Revision: 1.3 $
//   $Name: Build-33 $
//   $Author: aekamal $
//   $Date: 2010/09/03 21:06:11 $
// End CVS Header

// Copyright (C) 2010 by Pedro Mendes, Virginia Tech Intellectual
// Properties, Inc., University of Heidelberg, and The University
// of Manchester.
// All rights reserved.

// Copyright (C) 2008 by Pedro Mendes, Virginia Tech Intellectual
// Properties, Inc., EML Research, gGmbH, University of Heidelberg,
// and The University of Manchester.
// All rights reserved.

#ifndef CQPlotsWidget_h
#define CQPlotsWidget_h

#include "CQSortFilterProxyModel.h"
#include "ui_CQPlotsWidget.h"
#include "CQPlotDM.h"

class CQPlotsWidget : public CopasiWidget, public Ui::CQPlotsWidget
{
  Q_OBJECT

public:
  CQPlotsWidget(QWidget* parent = 0, const char* name = 0);
  ~CQPlotsWidget();

  virtual bool update(ListViews::ObjectType objectType, ListViews::Action action, const std::string & key);
  virtual bool leave();
  virtual void setFramework(int framework);

private:
  CQPlotDM* mpPlotDM;
  CQSortFilterProxyModel *mpProxyModel;
  void deleteSelectedPlots();
  void updateDeleteBtns();

protected:
  virtual void keyPressEvent(QKeyEvent* ev);
  virtual bool enterProtected();

protected slots:
  virtual void languageChange();
  virtual void slotBtnNewClicked();
  virtual void slotBtnDeleteClicked();
  virtual void slotBtnClearClicked();
  virtual void slotSelectionChanged(const QItemSelection& selected,
                                    const QItemSelection& deselected);
  virtual void slotDoubleClicked(const QModelIndex proxyIndex);
  virtual void dataChanged(const QModelIndex& topLeft,
                           const QModelIndex& bottomRight);
  virtual void slotFilterChanged();
};

#endif // CQPlotsWidget_h
