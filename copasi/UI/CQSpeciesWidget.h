// Begin CVS Header
//   $Source: /fs/turing/cvs/copasi_dev/copasi/UI/CQSpeciesWidget.h,v $
//   $Revision: 1.7 $
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

#ifndef CQSpeciesWidget_h
#define CQSpeciesWidget_h

#include "CQSortFilterProxyModel.h"
#include "CQComboDelegate.h"
#include "ui_CQSpeciesWidget.h"
#include "CQSpecieDM.h"

class CQSpeciesWidget : public CopasiWidget, public Ui::CQSpeciesWidget
{
  Q_OBJECT

public:
  CQSpeciesWidget(QWidget* parent = 0, const char* name = 0);
  ~CQSpeciesWidget();

  virtual bool update(ListViews::ObjectType objectType, ListViews::Action action, const std::string & key);
  virtual bool leave();
  virtual void setFramework(int framework);

private:
  CQSpecieDM* mpSpecieDM;
  CQSortFilterProxyModel *mpProxyModel;
  QStringList mCompartments;
  CQComboDelegate* mpCompartmentDelegate;
  CQIndexComboDelegate* mpTypeDelegate;
  void deleteSelectedSpecies();
  void updateDeleteBtns();

protected:
  virtual bool enterProtected();
  virtual void keyPressEvent(QKeyEvent* ev);
  void refreshCompartments();

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

#endif // CQSpeciesWidget_h
