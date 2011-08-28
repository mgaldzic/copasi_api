// Begin CVS Header
//   $Source: /fs/turing/cvs/copasi_dev/copasi/UI/CQComboDelegate.cpp,v $
//   $Revision: 1.1 $
//   $Name: Build-33 $
//   $Author: aekamal $
//   $Date: 2009/05/04 15:24:00 $
// End CVS Header

// Copyright (C) 2008 by Pedro Mendes, Virginia Tech Intellectual
// Properties, Inc., EML Research, gGmbH, University of Heidelberg,
// and The University of Manchester.
// All rights reserved.

#include <QComboBox>

#include "../copasi.h"
#include "CQComboDelegate.h"

CQComboDelegate::CQComboDelegate(const QStringList *pComboItems, QObject *parent)
    : QItemDelegate(parent)
{
  mpComboItems = pComboItems;
}

QWidget *CQComboDelegate::createEditor(QWidget *parent,
                                       const QStyleOptionViewItem & C_UNUSED(option),
                                       const QModelIndex & C_UNUSED(index)) const
{
  QComboBox *editor = new QComboBox(parent);
  editor->addItems(*mpComboItems);
  return editor;
}

void CQComboDelegate::setEditorData(QWidget *editor,
                                    const QModelIndex &index) const
{
  QString value = index.model()->data(index, Qt::EditRole).toString();
  QComboBox *comboBox = static_cast<QComboBox*>(editor);
  comboBox->setCurrentItem(comboBox->findText(value));
}

void CQComboDelegate::setModelData(QWidget *editor, QAbstractItemModel *model,
                                   const QModelIndex &index) const
{
  QComboBox *comboBox = static_cast<QComboBox*>(editor);
  QVariant value(comboBox->currentText());
  model->setData(index, value, Qt::EditRole);
}

void CQComboDelegate::updateEditorGeometry(QWidget *editor,
    const QStyleOptionViewItem &option, const QModelIndex & C_UNUSED(index)) const
{
  editor->setGeometry(option.rect);
}

CQIndexComboDelegate::CQIndexComboDelegate(const QStringList *pComboItems, QObject *parent)
    : CQComboDelegate(pComboItems, parent)
{
}

void CQIndexComboDelegate::setModelData(QWidget *editor, QAbstractItemModel *model,
                                        const QModelIndex &index) const
{
  QComboBox *comboBox = static_cast<QComboBox*>(editor);
  QVariant value(comboBox->currentIndex());
  model->setData(index, value, Qt::EditRole);
}
