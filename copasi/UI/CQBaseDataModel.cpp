// Begin CVS Header
//   $Source: /fs/turing/cvs/copasi_dev/copasi/UI/CQBaseDataModel.cpp,v $
//   $Revision: 1.7 $
//   $Name: Build-33 $
//   $Author: aekamal $
//   $Date: 2010/01/18 15:50:23 $
// End CVS Header

// Copyright (C) 2008 by Pedro Mendes, Virginia Tech Intellectual
// Properties, Inc., EML Research, gGmbH, University of Heidelberg,
// and The University of Manchester.
// All rights reserved.

#include "copasi.h"
#include "CQBaseDataModel.h"

CQBaseDataModel::CQBaseDataModel(QObject *parent)
    : QAbstractTableModel(parent)

{}

Qt::ItemFlags CQBaseDataModel::flags(const QModelIndex &index) const
{
  if (!index.isValid())
    return Qt::ItemIsEnabled;

  if (index.column() == COL_ROW_NUMBER)
    return QAbstractItemModel::flags(index);
  else
    return QAbstractItemModel::flags(index) | Qt::ItemIsEditable;
}

bool CQBaseDataModel::insertRow()
{
  return insertRows(rowCount() - 1, 1);
}

bool CQBaseDataModel::removeRow(int position)
{
  if (position >= 0 && (position < rowCount() - 1) && !isDefaultRow(index(position, 0)))
    return removeRows(position, 1);
  else
    return false;
}

bool CQBaseDataModel::clear()
{
  return removeRows(0, rowCount() - 1);
}

bool CQBaseDataModel::isDefaultRow(const QModelIndex& i) const
{
  //Index has to be from this model and should be valid.
  assert((i.model() == this) && i.isValid());

  return (i.row() == rowCount() - 1);
}

QString CQBaseDataModel::createNewName(const QString name, const int nameCol)
{
  QString nname = name;
  unsigned C_INT32 j, jmax = rowCount();

  for (unsigned C_INT32 i = 1;; ++i)
    {
      nname = name + "_" + QString::number(i);

      for (j = 0; j < jmax; ++j)
        if (index(j, nameCol).data() == nname) break;

      if (j == jmax) break;
    }

  return nname;
}
