// Begin CVS Header
//   $Source: /fs/turing/cvs/copasi_dev/copasi/UI/CQFunctionDM.h,v $
//   $Revision: 1.1.4.1 $
//   $Name: Build-33 $
//   $Author: aekamal $
//   $Date: 2010/09/27 13:44:56 $
// End CVS Header

// Copyright (C) 2010 by Pedro Mendes, Virginia Tech Intellectual
// Properties, Inc., University of Heidelberg, and The University
// of Manchester.
// All rights reserved.

// Copyright (C) 2008 by Pedro Mendes, Virginia Tech Intellectual
// Properties, Inc., EML Research, gGmbH, University of Heidelberg,
// and The University of Manchester.
// All rights reserved.

#ifndef CQFunctionDM_H
#define CQFunctionDM_H

#include "model/CModel.h"
#include "listviews.h"
#include "CQBaseDataModel.h"

#define COL_NAME_FUNCTIONS       1
#define COL_TYPE_FUNCTIONS       2
#define COL_MATH_DESC_FUNCTIONS  3
#define COL_SBML_ID_FUNCTIONS    4
#define TOTAL_COLS_FUNCTIONS     5

class CQFunctionDM : public CQBaseDataModel
{
  Q_OBJECT

public:
  CQFunctionDM(QObject *parent = 0);
  int rowCount(const QModelIndex &parent = QModelIndex()) const;
  int columnCount(const QModelIndex &parent = QModelIndex()) const;
  Qt::ItemFlags flags(const QModelIndex &index) const;
  QVariant data(const QModelIndex &index, int role) const;
  QVariant headerData(int section, Qt::Orientation orientation,
                      int role = Qt::DisplayRole) const;
  bool setData(const QModelIndex &index, const QVariant &value,
               int role = Qt::EditRole);
  bool removeRows(QModelIndexList rows, const QModelIndex &index = QModelIndex());

protected:
  bool insertRows(int position, int rows, const QModelIndex &index = QModelIndex());
  bool removeRows(int position, int rows, const QModelIndex &index = QModelIndex());

private:
  bool isFunctionReadOnly(const QModelIndex &index) const;
};

#endif //CQFunctionDM_H
