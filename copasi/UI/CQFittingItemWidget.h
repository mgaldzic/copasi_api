// Begin CVS Header
//   $Source: /fs/turing/cvs/copasi_dev/copasi/UI/CQFittingItemWidget.h,v $
//   $Revision: 1.26 $
//   $Name: Build-33 $
//   $Author: shoops $
//   $Date: 2009/07/20 19:31:34 $
// End CVS Header

// Copyright (C) 2008 by Pedro Mendes, Virginia Tech Intellectual
// Properties, Inc., EML Research, gGmbH, University of Heidelberg,
// and The University of Manchester.
// All rights reserved.

#ifndef CQFITTINGITEMWIDGET_H
#define CQFITTINGITEMWIDGET_H

#include <qvariant.h>

#include "ui_CQFittingItemWidget.h"

/*
//Added by qt3to4:
#include <Q3HBoxLayout>
#include <Q3GridLayout>
#include <Q3VBoxLayout>
 */

#include <QHBoxLayout>
#include <QGridLayout>
#include <QVBoxLayout>

#include <QLabel>
#include <QPixmap>

class CCopasiDataModel;
class CCopasiSelectionDialog;
class COptItem;
class CCopasiObject;
class CCopasiObjectName;
class CQValidatorBound;
class CQValidatorNotEmpty;
class QColor;
class CExperimentSet;
class CCopasiParameterGroup;
class CCrossValidationSet;

#ifndef COPASI_CROSSVALIDATION
# define pCrossValidationMap
#endif // not COPASI_CROSSVALIDATION

class CQFittingItemWidget : public QWidget, public Ui::CQFittingItemWidget
{
  Q_OBJECT

public:
  CQFittingItemWidget(QWidget* parent, const char* name = 0, Qt::WindowFlags fl = 0);
  ~CQFittingItemWidget();

  enum ItemType {OPT_ITEM = 0, OPT_CONSTRAINT, FIT_ITEM, FIT_CONSTRAINT};

  virtual bool load(CCopasiDataModel * pDataModel,
                    CCopasiParameterGroup * pItems,
                    const std::map<std::string, std::string> * pExperimentMap,
                    const std::map<std::string, std::string> * pCrossValidationMap);
  virtual bool save(const std::map<std::string, std::string> * pExperimentMap,
                    const std::map<std::string, std::string> * pCrossValidationMap);
  void setItemType(const ItemType & type);
  void setExperimentSet(const CExperimentSet * & pExperimentSet);
  void setCrossValidationSet(const CCrossValidationSet * & pCrossValidationSet);

signals:
  void numberChanged(int);

protected:
  const CCopasiDataModel * mpDataModel;
  const CCrossValidationSet **mppCrossValidationSet;
  std::set< unsigned int > mSelection;
  unsigned int mCurrentRow;
  std::vector< COptItem * > * mpItemsCopy;
  ItemType mItemType;
  QColor mChangedColor;
  QColor mSavedColor;
  bool mUpperInfChanged;
  bool mLowerInfChanged;
  CQValidatorBound * mpUpperValidator;
  CQValidatorBound * mpLowerValidator;
  CQValidatorNotEmpty * mpObjectValidator;
  const CCopasiObject* mpUpperObject;
  const CCopasiObject* mpLowerObject;
  CCopasiObjectName* mpObjectCN;
  const CExperimentSet ** mppExperimentSet;
  CCopasiParameterGroup * mpItems;

protected slots:
  virtual void languageChange();

private:
  void init();
  void destroy();
  void setTableText(const int & row, const COptItem * pItem);
  unsigned int currentRow();
  void loadSelection();
  void saveSelection();
  void selectRow(const unsigned int & row);
  void setItemSelection(const std::set<unsigned int> & selection);

private slots:
  void slotCheckLowerInf(bool checked);
  void slotCheckUpperInf(bool checked);
  void slotLowerEdit();
  void slotUpperEdit();
  void slotParamEdit();
  void slotExperiments();
  void slotExperimentChanged();
  void slotDelete();
  void slotCopy();
  void slotUp();
  void slotDown();
  void slotDuplicatePerExperiment();
  void slotNew();
  void slotSelectionChanged();
  void slotLowerLostFocus();
  void slotUpperLostFocus();
  void slotReset();
  void slotStartLostFocus();
  void slotCrossValidations();
  void slotCrossValidationChanged();
  void slotCheckAllCrossValidations(bool checked);
  void slotCheckAllExperiments(bool checked);

protected:
  enum IconID
  {
    image0_ID,
    image1_ID,
    image2_ID,
    image3_ID,
    image4_ID,
    image5_ID,
    image6_ID,
    unknown_ID
  };
  static QPixmap qt_get_icon(IconID id)
  {
    static const unsigned char image0_data[] =
    {
      0x89, 0x50, 0x4e, 0x47, 0x0d, 0x0a, 0x1a, 0x0a, 0x00, 0x00, 0x00, 0x0d,
      0x49, 0x48, 0x44, 0x52, 0x00, 0x00, 0x00, 0x10, 0x00, 0x00, 0x00, 0x10,
      0x08, 0x06, 0x00, 0x00, 0x00, 0x1f, 0xf3, 0xff, 0x61, 0x00, 0x00, 0x02,
      0x7c, 0x49, 0x44, 0x41, 0x54, 0x78, 0x9c, 0x7d, 0x93, 0x5d, 0x68, 0x8e,
      0x61, 0x18, 0xc7, 0x7f, 0xf7, 0xfd, 0x3c, 0xf7, 0xf3, 0xce, 0x5e, 0x66,
      0x63, 0xbe, 0x45, 0xbe, 0xb2, 0x48, 0xc9, 0x81, 0x14, 0x49, 0x42, 0x0e,
      0x70, 0xb4, 0x13, 0x07, 0x8e, 0x28, 0x25, 0xca, 0x28, 0x49, 0x51, 0x1c,
      0xcc, 0x47, 0x39, 0x21, 0x1f, 0x71, 0x20, 0x39, 0x51, 0x6a, 0x25, 0xc5,
      0x94, 0x1d, 0x88, 0xb1, 0x96, 0xcd, 0x47, 0xbe, 0xc6, 0x36, 0x9f, 0xb3,
      0x79, 0x37, 0xdb, 0xfb, 0xb1, 0x87, 0x77, 0xef, 0xf3, 0x3c, 0xf7, 0xe5,
      0xc0, 0xfb, 0xda, 0x8b, 0xe5, 0xaa, 0xab, 0xae, 0x83, 0xeb, 0xfa, 0xdf,
      0xff, 0x8f, 0x6e, 0xd5, 0xd8, 0xfc, 0xb2, 0xaf, 0xbc, 0x6c, 0x74, 0xcc,
      0x8a, 0x50, 0x28, 0x2b, 0x90, 0xcd, 0x0e, 0x11, 0xe4, 0x02, 0xfa, 0xfa,
      0x53, 0x3a, 0x0c, 0xec, 0x40, 0x32, 0xe5, 0xaf, 0xdf, 0xba, 0x65, 0xdd,
      0x73, 0xfe, 0xae, 0x47, 0x4f, 0xda, 0x7f, 0x58, 0x6b, 0xa5, 0xb8, 0xa3,
      0x28, 0x12, 0xff, 0x47, 0x20, 0xe9, 0x4c, 0x20, 0x8f, 0x9f, 0x75, 0xca,
      0xdb, 0x8e, 0x6e, 0xa9, 0x6f, 0x68, 0xed, 0xbd, 0x78, 0xe5, 0xd6, 0x2a,
      0x11, 0xa1, 0xb8, 0xf5, 0x3f, 0x88, 0x80, 0x52, 0x0a, 0xad, 0x41, 0x69,
      0x8b, 0x76, 0x34, 0x93, 0xa7, 0x54, 0xb2, 0x6a, 0xf9, 0xa2, 0xca, 0xaa,
      0x79, 0x33, 0x6f, 0x5f, 0xad, 0xbb, 0x57, 0x5d, 0xbc, 0x3b, 0x22, 0x80,
      0x00, 0x0a, 0xf0, 0x3c, 0x97, 0x77, 0xef, 0x3f, 0xd1, 0xd4, 0xd4, 0x4a,
      0xcb, 0xd3, 0x97, 0x84, 0x41, 0xd6, 0xcb, 0xa4, 0x93, 0xd7, 0x0e, 0x1e,
      0x39, 0x5b, 0xf3, 0x5f, 0x00, 0x00, 0xa5, 0x15, 0x5a, 0xc1, 0x84, 0xca,
      0x71, 0xe4, 0xac, 0xe2, 0x6b, 0x4f, 0x8a, 0xf6, 0x8e, 0x4f, 0x78, 0x25,
      0x1e, 0x99, 0xc1, 0x64, 0x6d, 0x61, 0xcf, 0x1d, 0xf1, 0x18, 0xd0, 0x4a,
      0x11, 0x06, 0x21, 0x8b, 0x16, 0x57, 0x91, 0xcb, 0x06, 0xf4, 0x65, 0x86,
      0xe8, 0xef, 0xee, 0xc6, 0x18, 0x8f, 0xe6, 0xa6, 0xc6, 0x11, 0x24, 0xd8,
      0x08, 0x5e, 0x35, 0x43, 0x5b, 0x0b, 0x16, 0x08, 0x05, 0x42, 0x63, 0x18,
      0xec, 0xec, 0x44, 0x9f, 0xd8, 0x83, 0x69, 0xa8, 0xc3, 0x94, 0xc6, 0xf1,
      0x3c, 0x83, 0xe3, 0x38, 0xbf, 0xcf, 0x86, 0x19, 0x7c, 0x7e, 0x03, 0x0f,
      0xeb, 0xc1, 0x4f, 0x22, 0x1f, 0x3b, 0xf1, 0x57, 0x54, 0x23, 0x1f, 0xda,
      0xa8, 0x38, 0x75, 0x80, 0x58, 0x43, 0x3d, 0xa3, 0x66, 0x3c, 0x20, 0x58,
      0xb8, 0x14, 0x3d, 0x71, 0x1a, 0x4a, 0x31, 0x02, 0xc0, 0x94, 0x39, 0x30,
      0x73, 0x3e, 0xdc, 0xa9, 0xc3, 0x69, 0xbb, 0x44, 0xfc, 0x51, 0x23, 0xba,
      0xbd, 0x1d, 0x73, 0xef, 0x2e, 0xc1, 0xd4, 0xd9, 0x0c, 0xec, 0x3c, 0x8c,
      0x53, 0x39, 0x09, 0xa3, 0x6c, 0x5e, 0x64, 0x01, 0x20, 0x3f, 0x67, 0xb5,
      0x8b, 0x5a, 0x59, 0x8d, 0xfa, 0x3e, 0x84, 0xb9, 0x7a, 0x9e, 0xd8, 0x8b,
      0x1b, 0x90, 0x18, 0x20, 0x58, 0xb0, 0x98, 0xfe, 0x43, 0x67, 0x08, 0x66,
      0x55, 0xe1, 0x0c, 0x0e, 0xe2, 0xb8, 0x2e, 0xaa, 0x88, 0x82, 0x5b, 0x30,
      0x21, 0x8c, 0x40, 0x19, 0x07, 0x57, 0x97, 0xa0, 0xba, 0x12, 0xf0, 0xa5,
      0x17, 0x32, 0x59, 0xe8, 0xfa, 0x86, 0xb2, 0x0a, 0xe5, 0x19, 0x1c, 0x05,
      0x8e, 0x2a, 0x04, 0x9d, 0x37, 0xd1, 0x16, 0x9c, 0x0f, 0x72, 0x78, 0x57,
      0x4e, 0x13, 0x3b, 0xbe, 0x1f, 0x7a, 0x53, 0x44, 0x4b, 0x96, 0x61, 0xe3,
      0xe3, 0x31, 0x2d, 0xaf, 0x19, 0xbf, 0x61, 0x2d, 0xb1, 0x9b, 0x37, 0x60,
      0x54, 0x69, 0xfe, 0xf5, 0x61, 0x06, 0xbf, 0x53, 0x70, 0xef, 0xd7, 0x63,
      0xce, 0x1c, 0x83, 0x9e, 0x3e, 0x72, 0x9b, 0xb7, 0x92, 0xb9, 0x7c, 0x1b,
      0xff, 0xe4, 0x05, 0xa4, 0x6c, 0x2c, 0x4e, 0x22, 0x41, 0xf9, 0xde, 0x5d,
      0xe8, 0x44, 0x0f, 0x62, 0xbc, 0x3f, 0x22, 0x1f, 0x96, 0xb0, 0x74, 0x35,
      0x6a, 0xf5, 0x26, 0x64, 0xfa, 0x5c, 0x72, 0xdb, 0x6a, 0x10, 0x11, 0x82,
      0x8d, 0x1b, 0x49, 0x5e, 0xbb, 0x4e, 0xe9, 0xbe, 0xdd, 0xf8, 0xdb, 0x77,
      0x60, 0xc7, 0x8c, 0xc5, 0x09, 0xc2, 0x3f, 0x24, 0xb8, 0x36, 0x3f, 0xc7,
      0x2b, 0xca, 0xe0, 0xe8, 0x39, 0x04, 0x30, 0x40, 0x14, 0x81, 0x08, 0x44,
      0x6b, 0x56, 0x92, 0x7b, 0xdc, 0x8a, 0x0a, 0x84, 0xb8, 0x9f, 0xc3, 0x75,
      0x0d, 0xc5, 0x39, 0xba, 0xa9, 0x74, 0xda, 0xaf, 0x3d, 0x7e, 0x3a, 0xd4,
      0xfa, 0x57, 0xa2, 0x82, 0x45, 0xa1, 0x11, 0xb1, 0x88, 0x80, 0xe4, 0x7f,
      0x86, 0x88, 0x80, 0x08, 0x4a, 0x69, 0xfc, 0x74, 0x66, 0xa8, 0x00, 0xf0,
      0x13, 0xf0, 0x32, 0x22, 0x79, 0x63, 0x5a, 0xb7, 0x6f, 0x00, 0x00, 0x00,
      0x00, 0x49, 0x45, 0x4e, 0x44, 0xae, 0x42, 0x60, 0x82
    };

    static const unsigned char image1_data[] =
    {
      0x89, 0x50, 0x4e, 0x47, 0x0d, 0x0a, 0x1a, 0x0a, 0x00, 0x00, 0x00, 0x0d,
      0x49, 0x48, 0x44, 0x52, 0x00, 0x00, 0x00, 0x10, 0x00, 0x00, 0x00, 0x10,
      0x08, 0x06, 0x00, 0x00, 0x00, 0x1f, 0xf3, 0xff, 0x61, 0x00, 0x00, 0x02,
      0xf3, 0x49, 0x44, 0x41, 0x54, 0x78, 0x9c, 0x65, 0x93, 0x4d, 0x68, 0x5c,
      0x65, 0x14, 0x86, 0x9f, 0xef, 0xbb, 0x73, 0xef, 0xcc, 0xdc, 0xcc, 0x4f,
      0x32, 0x49, 0x67, 0xd2, 0x80, 0x69, 0x9d, 0x0e, 0x8d, 0x96, 0xb6, 0x68,
      0xa1, 0x05, 0x63, 0x29, 0x51, 0xbb, 0x10, 0x09, 0x0a, 0x8a, 0xb3, 0x28,
      0x85, 0x56, 0x37, 0x51, 0x82, 0xe2, 0xce, 0x4d, 0x17, 0x01, 0x71, 0xe1,
      0x76, 0x28, 0x48, 0x40, 0x2b, 0x48, 0x28, 0x9a, 0x6d, 0x17, 0x52, 0x21,
      0x50, 0xdb, 0xba, 0xa8, 0x8d, 0xa6, 0x50, 0x4c, 0x6b, 0xd2, 0xe0, 0x24,
      0x8c, 0xe9, 0xfc, 0x64, 0x3a, 0x73, 0x6f, 0x93, 0x69, 0xef, 0xcc, 0xfd,
      0xee, 0xe7, 0x42, 0x3b, 0x9a, 0xf4, 0xc0, 0xd9, 0x9d, 0xf7, 0xe1, 0x3d,
      0xbc, 0xe7, 0x08, 0xad, 0x35, 0x3b, 0xeb, 0x95, 0xf7, 0x8a, 0x11, 0x37,
      0xea, 0x1e, 0x05, 0x3f, 0xa7, 0x75, 0x20, 0x10, 0x7a, 0x45, 0xaa, 0xca,
      0xcd, 0xf9, 0xe9, 0xf1, 0xd6, 0xce, 0x59, 0xf1, 0x7f, 0xc0, 0xc1, 0xc9,
      0xc5, 0x58, 0x44, 0xab, 0x73, 0xda, 0x30, 0x3f, 0x90, 0x56, 0x22, 0xf9,
      0x4c, 0x5f, 0x07, 0xdb, 0xf4, 0xf9, 0xa3, 0x1a, 0x21, 0xf0, 0x1c, 0x47,
      0xfa, 0x5b, 0xdf, 0x6b, 0xed, 0x7c, 0x3e, 0x3f, 0x7d, 0x72, 0xed, 0x29,
      0xc0, 0xb1, 0xc9, 0xc5, 0x41, 0x45, 0x70, 0x25, 0x93, 0xee, 0x1f, 0x39,
      0xf3, 0x6a, 0x92, 0xe3, 0xcf, 0x1b, 0x24, 0x6c, 0x03, 0x00, 0xaf, 0xa3,
      0xf9, 0xed, 0x5e, 0x8b, 0x0b, 0x73, 0x9b, 0x2c, 0x17, 0xab, 0x1b, 0xb6,
      0xae, 0xbc, 0x7b, 0xed, 0xcb, 0xd7, 0xaf, 0x74, 0x01, 0x07, 0x27, 0x17,
      0x63, 0x61, 0xd4, 0xfc, 0x91, 0x03, 0x43, 0x23, 0x9f, 0x8c, 0x47, 0x88,
      0x47, 0x43, 0x00, 0x5c, 0xbd, 0xed, 0x76, 0xdd, 0x9d, 0x38, 0x94, 0xa0,
      0xa3, 0x02, 0xbe, 0xbe, 0xdc, 0xe0, 0xf2, 0x2f, 0x95, 0x7a, 0x8f, 0xfa,
      0xfd, 0xd0, 0xd5, 0xe9, 0xd3, 0xf7, 0x43, 0x00, 0x16, 0x9d, 0xa9, 0x4c,
      0x3a, 0x3d, 0x72, 0x76, 0xcc, 0xc4, 0x6b, 0x83, 0xd7, 0xf6, 0x59, 0x5c,
      0x6b, 0xf1, 0xd6, 0x4b, 0x29, 0xe2, 0xff, 0xba, 0x98, 0x99, 0xab, 0x71,
      0x60, 0xd8, 0xe6, 0xed, 0xd1, 0x04, 0xab, 0xb5, 0xa0, 0x7f, 0x69, 0x25,
      0x3b, 0x95, 0xcf, 0xe7, 0x3f, 0x32, 0x7e, 0x5a, 0x1d, 0x8b, 0xb4, 0xad,
      0xc4, 0xc5, 0x77, 0x8e, 0x27, 0x23, 0x09, 0x1b, 0x5a, 0x9e, 0x62, 0xb5,
      0xe2, 0xf1, 0xc6, 0xb1, 0x5e, 0x7a, 0x63, 0xa1, 0xae, 0x83, 0xc3, 0xd9,
      0x1e, 0xe6, 0x6e, 0x39, 0x00, 0x24, 0xa3, 0x3e, 0xbf, 0xae, 0xe8, 0x11,
      0x53, 0x74, 0x66, 0x42, 0x6e, 0x34, 0x71, 0xd4, 0x8a, 0xf4, 0x24, 0xfb,
      0xa2, 0x8a, 0x5a, 0x33, 0x00, 0xe0, 0xe4, 0x91, 0x24, 0xa9, 0xb8, 0xf9,
      0x54, 0x3a, 0xa7, 0x5f, 0xdb, 0xc5, 0xa7, 0x17, 0xdb, 0xfc, 0xbc, 0xdc,
      0x87, 0x0c, 0xeb, 0x98, 0xe3, 0xe5, 0x46, 0x25, 0x90, 0x4b, 0xc5, 0x24,
      0x35, 0xd7, 0x67, 0xc3, 0xf1, 0x19, 0x3b, 0x9c, 0x60, 0xb0, 0xcf, 0x02,
      0xa0, 0x54, 0xf3, 0xba, 0xe2, 0xe5, 0xf5, 0x0e, 0x00, 0x5f, 0x9c, 0xb2,
      0x78, 0x71, 0x0f, 0x48, 0x33, 0x82, 0x94, 0x46, 0x4e, 0x6a, 0x1d, 0x88,
      0xb6, 0x32, 0x58, 0xab, 0x1b, 0xfc, 0xd5, 0x34, 0x68, 0x2b, 0x01, 0xc0,
      0xb9, 0xef, 0x3c, 0xae, 0x2f, 0xc9, 0x2e, 0xe0, 0xfd, 0xaf, 0x4c, 0x16,
      0x8a, 0x4f, 0xb2, 0x03, 0x01, 0x08, 0xa1, 0x65, 0x08, 0xa1, 0x57, 0x9c,
      0x96, 0xe6, 0x6e, 0x7d, 0x00, 0x80, 0x0f, 0xbf, 0x01, 0x8d, 0x22, 0xd0,
      0x26, 0xd9, 0xdd, 0xff, 0x01, 0x00, 0x3e, 0xfe, 0x16, 0x84, 0x00, 0x43,
      0x80, 0xf2, 0x3d, 0x2c, 0xd5, 0x2a, 0x49, 0xa9, 0x2a, 0x37, 0x95, 0xe7,
      0x38, 0x52, 0x80, 0x94, 0xff, 0xb4, 0x21, 0x41, 0xe9, 0xed, 0xe2, 0x90,
      0xd4, 0x44, 0x4c, 0x45, 0xd4, 0x54, 0x44, 0x42, 0x8a, 0xe0, 0x71, 0x63,
      0xab, 0xdf, 0xbf, 0xb3, 0x20, 0xe7, 0xa7, 0xc7, 0x5b, 0xa8, 0xcd, 0x0b,
      0xd2, 0x2b, 0x33, 0x9c, 0x70, 0xbb, 0xbd, 0x3f, 0xd5, 0x60, 0x97, 0xfd,
      0xa8, 0x0b, 0xd8, 0xdf, 0xdf, 0x24, 0xdb, 0xeb, 0xb2, 0x37, 0xe9, 0x22,
      0xbd, 0x32, 0xa6, 0x5f, 0xbd, 0x94, 0x0e, 0x97, 0xd6, 0x43, 0x00, 0xb6,
      0xa8, 0x7d, 0x56, 0xad, 0x18, 0x6f, 0x66, 0x07, 0x06, 0xf6, 0xa5, 0x92,
      0x56, 0x57, 0xb4, 0xd1, 0x0c, 0x80, 0x28, 0x7f, 0x96, 0x3d, 0x06, 0xe3,
      0x8f, 0x01, 0x78, 0xe0, 0xb4, 0xa9, 0x96, 0x6b, 0xa5, 0x21, 0x75, 0xe3,
      0x7c, 0xc3, 0x6d, 0xd4, 0xbb, 0xa7, 0xfc, 0xf2, 0xc4, 0x0f, 0xfb, 0x3a,
      0xa1, 0xf4, 0x8f, 0x2f, 0xe4, 0x92, 0xd9, 0xe7, 0x86, 0x63, 0x08, 0x21,
      0xb6, 0xad, 0xa0, 0xb5, 0xe6, 0xee, 0xda, 0x26, 0xb7, 0x96, 0x1a, 0xa5,
      0x8c, 0x7f, 0x6d, 0xe2, 0xd9, 0x68, 0xf1, 0x7a, 0xa1, 0x50, 0x70, 0xb7,
      0x3d, 0xd3, 0x89, 0x89, 0x99, 0xdd, 0x8f, 0x44, 0x76, 0xca, 0xb6, 0xed,
      0x53, 0x7b, 0x33, 0x56, 0xac, 0x37, 0x1e, 0x06, 0xa0, 0xf9, 0xd0, 0xa3,
      0x58, 0x6e, 0x6f, 0x79, 0xad, 0x07, 0x97, 0x32, 0xea, 0xc6, 0xf9, 0x3d,
      0xf6, 0xfd, 0xdb, 0x85, 0x42, 0xc1, 0x85, 0x1d, 0xdf, 0x08, 0x90, 0xcf,
      0xe7, 0xad, 0x87, 0x3d, 0xa3, 0x43, 0x8e, 0x91, 0x1b, 0x95, 0xd2, 0xc8,
      0x09, 0xa1, 0xa5, 0x50, 0xad, 0x52, 0xbf, 0x7f, 0x67, 0x21, 0x1d, 0x2e,
      0xad, 0x37, 0x1a, 0x8d, 0xfa, 0xec, 0xec, 0x6c, 0xfb, 0xc9, 0xfc, 0xdf,
      0xaa, 0x18, 0x4a, 0x98, 0xa0, 0x28, 0xcd, 0xf7, 0x00, 0x00, 0x00, 0x00,
      0x49, 0x45, 0x4e, 0x44, 0xae, 0x42, 0x60, 0x82
    };

    static const unsigned char image2_data[] =
    {
      0x89, 0x50, 0x4e, 0x47, 0x0d, 0x0a, 0x1a, 0x0a, 0x00, 0x00, 0x00, 0x0d,
      0x49, 0x48, 0x44, 0x52, 0x00, 0x00, 0x00, 0x10, 0x00, 0x00, 0x00, 0x10,
      0x08, 0x06, 0x00, 0x00, 0x00, 0x1f, 0xf3, 0xff, 0x61, 0x00, 0x00, 0x02,
      0x5c, 0x49, 0x44, 0x41, 0x54, 0x78, 0x9c, 0xa5, 0x93, 0x4d, 0x4f, 0x13,
      0x71, 0x10, 0x87, 0x9f, 0xff, 0xee, 0xd6, 0xd6, 0xb6, 0xdb, 0x52, 0x5a,
      0x40, 0x81, 0x16, 0x85, 0x5a, 0x5e, 0x04, 0xa2, 0xa9, 0x26, 0x26, 0x98,
      0x88, 0x89, 0x07, 0xaf, 0xc6, 0x0f, 0xe8, 0x47, 0xc0, 0x98, 0x78, 0x30,
      0x6a, 0x0c, 0x97, 0x62, 0xc4, 0x82, 0xc8, 0x6b, 0x63, 0x43, 0xaa, 0xa0,
      0xad, 0x5d, 0xe8, 0xf2, 0x62, 0x29, 0xdd, 0x6d, 0x3b, 0x1e, 0x10, 0x22,
      0x41, 0xf4, 0xe0, 0x9c, 0x26, 0x33, 0xf3, 0xfc, 0x0e, 0xbf, 0x99, 0x51,
      0x99, 0x4c, 0x86, 0xff, 0x09, 0xe3, 0x6f, 0x4d, 0x9f, 0xaf, 0x2c, 0xae,
      0xab, 0x61, 0xdb, 0x65, 0xb6, 0xb6, 0x2a, 0x78, 0xbd, 0x06, 0x3d, 0x3d,
      0x77, 0xd4, 0x3f, 0x05, 0x4c, 0xd3, 0x96, 0x4a, 0xe5, 0x80, 0x6c, 0xb6,
      0x40, 0xa1, 0x50, 0xc7, 0xeb, 0x15, 0x62, 0x31, 0x61, 0x60, 0x20, 0x4a,
      0x38, 0x5c, 0x91, 0xdd, 0xdd, 0x76, 0x75, 0xae, 0x80, 0xcf, 0x57, 0x97,
      0xf9, 0xf9, 0x0d, 0xe6, 0xe6, 0x2c, 0x82, 0x41, 0xc5, 0xc4, 0x44, 0x37,
      0x89, 0x44, 0x02, 0x91, 0x1f, 0x18, 0x86, 0xf0, 0x3b, 0x7c, 0x46, 0xc0,
      0xb6, 0xbf, 0xca, 0xcc, 0xcc, 0x67, 0x2c, 0xeb, 0x80, 0x64, 0xd2, 0xc3,
      0xe8, 0x68, 0x2f, 0xe1, 0xf0, 0x80, 0x72, 0x1c, 0x00, 0x3f, 0xae, 0xfb,
      0x17, 0x0f, 0xf2, 0xf9, 0x9c, 0x64, 0xe7, 0x4a, 0x04, 0xfc, 0x1a, 0x93,
      0x93, 0x7d, 0xf4, 0xf5, 0x5d, 0x53, 0x67, 0xc7, 0xcf, 0x11, 0xc8, 0xe7,
      0xd7, 0x64, 0x61, 0xe1, 0x3b, 0x83, 0xa9, 0x20, 0xe9, 0xf4, 0x75, 0x9a,
      0x4d, 0xdf, 0x29, 0x38, 0xac, 0xd6, 0x44, 0x73, 0x04, 0x6b, 0x7d, 0x95,
      0xda, 0x4e, 0x85, 0x50, 0x7b, 0x94, 0xee, 0xe1, 0x34, 0x5f, 0x1b, 0x09,
      0x65, 0x94, 0x4a, 0x39, 0xc9, 0x66, 0xbf, 0x91, 0x4c, 0x06, 0x19, 0x1b,
      0xeb, 0x3a, 0x05, 0x07, 0xdc, 0x05, 0xd9, 0x2f, 0x16, 0x59, 0x5a, 0x5a,
      0xa2, 0xb4, 0xd1, 0xc4, 0xda, 0x0f, 0x50, 0x6f, 0x18, 0xb4, 0x7b, 0x3e,
      0x71, 0x75, 0x3e, 0xc7, 0xed, 0xc7, 0x0f, 0xc4, 0x98, 0x9d, 0x9d, 0xc7,
      0x34, 0x03, 0xa4, 0x52, 0x09, 0x74, 0xbd, 0xf7, 0x04, 0xbe, 0xb0, 0xf1,
      0x42, 0x16, 0x3f, 0xac, 0x91, 0x2b, 0x18, 0x78, 0xc2, 0x31, 0xe2, 0x63,
      0x6d, 0xdc, 0x1c, 0xe9, 0xc7, 0x17, 0x08, 0x52, 0xdb, 0xb6, 0xc9, 0x3e,
      0x7f, 0x49, 0x7e, 0xe6, 0x15, 0x86, 0x65, 0xd9, 0x8c, 0x8f, 0xb7, 0x61,
      0x9a, 0xc9, 0x13, 0x58, 0xcf, 0x3f, 0x93, 0x77, 0xaf, 0x57, 0xd8, 0xb9,
      0x18, 0x65, 0xf0, 0x5e, 0x92, 0xc1, 0xf1, 0x1b, 0xb8, 0x44, 0x14, 0x40,
      0x03, 0xd0, 0xcd, 0x4d, 0xb1, 0x1b, 0x21, 0x76, 0xec, 0x6f, 0x67, 0xd7,
      0x98, 0x08, 0x94, 0xe5, 0xcd, 0xdb, 0x45, 0x0e, 0x82, 0x31, 0x26, 0x1f,
      0x3d, 0x04, 0x7f, 0x5c, 0xfd, 0x6e, 0xbe, 0x6e, 0xbd, 0x97, 0xcc, 0xf4,
      0x3a, 0x5e, 0xd1, 0xb8, 0x9c, 0x7e, 0x88, 0xd1, 0xd1, 0x11, 0xa1, 0x54,
      0x72, 0xb0, 0xac, 0x15, 0xe9, 0xe8, 0x18, 0x51, 0xf5, 0xba, 0xcd, 0x8f,
      0xea, 0x01, 0xa1, 0x78, 0x83, 0x20, 0x45, 0xcc, 0xe6, 0xb6, 0xb4, 0x9a,
      0x0d, 0xaa, 0xd5, 0x2a, 0x85, 0x8f, 0xab, 0x7c, 0x59, 0x2c, 0x72, 0xa8,
      0x5d, 0x62, 0xfc, 0xfe, 0x18, 0x46, 0xf7, 0x5d, 0xa5, 0xa6, 0xa6, 0x9e,
      0xc8, 0xf4, 0xf4, 0x17, 0xa2, 0x51, 0x97, 0x91, 0xe1, 0x5e, 0xa2, 0xa1,
      0x06, 0x9b, 0xb3, 0x4f, 0x29, 0xe7, 0xca, 0x38, 0xc5, 0x0b, 0x38, 0x7b,
      0xa0, 0xeb, 0x1a, 0x8e, 0x3f, 0x42, 0x35, 0x10, 0xa3, 0x2b, 0xd5, 0xc6,
      0xc8, 0xad, 0x21, 0x22, 0xf1, 0xb4, 0x02, 0x50, 0x99, 0x4c, 0x86, 0x5a,
      0xad, 0x20, 0xcb, 0xcb, 0x6b, 0x2c, 0x2f, 0x7b, 0xa8, 0x1f, 0x2a, 0xdc,
      0x9a, 0x8d, 0xb3, 0x5f, 0xa7, 0xbe, 0x27, 0xe8, 0x2d, 0xa1, 0x27, 0x2e,
      0xa4, 0x46, 0x4d, 0xae, 0x0c, 0xc5, 0x69, 0xef, 0x8a, 0xa0, 0x79, 0xfa,
      0x4f, 0xfc, 0x52, 0xc7, 0xdf, 0xe8, 0xf1, 0x6c, 0x8a, 0xeb, 0x7a, 0x68,
      0xb5, 0xe4, 0xa8, 0xd3, 0x82, 0x96, 0x80, 0xae, 0x0c, 0xfc, 0xc1, 0x10,
      0xd2, 0xda, 0xe1, 0xd0, 0xe9, 0x3c, 0x73, 0x5c, 0x27, 0x26, 0xba, 0xee,
      0xd1, 0x0a, 0x35, 0xed, 0x57, 0x41, 0x83, 0xe3, 0xb4, 0x76, 0x08, 0xd0,
      0xf9, 0xc7, 0x4b, 0xfc, 0x09, 0x52, 0xcb, 0x07, 0x62, 0x36, 0x43, 0x92,
      0xc6, 0x00, 0x00, 0x00, 0x00, 0x49, 0x45, 0x4e, 0x44, 0xae, 0x42, 0x60,
      0x82
    };

    static const unsigned char image3_data[] =
    {
      0x89, 0x50, 0x4e, 0x47, 0x0d, 0x0a, 0x1a, 0x0a, 0x00, 0x00, 0x00, 0x0d,
      0x49, 0x48, 0x44, 0x52, 0x00, 0x00, 0x00, 0x10, 0x00, 0x00, 0x00, 0x10,
      0x08, 0x06, 0x00, 0x00, 0x00, 0x1f, 0xf3, 0xff, 0x61, 0x00, 0x00, 0x01,
      0xdd, 0x49, 0x44, 0x41, 0x54, 0x78, 0x9c, 0x8d, 0x93, 0xbd, 0x6a, 0x55,
      0x41, 0x10, 0xc7, 0x7f, 0xb3, 0x67, 0xf7, 0x28, 0x18, 0x44, 0x24, 0x82,
      0x56, 0x56, 0xa2, 0x4d, 0x7c, 0x01, 0x9b, 0x54, 0xe2, 0x03, 0xd8, 0x5a,
      0xf9, 0x04, 0x3e, 0x80, 0x8d, 0x45, 0xc0, 0x5a, 0xb0, 0xb1, 0xb4, 0x11,
      0x04, 0x5b, 0xb5, 0xb0, 0x92, 0x40, 0x08, 0x1a, 0x3f, 0x90, 0x34, 0x26,
      0x51, 0x50, 0x50, 0xf1, 0x92, 0xfb, 0xed, 0xfd, 0xd8, 0xdd, 0x19, 0x8b,
      0x73, 0x4e, 0x3c, 0x68, 0x72, 0x71, 0x60, 0x61, 0x77, 0x76, 0xf7, 0xb7,
      0x33, 0xfb, 0x9f, 0x91, 0xf5, 0xcd, 0xed, 0xce, 0xa9, 0x93, 0x4b, 0xc7,
      0xd4, 0x8c, 0xc6, 0xd4, 0x60, 0x3a, 0x9d, 0x11, 0xe7, 0x91, 0xce, 0x7e,
      0xdf, 0xa5, 0xa8, 0xdd, 0x5e, 0x7f, 0x7c, 0xed, 0xe6, 0x8d, 0xab, 0x1f,
      0xf8, 0xdb, 0x5e, 0xbd, 0xdd, 0x99, 0xa8, 0xaa, 0xb5, 0x47, 0xce, 0xd9,
      0xc6, 0x93, 0x68, 0x83, 0x61, 0xb4, 0x37, 0xef, 0xf7, 0xec, 0xe3, 0xee,
      0x37, 0x7b, 0xf6, 0x62, 0xeb, 0xe7, 0x83, 0x87, 0x4f, 0x57, 0xcd, 0x8c,
      0xf6, 0x70, 0xff, 0x10, 0x01, 0x11, 0xc1, 0x39, 0x10, 0xa7, 0xb8, 0xc2,
      0x71, 0xf6, 0xdc, 0x32, 0xab, 0x57, 0x56, 0x96, 0x2f, 0x5d, 0x38, 0xff,
      0xfc, 0xd1, 0x93, 0x97, 0xd7, 0xdb, 0x67, 0x0f, 0x05, 0x18, 0x20, 0x40,
      0x59, 0x7a, 0x3e, 0x7d, 0xfe, 0xc2, 0xc6, 0xc6, 0x16, 0xaf, 0xdf, 0x6d,
      0x93, 0xe2, 0xb4, 0x1c, 0x0e, 0x7a, 0x8f, 0x6f, 0xdf, 0xb9, 0x7f, 0x6b,
      0x21, 0x00, 0x40, 0x9c, 0xe0, 0x04, 0xce, 0x2c, 0x9f, 0x66, 0x1e, 0x8d,
      0x1f, 0xdf, 0xbb, 0xec, 0xec, 0x7e, 0xa5, 0x3c, 0x5e, 0x32, 0x1c, 0xf5,
      0xd6, 0x9a, 0x73, 0xfe, 0xd0, 0xcb, 0x80, 0x13, 0x21, 0x45, 0xe5, 0xf2,
      0xca, 0x45, 0x66, 0xd3, 0xc8, 0x7e, 0xb7, 0x4f, 0x67, 0xbf, 0x47, 0x08,
      0x25, 0x9b, 0x1b, 0xeb, 0x2c, 0x04, 0x98, 0x19, 0x39, 0x2b, 0x89, 0x4c,
      0xaf, 0xf7, 0x8b, 0xd1, 0x68, 0xca, 0x6c, 0x36, 0xa7, 0x0c, 0x81, 0x10,
      0x3c, 0x45, 0x51, 0x2c, 0x06, 0x24, 0x35, 0x86, 0xa3, 0x39, 0x29, 0x25,
      0x34, 0x1b, 0x65, 0x19, 0x10, 0x07, 0xaa, 0x4a, 0xf0, 0x01, 0x11, 0x8e,
      0x06, 0xc4, 0xac, 0x4c, 0x26, 0x91, 0xac, 0x55, 0x32, 0xae, 0x10, 0x3c,
      0x0e, 0x91, 0x92, 0x9c, 0x13, 0x21, 0x84, 0x3a, 0xc9, 0x06, 0x50, 0xcf,
      0xa7, 0x31, 0x23, 0x40, 0x4c, 0xb9, 0x56, 0xc0, 0xa1, 0x19, 0xc4, 0x41,
      0x12, 0x45, 0xb2, 0x02, 0x9e, 0xc2, 0x7b, 0xa4, 0x15, 0x82, 0x6f, 0x64,
      0x48, 0x49, 0x0f, 0x42, 0x93, 0x42, 0x20, 0x19, 0x82, 0x54, 0x35, 0x21,
      0x82, 0x15, 0xe0, 0x14, 0x0a, 0x81, 0x4a, 0xe8, 0xca, 0x9c, 0xb6, 0x7f,
      0x1e, 0x2a, 0x48, 0xe3, 0xac, 0xe9, 0xae, 0x70, 0x08, 0x0e, 0x11, 0xa9,
      0x5f, 0xff, 0x13, 0xc1, 0x91, 0x75, 0xf0, 0xbf, 0xe6, 0x1a, 0x82, 0x69,
      0xd5, 0x44, 0x58, 0x15, 0x60, 0xab, 0xb7, 0x30, 0xad, 0x17, 0x56, 0x49,
      0xdc, 0x4e, 0xc1, 0x37, 0x7b, 0x4b, 0x27, 0xc2, 0x81, 0x53, 0xa1, 0x52,
      0xc1, 0x20, 0x25, 0x23, 0xab, 0x91, 0x72, 0x46, 0x63, 0xc6, 0x87, 0x92,
      0xb6, 0x8e, 0xbe, 0x3f, 0x18, 0x8c, 0xd7, 0xee, 0xde, 0x4b, 0xce, 0xf9,
      0xfa, 0x11, 0x45, 0x70, 0x98, 0x29, 0x66, 0x60, 0x75, 0x67, 0x98, 0x19,
      0x98, 0x21, 0xe2, 0x18, 0x0f, 0x86, 0xb3, 0x06, 0xf0, 0x1b, 0xda, 0xc4,
      0xfd, 0x8b, 0x61, 0x0d, 0x6a, 0xd5, 0x00, 0x00, 0x00, 0x00, 0x49, 0x45,
      0x4e, 0x44, 0xae, 0x42, 0x60, 0x82
    };

    static const unsigned char image4_data[] =
    {
      0x89, 0x50, 0x4e, 0x47, 0x0d, 0x0a, 0x1a, 0x0a, 0x00, 0x00, 0x00, 0x0d,
      0x49, 0x48, 0x44, 0x52, 0x00, 0x00, 0x00, 0x10, 0x00, 0x00, 0x00, 0x10,
      0x08, 0x06, 0x00, 0x00, 0x00, 0x1f, 0xf3, 0xff, 0x61, 0x00, 0x00, 0x02,
      0xff, 0x49, 0x44, 0x41, 0x54, 0x78, 0x9c, 0x65, 0x93, 0x4f, 0x68, 0x5c,
      0x75, 0x10, 0xc7, 0x3f, 0xbf, 0xdf, 0xdb, 0x7d, 0xfb, 0xf6, 0x65, 0xff,
      0x34, 0x49, 0xbb, 0x9b, 0xd6, 0xda, 0xc4, 0xb8, 0xb0, 0x5a, 0x2a, 0x54,
      0x21, 0x1e, 0xd6, 0xa2, 0x55, 0xa4, 0x78, 0x08, 0x28, 0x8a, 0x7b, 0x10,
      0xb1, 0xd5, 0x4b, 0x84, 0x50, 0xf0, 0xa2, 0x78, 0xf0, 0x10, 0x10, 0x8f,
      0x5e, 0x96, 0xa2, 0xe4, 0x60, 0xbd, 0x94, 0xa0, 0x39, 0x59, 0x2a, 0x4a,
      0x95, 0x62, 0xb5, 0x22, 0xd6, 0x44, 0xab, 0xa8, 0x49, 0x9b, 0x76, 0x69,
      0x12, 0xb6, 0xcd, 0xfe, 0xe9, 0xf6, 0xed, 0xbe, 0xfd, 0x93, 0x7d, 0xbb,
      0xef, 0xbd, 0x9f, 0x87, 0xe0, 0x56, 0x93, 0x81, 0x61, 0x60, 0x98, 0xf9,
      0x30, 0xcc, 0x77, 0x46, 0x28, 0xa5, 0xd8, 0x6e, 0x4f, 0xbf, 0xbe, 0x6a,
      0xd8, 0x61, 0x7b, 0x02, 0xdc, 0x94, 0x52, 0xbe, 0x40, 0xa8, 0xbc, 0xf4,
      0x4a, 0x0b, 0x8b, 0xb3, 0x93, 0xed, 0xed, 0xb5, 0xe2, 0xbf, 0x80, 0x43,
      0xd3, 0x4b, 0x11, 0x43, 0x79, 0xef, 0x29, 0x2d, 0xf8, 0xa6, 0xd4, 0x63,
      0xf1, 0xfb, 0x07, 0x7b, 0x98, 0x41, 0x97, 0x6b, 0x65, 0x03, 0xdf, 0xa9,
      0xd7, 0xa5, 0xdb, 0xfa, 0x5c, 0xa9, 0xfa, 0x07, 0x8b, 0xb3, 0xcf, 0xae,
      0xef, 0x00, 0x3c, 0x3e, 0xbd, 0x34, 0xe2, 0xe1, 0x5f, 0x4c, 0x26, 0x86,
      0xd3, 0xc7, 0x9f, 0x89, 0x73, 0xe4, 0x61, 0x8d, 0x98, 0xa9, 0x01, 0xe0,
      0xf4, 0x14, 0xbf, 0xdd, 0x68, 0x73, 0xfa, 0x42, 0x93, 0xeb, 0xab, 0xe5,
      0x3b, 0xa6, 0x2a, 0xbd, 0x7c, 0xe9, 0xe3, 0xe7, 0x2e, 0xf6, 0x01, 0x87,
      0xa6, 0x97, 0x22, 0x21, 0xbc, 0xc5, 0xc7, 0x0e, 0xee, 0x4b, 0xbf, 0x35,
      0x69, 0x10, 0x0d, 0x07, 0x00, 0x38, 0xfb, 0xd3, 0x5d, 0x00, 0x9e, 0xcf,
      0x0c, 0x01, 0xd0, 0xf3, 0x7c, 0x3e, 0x39, 0x6f, 0x71, 0xfe, 0x97, 0x52,
      0x75, 0xc0, 0xfb, 0xfb, 0x91, 0x1f, 0x66, 0x5f, 0xdd, 0x08, 0x00, 0xe8,
      0xf4, 0x66, 0x92, 0x89, 0x44, 0xfa, 0xc4, 0xd1, 0x20, 0x4e, 0x17, 0x9c,
      0xae, 0x0b, 0xc0, 0x89, 0x63, 0x09, 0xa2, 0xa6, 0xc6, 0x4a, 0x61, 0x93,
      0xe6, 0xa6, 0x0f, 0xc0, 0x8b, 0x99, 0x18, 0x6b, 0x15, 0x7f, 0x78, 0x25,
      0x3f, 0x3e, 0x93, 0xcd, 0x66, 0x4f, 0x6a, 0xdf, 0xaf, 0x1d, 0x35, 0xba,
      0x7a, 0x6c, 0xee, 0xa5, 0x23, 0x71, 0x23, 0x66, 0x42, 0xdb, 0xf1, 0xfa,
      0x7e, 0x70, 0xd4, 0x04, 0xc0, 0x6a, 0xba, 0x14, 0xaa, 0xce, 0x56, 0xbe,
      0xe3, 0x13, 0x0f, 0xbb, 0xfc, 0x9a, 0x57, 0xe9, 0xa0, 0xe8, 0x9d, 0x09,
      0xd8, 0xe1, 0xd8, 0x84, 0x6e, 0x0c, 0xc4, 0x07, 0xc3, 0x1e, 0x95, 0x9a,
      0xbf, 0x43, 0x11, 0x80, 0x76, 0xc7, 0xa7, 0x52, 0x73, 0xb1, 0x5b, 0x2e,
      0x5f, 0x2e, 0x0d, 0x61, 0x3b, 0x83, 0xc8, 0x90, 0x8a, 0xd4, 0x9d, 0x54,
      0x26, 0x00, 0xa4, 0x86, 0x22, 0x92, 0x8a, 0xed, 0xa2, 0x09, 0x81, 0x26,
      0x77, 0x02, 0x6a, 0x6d, 0xc5, 0xc2, 0xcd, 0x00, 0x6b, 0x56, 0x84, 0x46,
      0x37, 0x04, 0x80, 0x0c, 0x1a, 0x48, 0xa9, 0xa5, 0xa4, 0x52, 0xbe, 0xe8,
      0x7a, 0x1a, 0xeb, 0x55, 0x8d, 0x90, 0x61, 0x70, 0xf2, 0x85, 0xbd, 0x1c,
      0x9b, 0xd8, 0xc3, 0x85, 0xeb, 0xbb, 0xfa, 0x80, 0x1b, 0xd5, 0x08, 0x7e,
      0x70, 0x90, 0x2f, 0xde, 0x0e, 0x93, 0x7b, 0x0d, 0x84, 0x00, 0x01, 0x08,
      0xa1, 0x64, 0x00, 0xa1, 0xf2, 0xf5, 0xb6, 0xe2, 0x6a, 0x75, 0x37, 0xfa,
      0xc0, 0x56, 0x43, 0xfa, 0xbe, 0x00, 0x33, 0x59, 0xad, 0x0f, 0x78, 0x74,
      0x14, 0xde, 0x78, 0xea, 0xde, 0x44, 0x01, 0x09, 0x1d, 0xd7, 0x41, 0xf7,
      0xda, 0x05, 0x29, 0xbd, 0xd2, 0x82, 0xe7, 0xd4, 0xeb, 0x52, 0xc0, 0x5f,
      0x05, 0x78, 0x67, 0xce, 0x03, 0xe0, 0xf0, 0xa8, 0xb8, 0x07, 0x18, 0xdb,
      0x8a, 0x7f, 0xac, 0x29, 0xde, 0xfd, 0xcc, 0xc3, 0x08, 0x78, 0xf8, 0x1d,
      0xab, 0x35, 0xec, 0x2e, 0x5f, 0x91, 0x8b, 0xb3, 0x93, 0x6d, 0xbc, 0xe6,
      0x69, 0xe9, 0x14, 0x39, 0x10, 0xb3, 0xa9, 0xd7, 0x6d, 0x3e, 0x3c, 0xdb,
      0xd8, 0xb1, 0x87, 0x6b, 0xb7, 0x7a, 0x7c, 0xf4, 0x55, 0x8d, 0xb1, 0xb8,
      0x8d, 0x74, 0x8a, 0x04, 0xdd, 0xf2, 0xb9, 0x44, 0xa8, 0x70, 0x5b, 0x02,
      0x98, 0xa2, 0xf2, 0x7e, 0xb9, 0x54, 0xcc, 0xeb, 0xbe, 0xcd, 0x48, 0xb4,
      0x83, 0xdd, 0x68, 0xf2, 0xe9, 0xb7, 0xd5, 0x7e, 0xf3, 0xcd, 0xa2, 0xc3,
      0xdc, 0x77, 0x77, 0x18, 0x89, 0x76, 0xd0, 0x7d, 0x9b, 0x72, 0x71, 0xa3,
      0x90, 0xf4, 0x2e, 0x9f, 0xb2, 0x2c, 0xab, 0xda, 0x3f, 0xe5, 0x27, 0xa6,
      0xbe, 0x7e, 0xb0, 0x17, 0x48, 0x7c, 0x73, 0x38, 0x15, 0x1f, 0x7f, 0xe8,
      0x40, 0x04, 0x21, 0x04, 0xfb, 0x77, 0x07, 0xd9, 0xbf, 0x47, 0xe7, 0xe7,
      0xe5, 0x16, 0x4a, 0x29, 0xae, 0xae, 0x37, 0xf9, 0x7d, 0xc5, 0x2a, 0x24,
      0xdd, 0x4b, 0x53, 0x0f, 0x84, 0x57, 0x7f, 0xcc, 0xe5, 0x72, 0xf6, 0xff,
      0x9e, 0xe9, 0xc9, 0xa9, 0x33, 0x7b, 0x37, 0xc5, 0xf8, 0x8c, 0x69, 0x9a,
      0xaf, 0x8c, 0x25, 0xf5, 0xc8, 0xae, 0xe8, 0x96, 0x64, 0xb5, 0x86, 0xc3,
      0x6a, 0xb1, 0xdb, 0x72, 0xda, 0x77, 0xcf, 0x25, 0xbd, 0xcb, 0xa7, 0x46,
      0xcd, 0x8d, 0x3f, 0x73, 0xb9, 0x9c, 0x0d, 0xdb, 0xbe, 0x11, 0x20, 0x9b,
      0xcd, 0xea, 0x8d, 0x81, 0xcc, 0xbe, 0xba, 0x96, 0xca, 0x48, 0xa9, 0xa5,
      0x84, 0x50, 0x52, 0x78, 0xed, 0xc2, 0xb0, 0xbb, 0x7c, 0x25, 0x11, 0x2a,
      0xdc, 0xb6, 0x2c, 0xab, 0x3a, 0x3f, 0x3f, 0xdf, 0xfd, 0xb7, 0xfe, 0x1f,
      0x7d, 0x29, 0x67, 0x4b, 0xc6, 0x4f, 0x8f, 0x91, 0x00, 0x00, 0x00, 0x00,
      0x49, 0x45, 0x4e, 0x44, 0xae, 0x42, 0x60, 0x82
    };

    static const unsigned char image5_data[] =
    {
      0x89, 0x50, 0x4e, 0x47, 0x0d, 0x0a, 0x1a, 0x0a, 0x00, 0x00, 0x00, 0x0d,
      0x49, 0x48, 0x44, 0x52, 0x00, 0x00, 0x00, 0x10, 0x00, 0x00, 0x00, 0x10,
      0x08, 0x06, 0x00, 0x00, 0x00, 0x1f, 0xf3, 0xff, 0x61, 0x00, 0x00, 0x02,
      0x18, 0x49, 0x44, 0x41, 0x54, 0x78, 0x9c, 0x7d, 0x92, 0xcb, 0x4a, 0x5c,
      0x41, 0x10, 0x86, 0xbf, 0xea, 0xee, 0x39, 0x73, 0xc6, 0x11, 0x41, 0x72,
      0x11, 0x11, 0x02, 0x51, 0xb3, 0x90, 0x90, 0x59, 0x8a, 0xb8, 0xc8, 0x22,
      0x2f, 0x12, 0x12, 0x92, 0x07, 0x08, 0x24, 0x8b, 0x11, 0x02, 0x82, 0x1b,
      0x4d, 0xde, 0x20, 0x31, 0x97, 0x95, 0x1b, 0x1f, 0x20, 0xcc, 0xd2, 0x80,
      0x8a, 0x06, 0xc4, 0x45, 0x44, 0x54, 0x22, 0x51, 0x88, 0x59, 0xa8, 0xa8,
      0x73, 0xc6, 0xb9, 0x9e, 0xd3, 0x95, 0xc5, 0x78, 0xd7, 0x49, 0x41, 0x6d,
      0xaa, 0xe8, 0xbf, 0xff, 0xff, 0xeb, 0x96, 0xd9, 0xd9, 0xf5, 0x83, 0xce,
      0xce, 0x76, 0xa7, 0xea, 0x89, 0x63, 0x8f, 0xf7, 0x8a, 0x31, 0x96, 0x95,
      0x95, 0xcd, 0xf6, 0xe5, 0xe5, 0xb5, 0xe7, 0x13, 0x13, 0x2f, 0xbe, 0xf0,
      0xbf, 0x5a, 0x5a, 0xda, 0x6c, 0x78, 0xef, 0xf5, 0xb4, 0x6b, 0x35, 0xaf,
      0xfb, 0xfb, 0x15, 0x5d, 0x5d, 0xfd, 0xab, 0x53, 0x53, 0xdf, 0x75, 0x64,
      0xe4, 0xd3, 0x1b, 0x55, 0xa5, 0x55, 0x9b, 0xab, 0x82, 0xd6, 0x42, 0x10,
      0x04, 0x6c, 0x6c, 0xfc, 0x61, 0x6f, 0xef, 0x80, 0x30, 0x0c, 0xde, 0xe5,
      0xf3, 0x93, 0x6f, 0x5b, 0x19, 0xb8, 0x26, 0x20, 0x02, 0x49, 0x12, 0x93,
      0xcb, 0xf5, 0x33, 0x38, 0xf8, 0x90, 0x81, 0x81, 0x5e, 0x44, 0xcc, 0xfd,
      0x56, 0x02, 0xee, 0x26, 0x01, 0x6b, 0x2d, 0xc5, 0x62, 0x99, 0x4a, 0xa5,
      0x4e, 0x47, 0x47, 0x96, 0xe1, 0xe1, 0x47, 0x8f, 0x17, 0x16, 0x7e, 0x7d,
      0x88, 0xe3, 0x38, 0x53, 0x2c, 0x56, 0xb2, 0x33, 0x33, 0x3f, 0xbe, 0x8d,
      0x8f, 0xbf, 0xfc, 0x78, 0x23, 0x83, 0x28, 0xf2, 0xba, 0xbd, 0x7d, 0xa8,
      0x73, 0x73, 0xeb, 0x7a, 0x74, 0x54, 0xbd, 0x34, 0xdf, 0xd9, 0x89, 0x74,
      0x6b, 0xeb, 0x50, 0xa7, 0xa7, 0xe7, 0x74, 0x74, 0xf4, 0xeb, 0x6b, 0x55,
      0xbd, 0xec, 0xa0, 0x5a, 0x85, 0x72, 0xb9, 0x8a, 0x31, 0x16, 0x6b, 0x2d,
      0xe9, 0xf4, 0xf9, 0x3a, 0x9b, 0x05, 0xef, 0x33, 0xec, 0xee, 0x1e, 0x12,
      0x45, 0x11, 0xe9, 0x74, 0xf8, 0x3e, 0x9f, 0x9f, 0xbc, 0x7b, 0xc6, 0x40,
      0x15, 0xea, 0xf5, 0x06, 0xce, 0x19, 0xc2, 0xd0, 0x21, 0x22, 0xa8, 0x5e,
      0xc9, 0xeb, 0x84, 0x20, 0x48, 0xd1, 0xdd, 0x7d, 0x87, 0xfe, 0xfe, 0x7b,
      0xb4, 0xb5, 0x85, 0xaf, 0xce, 0x04, 0xea, 0x75, 0x10, 0xb1, 0x88, 0xb4,
      0xc2, 0x05, 0xc6, 0x08, 0xe9, 0x74, 0x8a, 0xae, 0xae, 0x5b, 0xf4, 0xf4,
      0xdc, 0xc6, 0xb9, 0x14, 0xce, 0x98, 0xf3, 0x87, 0x70, 0x4e, 0x30, 0xc6,
      0x91, 0x24, 0x70, 0x71, 0x7e, 0x5a, 0xe5, 0x72, 0x9d, 0x52, 0xa9, 0x8c,
      0xf7, 0x4a, 0x26, 0x93, 0x46, 0x44, 0x70, 0x51, 0x54, 0x21, 0x8a, 0x1a,
      0x44, 0x51, 0x15, 0xef, 0x15, 0xef, 0x3d, 0x49, 0xa2, 0xd4, 0xeb, 0x31,
      0x72, 0xc1, 0x4e, 0x14, 0x79, 0x8e, 0x8f, 0x6b, 0xa8, 0x0a, 0xd9, 0x6c,
      0x06, 0x6b, 0x2d, 0xce, 0x81, 0x2b, 0x14, 0xe6, 0xa7, 0x0a, 0x85, 0x79,
      0x73, 0xf9, 0xa6, 0x5a, 0x25, 0x97, 0x7b, 0xf0, 0x64, 0x68, 0xa8, 0xb7,
      0xef, 0x14, 0x6e, 0x1c, 0x37, 0x08, 0x02, 0x87, 0x31, 0x06, 0x63, 0xe4,
      0xc4, 0xa1, 0xe0, 0xc6, 0xc6, 0x9e, 0x3d, 0xbd, 0x29, 0x6f, 0xa1, 0xb0,
      0xfc, 0xd9, 0x7b, 0xfa, 0x8c, 0x81, 0x24, 0xf1, 0x58, 0x6b, 0x88, 0x63,
      0xc5, 0x5a, 0x83, 0xb5, 0x06, 0x11, 0xdb, 0xe4, 0xd2, 0x0a, 0x58, 0x18,
      0x86, 0x27, 0x51, 0x00, 0x9a, 0x51, 0xac, 0xbd, 0x4e, 0xf8, 0xda, 0x4f,
      0x3c, 0x27, 0xde, 0x3c, 0x98, 0x4a, 0x81, 0x88, 0xe0, 0x7d, 0x0a, 0x11,
      0x8f, 0x48, 0x0c, 0x80, 0x88, 0x69, 0x42, 0x6c, 0x25, 0xb0, 0xb8, 0xf8,
      0xb3, 0xb6, 0xb6, 0xf6, 0x1b, 0xa0, 0x04, 0xa0, 0xaa, 0x24, 0x49, 0x73,
      0xa7, 0x9a, 0x20, 0x02, 0xa5, 0x52, 0x25, 0xfe, 0x07, 0xd9, 0xbb, 0x1d,
      0xb5, 0x0b, 0x97, 0x87, 0x95, 0x00, 0x00, 0x00, 0x00, 0x49, 0x45, 0x4e,
      0x44, 0xae, 0x42, 0x60, 0x82
    };

    static const unsigned char image6_data[] =
    {
      0x89, 0x50, 0x4e, 0x47, 0x0d, 0x0a, 0x1a, 0x0a, 0x00, 0x00, 0x00, 0x0d,
      0x49, 0x48, 0x44, 0x52, 0x00, 0x00, 0x00, 0x10, 0x00, 0x00, 0x00, 0x10,
      0x08, 0x06, 0x00, 0x00, 0x00, 0x1f, 0xf3, 0xff, 0x61, 0x00, 0x00, 0x00,
      0x8f, 0x49, 0x44, 0x41, 0x54, 0x78, 0x9c, 0xa5, 0x52, 0x49, 0x12, 0xc0,
      0x20, 0x08, 0x8b, 0x4e, 0x1f, 0xc6, 0xd3, 0x7c, 0x1a, 0x3f, 0xb3, 0x07,
      0x6d, 0x59, 0x9a, 0xee, 0x78, 0x10, 0x42, 0x58, 0x26, 0x5a, 0x00, 0x85,
      0x99, 0xe0, 0x3a, 0xce, 0x26, 0x28, 0x80, 0xa2, 0x25, 0x52, 0x9b, 0x27,
      0x62, 0x42, 0x79, 0x35, 0x12, 0xb6, 0x22, 0x75, 0x8d, 0x84, 0x34, 0x32,
      0xac, 0x62, 0x0f, 0x3c, 0x91, 0xf9, 0xfe, 0xb6, 0xfc, 0x32, 0xa6, 0x69,
      0x5a, 0xf9, 0x29, 0xe6, 0xd6, 0xfd, 0x7a, 0x9f, 0x88, 0xc8, 0x04, 0xe3,
      0x22, 0x46, 0xc7, 0xf9, 0x47, 0x4c, 0x29, 0xaf, 0xc6, 0x95, 0xd8, 0x9a,
      0x71, 0x6a, 0xc6, 0x16, 0xeb, 0xe8, 0x89, 0x32, 0x67, 0xe4, 0xe2, 0xcc,
      0xd3, 0xa1, 0x41, 0x2c, 0x7c, 0x17, 0x17, 0x00, 0x87, 0xa7, 0xe1, 0x3f,
      0x91, 0x63, 0x70, 0xe2, 0x74, 0x00, 0x3d, 0x09, 0x76, 0x83, 0xc1, 0xbe,
      0xf2, 0x57, 0x73, 0x22, 0xb6, 0x92, 0x93, 0x4f, 0xb1, 0x5f, 0xb6, 0x02,
      0x6f, 0x3b, 0x53, 0x57, 0x71, 0xe6, 0x68, 0xdf, 0x00, 0x00, 0x00, 0x00,
      0x49, 0x45, 0x4e, 0x44, 0xae, 0x42, 0x60, 0x82
    };

    switch (id)
      {
        case image0_ID: {QImage img; img.loadFromData(image0_data, sizeof(image0_data), "PNG"); return QPixmap::fromImage(img);}
        case image1_ID: {QImage img; img.loadFromData(image1_data, sizeof(image1_data), "PNG"); return QPixmap::fromImage(img);}
        case image2_ID: {QImage img; img.loadFromData(image2_data, sizeof(image2_data), "PNG"); return QPixmap::fromImage(img);}
        case image3_ID: {QImage img; img.loadFromData(image3_data, sizeof(image3_data), "PNG"); return QPixmap::fromImage(img);}
        case image4_ID: {QImage img; img.loadFromData(image4_data, sizeof(image4_data), "PNG"); return QPixmap::fromImage(img);}
        case image5_ID: {QImage img; img.loadFromData(image5_data, sizeof(image5_data), "PNG"); return QPixmap::fromImage(img);}
        case image6_ID: {QImage img; img.loadFromData(image6_data, sizeof(image6_data), "PNG"); return QPixmap::fromImage(img);}
        default: return QPixmap();
      } // switch
  } // icon
};

#endif // CQFITTINGITEMWIDGET_H
