// Begin CVS Header
//   $Source: /fs/turing/cvs/copasi_dev/copasi/UI/CQTSSAWidget.h,v $
//   $Revision: 1.12 $
//   $Name: Build-33 $
//   $Author: aekamal $
//   $Date: 2010/05/10 16:12:15 $
// End CVS Header

// Copyright (C) 2010 by Pedro Mendes, Virginia Tech Intellectual
// Properties, Inc., University of Heidelberg, and The University
// of Manchester.
// All rights reserved.

// Copyright (C) 2008 by Pedro Mendes, Virginia Tech Intellectual
// Properties, Inc., EML Research, gGmbH, University of Heidelberg,
// and The University of Manchester.
// All rights reserved.

#ifndef CQTSSAWIDGET_H
#define CQTSSAWIDGET_H

#include <qvariant.h>
#include "TaskWidget.h"

#include "ui_CQTSSAWidget.h"

class CTSSAProblem;
class CQValidatorInt;
class CQValidatorDouble;
/*
#include <QtCore/QVariant>
#include <QtGui/QAction>
#include <QtGui/QApplication>
#include <QtGui/QButtonGroup>
#include <QtGui/QCheckBox>
#include <QtGui/QFrame>
#include <QtGui/QGridLayout>
#include <QtGui/QLabel>
#include <QtGui/QLineEdit>
#include <QtGui/QVBoxLayout>
#include "TaskWidget.h"
 */

#define TSSAMAX 10000000

class CQTSSAWidget : public TaskWidget, public Ui::CQTSSAWidget
{
  Q_OBJECT

public:
  CQTSSAWidget(QWidget* parent = 0, const char* name = 0);
  ~CQTSSAWidget();

  virtual bool runTask();


public slots:
  void slotDuration();
  void slotIntervalSize();
  void slotIntervals();

protected:
  virtual bool saveTask();
  virtual bool loadTask();
  virtual bool taskFinishedEvent();
  virtual CCopasiMethod * createMethod(const CCopasiMethod::SubType & type);

protected slots:
  virtual void languageChange();

private:
  CTSSAProblem * mpTSSAProblem;
  CQValidatorDouble * mpValidatorDuration;
  CQValidatorDouble * mpValidatorIntervalSize;
  CQValidatorDouble * mpValidatorDelay;
  CQValidatorInt * mpValidatorIntervals;

  void init();
  void destroy();
  void checkTimeSeries();
};

#endif // CQTSSAWIDGET_H
