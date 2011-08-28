// Begin CVS Header
//   $Source: /fs/turing/cvs/copasi_dev/copasi/UI/CQTaskThread.h,v $
//   $Revision: 1.3 $
//   $Name: Build-33 $
//   $Author: shoops $
//   $Date: 2010/04/12 19:26:39 $
// End CVS Header

// Copyright (C) 2010 by Pedro Mendes, Virginia Tech Intellectual
// Properties, Inc., University of Heidelberg, and The University
// of Manchester.
// All rights reserved.

#ifndef CQTaskThread_H
#define CQTaskThread_H

#include <QThread>

class TaskWidget;
class CCopasiException;

class CQTaskThread : public QThread
{
  Q_OBJECT

public:
  CQTaskThread(TaskWidget *tw);
  ~CQTaskThread();

  virtual void run();

protected:
  TaskWidget *mpTaskWidget;

signals:
  void exceptionOccured(CCopasiException* pException);
};

#endif //CQTaskThread_H

