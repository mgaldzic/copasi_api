// Begin CVS Header
//   $Source: /fs/turing/cvs/copasi_dev/copasi/UI/CQMathMatrixWidget.cpp,v $
//   $Revision: 1.8 $
//   $Name: Build-33 $
//   $Author: pwilly $
//   $Date: 2009/10/13 12:16:55 $
// End CVS Header

// Copyright (C) 2008 by Pedro Mendes, Virginia Tech Intellectual
// Properties, Inc., EML Research, gGmbH, University of Heidelberg,
// and The University of Manchester.
// All rights reserved.

// Copyright (C) 2001 - 2007 by Pedro Mendes, Virginia Tech Intellectual
// Properties, Inc. and EML Research, gGmbH.
// All rights reserved.

//#include <qpushbutton.h>
#include <qlayout.h>
#include <qlineedit.h>
#include <qlabel.h>
#include <q3table.h>
//Added by qt3to4:
#include <Q3GridLayout>

#include "copasi.h"

#include "CQMathMatrixWidget.h"
#include "qtUtilities.h"
//#include "sensitivities/CSensTask.h"
//#include "sensitivities/CSensProblem.h"

#include "CopasiDataModel/CCopasiDataModel.h"
#include "report/CCopasiRootContainer.h"
//#include "utilities/CCopasiVector.h"

#include "model/CModel.h"

#include <qtabwidget.h>

/**
 *  Constructs a CQMathMatrixWidget which is a child of 'parent', with the
 *  name 'name' and widget flags set to 'f'.
 */
CQMathMatrixWidget::CQMathMatrixWidget(QWidget* parent, const char* name, Qt::WFlags fl)
    : CopasiWidget(parent, name, fl)
{
  if (!name)
    setName("CQMathMatrixWidget");

  setCaption("CQMathMatrixWidget");

  mWidgetLayout = new Q3GridLayout(this, 1, 1, 11, 6, "CQMathMatrixWidgetLayout");

  // **********  Label **************
  //mLabelTitle = new QLabel(this, "MatrixLabel");
  //mLabelTitle->setText("Matrices");
  //mLabelTitle->setAlignment(int(Qt::AlignVCenter
  //                              | Qt::AlignLeft));
  //mLabelTitle->setSizePolicy(QSizePolicy::Fixed, QSizePolicy::Fixed);
  //mWidgetLayout->addWidget(mLabelTitle, 0, 0);

  // tab widget
  mpTab = new QTabWidget(this, "TabWidget");
  mWidgetLayout->addMultiCellWidget(mpTab, 0, 2, 0, 2);

  // 1
  mArrayWidget1 = new CQArrayAnnotationsWidget(mpTab, "ArrayWidget1");
  //mArrayWidget1->setColorCoding(new CColorScale1());
  //CColorScaleAverage * tcs = new CColorScaleAverage();
  CColorScaleSimple * tcs = new CColorScaleSimple();
  mArrayWidget1->setColorCoding(tcs);
  tcs->setMinMax(-1.5, 1.5);
  //tcs->setMinMax(-1,1);
  //tcs->setSymmetric(true);
  //tcs->setFactor(3.0);
  mArrayWidget1->setColorScalingAutomatic(false);
  mpTab->addTab(mArrayWidget1, "Stoichiometry Matrix");

  // 2
  mArrayWidget2 = new CQArrayAnnotationsWidget(mpTab, "ArrayWidget2");
  //mArrayWidge2->setColorCoding(new CColorScale1());
  tcs = new CColorScaleSimple();
  mArrayWidget2->setColorCoding(tcs);
  tcs->setMinMax(-1.5, 1.5);
  //tcs2->setMinMax(-1,1);
  //tcs2->setSymmetric(true);
  //tcs2->setFactor(3.0);
  mArrayWidget2->setColorScalingAutomatic(false);
  mpTab->addTab(mArrayWidget2, "Reduced Stoichiometry Matrix");

  // 3
  mArrayWidget3 = new CQArrayAnnotationsWidget(mpTab, "ArrayWidget3");
  //mArrayWidge2->setColorCoding(new CColorScale1());
  tcs = new CColorScaleSimple();
  mArrayWidget3->setColorCoding(tcs);
  tcs->setMinMax(-1.5, 1.5);
  //tcs2->setMinMax(-1,1);
  //tcs2->setSymmetric(true);
  //tcs2->setFactor(3.0);
  mArrayWidget3->setColorScalingAutomatic(false);
  mpTab->addTab(mArrayWidget3, "Link Matrix");
}

/*
 *  Destroys the object and frees any allocated resources
 */
CQMathMatrixWidget::~CQMathMatrixWidget()
{}

void CQMathMatrixWidget::loadMatrices()
{

  assert(CCopasiRootContainer::getDatamodelList()->size() > 0);
  const CModel* pModel = (*CCopasiRootContainer::getDatamodelList())[0]->getModel();

  //mpResult = pProblem->getResultAnnotated();
  //mpScaledResult = pProblem->getScaledResultAnnotated();

  const CArrayAnnotation * tmp;

  tmp = dynamic_cast<const CArrayAnnotation *>
        (pModel->getObject(CCopasiObjectName("Array=Stoichiometry(ann)")));
  mArrayWidget1->setArrayAnnotation(tmp);

  tmp = dynamic_cast<const CArrayAnnotation *>
        (pModel->getObject(CCopasiObjectName("Array=Reduced stoichiometry(ann)")));
  mArrayWidget2->setArrayAnnotation(tmp);

  tmp = dynamic_cast<const CArrayAnnotation *>
        (pModel->getObject(CCopasiObjectName("Array=Link matrix(ann)")));
  mArrayWidget3->setArrayAnnotation(tmp);
}

void CQMathMatrixWidget::clearArrays()
{
  mArrayWidget1->setArrayAnnotation(NULL);
  mArrayWidget2->setArrayAnnotation(NULL);
  mArrayWidget3->setArrayAnnotation(NULL);
}

//*************************************

bool CQMathMatrixWidget::update(ListViews::ObjectType C_UNUSED(objectType), ListViews::Action
                                C_UNUSED(action), const std::string & C_UNUSED(key))
{
  clearArrays();
  return true;
}

bool CQMathMatrixWidget::leave()
{
  return true;
}

bool CQMathMatrixWidget::enterProtected()
{
  loadMatrices();

  return true;
}
