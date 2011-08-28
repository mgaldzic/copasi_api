// Begin CVS Header
//   $Source: /fs/turing/cvs/copasi_dev/copasi/layoutUI/CQGLViewport.cpp,v $
//   $Revision: 1.8 $
//   $Name: Build-33 $
//   $Author: gauges $
//   $Date: 2010/02/03 13:53:00 $
// End CVS Header

// Copyright (C) 2008 by Pedro Mendes, Virginia Tech Intellectual
// Properties, Inc., EML Research, gGmbH, University of Heidelberg,
// and The University of Manchester.
// All rights reserved.

#include "CQGLViewport.h"

#include <qgl.h>
#include <QScrollBar>
#include <QRect>
#include <QHBoxLayout>
#include <QLayout>
#include <QResizeEvent>
#include <QVBoxLayout>
#include <QFrame>

#include <iostream>

#include "copasi/layout/CLBase.h"
#include "copasi/layout/CLayout.h"
#include "CQGLNetworkPainter.h"

/**
 * Constructor.
 */
CQGLViewport::CQGLViewport(QWidget* pParent, Qt::WFlags f):
    QFrame(pParent, f)
    , mpVerticalScrollbar(new QScrollBar(Qt::Vertical, NULL))
    , mpHorizontalScrollbar(new QScrollBar(Qt::Horizontal, NULL))
    , mpNetworkPainter(NULL)
{
  QVBoxLayout* pVBoxLayout = new QVBoxLayout();
  this->setLayout(pVBoxLayout);
  QFrame* pHBox = new QFrame(this);
  pHBox->setLayout(new QHBoxLayout());
  pVBoxLayout->addWidget(pHBox);
  QGLFormat format;
  format.setDoubleBuffer(TRUE);
  this->mpNetworkPainter = new CQGLNetworkPainter(format, pHBox);
  pHBox->layout()->addWidget(this->mpNetworkPainter);
  pHBox->layout()->addWidget(this->mpVerticalScrollbar);
  this->mpVerticalScrollbar->setSingleStep(1);
  pVBoxLayout->addWidget(this->mpHorizontalScrollbar);
  this->mpHorizontalScrollbar->setSingleStep(1);
  connect(this->mpVerticalScrollbar, SIGNAL(valueChanged(int)), this, SLOT(slotVValueChanged(int)));
  connect(this->mpHorizontalScrollbar, SIGNAL(valueChanged(int)), this, SLOT(slotHValueChanged(int)));
}

/**
 * Destructor.
 */
CQGLViewport::~CQGLViewport()
{}

void CQGLViewport::resizeEvent(QResizeEvent* e)
{
  this->updateScrollbars();
  QFrame::resizeEvent(e);
}

const CQGLNetworkPainter* CQGLViewport::getPainter() const
{
  return this->mpNetworkPainter;
}

CQGLNetworkPainter* CQGLViewport::getPainter()
{
  return this->mpNetworkPainter;
}

void CQGLViewport::createGraph(CLayout* IP)
{
  this->mpNetworkPainter->createGraph(IP);
  this->updateScrollbars();
}

void CQGLViewport::setZoomFactor(C_FLOAT64 zoom)
{
  // TODO make sure the GL window is only redrawn once
  this->mpNetworkPainter->setZoomFactor(zoom);
  this->updateScrollbars();
  this->mpNetworkPainter->update();
}

void CQGLViewport::updateScrollbars()
{
  // reset the scollbar range
  // TODO check te setting for the scroll range since there seem to be some
  // error messages
  // disconnect the scrollbar listeners and handle the update so that the GL
  // window is only redrawn once
  double zoom = this->mpNetworkPainter->getZoomFactor();
  CLPoint max = this->mpNetworkPainter->getGraphMax();
  CLPoint min = this->mpNetworkPainter->getGraphMin();
  double graphWidth = (max.getX() - min.getX()) * zoom;
  double graphHeight = (max.getY() - min.getY()) * zoom;
  double rectangleHeight = this->contentsRect().height();
  double rectangleWidth = this->contentsRect().width();

  if (graphHeight < rectangleHeight)
    {
      this->mpVerticalScrollbar->hide();
      this->mpVerticalScrollbar->setValue(0);
    }
  else
    {
      this->mpVerticalScrollbar->setPageStep(rectangleHeight);
      this->mpVerticalScrollbar->setRange(0, (unsigned int)(graphHeight - rectangleHeight));
      this->mpVerticalScrollbar->show();
      this->mpNetworkPainter->update();
    }

  if (graphWidth < rectangleWidth)
    {
      this->mpHorizontalScrollbar->hide();
      this->mpHorizontalScrollbar->setValue(0);
    }
  else
    {
      this->mpHorizontalScrollbar->setPageStep(rectangleWidth);
      this->mpHorizontalScrollbar->setRange(0, (unsigned int)(graphWidth - rectangleWidth));
      this->mpHorizontalScrollbar->show();
      this->mpNetworkPainter->update();
    }
}

void CQGLViewport::slotVValueChanged(int value)
{
  CLPoint p = this->mpNetworkPainter->getGraphMin();
  double zoom = this->mpNetworkPainter->getZoomFactor();
  this->mpNetworkPainter->setCurrentPositionY((double)(p.getY() + value / zoom));
}

void CQGLViewport::slotHValueChanged(int value)
{
  CLPoint p = this->mpNetworkPainter->getGraphMin();
  double zoom = this->mpNetworkPainter->getZoomFactor();
  this->mpNetworkPainter->setCurrentPositionX((double)(p.getX() + value / zoom));
}

void CQGLViewport::resetView()
{
  // TODO disconnect the scrollbar listeners
  // so that the display is only redrawn once
  this->mpNetworkPainter->resetView();
  this->updateScrollbars();
}

void CQGLViewport::updateWidget()
{
  this->mpNetworkPainter->update();
}

bool CQGLViewport::isCircleMode() const
{
  return this->mpNetworkPainter->isCircleMode();
}
