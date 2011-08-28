/* Begin CVS Header
 $Source: /fs/turing/cvs/copasi_dev/copasi/plotUI/scrollzoomer.cpp,v $
 $Revision: 1.5 $
 $Name: Build-33 $
 $Author: shoops $
 $Date: 2008/12/18 19:04:22 $
 End CVS Header */

// Copyright (C) 2008 by Pedro Mendes, Virginia Tech Intellectual
// Properties, Inc., EML Research, gGmbH, University of Heidelberg,
// and The University of Manchester.
// All rights reserved.

// Copyright � 1997   Josef Wilgen
// Copyright � 2002   Uwe Rathmann
//
// This file is published under the Qwt License, Version 1.0.
// You should have received a copy of this licence in the file
// QwtLicense.
//
// Modifications made to the original are
// Copyright � 2006 by Pedro Mendes, Virginia Tech Intellectual
// Properties, Inc. and EML Research, gGmbH.
// All rights reserved.

#include <qevent.h>
#include <qwt_plot_canvas.h>
#include <qwt_plot_layout.h>
#include <qwt_scale_engine.h>
//Added by qt3to4:
#include <QResizeEvent>
#include <QChildEvent>
#include "scrollbar.h"
#include "scrollzoomer.h"

LogPlotZoomer::LogPlotZoomer(QwtPlotCanvas *canvas):
    QwtPlotZoomer(canvas)
{}

QwtText LogPlotZoomer::trackerText(const QwtDoublePoint &pos) const
  {
    switch (rubberBand())
      {
      case HLineRubberBand:
        return QString().sprintf("%.4g", pos.y());
      case VLineRubberBand:
        return QString().sprintf("%.4g", pos.x());
      default:
        return QString().sprintf("%.4g, %.4g", pos.x(), pos.y());
      }
    return QwtText(); // make some dumb compilers happy
  }

/*void QwtPlotZoomer::move(double x, double y)
{
    x = qwtMax(x, zoomBase().left());
    x = qwtMin(x, zoomBase().right() - zoomRect().width());

    y = qwtMax(y, zoomBase().top());
    y = qwtMin(y, zoomBase().bottom() - zoomRect().height());

    if (x != zoomRect().left() || y != zoomRect().top())
    {
        d_data->zoomStack[d_data->zoomRectIndex].moveTo(x, y);
        rescale();
    }
}*/

void LogPlotZoomer::move(double x, double y)
{
  //QwtPlotZoomer::move(x,y);

  x = qwtMax(x, zoomBase().left());
  x = qwtMin(x, zoomBase().right() - zoomRect().width());

  y = qwtMax(y, zoomBase().top());
  y = qwtMin(y, zoomBase().bottom() - zoomRect().height());

  if (x != zoomRect().left() || y != zoomRect().top())
    {
      //zoomStack()[zoomRectIndex()].moveTo(x, y);
      QwtDoubleRect & rect = const_cast<QwtDoubleRect &>(zoomStack()[zoomRectIndex()]);

      //handle x axis
      const int xAxis = QwtPlotZoomer::xAxis();
      const QwtScaleEngine *sex = plot()->axisScaleEngine(xAxis);
      if (dynamic_cast<const QwtLog10ScaleEngine*>(sex))
        {
          //logarithmic
          double factor = rect.right() / rect.left();
          rect.setRight(x*factor);
          rect.setLeft(x);
        }
      else
        {
          rect.moveLeft(x);
        }

      const int yAxis = QwtPlotZoomer::yAxis();
      const QwtScaleEngine *sey = plot()->axisScaleEngine(yAxis);
      if (dynamic_cast<const QwtLog10ScaleEngine*>(sey))
        {
          //logarithmic
          double factor = rect.bottom() / rect.top();
          rect.setBottom(y*factor);
          rect.setTop(y);
        }
      else
        {
          rect.moveTop(y);
        }

      //zoomStack()[zoomRectIndex()].moveTo(x, y);
      rescale();
    }
}

//******************************************

class ScrollData
  {
  public:
    ScrollData():
        scrollBar(NULL),
        position(ScrollZoomer::OppositeToScale),
#if QT_VERSION < 0x040000
        mode(Q3ScrollView::Auto)
#else
        mode(Qt::ScrollBarAsNeeded)
#endif
    {
    }

    ~ScrollData()
    {
      delete scrollBar;
    }

    ScrollBar *scrollBar;
    ScrollZoomer::ScrollBarPosition position;
#if QT_VERSION < 0x040000
    Q3ScrollView::ScrollBarMode mode;
#else
    Qt::ScrollBarPolicy mode;
#endif
  };

//******************************************

ScrollZoomer::ScrollZoomer(QwtPlotCanvas *canvas):
    LogPlotZoomer(canvas),
    d_cornerWidget(NULL),
    d_hScrollData(NULL),
    d_vScrollData(NULL)
{
  if (!canvas)
    return;

  d_hScrollData = new ScrollData;
  d_vScrollData = new ScrollData;
}

ScrollZoomer::~ScrollZoomer()
{
  delete d_cornerWidget;
  delete d_vScrollData;
  delete d_hScrollData;
}

void ScrollZoomer::rescale()
{
  QwtPlotZoomer::rescale();
  updateScrollBars();
}

ScrollBar *ScrollZoomer::scrollBar(Qt::Orientation o)
{
  ScrollBar *&sb = (o == Qt::Vertical)
                   ? d_vScrollData->scrollBar : d_hScrollData->scrollBar;

  if (sb == NULL)
    {
      sb = new ScrollBar(o, canvas());
      sb->hide();
      connect(sb,
              SIGNAL(valueChanged(Qt::Orientation, double, double)),
              SLOT(scrollBarMoved(Qt::Orientation, double, double)));
    }
  return sb;
}

ScrollBar *ScrollZoomer::horizontalScrollBar() const
  {
    return d_hScrollData->scrollBar;
  }

ScrollBar *ScrollZoomer::verticalScrollBar() const
  {
    return d_vScrollData->scrollBar;
  }

#if QT_VERSION < 0x040000
void ScrollZoomer::setHScrollBarMode(Q3ScrollView::ScrollBarMode mode)
#else
void ScrollZoomer::setHScrollBarMode(Qt::ScrollBarPolicy mode)
#endif
{
  if (hScrollBarMode() != mode)
    {
      d_hScrollData->mode = mode;
      updateScrollBars();
    }
}

#if QT_VERSION < 0x040000
void ScrollZoomer::setVScrollBarMode(Q3ScrollView::ScrollBarMode mode)
#else
void ScrollZoomer::setVScrollBarMode(Qt::ScrollBarPolicy mode)
#endif
{
  if (vScrollBarMode() != mode)
    {
      d_vScrollData->mode = mode;
      updateScrollBars();
    }
}

#if QT_VERSION < 0x040000
Q3ScrollView::ScrollBarMode ScrollZoomer::hScrollBarMode() const
#else
Qt::ScrollBarPolicy ScrollZoomer::hScrollBarMode() const
#endif

  {
    return d_hScrollData->mode;
  }

#if QT_VERSION < 0x040000
Q3ScrollView::ScrollBarMode ScrollZoomer::vScrollBarMode() const
#else
Qt::ScrollBarPolicy ScrollZoomer::vScrollBarMode() const
#endif

  {
    return d_vScrollData->mode;
  }

void ScrollZoomer::setHScrollBarPosition(ScrollBarPosition pos)
{
  if (d_hScrollData->position != pos)
    {
      d_hScrollData->position = pos;
      updateScrollBars();
    }
}

void ScrollZoomer::setVScrollBarPosition(ScrollBarPosition pos)
{
  if (d_vScrollData->position != pos)
    {
      d_vScrollData->position = pos;
      updateScrollBars();
    }
}

ScrollZoomer::ScrollBarPosition ScrollZoomer::hScrollBarPosition() const
  {
    return d_hScrollData->position;
  }

ScrollZoomer::ScrollBarPosition ScrollZoomer::vScrollBarPosition() const
  {
    return d_vScrollData->position;
  }

void ScrollZoomer::setCornerWidget(QWidget *w)
{
  if (w != d_cornerWidget)
    {
      if (canvas())
        {
          delete d_cornerWidget;
          d_cornerWidget = w;
          if (d_cornerWidget->parent() != canvas())
            {
#if QT_VERSION < 0x040000
              d_cornerWidget->reparent(canvas(), QPoint(0, 0));
#else
              d_cornerWidget->setParent(canvas());
#endif
            }

          updateScrollBars();
        }
    }
}

QWidget *ScrollZoomer::cornerWidget() const
  {
    return d_cornerWidget;
  }

bool ScrollZoomer::eventFilter(QObject *o, QEvent *e)
{
  if (o == canvas())
    {
      switch (e->type())
        {
        case QEvent::Resize:
          {
            const int fw = ((QwtPlotCanvas *)canvas())->frameWidth();

            QRect rect;
            rect.setSize(((QResizeEvent *)e)->size());
            rect.setRect(rect.x() + fw, rect.y() + fw,
                         rect.width() - 2 * fw, rect.height() - 2 * fw);

            layoutScrollBars(rect);
            break;
          }
        case QEvent::ChildRemoved:
          {
            const QObject *child = ((QChildEvent *)e)->child();
            if (child == d_cornerWidget)
              d_cornerWidget = NULL;
            else if (child == d_hScrollData->scrollBar)
              d_hScrollData->scrollBar = NULL;
            else if (child == d_vScrollData->scrollBar)
              d_vScrollData->scrollBar = NULL;
            break;
          }
        default:
          break;
        }
    }
  return QwtPlotZoomer::eventFilter(o, e);
}

bool ScrollZoomer::needScrollBar(Qt::Orientation o) const
  {
#if QT_VERSION < 0x040000
    Q3ScrollView::ScrollBarMode mode;
#else
    Qt::ScrollBarPolicy mode;
#endif
    double zoomMin, zoomMax, baseMin, baseMax;

    if (o == Qt::Horizontal)
      {
        mode = d_hScrollData->mode;
        baseMin = zoomBase().left();
        baseMax = zoomBase().right();
        zoomMin = zoomRect().left();
        zoomMax = zoomRect().right();
      }
    else
      {
        mode = d_vScrollData->mode;
        baseMin = zoomBase().top();
        baseMax = zoomBase().bottom();
        zoomMin = zoomRect().top();
        zoomMax = zoomRect().bottom();
      }

    bool needed = false;
    switch (mode)
      {
#if QT_VERSION < 0x040000
      case Q3ScrollView::AlwaysOn:
#else
      case Qt::ScrollBarAlwaysOn:
#endif
        needed = true;
        break;
#if QT_VERSION < 0x040000
      case Q3ScrollView::AlwaysOff:
#else
      case Qt::ScrollBarAlwaysOff:
#endif
        needed = false;
        break;
      default:
        {
          if (baseMin < zoomMin || baseMax > zoomMax)
            needed = true;
          break;
        }
      }
    return needed;
  }

void ScrollZoomer::updateScrollBars()
{
  if (!canvas())
    return;

  const int xAxis = QwtPlotZoomer::xAxis();
  const int yAxis = QwtPlotZoomer::yAxis();

  int xScrollBarAxis = xAxis;
  if (hScrollBarPosition() == OppositeToScale)
    xScrollBarAxis = oppositeAxis(xScrollBarAxis);

  int yScrollBarAxis = yAxis;
  if (vScrollBarPosition() == OppositeToScale)
    yScrollBarAxis = oppositeAxis(yScrollBarAxis);

  QwtPlotLayout *layout = plot()->plotLayout();

  bool showHScrollBar = needScrollBar(Qt::Horizontal);
  if (showHScrollBar)
    {
      ScrollBar *sb = scrollBar(Qt::Horizontal);

      sb->setPalette(plot()->palette());

      const QwtScaleEngine *se = plot()->axisScaleEngine(xAxis);
      sb->setInverted(se->testAttribute(QwtScaleEngine::Inverted));
      sb->setLogScale(dynamic_cast<const QwtLog10ScaleEngine*>(se));

      sb->setBase(zoomBase().left(), zoomBase().right());
      sb->moveSlider(zoomRect().left(), zoomRect().right());

      if (!sb->isVisibleTo(canvas()))
        {
          sb->show();
          layout->setCanvasMargin(layout->canvasMargin(xScrollBarAxis)
                                  + sb->extent(), xScrollBarAxis);
        }
    }
  else
    {
      if (horizontalScrollBar())
        {
          horizontalScrollBar()->hide();
          layout->setCanvasMargin(layout->canvasMargin(xScrollBarAxis)
                                  - horizontalScrollBar()->extent(), xScrollBarAxis);
        }
    }

  bool showVScrollBar = needScrollBar(Qt::Vertical);
  if (showVScrollBar)
    {
      ScrollBar *sb = scrollBar(Qt::Vertical);

      sb->setPalette(plot()->palette());

      const QwtScaleEngine *se = plot()->axisScaleEngine(yAxis);
      sb->setInverted(!(se->testAttribute(QwtScaleEngine::Inverted)));
      sb->setLogScale(dynamic_cast<const QwtLog10ScaleEngine*>(se));

      sb->setBase(zoomBase().top(), zoomBase().bottom());
      sb->moveSlider(zoomRect().top(), zoomRect().bottom());

      if (!sb->isVisibleTo(canvas()))
        {
          sb->show();
          layout->setCanvasMargin(layout->canvasMargin(yScrollBarAxis)
                                  + sb->extent(), yScrollBarAxis);
        }
    }
  else
    {
      if (verticalScrollBar())
        {
          verticalScrollBar()->hide();
          layout->setCanvasMargin(layout->canvasMargin(yScrollBarAxis)
                                  - verticalScrollBar()->extent(), yScrollBarAxis);
        }
    }

  if (showHScrollBar && showVScrollBar)
    {
      if (d_cornerWidget == NULL)
        {
          d_cornerWidget = new QWidget(canvas());
          d_cornerWidget->setPalette(plot()->palette());
        }
      d_cornerWidget->show();
    }
  else
    {
      if (d_cornerWidget)
        d_cornerWidget->hide();
    }

  layoutScrollBars(((QwtPlotCanvas *)canvas())->contentsRect());
}

void ScrollZoomer::layoutScrollBars(const QRect &rect)
{
  int hPos = xAxis();
  if (hScrollBarPosition() == OppositeToScale)
    hPos = oppositeAxis(hPos);

  int vPos = yAxis();
  if (vScrollBarPosition() == OppositeToScale)
    vPos = oppositeAxis(vPos);

  ScrollBar *hScrollBar = horizontalScrollBar();
  ScrollBar *vScrollBar = verticalScrollBar();

  const int hdim = hScrollBar ? hScrollBar->extent() : 0;
  const int vdim = vScrollBar ? vScrollBar->extent() : 0;

  if (hScrollBar && hScrollBar->isVisible())
    {
      int x = rect.x();
      int y = (hPos == QwtPlot::xTop)
              ? rect.top() : rect.bottom() - hdim + 1;
      int w = rect.width();

      if (vScrollBar && vScrollBar->isVisible())
        {
          if (vPos == QwtPlot::yLeft)
            x += vdim;
          w -= vdim;
        }

      hScrollBar->setGeometry(x, y, w, hdim);
    }
  if (vScrollBar && vScrollBar->isVisible())
    {
      int pos = yAxis();
      if (vScrollBarPosition() == OppositeToScale)
        pos = oppositeAxis(pos);

      int x = (vPos == QwtPlot::yLeft)
              ? rect.left() : rect.right() - vdim + 1;
      int y = rect.y();

      int h = rect.height();

      if (hScrollBar && hScrollBar->isVisible())
        {
          if (hPos == QwtPlot::xTop)
            y += hdim;

          h -= hdim;
        }

      vScrollBar->setGeometry(x, y, vdim, h);
    }
  if (hScrollBar && hScrollBar->isVisible() &&
      vScrollBar && vScrollBar->isVisible())
    {
      if (d_cornerWidget)
        {
          QRect cornerRect(
            vScrollBar->pos().x(), hScrollBar->pos().y(),
            vdim, hdim);
          d_cornerWidget->setGeometry(cornerRect);
        }
    }
}

void ScrollZoomer::scrollBarMoved(Qt::Orientation o, double min, double)
{
  if (o == Qt::Horizontal)
    move(min, zoomRect().top());
  else
    move(zoomRect().left(), min);

  emit zoomed(zoomRect());
}

int ScrollZoomer::oppositeAxis(int axis) const
  {
    switch (axis)
      {
      case QwtPlot::xBottom:
        return QwtPlot::xTop;
      case QwtPlot::xTop:
        return QwtPlot::xBottom;
      case QwtPlot::yLeft:
        return QwtPlot::yRight;
      case QwtPlot::yRight:
        return QwtPlot::yLeft;
      default:
        break;
      }

    return axis;
  }
