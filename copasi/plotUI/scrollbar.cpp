/* Begin CVS Header
   $Source: /fs/turing/cvs/copasi_dev/copasi/plotUI/scrollbar.cpp,v $
   $Revision: 1.4 $
   $Name: Build-33 $
   $Author: shoops $
   $Date: 2006/06/20 13:19:33 $
   End CVS Header */

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

#include <qstyle.h>
#if QT_VERSION >= 0x040000
#include <qstyleoption.h>
#endif
#include "scrollbar.h"

#include <math.h>

ScrollBar::ScrollBar(QWidget * parent):
    QScrollBar(parent),
    mLogScale(false)
{
  init();
}

ScrollBar::ScrollBar(Qt::Orientation o,
                     QWidget *parent):
    QScrollBar(o, parent),
    mLogScale(false)
{
  init();
}

ScrollBar::ScrollBar(double minBase, double maxBase, bool logscale,
                     Qt::Orientation o, QWidget *parent):
    QScrollBar(o, parent),
    mLogScale(logscale)
{
  init();
  setBase(minBase, maxBase);
  moveSlider(minBase, maxBase);
}

void ScrollBar::init()
{
  d_inverted = orientation() == Qt::Vertical;
  d_baseTicks = 1000000;
  d_minBase = 0.0;
  d_maxBase = 1.0;
  moveSlider(d_minBase, d_maxBase);

  connect(this, SIGNAL(sliderMoved(int)), SLOT(catchSliderMoved(int)));
  connect(this, SIGNAL(valueChanged(int)), SLOT(catchValueChanged(int)));
}

void ScrollBar::setInverted(bool inverted)
{
  if (d_inverted != inverted)
    {
      d_inverted = inverted;
      moveSlider(minSliderValue(), maxSliderValue());
    }
}

bool ScrollBar::isInverted() const
  {
    return d_inverted;
  }

void ScrollBar::setBase(double min, double max)
{
  if (mLogScale)
    {
      min = log(min);
      max = log(max);
    }

  if (min != d_minBase || max != d_maxBase)
    {
      d_minBase = min;
      d_maxBase = max;

      moveSlider(minSliderValue(), maxSliderValue());
    }
}

void ScrollBar::moveSlider(double min, double max)
{
  if (mLogScale)
    {
      min = log(min);
      max = log(max);
    }

  int sliderTicks;
  sliderTicks = qRound((max - min) /
                       (d_maxBase - d_minBase) * d_baseTicks);

  // setRange initiates a valueChanged of the scrollbars
  // in some situations. So we block
  // and unblock the signals.

  blockSignals(true);

  setRange(sliderTicks / 2, d_baseTicks - sliderTicks / 2);
  int steps = sliderTicks / 200;
  if (steps <= 0)
    steps = 1;

#if QT_VERSION < 0x040000
  setSteps(steps, sliderTicks);
#else
  setSingleStep(steps);
  setPageStep(sliderTicks);
#endif

  int tick;
  tick = mapToTick(min + (max - min) / 2);

  if (isInverted())
    tick = d_baseTicks - tick;

#if QT_VERSION < 0x040000
  directSetValue(tick);
  rangeChange();
#else
  setSliderPosition(tick);
#endif
  blockSignals(false);
}

double ScrollBar::minBaseValue() const
  {
    if (mLogScale)
      return exp(d_minBase);
    else
      return d_minBase;
  }

double ScrollBar::maxBaseValue() const
  {
    if (mLogScale)
      return exp(d_maxBase);
    else
      return d_maxBase;
  }

void ScrollBar::sliderRange(int value, double &min, double &max) const
  {
    if (isInverted())
      value = d_baseTicks - value;

    const int visibleTicks = pageStep();

    min = mapFromTick(value - visibleTicks / 2);
    max = mapFromTick(value + visibleTicks / 2);

    if (mLogScale)
      {
        min = exp(min);
        max = exp(max);
      }
  }

double ScrollBar::minSliderValue() const
  {
    double min, dummy;
    sliderRange(value(), min, dummy);

    return min;
  }

double ScrollBar::maxSliderValue() const
  {
    double max, dummy;
    sliderRange(value(), dummy, max);

    return max;
  }

int ScrollBar::mapToTick(double v) const
  {
    return (int) ((v - d_minBase) / (d_maxBase - d_minBase) * d_baseTicks);
  }

double ScrollBar::mapFromTick(int tick) const
  {
    return d_minBase + (d_maxBase - d_minBase) * tick / d_baseTicks;
  }

void ScrollBar::catchValueChanged(int value)
{
  double min, max;
  sliderRange(value, min, max);
  emit valueChanged(orientation(), min, max);
}

void ScrollBar::catchSliderMoved(int value)
{
  double min, max;
  sliderRange(value, min, max);
  emit sliderMoved(orientation(), min, max);
}

int ScrollBar::extent() const
  {
#if QT_VERSION < 0x040000
    return style().pixelMetric(QStyle::PM_ScrollBarExtent, this);
#else
    QStyleOptionSlider opt;
    opt.init(this);
    opt.subControls = QStyle::SC_None;
    opt.activeSubControls = QStyle::SC_None;
    opt.orientation = orientation();
    opt.minimum = minimum();
    opt.maximum = maximum();
    opt.sliderPosition = sliderPosition();
    opt.sliderValue = value();
    opt.singleStep = singleStep();
    opt.pageStep = pageStep();
    opt.upsideDown = invertedAppearance();
    if (orientation() == Qt::Horizontal)
      opt.state |= QStyle::State_Horizontal;
    return style()->pixelMetric(QStyle::PM_ScrollBarExtent, &opt, this);
#endif
  }

void ScrollBar::setLogScale(bool l)
{
  mLogScale = l;
}

bool ScrollBar::isLogScale() const
  {
    return mLogScale;
  }
