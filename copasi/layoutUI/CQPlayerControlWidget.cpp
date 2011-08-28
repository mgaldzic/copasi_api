// Begin CVS Header
//   $Source: /fs/turing/cvs/copasi_dev/copasi/layoutUI/CQPlayerControlWidget.cpp,v $
//   $Revision: 1.10 $
//   $Name: Build-33 $
//   $Author: gauges $
//   $Date: 2010/03/10 12:33:51 $
// End CVS Header

// Copyright (C) 2010 by Pedro Mendes, Virginia Tech Intellectual
// Properties, Inc., University of Heidelberg, and The University
// of Manchester.
// All rights reserved.

// Copyright (C) 2008 by Pedro Mendes, Virginia Tech Intellectual
// Properties, Inc., EML Research, gGmbH, University of Heidelberg,
// and The University of Manchester.
// All rights reserved.

#include "CQPlayerControlWidget.h"

#include <QAction>
#include <QGridLayout>
#include <QIcon>
#include <QLayout>
#include <QPixmap>
#include <QPushButton>

#include "play.xpm"
#include "stop.xpm"
#include "pause.xpm"
#include "forward.xpm"
#include "backward.xpm"
#include "forward_single.xpm"
#include "backward_single.xpm"

CQPlayerControlWidget::CQPlayerControlWidget(QWidget* pParent): QWidget(pParent), mNumSteps(0), mCurrentStep(0), mPlaying(false)
{
  QGridLayout* pLayout = new QGridLayout;
  this->setLayout(pLayout);
  this->mpStepBackwardButton = new QPushButton;
  this->mpStepBackwardButton->setIcon(QPixmap(backward_single_xpm));
  pLayout->addWidget(this->mpStepBackwardButton, 0, 0, 1, 1);
  this->mpPlayButton = new QPushButton;
  this->mpPlayButton->setIcon(QPixmap(play_xpm));
  pLayout->addWidget(this->mpPlayButton, 0, 1, 1, 1);
  this->mpStepForwardButton = new QPushButton;
  pLayout->addWidget(this->mpStepForwardButton, 0, 2, 1, 1);
  this->mpStepForwardButton->setIcon(QPixmap(forward_single_xpm));
  this->mpBackwardButton = new QPushButton;
  this->mpBackwardButton->setIcon(QPixmap(backward_xpm));
  pLayout->addWidget(this->mpBackwardButton, 1, 0, 1, 1);
  this->mpStopButton = new QPushButton;
  this->mpStopButton->setIcon(QPixmap(stop_xpm));
  pLayout->addWidget(this->mpStopButton, 1, 1, 1, 1);
  this->mpForwardButton = new QPushButton;
  this->mpForwardButton->setIcon(QPixmap(forward_xpm));
  pLayout->addWidget(this->mpForwardButton, 1, 2, 1, 1);
  this->createActions();
  this->updateButtons();
  this->updateActions();
  this->setTabOrder(this->mpPlayButton, this->mpStopButton);
  this->setTabOrder(this->mpStopButton, this->mpStepForwardButton);
  this->setTabOrder(this->mpStepForwardButton, this->mpStepBackwardButton);
  this->setTabOrder(this->mpStepBackwardButton, this->mpForwardButton);
  this->setTabOrder(this->mpForwardButton, this->mpBackwardButton);
  connect(this->mpPlayButton, SIGNAL(clicked()), this, SLOT(slot_play_clicked()));
  connect(this->mpStopButton, SIGNAL(clicked()), this, SLOT(slot_stop_clicked()));
  connect(this->mpForwardButton, SIGNAL(clicked()), this, SLOT(slot_forward_clicked()));
  connect(this->mpBackwardButton, SIGNAL(clicked()), this, SLOT(slot_backward_clicked()));
  connect(this->mpStepForwardButton, SIGNAL(clicked()), this, SLOT(slot_step_forward_clicked()));
  connect(this->mpStepBackwardButton, SIGNAL(clicked()), this, SLOT(slot_step_backward_clicked()));
}

void CQPlayerControlWidget::createActions()
{
  this->mpPlayAction = new QAction(QPixmap(play_xpm), "play", this);
  connect(this->mpPlayAction, SIGNAL(triggered()), this, SLOT(slot_play_clicked()));
  this->mpPauseAction = new QAction(QPixmap(pause_xpm), "pause", this);
  connect(this->mpPauseAction, SIGNAL(triggered()), this, SLOT(slot_pause_clicked()));
  this->mpStopAction = new QAction(QPixmap(stop_xpm), "stop", this);
  connect(this->mpStopAction, SIGNAL(triggered()), this, SLOT(slot_stop_clicked()));
  this->mpForwardAction = new QAction(QPixmap(forward_xpm), "forward", this);
  connect(this->mpForwardAction, SIGNAL(triggered()), this, SLOT(slot_forward_clicked()));
  this->mpBackwardAction = new QAction(QPixmap(backward_xpm), "backward", this);
  connect(this->mpBackwardAction, SIGNAL(triggered()), this, SLOT(slot_backward_clicked()));
  this->mpStepForwardAction = new QAction(QPixmap(forward_single_xpm), "step forward", this);
  connect(this->mpStepForwardAction, SIGNAL(triggered()), this, SLOT(slot_step_forward_clicked()));
  this->mpStepBackwardAction = new QAction(QPixmap(backward_single_xpm), "step backward", this);
  connect(this->mpStepBackwardAction, SIGNAL(triggered()), this, SLOT(slot_step_backward_clicked()));
}

void CQPlayerControlWidget::setNumSteps(unsigned int numSteps)
{
  this->mNumSteps = numSteps;
  this->updateButtons();
  this->updateActions();
}

void CQPlayerControlWidget::setCurrentStep(unsigned int currentStep)
{
  if (currentStep < this->mNumSteps)
    {
      this->mCurrentStep = currentStep;
    }
  else
    {
      this->mCurrentStep = this->mNumSteps - 1;
    }

  this->updateButtons();
  this->updateActions();
}

void CQPlayerControlWidget::slot_pause_clicked()
{
  this->mpPlayButton->setIcon(QPixmap(play_xpm));
  this->mPlaying = !this->mPlaying;
  this->mpPauseAction->setEnabled(FALSE);
  emit pause();
}

void CQPlayerControlWidget::slot_play_clicked()
{
  if (this->mPlaying)
    {
      this->slot_pause_clicked();
    }
  else
    {
      this->mpPlayButton->setIcon(QPixmap(pause_xpm));
      this->mPlaying = !this->mPlaying;
      this->mpPauseAction->setEnabled(TRUE);
      emit play();
    }
}

void CQPlayerControlWidget::slot_stop_clicked()
{
  this->mCurrentStep = 0;
  this->mPlaying = false;
  this->updateButtons();
  this->updateActions();
  emit stop();
}

void CQPlayerControlWidget::slot_forward_clicked()
{
  this->mCurrentStep = this->mNumSteps - 1;
  this->updateButtons();
  this->updateActions();
  emit forward();
}

void CQPlayerControlWidget::updateActions()
{
  if (this->mNumSteps < 2)
    {
      // disable all actionsw
      if (this->mpPlayAction->isEnabled())
        {
          this->mpPlayAction->setEnabled(FALSE);
        }

      if (this->mpPauseAction->isEnabled())
        {
          this->mpPauseAction->setEnabled(FALSE);
        }

      if (this->mpStopAction->isEnabled())
        {
          this->mpStopAction->setEnabled(FALSE);
        }

      if (this->mpForwardAction->isEnabled())
        {
          this->mpForwardAction->setEnabled(FALSE);
        }

      if (this->mpStepForwardAction->isEnabled())
        {
          this->mpStepForwardAction->setEnabled(FALSE);
        }

      if (this->mpBackwardAction->isEnabled())
        {
          this->mpBackwardAction->setEnabled(FALSE);
        }

      if (this->mpStepBackwardAction->isEnabled())
        {
          this->mpStepBackwardAction->setEnabled(FALSE);
        }
    }
  else
    {
      if (this->mCurrentStep == 0)
        {
          // make sure backward and step_backward are disabled and all other
          // actions are enabled
          if (!this->mPlaying)
            {
              this->mpPauseAction->setEnabled(FALSE);
            }
          else
            {
              this->mpPauseAction->setEnabled(TRUE);
            }

          if (this->mpBackwardAction->isEnabled())
            {
              this->mpBackwardAction->setEnabled(FALSE);
            }

          if (this->mpStepBackwardAction->isEnabled())
            {
              this->mpStepBackwardAction->setEnabled(FALSE);
            }

          if (this->mpPlayAction->isEnabled() == FALSE)
            {
              this->mpPlayAction->setEnabled(TRUE);
            }

          if (this->mpStopAction->isEnabled() == FALSE)
            {
              this->mpStopAction->setEnabled(TRUE);
            }

          if (this->mpForwardAction->isEnabled() == FALSE)
            {
              this->mpForwardAction->setEnabled(TRUE);
            }

          if (this->mpStepForwardAction->isEnabled() == FALSE)
            {
              this->mpStepForwardAction->setEnabled(TRUE);
            }
        }
      else if (this->mCurrentStep == this->mNumSteps - 1)
        {
          // make sure play, forward and step_forward are disabled and all other
          // actions are enabled
          // reset the icon on the play button
          if (this->mpForwardAction->isEnabled())
            {
              this->mpForwardAction->setEnabled(FALSE);
            }

          if (this->mpStepForwardAction->isEnabled())
            {
              this->mpStepForwardAction->setEnabled(FALSE);
            }

          if (this->mpPlayAction->isEnabled())
            {
              this->mpPlayAction->setEnabled(FALSE);
            }

          if (this->mpPauseAction->isEnabled())
            {
              this->mpPauseAction->setEnabled(FALSE);
            }

          if (this->mpBackwardAction->isEnabled() == FALSE)
            {
              this->mpBackwardAction->setEnabled(TRUE);
            }

          if (this->mpStepBackwardAction->isEnabled() == FALSE)
            {
              this->mpStepBackwardAction->setEnabled(TRUE);
            }

          if (this->mpStopAction->isEnabled() == FALSE)
            {
              this->mpStopAction->setEnabled(TRUE);
            }
        }
      else
        {
          // make sure all buttons are enabled
          if (this->mpPlayAction->isEnabled() == FALSE)
            {
              this->mpPlayAction->setEnabled(TRUE);
            }

          if (this->mpStopAction->isEnabled() == FALSE)
            {
              this->mpStopAction->setEnabled(TRUE);
            }

          if (this->mpForwardAction->isEnabled() == FALSE)
            {
              this->mpForwardAction->setEnabled(TRUE);
            }

          if (this->mpBackwardAction->isEnabled() == FALSE)
            {
              this->mpBackwardAction->setEnabled(TRUE);
            }

          if (this->mpStepBackwardAction->isEnabled() == FALSE)
            {
              this->mpStepBackwardAction->setEnabled(TRUE);
            }

          if (this->mpStepForwardAction->isEnabled() == FALSE)
            {
              this->mpStepForwardAction->setEnabled(TRUE);
            }
        }
    }
}

void CQPlayerControlWidget::updateButtons()
{
  if (this->mNumSteps < 2)
    {
      // disable all buttons
      if (this->mpPlayButton->isEnabled())
        {
          this->mpPlayButton->setEnabled(FALSE);

          if (this->mPlaying)
            {
              // reset the icon
              this->mpPlayButton->setIcon(QPixmap(play_xpm));
            }
        }

      if (this->mpStopButton->isEnabled())
        {
          this->mpStopButton->setEnabled(FALSE);
        }

      if (this->mpForwardButton->isEnabled())
        {
          this->mpForwardButton->setEnabled(FALSE);
        }

      if (this->mpStepForwardButton->isEnabled())
        {
          this->mpStepForwardButton->setEnabled(FALSE);
        }

      if (this->mpBackwardButton->isEnabled())
        {
          this->mpBackwardButton->setEnabled(FALSE);
        }

      if (this->mpStepBackwardButton->isEnabled())
        {
          this->mpStepBackwardButton->setEnabled(FALSE);
        }
    }
  else
    {
      if (this->mCurrentStep == 0)
        {
          // make sure backward and step_backward are disabled and all other
          // buttons are enabled
          if (!this->mPlaying)
            {
              this->mpPlayButton->setIcon(QPixmap(play_xpm));
            }
          else
            {
              this->mpPlayButton->setIcon(QPixmap(pause_xpm));
            }

          if (this->mpBackwardButton->isEnabled())
            {
              this->mpBackwardButton->setEnabled(FALSE);
            }

          if (this->mpStepBackwardButton->isEnabled())
            {
              this->mpStepBackwardButton->setEnabled(FALSE);
            }

          if (this->mpPlayButton->isEnabled() == FALSE)
            {
              this->mpPlayButton->setEnabled(TRUE);
            }

          if (this->mpStopButton->isEnabled() == FALSE)
            {
              this->mpStopButton->setEnabled(TRUE);
            }

          if (this->mpForwardButton->isEnabled() == FALSE)
            {
              this->mpForwardButton->setEnabled(TRUE);
            }

          if (this->mpStepForwardButton->isEnabled() == FALSE)
            {
              this->mpStepForwardButton->setEnabled(TRUE);
            }
        }
      else if (this->mCurrentStep == this->mNumSteps - 1)
        {
          // make sure play, forward and step_forward are disabled and all other
          // buttons are enabled
          // reset the icon on the play button
          if (this->mpForwardButton->isEnabled())
            {
              this->mpForwardButton->setEnabled(FALSE);
            }

          if (this->mpStepForwardButton->isEnabled())
            {
              this->mpStepForwardButton->setEnabled(FALSE);
            }

          if (this->mpPlayButton->isEnabled())
            {
              this->mpPlayButton->setIcon(QPixmap(play_xpm));
              this->mpPlayButton->setEnabled(FALSE);
            }

          if (this->mpBackwardButton->isEnabled() == FALSE)
            {
              this->mpBackwardButton->setEnabled(TRUE);
            }

          if (this->mpStepBackwardButton->isEnabled() == FALSE)
            {
              this->mpStepBackwardButton->setEnabled(TRUE);
            }

          if (this->mpStopButton->isEnabled() == FALSE)
            {
              this->mpStopButton->setEnabled(TRUE);
            }
        }
      else
        {
          // make sure all buttons are enabled
          if (this->mpPlayButton->isEnabled() == FALSE)
            {
              this->mpPlayButton->setEnabled(TRUE);
            }

          if (this->mpStopButton->isEnabled() == FALSE)
            {
              this->mpStopButton->setEnabled(TRUE);
            }

          if (this->mpForwardButton->isEnabled() == FALSE)
            {
              this->mpForwardButton->setEnabled(TRUE);
            }

          if (this->mpBackwardButton->isEnabled() == FALSE)
            {
              this->mpBackwardButton->setEnabled(TRUE);
            }

          if (this->mpStepBackwardButton->isEnabled() == FALSE)
            {
              this->mpStepBackwardButton->setEnabled(TRUE);
            }

          if (this->mpStepForwardButton->isEnabled() == FALSE)
            {
              this->mpStepForwardButton->setEnabled(TRUE);
            }
        }
    }
}

void CQPlayerControlWidget::slot_backward_clicked()
{
  this->mCurrentStep = 0;
  this->updateButtons();
  this->updateActions();
  emit backward();
}

void CQPlayerControlWidget::slot_step_forward_clicked()
{
  if (this->mCurrentStep != this->mNumSteps - 1)
    {
      ++this->mCurrentStep;
    }

  this->updateButtons();
  this->updateActions();
  emit step_forward();
}

void CQPlayerControlWidget::slot_step_backward_clicked()
{
  if (this->mCurrentStep != 0)
    {
      --this->mCurrentStep;
    }

  this->updateButtons();
  this->updateActions();
  emit step_backward();
}

bool CQPlayerControlWidget::isPlaying() const
{
  return this->mPlaying;
}

QAction* CQPlayerControlWidget::getPlayAction()
{
  return this->mpPlayAction;
}

QAction* CQPlayerControlWidget::getPauseAction()
{
  return this->mpPauseAction;
}

QAction* CQPlayerControlWidget::getStopAction()
{
  return this->mpStopAction;
}

QAction* CQPlayerControlWidget::getForwardAction()
{
  return this->mpForwardAction;
}

QAction* CQPlayerControlWidget::getBackwardAction()
{
  return this->mpBackwardAction;
}

QAction* CQPlayerControlWidget::getStepForwardAction()
{
  return this->mpStepForwardAction;
}

QAction* CQPlayerControlWidget::getStepBackwardAction()
{
  return this->mpStepBackwardAction;
}
