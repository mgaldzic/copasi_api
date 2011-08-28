// Begin CVS Header
//   $Source: /fs/turing/cvs/copasi_dev/copasi/UI/CCopasiSelectionWidget.cpp,v $
//   $Revision: 1.15 $
//   $Name: Build-33 $
//   $Author: shoops $
//   $Date: 2009/04/21 16:20:31 $
// End CVS Header

// Copyright (C) 2008 by Pedro Mendes, Virginia Tech Intellectual
// Properties, Inc., EML Research, gGmbH, University of Heidelberg,
// and The University of Manchester.
// All rights reserved.

// Copyright (C) 2001 - 2007 by Pedro Mendes, Virginia Tech Intellectual
// Properties, Inc. and EML Research, gGmbH.
// All rights reserved.

#include <iostream>

#include "copasi.h"

#include "CCopasiSelectionWidget.h"
#include "ObjectBrowserWidget.h"
#include "CCopasiSimpleSelectionTree.h"

CCopasiSelectionWidget::CCopasiSelectionWidget(QWidget* parent , const char* name , Qt::WFlags fl): Q3WidgetStack(parent, name, fl), mpSimpleTree(NULL), mpObjectBrowser(NULL), mpOutputVector(NULL),
    mSingleSelect(false), mExpertMode(false)
{
  this->mpSimpleTree = new CCopasiSimpleSelectionTree(this);
  this->addWidget(this->mpSimpleTree);
  this->setSingleSelection(false);
  this->setExpertMode(false);
}

CCopasiSelectionWidget::~CCopasiSelectionWidget()
{
  pdelete(mpObjectBrowser);
  delete mpSimpleTree;
}

void CCopasiSelectionWidget::populateTree(const CModel * model,
    const CCopasiSimpleSelectionTree::ObjectClasses & classes)
{
  this->mpSimpleTree->populateTree(model, classes);
}

void CCopasiSelectionWidget::setOutputVector(std::vector< const CCopasiObject * > * outputVector)
{
  this->mpOutputVector = outputVector;

  if (this->mExpertMode)
    {
      this->mpObjectBrowser->setOutputVector(this->mpOutputVector);
    }
  else
    {
      this->mpSimpleTree->setOutputVector(this->mpOutputVector);
    }
}

void CCopasiSelectionWidget::setSingleSelection(bool singleSelection)
{
  if (this->mSingleSelect == singleSelection) return;

  this->mSingleSelect = singleSelection;
  this->mpSimpleTree->clearSelection();

  if (this->mSingleSelect)
    {
      this->mpSimpleTree->setSelectionMode(Q3ListView::Single);

      if (this->mpObjectBrowser)
        {
          /* this needs to be implemented first !

             this->mpObjectBrowser->setSingleSelection(true);

           */
        }
    }
  else
    {
      this->mpSimpleTree->setSelectionMode(Q3ListView::Extended);

      if (this->mpObjectBrowser)
        {
          /* this needs to be implemented first !

             this->mpObjectBrowser->setSingleSelection(true);

           */
        }
    }
}

bool CCopasiSelectionWidget::isSingleSelection() const
{
  return this->mSingleSelect;
}

void CCopasiSelectionWidget::setExpertMode(bool expertMode)
{
  if (this->mExpertMode == expertMode) return;

  this->mExpertMode = expertMode;

  if (this->mExpertMode)
    {
      if (!this->mpObjectBrowser)
        {
          // create the object browser only if needed
          // this saves some time if the model is large
          this->mpObjectBrowser = new ObjectBrowserWidget(this);
          this->addWidget(this->mpObjectBrowser);
        }

      this->mpSimpleTree->commitClicked();
      this->mpSimpleTree->setOutputVector(NULL);
      this->mpObjectBrowser->setOutputVector(this->mpOutputVector);
      this->raiseWidget(this->mpObjectBrowser);
    }
  else
    {
      // activate the SimpleTree and get the selection from the ObjectBrowser if there is one
      if (this->mpObjectBrowser)
        {
          this->mpObjectBrowser->commitClicked();
          this->mpObjectBrowser->setOutputVector(NULL);
        }

      this->mpSimpleTree->setOutputVector(this->mpOutputVector);
      this->raiseWidget(this->mpSimpleTree);
    }
}

bool CCopasiSelectionWidget::expertMode() const
{
  return this->mExpertMode;
}

void CCopasiSelectionWidget::commit()
{
  if (this->mpOutputVector)
    {
      if (this->mExpertMode)
        {
          this->mpObjectBrowser->commitClicked();

          if (this->mSingleSelect && this->mpOutputVector->size() > 1)
            {
              const CCopasiObject * object = this->mpOutputVector->at(0);
              this->mpOutputVector->clear();
              this->mpOutputVector->push_back(object);
            }
        }
      else
        {
          this->mpSimpleTree->commitClicked();
        }
    }
}
