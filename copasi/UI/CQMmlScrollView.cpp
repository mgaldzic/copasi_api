// Begin CVS Header
//   $Source: /fs/turing/cvs/copasi_dev/copasi/UI/CQMmlScrollView.cpp,v $
//   $Revision: 1.8.4.1 $
//   $Name: Build-33 $
//   $Author: shoops $
//   $Date: 2010/09/29 17:02:47 $
// End CVS Header

// Copyright (C) 2010 by Pedro Mendes, Virginia Tech Intellectual
// Properties, Inc., University of Heidelberg, and The University
// of Manchester.
// All rights reserved.

// Copyright (C) 2008 by Pedro Mendes, Virginia Tech Intellectual
// Properties, Inc., EML Research, gGmbH, University of Heidelberg,
// and The University of Manchester.
// All rights reserved.

/*!
 \file CQMmlScrollView.cpp
 \brief Implementation file of class CQMmlScrollView
 */

#include "CQMmlScrollView.h"

#include <QApplication>

#ifdef HAVE_MML
# include "mml/qtmmlwidget.h"
#endif // HAVE_MML

#include "qtUtilities.h"

///  Constructor
CQMmlScrollView::CQMmlScrollView(QWidget* parent, const char* /* name */, Qt::WFlags /* fl */)
    : QScrollArea(parent)
{
  setSizePolicy(QSizePolicy::Preferred, QSizePolicy::Preferred);

#ifdef HAVE_MML
  mpMmlWidget = new QtMmlWidget(viewport() /*, "mpMmlWidget"*/);
  mpMmlWidget->setMinimumSize(QSize(0, 0));
  setWidget(mpMmlWidget);
#endif // HAVE_MML

}

/// Destructor
CQMmlScrollView::~CQMmlScrollView()
{
  // no need to delete child widgets, Qt does it all for us
}

/// Function to update the MathML widget
void CQMmlScrollView::updateWidget(std::ostringstream &mml)
{
#ifdef HAVE_MML
  mpMmlWidget->setContent(FROM_UTF8(mml.str()));
  mpMmlWidget->setBaseFontPointSize(qApp->font().pointSize());
  mpMmlWidget->setFontName(QtMmlWidget::NormalFont, qApp->font().family());

  mpMmlWidget->resize(mpMmlWidget->sizeHint().width(), mpMmlWidget->sizeHint().height());

#endif // HAVE_MML
}
