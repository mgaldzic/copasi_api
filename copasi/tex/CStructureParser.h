// Begin CVS Header
//   $Source: /fs/turing/cvs/copasi_dev/copasi/tex/CStructureParser.h,v $
//   $Revision: 1.4 $
//   $Name: Build-33 $
//   $Author: shoops $
//   $Date: 2010/07/16 19:03:27 $
// End CVS Header

// Copyright (C) 2010 by Pedro Mendes, Virginia Tech Intellectual
// Properties, Inc., University of Heidelberg, and The University
// of Manchester.
// All rights reserved.

// Copyright (C) 2008 by Pedro Mendes, Virginia Tech Intellectual
// Properties, Inc., EML Research, gGmbH, University of Heidelberg,
// and The University of Manchester.
// All rights reserved.

// written by pwilly on 09.07.08

#ifndef STRUCTUREPARSER_H
#define STRUCTUREPARSER_H

#include <QXmlAttributes>
#include <QLinkedList>

class QString;

class CStructureParser : public QXmlDefaultHandler
{
public:
  CStructureParser(int);
  //  ~CStructureParser(){};

  bool startDocument();
  bool startElement(const QString&, const QString&, const QString& ,
                    const QXmlAttributes&);
  bool characters(const QString&);
  bool endElement(const QString&, const QString&, const QString&);

  bool ignorableWhitespace(const QString&);
  bool skippedEntity(const QString&);

  QString getTeX();

private:
  QString indent;

  // for writing formula in latex form
  QString texHead, tex, texTail;

  // tags
  QString tagName;

  // for column allignment
  bool needToWriteColumnAllignment;
  int sumColumns;
  int indexColumns;

  QLinkedList<QString> mListOfUncompletedTags;
};

#endif
