/* Begin CVS Header
$Source: /fs/turing/cvs/copasi_dev/copasi/utilities/CTableCell.cpp,v $
$Revision: 1.19 $
$Name: Build-33 $
$Author: shoops $
$Date: 2010/07/16 19:06:32 $
End CVS Header */

// Copyright (C) 2008 by Pedro Mendes, Virginia Tech Intellectual
// Properties, Inc., EML Research, gGmbH, University of Heidelberg,
// and The University of Manchester.
// All rights reserved.

// Copyright (C) 2001 - 2007 by Pedro Mendes, Virginia Tech Intellectual
// Properties, Inc. and EML Research, gGmbH.
// All rights reserved.

#include <limits>
#include <iostream>
#include <sstream>
#include <float.h>
#include <string.h>
#include <stdlib.h>

#include "copasi.h"

#include "CTableCell.h"

#include "utilities/utility.h"

CTableCell::CTableCell(const char & separator):
    mSeparator(separator),
    mName(""),
    mValue(std::numeric_limits<C_FLOAT64>::quiet_NaN()),
    mIsValue(false),
    mIsEmpty(true)
{}

CTableCell::CTableCell(const CTableCell & src):
    mSeparator(src.mSeparator),
    mName(src.mName),
    mValue(src.mValue),
    mIsValue(src.mIsValue),
    mIsEmpty(src.mIsEmpty)
{}

CTableCell::~CTableCell() {}

bool CTableCell::setSeparator(const char & separator)
{
  mSeparator = separator;
  return true;
}

const char & CTableCell::getSeparator() const {return mSeparator;}

const bool & CTableCell::isValue() const {return mIsValue;}

const bool & CTableCell::isEmpty() const {return mIsEmpty;}

const std::string & CTableCell::getName() const {return mName;}

const C_FLOAT64 & CTableCell::getValue() const {return mValue;}

std::istream & operator >> (std::istream &is, CTableCell & cell)
{
  static char buffer[256];

  cell.mName = "";

  do
    {
      is.clear();
      is.getline(buffer, 256, cell.mSeparator);
      cell.mName += buffer;
    }
  while (strlen(buffer) == 255 && !is.eof());

  /* Trim leading and trailing whitespaces from the string */
  std::string::size_type begin = cell.mName.find_first_not_of("\x20\x09\x0d\x0a");

  if (begin == std::string::npos)
    {
      cell.mName = "";
      cell.mIsValue = false;
      cell.mValue = std::numeric_limits<C_FLOAT64>::quiet_NaN();
      cell.mIsEmpty = true;

      return is;
    }

  std::string::size_type end = cell.mName.find_last_not_of("\x20\x09\x0d\x0a");

  if (end == std::string::npos)
    cell.mName = cell.mName.substr(begin);
  else
    cell.mName = cell.mName.substr(begin, end - begin + 1);

  cell.mIsEmpty = false;

  /* Try to convert the string into a number */
  const char * Tail = NULL;
  cell.mValue = strToDouble(cell.mName.c_str(), & Tail);

  if (Tail != NULL && *Tail == 0x0)
    {
      cell.mIsValue = true;
    }
  else if (cell.mName == "INF")
    {
      cell.mIsValue = true;
      cell.mValue = std::numeric_limits<C_FLOAT64>::infinity();
    }
  else if (cell.mName == "-INF")
    {
      cell.mIsValue = true;
      cell.mValue = - std::numeric_limits<C_FLOAT64>::infinity();
    }
  else
    {
      cell.mIsValue = false;
      cell.mValue = std::numeric_limits<C_FLOAT64>::quiet_NaN();
    }

  return is;
}

CTableRow::CTableRow(const unsigned C_INT32 & size,
                     const char & separator):
    mCells(0),
    mSeparator(separator),
    mIsEmpty(true)
{resize(size);}

CTableRow::CTableRow(const CTableRow & src):
    mCells(src.mCells),
    mSeparator(src.mSeparator),
    mIsEmpty(src.mIsEmpty)
{}

CTableRow::~CTableRow() {}

const std::vector< CTableCell > & CTableRow::getCells() const
{return mCells;}

bool CTableRow::resize(const unsigned C_INT32 & size)
{
  mCells.resize(size);

  std::vector< CTableCell >::iterator it = mCells.begin();
  std::vector< CTableCell >::iterator end = mCells.end();

  for (; it != end; ++it)
    it->setSeparator(mSeparator);

  return true;
}

unsigned C_INT32 CTableRow::size() const
{return mCells.size();}

const unsigned C_INT32 & CTableRow::getLastFilledCell() const
{return mLastFilledCell;}

unsigned C_INT32 CTableRow::guessColumnNumber(std::istream &is,
    const bool & rewind)
{
  std::istream::pos_type pos;

  if (rewind) pos = is.tellg();

  is >> *this;

  if (rewind) is.seekg(pos);

  unsigned C_INT32 count;

  for (count = mCells.size() - 1; count != C_INVALID_INDEX; count--)
    if (!mCells[count].isEmpty()) break;

  return count + 1;
}

const bool & CTableRow::isEmpty() const {return mIsEmpty;}

std::istream & CTableRow::readLine(std::istream & is)
{
  // Clear the line;
  std::stringstream line;

  char c;

  for (is.get(c); c != 0x0a && c != 0x0d; is.get(c))
    {
      if (is.fail() || is.eof()) break;

      line.put(c);
    }

  // Eat additional line break characters appearing on DOS and Mac text format;
  if ((c == 0x0d && is.peek() == 0x0a) || // DOS
      (c == 0x0a && is.peek() == 0x0d))   // Mac
    is.ignore(1);

  mIsEmpty = true;
  mLastFilledCell = C_INVALID_INDEX;

  std::vector< CTableCell >::iterator it = mCells.begin();
  std::vector< CTableCell >::iterator end = mCells.end();

  unsigned C_INT count;

  for (count = 0; it != end && !line.fail(); ++it, ++count)
    {
      line >> *it;

      if (!it->isEmpty())
        {
          mIsEmpty = false;
          mLastFilledCell = count;
        }
    }

  bool Finished = false;

  if (it == end)
    Finished = true;

  CTableCell Unread(mSeparator);

  while (!line.fail() && !line.eof())
    {
      mCells.push_back(Unread);
      line >> mCells.back();

      if (!mCells.back().isEmpty())
        {
          mIsEmpty = false;
          mLastFilledCell = count;
        }

      count++;
    }

  if (!Finished)
    {
      // Missing columns are filled with default
      for (; it != end; ++it)
        *it = Unread;
    }

  return is;
}

std::istream & operator >> (std::istream &is, CTableRow & row)
{return row.readLine(is);}
