// Begin CVS Header
//   $Source: /fs/turing/cvs/copasi_dev/copasi/elementaryFluxModes/CZeroSet.cpp,v $
//   $Revision: 1.7 $
//   $Name: Build-33 $
//   $Author: shoops $
//   $Date: 2010/01/29 21:59:25 $
// End CVS Header

// Copyright (C) 2008 by Pedro Mendes, Virginia Tech Intellectual
// Properties, Inc., EML Research, gGmbH, University of Heidelberg,
// and The University of Manchester.
// All rights reserved.

#include "copasi.h"

#include "CZeroSet.h"
#include "CStepMatrixColumn.h"

CZeroSet::CIndex::CIndex(const size_t & index):
    mIndex(index / (CHAR_BIT * sizeof(unsigned C_INT32))),
    mBit(1 << (index % (CHAR_BIT * sizeof(unsigned C_INT32)))),
    mNotBit(C_INVALID_INDEX - mBit)
{}

CZeroSet::CIndex::CIndex(const CZeroSet::CIndex & src):
    mIndex(src.mIndex),
    mBit(src.mBit),
    mNotBit(src.mNotBit)
{}

CZeroSet::CIndex::~CIndex()
{}

CZeroSet::CIndex & CZeroSet::CIndex::operator ++ ()
{
  mBit <<= 1;
  mNotBit = C_INVALID_INDEX - mBit;

  if (mBit == 0)
    {
      mIndex++;
      mBit = 1;
      mNotBit = C_INVALID_INDEX;
    }

  return *this;
}

CZeroSet::CIndex & CZeroSet::CIndex::operator = (const CZeroSet::CIndex & rhs)
{
  mIndex = rhs.mIndex;
  mBit = rhs.mBit;
  mNotBit = rhs.mNotBit;

  return *this;
}

CZeroSet::CIndex & CZeroSet::CIndex::operator = (const size_t & index)
{
  mIndex = index / (CHAR_BIT * sizeof(unsigned C_INT32));
  mBit = 1 << (index % (CHAR_BIT * sizeof(unsigned C_INT32)));
  mNotBit = C_INVALID_INDEX - mBit;

  return *this;
}

bool CZeroSet::CIndex::operator < (const CZeroSet::CIndex & rhs) const
{
  if (mIndex == rhs.mIndex)
    return mBit < rhs.mBit;

  return mIndex < rhs.mIndex;
}

CZeroSet::CZeroSet(const size_t & size):
    mBitSet(size / (CHAR_BIT * sizeof(unsigned C_INT32)) + 1),
    mIgnoredBits(mBitSet.size() * CHAR_BIT * sizeof(unsigned C_INT32) - size),
    mNumberSetBits(size)
{
  mBitSet = C_INVALID_INDEX;
}

CZeroSet::CZeroSet(const CZeroSet & src):
    mBitSet(src.mBitSet),
    mIgnoredBits(src.mIgnoredBits),
    mNumberSetBits(src.mNumberSetBits)
{}

CZeroSet::~CZeroSet()
{}

bool CZeroSet::isExtremeRay(const std::vector< CStepMatrixColumn * > & columns) const
{
  std::vector< CStepMatrixColumn * >::const_iterator it = columns.begin();
  std::vector< CStepMatrixColumn * >::const_iterator end = columns.end();

  for (; it != end; ++it)
    {
      if (*it != NULL &&
          (*it)->getZeroSet() >= *this)
        {
          return false;
        }
    }

  return true;
}

std::ostream & operator << (std::ostream & os, const CZeroSet & set)
{
  const unsigned C_INT32 * pIt = set.mBitSet.array();
  const unsigned C_INT32 * pEnd = pIt + set.mBitSet.size();

  size_t CurrentBit = 0;
  size_t LastBit = set.mBitSet.size() * CHAR_BIT * sizeof(unsigned C_INT32) - set.mIgnoredBits;

  os << ' ';

  for (; pIt != pEnd; ++pIt)
    {
      for (size_t i = 0;
           i < CHAR_BIT * sizeof(unsigned C_INT32) && CurrentBit < LastBit;
           ++i, ++CurrentBit)
        {
          if ((*pIt >> i) & 1)
            {
              os << 1;
            }
          else
            {
              os << 0;
            }
        }

      os << ' ';
    }

  os << " " << set.mNumberSetBits << std::endl;

  return os;
}
