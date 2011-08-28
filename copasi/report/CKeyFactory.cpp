/* Begin CVS Header
$Source: /fs/turing/cvs/copasi_dev/copasi/report/CKeyFactory.cpp,v $
$Revision: 1.20 $
$Name: Build-33 $
$Author: shoops $
$Date: 2010/02/15 18:17:01 $
End CVS Header */

// Copyright (C) 2008 by Pedro Mendes, Virginia Tech Intellectual
// Properties, Inc., EML Research, gGmbH, University of Heidelberg,
// and The University of Manchester.
// All rights reserved.

// Copyright (C) 2001 - 2007 by Pedro Mendes, Virginia Tech Intellectual
// Properties, Inc. and EML Research, gGmbH.
// All rights reserved.

/**
 * CKeyFactory class.
 * This class is used to create a unique key whithin COPASI. It also allows
 * retreival of the CCopasiObject the key is assigned to.
 *
 * Created for Copasi by Stefan Hoops 2003
 * Copyright Stefan Hoops
 */
#include <sstream>
#include <stdlib.h>

#include "copasi.h"

#include "CKeyFactory.h"

bool CKeyFactory::isValidKey(const std::string & key,
                             const std::string & prefix)
{
  if (key == "" && prefix == "") return true;

  unsigned C_INT32 digitsStart = key.length() - 1;

  while (isDigit(key[digitsStart]) && digitsStart) --digitsStart;

  if (digitsStart < 1 || digitsStart > key.length() - 2 || key[digitsStart] != '_') return false;

  if (prefix != "")
    {
      if (prefix != key.substr(0, digitsStart)) return false;
      else return true;
    }

  unsigned C_INT32 prefixEnd = 0;

  while (isPrefix(key[prefixEnd]) && prefixEnd < digitsStart) ++prefixEnd;

  return (prefixEnd == digitsStart);
}

CKeyFactory::CDecisionVector::CDecisionVector():
    CVector< bool >()
{}

CKeyFactory::CDecisionVector::CDecisionVector(const std::string & str):
    CVector< bool >(256)
{
  unsigned C_INT32 i, imax;

  for (i = 0, imax = size(); i < imax; i++)
    (*(CVector< bool > *) this)[i] = false;

  for (i = 0, imax = str.length(); i < imax; i++)
    (*(CVector< bool > *) this)[(unsigned C_INT32) str[i]] = true;
}

CKeyFactory::CDecisionVector::~CDecisionVector() {}

const bool & CKeyFactory::CDecisionVector::operator()(const unsigned char & c) const
{return (*(CVector< bool > *) this)[(unsigned C_INT32) c];}

CKeyFactory::HashTable::HashTable():
    mBeyond(0),
    mSize(128),
    mTable(128),
    mFree()
{memset(mTable.array(), 0, mSize * sizeof(CCopasiObject *));}

CKeyFactory::HashTable::HashTable(const CKeyFactory::HashTable & src):
    mBeyond(src.mBeyond),
    mSize(src.mSize),
    mTable(src.mTable),
    mFree(src.mFree)
{}

CKeyFactory::HashTable::~HashTable() {}

unsigned C_INT32 CKeyFactory::HashTable::add(CCopasiObject * pObject)
{
  unsigned C_INT32 index;

  if (!mFree.empty())
    {
      index = mFree.top();
      mFree.pop();
    }
  else
    {
      index = mBeyond;

      mBeyond++;

      if (mBeyond > mSize)
        {
          mTable.resize(mSize * 2, true);
          memset(mTable.array() + mSize, 0, mSize * sizeof(CCopasiObject *));
          mSize *= 2;
        }
    }

  mTable[index] = pObject;
  return index;
}

bool CKeyFactory::HashTable::addFix(const unsigned C_INT32 & index,
                                    CCopasiObject * pObject)
{
  while (index >= mSize)
    {
      mTable.resize(mSize * 2, true);
      memset(mTable.array() + mSize, 0,
             mSize * sizeof(CCopasiObject *));
      mSize *= 2;
    }

  if (mTable[index]) return false;

  mTable[index] = pObject;
  return true;
}

CCopasiObject * CKeyFactory::HashTable::get(const unsigned C_INT32 & index)
{
  if (index < mSize) return mTable[index];

  return NULL;
}

bool CKeyFactory::HashTable::remove(const unsigned C_INT32 & index)
{
  if (index < mSize)
    {
      if (!mTable[index]) return false;

      mTable[index] = NULL;
      mFree.push(index);

      return true;
    }

  return false;
}

CKeyFactory::CDecisionVector CKeyFactory::isDigit("0123456789");

CKeyFactory::CDecisionVector CKeyFactory::isPrefix("_abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789");

CKeyFactory::CKeyFactory():
    mKeyTable()
{}

CKeyFactory::~CKeyFactory() {}

std::string CKeyFactory::add(const std::string & prefix,
                             CCopasiObject * pObject)
{
  std::map< std::string, CKeyFactory::HashTable >::iterator it =
    mKeyTable.find(prefix);

  if (it == mKeyTable.end())
    {
      std::pair<std::map< std::string, CKeyFactory::HashTable >::iterator, bool> ret =
        mKeyTable.insert(std::map< std::string, CKeyFactory::HashTable >::value_type(prefix, CKeyFactory::HashTable()));

      it = ret.first;
    }

  std::stringstream key;
  key << prefix + "_" << it->second.add(pObject);

  return key.str();
}

bool CKeyFactory::addFix(const std::string & key, CCopasiObject * pObject)
{
  unsigned C_INT32 pos = key.length() - 1;

  while (isDigit(key[pos]) && pos) --pos;

  std::string Prefix = key.substr(0, pos);
  unsigned C_INT32 index = atoi(key.substr(pos + 1).c_str());

  std::map< std::string, CKeyFactory::HashTable >::iterator it =
    mKeyTable.find(Prefix);

  if (it == mKeyTable.end())
    {
      std::pair<std::map< std::string, CKeyFactory::HashTable >::iterator, bool> ret =
        mKeyTable.insert(std::map< std::string, CKeyFactory::HashTable >::value_type(Prefix, CKeyFactory::HashTable()));

      it = ret.first;
    }

  return it->second.addFix(index, pObject);
}

bool CKeyFactory::remove(const std::string & key)
{
  unsigned C_INT32 pos = key.length();

  if (pos == 0) return false;

  --pos;

  while (isDigit(key[pos]) && pos) --pos;

  std::string Prefix = key.substr(0, pos);

  unsigned C_INT32 index = 0;

  if (pos + 1 < key.length())
    index = atoi(key.substr(pos + 1).c_str());

  std::map< std::string, CKeyFactory::HashTable >::iterator it =
    mKeyTable.find(Prefix);

  if (it == mKeyTable.end()) return false;

  return it->second.remove(index);
}

CCopasiObject * CKeyFactory::get(const std::string & key)
{
  if (key.length() == 0) return NULL;

  unsigned C_INT32 pos = key.length() - 1; //TODO !!!pos can be invalid (-1); not anymore, but look for other errors like this

  while (isDigit(key[pos]) && pos) --pos;

  std::string Prefix = key.substr(0, pos);
  unsigned C_INT32 index = atoi(key.substr(pos + 1).c_str());

  std::map< std::string, CKeyFactory::HashTable >::iterator it =
    mKeyTable.find(Prefix);

  if (it == mKeyTable.end()) return NULL;

  return it->second.get(index);
}
