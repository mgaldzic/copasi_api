// Begin CVS Header
//   $Source: /fs/turing/cvs/copasi_dev/copasi/layoutUI/CSimSummaryInfo.cpp,v $
//   $Revision: 1.2 $
//   $Name: Build-33 $
//   $Author: urost $
//   $Date: 2007/05/11 10:39:29 $
// End CVS Header

// Copyright (C) 2007 by Pedro Mendes, Virginia Tech Intellectual
// Properties, Inc. and EML Research, gGmbH.
// All rights reserved.

#include "CSimSummaryInfo.h"

CSimSummaryInfo::CSimSummaryInfo()
{
  mNumberOfSteps = 0;
  mNumberOfReactants = 0;
  mTotalTime = 0.0;
  mMaxOverallConcentration = 0.0;
}

CSimSummaryInfo::CSimSummaryInfo(C_INT32 steps, C_INT32 numReactants, C_FLOAT32 time)
{
  mNumberOfSteps = steps;
  mNumberOfReactants = numReactants;
  mTotalTime = time;
  mMaxOverallConcentration = 0.0;
}

void CSimSummaryInfo::storeMax(std::string ndKey, C_FLOAT64 maxR)
{
  maxMap.insert(std::pair<std::string, C_FLOAT64>
                (ndKey, maxR));
}

void CSimSummaryInfo::storeMin(std::string ndKey, C_FLOAT64 minR)
{
  minMap.insert(std::pair<std::string, C_FLOAT64>
                (ndKey, minR));
}

C_FLOAT64 CSimSummaryInfo::getMaxForSpecies(std::string ndKey)
{
  std::map<std::string, C_FLOAT64>::const_iterator iter = maxMap.find(ndKey);
  if (iter != maxMap.end())
    return (*iter).second;
  else
    return 0.0;
}

C_FLOAT64 CSimSummaryInfo::getMinForSpecies(std::string ndKey)
{
  std::map<std::string, C_FLOAT64>::const_iterator iter = minMap.find(ndKey);
  if (iter != minMap.end())
    return (*iter).second;
  else
    return 0.0;
}

std::ostream & operator<<(std::ostream &os, const CSimSummaryInfo & inf)
{
  os << "number of steps: " << inf.mNumberOfSteps << std::endl;
  os << "number of reactants: " << inf.mNumberOfReactants << std::endl;
  os << "total time: " << inf.mTotalTime << std::endl;
  os << "maximum overall concentration: " << inf.mMaxOverallConcentration << std::endl;
  os << "range  of values for species: " << std::endl;

  std::map<std::string, C_FLOAT64>::const_iterator iter = inf.maxMap.begin();
  std::map<std::string, C_FLOAT64>::const_iterator iter2 = inf.minMap.begin();
  while ((iter != inf.maxMap.end()) && (iter2 != inf.minMap.end()))
    {
      os << " node: " << iter->first << " : [" << iter2->second;
      os << " - " << iter->second << "]" << std::endl;
      iter2++;
      iter++;
    }

  return os;
}
