// Begin CVS Header
//   $Source: /fs/turing/cvs/copasi_dev/copasi/report/COutputAssistant.cpp,v $
//   $Revision: 1.20 $
//   $Name: Build-33 $
//   $Author: shoops $
//   $Date: 2009/10/27 16:52:48 $
// End CVS Header

// Copyright (C) 2008 by Pedro Mendes, Virginia Tech Intellectual
// Properties, Inc., EML Research, gGmbH, University of Heidelberg,
// and The University of Manchester.
// All rights reserved.

// Copyright (C) 2001 - 2007 by Pedro Mendes, Virginia Tech Intellectual
// Properties, Inc. and EML Research, gGmbH.
// All rights reserved.

#include "COutputAssistant.h"
#include "CCopasiObject.h"
#include "CReportDefinition.h"
#include "CReportDefinitionVector.h"

#include "utilities/CCopasiProblem.h"
#include "utilities/CCopasiTask.h"
#include "trajectory/CTrajectoryProblem.h"
#include "steadystate/CSteadyStateProblem.h"
#include "model/CObjectLists.h"
#include "model/CModel.h"
#include "CopasiDataModel/CCopasiDataModel.h"
#include "report/CCopasiRootContainer.h"
#include "plot/COutputDefinitionVector.h"
#include "parameterFitting/CFitProblem.h"
#include "parameterFitting/CExperimentSet.h"
#include "parameterFitting/CExperiment.h"

//******* COutputAssistant **********************************

//static member variables
COutputAssistant::Map COutputAssistant::mMap;
const std::string COutputAssistant::emptyString("");
const CDefaultOutputDescription COutputAssistant::emptyItem;

//static
std::vector<C_INT32> COutputAssistant::getListOfDefaultOutputDescriptions(const CCopasiProblem * problem)
{
  //initializes the map on first call only
  initialize();

  std::vector<C_INT32> ret;

  problem = NULL; //DEBUG only!!!

  if (!problem) //generate full list
    {
      Map::const_iterator it, itEnd = mMap.end();

      for (it = mMap.begin(); it != itEnd; ++it)
        ret.push_back(it->first);

      return ret;
    }

  const CModel* pModel = problem->getModel();

  if (!pModel) return ret;

  //TODO use CObjectLists::existsFixedMetab()

  const CTrajectoryProblem* tp = dynamic_cast<const CTrajectoryProblem*>(problem);

  if (tp)
    {
      return ret;
    }

  const CSteadyStateProblem* ssp = dynamic_cast<const CSteadyStateProblem*>(problem);

  if (ssp)
    {
      return ret;
    }

  return ret;
}

//static
C_INT32 COutputAssistant::getDefaultReportIndex(const CCopasiProblem * problem)
{
  if (!problem) return - 1;

  switch (problem->getType())
    {
      case CCopasiTask::steadyState:
        return 1000;
      case CCopasiTask::timeCourse:
        return 1000;
      default:
        return - 1;
    }
}

//       steadyState = 0,
//       timeCourse,
//       scan,
//       fluxMode,
//       optimization,
//       parameterFitting,
//       mca,
//       lyap,
// #ifdef COPASI_DEBUG
//       tss,
// #endif // COPASI_DEBUG
//       sens,
// #ifdef COPASI_SSA
//       ssa,
// #endif // COPASI_SSA
//       unset,

//static
C_INT32 COutputAssistant::getDefaultPlotIndex(const CCopasiProblem * problem)
{
  if (!problem) return - 1;

  switch (problem->getType())
    {
      case CCopasiTask::steadyState:
        return 0;
      case CCopasiTask::timeCourse:
        return 0;
      default:
        return - 1;
    }
}

//static
const std::string & COutputAssistant::getItemName(C_INT32 id)
{
  Map::const_iterator it = mMap.find(id);

  if (it == mMap.end())
    return emptyString;
  else
    return it->second.name;
}

//static
const CDefaultOutputDescription & COutputAssistant::getItem(C_INT32 id)
{
  Map::const_iterator it = mMap.find(id);

  if (it == mMap.end())
    return emptyItem;
  else
    return it->second;
}

//static
bool COutputAssistant::initialize()
{
  //if map is already constructed do nothing
  if (mMap.size()) return true;

  std::pair<C_INT32, CDefaultOutputDescription> tmp;

  //first the plots
  tmp.first = -1;
  tmp.second.name = "-- Plots";
  tmp.second.description = "";
  tmp.second.isPlot = true;
  tmp.second.mTaskType = CCopasiTask::unset;
  mMap.insert(tmp);

  //concentrations plot
  tmp.first = 0;
  tmp.second.name = "Concentrations, Volumes, and Global Quantity Values";
  tmp.second.description = "A plot of the variable species concentrations, variable compartment volumes, and variable global quantity values vs. time.";
  tmp.second.isPlot = true;
  tmp.second.mTaskType = CCopasiTask::timeCourse;
  mMap.insert(tmp);

  //particle numbers plot
  tmp.first = 1;
  tmp.second.name = "Particle Numbers, Volumes, and Global Quantity Values";
  tmp.second.description = "A plot of the variable species particle numbers, variable compartment volumes, and variable global quantity values vs. time.";
  tmp.second.isPlot = true;
  tmp.second.mTaskType = CCopasiTask::timeCourse;
  mMap.insert(tmp);

  //complete concentrations plot
  tmp.first = 2;
  tmp.second.name = "Complete Concentrations, Volumes, and Global Quantity Values";
  tmp.second.description = "A plot of all the species concentrations, compartment volumes, and all global quantity values vs. time (includes fixed ones).";
  tmp.second.isPlot = true;
  tmp.second.mTaskType = CCopasiTask::timeCourse;
  mMap.insert(tmp);

  //complete particle numbers plot
  tmp.first = 3;
  tmp.second.name = "Complete Particle Numbers, Volumes, and Global Quantity Values";
  tmp.second.description = "A plot of all the species particle numbers, compartment volumes, and global quantity values vs. time (includes fixed ones).";
  tmp.second.isPlot = true;
  tmp.second.mTaskType = CCopasiTask::timeCourse;
  mMap.insert(tmp);

  //concentration rate plot
  tmp.first = 4;
  tmp.second.name = "Concentration Rates, Volume Rates, and Global Quantity Rates";
  tmp.second.description = "A plot of the rate of change of concentrations of species, compartment volume, and global quantities, which are determined by ODEs or reactions vs. time.";
  tmp.second.isPlot = true;
  tmp.second.mTaskType = CCopasiTask::timeCourse;
  mMap.insert(tmp);

  //particle rate plot
  tmp.first = 5;
  tmp.second.name = "Particle Number Rates, Volume Rates, and Global Quantity Rates";
  tmp.second.description = "A plot of the rate of change of particle numbers of all species, compartment volume, and global quantities, which are determined by ODEs or reactions vs. time.";
  tmp.second.isPlot = true;
  tmp.second.mTaskType = CCopasiTask::timeCourse;
  mMap.insert(tmp);

  //reaction particle flux
  tmp.first = 6;
  tmp.second.name = "Reaction Fluxes";
  tmp.second.description = "A plot of the fluxes of all reactions vs. time, in concentration/time unit.";
  tmp.second.isPlot = true;
  tmp.second.mTaskType = CCopasiTask::timeCourse;
  mMap.insert(tmp);

  //reaction particle flux
  tmp.first = 7;
  tmp.second.name = "Reaction Event Fluxes";
  tmp.second.description = "A plot of the fluxes of all reactions vs. time, in reaction events/time unit.";
  tmp.second.isPlot = true;
  tmp.second.mTaskType = CCopasiTask::timeCourse;
  mMap.insert(tmp);

  //fitting result plots
  tmp.first = 10;
  tmp.second.name = "Parameter Estimation Result";
  tmp.second.description = "Curves of all dependent values of all experiments are created in one plot. For each dependent value the experimental data, the fitted curve, and the weighted error are shown.";
  tmp.second.isPlot = true;
  tmp.second.mTaskType = CCopasiTask::parameterFitting;
  mMap.insert(tmp);

  //fitting result plots
  tmp.first = 11;
  tmp.second.name = "Plots of Parameter Estimation Results per Experiment";
  tmp.second.description = "For each experiment of the parameter estimation a plot is created. Each plot contains the experimental data, the fitted curve, and the weighted error for each dependent value.";
  tmp.second.isPlot = true;
  tmp.second.mTaskType = CCopasiTask::parameterFitting;
  mMap.insert(tmp);

  //fitting result plots
  tmp.first = 12;
  tmp.second.name = "Plots of Parameter Estimation Results per Dependent Value";
  tmp.second.description = "For each dependent value of the parameter estimation a plot is created. Each plot contains the experimental data, the fitted curves, and the weighted errors for each experiment a dependent value occurs.";
  tmp.second.isPlot = true;
  tmp.second.mTaskType = CCopasiTask::parameterFitting;
  mMap.insert(tmp);

  //fitting result plots
  tmp.first = 13;
  tmp.second.name = "Progress of Fit";
  tmp.second.description = "Plot of the sum of squares of residuals vs. number of function evaluations (for parameter estimation).";
  tmp.second.isPlot = true;
  tmp.second.mTaskType = CCopasiTask::parameterFitting;
  mMap.insert(tmp);

  //empty plot
  tmp.first = 99;
  tmp.second.name = "Empty";
  tmp.second.description = "A plot with nothing in it.";
  tmp.second.isPlot = true;
  tmp.second.mTaskType = CCopasiTask::unset;
  mMap.insert(tmp);

  //now the reports
  tmp.first = 999;
  tmp.second.name = "-- Reports";
  tmp.second.description = "";
  tmp.second.isPlot = true;
  tmp.second.mTaskType = CCopasiTask::unset;
  mMap.insert(tmp);

  //concentrations report
  tmp.first = 1000;
  tmp.second.name = "Time, Concentrations, Volumes, and Global Quantity Values";
  tmp.second.description = "A table of time, variable species concentrations, variable compartment volumes, and variable global quantity values.";
  tmp.second.isPlot = false; //report
  tmp.second.mTaskType = CCopasiTask::timeCourse;
  mMap.insert(tmp);

  tmp.first = 1001;
  tmp.second.name = "Time, Particle Numbers, Volumes, and Global Quantity Values";
  tmp.second.description = "A table of time, variable species particle numbers, variable compartment volumes, and variable global quantity values.";
  tmp.second.isPlot = false; //report
  tmp.second.mTaskType = CCopasiTask::timeCourse;
  mMap.insert(tmp);

  tmp.first = 1002;
  tmp.second.name = "Complete Time, Concentrations, Volumes, and Global Quantity Values";
  tmp.second.description = "A table of time, all species concentrations, all compartment volumes, and all global quantity values (includes fixed ones).";
  tmp.second.isPlot = false; //report
  tmp.second.mTaskType = CCopasiTask::timeCourse;
  mMap.insert(tmp);

  tmp.first = 1003;
  tmp.second.name = "Complete Time, Particle Numbers, Volumes, and Global Quantity Values";
  tmp.second.description = "A table of time, all species particle numbers, all compartment volumes, and all global quantity values (includes fixed ones).";
  tmp.second.isPlot = false; //report
  tmp.second.mTaskType = CCopasiTask::timeCourse;
  mMap.insert(tmp);

  tmp.first = 1004;
  tmp.second.name = "Time, Concentration Rates, Volume Rates, and Global Quantity Rates";
  tmp.second.description = "A table of time and the rate of change of concentrations of species, compartment volumes, and global quantities which are determined by reactions or ODEs.";
  tmp.second.isPlot = false; //report
  tmp.second.mTaskType = CCopasiTask::timeCourse;
  mMap.insert(tmp);  //not possible at the moment

  tmp.first = 1005;
  tmp.second.name = "Time, Particle Numbers Rates, Volume Rates, and Global Quantity Rates";
  tmp.second.description = "A table of time and the rate of change of particle numbers of species, compartment volumes, and global quantities which are determined by reactions or ODEs.";
  tmp.second.isPlot = false; //report
  tmp.second.mTaskType = CCopasiTask::timeCourse;
  mMap.insert(tmp);

  tmp.first = 1006;
  tmp.second.name = "Time and Reaction Fluxes";
  tmp.second.description = "A table of the fluxes of all reactions and time, in concentration/time unit.";
  tmp.second.isPlot = false; //report
  tmp.second.mTaskType = CCopasiTask::timeCourse;
  mMap.insert(tmp);

  tmp.first = 1007;
  tmp.second.name = "Time and Reaction Event Fluxes";
  tmp.second.description = "A table of the fluxes of all reactions and time, in reaction events/time unit.";
  tmp.second.isPlot = false; //report
  tmp.second.mTaskType = CCopasiTask::timeCourse;
  mMap.insert(tmp);

  tmp.first = 1008;
  tmp.second.name = "Time and all Variable Values (Concentration Units)";
  tmp.second.description = "This table includes all values which change over a time course. Species are measured in concentration unit and fluxes are in concentration/time unit.";
  tmp.second.isPlot = false; //report
  tmp.second.mTaskType = CCopasiTask::timeCourse;
  mMap.insert(tmp);

  tmp.first = 1009;
  tmp.second.name = "Time and all Variable Values (Particle Number Units)";
  tmp.second.description = "This table includes all values which change over a time course. Species are measured in particle numbers and fluxes are in events/time unit.";
  tmp.second.isPlot = false; //report
  tmp.second.mTaskType = CCopasiTask::timeCourse;
  mMap.insert(tmp);

  tmp.first = 1999;
  tmp.second.name = "Empty";
  tmp.second.description = "A table with nothing in it.";
  tmp.second.isPlot = false; //report
  tmp.second.mTaskType = CCopasiTask::unset;
  mMap.insert(tmp);

  return true;
}

//static
CCopasiObject* COutputAssistant::createDefaultOutput(C_INT32 id, CCopasiTask * task, CCopasiDataModel* pDataModel, bool activate)
{
  if (task == NULL)
    {
      return NULL;
    }

  if (task->getProblem() == NULL)
    {
      return NULL;
    }

  CModel* pModel = task->getProblem()->getModel();

  if (pModel == NULL)
    {
      return NULL;
    }

  bool isReport = (id >= 1000);
  C_INT32 idMod = id % 1000;

  const CCopasiObject* pTime = pModel->getObject(CCopasiObjectName("Reference=Time"));

  std::vector<const CCopasiObject*> data1, tmpdata;
  const CCopasiObject* data2 = NULL;

  switch (idMod)
    {
      case 0:
        data1 =
          CObjectLists::getListOfConstObjects(CObjectLists::NON_CONST_METAB_CONCENTRATIONS, pModel);
        tmpdata =
          CObjectLists::getListOfConstObjects(CObjectLists::NON_CONST_COMPARTMENT_VOLUMES, pModel);
        data1.insert(data1.end(), tmpdata.begin(), tmpdata.end());
        tmpdata =
          CObjectLists::getListOfConstObjects(CObjectLists::NON_CONST_GLOBAL_PARAMETER_VALUES, pModel);
        data1.insert(data1.end(), tmpdata.begin(), tmpdata.end());
        break;
      case 1:
        data1 =
          CObjectLists::getListOfConstObjects(CObjectLists::NON_CONST_METAB_NUMBERS, pModel);
        tmpdata =
          CObjectLists::getListOfConstObjects(CObjectLists::NON_CONST_COMPARTMENT_VOLUMES, pModel);
        data1.insert(data1.end(), tmpdata.begin(), tmpdata.end());
        tmpdata =
          CObjectLists::getListOfConstObjects(CObjectLists::NON_CONST_GLOBAL_PARAMETER_VALUES, pModel);
        data1.insert(data1.end(), tmpdata.begin(), tmpdata.end());
        break;
      case 2:
        data1 =
          CObjectLists::getListOfConstObjects(CObjectLists::METAB_CONCENTRATIONS, pModel);
        tmpdata =
          CObjectLists::getListOfConstObjects(CObjectLists::COMPARTMENT_VOLUMES, pModel);
        data1.insert(data1.end(), tmpdata.begin(), tmpdata.end());
        tmpdata =
          CObjectLists::getListOfConstObjects(CObjectLists::GLOBAL_PARAMETER_VALUES, pModel);
        data1.insert(data1.end(), tmpdata.begin(), tmpdata.end());
        break;
      case 3:
        data1 =
          CObjectLists::getListOfConstObjects(CObjectLists::METAB_NUMBERS, pModel);
        tmpdata =
          CObjectLists::getListOfConstObjects(CObjectLists::COMPARTMENT_VOLUMES, pModel);
        data1.insert(data1.end(), tmpdata.begin(), tmpdata.end());
        tmpdata =
          CObjectLists::getListOfConstObjects(CObjectLists::GLOBAL_PARAMETER_VALUES, pModel);
        data1.insert(data1.end(), tmpdata.begin(), tmpdata.end());
        break;
      case 4:
        data1 =
          CObjectLists::getListOfConstObjects(CObjectLists::METAB_CONC_RATES, pModel);
        tmpdata =
          CObjectLists::getListOfConstObjects(CObjectLists::COMPARTMENT_RATES, pModel);
        data1.insert(data1.end(), tmpdata.begin(), tmpdata.end());
        tmpdata =
          CObjectLists::getListOfConstObjects(CObjectLists::GLOBAL_PARAMETER_RATES, pModel);
        data1.insert(data1.end(), tmpdata.begin(), tmpdata.end());
        break;
      case 5:
        data1 =
          CObjectLists::getListOfConstObjects(CObjectLists::METAB_PART_RATES, pModel);
        tmpdata =
          CObjectLists::getListOfConstObjects(CObjectLists::COMPARTMENT_RATES, pModel);
        data1.insert(data1.end(), tmpdata.begin(), tmpdata.end());
        tmpdata =
          CObjectLists::getListOfConstObjects(CObjectLists::GLOBAL_PARAMETER_RATES, pModel);
        data1.insert(data1.end(), tmpdata.begin(), tmpdata.end());
        break;
      case 6:
        data1 =
          CObjectLists::getListOfConstObjects(CObjectLists::REACTION_CONC_FLUXES, pModel);
        break;
      case 7:
        data1 =
          CObjectLists::getListOfConstObjects(CObjectLists::REACTION_PART_FLUXES, pModel);
        break;
      case 8:
        data1 =
          CObjectLists::getListOfConstObjects(CObjectLists::NON_CONST_METAB_CONCENTRATIONS, pModel);
        tmpdata =
          CObjectLists::getListOfConstObjects(CObjectLists::NON_CONST_COMPARTMENT_VOLUMES, pModel);
        data1.insert(data1.end(), tmpdata.begin(), tmpdata.end());
        tmpdata =
          CObjectLists::getListOfConstObjects(CObjectLists::NON_CONST_GLOBAL_PARAMETER_VALUES, pModel);
        data1.insert(data1.end(), tmpdata.begin(), tmpdata.end());
        tmpdata =
          CObjectLists::getListOfConstObjects(CObjectLists::METAB_CONC_RATES, pModel);
        data1.insert(data1.end(), tmpdata.begin(), tmpdata.end());
        tmpdata =
          CObjectLists::getListOfConstObjects(CObjectLists::COMPARTMENT_RATES, pModel);
        data1.insert(data1.end(), tmpdata.begin(), tmpdata.end());
        tmpdata =
          CObjectLists::getListOfConstObjects(CObjectLists::GLOBAL_PARAMETER_RATES, pModel);
        data1.insert(data1.end(), tmpdata.begin(), tmpdata.end());
        tmpdata =
          CObjectLists::getListOfConstObjects(CObjectLists::REACTION_CONC_FLUXES, pModel);
        data1.insert(data1.end(), tmpdata.begin(), tmpdata.end());
        tmpdata =
          CObjectLists::getListOfConstObjects(CObjectLists::METAB_TRANSITION_TIME, pModel);
        data1.insert(data1.end(), tmpdata.begin(), tmpdata.end());
        break;
      case 9:
        data1 =
          CObjectLists::getListOfConstObjects(CObjectLists::NON_CONST_METAB_NUMBERS, pModel);
        tmpdata =
          CObjectLists::getListOfConstObjects(CObjectLists::NON_CONST_COMPARTMENT_VOLUMES, pModel);
        data1.insert(data1.end(), tmpdata.begin(), tmpdata.end());
        tmpdata =
          CObjectLists::getListOfConstObjects(CObjectLists::NON_CONST_GLOBAL_PARAMETER_VALUES, pModel);
        data1.insert(data1.end(), tmpdata.begin(), tmpdata.end());
        tmpdata =
          CObjectLists::getListOfConstObjects(CObjectLists::METAB_PART_RATES, pModel);
        data1.insert(data1.end(), tmpdata.begin(), tmpdata.end());
        tmpdata =
          CObjectLists::getListOfConstObjects(CObjectLists::COMPARTMENT_RATES, pModel);
        data1.insert(data1.end(), tmpdata.begin(), tmpdata.end());
        tmpdata =
          CObjectLists::getListOfConstObjects(CObjectLists::GLOBAL_PARAMETER_RATES, pModel);
        data1.insert(data1.end(), tmpdata.begin(), tmpdata.end());
        tmpdata =
          CObjectLists::getListOfConstObjects(CObjectLists::REACTION_PART_FLUXES, pModel);
        data1.insert(data1.end(), tmpdata.begin(), tmpdata.end());
        tmpdata =
          CObjectLists::getListOfConstObjects(CObjectLists::METAB_TRANSITION_TIME, pModel);
        data1.insert(data1.end(), tmpdata.begin(), tmpdata.end());
        break;
      case 10:  // :TODO: Implement me!
      {
        CPlotSpecification * pPlotSpecification = NULL;
        CCopasiTask * pTask = (*pDataModel->getTaskList())["Parameter Estimation"];

        if (pTask == NULL) return NULL;

        CFitProblem * pFitProblem = dynamic_cast< CFitProblem * >(pTask->getProblem());

        if (pFitProblem == NULL) return NULL;

        const CExperimentSet & ExperimentSet = pFitProblem->getExperiementSet();
        unsigned C_INT32 i, imax = ExperimentSet.getExperimentCount();

        std::vector< std::string > ChannelX;
        std::vector< std::string > Names;
        std::vector< unsigned C_INT32 > LineTypes;

        for (i = 0; i < imax; i++)
          {
            const CExperiment * pExperiment = ExperimentSet.getExperiment(i);
            const CCopasiVector< CFittingPoint > & FittingPoints = pExperiment->getFittingPoints();

            CCopasiVector< CFittingPoint >::const_iterator it = FittingPoints.begin();
            CCopasiVector< CFittingPoint >::const_iterator end = FittingPoints.end();

            if (it == end) continue;

            data2 = (*it)->getObject(CCopasiObjectName("Reference=Independent Value"));

            for (; it != end; ++it)
              {
                std::string Name = (*it)->getObjectName();
                const CCopasiObject * pObject = pDataModel->getObject(Name);

                if (pObject != NULL)
                  Name = pObject->getObjectDisplayName();

                Name = pExperiment->getObjectName() + "," + Name;

                data1.push_back((*it)->getObject(CCopasiObjectName("Reference=Measured Value")));
                ChannelX.push_back(data2->getCN());
                Names.push_back(Name + "(Measured Value)");
                LineTypes.push_back(2);

                data1.push_back((*it)->getObject(CCopasiObjectName("Reference=Fitted Value")));
                ChannelX.push_back(data2->getCN());
                Names.push_back(Name + "(Fitted Value)");

                if (pExperiment->getExperimentType() == CCopasiTask::timeCourse)
                  LineTypes.push_back(0);
                else
                  LineTypes.push_back(2);

                data1.push_back((*it)->getObject(CCopasiObjectName("Reference=Weighted Error")));
                ChannelX.push_back(data2->getCN());
                Names.push_back(Name + "(Weighted Error)");
                LineTypes.push_back(2);
              }
          }

        pPlotSpecification =
          createPlot(getItemName(id), data2, data1, getItem(id).mTaskType, pDataModel);

        if (pPlotSpecification != NULL)
          {
            CCopasiVector< CPlotItem > & Items = pPlotSpecification->getItems();
            CCopasiVector< CPlotItem >::const_iterator itItem = Items.begin();
            CCopasiVector< CPlotItem >::const_iterator endItem = Items.end();
            std::vector< std::string >::const_iterator itChannelX = ChannelX.begin();
            std::vector< std::string >::const_iterator itName = Names.begin();
            std::vector< unsigned C_INT32 >::const_iterator itLineType = LineTypes.begin();

            while (itItem != endItem)
              {
                (*itItem)->getChannels()[0] = CPlotDataChannelSpec(*itChannelX++);
                (*itItem)->setTitle(*itName++);
                (*itItem)->setActivity(COutputInterface::AFTER);
                (*itItem++)->setValue("Line type", *itLineType++);
              }
          }

        return pPlotSpecification;
      }
      break;
      case 11:
      {
        CPlotSpecification * pPlotSpecification = NULL;
        CCopasiTask * pTask = (*pDataModel->getTaskList())["Parameter Estimation"];

        if (pTask == NULL) return NULL;

        CFitProblem * pFitProblem = dynamic_cast< CFitProblem * >(pTask->getProblem());

        if (pFitProblem == NULL) return NULL;

        const CExperimentSet & ExperimentSet = pFitProblem->getExperiementSet();
        unsigned C_INT32 i, imax = ExperimentSet.getExperimentCount();

        for (i = 0; i < imax; i++)
          {
            const CExperiment * pExperiment = ExperimentSet.getExperiment(i);
            const CCopasiVector< CFittingPoint > & FittingPoints = pExperiment->getFittingPoints();

            CCopasiVector< CFittingPoint >::const_iterator it = FittingPoints.begin();
            CCopasiVector< CFittingPoint >::const_iterator end = FittingPoints.end();

            if (it == end) continue;

            data2 = (*it)->getObject(CCopasiObjectName("Reference=Independent Value"));
            data1.clear();

            for (; it != end; ++it)
              {
                data1.push_back((*it)->getObject(CCopasiObjectName("Reference=Measured Value")));
                data1.push_back((*it)->getObject(CCopasiObjectName("Reference=Fitted Value")));
                data1.push_back((*it)->getObject(CCopasiObjectName("Reference=Weighted Error")));
              }

            pPlotSpecification =
              createPlot(pExperiment->getObjectName(), data2, data1, getItem(id).mTaskType, pDataModel);

            if (pPlotSpecification != NULL)
              {
                CCopasiVector< CPlotItem > & Items = pPlotSpecification->getItems();
                CCopasiVector< CPlotItem >::const_iterator itItem = Items.begin();
                CCopasiVector< CPlotItem >::const_iterator endItem = Items.end();
                it = FittingPoints.begin();

                unsigned C_INT32 LineType;

                if (pExperiment->getExperimentType() == CCopasiTask::timeCourse)
                  LineType = 0;
                else
                  LineType = 2;

                while (itItem != endItem)
                  {
                    std::string Name = (*it++)->getObjectName();
                    const CCopasiObject * pObject = pDataModel->getObject(Name);

                    if (pObject != NULL)
                      Name = pObject->getObjectDisplayName();

                    (*itItem)->setTitle(Name + "(Measured Value)");
                    (*itItem)->setActivity(COutputInterface::AFTER);
                    (*itItem++)->setValue("Line type", (unsigned C_INT32) 2);

                    (*itItem)->setTitle(Name + "(Fitted Value)");
                    (*itItem)->setActivity(COutputInterface::AFTER);
                    (*itItem++)->setValue("Line type", (unsigned C_INT32) LineType);

                    (*itItem)->setTitle(Name + "(Weighted Error)");
                    (*itItem)->setActivity(COutputInterface::AFTER);
                    (*itItem++)->setValue("Line type", (unsigned C_INT32) 2);
                  }
              }
          }

        return pPlotSpecification;
      }
      break;
      case 12:
      {
        CPlotSpecification * pPlotSpecification = NULL;
        CCopasiTask * pTask = (*pDataModel->getTaskList())["Parameter Estimation"];

        if (pTask == NULL) return NULL;

        CFitProblem * pFitProblem = dynamic_cast< CFitProblem * >(pTask->getProblem());

        if (pFitProblem == NULL) return NULL;

        const CExperimentSet & ExperimentSet = pFitProblem->getExperiementSet();
        unsigned C_INT32 i, imax = ExperimentSet.getExperimentCount();

        std::map< const CCopasiObject *, CPlotSpecification * > PlotSpecMap;
        std::map< const CCopasiObject *, CPlotSpecification * >::iterator Found;

        for (i = 0; i < imax; i++)
          {
            const CExperiment * pExperiment = ExperimentSet.getExperiment(i);
            const CCopasiVector< CFittingPoint > & FittingPoints = pExperiment->getFittingPoints();

            CCopasiVector< CFittingPoint >::const_iterator it = FittingPoints.begin();
            CCopasiVector< CFittingPoint >::const_iterator end = FittingPoints.end();

            if (it == end) continue;

            std::string Name = pExperiment->getObjectName();
            CPlotDataChannelSpec ChannelX =
              (*it)->getObject(CCopasiObjectName("Reference=Independent Value"))->getCN();
            unsigned C_INT32 LineType;

            if (pExperiment->getExperimentType() == CCopasiTask::timeCourse)
              LineType = 0;
            else
              LineType = 2;

            for (; it != end; ++it)
              {
                const CCopasiObject * pObject =
                  pDataModel->getObject((*it)->getObjectName());

                if (pObject == NULL) continue;

                if ((Found = PlotSpecMap.find(pObject)) != PlotSpecMap.end())
                  pPlotSpecification = Found->second;
                else
                  {
                    unsigned C_INT32 i = 0;
                    std::ostringstream sname;
                    sname << pObject->getObjectDisplayName();

                    while (!(pPlotSpecification =
                               pDataModel->getPlotDefinitionList()->createPlotSpec(sname.str(),
                                   CPlotItem::plot2d)))
                      {
                        i++;
                        sname.str("");
                        sname << pObject->getObjectDisplayName() << "_" << i;
                      }

                    PlotSpecMap[pObject] = pPlotSpecification;
                  }

                if (pPlotSpecification != NULL)
                  {
                    CPlotItem * pItem =
                      pPlotSpecification->createItem(Name + "(Measured Value)", CPlotItem::curve2d);
                    pItem->setActivity(COutputInterface::AFTER);
                    pItem->setValue("Line type", (unsigned C_INT32) 2);
                    pItem->addChannel(ChannelX);
                    pItem->addChannel((*it)->getObject(CCopasiObjectName("Reference=Measured Value"))->getCN());

                    pItem =
                      pPlotSpecification->createItem(Name + "(Fitted Value)", CPlotItem::curve2d);
                    pItem->setActivity(COutputInterface::AFTER);
                    pItem->setValue("Line type", LineType);
                    pItem->addChannel(ChannelX);
                    pItem->addChannel((*it)->getObject(CCopasiObjectName("Reference=Fitted Value"))->getCN());

                    pItem =
                      pPlotSpecification->createItem(Name + "(Weighted Error)", CPlotItem::curve2d);
                    pItem->setActivity(COutputInterface::AFTER);
                    pItem->setValue("Line type", (unsigned C_INT32) 2);
                    pItem->addChannel(ChannelX);
                    pItem->addChannel((*it)->getObject(CCopasiObjectName("Reference=Weighted Error"))->getCN());
                  }
              }
          }

        return pPlotSpecification;
      }
      break;
      case 13:
      {
        CPlotSpecification * pPlotSpecification = NULL;
        CCopasiTask * pTask = (*pDataModel->getTaskList())["Parameter Estimation"];

        if (pTask == NULL) return NULL;

        CFitProblem * pFitProblem = dynamic_cast< CFitProblem * >(pTask->getProblem());

        if (pFitProblem == NULL) return NULL;

        //        const C_FLOAT64 & SolutionValue = pFitProblem->getSolutionValue();

        data2 = pFitProblem->getObject(CCopasiObjectName("Reference=Function Evaluations"));
        data1.clear();
        data1.push_back(pFitProblem->getObject(CCopasiObjectName("Reference=Best Value")));

        pPlotSpecification =
          createPlot("Progress of Fit" , data2, data1, getItem(id).mTaskType, pDataModel);

        if (pPlotSpecification != NULL)
          {
            pPlotSpecification->setLogY(true);
            CCopasiVector< CPlotItem > & Items = pPlotSpecification->getItems();
            CCopasiVector< CPlotItem >::const_iterator itItem = Items.begin();
            CCopasiVector< CPlotItem >::const_iterator endItem = Items.end();

            while (itItem != endItem)
              {
                (*itItem)->setTitle("sum of squares");
                (*itItem)->setActivity(COutputInterface::DURING);
                (*itItem++)->setValue("Line type", (unsigned C_INT32) 0);
              }
          }

        return pPlotSpecification;
      }
      break;
    }

  if (isReport)
    {
      data1.insert(data1.begin(), pTime);
      CReportDefinition* pReportDef = createTable(getItemName(id), data1, getItem(id).description, getItem(id).mTaskType, pDataModel);

      if (activate && pReportDef)
        {
          task->getReport().setReportDefinition(pReportDef);
          //TODO: also set a default filename?
        }

      return pReportDef;
    }
  else //plot
    {
      data2 = pTime;
      return createPlot(getItemName(id), data2, data1, getItem(id).mTaskType, pDataModel);
    }

  return NULL;
}

#include "plot/COutputDefinitionVector.h"

//static
CPlotSpecification* COutputAssistant::createPlot(const std::string & name,
    const CCopasiObject* x,
    const std::vector<const CCopasiObject*> & y,
    const CCopasiTask::Type & /* taskType */,
    CCopasiDataModel* pDataModel)
{
  if (!x) return NULL;

  std::vector<const CCopasiObject*>::const_iterator it, itEnd = y.end();

  //create plot with unique name
  unsigned C_INT32 i = 0;
  CPlotSpecification* pPl;
  std::ostringstream sname;
  sname << name;
  assert(pDataModel != NULL);

  while (!(pPl = pDataModel->getPlotDefinitionList()->createPlotSpec(sname.str(),
                 CPlotItem::plot2d)))
    {
      i++;
      sname.str("");
      sname << name << "_" << i;
    }

  // Set the task type
  // :TODO: This is currently not implemented for plots.

  //create curves

  CPlotDataChannelSpec name1 = x->getCN();
  CPlotDataChannelSpec name2;
  std::string itemTitle;
  CPlotItem * plItem;

  for (it = y.begin(); it != itEnd; ++it)
    {
      if (!(*it)) continue;

      name2 = (*it)->getCN();
      itemTitle = (*it)->getObjectDisplayName();

      plItem = pPl->createItem(itemTitle, CPlotItem::curve2d);
      plItem->addChannel(name1);
      plItem->addChannel(name2);
    }

  return pPl;
}

//static
CReportDefinition* COutputAssistant::createTable(const std::string & name,
    const std::vector<const CCopasiObject*> & d,
    const std::string & comment,
    const CCopasiTask::Type & taskType,
    CCopasiDataModel* pDataModel)
{
  std::vector<const CCopasiObject*>::const_iterator it, itEnd = d.end();

  //create plot with unique name
  unsigned C_INT32 i = 0;
  CReportDefinition * pReport = NULL;
  std::ostringstream sname;
  sname << name;
  assert(pDataModel != NULL);

  while (!(pReport = pDataModel->getReportDefinitionList()->createReportDefinition(sname.str(), comment)))
    {
      i++;
      sname.str("");
      sname << name << "_" << i;
    }

  // Set the task type
  pReport->setTaskType(taskType);

  pReport->setIsTable(true);
  pReport->setSeparator(CCopasiReportSeparator("\t"));

  for (it = d.begin(); it != itEnd; ++it)
    {
      if (!(*it)) continue;

      pReport->getTableAddr()->push_back((*it)->getCN());
    }

  return pReport;
}
