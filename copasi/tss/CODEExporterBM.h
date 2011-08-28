// Begin CVS Header
//   $Source: /fs/turing/cvs/copasi_dev/copasi/tss/CODEExporterBM.h,v $
//   $Revision: 1.6 $
//   $Name: Build-33 $
//   $Author: gauges $
//   $Date: 2009/02/18 20:55:35 $
// End CVS Header

// Copyright (C) 2008 by Pedro Mendes, Virginia Tech Intellectual
// Properties, Inc., EML Research, gGmbH, University of Heidelberg,
// and The University of Manchester.
// All rights reserved.

// Copyright (C) 2001 - 2007 by Pedro Mendes, Virginia Tech Intellectual
// Properties, Inc. and EML Research, gGmbH.
// All rights reserved.

#ifndef CODEExpoterBM_H__
#define CODEExpoterBM_H__

#include "CODEExporter.h"

class CODEExporterBM: public CODEExporter
  {
  protected:

  public:

    /**
     ** Constructor for the exporter.
     */
    CODEExporterBM();

    virtual bool exportTitleData(const CCopasiDataModel* pDataModel, std::ofstream & outFile);

    virtual void setReservedNames();

    virtual std::string translateTimeVariableName();

    std::string setODEName(const std::string & objName);

    std::string setConcentrationName(const std::string & objName);

    virtual std::string translateObjectName(const std::string & realName);

    std::string testName(const std::string & name);

    virtual bool exportSingleObject(std::ostringstream & which, std::string & name, std::string & expression, std::string & comments);

    virtual bool exportSingleMetabolite(const CMetab* metab, std::string & expression, std::string & comments);

    virtual bool exportSingleCompartment(const CCompartment* comp, std::string & expression, std::string & comments);

    virtual bool exportSingleModVal(const CModelValue* modval, std::string & expression, std::string & comments);

    virtual bool exportSingleModelEntity(const CModelEntity* tmp, std::string & expression, std::string & comments);

    virtual bool exportSingleParameter(const CCopasiParameter* param, std::string & expression, std::string & comments);

    virtual bool exportSingleODE(const CModelEntity* mentity, std::string & equation, std::string & comments);

    bool exportSingleFunction(const CFunction *func, std::set<std::string>& isExported);

    virtual std::string KineticFunction2ODEmember(const CReaction *reac);

    virtual std::string getDisplayFunctionString(CFunction * func);

    virtual std::string getDisplayExpressionString(CExpression * tmp);

    virtual std::string exportTitleString(const unsigned C_INT32 tmp);

    std::map< std::string, unsigned C_INT32 > Frequancy;
    std::set<std::string> NameSet;
  };
#endif
