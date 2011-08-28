// Begin CVS Header
//   $Source: /fs/turing/cvs/copasi_dev/copasi/xml/CCopasiXML.h,v $
//   $Revision: 1.23 $
//   $Name: Build-33 $
//   $Author: gauges $
//   $Date: 2010/03/10 12:51:27 $
// End CVS Header

// Copyright (C) 2010 by Pedro Mendes, Virginia Tech Intellectual
// Properties, Inc., University of Heidelberg, and The University
// of Manchester.
// All rights reserved.

// Copyright (C) 2008 by Pedro Mendes, Virginia Tech Intellectual
// Properties, Inc., EML Research, gGmbH, University of Heidelberg,
// and The University of Manchester.
// All rights reserved.

// Copyright (C) 2001 - 2007 by Pedro Mendes, Virginia Tech Intellectual
// Properties, Inc. and EML Research, gGmbH.
// All rights reserved.

/*!
 \file CCopasiXML.h
 \brief Header file of class CCopasiXML.
 */

/**
 * CCopasiXML class.
 * This class implements a CCopasiXMLInterface to the COPASI XML specified in
 * http://www.copasi.org/schema/copasi.xsd
 *
 * Created for Copasi by Stefan Hoops 2003
 * Copyright Stefan Hoops
 */

#ifndef COPASI_CCopasiXML
#define COPASI_CCopasiXML

#include "xml/CCopasiXMLInterface.h"
#include "utilities/CVersion.h"

class CModel;
class CEvaluationTree;
class CCopasiTask;
class CReportDefinitionVector;
class COutputDefinitionVector;
class CListOfLayouts;
class CLPoint;
class CLDimensions;
class CLBoundingBox;
class CLCurve;

class CCopasiParameter;
class CCopasiParameterGroup;
class CRegisteredObjectName;
class CCopasiDataModel;

#ifdef USE_CRENDER_EXTENSION

class CLLocalRenderInformation;
class CLGlobalRenderInformation;
class CLRenderInformationBase;
class CLLocalStyle;
class CLGlobalStyle;
class CLColorDefinition;
class CLGradientBase;
class CLLinearGradient;
class CLRadialGradient;
class CLLineEnding;
class CLRenderPoint;
class CLRenderCubicBezier;
class CLGroup;
class CLTransformation2D;
class CLImage;
class CLGraphicalPrimitive1D;
class CLText;
class CLRenderCurve;
class CLGraphicalPrimitive2D;
class CLRectangle;
class CLEllipse;
class CLPolygon;
class CLGradientStop;
class CLLineEnding;
class CLStyle;

#endif /* USE_CRENDER_EXTENSION */

class CCopasiXML : public CCopasiXMLInterface
{
  // Operations
public:
  /**
   * Constructor
   */
  CCopasiXML();

  /**
   * Destructor
   */
  ~CCopasiXML();

  /**
   * Save information to a given ostream.
   * @param std::ostream & os
   * @param const std::string & relativeTo
   * @return bool success
   */
  virtual bool save(std::ostream & os,
                    const std::string & relativeTo);

  /**
   * Load information from a given istream.
   * @param std::istream & is
   * @param const std::string & relativeTo
   * @return bool success
   */
  virtual bool load(std::istream & is,
                    const std::string & relativeTo);

  /**
   * Retrieve the version of the current XML file.
   * Before any load operation this contains the COPASI schema version supported by
   * the writer. After load it contains the schema version of the loaded file.
   * @return const CVersion & version
   */
  const CVersion & getVersion() const;

  /**
   * Set the model.
   * @param CModel * pModel
   * @return bool success
   */
  bool setModel(CModel * pModel);

  /**
   * Retreive the model.
   * @return CModel * pModel
   */
  CModel * getModel() const;

  /**
   * Retreive whether the XML contains a model.
   * @return bool have Model
   */
  bool haveModel() const;

  /**
   * Free the model.
   * @return bool success
   */
  bool freeModel();

  /**
   * Set the function list.
   * @param CopasiVectorN< CEvaluationTree > * pFunctionList
   * @return bool success
   */
  bool setFunctionList(CCopasiVectorN< CEvaluationTree > *pFunctionList);

  /**
   * Retreive the function list.
   * @return CCopasiVectorN< CEvaluationTree > * pFunctionList
   */
  CCopasiVectorN< CEvaluationTree > * getFunctionList() const;

  /**
   * Retreive whether the XML contains a function list.
   * @return bool haveFunctionList
   */
  bool haveFunctionList() const;

  /**
   * Free the function list.
   * @return bool success
   */
  bool freeFunctionList();

  /**
   * Set the task list.
   * @param CCopasiVectorN< CCopasiTask > *pTaskList
   * @return bool success
   */
  bool setTaskList(CCopasiVectorN< CCopasiTask > *pTaskList);

  /**
   * Set the datamodel.
   * @param CCopasiDataModel* pDataModel
   * @return bool success
   */
  bool setDatamodel(CCopasiDataModel* pDataModel);

  /**
   * Retreive the task list.
   * @return CCopasiVectorN< CCopasiTask > * taskList
   */
  CCopasiVectorN< CCopasiTask > * getTaskList() const;

  /**
   * Retreive whether the XML contains a task list.
   * @return bool haveTaskList
   */
  bool haveTaskList() const;

  /**
   * Free the task list.
   * @return bool success
   */
  bool freeTaskList();

  /**
   * Set the plot list.
   * @param COutputDefinitionVector * pPlotList
   * @return bool success
   */
  bool setPlotList(COutputDefinitionVector * pPlotList);

  /**
   * Retreive the plot list.
   * @return COutputDefinitionVector * plotList
   */
  COutputDefinitionVector * getPlotList() const;

  /**
   * Retreive whether the XML contains a plot list.
   * @return bool havePlotList
   */
  bool havePlotList() const;

  /**
   * Free the plot list.
   * @return bool success
   */
  bool freePlotList();

  /**
   * Set the report list.
   * @param CReportDefinitionVector *pReportList
   * @return bool success
   */
  bool setReportList(CReportDefinitionVector * pReportList);

  /**
   * Retreive the report list.
   * @return CReportDefinitionVector * reportList
   */
  CReportDefinitionVector * getReportList() const;

  /**
   * Retreive whether the XML contains a report list.
   * @return bool haveReportList
   */
  bool haveReportList() const;

  /**
   * Free the report list.
   * @return bool success
   */
  bool freeReportList();

  /**
   * Set the GUI.
   * @param SCopasiXMLGUI *pGUI
   * @return bool success
   */
  bool setGUI(SCopasiXMLGUI *pGUI);

  /**
   * Retreive the SCopasiXMLGUI.
   * @return SCopasiXMLGUI * pGUI
   */
  SCopasiXMLGUI * getGUI() const;

  /**
   * Retreive whether the XML contains a GUI.
   * @return bool have GUI
   */
  bool haveGUI() const;

  /**
   * Free the GUI.
   * @return bool success
   */
  bool freeGUI();

  /**
   * Set the layout list.
   * @param const CListOfLayouts & reportList
   * @return bool success
   */
  bool setLayoutList(const CListOfLayouts & reportList);

  /**
   * Retreive the layout list.
   * @return CListOfLayouts * layoutList
   */
  CListOfLayouts * getLayoutList() const;

  /**
   * Retreive whether the XML contains a layout list.
   * @return bool haveLayoutList
   */
  bool haveLayoutList() const;

  /**
   * Free the layout list.
   * @return bool success
   */
  bool freeLayoutList();

private:
  /**
   * Save the model.
   * @return bool success
   */
  bool saveModel();

  /**
   * Save the list of functions.
   * @return bool success
   */
  bool saveFunctionList();

  /**
   * Save the list of tasks.
   * @return bool success
   */
  bool saveTaskList();

  /**
   * Save the list of plots.
   * @return bool success
   */
  bool savePlotList();

  /**
   * Save the list of reports.
   * @return bool success
   */
  bool saveReportList();

  /**
   * Save GUI information
   * @return bool success
   */
  bool saveGUI();

  /**
   * Save the list of layout.
   * @return bool success
   */
  bool saveLayoutList();

  void savePosition(const CLPoint& p, const std::string & tag = "Position");

  void saveDimensions(const CLDimensions& d);

  void saveBoundingBox(const CLBoundingBox& bb);

  void saveCurve(const CLCurve& c);

  /**
   * Save the SBML reference information
   * @return bool success
   */
  bool saveSBMLReference();

  /**
   * Save a Report Section such as Header, Body or Footer.
   * @param const std::string & name
   * @param const std::vector <CCopasiObjectName> & section
   * @return bool success
   */
  bool saveReportSection(const std::string & name,
                         const std::vector <CRegisteredObjectName> & section);

  /**
   * Build a list of functions.
   * @return bool success
   */
  bool buildFunctionList();

#ifdef USE_CRENDER_EXTENSION

  /**
   * Saves the list of global render information objects.
   */
  void saveListOfGlobalRenderInformation(const CCopasiVector<CLGlobalRenderInformation>& list);

  /**
   * Saves the list of local render information objects.
   */
  void saveListOfLocalRenderInformation(const CCopasiVector<CLLocalRenderInformation>& list);

  /**
   * Saves a single global render information object.
   */
  void saveGlobalRenderInformation(const CLGlobalRenderInformation& renderInfo);

  /**
   * Saves a single local render information object.
   */
  void saveLocalRenderInformation(const CLLocalRenderInformation& renderInfo);

  /**
   * Saves the attributes that render information objects have in common.
   */
  void saveRenderInformationAttributes(const CLRenderInformationBase& renderInfo, CXMLAttributeList& attributes);

  /**
   * Saves color definitions, gradient definitions  and line endings.
   */
  void saveRenderInformationDefinitionElements(const CLRenderInformationBase& renderInfo);

  /**
   * Save a single color definition element.
   */
  void saveColorDefinition(const CLColorDefinition& color);

  /**
   * Saves a single linear gradient definition.
   */
  void saveLinearGradient(const CLLinearGradient& gradient);

  /**
   * Saves a single radial gradient definition.
   */
  void saveRadialGradient(const CLRadialGradient& gradient);

  /**
   * Adds the attributes common to radial and linear gradient.
   */
  void saveGradientAttributes(const CLGradientBase& gradient, CXMLAttributeList& attributes);

  /**
   * Saves the elements that are common to linear and radial gradients.
   */
  void saveGradientElements(const CLGradientBase& gradient);

  /**
   * Saves a single gradient stop element.
   */
  void saveGradientStop(const CLGradientStop& stop);

  /**
   * Saves a line ending definiton,
   */
  void saveLineEnding(const CLLineEnding& lineEnding);

  /**
   * Saves a single local style element.
   */
  void saveLocalStyle(const CLLocalStyle& style);

  /**
   * Saves a single local style element.
   */
  void saveGlobalStyle(const CLGlobalStyle& style);

  /**
   * Adds the attributes common to both style types.
   */
  void saveStyleAttributes(const CLStyle& style, CXMLAttributeList& attributes);

  /**
   * Saves the elements common to both style types.
   */
  void saveStyleElements(const CLStyle& style);

  /**
   * Saves a group element.
   */
  void saveGroupElement(const CLGroup& group);

  /**
   * Saves the attributes for a transformation.
   */
  void saveTransformationAttributes(const CLTransformation2D& transformation, CXMLAttributeList& attributes);

  /**
   * Saves the attributes for a 1D element
   */
  void save1DAttributes(const CLGraphicalPrimitive1D& primitive, CXMLAttributeList& attributes);

  /**
   * Saves the attributes for a 2D element
   */
  void save2DAttributes(const CLGraphicalPrimitive2D& primitive, CXMLAttributeList& attributes);

  /**
   * Saves the attributes for a text element.
   * We make this a template so that we can use it for a group as well as a text element.
   */
  template<typename TEXTELEMENT>
  void saveTextAttributes(const TEXTELEMENT& text, CXMLAttributeList& attributes);

  /**
   * Saves the startHead and endHead attribute as found in group and curves.
   * We write it as a template so that it can be used on curves and group elements.
   */
  template<typename HEADELEMENT>
  void saveArrowHeadAttributes(const HEADELEMENT& element, CXMLAttributeList& attributes);

  /**
   * Saves a class that is subclasses from Transformation2D.
   * This covers images, curves, rectangles, ellipses, polygons, text elements and groups.
   */
  void saveTransformation2D(const CLTransformation2D& transformation);

  /**
   * saves a single image element.
   */
  void saveImageElement(const CLImage& image);

  /**
   * saves a single rectangle element.
   */
  void saveRectangleElement(const CLRectangle& rectangle);

  /**
   * saves a single ellipse element.
   */
  void saveEllipseElement(const CLEllipse& ellipse);

  /**
   * saves a single text element.
   */
  void saveRenderTextElement(const CLText& text);

  /**
   * saves a single image element.
   */
  void savePolygonElement(const CLPolygon& polygon);

  /**
   * saves a single image element.
   */
  void saveRenderCurveElement(const CLRenderCurve& curve);

  /**
   * saves a vector of curve elements. This can be called from the polygon as well as the curve.
   */
  void saveCurveElements(const std::vector<CLRenderPoint*>& curveElements);

  /**
   * saves a single render point element.
   */
  void saveRenderPoint(const CLRenderPoint& point);

#endif /* USE_CRENDER_EXTENSION */

  // Attributes

  /**
   * The version of the COPASI XML Schema the current file adheres to.
   */
  CVersion mVersion;

  /**
   * Pointer to a model which has been loaded or is to be saved.
   * The ownership is handed to the user.
   */
  CModel * mpModel;

  /**
   * Pointer to a vector of functions which has been loaded or is to be saved.
   * The ownership is handed to the user.
   */
  CCopasiVectorN< CEvaluationTree > * mpFunctionList;

  /**
   * Pointer to a vector of tasks which has been loaded or is to be saved.
   * The ownership is handed to the user.
   */
  CCopasiVectorN< CCopasiTask > * mpTaskList;

  /**
   * Pointer to a vector of reports which has been loaded or is to be saved.
   * The ownership is handed to the user.
   */
  CReportDefinitionVector * mpReportList;

  /**
   * Pointer to a vector of plots which has been loaded or is to be saved.
   * The ownership is handed to the user.
   */
  COutputDefinitionVector * mpPlotList;

  /**
   * Pointer to a GUI related information, which has been loaded or is to be saved.
   * The ownership is handed to the user.
   */
  SCopasiXMLGUI * mpGUI;

  /**
   * Pointer to a vector of plots which has been loaded or is to be saved.
   * The ownership is handed to the user.
   */
  CListOfLayouts * mpLayoutList;

  /**
   * SBML Reference
   */
  std::map< std::string, std::string > mSBMLReference;

  /**
   * Pointer to the datamodel
   */
  CCopasiDataModel* mpDataModel;
};

#endif // COPASI_CCopasiXML
