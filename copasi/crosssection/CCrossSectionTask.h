// Begin CVS Header
//   $Source: /fs/turing/cvs/copasi_dev/copasi/crosssection/CCrossSectionTask.h,v $
//   $Revision: 1.1 $
//   $Name: Build-33 $
//   $Author: ssahle $
//   $Date: 2010/05/14 22:20:55 $
// End CVS Header

// Copyright (C) 2010 by Pedro Mendes, Virginia Tech Intellectual
// Properties, Inc., University of Heidelberg, and The University
// of Manchester.
// All rights reserved.

#ifndef CCROSSSECTIONTASK_H
#define CCROSSSECTIONTASK_H

#include <iostream>

#include "utilities/CCopasiTask.h"

class CCrossSectionTask : public CCopasiTask
{
public:

  /**
   * Default constructor
   * @param const CCopasiContainer * pParent (default: NULL)
   */
  CCrossSectionTask(const CCopasiContainer * pParent = NULL);

  /**
   * Copy constructor
   * @param const CCrossSectionTask & src
   * @param const CCopasiContainer * pParent (default: NULL)
   */
  CCrossSectionTask(const CCrossSectionTask & src,
                    const CCopasiContainer * pParent = NULL);

  /**
   * Destructor
   */
  virtual ~CCrossSectionTask();

  /**
   * Resizes result matrices and updates array annotations.
   * This is used when we need to know about the data structures of a task result
   * without actually performing the task, e.g. when selecting objects for output.
   * For now we assume that this functionality is also performed when
   * initialize() is called.
   */
  //virtual bool updateMatrices();

  /**
   * Initialize the task. If an ostream is given this ostream is used
   * instead of the target specified in the report. This allows nested
   * tasks to share the same output device.
   * @param const OutputFlag & of
   * @param COutputHandler * pOutputHandler
   * @param std::ostream * pOstream (default: NULL)
   * @return bool success
   */
  virtual bool initialize(const OutputFlag & of,
                          COutputHandler * pOutputHandler,
                          std::ostream * pOstream);

  /**
   * Process the task with or without initializing to the initial state.
   * @param const bool & useInitialValues
   * @return bool success
   */
  virtual bool process(const bool & useInitialValues);

  /**
   * This is the output method for any object. The default implementation
   * provided with CCopasiObject uses the ostream operator<< of the object
   * to print the object.To overide this default behaviour one needs to
   * reimplement the virtual print function.
   * @param std::ostream * ostream
   */
  virtual void print(std::ostream * ostream) const;

  // Friend functions
  friend std::ostream &operator<<(std::ostream &os,
                                  const CCrossSectionTask &A);

private:
  /**
   * cleanup()
   */
  void cleanup();
};

#endif // CCROSSSECTIONTASK_H
