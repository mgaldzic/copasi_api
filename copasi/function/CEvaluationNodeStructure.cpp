/* Begin CVS Header
$Source: /fs/turing/cvs/copasi_dev/copasi/function/CEvaluationNodeStructure.cpp,v $
$Revision: 1.9 $
$Name: Build-33 $
$Author: gauges $
$Date: 2009/02/19 15:37:57 $
End CVS Header */

// Copyright (C) 2008 by Pedro Mendes, Virginia Tech Intellectual
// Properties, Inc., EML Research, gGmbH, University of Heidelberg,
// and The University of Manchester.
// All rights reserved.

// Copyright (C) 2001 - 2007 by Pedro Mendes, Virginia Tech Intellectual
// Properties, Inc. and EML Research, gGmbH.
// All rights reserved.

#include <limits.h>

#include "copasi.h"
#include "CEvaluationNode.h"

#include "sbml/math/ASTNode.h"

CEvaluationNodeStructure::CEvaluationNodeStructure():
    CEvaluationNode(CEvaluationNode::INVALID, "")
{}

CEvaluationNodeStructure::CEvaluationNodeStructure(const SubType & subType,
    const Data & data):
    CEvaluationNode((Type) (CEvaluationNode::STRUCTURE | subType), data)
{
  switch (subType)
    {
    case OPEN:
    case VECTOR_OPEN:
      mPrecedence = PRECEDENCE_STRUCTURE_OPEN;
      break;

    case COMMA:
      mPrecedence = PRECEDENCE_STRUCTURE_COMMA;
      break;

    case CLOSE:
    case VECTOR_CLOSE:
      mPrecedence = PRECEDENCE_STRUCTURE_CLOSE;
      break;

    case INVALID:
      fatalError();
      break;
    }
}

CEvaluationNodeStructure::CEvaluationNodeStructure(const CEvaluationNodeStructure & src):
    CEvaluationNode(src)
{}

CEvaluationNodeStructure::~CEvaluationNodeStructure() {}

ASTNode* CEvaluationNodeStructure::toAST(const CCopasiDataModel* /*pDataModel*/) const
  {
    fatalError();
    return NULL;
  }
