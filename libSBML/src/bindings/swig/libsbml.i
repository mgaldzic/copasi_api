/**
 * \file    libsbml.i
 * \brief   Language-independent SWIG directives for wrapping libSBML
 * \author  Ben Bornstein and Ben Kovitz
 *
 * $Id: libsbml.i 11666 2010-08-05 05:12:53Z mhucka $
 * $HeadURL: https://sbml.svn.sourceforge.net/svnroot/sbml/trunk/libsbml/src/bindings/swig/libsbml.i $
 */
/* Copyright 2004 California Institute of Technology and Japan Science and
 * Technology Corporation.
 *
 * This library is free software; you can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation.  A copy of the license agreement is
 * provided in the file named "LICENSE.txt" included with this software
 * distribution.  It is also available online at
 * http://sbml.org/software/libsbml/license.html
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this library; if not, write to the Free Software Foundation,
 * Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA.
 */


%module libsbml

%pragma(java) moduleclassmodifiers="
/**
  * Wrapper class for global methods and constants defined by libSBML.
  * <p>
  * <em style='color: #555'>
  * This class of objects is defined by libSBML only and has no direct
  * equivalent in terms of SBML components.  This class is not prescribed by
  * the SBML specifications, although it is used to implement features
  * defined in SBML.
  * </em>
  * <p>
  * In the C++ and C versions of libSBML, there exists a small number of
  * methods that are global in scope; in addition, libSBML uses a number
  * of enum's to define such things as error codes in a way that can be
  * used by both C++ and C.  This poses a problem in languages such as
  * Java, where there is no concept of global method or global constant.
  * SWIG wraps these global identifiers in the class whose documentation
  * you see before you.
  */
public class"


%{
#include "libsbml.h"

LIBSBML_CPP_NAMESPACE_USE

#ifdef USE_LAYOUT
#include "../swig/layout.h"
#endif /* USE_LAYOUT */
#include "local.cpp"
%}

%import  sbml/common/libsbml-namespace.h
%import  sbml/common/extern.h
%import  sbml/common/sbmlfwd.h
%import  sbml/xml/XMLExtern.h

/**
 * Wraps List class by ListWrapper<TYPENAME> template class.
 * TYPENAME is replaced with a corresponding data type which is
 * stored in the List object (e.g. ModelCreator, CVTerm and ASTNode). 
 *
 * ListWrapper<TYPENAME> class is wrapped as TYPENAMEList class
 * (e.g. ListWrapper<CVTerm> -> CVTermList)
 *
 */

%include "ListWrapper.h"
%template(ModelCreatorList) ListWrapper<ModelCreator>;
%template(DateList)         ListWrapper<Date>;
%template(CVTermList)       ListWrapper<CVTerm>;
%template(ASTNodeList)      ListWrapper<ASTNode>;

/**
 *
 * Includes a language specific interface file.
 *
 */

%include local.i

/**
 * Unfortunately, SWIG makes no distinction between const and non-const
 * member functions (SWIG 1.3 Manual, Section 6.25), but in libSBML C++ we
 * have both const and non-const versions of most getter methods.  To avoid
 * a ton of warning messages about 'const' methods not being wrapped, we
 * disable Warning(516).
 */
#pragma SWIG nowarn=516


/**
 * SWIG is producing Warning(401): Nothing known about base class 
 * 'std::invalid_argument'. Ignored.
 * 
 * The SWIG documentation states:
 *   Note: For the best results, SWIG requires all base classes to be defined 
 *   in an interface. Otherwise, you may get a warning message like this.
 *
 * Suggested solutions are:
 *   1. Using %import to include the file that defines the unknown class
 *   2. Defining the unknown class as an empty class
 *   3. Ignoring the warning
 *
 * I (SK) tried the following:
 * 
 * 1) an import directive:   %import <exception.i>
 * 
 * this does not avoid the warning. It also breaks the exception handling - 
 * certainly in the csharp bindings. libSBML crashes with unhandled exceptions.
 * 
 * 2) Putting an empty class
 * 
 *    %inline
 *    %{
 *    class std::invalid_argument()
 *    {
 *    };
 *    %}
 * 
 * this does not avoid the warning either and produces a wrap file that 
 * will not build.
 *
 * Since this is the only class that produces the warning and the SWIG documentation
 * says "If any base class is undefined, SWIG still generates correct type 
 * relationships." we decided it was safe to ignore the warning for now and 
 * revisit the issue for future releases.
 */
#pragma SWIG nowarn=401


/**
 * Ignore the Visitor pattern accept() method (for now) on all SBML
 * objects.
 */
%ignore *::accept;

/**
 * Ignore internal implementation methods in ASTNode.h
 */
%ignore ASTNode(Token_t*);
%ignore ASTNode::setSemanticsFlag;
%ignore ASTNode::unsetSemanticsFlag;
%ignore ASTNode::getSemanticsFlag;
%ignore ASTNode::setDefinitionURL;
%ignore ASTNode::setParentSBMLObject;

/**
 * SWIG makes no distinction between int and long arguments.
 * (SWIG 1.3 Manual, Section 6.15.2)
 */
%ignore ASTNode::setValue(int);

/**
 * Ignore operator= and operator<< on all SBML objects.
 */
%ignore *::operator=;
%ignore *::operator<<;
%ignore operator==;
%ignore operator!=;

/**
 * Ignore certain internal implementation methods on all objects.
 */
%ignore *::writeElements;
%ignore *::getElementPosition;
%ignore *::setSBMLDocument;
%ignore *::setParentSBMLObject;
%ignore *::setInternalId;
%ignore *::getInternalId;

/**
 * Ignore internal implementation methods in MathML.h
 */
%ignore readMathML;
%ignore writeMathML;

/**
 * Ignore methods whose pointer argument serves as both input and output
 */
%ignore XMLAttributes::readInto;

/**
 * Ignore methods which receive List* or return List*.
 */
%ignore RDFAnnotationParser::parseRDFAnnotation(const XMLNode * annotation, List * CVTerms);
%ignore ASTNode::fillListOfNodes;
%ignore ASTNode::getListOfNodes(ASTNodePredicate predicate) const;

/**
 * Ignore methods which receive or return void*.
 */
%ignore ASTNode::setUserData;
%ignore ASTNode::getUserData;

/**
 * Ignore methods which receive std::list.
 */
%ignore XMLErrorLog::add(const std::list<XMLError>& errors);
%ignore SBMLErrorLog::add(const std::list<SBMLError>& errors);

/**
 * Ignore 'static ParentMap mParent;' in SBO.h
 */
%ignore mParent;

/**
 * Ignore 'struct xmlErrorTableEntry' in XMLError.h.
 */
%ignore xmlErrorTableEntry;

/**
 * Both "const std::string& SBase::getMetaId() const" and
 * "std:string& SBase::getMetaId()" are defined in SBase.cpp.
 * By default, SWIG doesn't convert non-const std:string& to and from
 * target language string.
 * So we ignore the non-const version.
 */
%ignore SBase::getMetaId();

/**
 * Ignore internal FormulaUnitsData methods on SBase
 */
%ignore SBase::removeDuplicateAnnotations;
%ignore SBase::setSBMLNamespaces;
%ignore SBase::getSBMLNamespaces;
%ignore SBase::read;
%ignore SBase::write;

/**
 * Ignore internal FormulaUnitsData methods on Model
 */
%ignore Model::addFormulaUnitsData;
%ignore Model::createFormulaUnitsData;
%ignore Model::getFormulaUnitsData;
%ignore Model::getListFormulaUnitsData;
%ignore Model::getNumFormulaUnitsData;
%ignore Model::isBoolean;
%ignore Model::removeMetaId;
%ignore Model::removeSBOTerms;
%ignore Model::removeHasOnlySubstanceUnits;
%ignore Model::removeSBOTermsNotInL2V2;
%ignore Model::removeDuplicateTopLevelAnnotations;
%ignore Model::convertToL1;
%ignore Model::convertToL2;
%ignore Model::convertToL2V1;
%ignore Model::convertToL2V2;
%ignore Model::convertToL2Strict;

/**
 * Ignore internal implementation methods in Rule
 */
%ignore Rule::setInternalIdOnly;
%ignore Rule::getInternalIdOnly;

/**
 * Ignore internal implementation methods in SpeciesReference
 */
%ignore SpeciesReference::sortMath;

/**
 * Ignore internal implementation methods in UnitDefinition
 */
%ignore UnitDefinition::areIdenticalSIUnits;

/**
 * Ignore internal implementation methods in XMLAttributes
 */
%ignore XMLAttributes::addResource;
%ignore XMLAttributes::write;
%ignore XMLAttributes::setErrorLog;

/**
 * Ignore internal implementation methods in Event
 */
%ignore Event::setInternalIdOnly;

/**
 * Ignore internal implementation methods in SBO
 */
%ignore SBO::readTerm;
%ignore SBO::writeTerm;


/**
 * Ignore internal implementation methods in SBMLErrorLog
 */
%ignore SBMLErrorLog::logError;
%ignore SBMLErrorLog::add;
%ignore SBMLErrorLog::remove;
%ignore SBMLErrorLog::SBMLErrorLog;

/**
 * Ignore internal implementation methods in XMLErrorLog
 */
%ignore XMLErrorLog::XMLErrorLog;
%ignore XMLErrorLog::add;
%ignore XMLErrorLog::setParser;


/**
 * Ignore internal implementation methods in ModelCreator
 */
%ignore ModelCreator::getAdditionalRDF;

/**
 * Ignore internal implementation methods in RDFAnnotationParser
 */
%ignore RDFAnnotationParser::hasRDFAnnotation;
%ignore RDFAnnotationParser::hasAdditionalRDFAnnotation;
%ignore RDFAnnotationParser::hasCVTermRDFAnnotation;
%ignore RDFAnnotationParser::hasHistoryRDFAnnotation;

/**
 * Ignore internal implementation methods in SyntaxChecer
 */
%ignore SyntaxChecker::isAllowedElement;
%ignore SyntaxChecker::hasDeclaredNS;
%ignore SyntaxChecker::isCorrectHTMLNode;

/**
 * Ignore internal implementation methods in SBMLNamespces
 */
%ignore SBMLNamespaces::setLevel;
%ignore SBMLNamespaces::setVersion;
%ignore SBMLNamespaces::setNamespaces;

/**
 * Ignore internal SBMLTransforms.
 */
%ignore SBMLTransforms::replaceFD;
%ignore SBMLTransforms::expandInitialAssignments;
%ignore SBMLTransforms::evaluateASTNode;
%ignore SBMLTransforms::mapComponentValues;

/**
 * Ignore internal implementation methods in XMLToken
 */
%ignore XMLToken::write;

/**
 * Ignore internal implementation methods in XMLNode
 */
%ignore XMLNode::XMLNode(XMLInputStream&);
%ignore XMLNode::write;

/**
 * Ignore internal implementation methods in XMLOutputStream
 */
%ignore XMLOutputStream::getStringStream;

/**
 * Ignore internal implementation classes
 */
%ignore XMLOutputStringStream;
%ignore XMLOutputFileStream;

/**
 * Ignore the unsigned int version of XMLOutputStream::writeAttribute method
 * in order to properly wrap the long version of XMLOutputStream::writeAttribute 
 * method which should be used instead of the unsigned int version.
 */
%ignore XMLOutputStream::writeAttribute(const std::string&, const unsigned int&);
%ignore XMLOutputStream::writeAttribute(const XMLTriple&,   const unsigned int&);

/**
 * The following methods will create new objects.  To prevent memory
 * leaks we must inform SWIG of this.
 */

%typemap(newfree) char * "free($1);";

%newobject *::clone;
%newobject SBase::toSBML;
%newobject SBMLReader::readSBMLFromString;
%newobject SBMLReader::readSBMLFromFile;
%newobject SBMLReader::readSBML;
%newobject readSBML(const char *);
%newobject readSBMLFromString(const char *);
%newobject readSBMLFromFile(const char *);
%newobject SBMLWriter::writeToString;
%newobject writeSBMLToString;
%newobject readMathMLFromString;
%newobject writeMathMLToString;
%newobject SBML_formulaToString;
%newobject SBML_parseFormula;
%newobject ASTNode::deepCopy;
%newobject ASTNode::getListOfNodes();
%newobject *::remove;
%newobject Model::removeFunctionDefinition;
%newobject Model::removeUnitDefinition;
%newobject Model::removeCompartmentType;
%newobject Model::removeSpeciesType;
%newobject Model::removeSpecies;
%newobject Model::removeCompartment;
%newobject Model::removeParameter;
%newobject Model::removeInitialAssignment;
%newobject Model::removeRule;
%newobject Model::removeConstraint;
%newobject Model::removeReaction;
%newobject Model::removeEvent;
%newobject Reaction::removeReactant;
%newobject Reaction::removeProduct;
%newobject Reaction::removeModifier;
%newobject Event::removeEventAssignment;
%newobject UnitDefinition::removeUnit;
%newobject KineticLaw::removeParameter;
%newobject KineticLaw::removeLocalParameter;
%newobject RDFAnnotationParser::parseRDFAnnotation(XMLNode *);
%newobject RDFAnnotationParser::deleteRDFAnnotation;
%newobject RDFAnnotationParser::parseCVTerms;
%newobject RDFAnnotationParser::parseModelHistory;
%newobject RDFAnnotationParser::createRDFAnnotation;
%newobject RDFAnnotationParser::createAnnotation;
%newobject RDFAnnotationParser::createRDFDescription;
%newobject RDFAnnotationParser::createCVTerms;
%newobject XMLNode::removeChild;
%newobject XMLNode::convertStringToXMLNode;
%newobject Unit::convertToSI;
%newobject UnitDefinition::convertToSI;
%newobject UnitDefinition::combine;

#ifdef USE_LAYOUT
%newobject Model::removeLayout;
#endif /* USE_LAYOUT */

/**
 * In the wrapped languages, these methods will appear as:
 *
 *  - libsbml.formulaToString()
 *  - libsbml.parseFormula()
 */
%rename(formulaToString) SBML_formulaToString;
%rename(parseFormula)    SBML_parseFormula;

/**
 * 
 * wraps "List* ASTNode::getListOfNodes(ASTNodePredicate)" function
 * as "ListWrapper<ASTNode>* ASTNode::getListOfNodes()" function 
 * which returns a list of all ASTNodes. 
 *
 */

%inline
%{
  int ASTNode_true(const ASTNode *node)
  {
    return 1;
  }
%}

%extend ASTNode
{
  ListWrapper<ASTNode>* getListOfNodes()
  {
    List *list = $self->getListOfNodes(ASTNode_true);
    return new ListWrapper<ASTNode>(list);
  }
}

/*
 * Wraps "static void RDFAnnotationParser::parseRDFAnnotation(const XMLNode *annotation, 
 * List *CVTerms)" function as 
 * "static void RDFAnnotationParser::parseRDFAnnotation(const XMLNode *annotation, 
 *  ListWrapper<CVTerm> *CVTerms);
 *
 */

%extend RDFAnnotationParser
{
  static void RDFAnnotationParser::parseRDFAnnotation(const XMLNode *annotation, 
                                                      ListWrapper<CVTerm> *CVTerms)
  {
    if (!CVTerms) return;

    List *list = CVTerms->getList();
    RDFAnnotationParser::parseRDFAnnotation(annotation,list);
  }
}



/**
 * Wrap these files.
 */

%include "std_string.i"

%include sbml/common/libsbml-version.h
%include sbml/common/operationReturnValues.h

%include sbml/SBMLReader.h
%include sbml/SBMLWriter.h
%include sbml/SBMLTypeCodes.h
%include sbml/SBase.h
%include sbml/ListOf.h
%include sbml/Model.h
%include sbml/SBMLDocument.h
%include sbml/FunctionDefinition.h
%include sbml/UnitKind.h
%include sbml/Unit.h
%include sbml/UnitDefinition.h
%include sbml/CompartmentType.h
%include sbml/SpeciesType.h
%include sbml/Compartment.h
%include sbml/Species.h
%include sbml/Parameter.h
%include sbml/LocalParameter.h
%include sbml/InitialAssignment.h
%include sbml/Rule.h
%include sbml/Constraint.h
%include sbml/Reaction.h
%include sbml/KineticLaw.h
%include sbml/SpeciesReference.h
%include sbml/Event.h
%include sbml/EventAssignment.h
%include sbml/Trigger.h
%include sbml/Delay.h
%include sbml/SBO.h
%include sbml/SyntaxChecker.h
%include sbml/StoichiometryMath.h
%include sbml/SBMLNamespaces.h
%include sbml/SBMLTransforms.h

%include sbml/math/MathML.h
%include sbml/math/ASTNode.h
%include sbml/math/FormulaParser.h
%include sbml/math/FormulaFormatter.h

%include sbml/xml/XMLAttributes.h
%include sbml/xml/XMLNamespaces.h
%include sbml/xml/XMLToken.h
%include sbml/xml/XMLNode.h
%include sbml/xml/XMLTriple.h
%include sbml/xml/XMLOutputStream.h
%include sbml/xml/XMLError.h
%include sbml/xml/XMLErrorLog.h

%include sbml/SBMLErrorLog.h
%include sbml/SBMLError.h

%include sbml/annotation/CVTerm.h
%include sbml/annotation/ModelHistory.h
%include sbml/annotation/RDFAnnotation.h

#ifdef USE_LAYOUT
%include ../swig/layout.i
#endif /* USE_LAYOUT */
