/**
 * @file    local.i
 * @brief   Python-specific SWIG directives for wrapping libSBML API
 * @author  Ben Bornstein
 * @author  Ben Kovitz
 * @author  Akiya Jouraku
 *
 * $Id: local.i 10866 2010-01-29 19:52:27Z mhucka $
 * $HeadURL: https://sbml.svn.sourceforge.net/svnroot/sbml/trunk/libsbml/src/bindings/python/local.i $
 *
 *<!---------------------------------------------------------------------------
 * This file is part of libSBML.  Please visit http://sbml.org for more
 * information about SBML, and the latest version of libSBML.
 *
 * Copyright 2005-2010 California Institute of Technology.
 * Copyright 2002-2005 California Institute of Technology and
 *                     Japan Science and Technology Corporation.
 * 
 * This library is free software; you can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation.  A copy of the license agreement is provided
 * in the file named "LICENSE.txt" included with this software distribution
 * and also available online as http://sbml.org/software/libsbml/license.html
 *----------------------------------------------------------------------- -->*/


/**
 * Turn on (minimal) Python docstrings and then append our own.
 */
%feature("autodoc", "1");
%include "pydoc.i"

/**
 *  Wraps std::cout, std::cerr, std::clog, std::ostream, and std::ostringstream, 
 *
 * (sample code) -----------------------------------------------------
 *
 * 1. wraps std::cout
 *
 *    xos = libsbml.XMLOutputStream(libsbml.cout)
 *
 * 2. wraps std::cerr
 *
 *    d = libsbml.readSBML("foo.xml")
 *    if d.getNumErrors() > 0 :
 *       d.printErrors(libsbml.cerr)
 *    
 *
 * 3. wraps std::ostringstream
 *
 *    oss = libsbml.ostringstream()
 *    xos = libsbml.XMLOutputStream(oss)
 *    ...
 *    libsbml.endl(oss)
 *    s = oss.str()
 * 
 */

// ignores C++ specific methods in std::string.
%ignore std::basic_string<char>::begin;
%ignore std::basic_string<char>::end;
%ignore std::basic_string<char>::rbegin;
%ignore std::basic_string<char>::rend;
%ignore std::basic_string<char>::get_allocator;
%ignore std::basic_string<char>::capacity;
%ignore std::basic_string<char>::reserve;

%include <std_alloc.i>
%include <std_basic_string.i>
%include <std_string.i>

#pragma SWIG nowarn=509
%warnfilter(401) basic_ios<char>;

namespace std
{
  // Template class basic ios
  template<typename _CharT, typename _Traits = char_traits<_CharT> >
  class basic_ios : public ios_base {};

  // Template class basic_ostream
  template<typename _CharT, typename _Traits = char_traits<_CharT> >
  class basic_ostream : virtual public basic_ios<_CharT, _Traits> 
  {
    public:
      explicit
      basic_ostream(std::basic_streambuf<_CharT, _Traits>* __sb);
      virtual 
      ~basic_ostream();
  };

  // Template class basic_ostringstream
  template<typename _CharT, typename _Traits = char_traits<_CharT>,
           typename _Alloc = allocator<_CharT> >
  class basic_ostringstream : public basic_ostream<_CharT, _Traits>
  {
    public:
      explicit
      basic_ostringstream(std::ios_base::openmode __mode = std::ios_base::out);
      ~basic_ostringstream();

      basic_string<_CharT, _Traits, _Alloc> 
      str() const;

      void
      str(const basic_string<_CharT, _Traits, _Alloc>& __s);
  };

  template<typename _CharT, typename _Traits = char_traits<_CharT> >
  basic_ostream<_CharT, _Traits>& 
  endl(basic_ostream<_CharT, _Traits>&);

  template<typename _CharT, typename _Traits = char_traits<_CharT> >
  basic_ostream<_CharT, _Traits>& 
  flush(basic_ostream<_CharT, _Traits>&);
}

namespace std
{
  /**
   *  std::ostream and std::ostringstream 
   *  (std::ios is not wrapped)
   */
  typedef basic_ios<char>           ios;
  typedef basic_ostream<char>       ostream ;
  typedef basic_ostringstream<char> ostringstream ;

  %template()              basic_ios<char>;
  %template(ostream)       basic_ostream<char>;
  %template(ostringstream) basic_ostringstream<char>;

  /**
   *  output manipulators
   */
  %template(endl)  endl<char, char_traits<char> >;
  %template(flush) flush<char, char_traits<char> >;

  /**
   *  std::cout, std::cerr, and std::clog.
   */
  %immutable;
  extern std::ostream cout;
  extern std::ostream cerr;
  extern std::ostream clog;
  %mutable;
}


/**
 * Convert an SBase object to a string.
 *
%extend SBase
{
  %pythoncode
  {
    def __str__(self):
      return self.toSBML()
  }
}*/


/**
 * Allows ListOf objects:
 *
 *   - To be indexed and sliced, e.g. lst[0].
 */
%extend ListOf
{
  int __len__()
  {
    return self->size();
  }

  %pythoncode
  {
    def __getitem__(self, key):

      try:
         keyIsSlice = isinstance(key, slice)
      except:
         keyIsSlice = 0

      if keyIsSlice:
        start = key.start
        if start is None:
          start = 0
        stop = key.stop
        if stop is None:
          stop = self.size()
        return [self[i] for i in range(
          self._fixNegativeIndex(start), self._fixNegativeIndex(stop)
        )]

      key = self._fixNegativeIndex(key)
      if key < 0 or key >= self.size():
        raise IndexError(key)
      return self.get(key)


    def _fixNegativeIndex(self, index):
      if index < 0:
        return index + self.size()
      else:
        return index


    def __iter__(self):
      for i in range(self.size()):
        yield self[i]


    def __repr__(self):
      return "[" + ", ".join([repr(self[i]) for i in range(len(self))]) + "]"


    def __str__(self):
      return repr(self)
  }
}



/**
 * Convert SBase, SimpleSpeciesReference and Rule objects into the most specific type possible.
 */
%typemap(out) SBase*, SimpleSpeciesReference*, Rule*
{
  $result = SWIG_NewPointerObj(SWIG_as_voidptr($1), GetDowncastSwigType($1),
                               $owner | %newpointer_flags);
}


/*
 * SWIG-generated wrapper code wrongly invokes 
 * XMLOutputStream::writeAttribute(.., const unsigned int& value) instead of
 * XMLOutputStream::writeAttribute(.., const bool& value) even if the writeAttribute 
 * function invoked with a bool value (True or False) in Python code.
 * It seems that a bool value could be casted to unsigned int, int, or long value 
 * in SWIG-generated internal type check code when these types are overloaded in the
 * wrapped function.
 *
 * To avoid this problem, XMLOutputStream::writeAttribute(.., const bool& value)
 * is internally wrapped as XMLOutputStream::writeAttributeBool(.., const bool&) 
 * and this function is properly invoked when the writeAttribute function is invoked
 * with a bool value in Python code.
 */

%extend XMLOutputStream
{
  void writeAttributeBool(const std::string& name, const bool& value)
  {
    $self->writeAttribute(name, value);
  }

  void writeAttributeBool(const XMLTriple& name, const bool& value)
  {
    $self->writeAttribute(name, value);
  }
}

#if SWIG_VERSION > 0x010336
%feature("pythonprepend")
XMLOutputStream::writeAttribute
%{
        if type(args[1]) == type(True): return _libsbml.XMLOutputStream_writeAttributeBool(self, *args)
%}
#else 
%feature("pythonprepend")
XMLOutputStream::writeAttribute
%{
        if type(args[2]) == type(True): return _libsbml.XMLOutputStream_writeAttributeBool(*args)
%}
#endif

/**
 * Add an equality operator to SBase.  All subclasses of SBase
 * will inherit this method.
 *
 * The %extend rewrites __cmp__ such that two objects of
 * disimilar type can be compared without throwing a TypeError.  For
 * example: the following will return false and not throw an exception:
 *
 *   c = libsbml.Compartment()
 *   n = 5
 *   c == n
 *
 * For more information, see testEquality() in test/TestPython.py
 */

%define SWIGPYTHON__CMP__(CLASS)
%extend CLASS
{
  %pythoncode
  {
    def __eq__(self, rhs):
      if ((self is None) and (rhs is None)): return True
      if ((self is None) or  (rhs is None)): return False
      if (hasattr(self, 'this') and hasattr(rhs, 'this')):
        if (self.this == rhs.this): return True
      return False

    def __ne__(self, rhs):
      if ((self is None) and (rhs is None)): return False
      if ((self is None) or  (rhs is None)): return True
      if (hasattr(self, 'this') and hasattr(rhs, 'this')):
        if (self.this == rhs.this): return False
      return True
  }
}
%enddef

SWIGPYTHON__CMP__(SBase)
SWIGPYTHON__CMP__(SBMLWriter)
SWIGPYTHON__CMP__(SBMLReader)
SWIGPYTHON__CMP__(ASTNode)
SWIGPYTHON__CMP__(CVTerm)
SWIGPYTHON__CMP__(Date)
SWIGPYTHON__CMP__(ModelHistory)
SWIGPYTHON__CMP__(ModelCreator)
SWIGPYTHON__CMP__(XMLNamespaces)
SWIGPYTHON__CMP__(SBMLNamespaces)
SWIGPYTHON__CMP__(XMLAttributes)
SWIGPYTHON__CMP__(XMLToken)
SWIGPYTHON__CMP__(XMLTriple)
SWIGPYTHON__CMP__(XMLError)
SWIGPYTHON__CMP__(XMLErrorLog)
SWIGPYTHON__CMP__(XMLOutputStream)

/**
 * The features directives below override the default SWIG generated
 * code for certain methods.  The idea is to tell SWIG to disown the
 * passed-in object.  The containing object will takeover ownership
 * and delete the object as appropriate.  This avoids a deadly
 * double-delete which can result in a segmentation fault.  For
 * example, each SBase that is appended to a ListOf is subsequently
 * owned by that ListOf.
 */

%define TAKEOVER_OWNERSHIP(METHOD_NAME,ARG_INDEX)
%feature("pythonprepend")
METHOD_NAME
%{
        if args[ARG_INDEX] is not None: args[ARG_INDEX].thisown = 0
%}
%enddef

// ----------------------------------------------------------------------
// ListOf
// ----------------------------------------------------------------------

#if SWIG_VERSION > 0x010336
TAKEOVER_OWNERSHIP(ListOf::appendAndOwn(SBase*),0)
#else
TAKEOVER_OWNERSHIP(ListOf::appendAndOwn(SBase*),1)
#endif

// ----------------------------------------------------------------------
// ASTNode
// ----------------------------------------------------------------------

#if SWIG_VERSION > 0x010336
TAKEOVER_OWNERSHIP(ASTNode::addChild(ASTNode*),0)
TAKEOVER_OWNERSHIP(ASTNode::prependChild(ASTNode*),0)
TAKEOVER_OWNERSHIP(ASTNode::insertChild(unsigned int, ASTNode*),1)
TAKEOVER_OWNERSHIP(ASTNode::replaceChild(unsigned int, ASTNode*),1)
TAKEOVER_OWNERSHIP(ASTNode::addSemanticsAnnotation(XMLNode*),0)
#else
TAKEOVER_OWNERSHIP(ASTNode::addChild(ASTNode*),1)
TAKEOVER_OWNERSHIP(ASTNode::prependChild(ASTNode*),1)
TAKEOVER_OWNERSHIP(ASTNode::insertChild(unsigned int, ASTNode*),2)
TAKEOVER_OWNERSHIP(ASTNode::replaceChild(unsigned int, ASTNode*),2)
TAKEOVER_OWNERSHIP(ASTNode::addSemanticsAnnotation(XMLNode*),1)
#endif

/**
 *
 * Wraps the SBMLConstructorException class (C++ exception defined by libSBML) 
 * as the VaueError class (Python built-in exception).
 *
 * For example, the exception can be catched in Python code as follows:
 *
 * --------------------------------------
 *  try:
 *    s = libsbml.CompartmentType(level,version)
 *  except ValueError, inst:
 *    errmsg = inst.args[0]
 * --------------------------------------
 */

%ignore SBMLConstructorException;

%define SBMLCONSTRUCTOR_EXCEPTION(SBASE_CLASS_NAME)
%exception SBASE_CLASS_NAME {
  try {
    $action
  }
  catch (SBMLConstructorException &e) {
    PyErr_SetString(PyExc_ValueError, const_cast<char*>(e.what()));
    return NULL;
  }
}
%enddef

SBMLCONSTRUCTOR_EXCEPTION(Compartment)
SBMLCONSTRUCTOR_EXCEPTION(CompartmentType)
SBMLCONSTRUCTOR_EXCEPTION(Constraint)
SBMLCONSTRUCTOR_EXCEPTION(Delay)
SBMLCONSTRUCTOR_EXCEPTION(Event)
SBMLCONSTRUCTOR_EXCEPTION(EventAssignment)
SBMLCONSTRUCTOR_EXCEPTION(FunctionDefinition)
SBMLCONSTRUCTOR_EXCEPTION(InitialAssignment)
SBMLCONSTRUCTOR_EXCEPTION(KineticLaw)
SBMLCONSTRUCTOR_EXCEPTION(Model)
SBMLCONSTRUCTOR_EXCEPTION(Parameter)
SBMLCONSTRUCTOR_EXCEPTION(LocalParameter)
SBMLCONSTRUCTOR_EXCEPTION(Reaction)
SBMLCONSTRUCTOR_EXCEPTION(AssignmentRule)
SBMLCONSTRUCTOR_EXCEPTION(AlgebraicRule)
SBMLCONSTRUCTOR_EXCEPTION(RateRule)
SBMLCONSTRUCTOR_EXCEPTION(Species)
SBMLCONSTRUCTOR_EXCEPTION(SpeciesReference)
SBMLCONSTRUCTOR_EXCEPTION(ModifierSpeciesReference)
SBMLCONSTRUCTOR_EXCEPTION(SpeciesType)
SBMLCONSTRUCTOR_EXCEPTION(StoichiometryMath)
SBMLCONSTRUCTOR_EXCEPTION(Trigger)
SBMLCONSTRUCTOR_EXCEPTION(Unit)
SBMLCONSTRUCTOR_EXCEPTION(UnitDefinition)

// ----------------------------------------------------------------------
// SBMLReader
// ----------------------------------------------------------------------


%pythoncode
%{
import sys
import os.path


def conditional_abspath (filename):
  """conditional_abspath (filename) -> filename

  Returns filename with an absolute path prepended, if necessary.
  Some combinations of platforms and underlying XML parsers *require*
  an absolute path to a filename while others do not.  This function
  encapsulates the appropriate logic.  It is used by readSBML() and
  SBMLReader.readSBML().
  """
  if sys.platform.find('cygwin') != -1:
    return filename
  else:
    return os.path.abspath(filename)
%}

%feature("shadow")
SBMLReader::readSBML(const std::string&)
%{
  def readSBML(*args):
    """readSBML(filename) -> SBMLDocument

    Reads an SBML document from the given file.  If filename does not exist
    or is not an SBML file, a fatal error will be logged.  Errors can be
    identified by their unique ids, e.g.:

      reader = libsbml.SBMLReader()
      d      = reader.readSBML(filename)

      if d.getNumErrors() > 0:
        pm = d.getError(0)
        if pm.getErrorId() == libsbml.XMLFileUnreadable
        if pm.getErrorId() == libsbml.XMLTagMismatch: 
    """
    args_copy    = list(args)
    args_copy[1] = conditional_abspath(args[1])
    return _libsbml.SBMLReader_readSBML(*args_copy)
%}


/**
 * Since we cannot seem to "shadow" readSBML() (maybe because it's
 * not a method of some object, but rather a top-level function, we
 * employ the following HACK: Tell SWIG to ignore readSBML and just
 * define it in terms of SBMLReader.readSBML().  This is less than
 * ideal, because the libSBML C/C++ core does essentially the same
 * thing, so now we're repeating ourselves.
 */

%ignore readSBML(const char*);

%pythoncode
%{
def readSBML(*args):
  """readSBML(filename) -> SBMLDocument

  Reads an SBML document from the given file.  If filename does not exist
  or is not an SBML file, a fatal error will be logged.  Errors can be
  identified by their unique ids, e.g.:

    d = readSBML(filename)

    if d.getNumErrors() > 0:
      pm = d.getError(0)
      if pm.getErrorId() == libsbml.XMLFileUnreadable
      if pm.getErrorId() == libsbml.XMLTagMismatch: 

  """
  reader = SBMLReader()
  return reader.readSBML(args[0])
%}


/**
 *  Wraps the following functions by using the corresponding 
 *  ListWrapper<TYPENAME> class.
 *
 *  - List* ModelHistory::getListCreators()
 *  - List* ModelHistory::getListModifiedDates()
 *  - List* SBase::getCVTerms()
 *
 *  ListWrapper<TYPENAME> class is wrapped as TYPENAMEList class.
 *  So, the above functions are wrapped as follows:
 *
 *  - ModelCreatorList ModelHistory.getListCreators()
 *  - DateList         ModelHistory.getListModifiedDates()
 *  - CVTermList       SBase.getCVTerms()
 *
 */

%feature("shadow")
ModelHistory::getListCreators
%{
  def getListCreators(self):
    """
    getListCreators(self) -> ModelCreatorList

    Get the ModelCreatorList of ModelCreator objects in this 
    ModelHistory.

    @return the ModelCreatorList for this ModelHistory.
          

    """
    return _libsbml.ModelHistory_getListCreators(self)
%}

%typemap(out) List* ModelHistory::getListCreators
{
  ListWrapper<ModelCreator> *listw = ($1 != 0) ? new ListWrapper<ModelCreator>($1) : 0;
  $result = SWIG_NewPointerObj(SWIG_as_voidptr(listw), 
#if SWIG_VERSION > 0x010333
                               SWIGTYPE_p_ListWrapperT_ModelCreator_t, 
#else
                               SWIGTYPE_p_ListWrapperTModelCreator_t, 
#endif
                               SWIG_POINTER_OWN |  0 );
}


%feature("shadow")
ModelHistory::getListModifiedDates
%{
  def getListModifiedDates(self):
    """
    getListModifiedDates(self) -> DateList

    Get the DateList of Date objects in this ModelHistory.

    @return the DateList for this ModelHistory.
          

    """
    return _libsbml.ModelHistory_getListModifiedDates(self)
%}

%typemap(out) List* ModelHistory::getListModifiedDates
{
  ListWrapper<Date> *listw = ($1 != 0) ? new ListWrapper<Date>($1) : 0;
  $result = SWIG_NewPointerObj(SWIG_as_voidptr(listw), 
#if SWIG_VERSION > 0x010333
                               SWIGTYPE_p_ListWrapperT_Date_t, 
#else
                               SWIGTYPE_p_ListWrapperTDate_t, 
#endif
                               SWIG_POINTER_OWN |  0 );
}

%feature("shadow")
SBase::getCVTerms
%{
  def getCVTerms(self):
    """
    getCVTerms(self) -> CVTermList

    Get the CVTermList of CVTerm objects in this SBase.

    @return the CVTermList for this SBase.


    """
    return _libsbml.SBase_getCVTerms(self)
%}

%typemap(out) List* SBase::getCVTerms
{
  ListWrapper<CVTerm> *listw = ($1 != 0) ? new ListWrapper<CVTerm>($1) : 0;
  $result = SWIG_NewPointerObj(SWIG_as_voidptr(listw), 
#if SWIG_VERSION > 0x010333
                               SWIGTYPE_p_ListWrapperT_CVTerm_t, 
#else
                               SWIGTYPE_p_ListWrapperTCVTerm_t, 
#endif
                               SWIG_POINTER_OWN |  0 );
}


// ----------------------------------------------------------------------
// Layout Extension
// ----------------------------------------------------------------------


#ifdef USE_LAYOUT
%include layout_local.i
#endif
