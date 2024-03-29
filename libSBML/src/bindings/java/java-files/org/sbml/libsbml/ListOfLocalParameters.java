/* ----------------------------------------------------------------------------
 * This file was automatically generated by SWIG (http://www.swig.org).
 * Version 1.3.40
 *
 * Do not make changes to this file unless you know what you are doing--modify
 * the SWIG interface file instead.
 * ----------------------------------------------------------------------------- */

package org.sbml.libsbml;

/** 
 * LibSBML implementation of SBML's ListOfLocalParameters construct.
 * <p>
 * The various ListOf___ classes in SBML are merely containers used for
 * organizing the main components of an SBML model.  All are derived from
 * the abstract class {@link SBase}, and inherit the various attributes and
 * subelements of {@link SBase}, such as 'metaid' as and 'annotation'.  The
 * ListOf___ classes do not add any attributes of their own.
 * <p>
 * {@link ListOfLocalParameters} is a subsidiary object class used only within
 * {@link KineticLaw}.  A {@link KineticLaw} object can have a single object of class
 * {@link ListOfLocalParameters} containing a set of local parameters used in that
 * kinetic law definition.
 * <p>
 * Readers may wonder about the motivations for using the ListOf___
 * containers.  A simpler approach in XML might be to place the components
 * all directly at the top level of the model definition.  We chose instead
 * to group them within XML elements named after {@link ListOf}<em>Classname</em>,
 * in part because we believe this helps organize the components and makes
 * visual reading of models in XML easier.  More importantly, the fact that
 * the container classes are derived from {@link SBase} means that software tools
 * can add information about the lists themselves into each list
 * container's 'annotation'.
 * <p>
 * @see ListOfFunctionDefinitions
 * @see ListOfUnitDefinitions
 * @see ListOfCompartmentTypes
 * @see ListOfSpeciesTypes
 * @see ListOfCompartments
 * @see ListOfSpecies
 * @see ListOfParameters
 * @see ListOfInitialAssignments
 * @see ListOfRules
 * @see ListOfConstraints
 * @see ListOfReactions
 * @see ListOfEvents
 */

public class ListOfLocalParameters extends ListOfParameters {
   private long swigCPtr;

   protected ListOfLocalParameters(long cPtr, boolean cMemoryOwn)
   {
     super(libsbmlJNI.SWIGListOfLocalParametersUpcast(cPtr), cMemoryOwn);
     swigCPtr = cPtr;
   }

   protected static long getCPtr(ListOfLocalParameters obj)
   {
     return (obj == null) ? 0 : obj.swigCPtr;
   }

   protected static long getCPtrAndDisown (ListOfLocalParameters obj)
   {
     long ptr = 0;

     if (obj != null)
     {
       ptr             = obj.swigCPtr;
       obj.swigCMemOwn = false;
     }

     return ptr;
   }

  protected void finalize() {
    delete();
  }

  public synchronized void delete() {
    if (swigCPtr != 0) {
      if (swigCMemOwn) {
        swigCMemOwn = false;
        libsbmlJNI.delete_ListOfLocalParameters(swigCPtr);
      }
      swigCPtr = 0;
    }
    super.delete();
  }

  
  /**
   * Creates and returns a deep copy of this {@link ListOfLocalParameters} object.
   * <p>
   * @return a (deep) copy of this {@link ListOfLocalParameters}.
   */
 public ListOfLocalParameters cloneObject() {
    long cPtr = libsbmlJNI.ListOfLocalParameters_cloneObject(swigCPtr, this);
    return (cPtr == 0) ? null : new ListOfLocalParameters(cPtr, true);
  }

  
  /**
   * Returns the libSBML type code for this SBML object.
   * <p>
   * LibSBML attaches an
   * identifying code to every kind of SBML object.  These are known as
   * <em>SBML type codes</em>.  In other languages, the set of type codes
   * is stored in an enumeration; in the Java language interface for
   * libSBML, the type codes are defined as static integer constants in
   * interface class {@link libsbmlConstants}.  The names of the type codes
   * all begin with the characters <code>SBML_</code>. 
   * <p>
   * @return the SBML type code for this object, or @link
   * SBMLTypeCode_t#SBML_UNKNOWN SBML_UNKNOWN@endlink (default).
   * <p>
   * @see #getElementName()
   */
 public int getTypeCode() {
    return libsbmlJNI.ListOfLocalParameters_getTypeCode(swigCPtr, this);
  }

  
  /**
   * Returns the libSBML type code for the objects contained in this {@link ListOf}
   * (i.e., {@link LocalParameter} objects, if the list is non-empty).
   * <p>
   * LibSBML attaches an
   * identifying code to every kind of SBML object.  These are known as
   * <em>SBML type codes</em>.  In other languages, the set of type codes
   * is stored in an enumeration; in the Java language interface for
   * libSBML, the type codes are defined as static integer constants in
   * interface class {@link libsbmlConstants}.  The names of the type codes
   * all begin with the characters <code>SBML_</code>. 
   * <p>
   * @return the SBML type code for the objects contained in this {@link ListOf}
   * instance, or @link SBMLTypeCode_t#SBML_UNKNOWN SBML_UNKNOWN@endlink
   * (default).
   * <p>
   * @see #getElementName()
   */
 public int getItemTypeCode() {
    return libsbmlJNI.ListOfLocalParameters_getItemTypeCode(swigCPtr, this);
  }

  
  /**
   * Returns the XML element name of this object.
   * <p>
   * For {@link ListOfLocalParameters}, the XML element name is <code>'listOfLocalParameters'</code>.
   * <p>
   * @return the name of this element, i.e., <code>'listOfLocalParameters'</code>.
   */
 public String getElementName() {
    return libsbmlJNI.ListOfLocalParameters_getElementName(swigCPtr, this);
  }

  
  /**
   * Returns the {@link LocalParameter} object located at position <code>n</code> within this
   * {@link ListOfLocalParameters} instance.
   * <p>
   * @param n the index number of the {@link LocalParameter} to get.
   * <p>
   * @return the nth {@link LocalParameter} in this {@link ListOfLocalParameters}.  If the
   * index <code>n</code> is out of bounds for the length of the list, then <code>NULL</code>
   * is returned.
   * <p>
   * @see #size()
   * @see #get(String sid)
   */
 public LocalParameter get(long n) {
    long cPtr = libsbmlJNI.ListOfLocalParameters_get__SWIG_0(swigCPtr, this, n);
    return (cPtr == 0) ? null : new LocalParameter(cPtr, false);
  }

  
  /**
   * Returns the first {@link LocalParameter} object matching the given identifier.
   * <p>
   * @param sid a string, the identifier of the {@link LocalParameter} to get.
   * <p>
   * @return the {@link LocalParameter} object found.  The caller owns the returned
   * object and is responsible for deleting it.  If none of the items have
   * an identifier matching <code>sid</code>, then <code>NULL</code> is returned.
   * <p>
   * @see #get(long n)
   * @see #size()
   */
 public LocalParameter get(String sid) {
    long cPtr = libsbmlJNI.ListOfLocalParameters_get__SWIG_2(swigCPtr, this, sid);
    return (cPtr == 0) ? null : new LocalParameter(cPtr, false);
  }

  
  /**
   * Removes the nth item from this {@link ListOfLocalParameters}, and returns a
   * pointer to it.
   * <p>
   * @param n the index of the item to remove.  
   * <p>
   * @return the item removed.  The caller owns the returned object and is
   * responsible for deleting it.  If the index number <code>n</code> is out of
   * bounds for the length of the list, then <code>NULL</code> is returned.
   * <p>
   * @see #size()
   * @see #remove(String sid)
   */
 public LocalParameter remove(long n) {
    long cPtr = libsbmlJNI.ListOfLocalParameters_remove__SWIG_0(swigCPtr, this, n);
    return (cPtr == 0) ? null : new LocalParameter(cPtr, true);
  }

  
  /**
   * Removes the first {@link LocalParameter} object in this {@link ListOfLocalParameters}
   * matching the given identifier, and returns a pointer to it.
   * <p>
   * @param sid the identifier of the item to remove.
   * <p>
   * @return the item removed.  The caller owns the returned object and is
   * responsible for deleting it.  If none of the items have an identifier
   * matching <code>sid</code>, then <code>NULL</code> is returned.
   */
 public LocalParameter remove(String sid) {
    long cPtr = libsbmlJNI.ListOfLocalParameters_remove__SWIG_1(swigCPtr, this, sid);
    return (cPtr == 0) ? null : new LocalParameter(cPtr, true);
  }

  public ListOfLocalParameters() {
    this(libsbmlJNI.new_ListOfLocalParameters(), true);
  }

}
