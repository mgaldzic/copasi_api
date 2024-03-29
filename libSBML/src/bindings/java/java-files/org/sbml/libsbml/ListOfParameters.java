/* ----------------------------------------------------------------------------
 * This file was automatically generated by SWIG (http://www.swig.org).
 * Version 1.3.40
 *
 * Do not make changes to this file unless you know what you are doing--modify
 * the SWIG interface file instead.
 * ----------------------------------------------------------------------------- */

package org.sbml.libsbml;

/** 
 * LibSBML implementation of SBML's ListOfParameters construct.
 * <p>
 * The various ListOf___ classes in SBML are merely containers used for
 * organizing the main components of an SBML model.  All are derived from
 * the abstract class {@link SBase}, and inherit the various attributes and
 * subelements of {@link SBase}, such as 'metaid' as and 'annotation'.  The
 * ListOf___ classes do not add any attributes of their own.
 * <p>
 * The relationship between the lists and the rest of an SBML model is
 * illustrated by the following (for SBML Level&nbsp;2 Version&nbsp;4):
 * <p>
 * <center><img src='listof-illustration.jpg'></center><br>
 * 
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

public class ListOfParameters extends ListOf {
   private long swigCPtr;

   protected ListOfParameters(long cPtr, boolean cMemoryOwn)
   {
     super(libsbmlJNI.SWIGListOfParametersUpcast(cPtr), cMemoryOwn);
     swigCPtr = cPtr;
   }

   protected static long getCPtr(ListOfParameters obj)
   {
     return (obj == null) ? 0 : obj.swigCPtr;
   }

   protected static long getCPtrAndDisown (ListOfParameters obj)
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
        libsbmlJNI.delete_ListOfParameters(swigCPtr);
      }
      swigCPtr = 0;
    }
    super.delete();
  }

  
  /**
   * Creates and returns a deep copy of this {@link ListOfParameters} instance.
   * <p>
   * @return a (deep) copy of this {@link ListOfParameters}.
   */
 public ListOfParameters cloneObject() {
    long cPtr = libsbmlJNI.ListOfParameters_cloneObject(swigCPtr, this);
    return (cPtr == 0) ? null : new ListOfParameters(cPtr, true);
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
    return libsbmlJNI.ListOfParameters_getTypeCode(swigCPtr, this);
  }

  
  /**
   * Returns the libSBML type code for the objects contained in this {@link ListOf}
   * (i.e., {@link Parameter} objects, if the list is non-empty).
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
    return libsbmlJNI.ListOfParameters_getItemTypeCode(swigCPtr, this);
  }

  
  /**
   * Returns the XML element name of this object.
   * <p>
   * For {@link ListOfParameters}, the XML element name is <code>'listOfParameters'</code>.
   * <p>
   * @return the name of this element, i.e., <code>'listOfParameters'</code>.
   */
 public String getElementName() {
    return libsbmlJNI.ListOfParameters_getElementName(swigCPtr, this);
  }

  
  /**
   * Returns the {@link Parameter} object located at position <code>n</code> within this
   * {@link ListOfParameters} instance.
   * <p>
   * @param n the index number of the {@link Parameter} to get.
   * <p>
   * @return the nth {@link Parameter} in this {@link ListOfParameters}.  If the index <code>n</code>
   * is out of bounds for the length of the list, then <code>NULL</code> is returned.
   * <p>
   * @see #size()
   * @see #get(String sid)
   */
 public Parameter get(long n) {
    long cPtr = libsbmlJNI.ListOfParameters_get__SWIG_0(swigCPtr, this, n);
    return (cPtr == 0) ? null : new Parameter(cPtr, false);
  }

  
  /**
   * Returns the first {@link Parameter} object matching the given identifier.
   * <p>
   * @param sid a string, the identifier of the {@link Parameter} to get.
   * <p>
   * @return the {@link Parameter} object found.  The caller owns the returned
   * object and is responsible for deleting it.  If none of the items have
   * an identifier matching <code>sid</code>, then <code>NULL</code> is returned.
   * <p>
   * @see #get(long n)
   * @see #size()
   */
 public Parameter get(String sid) {
    long cPtr = libsbmlJNI.ListOfParameters_get__SWIG_2(swigCPtr, this, sid);
    return (cPtr == 0) ? null : new Parameter(cPtr, false);
  }

  
  /**
   * Removes the nth item from this {@link ListOfParameters}, and returns a pointer
   * to it.
   * <p>
   * @param n the index of the item to remove
   * <p>
   * @return the item removed.  The caller owns the returned object and is
   * responsible for deleting it.  If the index number <code>n</code> is out of
   * bounds for the length of the list, then <code>NULL</code> is returned.
   * <p>
   * @see #size()
   */
 public Parameter remove(long n) {
    long cPtr = libsbmlJNI.ListOfParameters_remove__SWIG_0(swigCPtr, this, n);
    return (cPtr == 0) ? null : new Parameter(cPtr, true);
  }

  
  /**
   * Removes the first {@link Parameter} object in this {@link ListOfParameters}
   * matching the given identifier, and returns a pointer to it.
   * <p>
   * @param sid the identifier of the item to remove.
   * <p>
   * @return the item removed.  The caller owns the returned object and is
   * responsible for deleting it.  If none of the items have an identifier
   * matching <code>sid</code>, then <code>NULL</code> is returned.
   */
 public Parameter remove(String sid) {
    long cPtr = libsbmlJNI.ListOfParameters_remove__SWIG_1(swigCPtr, this, sid);
    return (cPtr == 0) ? null : new Parameter(cPtr, true);
  }

  public ListOfParameters() {
    this(libsbmlJNI.new_ListOfParameters(), true);
  }

}
