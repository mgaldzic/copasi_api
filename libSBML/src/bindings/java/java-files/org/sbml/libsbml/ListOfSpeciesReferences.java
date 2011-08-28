/* ----------------------------------------------------------------------------
 * This file was automatically generated by SWIG (http://www.swig.org).
 * Version 1.3.40
 *
 * Do not make changes to this file unless you know what you are doing--modify
 * the SWIG interface file instead.
 * ----------------------------------------------------------------------------- */

package org.sbml.libsbml;

/** 
 * LibSBML implementation of SBML's ListOfSpeciesReferences construct.
 * <p>
 * The {@link ListOfSpeciesReferences} class is used to store lists of reactants
 * and products in a {@link Reaction} object.
 * <p>
 * As with the various other ListOf___ classes in SBML, the
 * {@link ListOfSpeciesReferences} is merely a container used for organizing
 * instances of other objects, in this case {@link SpeciesReference} objects.
 * {@link ListOfSpeciesReferences} is derived from the abstract class {@link SBase}, and
 * inherit the various attributes and subelements of {@link SBase}, such as
 * 'metaid' as and 'annotation'.  The ListOf___ classes do not add any
 * attributes of their own.
 */

public class ListOfSpeciesReferences extends ListOf {
   private long swigCPtr;

   protected ListOfSpeciesReferences(long cPtr, boolean cMemoryOwn)
   {
     super(libsbmlJNI.SWIGListOfSpeciesReferencesUpcast(cPtr), cMemoryOwn);
     swigCPtr = cPtr;
   }

   protected static long getCPtr(ListOfSpeciesReferences obj)
   {
     return (obj == null) ? 0 : obj.swigCPtr;
   }

   protected static long getCPtrAndDisown (ListOfSpeciesReferences obj)
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
        libsbmlJNI.delete_ListOfSpeciesReferences(swigCPtr);
      }
      swigCPtr = 0;
    }
    super.delete();
  }

  
  /**
   * Creates a new, empty {@link ListOfSpeciesReferences}.
   */
 public ListOfSpeciesReferences() {
    this(libsbmlJNI.new_ListOfSpeciesReferences(), true);
  }

  
  /**
   * Creates and returns a deep copy of this {@link ListOfSpeciesReferences}
   * instance.
   * <p>
   * @return a (deep) copy of this {@link ListOfSpeciesReferences}.
   */
 public ListOfSpeciesReferences cloneObject() {
    long cPtr = libsbmlJNI.ListOfSpeciesReferences_cloneObject(swigCPtr, this);
    return (cPtr == 0) ? null : new ListOfSpeciesReferences(cPtr, true);
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
   * @return the SBML type code for this object, or @link SBMLTypeCode_t#SBML_UNKNOWN SBML_UNKNOWN@endlink (default).
   * <p>
   * @see #getElementName()
   */
 public int getTypeCode() {
    return libsbmlJNI.ListOfSpeciesReferences_getTypeCode(swigCPtr, this);
  }

  
  /**
   * Returns the libSBML type code for the objects contained in this {@link ListOf}
   * (i.e., {@link SpeciesReference} objects, if the list is non-empty).
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
   * instance, or @link SBMLTypeCode_t#SBML_UNKNOWN SBML_UNKNOWN@endlink (default).
   * <p>
   * @see #getElementName()
   */
 public int getItemTypeCode() {
    return libsbmlJNI.ListOfSpeciesReferences_getItemTypeCode(swigCPtr, this);
  }

  
  /**
   * Returns the XML element name of this object.
   * <p>
   * For {@link ListOfSpeciesReferences}, the XML element name is 
   * <code>'listOfSpeciesReferences'</code>.
   * <p>
   * @return the name of this element, i.e., <code>'listOfSpeciesReferences'</code>.
   */
 public String getElementName() {
    return libsbmlJNI.ListOfSpeciesReferences_getElementName(swigCPtr, this);
  }

  
  /**
   * Get a {@link SpeciesReference} from the {@link ListOfSpeciesReferences}.
   * <p>
   * @param n the index number of the {@link SpeciesReference} to get.
   * <p>
   * @return the nth {@link SpeciesReference} in this {@link ListOfSpeciesReferences}.
   * <p>
   * @see #size()
   */
 public SimpleSpeciesReference get(long n) {
  return (SimpleSpeciesReference) libsbml.DowncastSBase(libsbmlJNI.ListOfSpeciesReferences_get__SWIG_0(swigCPtr, this, n), false);
}

  
  /**
   * Get a {@link SpeciesReference} from the {@link ListOfSpeciesReferences}
   * based on its identifier.
   * <p>
   * @param sid a string representing the identifier 
   * of the {@link SpeciesReference} to get.
   * <p>
   * @return {@link SpeciesReference} in this {@link ListOfSpeciesReferences}
   * with the given id or <code>NULL</code> if no such
   * {@link SpeciesReference} exists.
   * <p>
   * @see #get(long n)
   * @see #size()
   */
 public SimpleSpeciesReference get(String sid) {
  return (SimpleSpeciesReference) libsbml.DowncastSBase(libsbmlJNI.ListOfSpeciesReferences_get__SWIG_2(swigCPtr, this, sid), false);
}

  
  /**
   * Removes the nth item from this {@link ListOfSpeciesReferences} items and returns a pointer to
   * it.
   * <p>
   * The caller owns the returned item and is responsible for deleting it.
   * <p>
   * @param n the index of the item to remove
   * <p>
   * @see #size()
   */
 public SimpleSpeciesReference remove(long n) {
  return (SimpleSpeciesReference) libsbml.DowncastSBase(libsbmlJNI.ListOfSpeciesReferences_remove__SWIG_0(swigCPtr, this, n), true);
}

  
  /**
   * Removes item in this {@link ListOfSpeciesReferences} items with the given identifier.
   * <p>
   * The caller owns the returned item and is responsible for deleting it.
   * If none of the items in this list have the identifier <code>sid</code>, then 
   * <code>NULL</code> is returned.
   * <p>
   * @param sid the identifier of the item to remove
   * <p>
   * @return the item removed.  As mentioned above, the caller owns the
   * returned item.
   */
 public SimpleSpeciesReference remove(String sid) {
  return (SimpleSpeciesReference) libsbml.DowncastSBase(libsbmlJNI.ListOfSpeciesReferences_remove__SWIG_1(swigCPtr, this, sid), true);
}

}
