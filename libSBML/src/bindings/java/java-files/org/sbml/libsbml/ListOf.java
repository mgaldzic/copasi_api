/* ----------------------------------------------------------------------------
 * This file was automatically generated by SWIG (http://www.swig.org).
 * Version 1.3.40
 *
 * Do not make changes to this file unless you know what you are doing--modify
 * the SWIG interface file instead.
 * ----------------------------------------------------------------------------- */

package org.sbml.libsbml;

/** 
 * Parent class for the various SBML 'ListOfXYZ' classes.
 * <p>
 * <em style='color: #555'>
This class of objects is defined by libSBML only and has no direct
equivalent in terms of SBML components.  This class is not prescribed by
the SBML specifications, although it is used to implement features
defined in SBML.
</em>

 * <p>
 */

public class ListOf extends SBase {
   private long swigCPtr;

   protected ListOf(long cPtr, boolean cMemoryOwn)
   {
     super(libsbmlJNI.SWIGListOfUpcast(cPtr), cMemoryOwn);
     swigCPtr = cPtr;
   }

   protected static long getCPtr(ListOf obj)
   {
     return (obj == null) ? 0 : obj.swigCPtr;
   }

   protected static long getCPtrAndDisown (ListOf obj)
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
        libsbmlJNI.delete_ListOf(swigCPtr);
      }
      swigCPtr = 0;
    }
    super.delete();
  }

  
  /**
   * Creates a new {@link ListOf}.
   */
 public ListOf() {
    this(libsbmlJNI.new_ListOf__SWIG_0(), true);
  }

  
  /**
   * Copy constructor.  Creates a copy of this {@link ListOf}.
   */
 public ListOf(ListOf orig) {
    this(libsbmlJNI.new_ListOf__SWIG_1(ListOf.getCPtr(orig), orig), true);
  }

  
  /**
   * Creates and returns a deep copy of this {@link ListOf}.
   * <p>
   * @return a (deep) copy of this {@link ListOf}.
   */
 public SBase cloneObject() {
  return libsbml.DowncastSBase(libsbmlJNI.ListOf_cloneObject(swigCPtr, this), true);
}

  
  /**
   * Adds item to the end of this {@link ListOf}.
   * <p>
   * This variant of the method makes a clone of the <code>item</code> handed to it.
   * This means that when the {@link ListOf} is destroyed, the original items will
   * not be destroyed.
   * <p>
   * @param item the item to be added to the list.
   * <p>
   * @see #appendAndOwn(SBase  item)
   */
 public void append(SBase item) {
    libsbmlJNI.ListOf_append(swigCPtr, this, SBase.getCPtr(item), item);
  }

  
  /**
   * Adds item to the end of this {@link ListOf}.
   * <p>
   * This variant of the method does not clone the <code>item</code> handed to it;
   * instead, it assumes ownership of it.  This means that when the {@link ListOf}
   * is destroyed, the item will be destroyed along with it.
   * <p>
   * @param item the item to be added to the list.
   * <p>
   * @see #append(SBase  item)
   */
 public void appendAndOwn(SBase item) {
    libsbmlJNI.ListOf_appendAndOwn(swigCPtr, this, SBase.getCPtrAndDisown(item), item);
  }

  
  /**
   * Get an item from the list.
   * <p>
   * @param n the index number of the item to get.
   * <p>
   * @return the nth item in this {@link ListOf} items.
   * <p>
   * @see #size()
   */
 public SBase get(long n) {
  return libsbml.DowncastSBase(libsbmlJNI.ListOf_get__SWIG_0(swigCPtr, this, n), false);
}

  
  /**
   * Removes all items in this {@link ListOf} object.
   * <p>
   * If doDelete is true (default), all items in this {@link ListOf} object are deleted
   * and cleared, and thus the caller doesn't have to delete those items.
   * Otherwise, all items are just cleared from this {@link ListOf} object and the caller 
   * is responsible for deleting all items (In this case, pointers to all items 
   * should be stored elsewhere before calling this function by the caller).
   * <p>
   * @param doDelete if true (default), all items are deleted and cleared.
   * Otherwise, all items are just cleared and not deleted. 
   */
 public void clear(boolean doDelete) {
    libsbmlJNI.ListOf_clear__SWIG_0(swigCPtr, this, doDelete);
  }

  
  /**
   * Removes all items in this {@link ListOf} object.
   * <p>
   * If doDelete is true (default), all items in this {@link ListOf} object are deleted
   * and cleared, and thus the caller doesn't have to delete those items.
   * Otherwise, all items are just cleared from this {@link ListOf} object and the caller 
   * is responsible for deleting all items (In this case, pointers to all items 
   * should be stored elsewhere before calling this function by the caller).
   * <p>
   * @param doDelete if true (default), all items are deleted and cleared.
   * Otherwise, all items are just cleared and not deleted. 
   */
 public void clear() {
    libsbmlJNI.ListOf_clear__SWIG_1(swigCPtr, this);
  }

  
  /**
   * Removes the nth item from this {@link ListOf} items and returns a pointer to
   * it.
   * <p>
   * The caller owns the returned item and is responsible for deleting it.
   * <p>
   * @param n the index of the item to remove
   * <p>
   * @see #size()
   */
 public SBase remove(long n) {
  return libsbml.DowncastSBase(libsbmlJNI.ListOf_remove(swigCPtr, this, n), true);
}

  
  /**
   * Get the size of this {@link ListOf}.
   * <p>
   * @return the number of items in this {@link ListOf} items.
   */
 public long size() {
    return libsbmlJNI.ListOf_size(swigCPtr, this);
  }

  
  /**
   * Returns the libSBML type code for this object, namely, 
   * <code>SBML_LIST_OF</code>.
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
    return libsbmlJNI.ListOf_getTypeCode(swigCPtr, this);
  }

  
  /**
   * Get the type code of the objects contained in this {@link ListOf}.
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
   */
 public int getItemTypeCode() {
    return libsbmlJNI.ListOf_getItemTypeCode(swigCPtr, this);
  }

  
  /**
   * Returns the XML element name of this object, which for {@link ListOf}, is
   * always <code>'listOf'</code>.
   * <p>
   * @return the XML name of this element.
   */
 public String getElementName() {
    return libsbmlJNI.ListOf_getElementName(swigCPtr, this);
  }

}