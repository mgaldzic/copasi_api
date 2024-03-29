/* ----------------------------------------------------------------------------
 * This file was automatically generated by SWIG (http://www.swig.org).
 * Version 1.3.40
 *
 * Do not make changes to this file unless you know what you are doing--modify
 * the SWIG interface file instead.
 * ----------------------------------------------------------------------------- */

package org.sbml.libsbml;

/** 
 * Class to store SBML level, version and namespace information.
 * <p>
 * <em style='color: #555'>
This class of objects is defined by libSBML only and has no direct
equivalent in terms of SBML components.  This class is not prescribed by
the SBML specifications, although it is used to implement features
defined in SBML.
</em>

 * <p>
 * There are differences in the definitions of components between different
 * SBML Levels, as well as Versions within Levels.  For example, the
 * 'sboTerm' attribute was not introduced until Level&nbsp;2
 * Version&nbsp;2, and then only on certain component classes; the SBML
 * Level&nbsp;2 Version&nbsp;3 specification moved the 'sboTerm' attribute
 * to the {@link SBase} class, thereby allowing nearly all components to have {@link SBO}
 * annotations.  As a result of differences such as those, libSBML needs to
 * track the SBML Level and Version of every object created.
 * <p>
 * The purpose of the {@link SBMLNamespaces} object class is to make it easier to
 * communicate SBML Level and Version data between libSBML constructors and
 * other methods.  The {@link SBMLNamespaces} object class tracks 3-tuples
 * (triples) consisting of SBML Level, Version, and the corresponding SBML
 * XML namespace.  (The plural name is not a mistake, because in SBML
 * Level&nbsp;3, objects may have extensions added by Level&nbsp;3 packages
 * used by a given model; however, until the introduction of SBML
 * Level&nbsp;3, the {@link SBMLNamespaces} object only records one SBML
 * Level/Version/namespace combination at a time.)  Most constructors for
 * SBML objects in libSBML take a {@link SBMLNamespaces} object as an argument,
 * thereby allowing the constructor to produce the proper combination of
 * attributes and other internal data structures for the given SBML
 * Level and Version.
 */

public class SBMLNamespaces {
   private long swigCPtr;
   protected boolean swigCMemOwn;

   protected SBMLNamespaces(long cPtr, boolean cMemoryOwn)
   {
     swigCMemOwn = cMemoryOwn;
     swigCPtr    = cPtr;
   }

   protected static long getCPtr(SBMLNamespaces obj)
   {
     return (obj == null) ? 0 : obj.swigCPtr;
   }

   protected static long getCPtrAndDisown (SBMLNamespaces obj)
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
        libsbmlJNI.delete_SBMLNamespaces(swigCPtr);
      }
      swigCPtr = 0;
    }
  }

  /**
   * Equality comparison method for SBMLNamespaces.
   * <p>
   * Because the Java methods for libSBML are actually wrappers around code
   * implemented in C++ and C, certain operations will not behave as
   * expected.  Equality comparison is one such case.  An instance of a
   * libSBML object class is actually a <em>proxy object</em>
   * wrapping the real underlying C/C++ object.  The normal <code>==</code>
   * equality operator in Java will <em>only compare the Java proxy objects</em>,
   * not the underlying native object.  The result is almost never what you
   * want in practical situations.  Unfortunately, Java does not provide a
   * way to override <code>==</code>.
   *  <p>
   * The alternative that must be followed is to use the
   * <code>equals()</code> method.  The <code>equals</code> method on this
   * class overrides the default java.lang.Object one, and performs an
   * intelligent comparison of instances of objects of this class.  The
   * result is an assessment of whether two libSBML Java objects are truly 
   * the same underlying native-code objects.
   *  <p>
   * The use of this method in practice is the same as the use of any other
   * Java <code>equals</code> method.  For example,
   * <em>a</em><code>.equals(</code><em>b</em><code>)</code> returns
   * <code>true</code> if <em>a</em> and <em>b</em> are references to the
   * same underlying object.
   *
   * @param sb a reference to an object to which the current object
   * instance will be compared
   *
   * @return <code>true</code> if <code>sb</code> refers to the same underlying 
   * native object as this one, <code>false</code> otherwise
   */
  public boolean equals(Object sb)
  {
    if ( this == sb ) 
    {
      return true;
    }
    return swigCPtr == getCPtr((SBMLNamespaces)(sb));
  }

  /**
   * Returns a hashcode for this SBMLNamespaces object.
   *
   * @return a hash code usable by Java methods that need them.
   */
  public int hashCode()
  {
    return (int)(swigCPtr^(swigCPtr>>>32));
  }

  
  /**
   * Creates a new {@link SBMLNamespaces} object corresponding to the given SBML
   * <code>level</code> and <code>version</code>.
   * <p>
   * {@link SBMLNamespaces} objects are used in libSBML to communicate SBML Level
   * and Version data between constructors and other methods.  The
   * {@link SBMLNamespaces} object class tracks 3-tuples (triples) consisting of
   * SBML Level, Version, and the corresponding SBML XML namespace.  Most
   * constructors for SBML objects in libSBML take a {@link SBMLNamespaces} object
   * as an argument, thereby allowing the constructor to produce the proper
   * combination of attributes and other internal data structures for the
   * given SBML Level and Version.
   * <p>
   * The plural name '{@link SBMLNamespaces}' is not a mistake, because in SBML
   * Level&nbsp;3, objects may have extensions added by Level&nbsp;3
   * packages used by a given model; however, until the introduction of
   * SBML Level&nbsp;3, the {@link SBMLNamespaces} object only records one SBML
   * Level/Version/namespace combination at a time.
   * <p>
   * @param level the SBML level
   * @param version the SBML version
   * <p>
   * @docnote The native C++ implementation of this method defines a default argument
value. In the documentation generated for different libSBML language
bindings, you may or may not see corresponding arguments in the method
declarations. For example, in Java, a default argument is handled by
declaring two separate methods, with one of them having the argument and
the other one lacking the argument. However, the libSBML documentation will
be <em>identical</em> for both methods. Consequently, if you are reading
this and do not see an argument even though one is described, please look
for descriptions of other variants of this method near where this one
appears in the documentation.

   */
 public SBMLNamespaces(long level, long version) {
    this(libsbmlJNI.new_SBMLNamespaces__SWIG_0(level, version), true);
  }

  
  /**
   * Creates a new {@link SBMLNamespaces} object corresponding to the given SBML
   * <code>level</code> and <code>version</code>.
   * <p>
   * {@link SBMLNamespaces} objects are used in libSBML to communicate SBML Level
   * and Version data between constructors and other methods.  The
   * {@link SBMLNamespaces} object class tracks 3-tuples (triples) consisting of
   * SBML Level, Version, and the corresponding SBML XML namespace.  Most
   * constructors for SBML objects in libSBML take a {@link SBMLNamespaces} object
   * as an argument, thereby allowing the constructor to produce the proper
   * combination of attributes and other internal data structures for the
   * given SBML Level and Version.
   * <p>
   * The plural name '{@link SBMLNamespaces}' is not a mistake, because in SBML
   * Level&nbsp;3, objects may have extensions added by Level&nbsp;3
   * packages used by a given model; however, until the introduction of
   * SBML Level&nbsp;3, the {@link SBMLNamespaces} object only records one SBML
   * Level/Version/namespace combination at a time.
   * <p>
   * @param level the SBML level
   * @param version the SBML version
   * <p>
   * @docnote The native C++ implementation of this method defines a default argument
value. In the documentation generated for different libSBML language
bindings, you may or may not see corresponding arguments in the method
declarations. For example, in Java, a default argument is handled by
declaring two separate methods, with one of them having the argument and
the other one lacking the argument. However, the libSBML documentation will
be <em>identical</em> for both methods. Consequently, if you are reading
this and do not see an argument even though one is described, please look
for descriptions of other variants of this method near where this one
appears in the documentation.

   */
 public SBMLNamespaces(long level) {
    this(libsbmlJNI.new_SBMLNamespaces__SWIG_1(level), true);
  }

  
  /**
   * Creates a new {@link SBMLNamespaces} object corresponding to the given SBML
   * <code>level</code> and <code>version</code>.
   * <p>
   * {@link SBMLNamespaces} objects are used in libSBML to communicate SBML Level
   * and Version data between constructors and other methods.  The
   * {@link SBMLNamespaces} object class tracks 3-tuples (triples) consisting of
   * SBML Level, Version, and the corresponding SBML XML namespace.  Most
   * constructors for SBML objects in libSBML take a {@link SBMLNamespaces} object
   * as an argument, thereby allowing the constructor to produce the proper
   * combination of attributes and other internal data structures for the
   * given SBML Level and Version.
   * <p>
   * The plural name '{@link SBMLNamespaces}' is not a mistake, because in SBML
   * Level&nbsp;3, objects may have extensions added by Level&nbsp;3
   * packages used by a given model; however, until the introduction of
   * SBML Level&nbsp;3, the {@link SBMLNamespaces} object only records one SBML
   * Level/Version/namespace combination at a time.
   * <p>
   * @param level the SBML level
   * @param version the SBML version
   * <p>
   * @docnote The native C++ implementation of this method defines a default argument
value. In the documentation generated for different libSBML language
bindings, you may or may not see corresponding arguments in the method
declarations. For example, in Java, a default argument is handled by
declaring two separate methods, with one of them having the argument and
the other one lacking the argument. However, the libSBML documentation will
be <em>identical</em> for both methods. Consequently, if you are reading
this and do not see an argument even though one is described, please look
for descriptions of other variants of this method near where this one
appears in the documentation.

   */
 public SBMLNamespaces() {
    this(libsbmlJNI.new_SBMLNamespaces__SWIG_2(), true);
  }

  
  /**
   * Copy constructor; creates a copy of a {@link SBMLNamespaces}.
   * <p>
   * @param orig the {@link SBMLNamespaces} instance to copy.
   */
 public SBMLNamespaces(SBMLNamespaces orig) {
    this(libsbmlJNI.new_SBMLNamespaces__SWIG_3(SBMLNamespaces.getCPtr(orig), orig), true);
  }

  
  /**
   * Creates and returns a deep copy of this {@link SBMLNamespaces}.
   * <p>
   * @return a (deep) copy of this {@link SBMLNamespaces}.
   */
 public SBMLNamespaces cloneObject() {
    long cPtr = libsbmlJNI.SBMLNamespaces_cloneObject(swigCPtr, this);
    return (cPtr == 0) ? null : new SBMLNamespaces(cPtr, true);
  }

  
  /**
   * Returns a string representing the SBML XML namespace for the 
   * given <code>level</code> and <code>version</code> of SBML.
   * <p>
   * @param level the SBML level
   * @param version the SBML version
   * <p>
   * @return a string representing the SBML namespace that reflects the
   * SBML Level and Version specified.
   */
 public static String getSBMLNamespaceURI(long level, long version) {
    return libsbmlJNI.SBMLNamespaces_getSBMLNamespaceURI(level, version);
  }

  
  /**
   * Get the SBML Level of this {@link SBMLNamespaces} object.
   * <p>
   * @return the SBML Level of this {@link SBMLNamespaces} object.
   */
 public long getLevel() {
    return libsbmlJNI.SBMLNamespaces_getLevel__SWIG_0(swigCPtr, this);
  }

  
  /**
   * Get the SBML Version of this {@link SBMLNamespaces} object.
   * <p>
   * @return the SBML Version of this {@link SBMLNamespaces} object.
   */
 public long getVersion() {
    return libsbmlJNI.SBMLNamespaces_getVersion__SWIG_0(swigCPtr, this);
  }

  
  /**
   * Get the XML namespaces list for this {@link SBMLNamespaces} object.
   * <p>
   * The plural is not a mistake, because in SBML Level&nbsp;3, objects may
   * have extensions added by Level&nbsp;3 packages used by a given model,
   * and therefore there may be multiple XML namespaces involved too.
   * However, until the introduction of SBML Level&nbsp;3, the
   * {@link SBMLNamespaces} object only records one SBML Level/Version/namespace
   * combination at a time, and so this method will also only return
   * a list of one item.
   * <p>
   * @return the XML namespaces of this {@link SBMLNamespaces} object.
   */
 public XMLNamespaces getNamespaces() {
    long cPtr = libsbmlJNI.SBMLNamespaces_getNamespaces__SWIG_0(swigCPtr, this);
    return (cPtr == 0) ? null : new XMLNamespaces(cPtr, false);
  }

  
  /**
   * Add the XML namespaces list to the set of namespaces
   * within this {@link SBMLNamespaces} object.
   * <p>
   * @param xmlns the XML namespaces to be added.
   */
 public void addNamespaces(XMLNamespaces xmlns) {
    libsbmlJNI.SBMLNamespaces_addNamespaces(swigCPtr, this, XMLNamespaces.getCPtr(xmlns), xmlns);
  }

}
