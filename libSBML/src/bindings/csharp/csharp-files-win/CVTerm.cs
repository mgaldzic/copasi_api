/* ----------------------------------------------------------------------------
 * This file was automatically generated by SWIG (http://www.swig.org).
 * Version 1.3.40
 *
 * Do not make changes to this file unless you know what you are doing--modify
 * the SWIG interface file instead.
 * ----------------------------------------------------------------------------- */

namespace libsbml {

using System;
using System.Runtime.InteropServices;

public class CVTerm : IDisposable {
	private HandleRef swigCPtr;
	protected bool swigCMemOwn;
	
	internal CVTerm(IntPtr cPtr, bool cMemoryOwn)
	{
		swigCMemOwn = cMemoryOwn;
		swigCPtr    = new HandleRef(this, cPtr);
	}
	
	internal static HandleRef getCPtr(CVTerm obj)
	{
		return (obj == null) ? new HandleRef(null, IntPtr.Zero) : obj.swigCPtr;
	}
	
	internal static HandleRef getCPtrAndDisown (CVTerm obj)
	{
		HandleRef ptr = new HandleRef(null, IntPtr.Zero);
		
		if (obj != null)
		{
			ptr             = obj.swigCPtr;
			obj.swigCMemOwn = false;
		}
		
		return ptr;
	}

  ~CVTerm() {
    Dispose();
  }

  public virtual void Dispose() {
    lock(this) {
      if (swigCPtr.Handle != IntPtr.Zero) {
        if (swigCMemOwn) {
          swigCMemOwn = false;
          libsbmlPINVOKE.delete_CVTerm(swigCPtr);
        }
        swigCPtr = new HandleRef(null, IntPtr.Zero);
      }
      GC.SuppressFinalize(this);
    }
  }

  public static bool operator==(CVTerm lhs, CVTerm rhs)
  {
    if((Object)lhs == (Object)rhs)
    {
      return true;
    }

    if( ((Object)lhs == null) || ((Object)rhs == null) )
    {
      return false;
    }

    return (getCPtr(lhs).Handle.ToString() == getCPtr(rhs).Handle.ToString());
  }

  public static bool operator!=(CVTerm lhs, CVTerm rhs)
  {
    return !(lhs == rhs);
  }

  public override bool Equals(Object sb)
  {
    if ( ! (sb is CVTerm) )
    {
      return false;
    }

    return this == (CVTerm)sb;
  }

  public override int GetHashCode()
  {
    return swigCPtr.Handle.ToInt32();
  }

  public CVTerm(int type) : this(libsbmlPINVOKE.new_CVTerm__SWIG_0(type), true) {
  }

  public CVTerm() : this(libsbmlPINVOKE.new_CVTerm__SWIG_1(), true) {
  }

  public CVTerm(XMLNode node) : this(libsbmlPINVOKE.new_CVTerm__SWIG_2(XMLNode.getCPtr(node)), true) {
    if (libsbmlPINVOKE.SWIGPendingException.Pending) throw libsbmlPINVOKE.SWIGPendingException.Retrieve();
  }

  public CVTerm(CVTerm orig) : this(libsbmlPINVOKE.new_CVTerm__SWIG_3(CVTerm.getCPtr(orig)), true) {
    if (libsbmlPINVOKE.SWIGPendingException.Pending) throw libsbmlPINVOKE.SWIGPendingException.Retrieve();
  }

  public CVTerm clone() {
    IntPtr cPtr = libsbmlPINVOKE.CVTerm_clone(swigCPtr);
    CVTerm ret = (cPtr == IntPtr.Zero) ? null : new CVTerm(cPtr, true);
    return ret;
  }

  public int getQualifierType() {
    int ret = libsbmlPINVOKE.CVTerm_getQualifierType(swigCPtr);
    return ret;
  }

  public int getModelQualifierType() {
    int ret = libsbmlPINVOKE.CVTerm_getModelQualifierType(swigCPtr);
    return ret;
  }

  public int getBiologicalQualifierType() {
    int ret = libsbmlPINVOKE.CVTerm_getBiologicalQualifierType(swigCPtr);
    return ret;
  }

  public XMLAttributes getResources() {
    IntPtr cPtr = libsbmlPINVOKE.CVTerm_getResources__SWIG_0(swigCPtr);
    XMLAttributes ret = (cPtr == IntPtr.Zero) ? null : new XMLAttributes(cPtr, false);
    return ret;
  }

  public long getNumResources() { return (long)libsbmlPINVOKE.CVTerm_getNumResources(swigCPtr); }

  public string getResourceURI(long n) {
    string ret = libsbmlPINVOKE.CVTerm_getResourceURI(swigCPtr, n);
    return ret;
  }

  public int setQualifierType(int type) {
    int ret = libsbmlPINVOKE.CVTerm_setQualifierType(swigCPtr, type);
    return ret;
  }

  public int setModelQualifierType(int type) {
    int ret = libsbmlPINVOKE.CVTerm_setModelQualifierType(swigCPtr, type);
    return ret;
  }

  public int setBiologicalQualifierType(int type) {
    int ret = libsbmlPINVOKE.CVTerm_setBiologicalQualifierType(swigCPtr, type);
    return ret;
  }

  public int addResource(string resource) {
    int ret = libsbmlPINVOKE.CVTerm_addResource(swigCPtr, resource);
    return ret;
  }

  public int removeResource(string resource) {
    int ret = libsbmlPINVOKE.CVTerm_removeResource(swigCPtr, resource);
    return ret;
  }

  public bool hasRequiredAttributes() {
    bool ret = libsbmlPINVOKE.CVTerm_hasRequiredAttributes(swigCPtr);
    return ret;
  }

}

}
