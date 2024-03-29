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

public class ListOfReactions : ListOf {
	private HandleRef swigCPtr;
	
	internal ListOfReactions(IntPtr cPtr, bool cMemoryOwn) : base(libsbmlPINVOKE.ListOfReactionsUpcast(cPtr), cMemoryOwn)
	{
		//super(libsbmlPINVOKE.ListOfReactionsUpcast(cPtr), cMemoryOwn);
		swigCPtr = new HandleRef(this, cPtr);
	}
	
	internal static HandleRef getCPtr(ListOfReactions obj)
	{
		return (obj == null) ? new HandleRef(null, IntPtr.Zero) : obj.swigCPtr;
	}
	
	internal static HandleRef getCPtrAndDisown (ListOfReactions obj)
	{
		HandleRef ptr = new HandleRef(null, IntPtr.Zero);
		
		if (obj != null)
		{
			ptr             = obj.swigCPtr;
			obj.swigCMemOwn = false;
		}
		
		return ptr;
	}

  ~ListOfReactions() {
    Dispose();
  }

  public override void Dispose() {
    lock(this) {
      if (swigCPtr.Handle != IntPtr.Zero) {
        if (swigCMemOwn) {
          swigCMemOwn = false;
          libsbmlPINVOKE.delete_ListOfReactions(swigCPtr);
        }
        swigCPtr = new HandleRef(null, IntPtr.Zero);
      }
      GC.SuppressFinalize(this);
      base.Dispose();
    }
  }

  public new ListOfReactions clone() {
    IntPtr cPtr = libsbmlPINVOKE.ListOfReactions_clone(swigCPtr);
    ListOfReactions ret = (cPtr == IntPtr.Zero) ? null : new ListOfReactions(cPtr, true);
    return ret;
  }

  public override int getTypeCode() {
    int ret = libsbmlPINVOKE.ListOfReactions_getTypeCode(swigCPtr);
    return ret;
  }

  public override int getItemTypeCode() {
    int ret = libsbmlPINVOKE.ListOfReactions_getItemTypeCode(swigCPtr);
    return ret;
  }

  public override string getElementName() {
    string ret = libsbmlPINVOKE.ListOfReactions_getElementName(swigCPtr);
    return ret;
  }

  public new Reaction get(long n) {
    IntPtr cPtr = libsbmlPINVOKE.ListOfReactions_get__SWIG_0(swigCPtr, n);
    Reaction ret = (cPtr == IntPtr.Zero) ? null : new Reaction(cPtr, false);
    return ret;
  }

  public virtual Reaction get(string sid) {
    IntPtr cPtr = libsbmlPINVOKE.ListOfReactions_get__SWIG_2(swigCPtr, sid);
    Reaction ret = (cPtr == IntPtr.Zero) ? null : new Reaction(cPtr, false);
    return ret;
  }

  public new Reaction remove(long n) {
    IntPtr cPtr = libsbmlPINVOKE.ListOfReactions_remove__SWIG_0(swigCPtr, n);
    Reaction ret = (cPtr == IntPtr.Zero) ? null : new Reaction(cPtr, true);
    return ret;
  }

  public virtual Reaction remove(string sid) {
    IntPtr cPtr = libsbmlPINVOKE.ListOfReactions_remove__SWIG_1(swigCPtr, sid);
    Reaction ret = (cPtr == IntPtr.Zero) ? null : new Reaction(cPtr, true);
    return ret;
  }

  public ListOfReactions() : this(libsbmlPINVOKE.new_ListOfReactions(), true) {
  }

}

}
