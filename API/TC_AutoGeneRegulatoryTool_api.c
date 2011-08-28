#include "TC_AutoGeneRegulatoryTool_api.h"
#include "TC_BasicInformationTool_api.h"

tc_items (*_tc_partsIn)(long) = 0;
/*! 
 \brief Get all DNA parts inside the given container or module
 \ingroup Get items
*/ TCAPIEXPORT 
tc_items tc_partsIn(long o)
{
	if (_tc_partsIn)
		return _tc_partsIn(o);
	return tc_createItemsArray(0);
}

tc_items (*_tc_partsUpstream)(long) = 0;
/*! 
 \brief Get all DNA parts upstream of the given part
 \ingroup Get items
*/ TCAPIEXPORT 
tc_items tc_partsUpstream(long o)
{
	if (_tc_partsUpstream)
		return _tc_partsUpstream(o);
	return tc_createItemsArray(0);
}

tc_items (*_tc_partsDownstream)(long) = 0;
/*! 
 \brief Get all DNA parts downstream of the given part
 \ingroup Get items
*/ TCAPIEXPORT 
tc_items tc_partsDownstream(long o)
{
	if (_tc_partsDownstream)
		return _tc_partsDownstream(o);
	return tc_createItemsArray(0);
}

void (*_tc_alignParts)(tc_items) = 0;
/*! 
 \brief Align the given DNA parts in the order given
 \ingroup Get and set position
*/ TCAPIEXPORT 
void tc_alignParts(tc_items a)
{
	if (_tc_alignParts)
		_tc_alignParts(a);
}

void (*_tc_alignPartsOnPlasmid)(long,tc_items) = 0;
/*! 
 \brief Align the given DNA parts in the order given
 \ingroup Get and set position
*/ TCAPIEXPORT 
void tc_alignPartsOnPlasmid(long o, tc_items a)
{
	if (_tc_alignPartsOnPlasmid)
		_tc_alignPartsOnPlasmid(o, a);
}

/*! 
 \brief Assign DNA sequence to a part
 \ingroup Get and set position
*/ TCAPIEXPORT 
void tc_setSequence(long o, const char * s)
{
	tc_setTextAttribute(o,"sequence",s);
}

/*! 
 \brief initialize grouping
 \ingroup init
*/ TCAPIEXPORT 
void tc_AutoGeneRegulatoryTool_api(
		tc_items (*f1)(long), 
		tc_items (*f2)(long), 
		tc_items (*f3)(long), 
		void (*f4)(tc_items),
		void (*f5)(long,tc_items)
	)
{
	_tc_partsIn = f1;
	_tc_partsUpstream = f2;
	_tc_partsDownstream = f3;
	_tc_alignParts = f4;
	_tc_alignPartsOnPlasmid = f5;
}

