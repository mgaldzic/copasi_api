/* ----------------------------------------------------------------------------
 * This file was automatically generated by SWIG (http://www.swig.org).
 * Version 1.3.40
 *
 * Do not make changes to this file unless you know what you are doing--modify
 * the SWIG interface file instead.
 * ----------------------------------------------------------------------------- */

package org.sbml.libsbml;

/** 
 * LibSBML implementation of SBML's UnitDefinition construct.
 * <p>
 * Units of measurement may be supplied in a number of contexts in an SBML
 * model.  The SBML unit definition facility uses two classes of objects,
 * {@link UnitDefinition} and {@link Unit}.  The approach to defining units in SBML is
 * compositional; for example, <em>meter second<sup> &ndash;2</sup></em> is
 * constructed by combining a {@link Unit} object representing <em>meter</em> with
 * another {@link Unit} object representing <em>second<sup> &ndash;2</sup></em>.
 * The combination is wrapped inside a {@link UnitDefinition}, which provides for
 * assigning an identifier and optional name to the combination.  The
 * identifier can then be referenced from elsewhere in a model.  Thus, the
 * {@link UnitDefinition} class is the container, and {@link Unit} instances are placed
 * inside {@link UnitDefinition} instances.
 * <p>
 * Two points are worth discussing in the context of SBML units.  First,
 * unit declarations in SBML models are \<em>optional</em>.  The consequence of
 * this is that a model must be numerically self-consistent independently
 * of unit declarations, for the benefit of software tools that cannot
 * interpret or manipulate units.  {@link Unit} declarations in SBML are thus more
 * akin to a type of annotation; they can indicate intentions, and can be
 * used by model readers for checking the consistency of the model,
 * labeling simulation output, etc., but any transformations of values
 * implied by different units must be incorporated \<em>explicitly</em> into a
 * model.
 * <p>
 * Second, the vast majority of situations that require new SBML unit
 * definitions involve simple multiplicative combinations of base units and
 * factors.  An example is <em>moles per litre per second</em>.  What
 * distinguishes these sorts of unit definitions from more complex ones is
 * that they may be expressed without the use of an additive offset from a
 * zero point.  The use of offsets complicates all unit definition systems,
 * yet in the domain of SBML, the real-life cases requiring offsets are few
 * (and in fact, to the best of our knowledge, only involve temperature).
 * Consequently, the SBML unit system has been consciously designed to
 * simplify implementation of unit support for the most common cases in
 * systems biology.  The cost of this simplification is to require units
 * with offsets to be handled explicitly by the modeler.
 * <p>
 * <h2>Summary of the {@link UnitDefinition} construct</h2>
 * <p>
 * {@link UnitDefinition} has two attributes and one subelement.  The two
 * attributes are 'id' and 'name', and the subelement is {@link ListOfUnits}.
 * <p>
 * The required attribute 'id' and optional attribute 'name' are both
 * strings.  The 'id' attribute is used to give the defined unit a unique
 * identifier by which other parts of an SBML model definition can refer to
 * it.  The 'name' attribute is intended to be used for giving the unit
 * definition an optional human-readable name.  Please see the <a
 * href='#unitdef-id'>next section</a> for information about the values
 * permitted for 'id'.
 * <p>
 * A {@link UnitDefinition} must contain exactly one {@link ListOfUnits}, and this list
 * must contain one or more {@link Unit} definitions; see the definitions of these
 * other object classes for more information about them.  The following
 * example illustrates a complete unit definition (when written in XML)
 * when they all the pieces are combined together.  This defines 'mmls'
 * to be millimoles per litre per second.
 * <div class='fragment'><pre>
 * &lt;listOfUnitDefinitions&gt;
 *     &lt;unitDefinition id='mmls'&gt;
 *         &lt;listOfUnits&gt;
 *             &lt;unit kind='mole'   scale='-3'/&gt;
 *             &lt;unit kind='litre'  exponent='-1'/&gt;
 *             &lt;unit kind='second' exponent='-1'/&gt;
 *         &lt;/listOfUnits&gt;
 *     &lt;/unitDefinition&gt;
 * &lt;/listOfUnitDefinitions&gt;</pre></div>
 * <p>
 * <h2>Special considerations for {@link Unit} object identifiers</h2>
 * <p>
 * The attribute 'id' in {@link UnitDefinition} cannot be given simply any value,
 * and the precise details of the values permitted differ slightly between
 * Levels of SBML:
 * <ul>
 * <p>
 * <li> The 'id' of a {@link UnitDefinition} must <em>not</em> contain a value from the
 * list of SBML's predefined base unit names (i.e., the strings <code>gram</code>, 
 * <code>litre</code>, etc.).  In SBML Level&nbsp;3, this list consists of the
 * following:
 * <p>
 * <div>
<ul style='list-style-type: none'>
<li>ampere</li>
<li>avogadro</li>
<li>becquerel</li>
<li>candela</li>
<li>coulomb</li>
<li>dimensionless</li>
</ul>
</div>
<div style='margin-left: 10em; margin-top: -8em'>
<ul style='list-style-type: none'>
<li>farad</li>
<li>gram</li>
<li>gray</li>
<li>henry</li>
<li>hertz</li>
<li>item</li>
</ul>
</div>
<div style='margin-left: 19em; margin-top: -8em'>
<ul style='list-style-type: none'>
<li>joule</li>
<li>katal</li>
<li>kelvin</li>
<li>kilogram</li>
<li>litre</li>
<li>lumen</li>
</ul>
</div>
<div style='margin-left: 28em; margin-top: -8em'>
<ul style='list-style-type: none'>
<li>lux</li>
<li>metre</li>
<li>mole</li>
<li>newton</li>
<li>ohm</li>
<li>pascal</li>
</ul>
</div>
<div style='margin-left: 37em; margin-top: -8em'>
<ul style='list-style-type: none'>
<li>radian</li>
<li>second</li>
<li>siemens</li>
<li>sievert</li>
<li>steradian</li>
<li>tesla</li>
</ul>
</div>
<div style='margin-left: 46em; margin-top: -8em'>
<ul style='list-style-type: none'>
<li>volt</li>
<li>watt</li>
<li>weber</li>
</ul>
</div>
<br style='clear: both'>
<p style='padding-bottom: 1.5em'>

 * <p>
 * This list of predefined base units is nearly identical in SBML
 * Level&nbsp;2 Version&nbsp;4, the exception being that Level&nbsp;2 does
 * not define <code>avogadro</code>.  SBML Level&nbsp;2 Version&nbsp;1 (and only this
 * Level+Version combination) provides an additional predefined unit name,
 * <code>Celsius</code>.  SBML Level&nbsp;1 Versions&nbsp;2&ndash;3 provide two more
 * additional predefined unit names, <code>meter</code> and <code>liter</code>.
 * <p>
 * <li> In SBML Level&nbsp;2 (all Versions), there is an additional set of
 * reserved identifiers: <code>substance</code>, <code>volume</code>, <code>area</code>, <code>length</code>, and
 * <code>time</code>.  Using one of these values for the attribute 'id' of a
 * {@link UnitDefinition} has the effect of redefining the model-wide default units
 * for the corresponding quantities.  The list of special unit names in
 * SBML Level&nbsp;2 is given in the table below:
 * <p>
 *   <center>
<table border='0' class='text-table width80 normal-font alt-row-colors'>
 <tr>
     <th align='left'>Identifier</th>
     <th align='left'>Possible scalable units</th>
     <th align='left'>Default units</th>
 </tr>
<tr><td><code>substance</code></td><td>mole, item, gram, kilogram, dimensionless</td><td>mole</td></tr>
<tr><td><code>volume</code></td><td>litre, cubic metre, dimensionless</td><td>litre</td></tr>
<tr><td><code>area</code></td><td>square metre, dimensionless</td><td>square metre</td></tr>
<tr><td><code>length</code></td><td>metre, dimensionless</td><td>metre</td></tr>
<tr><td><code>time</code></td><td>second, dimensionless</td><td>second</td></tr>
</table>
</center>

 * <p>
 * Also, SBML Level&nbsp;2 imposes two limitations on redefining the
 * predefined unit <code>substance</code>, <code>volume</code>, <code>area</code>, <code>length</code>, and 
 * <code>time</code>: (1) The {@link UnitDefinition} of a predefined SBML unit can only contain
 * a single {@link Unit} object within it.  (2) The value of the 'kind' attribute
 * in a {@link Unit} instance must be drawn from one of the values in the second
 * column of the table above.
 * <p>
 * The special unit names <code>substance</code>, <code>volume</code>, <code>area</code>, <code>length</code>, and
 * <code>time</code> are not defined by SBML Level&nbsp;3, which uses a different
 * approach to setting model-wide inherited units.
 * <p>
 * </ul>
 * <p>
 * <p>
 * <h2>Further comments about SBML's unit definition system</h2>
 * <p>
 * The vast majority of modeling situations requiring new SBML unit
 * definitions involve simple multiplicative combinations of base units and
 * factors.  An example of this might be <em>moles per litre per
 * second</em>.  What distinguishes these sorts of simpler unit definitions
 * from more complex ones is that they may be expressed without the use of
 * an additive offset from a zero point.  The use of offsets complicates
 * all unit definition systems, yet in the domain of SBML the real-life
 * cases requiring offsets are few (and in fact, to the best of our
 * knowledge, only involve temperature).  Consequently, the SBML unit
 * system has been consciously designed in a way that attempts to simplify
 * implementation of unit support for the most common cases in systems
 * biology.
 * <p>
 * As of SBML Level&nbsp;2 Version&nbsp;2, {@link Unit} no longer has the
 * attribute called 'offset' introduced in SBML Level&nbsp;2
 * Version&nbsp;1.  It turned out that the general case involving units
 * with offsets was incorrectly defined, and few (if any) developers even
 * attempted to support offset-based units in their software.  In the
 * development of Level&nbsp;2 Version&nbsp;2, a consensus among SBML
 * developers emerged that a fully generalized unit scheme is <em>so</em>
 * confusing and complicated that it actually <em>impedes</em> interoperability.
 * SBML Level&nbsp;2 Version&nbsp;2, Version&nbsp;3 and Version&nbsp;4 acknowledge this
 * reality by reducing and simplifying the unit system, specifically by
 * removing the 'offset' attribute on {@link Unit} and <code>Celsius</code> as a pre-defined
 * unit.
 * <p>
 * The following guidelines suggest methods for handling units that do
 * require the use of zero offsets for their definitions:
 * <ul>
 * <li> <em>Handling Celsius</em>.  A model in which certain quantities are
 *   temperatures measured in degrees Celsius can be converted
 *   straightforwardly to a model in which those temperatures are in
 *   kelvin.  A software tool could do this by performing a straightforward
 *   substitution using the following relationship: T<sub> kelvin</sub> =
 *   T<sub> Celsius</sub> + 273.15.  In every mathematical formula of the
 *   model where a quantity (call it <em>x</em>) in degrees Celsius appears,
 *   replace <em>x</em> with x<sub> k</sub>+ 273.15, where x<sub> k</sub> is now
 *   in kelvin.  An alternative approach would be to use a
 *   {@link FunctionDefinition} to define a function encapsulating this
 *   relationship above and then using that in the rest of the model as
 *   needed.  Since Celsius is a commonly-used unit, software tools could
 *   help users by providing users with the ability to express temperatures
 *   in Celsius in the tools' interfaces, and making substitutions
 *   automatically when writing out the SBML.
 * <p>
 * <li> <em>Other units requiring offsets</em>.  One approach to handling
 *   other kinds of units is to use a {@link FunctionDefinition} to define a function
 *   encapsulating the necessary mathematical relationship, then
 *   substituting a call to this function wherever the original quantity
 *   appeared in the model.  For example, here is a possible definition for
 *   converting Fahrenheit to Celsius degrees:
 *   <div class='fragment'><pre>
 * &lt;functionDefinition id='Fahrenheit_to_kelvin'&gt;
 *     &lt;math xmlns='http://www.w3.org/1998/Math/MathML'&gt;
 *         &lt;lambda&gt;
 *             &lt;bvar&gt;&lt;ci&gt; temp_in_fahrenheit &lt;/ci&gt;&lt;/bvar&gt;
 *             &lt;apply&gt;
 *                 &lt;divide/&gt;
 *                 &lt;apply&gt;
 *                     &lt;plus/&gt;
 *                     &lt;ci&gt; temp_in_fahrenheit &lt;/ci&gt;
 *                     &lt;cn&gt; 459.67 &lt;/cn&gt;
 *                 &lt;/apply&gt;
 *                 &lt;cn&gt; 1.8 &lt;/cn&gt;
 *             &lt;/apply&gt;
 *         &lt;/lambda&gt;
 *     &lt;/math&gt;
 * &lt;/functionDefinition&gt;</pre></div>
 * <p>
 * <li> An alternative approach not requiring the use of function definitions
 *   is to use an {@link AssignmentRule} for each variable in Fahrenheit units.
 *   The {@link AssignmentRule} could compute the conversion from Fahrenheit to
 *   (say) kelvin, assign its value to a variable (in Kelvin units), and
 *   then that variable could be used elsewhere in the model.
 * <p>
 * <li> Still another approach is to rewrite the mathematical formulas of a
 *   model to directly incorporate the conversion formula wherever the
 *   original quantity appeared.
 * </ul>
 * <p>
 * Please consult the SBML specifications for more information about this
 * and other issues involving units.
 * <p>
 * <!-- leave this next break as-is to work around some doxygen bug -->
 */

public class UnitDefinition extends SBase {
   private long swigCPtr;

   protected UnitDefinition(long cPtr, boolean cMemoryOwn)
   {
     super(libsbmlJNI.SWIGUnitDefinitionUpcast(cPtr), cMemoryOwn);
     swigCPtr = cPtr;
   }

   protected static long getCPtr(UnitDefinition obj)
   {
     return (obj == null) ? 0 : obj.swigCPtr;
   }

   protected static long getCPtrAndDisown (UnitDefinition obj)
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
        libsbmlJNI.delete_UnitDefinition(swigCPtr);
      }
      swigCPtr = 0;
    }
    super.delete();
  }

  
  /**
   * Creates a new {@link UnitDefinition} using the given SBML <code>level</code> and <code>version</code>
   * values.
   * <p>
   * @param level a long integer, the SBML Level to assign to this {@link UnitDefinition}
   * <p>
   * @param version a long integer, the SBML Version to assign to this
   * {@link UnitDefinition}
   * <p>
   * @note Upon the addition of a {@link UnitDefinition} object to an {@link SBMLDocument}
   * (e.g., using Model.addUnitDefinition()), the SBML Level, SBML Version
   * and XML namespace of the document <em>override</em> the values used
   * when creating the {@link UnitDefinition} object via this constructor.  This is
   * necessary to ensure that an SBML document is a consistent structure.
   * Nevertheless, the ability to supply the values at the time of creation
   * of a {@link UnitDefinition} is an important aid to producing valid SBML.
   * Knowledge of the intented SBML Level and Version determine whether it
   * is valid to assign a particular value to an attribute, or whether it
   * is valid to add an object to an existing {@link SBMLDocument}.
   */
 public UnitDefinition(long level, long version) throws org.sbml.libsbml.SBMLConstructorException {
    this(libsbmlJNI.new_UnitDefinition__SWIG_0(level, version), true);
  }

  
  /**
   * Creates a new {@link UnitDefinition} using the given {@link SBMLNamespaces} object
   * <code>sbmlns</code>.
   * <p>
   * The {@link SBMLNamespaces} object encapsulates SBML Level/Version/namespaces
   * information.  It is used to communicate the SBML Level, Version, and
   * (in Level&nbsp;3) packages used in addition to SBML Level&nbsp;3 Core.
   * A common approach to using this class constructor is to create an
   * {@link SBMLNamespaces} object somewhere in a program, once, then pass it to
   * object constructors such as this one when needed.
   * <p>
   * @param sbmlns an {@link SBMLNamespaces} object.
   * <p>
   * @note Upon the addition of a {@link UnitDefinition} object to an {@link SBMLDocument}
   * (e.g., using Model.addUnitDefinition()), the SBML XML namespace of
   * the document <em>overrides</em> the value used when creating the
   * {@link UnitDefinition} object via this constructor.  This is necessary to
   * ensure that an SBML document is a consistent structure.  Nevertheless,
   * the ability to supply the values at the time of creation of a
   * {@link UnitDefinition} is an important aid to producing valid SBML.  Knowledge
   * of the intented SBML Level and Version determine whether it is valid
   * to assign a particular value to an attribute, or whether it is valid
   * to add an object to an existing {@link SBMLDocument}.
   */
 public UnitDefinition(SBMLNamespaces sbmlns) throws org.sbml.libsbml.SBMLConstructorException {
    this(libsbmlJNI.new_UnitDefinition__SWIG_1(SBMLNamespaces.getCPtr(sbmlns), sbmlns), true);
  }

  
  /**
  * Copy constructor; creates a copy of this {@link UnitDefinition}.
  */
 public UnitDefinition(UnitDefinition orig) throws org.sbml.libsbml.SBMLConstructorException {
    this(libsbmlJNI.new_UnitDefinition__SWIG_2(UnitDefinition.getCPtr(orig), orig), true);
  }

  
  /**
   * Creates and returns a deep copy of this {@link UnitDefinition}.
   * <p>
   * @return a (deep) copy of this {@link UnitDefinition}.
   */
 public UnitDefinition cloneObject() {
    long cPtr = libsbmlJNI.UnitDefinition_cloneObject(swigCPtr, this);
    return (cPtr == 0) ? null : new UnitDefinition(cPtr, true);
  }

  
  /**
   * Returns the value of the 'id' attribute of this {@link UnitDefinition}.
   * <p>
   * @return the id of this {@link UnitDefinition}.
   */
 public String getId() {
    return libsbmlJNI.UnitDefinition_getId(swigCPtr, this);
  }

  
  /**
   * Returns the value of the 'name' attribute of this {@link UnitDefinition}.
   * <p>
   * @return the name of this {@link UnitDefinition}.
   */
 public String getName() {
    return libsbmlJNI.UnitDefinition_getName(swigCPtr, this);
  }

  
  /**
   * Predicate returning <code>true</code> if this
   * {@link UnitDefinition}'s 'id' attribute has been set.
   * <p>
   * @return <code>true</code> if the 'id' attribute of this {@link UnitDefinition} has been
   * set, <code>false</code> otherwise.
   */
 public boolean isSetId() {
    return libsbmlJNI.UnitDefinition_isSetId(swigCPtr, this);
  }

  
  /**
   * Predicate returning <code>true</code> if this
   * {@link UnitDefinition}'s 'name' attribute has been set.
   * <p>
   * @return <code>true</code> if the 'name' attribute of this {@link UnitDefinition} has been
   * set, <code>false</code> otherwise.
   */
 public boolean isSetName() {
    return libsbmlJNI.UnitDefinition_isSetName(swigCPtr, this);
  }

  
  /**
   * Sets the value of the 'id' attribute of this {@link UnitDefinition}.
   * <p>
   * The string <code>sid</code> is copied.  Note that SBML has strict requirements
   * for the syntax of identifiers.  The following is a summary of the definition of the SBML identifier type 
<code>SId</code>, which defines the permitted syntax of identifiers.  We
express the syntax using an extended form of BNF notation: 
<p>
<code style='margin-left: 2em'>letter .= 'a'..'z','A'..'Z'</code><br>
<code style='margin-left: 2em'>digit  .= '0'..'9'</code><br>
<code style='margin-left: 2em'>idChar .= letter | digit | '_'</code><br>
<code style='margin-left: 2em'>SId    .= ( letter | '_' ) idChar*</code><br>
<p>
The characters <code>(</code> and <code>)</code> are used for grouping, the
character <code>*</code> 'zero or more times', and the character
<code>|</code> indicates logical 'or'.  The equality of SBML identifiers is
determined by an exact character sequence match; i.e., comparisons must be
performed in a case-sensitive manner.  In addition, there are a few
conditions for the uniqueness of identifiers in an SBML model.  Please
consult the SBML specifications for the exact formulations.
<p>

   * <p>
   * @param sid the string to use as the identifier of this {@link UnitDefinition}
   * <p>
   * @return integer value indicating success/failure of the
   * function.   The possible values
   * returned by this function are:
   * <li> {@link  libsbmlConstants#LIBSBML_OPERATION_SUCCESS LIBSBML_OPERATION_SUCCESS }
   * <li> {@link  libsbmlConstants#LIBSBML_INVALID_ATTRIBUTE_VALUE LIBSBML_INVALID_ATTRIBUTE_VALUE }
   */
 public int setId(String sid) {
    return libsbmlJNI.UnitDefinition_setId(swigCPtr, this, sid);
  }

  
  /**
   * Sets the value of the 'name' attribute of this {@link UnitDefinition}.
   * <p>
   * The string in <code>name</code> is copied.
   * <p>
   * @param name the new name for the {@link UnitDefinition}
   * <p>
   * @return integer value indicating success/failure of the
   * function.   The possible values
   * returned by this function are:
   * <li> {@link  libsbmlConstants#LIBSBML_OPERATION_SUCCESS LIBSBML_OPERATION_SUCCESS }
   * <li> {@link  libsbmlConstants#LIBSBML_INVALID_ATTRIBUTE_VALUE LIBSBML_INVALID_ATTRIBUTE_VALUE }
   */
 public int setName(String name) {
    return libsbmlJNI.UnitDefinition_setName(swigCPtr, this, name);
  }

  
  /**
   * Unsets the value of the 'name' attribute of this {@link UnitDefinition}.
   * <p>
   * @return integer value indicating success/failure of the
   * function.   The possible values
   * returned by this function are:
   * <li> {@link  libsbmlConstants#LIBSBML_OPERATION_SUCCESS LIBSBML_OPERATION_SUCCESS }
   * <li> {@link  libsbmlConstants#LIBSBML_OPERATION_FAILED LIBSBML_OPERATION_FAILED }
   */
 public int unsetName() {
    return libsbmlJNI.UnitDefinition_unsetName(swigCPtr, this);
  }

  
  /**
   * Convenience function for testing if a given unit definition is a
   * variant of the predefined unit identifier <code>'area'</code>.
   * <p>
   * @return <code>true</code> if this {@link UnitDefinition} is a variant of the predefined
   * unit <code>area</code>, meaning square metres with only abritrary variations
   * in scale or multiplier values; <code>false</code> otherwise.
   */
 public boolean isVariantOfArea() {
    return libsbmlJNI.UnitDefinition_isVariantOfArea(swigCPtr, this);
  }

  
  /**
   * Convenience function for testing if a given unit definition is a
   * variant of the predefined unit identifier <code>'length'</code>.
   * <p>
   * @return <code>true</code> if this {@link UnitDefinition} is a variant of the predefined
   * unit <code>length</code>, meaning metres with only abritrary variations in scale
   * or multiplier values; <code>false</code> otherwise.
   */
 public boolean isVariantOfLength() {
    return libsbmlJNI.UnitDefinition_isVariantOfLength(swigCPtr, this);
  }

  
  /**
   * Convenience function for testing if a given unit definition is a
   * variant of the predefined unit identifier <code>'substance'</code>.
   * <p>
   * @return <code>true</code> if this {@link UnitDefinition} is a variant of the predefined
   * unit <code>substance</code>, meaning moles or items (and grams or kilograms from
   * SBML Level&nbsp;2 Version&nbsp;2 onwards) with only abritrary variations
   * in scale or multiplier values; <code>false</code> otherwise.
   */
 public boolean isVariantOfSubstance() {
    return libsbmlJNI.UnitDefinition_isVariantOfSubstance(swigCPtr, this);
  }

  
  /**
   * Convenience function for testing if a given unit definition is a
   * variant of the predefined unit identifier <code>'time'</code>.
   * <p>
   * @return <code>true</code> if this {@link UnitDefinition} is a variant of the predefined
   * unit <code>time</code>, meaning seconds with only abritrary variations in scale or
   * multiplier values; <code>false</code> otherwise.
   */
 public boolean isVariantOfTime() {
    return libsbmlJNI.UnitDefinition_isVariantOfTime(swigCPtr, this);
  }

  
  /**
   * Convenience function for testing if a given unit definition is a
   * variant of the predefined unit identifier <code>'volume'</code>.
   * <p>
   * @return <code>true</code> if this {@link UnitDefinition} is a variant of the predefined
   * unit <code>volume</code>, meaning litre or cubic metre with only abritrary
   * variations in scale or multiplier values; <code>false</code> otherwise.
   */
 public boolean isVariantOfVolume() {
    return libsbmlJNI.UnitDefinition_isVariantOfVolume(swigCPtr, this);
  }

  
  /**
   * Convenience function for testing if a given unit definition is a
   * variant of the unit <code>'dimensionless'</code>.
   * <p>
   * @return <code>true</code> if this {@link UnitDefinition} is a variant of 
   * <code>dimensionless</code>, meaning dimensionless with only abritrary variations in
   * scale or multiplier values; <code>false</code> otherwise.
   */
 public boolean isVariantOfDimensionless() {
    return libsbmlJNI.UnitDefinition_isVariantOfDimensionless(swigCPtr, this);
  }

  
  /**
   * Convenience function for testing if a given unit definition is a
   * variant of the predefined unit identifier <code>'mass'</code>.
   * <p>
   * @return <code>true</code> if this {@link UnitDefinition} is a variant of mass units,
   * meaning gram or kilogram with only abritrary variations in scale or
   * multiplier values; <code>false</code> otherwise.
   */
 public boolean isVariantOfMass() {
    return libsbmlJNI.UnitDefinition_isVariantOfMass(swigCPtr, this);
  }

  
  /**
   * Convenience function for testing if a given unit definition is a
   * variant of the predefined unit <code>'substance'</code> divided by the predefined
   * unit <code>'time'</code>.
   * <p>
   * @return <code>true</code> if this {@link UnitDefinition} is a variant of the predefined
   * unit <code>substance</code> per predefined unit <code>time</code>, meaning it contains two
   * units one of which is a variant of substance and the other is a
   * variant of time which an exponent of -1; <code>false</code> otherwise.
   */
 public boolean isVariantOfSubstancePerTime() {
    return libsbmlJNI.UnitDefinition_isVariantOfSubstancePerTime(swigCPtr, this);
  }

  
  /**
   * Adds a copy of the given {@link Unit} to this {@link UnitDefinition}.
   * <p>
   * @param u the {@link Unit} instance to add to this {@link UnitDefinition}.
   * <p>
   * @return integer value indicating success/failure of the
   * function.   The possible values
   * returned by this function are:
   * <li> {@link  libsbmlConstants#LIBSBML_OPERATION_SUCCESS LIBSBML_OPERATION_SUCCESS }
   * <li> {@link  libsbmlConstants#LIBSBML_LEVEL_MISMATCH LIBSBML_LEVEL_MISMATCH }
   * <li> {@link  libsbmlConstants#LIBSBML_VERSION_MISMATCH LIBSBML_VERSION_MISMATCH }
   * <li> {@link  libsbmlConstants#LIBSBML_DUPLICATE_OBJECT_ID LIBSBML_DUPLICATE_OBJECT_ID }
   * <li> {@link  libsbmlConstants#LIBSBML_OPERATION_FAILED LIBSBML_OPERATION_FAILED }
   * <p>
   * @note This method should be used with some caution.  The fact that
   * this method <em>copies</em> the object passed to it means that the caller
   * will be left holding a physically different object instance than the
   * one contained in this {@link UnitDefinition}.  Changes made to the original
   * object instance (such as resetting attribute values) will <em>not
   * affect the instance in the {@link UnitDefinition}</em>.  In addition, the
   * caller should make sure to free the original object if it is no longer
   * being used, or else a memory leak will result.  Please see
   * UnitDefinition.createUnit() for a method that does not lead to these
   * issues.
   * <p>
   * @see #createUnit()
   */
 public int addUnit(Unit u) {
    return libsbmlJNI.UnitDefinition_addUnit(swigCPtr, this, Unit.getCPtr(u), u);
  }

  
  /**
   * Creates a new and empty {@link Unit}, adds it to this {@link UnitDefinition}'s list of
   * units, and returns it.
   * <p>
   * @return a newly constructed (and empty) {@link Unit} instance.
   * <p>
   * @note It is worth emphasizing that the attribute 'kind' value of a
   * {@link Unit} is a required attribute for a valid {@link Unit} definition.  The
   * createUnit() method does not assign a valid kind to the constructed
   * unit (instead, it sets the 'kind' to <code>UNIT_KIND_INVALID</code>).  Callers
   * are cautioned to set the newly-constructed {@link Unit}'s kind using
   * Unit.setKind() soon after calling this method.
   * <p>
   * @see #addUnit(Unit  u)
   */
 public Unit createUnit() {
    long cPtr = libsbmlJNI.UnitDefinition_createUnit(swigCPtr, this);
    return (cPtr == 0) ? null : new Unit(cPtr, false);
  }

  
  /**
   * Returns the list of Units for this {@link UnitDefinition} instance.
   * @return the {@link ListOfUnits} value for this {@link UnitDefinition}.
   */
 public ListOfUnits getListOfUnits() {
    long cPtr = libsbmlJNI.UnitDefinition_getListOfUnits__SWIG_0(swigCPtr, this);
    return (cPtr == 0) ? null : new ListOfUnits(cPtr, false);
  }

  
  /**
   * Returns a specific {@link Unit} instance belonging to this {@link UnitDefinition}.
   * <p>
   * @param n an integer, the index of the {@link Unit} to be returned.
   * <p>
   * @return the nth {@link Unit} of this {@link UnitDefinition}.
   * <p>
   * @see #getNumUnits()
   */
 public Unit getUnit(long n) {
    long cPtr = libsbmlJNI.UnitDefinition_getUnit__SWIG_0(swigCPtr, this, n);
    return (cPtr == 0) ? null : new Unit(cPtr, false);
  }

  
  /**
   * Returns the number of {@link Unit} objects contained within this
   * {@link UnitDefinition}.
   * <p>
   * @return an integer representing the number of Units in this
   * {@link UnitDefinition}.
   */
 public long getNumUnits() {
    return libsbmlJNI.UnitDefinition_getNumUnits(swigCPtr, this);
  }

  
  /**
   * Removes the nth {@link Unit} object from this {@link UnitDefinition} object and
   * returns a pointer to it.
   * <p>
   * The caller owns the returned object and is responsible for deleting it.
   * <p>
   * @param n the index of the {@link Unit} object to remove
   * <p>
   * @return the {@link Unit} object removed, or <code>NULL</code> if the given index 
   * is out of range.
   * <p>
   */
 public Unit removeUnit(long n) {
    long cPtr = libsbmlJNI.UnitDefinition_removeUnit(swigCPtr, this, n);
    return (cPtr == 0) ? null : new Unit(cPtr, true);
  }

  
  /**
   * Returns the libSBML type code for this object instance.
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
    return libsbmlJNI.UnitDefinition_getTypeCode(swigCPtr, this);
  }

  
  /**
   * Returns the XML element name of this object, which for {@link UnitDefinition},
   * is always <code>'unitDefinition'</code>.
   * <p>
   * @return the name of this element, i.e., <code>'unitDefinition'</code>.
   */
 public String getElementName() {
    return libsbmlJNI.UnitDefinition_getElementName(swigCPtr, this);
  }

  
  /** 
   * Simplifies the {@link UnitDefinition} such that any given kind of {@link Unit} object
   * occurs only once in the {@link ListOfUnits}.
   * <p>
   * For example, the following definition,
   * <div class='fragment'><pre>
   * &lt;unitDefinition&gt;
   *  &lt;listOfUnits&gt;
   *    &lt;unit kind='metre' exponent='1'/&gt;
   *    &lt;unit kind='metre' exponent='2'/&gt;
   *  &lt;/listOfUnits&gt;
   * &lt;unitDefinition&gt;</pre></div>
   * will be simplified to 
   * <div class='fragment'><pre>
   * &lt;unitDefinition&gt;
   *   &lt;listOfUnits&gt;
   *     &lt;unit kind='metre' exponent='3'/&gt;
   *   &lt;/listOfUnits&gt;
   * &lt;unitDefinition&gt;</pre></div>
   * <p>
   * @param ud the {@link UnitDefinition} object to be simplified.
   */
 public static void simplify(UnitDefinition ud) {
    libsbmlJNI.UnitDefinition_simplify(UnitDefinition.getCPtr(ud), ud);
  }

  
  /** 
   * Alphabetically orders the {@link Unit} objects within the {@link ListOfUnits} of a
   * {@link UnitDefinition}.
   * <p>
   * @param ud the {@link UnitDefinition} object whose units are to be reordered.
   */
 public static void reorder(UnitDefinition ud) {
    libsbmlJNI.UnitDefinition_reorder(UnitDefinition.getCPtr(ud), ud);
  }

  
  /**
   * Convert a given {@link UnitDefinition} into a new {@link UnitDefinition} object
   * that uses SI units.
   * <p>
   * @param ud the {@link UnitDefinition} object to convert to SI
   * <p>
   * @return a new {@link UnitDefinition} object representing the results of the
   * conversion.
   */
 public static UnitDefinition convertToSI(UnitDefinition arg0) {
    long cPtr = libsbmlJNI.UnitDefinition_convertToSI(UnitDefinition.getCPtr(arg0), arg0);
    return (cPtr == 0) ? null : new UnitDefinition(cPtr, true);
  }

  
  /** 
   * Predicate returning <code>true</code> if two
   * {@link UnitDefinition} objects are identical.
   * <p>
   * For the purposes of performing this comparison, two {@link UnitDefinition}
   * objects are considered identical when they contain identical lists of
   * {@link Unit} objects.  Pairs of {@link Unit} objects in the lists are in turn
   * considered identical if they satisfy the predicate
   * Unit.areIdentical().  The predicate compares every attribute of the
   * {@link Unit} objects.
   * <p>
   * @param ud1 the first {@link UnitDefinition} object to compare
   * @param ud2 the second {@link UnitDefinition} object to compare
   * <p>
   * @return <code>true</code> if all the {@link Unit} objects in ud1 are identical to the
   * {@link Unit} objects of ud2, <code>false</code> otherwise.
   * <p>
   * @see #areEquivalent(UnitDefinition  ud1, UnitDefinition  ud2)
   * @see Unit#areIdentical(Unit  unit1, Unit  unit2)
   */
 public static boolean areIdentical(UnitDefinition ud1, UnitDefinition ud2) {
    return libsbmlJNI.UnitDefinition_areIdentical(UnitDefinition.getCPtr(ud1), ud1, UnitDefinition.getCPtr(ud2), ud2);
  }

  
  /** 
   * Predicate returning <code>true</code> if two
   * {@link UnitDefinition} objects are equivalent.
   * <p>
   * For the purposes of performing this comparison, two {@link UnitDefinition}
   * objects are considered equivalent when they contain <em>equivalent</em>
   * list of {@link Unit} objects.  {@link Unit} objects are in turn considered equivalent
   * if they satisfy the predicate Unit.areEquivalent().  The predicate
   * tests a subset of the objects's attributes.
   * <p>
   * @param ud1 the first {@link UnitDefinition} object to compare
   * <p>
   * @param ud2 the second {@link UnitDefinition} object to compare
   * <p>
   * @return <code>true</code> if all the {@link Unit} objects in ud1 are equivalent
   * to the {@link Unit} objects in ud2, <code>false</code> otherwise.
   * <p>
   * @see #areIdentical(UnitDefinition  ud1, UnitDefinition  ud2)
   * @see Unit#areEquivalent(Unit  unit1, Unit  unit2)
   */
 public static boolean areEquivalent(UnitDefinition ud1, UnitDefinition ud2) {
    return libsbmlJNI.UnitDefinition_areEquivalent(UnitDefinition.getCPtr(ud1), ud1, UnitDefinition.getCPtr(ud2), ud2);
  }

  
  /** 
   * Combines two {@link UnitDefinition} objects into a single {@link UnitDefinition}.
   * <p>
   * This takes {@link UnitDefinition} objects <code>ud1</code> and <code>ud2</code>, and creates a
   * {@link UnitDefinition} object that expresses the product of the units of 
   * <code>ud1</code> and <code>ud2</code>.
   * <p>
   * @param ud1 the first {@link UnitDefinition} object 
   * @param ud2 the second {@link UnitDefinition} object
   * <p>
   * @return a {@link UnitDefinition} which represents the product of the 
   * units of the two argument UnitDefinitions.
   */
 public static UnitDefinition combine(UnitDefinition ud1, UnitDefinition ud2) {
    long cPtr = libsbmlJNI.UnitDefinition_combine(UnitDefinition.getCPtr(ud1), ud1, UnitDefinition.getCPtr(ud2), ud2);
    return (cPtr == 0) ? null : new UnitDefinition(cPtr, true);
  }

  
  /** 
   * Expresses the given definition in a plain-text form.
   * <p>
   * For example, printUnits() applied to
   * <div class='fragment'><pre>
   * &lt;unitDefinition&gt;
   *  &lt;listOfUnits&gt;
   *    &lt;unit kind='metre' exponent='1'/&gt;
   *    &lt;unit kind='second' exponent='-2'/&gt;
   *  &lt;/listOfUnits&gt;
   * &lt;unitDefinition&gt;</pre></div>
   * will return the string <code>'metre (exponent = 1, multiplier = 1,
   * scale = 0) second (exponent = -2, multiplier = 1, scale = 0)'</code>
   * or, if the optional parameter <code>compact</code> is given the value <code>true</code>,
   * the string <code>'(1 metre)^1 (1 second)^-2'</code>.  This method may
   * be useful for printing unit information to human users, or in
   * debugging software, or other situations.
   * <p>
   * @param ud the {@link UnitDefinition} object
   * @param compact boolean indicating whether the compact form
   * should be used (defaults to false)
   * <p>
   * @return a string expressing the unit definition defined by the given
   * {@link UnitDefinition} object <code>ud</code>.
   */
 public static String printUnits(UnitDefinition ud, boolean compact) {
    return libsbmlJNI.UnitDefinition_printUnits__SWIG_0(UnitDefinition.getCPtr(ud), ud, compact);
  }

  
  /** 
   * Expresses the given definition in a plain-text form.
   * <p>
   * For example, printUnits() applied to
   * <div class='fragment'><pre>
   * &lt;unitDefinition&gt;
   *  &lt;listOfUnits&gt;
   *    &lt;unit kind='metre' exponent='1'/&gt;
   *    &lt;unit kind='second' exponent='-2'/&gt;
   *  &lt;/listOfUnits&gt;
   * &lt;unitDefinition&gt;</pre></div>
   * will return the string <code>'metre (exponent = 1, multiplier = 1,
   * scale = 0) second (exponent = -2, multiplier = 1, scale = 0)'</code>
   * or, if the optional parameter <code>compact</code> is given the value <code>true</code>,
   * the string <code>'(1 metre)^1 (1 second)^-2'</code>.  This method may
   * be useful for printing unit information to human users, or in
   * debugging software, or other situations.
   * <p>
   * @param ud the {@link UnitDefinition} object
   * @param compact boolean indicating whether the compact form
   * should be used (defaults to false)
   * <p>
   * @return a string expressing the unit definition defined by the given
   * {@link UnitDefinition} object <code>ud</code>.
   */
 public static String printUnits(UnitDefinition ud) {
    return libsbmlJNI.UnitDefinition_printUnits__SWIG_1(UnitDefinition.getCPtr(ud), ud);
  }

  
  /**
   * Predicate returning <code>true</code> if
   * all the required attributes for this {@link UnitDefinition} object
   * have been set.
   * <p>
   * @note The required attributes for a {@link UnitDefinition} object are:
   * <li> 'id'
   * <p>
   * @return a boolean value indicating whether all the required
   * attributes for this object have been defined.
   */
 public boolean hasRequiredAttributes() {
    return libsbmlJNI.UnitDefinition_hasRequiredAttributes(swigCPtr, this);
  }

  
  /**
   * Predicate returning <code>true</code> if
   * all the required elements for this {@link UnitDefinition} object
   * have been set.
   * <p>
   * @note The required elements for a {@link Constraint} object are:
   * <li> 'listOfUnits' (required in SBML Level&nbsp;2 only, optional in Level&nbsp;3)
   * <p>
   * @return a boolean value indicating whether all the required
   * elements for this object have been defined.
   */
 public boolean hasRequiredElements() {
    return libsbmlJNI.UnitDefinition_hasRequiredElements(swigCPtr, this);
  }

}
