// Begin CVS Header
//   $Source: /fs/turing/cvs/copasi_dev/copasi/UI/parametertable.cpp,v $
//   $Revision: 1.30 $
//   $Name: Build-33 $
//   $Author: ssahle $
//   $Date: 2009/05/08 22:38:18 $
// End CVS Header

// Copyright (C) 2008 by Pedro Mendes, Virginia Tech Intellectual
// Properties, Inc., EML Research, gGmbH, University of Heidelberg,
// and The University of Manchester.
// All rights reserved.

// Copyright (C) 2001 - 2007 by Pedro Mendes, Virginia Tech Intellectual
// Properties, Inc. and EML Research, gGmbH.
// All rights reserved.

#include <qstringlist.h>
#include <q3header.h>
#include <qcombobox.h>
#include <qlineedit.h>
//Added by qt3to4:
#include <QPixmap>
#include <stdlib.h>

#include "parametertable.h"
#include "model/CReactionInterface.h"
#include "model/CModel.h"
#include "model/CMetabNameInterface.h"
#include "qtUtilities.h"
#include "utilities/CDimension.h"
#include "copasi/report/CCopasiRootContainer.h"

#include "./icons/product.xpm"
#include "./icons/substrate.xpm"
#include "./icons/modifier.xpm"
#include "./icons/locked.xpm"
#include "./icons/unlocked.xpm"

ParameterTable::ParameterTable(QWidget * parent, const char * name)
    : Q3Table(parent, name),
    mOldRow(0)
{
  initTable();

  connect(this, SIGNAL(currentChanged(int, int)),
          this, SLOT(handleCurrentCell(int, int)));

  connect(this, SIGNAL(valueChanged(int, int)),
          this, SLOT(slotCellChanged(int, int)));

  connect(this, SIGNAL(signalChanged(int, int, QString)),
          parent, SLOT(slotTableChanged(int, int, QString)));

  connect(this, SIGNAL(parameterStatusChanged(int, bool)),
          parent, SLOT(slotParameterStatusChanged(int, bool)));
}

void ParameterTable::initTable()
{
  setNumRows(0);
  verticalHeader()->hide();
  setLeftMargin(0);
  //setSelectionMode(QTable::Single);
  setSelectionMode(Q3Table::NoSelection);
  //setFocusStyle(QTable::FollowStyle);

  //setNumRows(3);
  setNumCols(5);
  horizontalHeader()->setLabel(0, "Description");
  horizontalHeader()->setLabel(1, "Name");
  horizontalHeader()->setLabel(2, "");
  horizontalHeader()->setLabel(3, "Value");
  horizontalHeader()->setLabel(4, "Unit");
  setColumnStretchable(0, false);
  setColumnStretchable(1, true);
  setColumnStretchable(2, false);
  setColumnStretchable(3, true);
  setColumnStretchable(4, false);

  setShowGrid(false);
}

const std::vector<std::string> ParameterTable::getListOfAllMetabNames(const CModel & model,
    const CReactionInterface & ri)
{
  std::vector<std::string> ret;

  ret.push_back("unknown");

  //first all the metabs in the model
  unsigned C_INT32 i, imax = model.getMetabolites().size();

  for (i = 0; i < imax; ++i)
    ret.push_back(CMetabNameInterface::getDisplayName(&model, *model.getMetabolites()[i]));

  //now all the metabs from the ReactionInterface (that have not yet been committed to the model)
  std::vector<std::string> lll;
  std::vector<std::string>::const_iterator sourceIt, sourceItEnd;
  std::vector<std::string>::const_iterator searchIt, searchItEnd;

  lll = ri.getListOfMetabs(CFunctionParameter::SUBSTRATE);
  sourceItEnd = lll.end();

  for (sourceIt = lll.begin(); sourceIt != sourceItEnd; ++sourceIt)
    {
      searchItEnd = ret.end();

      for (searchIt = ret.begin(); searchIt != searchItEnd; ++searchIt)
        if (*searchIt == *sourceIt)
          break;

      if (searchIt == searchItEnd) //that means new metab name is not in model yet
        ret.push_back(*sourceIt);
    }

  lll = ri.getListOfMetabs(CFunctionParameter::PRODUCT);
  sourceItEnd = lll.end();

  for (sourceIt = lll.begin(); sourceIt != sourceItEnd; ++sourceIt)
    {
      searchItEnd = ret.end();

      for (searchIt = ret.begin(); searchIt != searchItEnd; ++searchIt)
        if (*searchIt == *sourceIt)
          break;

      if (searchIt == searchItEnd) //that means new metab name is not in model yet
        ret.push_back(*sourceIt);
    }

  lll = ri.getListOfMetabs(CFunctionParameter::MODIFIER);
  sourceItEnd = lll.end();

  for (sourceIt = lll.begin(); sourceIt != sourceItEnd; ++sourceIt)
    {
      searchItEnd = ret.end();

      for (searchIt = ret.begin(); searchIt != searchItEnd; ++searchIt)
        if (*searchIt == *sourceIt)
          break;

      if (searchIt == searchItEnd) //that means new metab name is not in model yet
        ret.push_back(*sourceIt);
    }

  return ret;
}

//static
QStringList ParameterTable::getListOfAllGlobalParameterNames(const CModel & model)
{
  QStringList ret;

  ret += "unknown";

  //all the global paramters  in the model
  unsigned C_INT32 i, imax = model.getNumModelValues();

  for (i = 0; i < imax; ++i)
    ret += FROM_UTF8(model.getModelValues()[i]->getObjectName());

  return ret;
}

//static
QStringList ParameterTable::getListOfAllCompartmentNames(const CModel & model)
{
  QStringList ret;

  ret += "unknown";

  //all the global paramters  in the model
  unsigned C_INT32 i, imax = model.getCompartments().size();

  for (i = 0; i < imax; ++i)
    ret += FROM_UTF8(model.getCompartments()[i]->getObjectName());

  return ret;
}

void ParameterTable::updateTable(const CReactionInterface & ri, const CModel & model)
{
  //first get the units strings
  CFindDimensions units(ri.getFunction(), model.getQuantityUnitEnum() == CModel::dimensionlessQuantity,
                        model.getVolumeUnitEnum() == CModel::dimensionlessVolume,
                        model.getTimeUnitEnum() == CModel::dimensionlessTime,
                        model.getAreaUnitEnum() == CModel::dimensionlessArea,
                        model.getLengthUnitEnum() == CModel::dimensionlessLength
                       );
  units.setUseHeuristics(true);
  units.setMolecularitiesForMassAction(ri.getChemEqInterface().getMolecularity(CFunctionParameter::SUBSTRATE),
                                       ri.getChemEqInterface().getMolecularity(CFunctionParameter::PRODUCT));
  units.findDimensions(ri.isMulticompartment());

  C_INT32 i, imax = ri.size();
  C_INT32 j, jmax;
  C_INT32 rowCounter = 0;

  //ColorTableItem *item;
  Q3TableItem *item;
  //ComboItem *combo;
  Q3ComboTableItem *combo;
  QStringList qsl;

  QColor subsColor(255, 210, 210);
  QColor prodColor(210, 255, 210);
  QColor modiColor(250, 250, 190);
  QColor paraColor(210, 210, 255);
  QColor volColor(210, 210, 255);
  QColor timeColor(210, 210, 210);

  QPixmap * pProduct = new QPixmap((const char**)product_xpm);
  QPixmap * pSubstrate = new QPixmap((const char**)substrate_xpm);
  QPixmap * pModifier = new QPixmap((const char**)modifier_xpm);
  QPixmap * pLocked = new QPixmap((const char**)locked_xpm);
  QPixmap * pUnlocked = new QPixmap((const char**)unlocked_xpm);

  CFunctionParameter::Role usage;
  QString qUsage;
  QColor color;
  const std::vector<std::string> * metabNames;

  mIndex2Line.resize(imax);
  mLine2Index.clear();

  setNumRows(0); // this is a hack to clear the table.
  setNumRows(imax*2);

  for (i = 0; i < imax; ++i)
    {
      // add additional space
      clearCell(rowCounter, 0); clearCell(rowCounter, 1); clearCell(rowCounter, 2);
      setRowReadOnly(rowCounter, true);
      setRowHeight(rowCounter++, 8);

      //
      setRowReadOnly(rowCounter, false);
      mIndex2Line[i] = rowCounter;

      // set the stuff that is different for the specific usages
      usage = ri.getUsage(i);
      qUsage = FROM_UTF8(CFunctionParameter::RoleNameDisplay[usage]);

      switch (usage)
        {
          case CFunctionParameter::SUBSTRATE:
            color = subsColor;
            break;
          case CFunctionParameter::PRODUCT:
            color = prodColor;
            break;
          case CFunctionParameter::MODIFIER:
            color = modiColor;
            break;
          case CFunctionParameter::PARAMETER:
            color = paraColor;
            break;
          case CFunctionParameter::VOLUME:
            color = volColor;
            break;
          case CFunctionParameter::TIME:
            color = timeColor;
            break;
          case CFunctionParameter::VARIABLE:
            color = QColor(255, 20, 20);
            break;
          default :
            qUsage = "unknown";
            color = QColor(255, 20, 20);
        }

      // add first column
      item = new ColorTableItem(this, Q3TableItem::Never, color, qUsage);

      if (usage == CFunctionParameter::SUBSTRATE) item->setPixmap(*pSubstrate);

      if (usage == CFunctionParameter::PRODUCT) item->setPixmap(*pProduct);

      if (usage == CFunctionParameter::MODIFIER) item->setPixmap(*pModifier);

      setItem(rowCounter, 0, item);

      // add second column
      item = new ColorTableItem(this, Q3TableItem::Never, color, FROM_UTF8(ri.getParameterName(i)));

      if ((usage != CFunctionParameter::PARAMETER)
          && (usage != CFunctionParameter::VOLUME)
          && (usage != CFunctionParameter::TIME))
        {
          if (ri.isLocked(i)) item->setPixmap(*pLocked); else item->setPixmap(*pUnlocked);
        }

      setItem(rowCounter, 1, item);

      // add third column
      if (usage == CFunctionParameter::PARAMETER)
        {
          item = new ColorCheckTableItem(this, color, "global");
          dynamic_cast<ColorCheckTableItem*>(item)->setChecked(!ri.isLocalValue(i));
        }
      else
        {
          item = new ColorTableItem(this, Q3TableItem::Never, color, "");
        }

      setItem(rowCounter, 2, item);

      // add units column
      assert(CCopasiRootContainer::getDatamodelList()->size() > 0);
      CCopasiDataModel* pDataModel = (*CCopasiRootContainer::getDatamodelList())[0];
      assert(pDataModel != NULL);
      item = new ColorTableItem(this, Q3TableItem::Never, color,
                                FROM_UTF8(" " + units.getDimensions()[i].getDisplayString(pDataModel)));
      setItem(rowCounter, 4, item);

      // add a line for a metabolite Parameter
      if ((usage == CFunctionParameter::SUBSTRATE)
          || (usage == CFunctionParameter::PRODUCT)
          || (usage == CFunctionParameter::MODIFIER))
        {
          // get the list of possible metabs (for the combo box)
          if (usage == CFunctionParameter::MODIFIER) //get all metabs; modifiers are never locked
            vectorOfStrings2QStringList(getListOfAllMetabNames(model, ri), qsl);
          else //only get the modifiers from the ChemEq
            {
              if (!ri.isLocked(i))
                vectorOfStrings2QStringList(ri.getListOfMetabs(usage), qsl);
            }

          metabNames = &(ri.getMappings(i));

          if (!ri.isVector(i))
            {
              if (ri.isLocked(i))
                {
                  item = new ColorTableItem(this, Q3TableItem::Never, color, FROM_UTF8((*metabNames)[0]));
                  setItem(rowCounter, 3, item);
                }
              else
                {
                  //combo = new ComboItem(this, QTableItem::WhenCurrent, color, qsl);
                  combo = new Q3ComboTableItem(this, qsl);
                  //combo->setText(FROM_UTF8((*metabNames)[0]));
                  combo->setCurrentItem(FROM_UTF8((*metabNames)[0]));
                  setItem(rowCounter, 3, combo);
                }
            }
          else
            {
              if (ri.isLocked(i))
                {
                  item = new ColorTableItem(this, Q3TableItem::Never, color, "");
                  setItem(rowCounter, 3, item);
                }
              else // this should not happen
                {
                  //combo = new ComboItem(this, QTableItem::WhenCurrent, color, qsl);
                  combo = new Q3ComboTableItem(this, qsl);
                  //combo->setText("add species");
                  combo->setCurrentItem("add species");
                  setItem(rowCounter, 3, combo);
                }

              // add lines for vector parameters
              jmax = metabNames->size();
              setNumRows(numRows() + jmax);

              for (j = 0; j < jmax; ++j)
                {
                  ++rowCounter;
                  item = new ColorTableItem(this, Q3TableItem::Never, color, FROM_UTF8((*metabNames)[j]));
                  setItem(rowCounter, 3, item);
                }
            }
        }
      // add a line for a kinetic parameter
      else if (usage == CFunctionParameter::PARAMETER)
        {
          if (ri.isLocalValue(i))
            {
              item = new ColorTableItem(this, Q3TableItem::OnTyping, color, QString::number(ri.getLocalValue(i)));
              setItem(rowCounter, 3, item);
            }
          else //global parameter
            {
              combo = new Q3ComboTableItem(this, getListOfAllGlobalParameterNames(model));
              combo->setCurrentItem(FROM_UTF8(ri.getMapping(i)));
              setItem(rowCounter, 3, combo);
            }
        }
      // add a line for a kinetic parameter
      else if (usage == CFunctionParameter::VOLUME)
        {
          combo = new Q3ComboTableItem(this, getListOfAllCompartmentNames(model));
          combo->setCurrentItem(FROM_UTF8(ri.getMapping(i)));
          setItem(rowCounter, 3, combo);
        }
      // add a line for time
      else if (usage == CFunctionParameter::TIME)
        {
          item = new ColorTableItem(this, Q3TableItem::OnTyping, color, "");
          setItem(rowCounter, 3, item);
        }
      // add a line for an unknown role
      else
        {
          item = new ColorTableItem(this, Q3TableItem::OnTyping, color, QString::number(ri.getLocalValue(i)));
          setItem(rowCounter, 3, item);
        }

      adjustRow(rowCounter);

      //mLine2Index

      ++rowCounter;
    }

  adjustColumn(0);
  adjustColumn(2);
  adjustColumn(4);
}

void ParameterTable::handleCurrentCell(int row, int col)
{
  bool changed = false;

  int i, imax = mIndex2Line.size();

  for (i = 0; i < imax; ++i)
    if (mIndex2Line[i] - 1 == row)
      {
        changed = true;
        if ((mOldRow < row) || (row == 0)) ++row; else --row;

        break;
      }

  if (col > 3) {changed = true; col = 3;}

  if (col < 2) {changed = true; col = 2;}

  mOldRow = row;

  if (changed) setCurrentCell(row, col);

  //TODO: allow keyboard editing of col 2
}

void ParameterTable::slotCellChanged(int row, int col)
{
  // find the index of the parameter
  C_INT32 i, imax = mIndex2Line.size();

  for (i = imax - 1; i >= 0; --i)
    if (mIndex2Line[i] <= row) break;

  //handle the check boxes
  if (col == 2) //only checkboxes is this column
    {
      Q3CheckTableItem *tmp = dynamic_cast<Q3CheckTableItem*>(this->item(row, col));

      if (!tmp) return;

      /*if (tmp->isChecked())
        {
        }
      else
        {
        }*/
      emit parameterStatusChanged(i, !tmp->isChecked());
    }
  else
    {
      emit signalChanged(i, row - mIndex2Line[i], text(row, col));
    }
}

//**************************************************************************

ComboItem::ComboItem(Q3Table *t, EditType et, QColor c, const QStringList & sl)
    : ColorTableItem(t, et, c, "Yes"), cb(0)
{
  // we do not want this item to be replaced
  setReplaceable(false);
  mSL = sl;
}

QWidget *ComboItem::createEditor() const
{
  // create an editor - a combobox in our case
  ((ComboItem*)this)->cb = new QComboBox(false, table()->viewport());
  QObject::connect(cb, SIGNAL(activated(int)), table(), SLOT(doValueChanged()));
  cb->insertStringList(mSL);
  cb->setCurrentText(text());
  return cb;
}

void ComboItem::setContentFromEditor(QWidget *w)
{
  // the user changed the value of the combobox, so synchronize the
  // value of the item (its text), with the value of the combobox
  if (w->inherits("QComboBox"))
    setText(((QComboBox*)w)->currentText());
  else
    Q3TableItem::setContentFromEditor(w);
}

void ComboItem::setText(const QString &s)
{
  Q3TableItem::setText(s);
}

//**********************************************************************

ColorTableItem::ColorTableItem(Q3Table *t, EditType et, QColor c, const QString txt)
    : Q3TableItem(t, et, txt)
{
  color = c;
}

ColorTableItem::~ColorTableItem()
{}

void ColorTableItem::paint(QPainter *p, const QColorGroup &cg,
                           const QRect &cr, bool selected)
{
  QColorGroup g(cg);
  g.setColor(QColorGroup::Base, color);
  Q3TableItem::paint(p, g, cr, selected);
}

//**********************************************************************

ColorCheckTableItem::ColorCheckTableItem(Q3Table *t, QColor c, const QString txt)
    : Q3CheckTableItem(t, txt)
{
  color = c;
}

ColorCheckTableItem::~ColorCheckTableItem()
{}

void ColorCheckTableItem::paint(QPainter *p, const QColorGroup &cg,
                                const QRect &cr, bool selected)
{
  QColorGroup g(cg);
  g.setColor(QColorGroup::Base, color);
  Q3CheckTableItem::paint(p, g, cr, selected);
}
