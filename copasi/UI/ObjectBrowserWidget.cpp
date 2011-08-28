// Begin CVS Header
//   $Source: /fs/turing/cvs/copasi_dev/copasi/UI/ObjectBrowserWidget.cpp,v $
//   $Revision: 1.23 $
//   $Name: Build-33 $
//   $Author: aekamal $
//   $Date: 2009/08/10 15:14:44 $
// End CVS Header

// Copyright (C) 2008 by Pedro Mendes, Virginia Tech Intellectual
// Properties, Inc., EML Research, gGmbH, University of Heidelberg,
// and The University of Manchester.
// All rights reserved.

// Copyright (C) 2001 - 2007 by Pedro Mendes, Virginia Tech Intellectual
// Properties, Inc. and EML Research, gGmbH.
// All rights reserved.

/********************************************************
Author: Liang Xu
Version : 1.xx  <first>
Description:
Date: 04/03
Comment : Copasi Object Browser:
Contact: Please contact lixu1@vt.edu.
 *********************************************************/

#include <qmessagebox.h>
#include <qvariant.h>
#include <q3header.h>
#include <q3listview.h>
#include <q3textedit.h>
#include <q3listbox.h>
#include <qpushbutton.h>
#include <qlayout.h>
#include <qtooltip.h>
#include <q3whatsthis.h>
#include <qimage.h>
#include <qpixmap.h>
#include <q3simplerichtext.h>
#include <q3frame.h>
//Added by qt3to4:
#include <Q3GridLayout>

#include "copasi.h"

#include "ObjectBrowserWidget.h"
#include "ObjectBrowserItem.h"
#include "copasiui3window.h"
#include "qtUtilities.h"
#include "report/CCopasiObject.h"
#include "report/CCopasiObjectName.h"
#include "report/CCopasiContainer.h"
#include "report/CCopasiStaticString.h"
#include "utilities/CCopasiVector.h"
#include "report/CCopasiRootContainer.h"

#include "./icons/objectAll.xpm"
#include "./icons/objectParts.xpm"
#include "./icons/objectNone.xpm"

QPixmap *pObjectAll = 0;
QPixmap *pObjectParts = 0;
QPixmap *pObjectNone = 0;

/*
 *  Constructs a ObjectBrowser which is a child of 'parent', with the
 *  name 'name' and widget flags set to 'f'.
 */
ObjectBrowserWidget::ObjectBrowserWidget(QWidget* parent, const char* name, Qt::WFlags fl, int state)
    : QWidget(parent, name, fl),
    objectItemList(NULL),
    refreshList(NULL),
    ObjectBrowserLayout(NULL),
    clearButton(NULL),
    toggleViewButton(NULL),
    commitButton(NULL),
    ObjectListView(NULL),
    Line1(NULL),
    spacer(NULL),
    ObjectItemText(NULL),
    mOutputObjectVector(NULL),
    currentPage(LISTVIEWPAGE)
{
  if (!name)
    setName("ObjectBrowser");

  if (state == 0)
    ObjectBrowserLayout = new Q3GridLayout(this, 1, 1, 0, -1, "ObjectBrowserLayout");
  else
    ObjectBrowserLayout = new Q3GridLayout(this, 2, 4, 0, 6, "ObjectBrowserLayout");

  ObjectBrowserLayout->setAutoAdd(false);
  ObjectListView = new Q3ListView(this, "ObjectListView");
  //  ObjectListView->addColumn(trUtf8("Object Browser"));
  ObjectListView->addColumn(trUtf8("Name"));
  ObjectListView->addColumn(trUtf8("Type"));
  ObjectListView->header()->setClickEnabled(false, ObjectListView->header()->count() - 1);
  ObjectListView->setAcceptDrops(false);
  ObjectListView->setResizeMode(Q3ListView::LastColumn);
  ObjectListView->setTreeStepSize(19);
  //  ObjectListView->setSorting(-1);
  ObjectListView->setSelectionMode(Q3ListView::NoSelection);

  ObjectItemText = new Q3TextEdit(this, "ObjectItemText");
  ObjectItemText ->hide();

  ObjectBrowserLayout->addMultiCellWidget(ObjectListView, 0, 0, 0, 3);
  ObjectBrowserLayout->addMultiCellWidget(ObjectItemText, 0, 0, 0, 3);

  if (state != 0)
    {
      Line1 = new Q3Frame(this, "Line1");
      Line1->setFrameShape(Q3Frame::HLine);
      Line1->setFrameShadow(Q3Frame::Sunken);
      Line1->setFrameShape(Q3Frame::HLine);

      ObjectBrowserLayout->addMultiCellWidget(Line1, 1, 1, 0, 3);

      clearButton = new QPushButton(this, "clearButton");
      clearButton->setText(trUtf8("Clear"));
      ObjectBrowserLayout->addWidget(clearButton, 2, 0);

      commitButton = new QPushButton(this, "commitButton");
      commitButton->setText(trUtf8("Commit"));
      ObjectBrowserLayout->addWidget(commitButton, 2, 3);

      toggleViewButton = new QPushButton(this, "toggleViewButton");
      toggleViewButton->setText(trUtf8("Selected Items"));
      ObjectBrowserLayout->addWidget(toggleViewButton, 2, 2);

      spacer = new QSpacerItem(131, 31, QSizePolicy::Expanding, QSizePolicy::Minimum);
      ObjectBrowserLayout->addItem(spacer, 2, 1);

      // signals and slots connections
      connect(clearButton, SIGNAL(clicked()), this, SLOT(clearClicked()));
      connect(toggleViewButton, SIGNAL(clicked()), this, SLOT(toggleViewClicked()));
      connect(commitButton, SIGNAL(clicked()), this, SLOT(commitClicked()));

      // tab order
      setTabOrder(ObjectListView, clearButton);
      setTabOrder(clearButton, toggleViewButton);
      setTabOrder(toggleViewButton, commitButton);
    }

  connect(ObjectListView, SIGNAL(clicked(Q3ListViewItem*)), this, SLOT(listviewChecked(Q3ListViewItem*)));

  pObjectAll = new QPixmap((const char**)ptrObjectAll);
  pObjectNone = new QPixmap((const char**)ptrObjectNone);
  pObjectParts = new QPixmap((const char**)ptrObjectParts);

  objectItemList = new ObjectList();
  refreshList = new ObjectList();
  ObjectBrowserItem::resetKeySpace();
  loadData();
  currentPage = LISTVIEWPAGE;
}

ObjectBrowserWidget::~ObjectBrowserWidget()
{
  cleanup();
}

void ObjectBrowserWidget::cleanup()
{
  delete objectItemList;
  delete refreshList;
}

void ObjectBrowserWidget::clearClicked()
{
  if (!ObjectListView->firstChild()) return;

  ObjectBrowserItem* root = static_cast< ObjectBrowserItem * >(ObjectListView->firstChild());
  setUncheck(root);
  updateUI();
  //  qWarning("ObjectBrowser::clearClicked(): Not implemented yet!");
}

void ObjectBrowserWidget::listviewChecked(Q3ListViewItem* pCurrent)
{
  //  if ((pCurrent == NULL)||((static_cast< ObjectBrowserItem * >(pCurrent))->getType() == FIELDATTR))
  if ((pCurrent == NULL) || (pCurrent->text(0) == "Select by attribute"))
    return;

  clickToReverseCheck(static_cast< ObjectBrowserItem * >(pCurrent));
  updateUI();
  //loadUI();
}

void ObjectBrowserWidget::clickToReverseCheck(ObjectBrowserItem* pCurrent)
{
  refreshList->insert(pCurrent);

  ObjectBrowserItem* pTmp = pCurrent;

  while (pTmp->parent() != NULL)
    {
      pTmp = static_cast< ObjectBrowserItem * >(pTmp->parent());
      refreshList->insert(pTmp);
    }

  if (pCurrent->isChecked() == ALLCHECKED)
    {
      if (pCurrent->isChecked())
        pCurrent->reverseChecked();

      setUncheck(static_cast< ObjectBrowserItem * >(pCurrent->firstChild()));
      return;
    }

  //else no check or partly checked
  if (!pCurrent->isChecked())
    pCurrent->reverseChecked();

  setCheck(static_cast< ObjectBrowserItem * >(pCurrent->firstChild()));
  return;
}

void ObjectBrowserWidget::setUncheck(ObjectBrowserItem* pCurrent)
{
  if (pCurrent == NULL)
    return;

  refreshList->insert(pCurrent);

  if (pCurrent->isChecked())
    pCurrent->reverseChecked();

  if (pCurrent->firstChild() != NULL)
    setUncheck(static_cast< ObjectBrowserItem * >(pCurrent->firstChild()));

  if (pCurrent->nextSibling() != NULL)
    setUncheck(static_cast< ObjectBrowserItem * >(pCurrent->nextSibling()));
}

void ObjectBrowserWidget::setCheck(ObjectBrowserItem* pCurrent)
{
  if (pCurrent == NULL)
    return;

  refreshList->insert(pCurrent);

  if (!pCurrent->isChecked())
    pCurrent->reverseChecked();

  if (pCurrent->firstChild() != NULL)
    setCheck(static_cast< ObjectBrowserItem * >(pCurrent->firstChild()));

  if (pCurrent->nextSibling() != NULL)
    setCheck(static_cast< ObjectBrowserItem * >(pCurrent->nextSibling()));
}

void ObjectBrowserWidget::toggleViewClicked()
{
  switch (currentPage)
    {
      case LISTVIEWPAGE:
        updateSelectedItemsView();
        currentPage = SELECTEDITEMPAGE;

        if (toggleViewButton) toggleViewButton->setText(trUtf8("Tree View"));

        ObjectListView->hide();
        ObjectItemText->show();
        break;
      case SELECTEDITEMPAGE:
        currentPage = LISTVIEWPAGE;

        if (toggleViewButton) toggleViewButton->setText(trUtf8("Selected Items"));

        ObjectListView->show();
        ObjectItemText->hide();
        break;
    }
}

void ObjectBrowserWidget::updateSelectedItemsView()
{
  std::vector< const CCopasiObject * > * outputVector;
  ObjectBrowserItem* rootItem;
  unsigned C_INT32 i;
  rootItem = objectItemList->getRoot()->pItem;
  outputVector = new std::vector< const CCopasiObject * >();
  eXport(rootItem, outputVector);
  //      CQMessageBox::information(this, "Output object list done!", "Selected CopasiObject list done!");
  //   ObjectListItem* pHead;
  ObjectItemText->clear();

  for (i = 0; i < outputVector->size(); i++)
    {
      //   if (double(i) / 2 == int(i / 2))
      //     ObjectItemText->setColor(red);
      //   else
      //     ObjectItemText->setColor(blue);
      if ((*outputVector)[i])
        {
          ObjectItemText->insertParagraph(FROM_UTF8((*outputVector)[i]->getObjectDisplayName()), -1);
          //ObjectItemText->insertParagraph(pHead->pItem->getObject()->pCopasiObject->getCN()., -1);
          //ObjectItemText->insertParagraph(FROM_UTF8((*outputVector)[i]->getObjectType()), -1);
        }
    }

  /*
  for (pHead = outputVector->getRoot(), i = 1; pHead != NULL; pHead = pHead->pNext)
  {
  if (double(i) / 2 == int(i / 2))
  ObjectItemText->setColor(red);
  else
  ObjectItemText->setColor(blue);
  if (pHead->pItem->getObject()->pCopasiObject)
  {
  //ObjectItemText->insertParagraph(pHead->pItem->getObject()->pCopasiObject->getCN()., -1);
  ObjectItemText->insertParagraph(pHead->pItem->getObject()->pCopasiObject->getObjectType()., -1);
  i++;
  }
  }
  */
  pdelete(outputVector);
}

void ObjectBrowserWidget::setOutputVector(std::vector< const CCopasiObject * > * pObjectVector)
{
  mOutputObjectVector = pObjectVector;
  this->clearClicked();

  if (mOutputObjectVector) selectObjects(mOutputObjectVector);
}

void ObjectBrowserWidget::commitClicked()
{
  ObjectBrowserItem* rootItem;
  rootItem = objectItemList->getRoot()->pItem;

  if (mOutputObjectVector) mOutputObjectVector->clear();

  eXport(rootItem, mOutputObjectVector);
  return;
}

void ObjectBrowserWidget::eXport(ObjectBrowserItem* pCurrent, std::vector< const CCopasiObject * > * outputVector)
{
  if (!outputVector) return;

  if (pCurrent->firstChild())
    {
      ObjectBrowserItem* pChild = static_cast< ObjectBrowserItem * >(pCurrent->firstChild());

      for (; pChild != NULL; pChild = static_cast< ObjectBrowserItem * >(pChild->nextSibling()))
        if (pChild->getType() != FIELDATTR)
          eXport(pChild, outputVector);
    }
  else //it has no child
    {
      if (pCurrent->isChecked() && (pCurrent->getType() != FIELDATTR))
        {
          if (pCurrent->getObject())
            {
              bool insertFlag = true;
              unsigned C_INT32 i = 0;

              for (; i < outputVector->size(); i++)
                if (pCurrent->getObject()->pCopasiObject->getCN() == (*outputVector)[i]->getCN()) insertFlag = false;

              if (insertFlag)
                outputVector->push_back(pCurrent->getObject()->pCopasiObject);
            }
        }

      // else skip current item
    }
}

void ObjectBrowserWidget::swap(int first_pos, int second_pos, ObjectBrowserItem** array)
{
  ObjectBrowserItem* tmp = array[first_pos];
  array[first_pos] = array[second_pos];
  array[second_pos] = tmp;
}

int ObjectBrowserWidget::partition(int split_pos, int start_pos, int end_pos, ObjectBrowserItem** quick_sort_array)
{
  const void * tmp = quick_sort_array[split_pos]->getObject()->pCopasiObject;

  if (split_pos != start_pos)
    swap(split_pos, start_pos, quick_sort_array);

  split_pos = start_pos;
  end_pos++;

  while (start_pos < end_pos)
    {
      while ((++start_pos < end_pos) && (quick_sort_array[start_pos]->getObject()->pCopasiObject <= tmp));

      while ((--end_pos > start_pos) && (quick_sort_array[end_pos]->getObject()->pCopasiObject > tmp));

      if (start_pos < end_pos)
        swap(start_pos, end_pos, quick_sort_array);
      else
        break;
    }

  start_pos--;

  if (split_pos != start_pos)
    swap(split_pos, start_pos, quick_sort_array);

  return start_pos;
}

void ObjectBrowserWidget::quick_sort(int m, int n, ObjectBrowserItem** quick_sort_array)
{
  if (m < n)
    {
      int medium = (m + n) / 2;
      medium = partition(medium, m, n, quick_sort_array);

      if (m < medium - 1)
        quick_sort(m, medium - 1, quick_sort_array);

      if (medium + 1 < n)
        quick_sort(medium + 1, n, quick_sort_array);
    }
}

void ObjectBrowserWidget::loadData()
{
  const CCopasiContainer * root = CCopasiRootContainer::getRoot();
  ObjectBrowserItem * itemRoot = new ObjectBrowserItem(ObjectListView, NULL, root, objectItemList);
  itemRoot->attachKey();
  itemRoot->setObjectType(CONTAINERATTR);
  //  itemRoot->setText(0, FROM_UTF8(root->getObjectName()));
  itemRoot->setText(0, QString("COPASI"));
  itemRoot->setOpen(true);
  loadChild(itemRoot, root, true);
  removeDuplicate(objectItemList);
  loadUI();
}

void ObjectBrowserWidget::removeDuplicate(ObjectList* objectItemList)
{
  int length = objectItemList->len();
  ObjectBrowserItem** bufferVector = new ObjectBrowserItem * [length];

  // copy out the list
  ObjectListItem* pCurrent = objectItemList->getRoot();

  for (int i = 0; i < length; i++)
    {
      bufferVector[i] = pCurrent->pItem;
      pCurrent = pCurrent->pNext;
    }

  // sort list according to CCopasiObject
  quick_sort(0, length - 1, bufferVector);

  int index;
  CBrowserObject* pBrowserObject = NULL;
  ObjectBrowserItem* pBrowserItem = NULL;

  for (index = 0; (index < length) && (!bufferVector[index]->getObject()->pCopasiObject); index++);

  pBrowserObject = bufferVector[index++]->getObject();

  for (; index < length; index++)
    {
      if (bufferVector[index]->getObject()->pCopasiObject == pBrowserObject->pCopasiObject) //duplicate point
        {
          pBrowserItem = bufferVector[index]->getObject()->referenceList->pop();
          pBrowserObject->referenceList->insert(pBrowserItem);
          pBrowserItem->setBrowserObject(pBrowserObject);
        }
      else
        pBrowserObject = bufferVector[index]->getObject();
    }
}

void ObjectBrowserWidget::loadChild(ObjectBrowserItem* parent,
                                    const CCopasiContainer* copaParent,
                                    bool nField)
{
  unsigned int i;
  ObjectBrowserItem* last = NULL;
  CCopasiObject* current = NULL;

  ObjectList* childStack = new ObjectList();

  const CCopasiContainer::objectMap * pObjectList = & copaParent->getObjects();
  CCopasiContainer::objectMap::const_iterator it = pObjectList->begin();
  CCopasiContainer::objectMap::const_iterator end = pObjectList->end();

  if ((copaParent->isVector()) && (nField))
    {
      if ((static_cast< const CCopasiVector < CCopasiObject > * >(copaParent)->size() >= 1) &&
          ((*static_cast< const CCopasiVector < CCopasiObject > * >(copaParent))[0]->isContainer()))
        {//add attribute list
          ObjectBrowserItem* fieldChild = new ObjectBrowserItem(parent, NULL, NULL, objectItemList);
          fieldChild->setObjectType(FIELDATTR);
          fieldChild->setText(0, "Select by attribute");
          fieldChild->setSelectable(false);
          loadField(fieldChild, const_cast<CCopasiVector < CCopasiObject > *>(static_cast< const CCopasiVector < CCopasiObject > * >(copaParent)));
          fieldChild->attachKey();
          last = fieldChild;
        }
    }

  if (copaParent->isVector())
    {
      for (i = 0; i < static_cast< const CCopasiVector < CCopasiObject > * >(copaParent)->size(); i++)
        {
          current = (*static_cast< const CCopasiVector < CCopasiObject > * >(copaParent))[i];
          ObjectBrowserItem* currentItem = new ObjectBrowserItem(parent, last, current, objectItemList);
          last = currentItem;
          currentItem->setText(0, FROM_UTF8(current->getObjectName()));

          if (current->isContainer())
            {
              currentItem->setObjectType(CONTAINERATTR);
              currentItem->attachKey();

              if (current->isVector())
                currentItem->setText(0, currentItem->text(0) + "[]");

              loadChild(currentItem, static_cast< CCopasiContainer * >(current), nField);
            }
          else
            {
              currentItem->setObjectType(OBJECTATTR);
              childStack->insert(currentItem);  //attach the key later
            }
        }
    }
  else
    {
      while (it != end)
        {
          current = it->second;

          // Skip all strings
          if (dynamic_cast<CCopasiStaticString *>(current))
            {
              it++;
              continue;
            }

          ObjectBrowserItem* currentItem = new ObjectBrowserItem(parent, last, current, objectItemList);
          last = currentItem;
          currentItem->setText(0, FROM_UTF8(current->getObjectName()));

          if (current->isContainer())
            {
              currentItem->setObjectType(CONTAINERATTR);
              currentItem->attachKey();

              if (current->isVector())
                currentItem->setText(0, currentItem->text(0) + "[]");

              loadChild(currentItem, static_cast< CCopasiContainer * >(current), nField);
            }
          else
            {
              currentItem->setObjectType(OBJECTATTR);
              childStack->insert(currentItem);  //attach the key later
            }

          it++;
        }
    }

  ObjectBrowserItem* pCurrent;

  while (childStack->len() > 0)
    {
      pCurrent = childStack->pop();
      pCurrent->attachKey();
    }

  pdelete(childStack);
}

void ObjectBrowserWidget::loadField(ObjectBrowserItem* parent, CCopasiVector <CCopasiObject>* copaParent)
{
  unsigned int i;
  ObjectBrowserItem* lastFieldItem = NULL;
  CCopasiObject* currentFieldObject = NULL;
  ObjectBrowserItem* lastObjectItem = NULL;
  CCopasiObject* currentObject = NULL;

  if ((copaParent->size() < 1) || (!(*copaParent)[0]->isContainer())) return; //empty list

  const CCopasiContainer::objectMap * pFieldList =
    &(static_cast< CCopasiContainer * >((*copaParent)[0])->getObjects());

  CCopasiContainer::objectMap::const_iterator fieldIt = pFieldList->begin();
  CCopasiContainer::objectMap::const_iterator fieldEnd = pFieldList->end();

  while (fieldIt != fieldEnd)
    {
      currentFieldObject = fieldIt->second;
      ObjectBrowserItem* currentFieldItem = new ObjectBrowserItem(parent, lastFieldItem, NULL, objectItemList);
      currentFieldItem->attachKey();
      currentFieldItem->setObjectType(FIELDATTR);
      currentFieldItem->setText(0, FROM_UTF8(currentFieldObject->getObjectName()));
      lastFieldItem = currentFieldItem;
      lastObjectItem = NULL;

      for (i = 0; i < copaParent->size(); i++)
        {
          CCopasiObject* pSubField;
          currentObject = (*static_cast< CCopasiVector < CCopasiObject > * >(copaParent))[i];

          if (currentObject->isContainer())
            pSubField =
              getFieldCopasiObject(static_cast< CCopasiContainer * >(currentObject),
                                   FROM_UTF8(currentFieldObject->getObjectName()));
          else
            {
              pSubField = NULL; // this shall be an exception error
            }

          ObjectBrowserItem* currentItem = new ObjectBrowserItem(currentFieldItem, lastObjectItem, pSubField, objectItemList);
          currentItem->setText(0, FROM_UTF8(currentObject->getObjectName()));
          //   if ((pSubField)&&(pSubField->isVector()))
          //            currentItem->setText(0, currentItem->text(0) + "[]");
          currentItem->setObjectType(FIELDATTR);
          currentItem->attachKey();

          if (pSubField)
            if (pSubField->isContainer()) //container field recursively to find all possible items
              {
                loadChild(currentItem, static_cast< CCopasiContainer * >(pSubField), false); // wont show the attribute and field list
              }

          lastObjectItem = currentItem;
        }

      fieldIt++;
    }
}

void ObjectBrowserWidget::updateUI()
{
  //refresh List stores all affected items,

  refreshList->createBucketIndex((int)ObjectBrowserItem::getKeySpace()); //construct index to do binary search

  for (ObjectListItem* pCurrent = refreshList->getRoot(); pCurrent != NULL; pCurrent = pCurrent->pNext)
    {
      ObjectListItem * pHead = pCurrent->pItem->getObject()->referenceList->getRoot();

      while (pHead != NULL)
        {
          ObjectBrowserItem * pCurrentLevel = pHead->pItem;

          if (pCurrent != pHead)
            for (; (pCurrentLevel != NULL); pCurrentLevel =
                   static_cast< ObjectBrowserItem * >(pCurrentLevel->parent()))
              refreshList->insertBucket(pCurrentLevel);

          pHead = pHead->pNext;
        }
    }

  int cursor = 0;

  for (ObjectBrowserItem* pUpdate = refreshList->bucketPop(cursor); pUpdate != NULL; pUpdate = refreshList->bucketPop(cursor))
    setCheckMark(pUpdate);

  refreshList->destroyBucket();

  if (currentPage == SELECTEDITEMPAGE) updateSelectedItemsView();
}

void ObjectBrowserWidget::setCheckMark(ObjectBrowserItem* pCurrent)
{
  if (pCurrent->text(0) == "Select by attribute")
    {
      pCurrent->setPixmap(0, NULL);
      return;
    }

  switch (pCurrent->nUserChecked())
    {
      case NOCHECKED:
        pCurrent->setPixmap(0, *pObjectNone);
        break;
      case ALLCHECKED:
        pCurrent->setPixmap(0, *pObjectAll);
        break;
      case PARTCHECKED:
        pCurrent->setPixmap(0, *pObjectParts);
        break;
    }
}

void ObjectBrowserWidget::loadUI()
{
  ObjectListItem* pCurrent = objectItemList->getRoot();
  setCheckMark(pCurrent->pItem);

  for (; pCurrent != NULL; pCurrent = pCurrent->pNext)
    setCheckMark(pCurrent->pItem);
}

CCopasiObject* ObjectBrowserWidget::getFieldCopasiObject(CCopasiContainer * pCurrent, const char* name)
{
  const CCopasiContainer::objectMap * pObjectList = & pCurrent->getObjects();
  CCopasiContainer::objectMap::const_iterator it = pObjectList->begin();
  CCopasiContainer::objectMap::const_iterator end = pObjectList->end();

  //  CCopasiObject* pResult;
  while (it != end)
    {
      //      if (it->second->isContainer())
      //        {
      //          pResult = getFieldCopasiObject((CCopasiContainer *)it->second, name);
      //          if (pResult)
      //            return pResult;
      //}
      if (FROM_UTF8(it->second->getObjectName()) == name)
        return it->second;

      it++;
    }

  return NULL;
}

void ObjectBrowserWidget::selectObjects(std::vector< const CCopasiObject * > * pObjectVector)
{
  unsigned int i;
  ObjectBrowserItem* rootItem;

  rootItem = objectItemList->getRoot()->pItem;

  for (i = 0; i < pObjectVector->size(); i++)
    {
      if ((*pObjectVector)[i]) selectObjects(rootItem, (*pObjectVector)[i]);
    }

  updateUI();
}

void ObjectBrowserWidget::selectObjects(ObjectBrowserItem* browserItem,
                                        const CCopasiObject * selectObject)
{
  ObjectBrowserItem* pCurrent;
  pCurrent = browserItem;

  if (pCurrent->firstChild())
    {
      ObjectBrowserItem* pChild = static_cast< ObjectBrowserItem * >(pCurrent->firstChild());

      for (; pChild != NULL; pChild = static_cast< ObjectBrowserItem * >(pChild->nextSibling()))
        if (pChild->getType() != FIELDATTR)
          selectObjects(pChild, selectObject);
    }
  else
    {//it has no child
      if (pCurrent->getType() != FIELDATTR)
        {
          if (pCurrent->getObject() && (pCurrent->getObject()->pCopasiObject == selectObject))
            {
              refreshList->insert(pCurrent);

              ObjectBrowserItem* pTmp = pCurrent;

              while (pTmp->parent() != NULL)
                {
                  pTmp = static_cast< ObjectBrowserItem * >(pTmp->parent());
                  refreshList->insert(pTmp);
                }

              pCurrent->getObject()->mChecked = ALLCHECKED;
              setCheck(static_cast< ObjectBrowserItem * >(pCurrent->firstChild()));
            }

          // else skip current item
        }
    }
}
