/**
 * @file    List.cpp
 * @brief   Simple, generic list class.
 * @author  Ben Bornstein
 *
 * $Id: List.cpp 11634 2010-08-03 03:57:18Z mhucka $
 * $HeadURL: https://sbml.svn.sourceforge.net/svnroot/sbml/trunk/libsbml/src/util/List.cpp $
 *
 *<!---------------------------------------------------------------------------
 * This file is part of libSBML.  Please visit http://sbml.org for more
 * information about SBML, and the latest version of libSBML.
 *
 * Copyright 2005-2010 California Institute of Technology.
 * Copyright 2002-2005 California Institute of Technology and
 *                     Japan Science and Technology Corporation.
 * 
 * This library is free software; you can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation.  A copy of the license agreement is provided
 * in the file named "LICENSE.txt" included with this software distribution and
 * also available online as http://sbml.org/software/libsbml/license.html
 *----------------------------------------------------------------------- -->*/

#include <sbml/util/List.h>

LIBSBML_CPP_NAMESPACE_BEGIN

/*
 * Creates a new List.
 */
List::List ():
    size(0)
  , head(0)
  , tail(0)
{
}


/*
 * Destroys the given List.
 *
 * This function does not delete List items.  It destroys only the List and
 * its constituent ListNodes (if any).
 *
 * Presumably, you either i) have pointers to the individual list items
 * elsewhere in your program and you want to keep them around for awhile
 * longer or ii) the list has no items (List.size() == 0).  If neither are
 * true, try List_freeItems() instead.
 */
List::~List ()
{
  ListNode *node;
  ListNode *temp;


  node = head;

  while (node != 0)
  {
    temp = node;
    node = node->next;

    delete temp;
  }
}


/*
 * Adds item to the end of this List.
 */
void
List::add (void *item)
{
  ListNode* node = new ListNode(item);


  if (head == 0)
  {
    head = node;
    tail = node;
  }
  else
  {
    tail->next = node;
    tail       = node;
  }

  size++;
}


/*
 * @return the number of items in this List for which predicate(item)
 * returns true.
 *
 * The typedef for ListItemPredicate is:
 *
 *   int (*ListItemPredicate) (const void *item);
 *
 * where a return value of non-zero represents true and zero represents
 * false.
 */
unsigned int
List::countIf (ListItemPredicate predicate) const
{
  unsigned int count = 0;
  ListNode     *node = head;


  while (node != 0)
  {
    if (predicate(node->item) != 0)
    {
      count++;
    }

    node = node->next;
  }

  return count;
}


/*
 * @return the first occurrence of item1 in this List or NULL if item was
 * not found.  ListItemComparator is a pointer to a function used to find
 * item.  The typedef for ListItemComparator is:
 *
 *   int (*ListItemComparator) (void *item1, void *item2);
 *
 * The return value semantics are the same as for strcmp:
 *
 *   -1    item1 <  item2,
 *    0    item1 == item 2
 *    1    item1 >  item2
 */
void *
List::find (const void *item1, ListItemComparator comparator) const
{
  void     *item2  = 0;
  ListNode *node   = head;


  while (node != 0)
  {
    if (comparator(item1, node->item) == 0)
    {
      item2 = node->item;
      break;
    }

    node = node->next;
  }

  return item2;
}


/*
 * @return a new List containing (pointers to) all items in this List for
 * which predicate(item) was true.
 *
 * The returned list may be empty.
 *
 * The caller owns the returned list (but not its constituent items) and is
 * responsible for deleting it.
 */
List *
List::findIf (ListItemPredicate predicate) const
{
  List     *result = new List();
  ListNode *node   = head;


  while (node != 0)
  {
    if (predicate(node->item) != 0)
    {
      result->add(node->item);
    }

    node = node->next;
  }

  return result;
}


/*
 * Returns the nth item in this List.  If n > List.size() returns 0.
 */
void *
List::get (unsigned int n) const
{
  ListNode* node = head;


  if (n >= size)
  {
    return 0;
  }

  /**
   * Special case to retreive last item in the list without a full list
   * traversal.
   */
  if (n == (size - 1))
  {
    node = tail;
  }
  else
  {
    /* Point node to the nth item. */
    while (n-- > 0)
    {
      node = node->next;
    }
  }

  return node->item;
}


/*
 * Adds item to the beginning of this List.
 */
void
List::prepend (void *item)
{
  ListNode* node = new ListNode(item);


  if (head == 0)
  {
    head = node;
    tail = node;
  }
  else
  {
    node->next = head;
    head       = node;
  }

  size++;
}


/*
 * Removes the nth item from this List and returns a pointer to it.  If n >
 * List.size() returns 0.
 */
void *
List::remove (unsigned int n)
{
  void*     item;
  ListNode* prev;
  ListNode* temp;
  ListNode* next;


  if (n >= size)
  {
    return 0;
  }

  /**
   * temp = node to be removed
   * prev = node before temp (or NULL if temp == list->head)
   * next = node after  temp (or NULL if temp == list->tail)
   */
  prev = 0;
  temp = head;
  next = temp->next;

  /**
   * Point temp to nth item.
   */
  while (n-- > 0)
  {
    prev = temp;
    temp = temp->next;
    next = temp->next;
  }

  /**
   * If the first item in the list is being removed, only list->head needs
   * to be updated to remove temp.  Otherwise, prev->next must "forget"
   * about temp and point to next instead.
   */
  if (head == temp)
  {
    head = next;
  }
  else
  {
    prev->next = next;
  }

  /**
   * Regardless of the restructuring above, if the last item in the list
   * has been removed, update list->tail.
   */
  if (tail == temp)
  {
    tail = prev;
  }

  item = temp->item;
  delete temp;

  size--;

  return item;
}


/*
 * Returns the number of elements in this List.
 */
unsigned int
List::getSize () const
{
  return size;
}


/** @cond doxygen-c-only */


/**
 * Creates a new List and returns a pointer to it.
 */
LIBSBML_EXTERN
List_t *
List_create (void)
{
  return new List;
}


/**
 * Creates a new ListNode (with item) and returns a pointer to it.
 */
LIBSBML_EXTERN
ListNode_t *
ListNode_create (void *item)
{
  return new ListNode(item);
}


/**
 * Frees the given List.
 *
 * This function does not free List items.  It frees only the List_t
 * structure and its constituent ListNode_t structures (if any).
 *
 * Presumably, you either i) have pointers to the individual list items
 * elsewhere in your program and you want to keep them around for awhile
 * longer or ii) the list has no items (List_size(lst) == 0).  If neither
 * are true, try List_freeItems() instead.
 */
LIBSBML_EXTERN
void
List_free (List_t *lst)
{
  delete static_cast<List*>(lst);
}


/**
 * Frees the given ListNode.
 */
void
ListNode_free (ListNode_t *node)
{
  delete static_cast<ListNode*>(node);
}


/**
 * Adds item to the end of this List.
 */
LIBSBML_EXTERN
void
List_add (List_t *lst, void *item)
{
  static_cast<List*>(lst)->add(item);
}


/**
 * @return the number of items in this List for which predicate(item)
 * returns true.
 *
 * The typedef for ListItemPredicate is:
 *
 *   int (*ListItemPredicate) (const void *item);
 *
 * where a return value of non-zero represents true and zero represents
 * false.
 */
LIBSBML_EXTERN
unsigned int
List_countIf (const List_t *lst, ListItemPredicate predicate)
{
  return static_cast<const List*>(lst)->countIf(predicate);
}


/**
 * @return the first occurrence of item1 in this List or NULL if item was
 * not found.  ListItemComparator is a pointer to a function used to find
 * item.  The typedef for ListItemComparator is:
 *
 *   int (*ListItemComparator) (void *item1, void *item2);
 *
 * The return value semantics are the same as for strcmp:
 *
 *   -1    item1 <  item2,
 *    0    item1 == item 2
 *    1    item1 >  item2
 */
LIBSBML_EXTERN
void *
List_find ( const List_t *lst,
            const void   *item1,
            ListItemComparator comparator )
{
  return static_cast<const List*>(lst)->find(item1, comparator);
}


/**
 * @return a new List containing (pointers to) all items in this List for
 * which predicate(item) was true.
 *
 * The returned list may be empty.
 *
 * The caller owns the returned list (but not its constituent items) and is
 * responsible for freeing it with List_free().
 */
LIBSBML_EXTERN
List_t *
List_findIf (const List_t *lst, ListItemPredicate predicate)
{
  return static_cast<const List*>(lst)->findIf(predicate);
}


/**
 * Returns the nth item in this List.  If n > List_size(lst) returns NULL.
 */
LIBSBML_EXTERN
void *
List_get (const List_t *lst, unsigned int n)
{
  return static_cast<const List*>(lst)->get(n);
}


/**
 * Adds item to the beginning of this List.
 */
LIBSBML_EXTERN
void
List_prepend (List_t *lst, void *item)
{
  static_cast<List*>(lst)->prepend(item);
}


/**
 * Removes the nth item from this List and returns a pointer to it.  If n >
 * List_size(lst) returns NULL.
 */
LIBSBML_EXTERN
void *
List_remove (List_t *lst, unsigned int n)
{
  return static_cast<List*>(lst)->remove(n);
}


/**
 * @return the number of elements in this List.
 */
LIBSBML_EXTERN
unsigned int
List_size (const List_t *lst)
{
  return static_cast<const List*>(lst)->getSize();
}

LIBSBML_CPP_NAMESPACE_END

/** @endcond */
