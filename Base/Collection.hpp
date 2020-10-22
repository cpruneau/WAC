// Author: Claude Pruneau   09/25/2019

/***********************************************************************
 * Copyright (C) 2019, Claude Pruneau.                                 *
 * All rights reserved.                                                *
 * Based on the ROOT package and environment                           *
 *                                                                     *
 * For the licensing terms see LICENSE.                                *
 **********************************************************************/

#ifndef WAC_Collection
#define WAC_Collection

// ====================================================
//
// Collection
//
// Generic base class for Collection services.
// Does not own a store
//
// ====================================================
#include <iostream>
#include "TObject.h"

using namespace std;

template < class T >
class Collection : public TObject
{
protected:
  /////////////////////////////////////////////////////////////
  // Data Members
  /////////////////////////////////////////////////////////////
  long nObjectCapacity;
  long nObject;
  bool bOwnObjects; // if the instance owns its own objects
                    // use to be the case unless the objects are incorporated to an external list
                    // when true at instance destruction time the objects are deleted
  T ** objects;

  public:
  /////////////////////////////////////////////////////////////
  // Member functions
  /////////////////////////////////////////////////////////////
  Collection(int initialCapacity=10)
  :
  TObject(),
  nObjectCapacity(0),
  nObject(0),
  bOwnObjects(true),
  objects(0)
  {
  initialize(initialCapacity);
  }

  Collection(const Collection & source)
  :
  TObject(source),
  nObjectCapacity(0),
  nObject(0),
  objects(0)
  {
  initialize(source.nObjectCapacity);
  bOwnObjects = true;
  deepcopy(source.nObject,source.objects);
  }


  virtual ~Collection()
  {
    if (bOwnObjects) {
      deleteObjects();
    }
    delete[] objects;
  }

  void SetOwn(bool ownobjects = true) {
    bOwnObjects = ownobjects;
  }

  Collection & operator=(const Collection source)
  {
  if (this==&source) return *this;
  if (bOwnObjects) {
    deleteObjects();
  }
  delete[] objects;
  bOwnObjects = source.bOwnObjects;
  initialize(source.nObjectCapacity);
  copy(source.nObject,source.objects);
  return *this;
  }

  T * operator[](long index)
  {
  if (index<0) return objects[0];
  if (index>=nObject) return objects[nObject-1];
  return objects[index];
  }

  // Returns the nObjectCapacity of this store/Collection
  long getCollectionCapacity() const
  {
    return nObjectCapacity;
  }

  // Returns the number of object currently
  // in use.
  long getCollectionSize() const
  {
    return nObject;
  }

T * add(T * object)
  {
  if (nObject >= nObjectCapacity)
    {
    T ** work = objects;
    initialize(2*nObjectCapacity);
    copy(nObject,work);
    delete[] work;
    }
  objects[nObject] = object;
  nObject++;
  return object;
  }

  T * getObjectAt(long index)
  {
  if (index<0) return objects[0];
  if (index>=nObject) return objects[nObject-1];
  return objects[index];
  }

protected:

  void initialize(int initialCapacity)
  {
    nObjectCapacity = initialCapacity;
    objects = new T*[nObjectCapacity];
  }

  void deleteObjects() {
    for(int i = 0; i < nObject; ++i) {
      delete objects[i];
    }
  }


void copy(long n, T ** source)
  {
  for (int iObject=0; iObject<n; iObject++)
    {
    //objects[iObject] = new T(*source[iObject]);
    objects[iObject] = source[iObject];
    }
  nObject = n;
  }

void deepcopy(long n, T ** source)
  {
  for (int iObject=0; iObject<n; iObject++)
    {
    objects[iObject] = (T*) source[iObject]->Clone();
    }
  nObject = n;
  }

  ClassDef(Collection,1)
};

#endif /* WAC_Collection */
