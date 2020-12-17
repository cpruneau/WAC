// Author: Claude Pruneau   09/25/2019

/*************************************************************************
 * Copyright (C) 2019, Claude Pruneau.
 * All rights reserved.
 * Based on the ROOT package and environment
 *
 * For the licensing terms see LICENSE.
 *************************************************************************/

#ifndef WAC_TaskCollection
#define WAC_TaskCollection

/////////////////////////////////////////////////////////////////////////
//
// TaskCollection
//
// Class defining a nested task
//
/////////////////////////////////////////////////////////////////////////
#include <ostream>
#include "Task.hpp"

using namespace std;

class TaskCollection : public Task
{
public:


  // ====================
  // Member functions
  // ====================

  TaskCollection(const TString & name,
                 TaskConfiguration * configuration,
                 int             nTasksMax = 100,
                 LogLevel selectedLevel=Info);
  virtual ~TaskCollection();
  virtual void initialize();
  virtual void execute();
  virtual void finalize();
  virtual void reset();
  virtual void clear();
  virtual void savePartialResults();
  virtual void printTaskConfiguration(ostream & output);
  Task *  addTask(Task * task);

  // ====================
  // Member Data
  // ====================

protected:

  int     nTasksMax;
  int     nTasks;
  Task ** tasks;

  ClassDef(TaskCollection,0)
};


#endif /* TaskCollection_hpp */
