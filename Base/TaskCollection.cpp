// Author: Claude Pruneau   09/25/2019

/*************************************************************************
 * Copyright (C) 2019, Claude Pruneau.
 * All rights reserved.
 * Based on the ROOT package and environment
 *
 * For the licensing terms see LICENSE.  
 *************************************************************************/
/**
 \class TaskCollection
 \ingroup WAC

 Class defining Nested Tasks
 */

#include "TaskCollection.hpp"

ClassImp(TaskCollection);


///////////////////////////////////////////////////////////////////////////////////////////////////////////////
// CTOR
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
TaskCollection::TaskCollection(const TString & name,
                               TaskConfiguration * configuration,
                               int nTasksMaxSelected,
                               LogLevel selectedLevel)
:
Task( name, configuration, nullptr,selectedLevel),
nTasksMax(nTasksMaxSelected>1?nTasksMaxSelected:1),
nTasks(0),
tasks(0)
{
  if (reportStart("TaskCollection",getTaskName(),"CTOR"))
    ;
  tasks = new Task* [nTasksMax];
  if (reportEnd("TaskCollection",getTaskName(),"CTOR"))
    ;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
// DTOR
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
TaskCollection::~TaskCollection()
{
  if (reportStart("TaskCollection",getTaskName(),"DTOR"))
    ;
  delete[] tasks;
  if (reportEnd("TaskCollection",getTaskName(),"DTOR"))
    ;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Initialize task
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void TaskCollection::initialize()
{
  if (reportInfo("TaskCollection",getTaskName(),"initialize()"))
    cout << "Initializing " << nTasks << " tasks" << endl;
  for (int iTask=0; iTask<nTasks; iTask++)
    {
      if (isTaskOk())
        {
        if (!tasks[iTask])
          {
          if (reportFatal("TaskCollection",getTaskName(),"initialize()")) cout << "Null pointer given for iTask=" << iTask << endl;
          postTaskFatal();
          return;
          }
        }
    tasks[iTask] ->initialize();
    }
  if (reportEnd("TaskCollection",getTaskName(),"initialize()"))
    ;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Execute task
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void TaskCollection::execute()
{
  //if (reportStart("TaskCollection",getTaskName(),"execute()"));
  for (int iTask=0; iTask<nTasks; iTask++) if (isTaskOk()) tasks[iTask]->execute();
  //if (reportEnd("TaskCollection",getTaskName(),"execute()"));
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Finalize task
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void TaskCollection::finalize()
{
  if (reportStart("TaskCollection",getTaskName(),"finalize()"))
    ;
   for (int iTask=0; iTask<nTasks; iTask++) if (isTaskOk()) tasks[iTask]->finalize();
  if (reportEnd("TaskCollection",getTaskName(),"finalize()"))
    ;
}

////////////////////////////////////////////
// Reset this task
////////////////////////////////////////////
void TaskCollection::reset()
{
  if (reportStart("TaskCollection",getTaskName(),"reset()"))
    ;
  for (int iTask=0; iTask<nTasks; iTask++) if (isTaskOk()) tasks[iTask]->reset();
  if (reportEnd("TaskCollection",getTaskName(),"reset()"))
    ;
}

void TaskCollection::clear()
{
  if (reportStart("TaskCollection",getTaskName(),"clear()"))
    ;
  for (int iTask=0; iTask<nTasks; iTask++) if (isTaskOk()) tasks[iTask]->clear();
  if (reportEnd("TaskCollection",getTaskName(),"clear()"))
    ;
}


void TaskCollection::savePartialResults()
{
  if (reportStart("TaskCollection",getTaskName(),"savePartialResults()"))
    ;
  for (int iTask=0; iTask<nTasks; iTask++) if (isTaskOk()) tasks[iTask]->savePartialResults();
  if (reportEnd("TaskCollection",getTaskName(),"savePartialResults()()"))
    ;
}

////////////////////////////////////////////
// Print task configuration
////////////////////////////////////////////
void TaskCollection::printTaskConfiguration(ostream & output)
{
  if (reportInfo("TaskCollection",getTaskName(),"printConfiguration(ostream & output)"))
    {
    for (int iTask=0; iTask<nTasks; iTask++) if (isTaskOk()) tasks[iTask]->printTaskConfiguration(output);
    }
  if (reportEnd("TaskCollection",getTaskName(),"printConfiguration(ostream & output)"))
    ;
}

Task *  TaskCollection::addTask(Task * task)
{
  if (!task)
    {
    if (reportFatal("TaskCollection",getTaskName(),"addTask(Task * task)")) cout << "Given task pointer is null. Abort." << endl;
    postTaskFatal();
    return task;
    }

  if (nTasks>=nTasksMax)
    {
    if (reportFatal("TaskCollection",getTaskName(),"addTask(Task * task)")) cout << "Cannot add subtask " << task->getTaskName() << " to task " << getTaskName()
      << ". Maximum of " << nTasksMax << " tasks exceeded" << endl;
    postTaskFatal();
    }
  tasks[nTasks++] = task;
  if (reportInfo("TaskCollection",getTaskName(),"addTask(Task * task)")) cout << "Added task " << task->getTaskName() << " to task " << getTaskName() << endl;
  return task;
}
