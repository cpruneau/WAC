// Author: Claude Pruneau   09/25/2019

/***********************************************************************
 * Copyright (C) 2019, Claude Pruneau.
 * All rights reserved.
 * Based on the ROOT package and environment
 *
 * For the licensing terms see LICENSE.
 **********************************************************************/

#ifndef WAC_MessageLogger
#define WAC_MessageLogger
#include "TClass.h"
#include <iostream>
using namespace std;
using std::cout;
using std::endl;

class MessageLogger
{
public:

  enum LogLevel   { Unknown, Debug, Info, Warning, Error, Fatal};

  LogLevel  reportLevel;

  MessageLogger(LogLevel selectedLevel=Info) : reportLevel( selectedLevel )
  {   }

  virtual ~MessageLogger()
  { }

  void setReportLevel(LogLevel requiredLevel)
  {
    reportLevel = requiredLevel;
  }

  void setLogLevel(LogLevel requiredLevel)
  {
    reportLevel = requiredLevel;
  }

  LogLevel getLogLevel() const
  {
    return reportLevel;
  }

  LogLevel getReportLevel() const
  {
    return reportLevel;
  }

  bool reportDebug(ostream & output=cout) const;
  bool reportInfo(ostream & output=cout) const;
  bool reportWarning(ostream & output=cout) const;
  bool reportError(ostream & output=cout) const;
  bool reportFatal(ostream & output=cout) const;

  bool reportDebug(TString  className, TString  fctName, TString  taskName, ostream & output=cout) const;
  bool reportNoOps(TString  className, TString  fctName, TString  taskName, ostream & output=cout) const;
  bool reportStart(TString  className, TString  fctName, TString  taskName, ostream & output=cout) const;
  bool reportEnd(TString  className, TString  fctName, TString  taskName, ostream & output=cout) const;
  bool reportInfo(TString  className, TString  fctName, TString  taskName, ostream & output=cout) const;
  bool reportWarning(TString  className, TString  fctName, TString  taskName, ostream & output=cout) const;
  bool reportError(TString  className, TString  fctName, TString  taskName, ostream & output=cout) const;
  bool reportFatal(TString  className, TString  fctName, TString  taskName, ostream & output=cout) const;


  ClassDef(MessageLogger,0)
};

#endif /* WAC_MessageLogger */
