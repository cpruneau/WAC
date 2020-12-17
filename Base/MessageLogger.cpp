// Author: Claude Pruneau   09/25/2019

/***********************************************************************
 * Copyright (C) 2019, Claude Pruneau.
 * All rights reserved.
 * Based on the ROOT package and environment
 *
 * For the licensing terms see LICENSE.
 **********************************************************************/
/**
 \class MessageLogger
 \ingroup WAC
 */

#include "MessageLogger.hpp"

ClassImp(MessageLogger);


bool MessageLogger::reportDebug(TString className, TString taskName, TString fctName, ostream & output) const
{
 if (reportLevel<=Debug)
   {
   output << " <DEBUG> " << className << "[" << taskName << "]::" << fctName << "  ";
   return true;
   }
 else
   return false;
}

bool MessageLogger::reportNoOps(TString  className, TString taskName, TString fctName, ostream & output) const
{
 if (reportLevel<=Debug)
   {
   output << " <DEBUG> " << className << "[" << taskName << "]::" << fctName << "  No ops." << endl;
   return true;
   }
 else
   return false;
}

bool MessageLogger::reportStart(TString  className, TString taskName, TString fctName, ostream & output) const
{
 if (reportLevel<=Debug)
   {
   output << " <DEBUG> " << className << "[" << taskName << "]::" << fctName << "  Started." << endl;
   return true;
   }
 else
   return false;
}

bool MessageLogger::reportEnd(TString  className, TString taskName, TString fctName, ostream & output) const
{
 if (reportLevel<=Debug)
   {
   output << " <DEBUG> " << className << "[" << taskName << "]::" << fctName << "  Completed." << endl;
   return true;
   }
 else
   return false;
}


bool MessageLogger::reportInfo(TString  className, TString taskName, TString fctName, ostream & output) const
{
 if (reportLevel<=Info)
   {
   output << " <INFO> " << className << "[" << taskName << "]::" << fctName << "  ";
   return true;
   }
 else
   return false;
}


bool MessageLogger::reportWarning(TString  className, TString taskName, TString fctName, ostream & output) const
{
 if (reportLevel<=Warning)
   {
   output << " <WARNING> " << className << "(" << taskName << ")::" << fctName << "  ";
   return true;
   }
 else
   return false;
}


bool MessageLogger::reportError(TString  className, TString taskName, TString fctName, ostream & output) const
{
 if (reportLevel<=Error)
   {
   output << " <ERROR> " << className << "(" << taskName << ")::" << fctName << "  ";
   return true;
   }
 else
   return false;
}

bool MessageLogger::reportFatal(TString  className, TString taskName, TString fctName, ostream & output) const
{
 if (reportLevel<=Fatal)
   {
   output << " <FATAL> " << className << "(" << taskName << ")::" << fctName << "  ";
   return true;
   }
 else
   return false;
}



bool MessageLogger::reportDebug(ostream & output) const
{
  if (reportLevel<=Debug)
    {
    output << " <DEBUG> ";
    return true;
    }
  else
    return false;
}

bool MessageLogger::reportInfo(ostream & output) const
{
  if (reportLevel<=Info)
    {
    output << " <INFO> ";
    return true;
    }
  else
    return false;
}

bool MessageLogger::reportWarning(ostream & output) const
{
  if (reportLevel<=Warning)
    {
    output << " <WARNING> ";
    return true;
    }
  else
    return false;
}

bool MessageLogger::reportError(ostream & output) const
{
  if (reportLevel<=Error)
    {
    output << " <ERROR> ";
    return true;
    }
  else
    return false;
}

bool MessageLogger::reportFatal(ostream & output) const
{
  if (reportLevel<=Fatal)
    {
    output << " <FATAL> ";
    return true;
    }
  else
    return false;
}
