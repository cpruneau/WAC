// Author: Claude Pruneau   09/25/2019

/*************************************************************************
 * Copyright (C) 2019, Claude Pruneau.                                   *
 * All rights reserved.                                                  *
 * Based on the ROOT package and environment                             *
 *                                                                       *
 * For the licensing terms see LICENSE.                                  *
 *************************************************************************/
/**
 \class EccentricityCalc
 \ingroup WAC

 Class defining EccentricityCalc
 */

#include "EccentricityCalc.hpp"

ClassImp(EccentricityCalc);

///////////////////////////////////////////////////////////////////////////////////////////////////
// Constructor
///////////////////////////////////////////////////////////////////////////////////////////////////
EccentricityCalc::EccentricityCalc(const TString & _systemName,
                                   int nX, double minX, double maxX,
                                   int nY, double minY, double maxY,
                                   int nEcc, double minEcc, double maxEcc,
                                   bool _doHistograms=false)
:
systemName(_systemName),
doHistograms(_doHistograms)
{
  cout << "-INFO- EccentricityCalc::EccentricityCalc() ctor Started" << endl;
  if (doHistograms)
    {
    TString baseName = systemName;
    baseName += "_";
    h_x         = new TH1D(baseName+"x",        baseName+"x",nX,minX,maxX);
    h_y         = new TH1D(baseName+"y",        baseName+"y",nY,minY,maxY);
    h2_yVsX     = new TH2D(baseName+"yVsX",     baseName+"yVsX",nX,minX,maxX,nY,minY,maxY);

    h_avg_x     = new TH1D(baseName+"avg_x",    baseName+"avg_x",nX,minX,maxX);
    h_avg_y     = new TH1D(baseName+"avg_y",    baseName+"avg_x",nY,minY,maxY);
    h_var_x     = new TH1D(baseName+"var_x",    baseName+"var_x", nX,0.0,maxX*maxX);
    h_var_y     = new TH1D(baseName+"var_y",    baseName+"var_y", nY,0.0,maxY*maxY);
    h_var_xy    = new TH1D(baseName+"var_xy",   baseName+"var_xy",nX,0.0,maxX*maxY);
    h_eps_denom = new TH1D(baseName+"eps_denom",baseName+"eps_denom",nX,0.0,maxX*maxX);
    h_eps_x     = new TH1D(baseName+"eps_x",    baseName+"eps_x",    nEcc,minEcc, maxEcc);
    h_eps_y     = new TH1D(baseName+"eps_y",    baseName+"eps_y",    nEcc,minEcc, maxEcc);
    h_eps_mod   = new TH1D(baseName+"eps_mod",  baseName+"eps_mod",  nEcc, 0.0, maxEcc);
    h2_eps_xy   = new TH2D(baseName+"eps_xy",   baseName+"eps_xy",   nEcc,minEcc, maxEcc,nEcc,minEcc, maxEcc);
    }
  clear();
  cout << "-INFO- EccentricityCalc::EccentricityCalc() Completed" << endl;
}

void EccentricityCalc::fill(TH2 * h)
{
  int nX = h->GetXaxis()->GetNbins();
  double minX = h->GetXaxis()->GetXmin();
  double maxX = h->GetXaxis()->GetXmax();
  int nY = h->GetYaxis()->GetNbins();
  double minY = h->GetYaxis()->GetXmin();
  double maxY = h->GetYaxis()->GetXmax();
  double x, y, w;
  TAxis * xAxis = h->GetXaxis();
  TAxis * yAxis = h->GetYaxis();

  for (int iX=1; iX<=nX; iX++)
    {
    x = xAxis->GetBinCenter(iX);
    for (int iY=0; iY<nY; iY++)
      {
      y = yAxis->GetBinCenter(iY);
      w = h->GetBinContent(iX,iY);
      fill(x,y,w);
      }
    }

}


// b is the impact parameter in fm.
void EccentricityCalc::fillWithNuclei(double b, GeometryGenerator * nucleus1, GeometryGenerator * nucleus2)
{
  double xAvg, yAvg, dx, dy, dSq;

  for (int i1=0; i1<nucleus1->nucleusSize; i1++)
    {
    for (int i2=0; i2<nucleus2->nucleusSize; i2++)
       {
       dx = nucleus2->points[i2]->x - nucleus1->points[i1]->x + b;
       dy = nucleus2->points[i2]->y - nucleus1->points[i1]->y;
       dSq = dx*dx + dy*dy;
       if (dSq<nucleus1->maxDistanceSq)
         {
         // this is a collision
         xAvg = (nucleus1->points[i1]->x + nucleus2->points[i2]->x)/2.0;
         yAvg = (nucleus1->points[i1]->y + nucleus2->points[i2]->y)/2.0;
         fill(xAvg,yAvg,1.0);
         }
       }
    }
}

void EccentricityCalc::fill(double x, double y, double w=1.0)
{
  double x2 = x*x;
  double x3 = x2*x;
  double x4 = x3*x;
  double x6 = x3*x3;
  double x8 = x4*x4;
  double y2 = y*y;
  double y3 = y2*y;
  double y4 = y3*y;
  double y6 = y3*y3;
  double y8 = y4*y4;
  double xy = x*y;
  counts  += w;
  mom_x   += w*x;
  mom_y   += w*y;
  mom_x2  += w*x2;
  mom_y2  += w*y2;
  mom_x3  += w*x3;
  mom_y3  += w*y3;
  mom_x4  += w*x4;
  mom_y4  += w*y4;
  mom_x6  += w*x6;
  mom_y6  += w*y6;
  mom_x8  += w*x8;
  mom_y8  += w*y8;
  mom_xy  += w*xy;

  if (doHistograms)
    {
    // cout << "x:" << x << " y:" << y << " w:" << w << endl;
    h_x->Fill(x,w);
    h_y->Fill(y,w);
    h2_yVsX->Fill(x,y,w);
    }

}

/// Fill average histograms -- call after calculateAverages()
void EccentricityCalc::fillAvg()
{
  if (doHistograms && counts>0)
    {
    h_avg_x->Fill(avg_x);
    h_avg_y->Fill(avg_y);
    h_var_x->Fill(var_x);
    h_var_y->Fill(var_y);
    h_var_xy->Fill(var_xy);
    h_eps_denom->Fill(eps_denom);
    h_eps_x->Fill(eps_x);
    h_eps_y->Fill(eps_y);
    h_eps_mod->Fill(eps_mod);
    h2_eps_xy->Fill(eps_x,eps_y);
    }
}


/// Clear container to fresh start
void EccentricityCalc::clear()
{
  if (doHistograms)
    {
    h_x->Reset();
    h_y->Reset();
    h2_yVsX->Reset();
    h_avg_x->Reset();
    h_avg_y->Reset();
    h_var_x->Reset();
    h_var_y->Reset();
    h_var_xy->Reset();
    h_eps_denom->Reset();
    h_eps_x->Reset();
    h_eps_y->Reset();
    h_eps_mod->Reset();
    h2_eps_xy->Reset();
    }
  reset();
}

/// Reset event to zero but not the histograms.
void EccentricityCalc::reset()
{
  counts  = 0;
  mom_x   = 0;
  mom_y   = 0;
  mom_x2  = 0;
  mom_y2  = 0;
  mom_x3  = 0;
  mom_y3  = 0;
  mom_x4  = 0;
  mom_y4  = 0;
  mom_xy  = 0;

  avg_x  = 0;
  avg_y  = 0;
  avg_x2  = 0;
  avg_y2  = 0;
  avg_x3  = 0;
  avg_y3  = 0;
  avg_x4  = 0;
  avg_y4  = 0;
  avg_xy  = 0;

  var_x   = 0;
  var_y   = 0;
  var_xy  = 0;

  eps_denom = 0;
  eps_x     = 0;
  eps_y     = 0;
  eps_mod   = 0;
}

/// Save histograms.
 


void EccentricityCalc::calculateAverages()
{
  //cout << "-INFO- EccentricityCalc::calculateAverages() Started" << endl;

  if (counts<1.0E-10)
    {
    avg_x    = 0.0;
    avg_y    = 0.0;
    avg_x2   = 0.0;
    avg_y2   = 0.0;
    avg_x3   = 0.0;
    avg_y3   = 0.0;
    avg_x4   = 0.0;
    avg_y4   = 0.0;
    avg_x6   = 0.0;
    avg_y6   = 0.0;
    avg_x8   = 0.0;
    avg_y8   = 0.0;
    avg_xy   = 0.0;

    var_x    = 0.0;
    var_y    = 0.0;
    var_xy   = 0.0;

    eps_denom = 0.0;
    eps_x     = 0.0;
    eps_y     = 0.0;
    eps_mod   = 0.0;

    eps_22   = 0.0;
    eps_24   = 0.0;
    eps_26   = 0.0;
    eps_28   = 0.0;

    psi2     = 0.0;
    }
  else
    {
    avg_x   = mom_x/counts;
    avg_y   = mom_y/counts;
    avg_x2   = mom_x2/counts;
    avg_y2   = mom_y2/counts;
    avg_x3   = mom_x3/counts;
    avg_y3   = mom_y3/counts;
    avg_x4   = mom_x4/counts;
    avg_y4   = mom_y4/counts;
    avg_x6   = mom_x6/counts;
    avg_y6   = mom_y6/counts;
    avg_x8   = mom_x8/counts;
    avg_y8   = mom_y8/counts;
    avg_xy   = mom_xy/counts;

    var_x  = avg_x2 - avg_x*avg_x;   //cout << " var_x:" << var_x << endl;
    var_y  = avg_y2 - avg_x*avg_y;   //cout << " var_y:" << var_y << endl;
    var_xy = avg_xy - avg_x*avg_y;   //cout << " var_xy:" << var_xy << endl;

    eps_denom = var_x + var_y;      //cout << " eps_denom:" << eps_denom << endl;
    eps_x     = (var_y - var_x)/eps_denom;      //cout << " eps_x:" << eps_x << endl;
    eps_y     = 2*var_xy/eps_denom;             //cout << " eps_y:" << eps_y << endl;
    eps_mod   = sqrt(eps_x*eps_x + eps_y*eps_y); //cout << " eps_mod:" << eps_mod << endl;

    eps_22   = eps_mod*eps_mod;
    eps_24   = eps_22*eps_22;
    eps_26   = eps_22*eps_24;
    eps_28   = eps_24*eps_24;

    psi2 = atan2(eps_y,eps_x);
    }

  //cout << "-INFO- EccentricityCalc::calculateAverages() Completed" << endl;
}


