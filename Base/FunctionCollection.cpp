// Author: Claude Pruneau   09/25/2019

/*************************************************************************
 * Copyright (C) 2019, Claude Pruneau.                                   *
 * All rights reserved.                                                  *
 * Based on the ROOT package and environment                             *
 *                                                                       *
 * For the licensing terms see LICENSE.                                  *
 *************************************************************************/
/**
 \class FunctionCollection
 \ingroup WAC

 Container class for root functions TF1
 */

#include "FunctionCollection.hpp"

ClassImp(FunctionCollection);

FunctionCollection::FunctionCollection(long initialCapacity)
:
Collection<TF1>(initialCapacity)
{
  /* no ops*/
}


FunctionCollection::FunctionCollection(const FunctionCollection & source)
:
Collection<TF1>(source)
{
  /* no ops*/
}

FunctionCollection::~FunctionCollection()
{
  //deletefunctions();
}

FunctionCollection & FunctionCollection::operator=(const FunctionCollection & source)
{
  Collection<TF1>::operator=(source);
  return *this;
}



////////////////////////////////////////////////////
// Set Default Style for Plots
////////////////////////////////////////////////////
void FunctionCollection::setDefaultOptions(bool color)
{
  if (color)
    gStyle->SetPalette(1,0);
  else
    gStyle->SetPalette(7,0);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  gStyle->SetOptDate(0);
  gStyle->SetOptTitle(0);
  gStyle->SetPadBorderMode(0);

}

////////////////////////////////////////////////////////////////////////
// Setting Function Properties
////////////////////////////////////////////////////////////////////////


void FunctionCollection::setFunctionProperties(TF1 * f, const GraphConfiguration & graphConfiguration)
{
  cout << "Setting properties of histo: " << f->GetTitle() << endl;
  f->SetLineColor(graphConfiguration.lineColor);
  f->SetLineStyle(graphConfiguration.lineStyle);
  f->SetLineWidth(graphConfiguration.lineWidth);
  TAxis * xAxis = (TAxis *) f->GetXaxis();
  xAxis->SetNdivisions(graphConfiguration.nXDivisions);
  xAxis->SetTitleSize(graphConfiguration.xTitleSize);
  xAxis->SetTitleOffset(graphConfiguration.xTitleOffset);
  //xAxis->SetTitle(graphConfiguration.xTitle);
  xAxis->SetLabelSize(graphConfiguration.xLabelSize);
  xAxis->SetLabelOffset(graphConfiguration.xLabelOffset);
  TAxis * yAxis = (TAxis *) f->GetYaxis();
  yAxis->SetNdivisions(graphConfiguration.nYDivisions);
  yAxis->SetTitleSize(graphConfiguration.yTitleSize);
  yAxis->SetTitleOffset(graphConfiguration.yTitleOffset);
  yAxis->SetLabelSize(graphConfiguration.yLabelSize);
  yAxis->SetLabelOffset(graphConfiguration.yLabelOffset);
  //yAxis->SetTitle(graphConfiguration.yTitle);
  if (f->IsA() == TF2::Class() )
    {
    cout << "Setting properties as 2D function: " << f->GetTitle() << endl;
    TAxis * zAxis = (TAxis *) f->GetZaxis();
    zAxis->SetNdivisions(graphConfiguration.nZDivisions);
    zAxis->SetTitleSize(graphConfiguration.zTitleSize);
    zAxis->SetTitleOffset(graphConfiguration.zTitleOffset);
    zAxis->SetLabelSize(graphConfiguration.zLabelSize);
    zAxis->SetLabelOffset(graphConfiguration.zLabelOffset);
    }
}


void FunctionCollection::plotAllFunctions(const TString & outputPath, bool doPrint)
{
  //cout << "-INFO- plotAllHistos(...) Ploting all histos of the collection." << endl;
  cout << "-INFO- plotAllFunctions(...) n functions :" << getNFunction() << endl;

  GraphConfiguration * gc1D = new GraphConfiguration(1,0);
  GraphConfiguration * gc2D = new GraphConfiguration(2,0);
  CanvasConfiguration * cc1D = new CanvasConfiguration(1);
  CanvasConfiguration * cc2D = new CanvasConfiguration(0);
  CanvasCollection * canvasCollection = new CanvasCollection();
  //cout << "-INFO- plotAllHistos(...) Creating directory (unless it already exists):" << outputPath << endl;
  canvasCollection->createDirectory(outputPath);
  TString name;
  //cout << "-INFO- plotAllHistos(...) getNFunction() :" << getNFunction()  << " with capacity:" << nFunctionCapacity << endl;

  for (int iFunc=0; iFunc<getNFunction() ; iFunc++)
    {
    TF1* f = getObjectAt(iFunc);
    name = f->GetName();
    if (f->IsA() == TF1::Class())
      {
      cout << "-INFO- plotAllHistos(...) Plotting 1D histo #" << iFunc << " named " << f->GetTitle() << endl;
      canvasCollection->createCanvas(name, *cc1D, 30);
      setFunctionProperties(f, *gc1D);
      f->Draw();
      }
    else if (f->IsA() == TF2::Class())
      {
      cout << "-INFO- plotAllHistos(...) Plotting 2D histo #" << iFunc << " named " << f->GetTitle() << endl;
      canvasCollection->createCanvas(name, *cc2D, 30);
      setFunctionProperties(f, *gc2D);
      f->Draw("SURF3");
      }
    }
  if (doPrint) canvasCollection->printAllCanvas(outputPath);
  }
