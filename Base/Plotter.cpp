// Author: Claude Pruneau   09/25/2019

/*************************************************************************
 * Copyright (C) 2019, Claude Pruneau.                                   *
 * All rights reserved.                                                  *
 * Based on the ROOT package and environment                             *
 *                                                                       *
 * For the licensing terms see LICENSE.                                  *
 *************************************************************************/
/**
 \class Task
 \ingroup WAC

 Class defining Event Plotter
 */

#include "Plotter.hpp"
ClassImp(Plotter);


Plotter::Plotter()
:
CanvasCollection(),
MessageLogger()
{
  if (reportDebug()) cout << "Plotter::CTOR() No OPS." << endl;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
// DTOR
////////////////////////////////////////////////////////////////////////////////////////////////////////
Plotter::~Plotter()
{
  if (reportDebug()) cout << "Plotter::DTOR() No ops." << endl;
}

TCanvas *  Plotter::plot(TString  canvasName, CanvasConfiguration * cc, GraphConfiguration * gc,
                          TString  xTitle,  double xMin, double xMax,
                          TString  yTitle,  double yMin, double yMax,
                          TH1 * h,
                          TString  legendText,
                          double xMinLeg, double yMinLeg, double xMaxLeg, double yMaxLeg,
                          double legendSize)
{
  TCanvas * canvas  = createCanvas(canvasName,*cc);
  setProperties(h,*gc);
  h->GetXaxis()->SetTitle(xTitle);
  h->GetYaxis()->SetTitle(yTitle);
  h->SetMinimum(yMin);
  h->SetMaximum(yMax);
  h->GetXaxis()->SetRangeUser(xMin,xMax);
  h->DrawCopy(gc->plotOption);
  createLegend(h,legendText,xMinLeg, yMinLeg, xMaxLeg, yMaxLeg,0, legendSize);
  return canvas;
}

TCanvas *  Plotter::plot(TString   canvasName, CanvasConfiguration * cc, GraphConfiguration * gc,
                         TString   xTitle,  double xMin, double xMax,
                         TString   yTitle,  double yMin, double yMax,
                         TString   zTitle,  double zMin, double zMax,
                         TH2 * h)
{
  TCanvas * canvas  = createCanvas(canvasName,*cc);
  setProperties(h,*gc);
  h->GetXaxis()->SetTitle(xTitle);
  h->GetYaxis()->SetTitle(yTitle);
  h->GetZaxis()->SetTitle(zTitle);

  if (xMin < xMax)
    {
    //cout << "Setting x range: " << xMin << " to: " << xMax << endl;
    h->GetXaxis()->SetRangeUser(xMin,xMax);
    }
  if (yMin < yMax)
    {
    //cout << "Setting y range: " << yMin << " to: " << yMax << endl;
    h->GetYaxis()->SetRangeUser(yMin,yMax);
    }

  if (zMin < zMax)
    {
    //cout << "Setting z range: " << zMin << " to: " << zMax << endl;
    h->SetMinimum(zMin);
    h->SetMaximum(zMax);
    }
  h->DrawCopy(gc->plotOption);
  return canvas;
}


TCanvas *  Plotter::plot(int nGraphs, TString  canvasName, CanvasConfiguration * cc, GraphConfiguration ** gc,
                          TString  xTitle,  double xMin, double xMax,
                          TString  yTitle,  double yMin, double yMax,
                          TH1 ** h,
                          TString ** legendTexts,
                          double xMinLeg, double yMinLeg, double xMaxLeg, double yMaxLeg,
                          double legendSize
                          )
{
  TCanvas * canvas = createCanvas(canvasName,*cc);
  setProperties(h[0],*gc[0]);
  h[0]->GetXaxis()->SetTitle(xTitle);
  h[0]->GetYaxis()->SetTitle(yTitle);
  double min = -1.0;
  double max =  1.0;
  if (yMin < yMax)
    {
    min = yMin;
    max = yMax;
    }
  else if (yMin >= yMax)
    {
    min =  1.0E100;
    max = -1.0E100;
    //cout << "Histo: " << h[0]->GetTitle() << endl;
    for (int iGraph=0; iGraph<nGraphs; iGraph++)
      {
      double vMax = h[iGraph]->GetBinContent(h[iGraph]->GetMaximumBin());
      if (vMax>max){max = vMax;}
      double vMin = h[iGraph]->GetBinContent(h[iGraph]->GetMinimumBin());
      if (vMin<min){min = vMin;}
      }
    //cout << "Found min : " << min << " max : " << max << endl;
    if (max>min)
      {
      if (yMax<0.0)
        {
        if (max>0)
          {
          max  = fabs(yMax)*max;
          }
        else if (max<0)
          {
          max  = max/fabs(yMax);
          }
        }
      if (yMin>0.0)
        {
        if (min<0)
          {
          min = fabs(yMin)*min;
          }
        else if (min>0)
          {
          min = min/fabs(yMin);
          }
        }
      }
    else
      {
      min = -1.0;
      max =  1.0;
      }
    }
 // cout << "Setting min to : " << min << " max to: " << max << endl;
  h[0]->SetMinimum(min);
  h[0]->SetMaximum(max);


  if (xMin < xMax)
    {
    h[0]->GetXaxis()->SetRangeUser(xMin,xMax);
    }

  TString plotOption;
  setProperties(h[0],*gc[0]);
  plotOption = gc[0]->plotOption;
  h[0]->DrawCopy(plotOption);
  for (int iGraph=1; iGraph<nGraphs; iGraph++)
    {
    setProperties(h[iGraph],*gc[iGraph]);
    plotOption = gc[iGraph]->plotOption;
    h[iGraph]->DrawCopy(plotOption+" SAME");
    }
  if (nGraphs<6)
    createLegend(nGraphs,h,legendTexts,xMinLeg, yMinLeg, xMaxLeg, yMaxLeg,0, legendSize);
  else
    {
    int n1 = nGraphs/2;
    int n2 = nGraphs - n1;
    createLegend(n1,h,legendTexts,xMinLeg, yMinLeg, xMaxLeg, yMaxLeg,0, legendSize);
    createLegend(n2,h+n1,legendTexts+n1,xMinLeg+0.2, yMinLeg, xMaxLeg+0.2, yMaxLeg,0, legendSize);
    }
  return canvas;
}


// =======================================
// Create a label within an existing pad
// =======================================
TLatex * Plotter::createLabel(double x, double y, int color, int fontType, double fontSize, const TString & text, bool doDraw)
{
  TLatex * label;
  label = new TLatex(x,y,text);
  label->SetTextColor(color);
  fontType=0; // to avoid warnings.
              //label->SetTextFont(fontType);
  label->SetTextSize(fontSize);
  if (doDraw) label->Draw();
  return label;
}

TLatex * Plotter::createLabel(double x, double y, double angle, int color, int fontType, double fontSize, const TString & text, bool doDraw)
{
  TLatex * label;
  label = new TLatex(x,y,text);
  label->SetTextColor(color);
  label->SetTextAngle(angle);
  fontType=0; // to avoid warnings.
              //label->SetTextFont(fontType);
  label->SetTextSize(fontSize);
  if (doDraw) label->Draw();
  return label;
}


////////////////////////////////////////////////////
// Create standard legend
////////////////////////////////////////////////////
TLegend * Plotter::createLegend(float x1, float y1, float x2, float y2, int fontType, float fontSize)
{
  TLegend *legend = new TLegend(x1,y1,x2,y2);
  //legend->SetTextFont(fontType);
  fontType=0; // to avoid warnings.
  legend->SetTextSize(fontSize);
  legend->SetFillColor(0);
  legend->SetBorderSize(0);
  return legend;
}

TLegend * Plotter::createLegend(TH1 * h, TString & legendText,
                                 float x1, float y1, float x2, float y2, int fontType, float fontSize, bool doDraw)
{
  TLegend *legend = createLegend(x1,y1,x2,y2,fontType,fontSize);
  legend->AddEntry(h,legendText);
  if (doDraw) legend->Draw();
  return legend;
}

TLegend * Plotter::createLegend(int nGraphs, TH1 ** h, TString ** legendTexts,
                                 float x1, float y1, float x2, float y2, int fontType, float fontSize, bool doDraw)
{
  TLegend *legend = createLegend(x1,y1,x2,y2,fontType,fontSize);
  for (int iGraph=0; iGraph<nGraphs; iGraph++)
    {
    legend->AddEntry(h[iGraph],*legendTexts[iGraph]);
    }
  if (doDraw) legend->Draw();
  return legend;
}

////////////////////////////////////////////////////
// Create simple line
////////////////////////////////////////////////////
TLine * Plotter::createLine(float x1, float y1, float x2, float y2, int style, int color, int width, bool doDraw)
{
  TLine *line = new TLine(x1,y1,x2,y2);
  line->SetLineStyle(style);
  line->SetLineColor(color);
  line->SetLineWidth(width);
  if (doDraw) line->Draw();
  return line;
}

////////////////////////////////////////////////////
// Create Arrow Line
////////////////////////////////////////////////////
TArrow * Plotter::createArrow(float x1, float y1, float x2, float y2, float arrowSize, Option_t* option, int style, int color, int width, bool doDraw)
{
  TArrow *line = new TArrow(x1,y1,x2,y2,arrowSize,option);
  line->SetLineStyle(style);
  line->SetLineColor(color);
  line->SetLineWidth(width);
  if (doDraw) line->Draw();
  return line;
}


////////////////////////////////////////////////////////////////////////
// Setting Histogram Properties
////////////////////////////////////////////////////////////////////////
void Plotter::setProperties(TH1 * h, const GraphConfiguration & graphConfiguration)
{
  if (reportDebug()) cout << "Setting properties of histo: " << h->GetTitle() << endl;
  h->SetLineColor(graphConfiguration.lineColor);
  h->SetLineStyle(graphConfiguration.lineStyle);
  h->SetLineWidth(graphConfiguration.lineWidth);
  h->SetMarkerColor(graphConfiguration.markerColor);
  h->SetMarkerStyle(graphConfiguration.markerStyle);
  h->SetMarkerSize (graphConfiguration.markerSize);
  TAxis * xAxis = (TAxis *) h->GetXaxis();
  xAxis->SetNdivisions(graphConfiguration.nXDivisions);
  xAxis->SetTitleSize(graphConfiguration.xTitleSize);
  xAxis->SetTitleOffset(graphConfiguration.xTitleOffset);
  xAxis->SetLabelSize(graphConfiguration.xLabelSize);
  xAxis->SetLabelOffset(graphConfiguration.xLabelOffset);
  TAxis * yAxis = (TAxis *) h->GetYaxis();
  yAxis->SetNdivisions(graphConfiguration.nYDivisions);
  yAxis->SetTitleSize(graphConfiguration.yTitleSize);
  yAxis->SetTitleOffset(graphConfiguration.yTitleOffset);
  yAxis->SetLabelSize(graphConfiguration.yLabelSize);
  yAxis->SetLabelOffset(graphConfiguration.yLabelOffset);
  if (h->IsA() == TH2::Class()  || h->IsA() == TH2F::Class() || h->IsA() == TH2F::Class() )
    {
    if (reportDebug()) cout << "Setting properties as 2D histo: " << h->GetTitle() << endl;
    TAxis * zAxis = (TAxis *) h->GetZaxis();
    zAxis->SetNdivisions(graphConfiguration.nZDivisions);
    zAxis->SetTitleSize(graphConfiguration.zTitleSize);
    zAxis->SetTitleOffset(graphConfiguration.zTitleOffset);
    zAxis->SetLabelSize(graphConfiguration.zLabelSize);
    zAxis->SetLabelOffset(graphConfiguration.zLabelOffset);
    }
}

void Plotter::setProperties(TH1 * h, const GraphConfiguration & graphConfiguration,
                            const TString & xTitle, const TString & yTitle, const TString & zTitle)
{
  setProperties(h, graphConfiguration);
  TAxis * xAxis = (TAxis *) h->GetXaxis();
  xAxis->SetTitle(xTitle);
  TAxis * yAxis = (TAxis *) h->GetYaxis();
  yAxis->SetTitle(yTitle);
  if (h->IsA() == TH2::Class()  || h->IsA() == TH2F::Class() || h->IsA() == TH2F::Class() )
    {
    TAxis * zAxis = (TAxis *) h->GetZaxis();
    zAxis->SetTitle(zTitle);
    }
}

void Plotter::setProperties(TGraph * g, const GraphConfiguration & graphConfiguration)
{
  if (reportDebug()) cout << "Setting properties of histo: " << g->GetTitle() << endl;
  g->SetLineColor(graphConfiguration.lineColor);
  g->SetLineStyle(graphConfiguration.lineStyle);
  g->SetLineWidth(graphConfiguration.lineWidth);
  g->SetMarkerColor(graphConfiguration.markerColor);
  g->SetMarkerStyle(graphConfiguration.markerStyle);
  g->SetMarkerSize (graphConfiguration.markerSize);
  TAxis * xAxis = (TAxis *) g->GetXaxis();
  xAxis->SetNdivisions(graphConfiguration.nXDivisions);
  xAxis->SetTitleSize(graphConfiguration.xTitleSize);
  xAxis->SetTitleOffset(graphConfiguration.xTitleOffset);
  xAxis->SetLabelSize(graphConfiguration.xLabelSize);
  xAxis->SetLabelOffset(graphConfiguration.xLabelOffset);
  TAxis * yAxis = (TAxis *) g->GetYaxis();
  yAxis->SetNdivisions(graphConfiguration.nYDivisions);
  yAxis->SetTitleSize(graphConfiguration.yTitleSize);
  yAxis->SetTitleOffset(graphConfiguration.yTitleOffset);
  yAxis->SetLabelSize(graphConfiguration.yLabelSize);
  yAxis->SetLabelOffset(graphConfiguration.yLabelOffset);
}


 void Plotter::setProperties(TGraph * g, const GraphConfiguration & graphConfiguration, const TString & xTitle, const TString & yTitle)
{
  setProperties(g,graphConfiguration);
  if (reportDebug()) cout << "Setting properties of histo: " << g->GetTitle() << endl;
  TAxis * xAxis = (TAxis *) g->GetXaxis();
  xAxis->SetTitle(xTitle);
  TAxis * yAxis = (TAxis *) g->GetYaxis();
  yAxis->SetTitle(yTitle);
}


// ======================================================================
// Set Default Style for Plots
// ======================================================================
void Plotter::setDefaultOptions(bool color)
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
