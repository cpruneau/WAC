// Author: Claude Pruneau   09/25/2019

/*************************************************************************
 * Copyright (C) 2019, Claude Pruneau.
 * All rights reserved.
 * Based on the ROOT package and environment
 *
 * For the licensing terms see LICENSE.
 *************************************************************************/
#ifndef WAC_Plotter
#define WAC_Plotter
#include "CanvasCollection.hpp"
#include "CanvasConfiguration.hpp"
#include "GraphConfiguration.hpp"
#include "MessageLogger.hpp"
#include "TH1.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TH2.h"
#include "TH2F.h"
#include "TH2D.h"
#include "TGraph.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TArrow.h"

class Plotter : public CanvasCollection, public MessageLogger
{
public:
  Plotter();
  virtual ~Plotter();

  // ================================================================================================
  // Function to plot a single 1D histogram
  // h       : pointer  to one histograms
  // legends : label  used as legend in the body of the plot
  // ================================================================================================
  TCanvas *  plot(TString   canvasName, CanvasConfiguration * cc, GraphConfiguration * gc,
                  TString   xTitle,  double xMin, double xMax,
                  TString   yTitle,  double yMin, double yMax,
                  TH1 * h,
                  TString   hLegend,
                  double xMinLeg, double yMinLeg, double xMaxLeg, double yMaxLeg,
                  double legendSize);

  // ================================================================================================
  // Function to plot a single 2D histogram
  // ================================================================================================
  TCanvas *  plot(TString   canvasName, CanvasConfiguration * cc, GraphConfiguration * gc,
                  TString   xTitle,  double xMin, double xMax,
                  TString   yTitle,  double yMin, double yMax,
                  TString   zTitle,  double zMin, double zMax,
                  TH2 * h);

  // ================================================================================================
  // Function to plot nHists 1D histogram
  // h       : array of nHists pointers to histograms
  // legends : array of nHists pointers to labels used as legends in the body of the plot
  // ================================================================================================
  TCanvas *  plot(int nHists, TString  canvasName, CanvasConfiguration * cc, GraphConfiguration ** gc,
                  TString   xTitle,  double xMin, double xMax,
                  TString   yTitle,  double yMin, double yMax,
                  TH1 ** h,
                  TString ** legends,
                  double xMinLeg, double yMinLeg, double xMaxLeg, double yMaxLeg,
                  double legendSize);

  ////////////////////////////////////////////////////////////////////////
  // Setting Histogram Properties
  ////////////////////////////////////////////////////////////////////////
  void setProperties(TH1 * h, const GraphConfiguration & graphConfiguration);
  void setProperties(TH1 * h, const GraphConfiguration & graphConfiguration, const TString & xTitle, const TString & yTitle, const TString & zTitle="");

  void setProperties(TGraph * g, const GraphConfiguration & graphConfiguration);
  void setProperties(TGraph * g, const GraphConfiguration & graphConfiguration, const TString & xTitle, const TString & yTitle);

  TLatex  * createLabel(double x, double y, int color, int fontType, double fontSize, const TString & text, bool doDraw=true);
  TLatex  * createLabel(double x, double y, double angle, int color, int fontType, double fontSize, const TString & text, bool doDraw=true);

  TLegend * createLegend(float x1, float y1, float x2, float y2, int fontType, float fontSize);
  TLegend * createLegend(TH1 * h, TString & legendText,
                         float x1, float y1, float x2, float y2, int fontType, float fontSize, bool doDraw=true);
  TLegend * createLegend(int nGraphs, TH1 ** h, TString ** legendTexts,
                         float x1, float y1, float x2, float y2, int fontType, float fontSize, bool doDraw=true);
  TLine   * createLine(float x1, float y1, float x2, float y2, int style, int color, int width, bool doDraw=true);
  TArrow  * createArrow(float x1, float y1, float x2, float y2, float arrowSize, Option_t* option, int style, int color, int width, bool doDraw=true);

    ClassDef(Plotter,0)
};

#endif /* WAC_Plotter */
