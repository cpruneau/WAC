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
  
  
  ////////////////////////////////////////////////////////////////////////////////////////////////////////
  // CTOR
  ////////////////////////////////////////////////////////////////////////////////////////////////////////
  Plotter();
  virtual ~Plotter();

  TCanvas *  plot(TString   canvasName, CanvasConfiguration * cc, GraphConfiguration * gc,
                  TString   xTitle,  double xMin, double xMax,
                  TString   yTitle,  double yMin, double yMax,
                  TH1 * h,
                  TString   hLegend,
                  double xMinLeg, double yMinLeg, double xMaxLeg, double yMaxLeg,
                  double legendSize);

  TCanvas *  plot(int nGraphs, TString  canvasName, CanvasConfiguration * cc, GraphConfiguration ** gc,
                   TString   xTitle,  double xMin, double xMax,
                   TString   yTitle,  double yMin, double yMax,
                   TH1 ** h,
                   TString ** hLegend,
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
