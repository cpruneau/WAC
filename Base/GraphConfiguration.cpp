// Author: Claude Pruneau   09/25/2019

/*************************************************************************
 * Copyright (C) 2019, Claude Pruneau.                                   *
 * All rights reserved.                                                  *
 * Based on the ROOT package and environment                             *
 *                                                                       *
 * For the licensing terms see LICENSE.                                  *
 *************************************************************************/
/**
 \class GenericEventAnalyzer
 \ingroup WAC

 Utility class used to define graph parameters
 */

#include "GraphConfiguration.hpp"

ClassImp(GraphConfiguration);

GraphConfiguration::GraphConfiguration(int dim, int type)
:
npx(100),
npy(100),
lineColor(1),
lineStyle(1),
lineWidth(2),
markerColor(1),
markerStyle(20),
markerSize(1.1),
nXDivisions(5),
xTitleSize(0.06),
xTitleOffset(1.1),
xLabelSize(0.05),
xLabelOffset(0.01),
xTitle("x"),
nYDivisions(5),
yTitleSize(0.06),
yTitleOffset(1.2),
yLabelSize(0.05),
yLabelOffset(0.01),
yTitle("y"),
nZDivisions(5),
zTitleSize(0.06),
zTitleOffset(1.5),
zLabelSize(0.05),
zLabelOffset(0.01),
zTitle("z"),
plotOption("E0")
{
  if (dim==1 || dim==0)
    {
    switch (type)
      {
        default:
        case  0:
        lineColor = kBlack;
        lineStyle = kSolid;
        lineWidth = 2;
        markerColor = kBlack;
        markerStyle = kFullSquare;
        markerSize  = 1.5;
        break;

        case  1:
        lineColor = kBlue;
        lineStyle = kSolid;
        lineWidth = 2;
        markerColor = kBlue;
        markerStyle = kFullCircle;
        markerSize  = 1.5;
        break;
        case  2:
        lineColor = kRed;
        lineStyle = kSolid;
        lineWidth = 2;
        markerColor = kRed;
        markerStyle = kFullSquare;
        markerSize  = 1.5;
        break;
        case  3:
        lineColor = kMagenta;
        lineStyle = kSolid;
        lineWidth = 2;
        markerColor = kMagenta;
        markerStyle = kFullDiamond;
        markerSize  = 1.5;
        break;
        case  4:
        lineColor = kBlue;
        lineStyle = kDashed;
        lineWidth = 2;
        markerColor = kBlue;
        markerStyle = kOpenSquare;
        markerSize  = 1.5;
        break;

        case  5:
        lineColor = kRed;
        lineStyle = kDashed;
        lineWidth = 1;
        markerColor = kRed;
        markerStyle = kOpenCircle;
        markerSize  = 1.5;
        break;
        case  6:
        lineColor = kBlack;
        lineStyle = kDashed;
        lineWidth = 1;
        markerColor = kBlack;
        markerStyle = kOpenSquare;
        markerSize  = 1.5;
        break;
        case  7:
        lineColor = kMagenta;
        lineStyle = kDashed;
        lineWidth = 1;
        markerColor = kMagenta;
        markerStyle = kOpenDiamond;
        markerSize  = 1.5;
        case  8:
        lineColor = kBlue;
        lineStyle = kDashed;
        lineWidth = 1;
        markerColor = kBlue;
        markerStyle = kFullDiamond;
        markerSize  = 2.1;
        break;
        case  9:
        lineColor = kBlack;
        lineStyle = 3;
        lineWidth = 1;
        markerColor = kBlack;
        markerStyle = kOpenSquare;
        markerSize  = 1.5;
        break;

        case  10:
        lineColor = kBlack;
        lineStyle = 1;
        lineWidth = 4;
        markerColor = kBlack;
        markerStyle = 20;
        markerSize  = 1.5;
        break;
        case  100:
           lineColor = kBlack;
           lineStyle = 1;
           lineWidth = 4;
           markerColor = kBlack;
           markerStyle = 20;
           markerSize  = 1.5;

           nXDivisions  = 4;
           xTitleSize   = 0.07;
           xTitleOffset = 0.7;
           xLabelSize   = 0.07;
           xLabelOffset = 0.02;
           nYDivisions  = 4;
           yTitleSize   = 0.07;
           yTitleOffset = 0.8;
           yLabelSize   = 0.07;
           yLabelOffset = 0.02;
           break;
      }
    }
  else if (dim==2)
    {
    nYDivisions  = 4;
    xTitleSize   = 0.05;
    xTitleOffset = 1.3;
    xLabelSize   = 0.05;
    xLabelOffset = 0.002;
    nYDivisions  = 4;
    yTitleSize   = 0.05;
    yTitleOffset = 1.5;
    yLabelSize   = 0.05;
    yLabelOffset = 0.002;
    nZDivisions  = 4;
    zTitleSize   = 0.05;
    zTitleOffset = 1.7;
    zLabelSize   = 0.05;
    zLabelOffset = 0.005;
    plotOption   = "SURF3";
    }
}

GraphConfiguration::GraphConfiguration(const GraphConfiguration & source)
:
npx(source.npx),
npy(source.npy),
lineColor(source.lineColor),
lineStyle(source.lineStyle),
lineWidth(source.lineWidth),
markerColor(source.markerColor),
markerStyle(source.markerStyle),
markerSize(source.markerSize),
nXDivisions(source.nXDivisions),
xTitleSize(source.xTitleSize),
xTitleOffset(source.xTitleOffset),
xLabelSize(source.xLabelSize),
xLabelOffset(source.xLabelOffset),
xTitle(source.xTitle),
nYDivisions(source.nYDivisions),
yTitleSize(source.yTitleSize),
yTitleOffset(source.yTitleOffset),
yLabelSize(source.yLabelSize),
yLabelOffset(source.yLabelOffset),
yTitle(source.yTitle),
nZDivisions(source.nZDivisions),
zTitleSize(source.zTitleSize),
zTitleOffset(source.zTitleOffset),
zLabelSize(source.zLabelSize),
zLabelOffset(source.zLabelOffset),
zTitle(source.zTitle),
plotOption(source.plotOption)
{

}

