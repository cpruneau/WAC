// Author: Claude Pruneau   09/25/2019

/*************************************************************************
 * Copyright (C) 2019, Claude Pruneau.                                   *
 * All rights reserved.                                                  *
 * Based on the ROOT package and environment                             *
 *                                                                       *
 * For the licensing terms see LICENSE.                                  *
 *************************************************************************/
/**
 \class CanvasConfiguration
 \ingroup WAC

 Utility class used to define the parameters of a root canvas
 */

#include "CanvasCollection.hpp"

ClassImp(CanvasCollection);

CanvasConfiguration::CanvasConfiguration(CanvasFormat format, CanvasAxes axes)
:
x(20),
y(20),
width(700),
height(700),
leftMargin(0.15),
topMargin(0.03),
rightMargin(0.03),
bottomMargin(0.15),
logx(0),
logy(0),
logz(0),
theta(35.0),
phi(25.0)
{
  // chose the aspect ratio and the margin style (tight, regular, wide)
  setFormat(format);
  // chose the axes settings
  setAxes(axes);
}

// deprecated...
CanvasConfiguration::CanvasConfiguration(bool log, int style, int dimension)
:
x(20),
y(20),
width(700),
height(700),
leftMargin(0.15),
topMargin(0.03),
rightMargin(0.03),
bottomMargin(0.15),
logx(0),
logy(0),
logz(0),
theta(35.0),
phi(25.0)
{
  switch (dimension)
    {
      default:
      case 0: break;
      case 1:
      switch (style)
        {
          default:
          case 0:
          leftMargin   = 0.15;
          topMargin    = 0.03;
          rightMargin  = 0.14;
          bottomMargin = 0.15;
          break;
          case 1:
          leftMargin   = 0.15;
          topMargin    = 0.03;
          rightMargin  = 0.05;
          bottomMargin = 0.15;
          break;
        }
      if (log)
        {
        logy = 1;
        }
      break;
      case 2:

      switch (style)
        {
          default:
          case 0:
          leftMargin   = 0.15;
          topMargin    = 0.03;
          rightMargin  = 0.14;
          bottomMargin = 0.15;
          break;
          case 1:
          leftMargin   = 0.15;
          topMargin    = 0.03;
          rightMargin  = 0.05;
          bottomMargin = 0.15;
          break;
        }
      if (log)
        {
        logz = 1;
        }
      break;
    }
}

////////////////////////////////////////////////////
// COPY CTOR
////////////////////////////////////////////////////
CanvasConfiguration::CanvasConfiguration(const CanvasConfiguration & canvasConfig)
:
x(canvasConfig.x),
y(canvasConfig.y),
width(canvasConfig.width),
height(canvasConfig.height),
leftMargin(canvasConfig.leftMargin),
topMargin(canvasConfig.topMargin),
rightMargin(canvasConfig.rightMargin),
bottomMargin(canvasConfig.bottomMargin),
logx(canvasConfig.logx),
logy(canvasConfig.logy),
logz(canvasConfig.logz),
theta(canvasConfig.theta),
phi(canvasConfig.phi)
{}

////////////////////////////////////////////////////
// Assignment operator
////////////////////////////////////////////////////
CanvasConfiguration & CanvasConfiguration::operator=(const CanvasConfiguration & canvasConfig)
{
  x = canvasConfig.x;
  y = canvasConfig.y;
  width  = canvasConfig.width;
  height = canvasConfig.height;
  leftMargin = canvasConfig.leftMargin;
  topMargin  = canvasConfig.topMargin;
  rightMargin = canvasConfig.rightMargin;
  bottomMargin = canvasConfig.bottomMargin;
  logx  = canvasConfig.logx;
  logy  = canvasConfig.logy;
  logz  = canvasConfig.logz;
  theta = canvasConfig.theta;
  phi   = canvasConfig.phi;

  return *this;
}

CanvasConfiguration::~CanvasConfiguration()
{

}

void CanvasConfiguration::setFormat(CanvasFormat format)
{
  switch (format)
    {
      case PortraitTight:
      width    = 500;
      height   = 800;
      leftMargin   = 0.10;
      topMargin    = 0.03;
      rightMargin  = 0.03;
      bottomMargin = 0.10;
      break;
      case Portrait:
      width    = 500;
      height   = 800;
      leftMargin   = 0.15;
      topMargin    = 0.05;
      rightMargin  = 0.05;
      bottomMargin = 0.15;
      break;
      case PortraitWide:
      width    = 500;
      height   = 800;
      leftMargin   = 0.2;
      topMargin    = 0.1;
      rightMargin  = 0.1;
      bottomMargin = 0.2;
      break;
      case SquareTight:
      width    = 700;
      height   = 700;
      leftMargin   = 0.10;
      topMargin    = 0.03;
      rightMargin  = 0.03;
      bottomMargin = 0.10;
      break;
      case Square:
      width    = 700;
      height   = 700;
      leftMargin   = 0.15;
      topMargin    = 0.05;
      rightMargin  = 0.05;
      bottomMargin = 0.15;
      break;
      case SquareWide:
      width    = 700;
      height   = 700;
      leftMargin   = 0.2;
      topMargin    = 0.1;
      rightMargin  = 0.1;
      bottomMargin = 0.2;
      break;
      case LandscapeTight:
      width    = 800;
      height   = 500;
      leftMargin   = 0.10;
      topMargin    = 0.03;
      rightMargin  = 0.03;
      bottomMargin = 0.10;
      break;
      case Landscape:
      width    = 800;
      height   = 500;
      leftMargin   = 0.15;
      topMargin    = 0.05;
      rightMargin  = 0.1;
      bottomMargin = 0.15;
      break;
      case LandscapeWide:
      width    = 800;
      height   = 500;
      leftMargin   = 0.2;
      topMargin    = 0.1;
      rightMargin  = 0.15;
      bottomMargin = 0.2;
      break;
    };
}


void CanvasConfiguration::setAxes(CanvasAxes axes)
{
  switch (axes)
    {
      case Linear:    logx = 0;  logy = 0;  logz = 0; break;
      case LogX:      logx = 1;  logy = 0;  logz = 0; break;
      case LogY:      logx = 0;  logy = 1;  logz = 0; break;
      case LogZ:      logx = 0;  logy = 0;  logz = 1; break;
      case LogXY:     logx = 1;  logy = 1;  logz = 0; break;
      case LogXZ:     logx = 1;  logy = 0;  logz = 1; break;
      case LogYZ:     logx = 0;  logy = 1;  logz = 1; break;
      case LogXYZ:    logx = 1;  logy = 1;  logz = 1; break;
    };
}

void CanvasConfiguration::setMargins(float left, float top, float right, float bottom)
{
  leftMargin   = left;
  topMargin    = top;
  rightMargin  = right;
  bottomMargin = bottom;
}

void CanvasConfiguration::setSize(int widthReq, int heightReq)
{
  width  = widthReq;
  height = heightReq;
}
