/*
  ==============================================================================

    RippleDisplayCanvas.h
    Created: 20 Jun 2013 3:32:29pm
    Author:  sethi-san

  ==============================================================================
*/
/*
    ------------------------------------------------------------------

    This file is part of the Open Ephys GUI
    Copyright (C) 2013 Open Ephys

    ------------------------------------------------------------------

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

*/

#ifndef __RIPPLEDISPLAYCANVAS_H_D58AF97A__
#define __RIPPLEDISPLAYCANVAS_H_D58AF97A__

#include "../../../JuceLibraryCode/JuceHeader.h"

#include "../RippleDisplayNode.h"
#include "RippleObject.h"

#include "Visualizer.h"
#include <vector>

#define WAVE1 0
#define WAVE2 1
#define WAVE3 2
#define WAVE4 3
#define PROJ1x2 4
#define PROJ1x3 5
#define PROJ1x4 6
#define PROJ2x3 7
#define PROJ2x4 8
#define PROJ3x4 9

#define TETRODE_PLOT 1004
#define STEREO_PLOT  1002
#define SINGLE_PLOT  1001

#define MAX_NUMBER_OF_RIPPLE_SOURCES 128
#define MAX_N_CHAN 4

class RippleDisplayNode;
class RippleDisplay;
class RippleGenericAxes;
class RippleProjectionAxes;
class RippleWaveAxes;
class RipplePlot;

/**

  Displays Ripple waveforms and projections.

  @see RippleDisplayNode, RippleDisplayEditor, Visualizer

*/

class RippleDisplayCanvas : public Visualizer, public Button::Listener

{
public:
    RippleDisplayCanvas(RippleDisplayNode* n);
    ~RippleDisplayCanvas();

    void paint(Graphics& g);

    void refresh();

    void processRippleEvents();

    void beginAnimation();
    void endAnimation();

    void refreshState();

    void setParameter(int, float) {}
    void setParameter(int, int, int, float) {}

    void update();

    void resized();

    bool keyPressed(const KeyPress& key);

    void buttonClicked(Button* button);

private:

    RippleDisplayNode* processor;
    MidiBuffer* RippleBuffer;

    ScopedPointer<RippleDisplay> rippleDisplay;
    ScopedPointer<Viewport> viewport;

    ScopedPointer<UtilityButton> clearButton;

    bool newRipple;
    RippleObject Ripple;

    int scrollBarThickness;

    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(RippleDisplayCanvas);

};

class RippleDisplay : public Component
{
public:
    RippleDisplay(RippleDisplayCanvas*, Viewport*);
    ~RippleDisplay();

    void removePlots();
    void clear();
    void addRipplePlot(int numChannels, int electrodeNum);

    void paint(Graphics& g);

    void resized();

    void mouseDown(const MouseEvent& event);

    void plotRipple(const RippleObject& Ripple, int electrodeNum);

    int getTotalHeight()
    {
        return totalHeight;
    }

private:

    //void computeColumnLayout();
    //void initializeRipplePlots();
    //void repositionRipplePlots();

    int numColumns;

    int totalHeight;

    RippleDisplayCanvas* canvas;
    Viewport* viewport;

    OwnedArray<RipplePlot> ripplePlots;

    // float tetrodePlotMinWidth, stereotrodePlotMinWidth, singleElectrodePlotMinWidth;
    // float tetrodePlotRatio, stereotrodePlotRatio, singleElectrodePlotRatio;

};

/**

  Class for drawing the waveforms and projections of incoming Ripples.

*/

class RipplePlot : public Component, Button::Listener
{
public:
    RipplePlot(RippleDisplayCanvas*, int elecNum, int plotType);
    virtual ~RipplePlot();

    void paint(Graphics& g);
    void resized();

    void select();
    void deselect();

    void processRippleObject(const RippleObject& s);

    RippleDisplayCanvas* canvas;

    bool isSelected;

    int electrodeNumber;

    int nChannels;

    void initAxes();
    void getBestDimensions(int*, int*);

    void clear();

    float minWidth;
    float aspectRatio;

    void buttonClicked(Button* button);

private:


    int plotType;
    int nWaveAx;
    int nProjAx;

    bool limitsChanged;

    double limits[MAX_N_CHAN][2];

    OwnedArray<RippleProjectionAxes> pAxes;
    OwnedArray<RippleWaveAxes> wAxes;
    OwnedArray<UtilityButton> rangeButtons;
    Array<float> ranges;

    void initLimits();
    void setLimitsOnAxes();
    void updateAxesPositions();


    Font font;

};

/**

  Base class for drawing axes for Ripple visualization.

  @see RippleDisplayCanvas

*/


class RippleGenericAxes : public Component
{
public:

    RippleGenericAxes(int t);

    virtual ~RippleGenericAxes();

    virtual void updateRippleData(const RippleObject& s);

    void setXLims(double xmin, double xmax);
    void getXLims(double* xmin, double* xmax);
    void setYLims(double ymin, double ymax);
    void getYLims(double* ymin, double* ymax);

    void setType(int type);
    int getType();

    virtual void paint(Graphics& g) = 0;

    int roundUp(int, int);
    void makeLabel(int val, int gain, bool convert, char* s);

protected:
    double xlims[2];
    double ylims[2];

    RippleObject s;

    bool gotFirstRipple;

    int type;

    Font font;

    double ad16ToUv(int x, int gain);

};


/*

  Class for drawing Ripple waveforms.
*/

class RippleWaveAxes : public RippleGenericAxes
{
public:
    RippleWaveAxes(int channel);
    ~RippleWaveAxes() {}

    void updateRippleData(const RippleObject& s);

    void paint(Graphics& g);

    void plotRipple(const RippleObject& s, Graphics& g);

    void clear();

    void mouseMove(const MouseEvent& event);
    void mouseExit(const MouseEvent& event);
    void mouseDown(const MouseEvent& event);
    void mouseDrag(const MouseEvent& event);

    void setRange(float);
    float getRange() {return range;}

    //MouseCursor getMouseCursor();

private:

    Colour waveColour;
    Colour thresholdColour;
    Colour gridColour;

    bool drawGrid;

    float thresholdLevel;
    float x;

    void drawWaveformGrid(Graphics& g);

    void drawThresholdSlider(Graphics& g);

    Font font;

    Array<RippleObject> RippleBuffer;

    int RippleIndex;
    int bufferSize;

    float range;

    bool isOverThresholdSlider;
    bool isDraggingThresholdSlider;

    MouseCursor::StandardCursorType cursorType;

};





/*

  Class for drawing the peak projections of Ripple waveforms.
*/

class RippleProjectionAxes : public RippleGenericAxes
{
public:
    RippleProjectionAxes(int projectionNum);
    ~RippleProjectionAxes() {}

    void updateRippleData(const RippleObject& s);

    void paint(Graphics& g);

    void clear();

    void setRange(float, float);

    static void n2ProjIdx(int i, int* p1, int* p2);

private:

    void updateProjectionImage(uint16_t, uint16_t, uint16_t);

    void calcWaveformPeakIdx(const RippleObject&, int, int, int*, int*);

    int ampDim1, ampDim2;

    Image projectionImage;

    Colour pointColour;
    Colour gridColour;

    int imageDim;

    int rangeX;
    int rangeY;


};




#endif  // __RIPPLEDISPLAYCANVAS_H_D58AF97A__
