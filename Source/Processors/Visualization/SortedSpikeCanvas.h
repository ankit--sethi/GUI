/*
  ==============================================================================

    SortedSpikeCanvas.h
    Created: 18 Aug 2014 1:35:47pm
    Author:  sethisan

  ==============================================================================
*/

#ifndef __SORTEDSPIKECANVAS_H_4FA3496D__
#define __SORTEDSPIKECANVAS_H_4FA3496D__

#include "../../../JuceLibraryCode/JuceHeader.h"

#include "../SortedSpikeViewer.h"
#include "SpikeObject.h"
#include "SpikeDisplayCanvas.h"
#include "../../AccessClass.h"
#include "Visualizer.h"
#include <vector>

#define K_IS_TWO 7002
#define K_IS_THREE 7003
#define K_IS_FOUR 7004
#define K_IS_FIVE 7005

#define PC1x2 1
#define PC1x3 2
#define PC1x4 3
#define PC1x5 4
#define PC2x3 5
#define PC2x4 6
#define PC2x5 7
#define PC3x4 8
#define PC3x5 9
#define PC4x5 10

#define MAX_N_PROJ_COMBS 10

class SpikeSortingDisplay;
class GenericAxes;
class PCAxes;
class PCPlot;
class PCDisplay;
struct RasterArray;
struct RasterData;

class SortedSpikeCanvas : public Visualizer, public Button::Listener, public AccessClass
{
public:
    SortedSpikeCanvas(SpikeSortingDisplay* n);
    ~SortedSpikeCanvas();

    void paint(Graphics& g);

    void refresh();

    void beginAnimation();
    void endAnimation();

    void refreshState();

    void processSpikeEvents();

    void setParameter(int, float) {}
    void setParameter(int, int, int, float) {}

    void update();

    void resized();

    bool keyPressed(const KeyPress& key);

    void buttonClicked(Button* button);

    SpikeSortingDisplay* processor;

private:

    ScopedPointer<PCDisplay> pcDisplay;
    ScopedPointer<Viewport> viewport;
    ScopedPointer<UtilityButton> clearButton;

    bool newSpike;
    SpikeObject spike;

    int scrollBarThickness;

    ScopedPointer<UtilityButton> lockThresholdsButton;

    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(SortedSpikeCanvas);

};

class PCDisplay : public Component
{

public:
    PCDisplay(SortedSpikeCanvas*, Viewport*);
    ~PCDisplay();

    void removePlots();
    void clear();
    PCPlot* addPCPlot(int numChannels, int electrodeNum, String name);

    void paint(Graphics& g);

    void resized();

    void mouseDown(const MouseEvent& event);

    void plotSpike(const SortedSpikeObject& spike, int electrodeNum);

    int getTotalHeight()
    {
        return totalHeight;
    }

    int svdCols;

private:

    int numColumns;

    int totalHeight;

    SortedSpikeCanvas* canvas;
    Viewport* viewport;

    OwnedArray<PCPlot> pcPlots;

};
class RasterPlot;
class PCPlot : public Component, Button::Listener
{
    public:
    PCPlot(SortedSpikeCanvas*, int elecNum, int dim, String name_);
    virtual ~PCPlot();

    void paint(Graphics& g);
    void resized();

    void select();
    void deselect();

    void processSortedSpikeObject(const SortedSpikeObject& s);
    void processRasterPlot(const RasterArray& r, int electrodeNum);

    SortedSpikeCanvas* canvas;

    bool isSelected;

    int electrodeNumber;

    int nChannels;

    void initAxes();

    void clear();

    float minWidth;
    float aspectRatio;

    void buttonClicked(Button* button);

private:

    bool limitsChanged;

    int svdCol;

    double limits[MAX_N_PROJ_COMBS][2];

    OwnedArray<PCAxes> pAxes;
    OwnedArray<UtilityButton> rangeButtons;
    Array<float> ranges;

    int nPCAx;
    void initLimits();
    void setLimitsOnAxes();
    void updateAxesPositions();

    ScopedPointer <RasterPlot> raster;

    String name;

    Font font;
};

class RasterPlot : public Component
{
public:
    RasterPlot();
    ~RasterPlot() {}

    bool updateRasterData(const RasterArray& s, int electrodeNum);
    void situateRasterData(const RasterData& s, unsigned long int startTimestamp, unsigned long int endTimeStamp);
    Image ProjectionImage;
    //float startTimeStampPixel;
    //float endTimeStampPixel;
    void clear();
    void paint(Graphics& g);
};

class PCAxes : public GenericAxes
{
public:
    PCAxes(int projectionNum, int dictlength);
    ~PCAxes() {}

    bool updateSpikeData(const SortedSpikeObject& s);
    bool updateSpikeData(const SortedSpikeObject& s, int index);

    void paint(Graphics& g);

    void clear();

    void setRange(float, float);

    void drawGrid(Graphics& g);
    void updatePrincipalComponents(const SortedSpikeObject& s, int index);

    Image projectionImage;
private:

    void updateProjectionImage(uint16_t, uint16_t, uint16_t);

    int ampDim1, ampDim2;



    Colour pointColour;
    Colour gridColour;

    int svdCols;

    float rangeX;
    float rangeY;

    int spikesReceivedSinceLastRedraw;

};

#endif  // __SORTEDSPIKECANVAS_H_4FA3496D__
