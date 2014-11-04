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
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.
You should have received a copy of the GNU General Public License
along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef __SORTEDSPIKEVIEWER_H_404FBFF2__
#define __SORTEDSPIKEVIEWER_H_404FBFF2__



#ifdef _WIN32
#include <Windows.h>
#endif
#include "../../JuceLibraryCode/JuceHeader.h"
#include "GenericProcessor.h"
#include "Editors/SortedSpikeViewerEditor.h"
#include "Editors/VisualizerEditor.h"
#include "Visualization/SpikeObject.h"
#include "Visualization/RippleObject.h"
#include "ParameterEstimation.h"
#include "RippleDetector.h"

class DataViewport;
class PCPlot;

struct RasterData
{
    int neuronID;
    int electrodeNum;
    unsigned long int timestamp;
};

struct RasterArray
{
    Array<RasterData> accruedRasterMarks;

    unsigned long int startRipple;
    unsigned long int stopRipple;
};

class SpikeSortingDisplay : public GenericProcessor
{
public:

    SpikeSortingDisplay();

    ~SpikeSortingDisplay();

    bool isSource()
    {
        return false;
    }

    bool isSink()
    {
        return true;
    }
    AudioProcessorEditor* createEditor();

    void updateSettings();
    bool enable();
    bool disable();
    int getNumberOfChannelsForElectrode(int i);
    String getNameForElectrode(int i);
    int getNumElectrodes();
    void handleEvent(int eventType, MidiMessage& event, int samplePosition);
    void getRippleDetectorPointer();

    struct Electrode
    {
        String name;

        int numChannels;

        PCPlot* pcPlot;

        int currentSpikeIndex;

        Array<SortedSpikeObject> mostRecentSortedSpikes;

        RasterArray forRasterPlot;
    };


    int svdCols;

    Array<Electrode> electrodes;

    bool ripplePresent;

    unsigned long int rippleTimestamp;
    unsigned long int ripplesDisplayed;
    bool pushRaster;
    bool canDrawOne;
    RippleDetector* node;
    RippleObject newRipple;

    void addPCPlotForElectrode(PCPlot* sp, int i);

    void process(AudioSampleBuffer& buffer, MidiBuffer& events, int& nSamples);

    void setParameter(int parameterIndex, float newValue);

    void removePCPlots();



    int displayBufferSize;
    bool redraw;

private:
    // private members and methods go here
    //
    // e.g.:
    //
    // float threshold;
    // bool state;
    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(SpikeSortingDisplay);
};



#endif  // __SORTEDSPIKEVIEWER_H_404FBFF2__
