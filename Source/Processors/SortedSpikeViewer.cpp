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
#include <stdio.h>
#include "SortedSpikeViewer.h"
#include "Visualization/SortedSpikeCanvas.h"
#include "Channel.h"

SpikeSortingDisplay::SpikeSortingDisplay()
    : GenericProcessor("Sorted Spike Viewer"), displayBufferSize(5)
{
    redraw = false;
    ripplePresent = false;
    forRasterPlot.rippleID = 0;
    ripplesDisplayed = 0;
}

SpikeSortingDisplay::~SpikeSortingDisplay()
{
}

AudioProcessorEditor* SpikeSortingDisplay::createEditor()
{
    //std::cout<<"Creating SpikeDisplayCanvas."<<std::endl;

    editor = new SortedSpikeDisplayEditor(this);
    return editor;

}



void SpikeSortingDisplay::updateSettings()
{
    electrodes.clear();

    for (int i = 0; i < eventChannels.size(); i++)
    {
        if ((eventChannels[i]->eventType < 999) && (eventChannels[i]->eventType > SPIKE_BASE_CODE))
        {

            Electrode elec;
            elec.numChannels = eventChannels[i]->eventType - 100;
            elec.name = eventChannels[i]->name;
            elec.currentSpikeIndex = 0;
            elec.mostRecentSortedSpikes.ensureStorageAllocated(displayBufferSize);
            electrodes.add(elec);
        }
    }
    // Set K
    ProcessorGraph* gr = getProcessorGraph();
    Array<GenericProcessor*> p = gr->getListOfProcessors();
    ParameterEstimator* node;
    bool flag = false;
    for (int k=0;k<p.size();k++)
    {
        if (p[k]->getName() == "Parameter Estimator")
        {
            node = (ParameterEstimator*)p[k];
            flag = true;
        }

    }
    if (!flag)
    {
        std::cout << "Could not find a Parameter Estimator." << std::endl;
        svdCols = 3;
    }
    else
    {
        svdCols = node->SVDCols;
        std::cout << "SVDCols is now " << node->SVDCols <<  "  /  ";
    }
    // endof Set K
}

bool SpikeSortingDisplay::enable()
{
    std::cout << "SpikeDisplayNode::enable()" << std::endl;
    SortedSpikeDisplayEditor * editor = (SortedSpikeDisplayEditor*) getEditor();
    editor->enable();
    return true;
}


bool SpikeSortingDisplay::disable()
{
    std::cout << "SpikeDisplayNode disabled!" << std::endl;
    SortedSpikeDisplayEditor* editor = (SortedSpikeDisplayEditor*) getEditor();
    editor->disable();
    return true;
}

int SpikeSortingDisplay::getNumberOfChannelsForElectrode(int i)
{
    if (i > -1 && i < electrodes.size())
    {
        return electrodes[i].numChannels;
    } else {
        return 0;
    }
}

String SpikeSortingDisplay::getNameForElectrode(int i)
{

    if (i > -1 && i < electrodes.size())
    {
        return electrodes[i].name;
    } else {
        return " ";
    }
}

int SpikeSortingDisplay::getNumElectrodes()
{
    return electrodes.size();

}
void SpikeSortingDisplay::addPCPlotForElectrode(PCPlot* sp, int i)
{
    Electrode& e = electrodes.getReference(i);
    e.pcPlot = sp;

}

void SpikeSortingDisplay::removePCPlots()
{
    for (int i = 0; i < getNumElectrodes(); i++)
    {
        Electrode& e = electrodes.getReference(i);
        e.pcPlot = nullptr;
    }
}

void SpikeSortingDisplay::setParameter(int param, float val)
{
    redraw = true;
}

void SpikeSortingDisplay::process(AudioSampleBuffer& buffer,
                                  MidiBuffer& events,
                                  int& nSamples)
{
    checkForEvents(events);

    if (redraw)
    {
        for (int i = 0; i < getNumElectrodes(); i++)
        {

            Electrode& e = electrodes.getReference(i);

            // transfer buffered spikes to spike plot
            for (int j = 0; j < e.currentSpikeIndex; j++)
            {
                //std::cout << "Transferring spikes." << std::endl;
                e.pcPlot->processSortedSpikeObject(e.mostRecentSortedSpikes[j]);
                e.currentSpikeIndex = 0;
            }

        }
        redraw = false;
    }

    if(canDrawOne)
    {
        for (int i = 0; i < getNumElectrodes(); i++)
        {

            Electrode& e = electrodes.getReference(i);
            e.pcPlot->processRasterPlot(forRasterPlot, i);
        }
        canDrawOne = false;
        forRasterPlot.accruedRasterMarks.clearQuick();
    }

}

void SpikeSortingDisplay::handleEvent(int eventType, MidiMessage& event, int samplePosition)
{
    if (eventType == SORTEDSPIKE)
    {
        //std::cout << "reached handle event";
        const uint8_t* dataptr = event.getRawData();
        int bufferSize = event.getRawDataSize();

        if (bufferSize > 0)
        {
            //std::cout << "buffer is positive";
            SortedSpikeObject newSpike;

            bool isValid = unpackSortedSpike(&newSpike, dataptr, bufferSize);

            if (isValid)
            {
                //std::cout << "buffer is valid";
                int electrodeNum = newSpike.source;

                if (ripplePresent)
                {
                    if(newSpike.timestamp > forRasterPlot.startRipple)
                    {
                    RasterData newMark;
                    newMark.neuronID = newSpike.neuronID;
                    newMark.timestamp = newSpike.timestamp;
                    newMark.electrodeNum = electrodeNum;

                    forRasterPlot.accruedRasterMarks.add(newMark);
                    }
                }
                else
                {
                    if(newSpike.timestamp < forRasterPlot.stopRipple)
                    {
                    RasterData newMark;
                    newMark.neuronID = newSpike.neuronID;
                    newMark.timestamp = newSpike.timestamp;
                    newMark.electrodeNum = electrodeNum;

                    forRasterPlot.accruedRasterMarks.add(newMark);
                    }
                }


                Electrode& e = electrodes.getReference(electrodeNum);

                if (e.currentSpikeIndex < displayBufferSize)
                {
                    //  std::cout << "Adding spike " << e.currentSpikeIndex + 1 << std::endl;
                    e.mostRecentSortedSpikes.set(e.currentSpikeIndex, newSpike);
                    e.currentSpikeIndex++;
                }
            }
        }
    }
    if (eventType == RIPPLE)
    {
        const uint8_t* dataptr = event.getRawData();
        int bufferSize = event.getRawDataSize();

        if (bufferSize > 0)
        {


            bool isValid = unpackRipple(&newRipple, dataptr, bufferSize);

            if (isValid)
            {
                if (newRipple.start == 1)
                {
                    ripplePresent = true;
                    forRasterPlot.startRipple = newRipple.timestamp;
                    forRasterPlot.rippleID = newRipple.eventId;
                }
                if (newRipple.start == 0)
                {
                    ripplePresent = false;
                    forRasterPlot.stopRipple = newRipple.timestamp;
                    canDrawOne = true;
                }
            }
        }

        ripplePresent = true;
    }
}
