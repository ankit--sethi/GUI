/*
  ==============================================================================

    RippleDisplayNode.cpp
    Created: 20 Jun 2013 3:20:26pm
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
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program. If not, see <http://www.gnu.org/licenses/>.

*/

#include "RippleDisplayNode.h"
#include <stdio.h>

#include "Channel.h"

RippleDisplayNode::RippleDisplayNode()
    : GenericProcessor("Ripple Viewer"),
      bufferSize(0)

{
    // displayBuffer = new AudioSampleBuffer(8, 100);
    eventBuffer = new MidiBuffer();
}

RippleDisplayNode::~RippleDisplayNode()
{
    //deleteAndZero(displayBuffer);
    //deleteAndZero(eventBuffer);
}

AudioProcessorEditor* RippleDisplayNode::createEditor()
{
    std::cout<<"Creating RippleDisplayCanvas."<<std::endl;

    editor = new RippleDisplayEditor(this);
    return editor;

}

// void RippleDisplayNode::updateSettings()
// {
// //std::cout << "Setting num inputs on RippleDisplayNode to " << getNumInputs() << std::endl;

// }

// void RippleDisplayNode::updateVisualizer()
// {

// }

bool RippleDisplayNode::enable()
{
    std::cout<<"RippleDisplayNode::enable()"<<std::endl;
    RippleDisplayEditor* editor = (RippleDisplayEditor*) getEditor();
    editor->enable();
    return true;

}

bool RippleDisplayNode::disable()
{
    std::cout<<"RippleDisplayNode disabled!"<<std::endl;
    RippleDisplayEditor* editor = (RippleDisplayEditor*) getEditor();
    editor->disable();
    return true;
}

int RippleDisplayNode::getNumberOfChannelsForElectrode(int elec)
{
    std::cout<<"RippleDisplayNode::getNumberOfChannelsForInput(" << elec << ")"<<std::endl;

    int electrodeIndex = -1;

    for (int i = 0; i < eventChannels.size(); i++)
    {
        if (eventChannels[i]->eventType < 999)
        {
            electrodeIndex++;

            if (electrodeIndex == elec)
            {
                std::cout << "Electrode " << elec << " has " << eventChannels[i]->eventType << " channels" << std::endl;
                return eventChannels[i]->eventType;
            }
        }
    }

    return 0;
}

int RippleDisplayNode::getNumElectrodes()
{
    int nElectrodes = 0;

    for (int i = 0; i < eventChannels.size(); i++)
    {
        if (eventChannels[i]->eventType < 999)
        {
            nElectrodes++;
        }
    }

    return nElectrodes;

}


void RippleDisplayNode::setParameter(int param, float val)
{
    std::cout<<"Got Param:"<< param<< " with value:"<<val<<std::endl;
}



void RippleDisplayNode::process(AudioSampleBuffer& buffer, MidiBuffer& events, int& nSamples)
{

    checkForEvents(events); // automatically calls 'handleEvent

}

void RippleDisplayNode::handleEvent(int eventType, MidiMessage& event, int samplePosition)
{

    //std::cout << "Received event of type " << eventType << std::endl;

    if (eventType == RIPPLE)
    {

        eventBuffer->addEvent(event, 0);

    }

}
