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

#include "TTLTrigger.h"
#include "FPGAOutput.h"

TTLTrigger::TTLTrigger()
    : GenericProcessor("TTLTrigger") , firstTime(true) , enabled(true)

{

    //parameters.add(Parameter("thresh", 0.0, 500.0, 200.0, 0));

}

TTLTrigger::~TTLTrigger()
{

}



void TTLTrigger::setParameter(int parameterIndex, float newValue)
{
    editor->updateParameterButtons(parameterIndex);
}

void TTLTrigger::process(AudioSampleBuffer& buffer,
                               MidiBuffer& events,
                               int& nSamples)
{



    if (dataThread != 0)
    {


        {
            checkForEvents(events);
        }

    }



}


void TTLTrigger::handleEvent(int eventType, MidiMessage& event, int sampleNum)
{

    //

        if (eventType == RIPPLE )
        {
            //std::cout << "reached here" << std::endl;
            //std::cout << "Get Bitvolts" << std::endl;
            //dataThread = (RHD2000Thread*) s->getThread();
            //std::cout<<"Testing"<<"End of Testing"<<std::endl;
            if (enabled == true)
            {
            int ttlOutArray[16] = {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};
            dataThread->setTTLPins(ttlOutArray);
            //std::cout << dataThread->isAcquisitionActive() << std::endl;
            enabled = false;
            }
            else
            {
                int ttlOutArray[16] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
            dataThread->setTTLPins(ttlOutArray);
            //std::cout << dataThread->isAcquisitionActive() << std::endl;
            enabled = true;
            }



        }

}

void TTLTrigger::updateSettings()
{
    removeAllActionListeners();

    GenericProcessor* src;
    GenericProcessor* lastSrc;

    lastSrc = getSourceNode();
    src = getSourceNode();

    while (src != 0)
    {
        lastSrc = src;
        src = lastSrc->getSourceNode();
    }

    if (lastSrc != 0)
    {
        SourceNode* s = (SourceNode*) lastSrc;
        addActionListener(s);
        std::cout << "FPGA Output node communicating with " << lastSrc->getName() << std::endl;
        dataThread = (RHD2000Thread*) s->getThread();
    }
    else
    {
        std::cout << "FPGA Output couldn't find a source" << std::endl;
    }


    //dataThread = (FPGAThread*) s->getThread();
}
