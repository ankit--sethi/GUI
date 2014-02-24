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
/*
==============================================================================

TTLTrigger.h
Created: 20 Jan 2014 1:26:00pm
Author: sethisan

==============================================================================
*/

#ifndef __TTLTRIGGER_H_581AD2__
#define __TTLTRIGGER_H_581AD2__


#ifdef _WIN32
#include <Windows.h>
#endif

#include "../../JuceLibraryCode/JuceHeader.h"
#include "GenericProcessor.h"
#include "SourceNode.h"
#include "DataThreads/RHD2000Thread.h"
#include "FPGAOutput.h"


class TTLTrigger : public GenericProcessor
{
public:

    TTLTrigger();
    ~TTLTrigger();

    bool isSource()
    {
        return false;
    }


    bool isSink()
    {
        return true;
    }

    void updateSettings();
    void process(AudioSampleBuffer& buffer, MidiBuffer& events, int& nSamples);
    void setParameter(int parameterIndex, float newValue);
    void handleEvent(int eventType, MidiMessage& event, int sampleNum);

    bool firstTime;
    bool enabled;




private:

    RHD2000Thread* dataThread;

    void timerCallback();

    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(TTLTrigger);

};

#endif // __TTLTRIGGER_H_581AD2__
