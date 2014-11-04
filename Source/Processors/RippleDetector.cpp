/*
  ==============================================================================

    RippleDetector.cpp
    Created: 24 May 2013 12:09:14pm
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

#include "RippleDetector.h"
#include <stdio.h>
#include <math.h>
#include <complex>
//#include "../../../../../../../usr/include/fftw3.h"


RippleDetector::RippleDetector()
    : GenericProcessor("Ripple Detector"),
    overflowBuffer(2,1), dataBuffer(overflowBuffer),
      overflowBufferSize(1), currentElectrode(-1), binsActive(false), ripplePresent(false), delta(0), term(0) , sinp(0), np(0), numOfEvents(1) //rippleContinues(0)
{

    //parameters.add(Parameter("thresh", 0.0, 500.0, 200.0, 0));
	electrodeTypes.add("single electrode");
    	electrodeTypes.add("stereotrode");
    	electrodeTypes.add("tetrode");
	
	for (int i = 0; i < electrodeTypes.size()+1; i++)
	    {
		electrodeCounter.add(0);
	    }
    rippleBuffer = new uint8_t[MAX_RIPPLE_BUFFER_LEN];
}

RippleDetector::~RippleDetector()
{

}

AudioProcessorEditor* RippleDetector::createEditor()
{
    editor = new RippleDetectorEditor(this, true);
    return editor;
}

void RippleDetector::updateSettings()
{

    if (getNumInputs() > 0)
        overflowBuffer.setSize(getNumInputs(), overflowBufferSize);

    for (int i = 0; i < electrodes.size(); i++)
    {

        Channel* ch = new Channel(this, i);
        ch->isEventChannel = true;
        ch->eventType = electrodes[i]->numChannels;
        ch->name = electrodes[i]->name;

        eventChannels.add(ch);
    }



}

bool RippleDetector::addElectrode(int nChans)
{

    std::cout << "Adding electrode with " << nChans << " channels." << std::endl;

    int firstChan;

    if (electrodes.size() == 0)
    {
        firstChan = 0;
    }
    else
    {
        Electrode* e = electrodes.getLast();
        firstChan = *(e->channels+(e->numChannels-1))+1;
    }

    if (firstChan + nChans > getNumInputs())
    {
        return false;
    }

    int currentVal = electrodeCounter[nChans];
    electrodeCounter.set(nChans,++currentVal);

    String electrodeName;

    // hard-coded for tetrode configuration
    if (nChans < 3)
        electrodeName = electrodeTypes[nChans-1];
    else
        electrodeName = electrodeTypes[nChans-2];

    String newName = electrodeName.substring(0,1);
    newName = newName.toUpperCase();
    electrodeName = electrodeName.substring(1,electrodeName.length());
    newName += electrodeName;
    newName += " ";
    newName += electrodeCounter[nChans];

    Electrode* newElectrode = new Electrode();

    newElectrode->name = newName;
    newElectrode->numChannels = nChans;
    newElectrode->prePeakSamples = 0;
    newElectrode->postPeakSamples = 256;
    newElectrode->thresholds = new double[nChans];
    newElectrode->isActive = new bool[nChans];
    newElectrode->channels = new int[nChans];
    newElectrode->fftValidCount = 0;
    newElectrode->paramAveragingCount = 0;
    newElectrode->mu =0;
    newElectrode->sigma = 0;
    newElectrode->amplitudeCount = 0;
    newElectrode->sumOfSquaresOfDifferences = 0;
    newElectrode->partialSum = 0;
    newElectrode->maxAmplitudeFound = -10;
    newElectrode->rippleAmplitude = 1400;

    for (int i = 0; i < nChans; i++)
    {
        *(newElectrode->channels+i) = firstChan+i;
        *(newElectrode->thresholds+i) = getDefaultThreshold();
        *(newElectrode->isActive+i) = true;
    }

    resetElectrode(newElectrode);

    electrodes.add(newElectrode);

    return true;

}
float RippleDetector::getDefaultThreshold()
{
    return 200.0f;
}
void RippleDetector::resetElectrode(Electrode* e)
{
    e->lastBufferIndex = 0;
}

bool RippleDetector::isChannelActive(int electrodeIndex, int i)
{
    return *(electrodes[electrodeIndex]->isActive+i);
}

bool RippleDetector::enable()
{

    useOverflowBuffer = false;
    return true;
}

bool RippleDetector::disable()
{

    for (int n = 0; n < electrodes.size(); n++)
    {
        resetElectrode(electrodes[n]);
    }

    return true;
}
void RippleDetector::setParameter(int parameterIndex, float newValue)
{
    electrodes[newValue]->rippleStatus.removeRange(0,2);
}

bool RippleDetector::samplesAvailable(int& nSamples)
{

    if (sampleIndex >= nSamples)
    {
        return false;
    }
    else
    {
        return true;
    }

}
void RippleDetector::setChannel(int electrodeIndex, int channelNum, int newChannel)
{

    std::cout << "Setting electrode " << electrodeIndex << " channel " << channelNum <<
               " to " << newChannel << std::endl;

    *(electrodes[electrodeIndex]->channels+channelNum) = newChannel;
}

int RippleDetector::getNumChannels(int index)
{
    return electrodes[index]->numChannels;
}

int RippleDetector::getChannel(int index, int i)
{
    return *(electrodes[index]->channels+i);
}
StringArray RippleDetector::getElectrodeNames()
{
    StringArray names;

    for (int i = 0; i < electrodes.size(); i++)
    {
        names.add(electrodes[i]->name);
    }

    return names;
}
void RippleDetector::setElectrodeName(int index, String newName)
{
    electrodes[index-1]->name = newName;
}
bool RippleDetector::removeElectrode(int index)
{

    // std::cout << "Spike detector removing electrode" << std::endl;

    if (index > electrodes.size() || index < 0)
        return false;

    electrodes.remove(index);
    return true;
}
double RippleDetector::getChannelThreshold(int electrodeNum, int channelNum)
{
    return *(electrodes[electrodeNum]->thresholds+channelNum);
}
void RippleDetector::setChannelThreshold(int electrodeNum, int channelNum, float thresh)
{
    currentElectrode = electrodeNum;
    currentChannelIndex = channelNum;
    setParameter(99, thresh);
}
void RippleDetector::setChannelActive(int electrodeIndex, int i, bool active)
{
    *(electrodes[electrodeIndex]->isActive+i) = active;
}

void RippleDetector::addWaveformToRippleObject(RippleObject* s,
                                             int& startIndex, int& nSamples,
                                             int& electrodeNumber,
                                             int& currentChannel)
{
  /*  int spikeLength = 150;
   // std::cout<<"The length of the spike stored in data is "<<nSamples<<"   "<<startIndex<<"     "<<transformLength/20<<"    "<<spikeLength;
    s->nSamples = spikeLength;

    int chan = *(electrodes[electrodeNumber]->channels+currentChannel);

    s->gain[currentChannel] = (int)(1.0f / channels[chan]->bitVolts)*1000;
    s->threshold[currentChannel] = (int) *(electrodes[electrodeNumber]->thresholds+currentChannel) / channels[chan]->bitVolts * 1000;

    // cycle through buffer
    //if (isChannelActive(electrodeNumber, currentChannel))
  /*  if (currentChannel < 1)
    {
        for (int sample = 0; sample < spikeLength; sample++)
        {

            // warning -- be careful of bitvolts conversion
            s->data[currentIndex] = uint16(getNextSample(*(electrodes[electrodeNumber]->channels+currentChannel)) / channels[chan]->bitVolts + 32768);

            currentIndex++;
            sampleIndex++;

            //std::cout << currentIndex << std::endl;

        }
    }

    sampleIndex -= spikeLength; // reset sample index
*/

}
void RippleDetector::addRippleEvent(RippleObject* s, MidiBuffer& eventBuffer, int startIndex)
{


    int numBytes = packRipple(s, rippleBuffer, 256);
    eventBuffer.addEvent(rippleBuffer, numBytes, startIndex);

}
void RippleDetector::addFrequencyBins(int& fRes)
{
    float fs = RippleDetector::getSampleRate();

    int N = ceil(fs/float(fRes));
    if ( (N%2) == 1)
    {
        N=N+1;
    }
    std::cout<<N<<" is the number of sample points in the FFT"<<std::endl;

    binCount=0;
    for (int i = 0; i <= (N/2)-1 ; i++)
    {
        if (i >= 150*N/fs && i<=250*N/fs)
        {
            binCount += 1;
            binNumbers.add(i);
        }

    }
        std::cout<<binCount<<" ris the number of bins between 150 & 250 Hz"<<std::endl;

    if (binCount == 0)
    {
        return;
    }

    for (int i = 0; i < electrodes.size(); i++)
    {
        Electrode* electrode= electrodes[i];

        for (int j = 0; j < electrode->numChannels; j++)
        {
            if (*(electrode->isActive+j))
            {
                for (int k = 0; k < binCount; k++)
                {
                    FrequencyBin* newBin = new FrequencyBin();
                    newBin->binValue = std::polar(0,0);
                    fftBins.add(newBin);

                }
            }

        }
    }
    binsActive = true;
    transformLength = N;
}

/*float RippleDetector::normal_pdf(float x, float m, float s)
{
    float a = (x - m) / s;

    return (inv_sqrt_2pi / s) * std::exp(-0.5f * a * a);
}*/

/*void RippleDetector::timerCallback()
{
    //testThread->setOutputLow();

    sendActionMessage("LO");
    int ttlOutArray[16] = {0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1};
    dataThread->setTTLPins(ttlOutArray);
    //std::cout << "reached callback" << std::endl;

    stopTimer();


}*/

void RippleDetector::handleEvent(int eventType, MidiMessage& event, int sampleNum)
{

    if (eventType == TIMESTAMP)
    {
        const uint8* dataptr = event.getRawData();

        memcpy(&timestamp, dataptr + 4, 8); // remember to skip first four bytes
    }


}
void RippleDetector::process(AudioSampleBuffer& buffer,
                               MidiBuffer& events,
                               int& nSamples)
{

    Electrode* electrode;
    dataBuffer = buffer;
    checkForEvents(events);


    for(int i = 0; i < electrodes.size(); i++)
    {
        electrode = electrodes[i];
        sampleIndex = electrode->lastBufferIndex -1;// subtract 1 to account for
        // increment at start of getNextSample()

       while (samplesAvailable(nSamples))
        {
            sampleIndex++;

            for(int chan=0; chan < 1 ; chan++) //detecting on 1 channel only currently. to be expanded. CUSUM -> multichannel/distributed CUSUM
            {
                if (*(electrode->isActive+chan))
                {
                    int currentChannel = *(electrode->channels+chan);

                    double var = getNextSample(currentChannel);

                        if (electrode->paramAveragingCount < 30000) //averaging for 1 s to find noise mu and sigma (needs to be changed to that each electrode gets its own noise measurement)
                        {
                            electrode->paramAveragingCount++;
                            delta = (var - electrode->mu);
                            electrode->mu = electrode->mu + (delta/electrode->paramAveragingCount);
                            electrode->sumOfSquaresOfDifferences = electrode->sumOfSquaresOfDifferences + delta*(var - electrode->mu);
                        }
                        else  //Begin the hunt for ripples!
                        {
                            if (electrode->paramAveragingCount == 30000)
                            {
                                electrode->sigma = sqrt(electrode->sumOfSquaresOfDifferences/electrode->paramAveragingCount);
                                electrode->paramAveragingCount++;
                                electrode->rippleAmplitude = electrode->mu + 5*electrode->sigma;
                                std::cout << "Mean = " << electrode->mu << " Sigma is = " << electrode->sigma << std::endl;

                            }

                            np = ((var-electrode->mu)/electrode->sigma);
                            np = 0.5*np*np - 8.0; //norm pdf (var) > 8 counts towards ripple probability

                            logLikelihood = np;

                            if (ripplePresent == false)
                            {

                                electrode->partialSum += logLikelihood;
                                electrode->partialSum = std::max(electrode->partialSum,float(0.0)); //CUSUM

                                if (electrode->partialSum > 5)//*(electrode->thresholds+chan)) //ripple detected
                                {

                                    //std::cout << "PS is " << electrode->partialSum << "LL is " << logLikelihood << " sinp " << sinp << " np " << np << std::endl;
                                    electrode->partialSum = 0.0; //reset algorithm
                                    std::cout << "Ripple Detected!" << std::endl;
                                    ripplePresent = true;

                                        electrode->rippleStatus.add(timestamp + sampleIndex);
                                        RippleObject newRipple;
                                        newRipple.start = 1;
                                        newRipple.eventType = RIPPLE;
                                        newRipple.eventId = 1;
                                        newRipple.eventChannel = 0;


                                        addRippleEvent(&newRipple, events, sampleIndex);
                                }
                                break;
                            }


                            else // if ripple is still present
                            {
                                electrode->partialSum = electrode->partialSum + logLikelihood;
                                electrode->partialSum = std::min(electrode->partialSum, float(0.0));



                                if (electrode->partialSum < -150) // need to play around with this a bit. the PS needs to hit this thresh
                                    //after or close to when the ripple ends to prevent re-detection of the same ripple. Not an issue when disrupting.
                                {
                                    //std::cout << electrode->partialSum << std::endl;
                                    ripplePresent = false;
                                    electrode->rippleStatus.add(timestamp + sampleIndex);
                                    electrode->partialSum = 0.0; //reset again
                                    RippleObject newRipple;
                                    newRipple.start = 0;
                                    newRipple.eventType = RIPPLE;
                                    newRipple.eventId = 1;
                                    newRipple.eventChannel = 0;


                                    addRippleEvent(&newRipple, events, sampleIndex);

                                }


                            }



                            }


                     }
                }
               }

        electrode->lastBufferIndex = 0;
    }




}

void RippleDetector::saveCustomParametersToXml(XmlElement* parentElement)
{

    for (int i = 0; i < electrodes.size(); i++)
    {
         XmlElement* electrodeNode = parentElement->createNewChildElement("ELECTRODE");
         electrodeNode->setAttribute("name", electrodes[i]->name);
         electrodeNode->setAttribute("numChannels", electrodes[i]->numChannels);
         electrodeNode->setAttribute("prePeakSamples", electrodes[i]->prePeakSamples);
         electrodeNode->setAttribute("postPeakSamples", electrodes[i]->postPeakSamples);

         for (int j = 0; j < electrodes[i]->numChannels; j++)
        {
            XmlElement* channelNode = electrodeNode->createNewChildElement("SUBCHANNEL");
            channelNode->setAttribute("ch",*(electrodes[i]->channels+j));
            channelNode->setAttribute("thresh",*(electrodes[i]->thresholds+j));
            channelNode->setAttribute("isActive",*(electrodes[i]->isActive+j));

        }
    }


}

void RippleDetector::loadCustomParametersFromXml()
{

    if (parametersAsXml != nullptr)
    {
        // use parametersAsXml to restore state

        int electrodeIndex = -1;

        forEachXmlChildElement(*parametersAsXml, xmlNode)
        {
           if (xmlNode->hasTagName("ELECTRODE"))
            {

                electrodeIndex++;

                int channelsPerElectrode = xmlNode->getIntAttribute("numChannels");

                RippleDetectorEditor* sde = (RippleDetectorEditor*) getEditor();
                sde->addElectrode(channelsPerElectrode);

                setElectrodeName(electrodeIndex+1, xmlNode->getStringAttribute("name"));

                int channelIndex = -1;

                forEachXmlChildElement(*parametersAsXml, channelNode)
                {
                    if (channelNode->hasTagName("SUBCHANNEL"))
                    {
                        channelIndex++;

                        setChannel(electrodeIndex, channelIndex, channelNode->getIntAttribute("ch"));
                        setChannelThreshold(electrodeIndex, channelIndex, channelNode->getDoubleAttribute("thresh"));
                        setChannelActive(electrodeIndex, channelIndex, channelNode->getBoolAttribute("isActive"));
                    }
                }

            }
        }
    }
}

float RippleDetector::getNextSample(int& chan)
{


/*
    //if (useOverflowBuffer)
    //{
    if (sampleIndex < 0)
    {
        // std::cout << "  sample index " << sampleIndex << "from overflowBuffer" << std::endl;
        int ind = overflowBufferSize + sampleIndex;

        if (ind < overflowBuffer.getNumSamples())
        {
            //std::cout<<ind<<std::endl;
            return *overflowBuffer.getSampleData(chan, ind);
        }
        else
            return 0;

    }
    else*/
    {
        //  useOverflowBuffer = false;
        // std::cout << "  sample index " << sampleIndex << "from regular buffer" << std::endl;

        if (sampleIndex < dataBuffer.getNumSamples())
        {
            //std::cout<<sampleIndex<<std::endl;
            return *dataBuffer.getSampleData(chan, sampleIndex);
        }
        else
            return 0;
    }
    //} else {
    //    std::cout << "  sample index " << sampleIndex << "from regular buffer" << std::endl;
    //     return *dataBuffer.getSampleData(chan, sampleIndex);
    //}

}

float RippleDetector::getCurrentSample(int& chan)
{

    // if (useOverflowBuffer)
    // {
    //     return *overflowBuffer.getSampleData(chan, overflowBufferSize + sampleIndex - 1);
    // } else {
    //     return *dataBuffer.getSampleData(chan, sampleIndex - 1);
    // }

    /*if (sampleIndex < 1)
    {
        //std::cout << "  sample index " << sampleIndex << "from overflowBuffer" << std::endl;
        //std::cout<<overflowBufferSize + sampleIndex - 1<<std::endl;
        return *overflowBuffer.getSampleData(chan, overflowBufferSize + sampleIndex - 1);
    }
    else*/
    if (sampleIndex > 0)
    {
        //  useOverflowBuffer = false;
        // std::cout << "  sample index " << sampleIndex << "from regular buffer" << std::endl;
        //std::cout<<overflowBufferSize + sampleIndex - 1<<std::endl;
        return *dataBuffer.getSampleData(chan, sampleIndex - 1);
    }
    else
        return 0;
    //} else {

}
