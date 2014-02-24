/*
  ==============================================================================

    RippleDetector.h
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

#ifndef __RIPPLEDETECTOR_H_8D9668F1__
#define __RIPPLEDETECTOR_H_8D9668F1__

#ifdef _WIN32
#include <Windows.h>
#endif
#include "Editors/RippleDetectorEditor.h"
#include "../../JuceLibraryCode/JuceHeader.h"
#include "GenericProcessor.h"
#include "Visualization/RippleObject.h"
#include <math.h>
#include <complex>
#include "DataThreads/RHD2000Thread.h"
#include "FPGAOutput.h"
#include "SourceNode.h"

class RippleDetector : public GenericProcessor

{
public:

    /** The class constructor, used to initialize any members. */
    RippleDetector();

    /** The class destructor, used to deallocate memory */
    ~RippleDetector();

    /** Determines whether the processor is treated as a source. */
    bool isSource()
    {
        return false;
    }

    /** Determines whether the processor is treated as a sink. */
    bool isSink()
    {
        return false;
    }

    /** Defines the functionality of the processor.

        The process method is called every time a new data buffer is available.

        Processors can either use this method to add new data, manipulate existing
        data, or send data to an external target (such as a display or other hardware).

        Continuous signals arrive in the "buffer" variable, event data (such as TTLs
        and spikes) is contained in the "events" variable, and "nSamples" holds the
        number of continous samples in the current buffer (which may differ from the
        size of the buffer).
         */
    void process(AudioSampleBuffer& buffer, MidiBuffer& events, int& nSamples);

    /** Any variables used by the "process" function _must_ be modified only through
        this method while data acquisition is active. If they are modified in any
        other way, the application will crash.  */
    void setParameter(int parameterIndex, float newValue);
	

    /** Called whenever the signal chain is altered. */
    void updateSettings();

    /** Called prior to start of acquisition. */
    bool enable();

    /** Called after acquisition is finished. */
    bool disable();

    /** Creates the RippleDetectorEditor. */
    AudioProcessorEditor* createEditor();

	// INTERNAL BUFFERS //

    /** Extra samples are placed in this buffer to allow seamless
        transitions between callbacks. */
    AudioSampleBuffer overflowBuffer;

	// CREATE AND DELETE ELECTRODES //

    /** Adds an electrode with n channels to be processed. */
    bool addElectrode(int nChans);

    /** Removes an electrode with a given index. */
    bool removeElectrode(int index);


    // EDIT AND QUERY ELECTRODE SETTINGS //

    /** Returns the number of channels for a given electrode. */
    int getNumChannels(int index);

    /** Edits the mapping between input channels and electrode channels. */
    void setChannel(int electrodeIndex, int channelNum, int newChannel);

    /** Returns the continuous channel that maps to a given
    	electrode channel. */
    int getChannel(int index, int chan);

    /** Sets the name of a given electrode. */
    void setElectrodeName(int index, String newName);

    /** */
    void setChannelActive(int electrodeIndex, int channelNum, bool active);

    /** */
    bool isChannelActive(int electrodeIndex, int channelNum);

    // RETURN STRING ARRAYS //

    /** Returns a StringArray containing the names of all electrodes */
    StringArray getElectrodeNames();

    /** Returns a list of possible electrode types (e.g., stereotrode, tetrode). */
    StringArray electrodeTypes;

    void setChannelThreshold(int electrodeNum, int channelNum, float threshold);

    double getChannelThreshold(int electrodeNum, int channelNum);

    void addFrequencyBins(int &fRes);

    //float normal_pdf(float x, float m, float s);

    //static const float inv_sqrt_2pi = 0.3989422804014327;




    void saveCustomParametersToXml(XmlElement* parentElement);
    void loadCustomParametersFromXml();

private:

    // private members and methods go here
    //
    // e.g.:
    //
    // float threshold;
    // bool state;

	/** Reference to a continuous buffer. */
    AudioSampleBuffer dataBuffer;

    float getDefaultThreshold();

    double logLikelihood;

    float delta;
    float term;
    float np;
    float sinp;
    int numOfEvents;


    int overflowBufferSize;

    int sampleIndex;

    Array<int> electrodeCounter;

    float getNextSample(int& chan);
    float getCurrentSample(int& chan);
    bool samplesAvailable(int& nSamples);


    struct FrequencyBin
    {

       std::complex<double> binValue;

    };

    Array<FrequencyBin*> fftBins;


    bool useOverflowBuffer;
    bool binsActive;
    int binCount;
    int transformLength;
    Array<int> binNumbers;

    int currentElectrode;
    int currentChannelIndex;
    int currentIndex;

    int startIndex;
    int stopIndex;
    bool ripplePresent;
    //int rippleContinues;



    struct Electrode
    {

        String name;

        int numChannels;
        int prePeakSamples, postPeakSamples;
        int lastBufferIndex;

        int* channels;
        double* thresholds;
        bool* isActive;

        int fftValidCount;
        double paramAveragingCount;
        float mu, sigma;
        double rippleAmplitude;
        float amplitudeCount;
        float sumOfSquaresOfDifferences;
        float partialSum;
        double maxAmplitudeFound;

    };

    uint8_t* rippleBuffer;///[256];

    Array<Electrode*> electrodes;

    // void createSpikeEvent(int& peakIndex,
    // 					  int& electrodeNumber,
    // 					  int& currentChannel,
    // 					  MidiBuffer& eventBuffer);

    void addRippleEvent(RippleObject* s, MidiBuffer& eventBuffer, int startIndex);
    void addWaveformToRippleObject(RippleObject* s,
                                  int& startIndex, int& nSamples,
                                  int& electrodeNumber,
                                  int& currentChannel);

    void resetElectrode(Electrode*);

    RHD2000Thread* dataThread;
    FPGAThread* testThread;

    void timerCallback();

    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(RippleDetector);

};


#endif  // __RIPPLEDETECTOR_H_8D9668F1__
