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

    ParameterEstimation.h
    Created: 15 May 2014 9:29:25am
    Author:  sethisan

  ==============================================================================
*/


#ifndef __PARAMETERESTIMATION_H_1B5A8C4B__
#define __PARAMETERESTIMATION_H_1B5A8C4B__

#ifdef _WIN32
#include <Windows.h>
#endif

#include "../../JuceLibraryCode/JuceHeader.h"
#include "GenericProcessor.h"
#include "Editors/ParameterEstimatorEditor.h"
#include "Visualization/SpikeObject.h"
#include <algorithm>    // std::sort
#include <list>
#include <queue>

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
    float mu, sigma, musqrd;
    float sumOfSquaresOfDifferences;


};


class TwoDimMatrix
{
    public:
    TwoDimMatrix()
        {
            xSize = 0; ySize = 0;
        }
       ~TwoDimMatrix()
        {

        }

        float get( unsigned int x, unsigned int y ) const
        {
            return data[ x + y * ySize ];
        }
        void setdim( unsigned int x, unsigned int y )
        {
            xSize = x;
            ySize = y;
        }
        int getrowdim()
        {
            return xSize;
        }
        int getcoldim()
        {
            return ySize;
        }


        void set( unsigned int x, unsigned int y, float element )
        {
            data[x + y * ySize] = element;

        }
        void addnew( float element )
        {
            data.push_back(element);
        }


    private:
        unsigned int xSize;
        unsigned int ySize;
        std::vector<float> data;


};





class SVDjob
{
public:
    //SVDjob(Array<SpikeObject> _spikes, bool _reportDone);
    SVDjob() {}
    ~SVDjob();
    void SVDsetdim(Array<SpikeObject> _spikes, bool _reportDone);
    float **U;
    Array<SpikeObject> spikes;
    int dim;
    TwoDimMatrix dict;
    bool reportDone;
    int svdcmp(float **a, int nRows, int nCols, float *w, float **v);
private:

    float pythag(float a, float b);

};

class SVDcomputingThread : juce::Thread
{
public:
    SVDcomputingThread();
    void run(); // computes SVD on spike waveforms
    void addSVDjob(SVDjob job);

    std::queue<SVDjob> jobs;

    SVDjob J;

};



class ParameterEstimator : public GenericProcessor

{
public:

    /** The class constructor, used to initialize any members. */
    ParameterEstimator();

    /** The class destructor, used to deallocate memory */
    ~ParameterEstimator();

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

    bool samplesAvailable(int& nSamples);
    float getNextSample(int& chan);
    float getCurrentSample(int& chan);

    void updateSettings();

    bool isChannelActive(int electrodeIndex, int channelNum);

    AudioProcessorEditor* createEditor();

    bool addElectrode(int nChans);

    /** Used to alter parameters of data acquisition. */

    void setChannelThreshold(int electrodeNum, int channelNum, float threshold);

    double getChannelThreshold(int electrodeNum, int channelNum);
    /** */


    void setChannelActive(int electrodeIndex, int channelNum, bool active);

    /** Returns the number of channels for a given electrode. */
    int getNumChannels(int index);

    /** Edits the mapping between input channels and electrode channels. */
    void setChannel(int electrodeIndex, int channelNum, int newChannel);

    /** Returns the continuous channel that maps to a given
        electrode channel. */
    int getChannel(int index, int chan);

    /** Sets the name of a given electrode. */
    void setElectrodeName(int index, String newName);

    /** Removes an electrode with a given index. */
    bool removeElectrode(int index);

    void setSVDColumnsToUse(int dim);


    void process(AudioSampleBuffer& buffer, MidiBuffer& events, int& nSamples);

    /** Any variables used by the "process" function _must_ be modified only through
this method while data acquisition is active. If they are modified in any
other way, the application will crash. */
    void setParameter(int parameterIndex, float newValue);

    /** Called prior to start of acquisition. */
    bool enable();

    /** Called after acquisition is finished. */
    bool disable();

    float getElectrodeNoiseVar(int index);

    StringArray electrodeTypes;
    StringArray getElectrodeNames();

    bool allParametersEstimated;

    void saveCustomParametersToXml(XmlElement* parentElement);
    void loadCustomParametersFromXml();

    float acfLag1;

    Array<Electrode*> electrodes;

    TwoDimMatrix dictionary;

    SVDjob job;
    int SVDCols;

private:

    // private members and methods go here
    //
    // e.g.:
    //
    // float threshold;
    // bool state;

    AudioSampleBuffer& dataBuffer;
    AudioSampleBuffer overflowBuffer;

    int currentElectrode;
    int currentChannelIndex;
    int currentIndex;

    int sampleIndex;

    int overflowBufferSize;

    float delta, cssp;
    float sampleRate;

    int dictionarySize;



    Array<int> electrodeCounter;

    SVDcomputingThread dictionaryThread;

    bool useOverflowBuffer;

    float getDefaultThreshold();

    bool SVDmethod;


    int64 timestamp;

    Array<SpikeObject> detectedSpikesAllElectrodes;
    //Array<Array<SpikeObject>> detectedSpikesPerElectrode;

    void addWaveformToSpikeObject(SpikeObject* s,
                                  int& peakIndex,
                                  int& electrodeNumber,
                                  int& currentChannel);
    void handleEvent(int eventType, MidiMessage& event, int sampleNum);

    void resetElectrode(Electrode*);


    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(ParameterEstimator);

};







#endif  // __PARAMETERESTIMATION_H_1B5A8C4B__
