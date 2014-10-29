/*
  ==============================================================================

    SpikeSorter.h
    Created: 27 May 2014 1:34:38pm
    Author:  sethisan

  ==============================================================================
*/

#ifndef __SPIKESORTER_H_4BC9FE8__
#define __SPIKESORTER_H_4BC9FE8__


#ifdef _WIN32
#include <Windows.h>
#endif

#include "../../JuceLibraryCode/JuceHeader.h"
#include "GenericProcessor.h"
#include "ParameterEstimation.h"
#include "SourceNode.h"
//#include "SpikeDetector.h"
#include "Editors/SpikeSorterEditor.h"
#include "../../Resources/eigen-eigen-6b38706d90a9/Eigen/Dense"
#include "../../Resources/eigen-eigen-6b38706d90a9/Eigen/QR"

#define MAX_CHANNELS 4
#define MAX_ELECTRODES 16
//using namespace Eigen;
//using namespace juce;

class expandedfloat
{
public:
    expandedfloat()
    {
    }
    ~expandedfloat()
    {

    }
    bool isInf()
    {
        return infFlag;
    }
    void setInf()
    {
        infFlag = true;
    }
    void setNonInf()
    {
        infFlag = false;
    }
    void setVar(float element)
    {
        var = element;
    }
    float getVar()
    {
        if (infFlag == false)
        {
            return var;
        }
        else
        {
            return 0.0;
        }
    }

private:
    float var;
    bool infFlag;
};


class ContinuousCircularBuffer
{
public:
    ContinuousCircularBuffer(int NumCh, float SamplingRate, int SubSampling, float NumSecInBuffer);
    void reallocate(int N);
    void update(AudioSampleBuffer& buffer, int64 hardware_ts, int64 software_ts, int numpts);
    int GetPtr();
    void addTrialStartToSmartBuffer(int trialID);
    int numCh;
    int subSampling;
    float samplingRate;
    CriticalSection mut;
    int numSamplesInBuf;
    float numTicksPerSecond;
    int ptr;
    int bufLen;
    int leftover_k;
    std::vector<std::vector<float> > Buf;
    std::vector<bool> valid;
    std::vector<int64> hardwareTS,softwareTS;
};


class CircularQueue
{
    public:
        CircularQueue() { size = 0;}
        ~CircularQueue() { }
        void setsize(int siz);
        void enqueue(float item);
        void show(int index, Eigen::Matrix<float,30,1>& returnedX);
        void dequeue();
        bool isBufferPlush(int minsize);
        int getsize();
    private:
        int front, back, count, size;
        float *array;

};
class SpikeSorter;
class ParameterEstimator;

class SortingElectrodeChannel
{
    public:
    juce::Array <Eigen::Matrix3f> phi;
    juce::Array <Eigen::Matrix3f> lamclus;
    juce::Array <Eigen::Matrix3f> lamclusinv;
    juce::Array <Eigen::Matrix3f> R;
    juce::Array <Eigen::Matrix3f> Rinv;
    juce::Array<Eigen::Matrix<float, 30, 30>> Q;
    juce::Array<float> detQ;
    Eigen::MatrixXf muu0;
    Eigen::MatrixXf muu;
    CircularQueue circbuffer;
    Eigen::Matrix<float, 30, 1> xwind; //remember, xwind is "of the moment", xwind lacks history
    Eigen::Matrix<float, 30, 1> xwindloop;
    Eigen::VectorXf xwindLonger;
    Eigen::HouseholderQR<Eigen::MatrixXf> lambdaQR;
    juce::Array<Eigen::HouseholderQR<Eigen::MatrixXf>> QQR;
    juce::Array<Eigen::LLT<Eigen::Matrix<float,30,30>>> QLLT;
    Eigen::VectorXf QLogAbsDeterminant;
    Eigen::VectorXf cLL;
};

class SortingElectrode
{
    public:
    Eigen::RowVectorXf nu;
    Eigen::RowVectorXf kappa;
    Eigen::Matrix<float, 30, 30> lambda;
    Eigen::Matrix<float, 30, 30> sigma;
    Eigen::Matrix<float, 30, 3> ReducedDictionary;
    Eigen::Matrix<float, 3, 30> ReducedDictionaryTranspose;
    Eigen::Matrix3f Qmat;
    Eigen::VectorXf yhat;
    Eigen::VectorXf mhat;
    Eigen::MatrixXf Qhat;
    unsigned int neuronCount;

    // --------------- temporary variables

    Eigen::Matrix3f RTLR;
    Eigen::Matrix<float, 3, 30> RTL;
    Eigen::Matrix<float, 30, 3> LR;
    Eigen::Matrix3f Qmatinv;
    Eigen::Matrix3f C;
    Eigen::Matrix3f bracket;
    Eigen::Matrix<float, MAX_CHANNELS, 1> sample;
    Eigen::HouseholderQR<Eigen::MatrixXf> lambdaQR;
    float detSigma;
    bool flag;
    // ---------------

    double Hadj;
    float logPlusDetTermForNoiseLL;
    int idx;
    float likelihoodThreshold;
    bool suppress;
    bool test;
    Eigen::RowVectorXf tlastspike;
    Eigen::RowVectorXf ltheta;
    Eigen::RowVectorXf ngamma;
    float threshold;
    float totalSpikesFound;
    int deltaT;
    bool setLambda;
    Eigen::VectorXf lthr;
    Eigen::MatrixXf lon;
    SortedSpikeObject currentSpike;
    int currentIndex;
    bool spikeDetected;
    bool samplesBeingCollected;
    float likelihoodNoSpike;
    float maximumLikelihoodPerNeuron;
    Eigen::RowVectorXf likelihoodPerNeuron;
    double Hadjsum;
    bool thingsHaveChanged;
    uint8_t* spikeBuffer;

    SortingElectrodeChannel sec[MAX_CHANNELS];
};


class SpikeSorter : public GenericProcessor

{
public:

    SpikeSorter();
    ~SpikeSorter();
    bool checkIfLogIsInf(float  num, float den);
    int P;
    float alpha;
    float kappa0;
    float nu0;
    unsigned int K;
    Eigen::Matrix3f phi0;
    float apii;
    float bpii;
    float beta;
    float tau;
    unsigned int Cmax;
    int curndx;
    int lookahead;
    int range;
    float pii;
    float suppresslikelihood;
    SortingElectrode se[MAX_ELECTRODES];
    int64 timestamp;
    bool flag;

    bool isSource()
    {
        return false;
    }


    bool isSink()
    {
        return false;
    }

    bool checkIfAllParametersEstimated;
    void collectSamplesForSpikeObject(int electrodeIndex, int index);
    void addNewSampleAndLikelihoodsToCurrentSpikeObject(Eigen::Matrix<float, 4, 1> sample, MidiBuffer &eventBuffer, int electrodeIndex, int channelCount);
    void PackageCurrentSortedSpikeIntoBuffer(MidiBuffer &eventBuffer1, int electrodeIndex);
    void updateAllSortingParameters(int electrodeIndex, int chan);
    void process(AudioSampleBuffer& buffer, MidiBuffer& events, int& nSamples);
    float getNextSample(int& chan);

    void setParameter(int parameterIndex, float newValue);
    void handleEvent(int eventType, MidiMessage& event, int sampleNum);
    bool enable();
    bool disable();
    AudioProcessorEditor* createEditor();
    void updateSettings();
    bool allParametersEstimated;

private:

    AudioSampleBuffer& dataBuffer;
    AudioSampleBuffer nullbuffer;
    //Array<Electrode*> electrodes;
    ParameterEstimator* node;
    unsigned long int sampleIndex;
    float masterSampleIndex;
    bool timeflag;
    float timetaken;
    float n;
    bool buffersArePlush;
    int number;
    double startTime;
    double stopTime;
    int downsamplingFactor;

    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(SpikeSorter);
};


#endif  // __SPIKESORTER_H_4BC9FE8__
