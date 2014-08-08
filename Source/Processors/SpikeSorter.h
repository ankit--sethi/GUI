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
        void show(int index, Eigen::VectorXd& returnedX);
        void dequeue();
        bool isBufferPlush(int minsize);
        int getsize();
    private:
        int front, back, count, size;
        float *array;

};


class SpikeSorter : public GenericProcessor

{
public:

    SpikeSorter();
    ~SpikeSorter(){ }
    bool checkIfLogIsInf(float  num, float den);

    float P; //column span of the dictionary
    float alpha;
    float kappa0;
    float nu0;
    unsigned int K;
    //Array<Array<float>> phi0;
    Eigen::MatrixXd phi0;
    float apii;
    float bpii;
    float beta;
    float tau;
    unsigned int Cmax;  //maximum possible number of neurons present
    int curndx; //used to index current location in buffer
    int lookahead; //functionally, it is the size of the buffer
    int range;
    float pii;
    //Array <float> nu;
    Eigen::RowVectorXd nu;
    juce::Array <Eigen::MatrixXd> phi;
    juce::Array <Eigen::MatrixXd> lamclus;
    juce::Array <Eigen::MatrixXd> lamclusinv;
    juce::Array <Eigen::MatrixXd> R;
    juce::Array <Eigen::MatrixXd> Rinv;
    Eigen::MatrixXd Qt;
    Eigen::MatrixXd muu0;
    Eigen::RowVectorXd kappa;
    Eigen::MatrixXd muu;
    Eigen::MatrixXd lambda;
    double logDeterminantOfLambda;
    Eigen::MatrixXd sigma;
    juce::Array<Eigen::MatrixXd> Q;
    //Eigen::MatrixXd Q;
    Eigen::MatrixXd Qinv;
    Eigen::MatrixXd ReducedDictionary; //This is A (PxK)
    Eigen::MatrixXd ReducedDictionaryTranspose;
    Eigen::MatrixXd Qmat;
    Eigen::VectorXd yhat;
    Eigen::VectorXd mhat;
    Eigen::MatrixXd Qhat;
    Eigen::VectorXd xwindLonger;
    unsigned int neuronCount; // this is C in the code

    double Hadj;
    double logdetQ;
    double logPlusDetTermForNoiseLL;
    double logPlusDetTermForNewNeuronLL;
    int idx; //this is Q in the second half
    float likelihoodThreshold;

    bool suppress;
    bool test;
    Eigen::RowVectorXd tlastspike;
    Eigen::RowVectorXd ltheta;
    Eigen::RowVectorXd ngamma;

    float threshold;
    int totalSpikesFound;


    //Array <TwoDimMatrix> lamclusS;

    float suppresslikelihood;
    int deltaT;

    int64 timestamp;

    juce::Array <CircularQueue*> circbuffer;

    Eigen::VectorXd xwind; //remember, xwind is "of the moment", xwind lacks history
    Eigen::VectorXd xwindloop;

    bool setLambda;

    bool isSource()
    {
        return false;
    }


    bool isSink()
    {
        return false;
    }

    bool checkIfAllParametersEstimated;
    void collectSamplesForSpikeObject(int electrodeIndex, float trigSample);
    void addNewSampleAndLikelihoodsToCurrentSpikeObject(float sample, MidiBuffer &eventBuffer, int chan);
    void PackageCurrentSortedSpikeIntoBuffer(MidiBuffer &eventBuffer1);
    void updateAllSortingParameters();
    void process(AudioSampleBuffer& buffer, MidiBuffer& events, int& nSamples);
    float getNextSample(int& chan);

    void setParameter(int parameterIndex, float newValue);
    void handleEvent(int eventType, MidiMessage& event, int sampleNum);
    bool enable();
    bool disable();
    AudioProcessorEditor* createEditor();
    void updateSettings();

private:

    AudioSampleBuffer& dataBuffer;
    AudioSampleBuffer nullbuffer;
    //Array<Electrode*> electrodes;
    ParameterEstimator* node;
    unsigned long int sampleIndex;
    unsigned long int masterSampleIndex;

    bool paramsCopied;

    SortedSpikeObject currentSpike;
    int currentIndex;  // this is for keeping track of indices after threshold crossing

    Eigen::VectorXd lthr;
    Eigen::MatrixXd lon;

    bool buffersArePlush;

    bool spikeDetected;
    bool samplesBeingCollected;

    int number;
    double startTime;
    double stopTime;

    float likelihoodNoSpike;
    float maximumLikelihoodPerNeuron;
    Eigen::RowVectorXd likelihoodPerNeuron;
    double Hadjsum;
    bool thingsHaveChanged;

    Eigen::HouseholderQR<Eigen::MatrixXd> lambdaQR;
    Eigen::HouseholderQR<Eigen::MatrixXd> QQR;
    double Quad;

    uint8_t* spikeBuffer;
    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(SpikeSorter);

};


#endif  // __SPIKESORTER_H_4BC9FE8__
