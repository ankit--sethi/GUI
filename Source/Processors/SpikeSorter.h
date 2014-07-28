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
#include "../../Resources/eigen-eigen-6b38706d90a9/Eigen/LU"

using namespace Eigen;

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
        Array<float> show(int index);
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
    MatrixXd phi0;
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
    VectorXd nu;
    juce:Array <MatrixXd> phi;
    juce:Array <MatrixXd> lamclus;
    juce:Array <MatrixXd> lamclusinv;
    juce:Array <MatrixXd> R;
    juce:Array <MatrixXd> Rinv;
    MatrixXd Qt;
    MatrixXd muu0;
    VectorXd kappa;
    MatrixXd muu;
    MatrixXd lambda;
    double logDeterminantOfLambda;
    MatrixXd sigma;
    MatrixXd Q;
    MatrixXd Qinv;
    MatrixXd ReducedDictionary; //This is A (PxK)
    MatrixXd ReducedDictionaryTranspose;
    MatrixXd Qmat;
    MatrixXd yhat;
    MatrixXd mhat;
    MatrixXd Qhat;
    VectorXd xwindLonger;
    unsigned int neuronCount; // this is C in the code

    double Hadj;
    int idx; //this is Q in the second half
    float likelihoodThreshold;

    bool suppress;
    bool test;
    juce:Array <> tlastspike;
    juce:Array <expandedfloat> ltheta;
    VectorXd ngamma;

    float threshold;
    int totalSpikesFound;


    //Array <TwoDimMatrix> lamclusS;

    float suppresslikelihood;
    int deltaT;

    int64 timestamp;

    Array <CircularQueue*> circbuffer;

    float getaBaInnerProductSum(Array<float>& xwindm, Array<Array<float>>& InnerMatrix);

    juce:Array<float> xwind; //remember, xwind is "of the moment", xwind lacks history
    juce:Array<float> xwindloop;

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
    void invertSquareTwoDimMatrix(const Array<Array<float>> &A, Array<Array<float>>& result);
    void invertSquareTwoDimMatrix(const Array<Array<float> > &A, Array<Array<long double>>& result);
    void MatrixMultiplication(const Array<Array<float>>& A, const Array<Array<float>>& B, Array<Array<float> > &result);
    void MatrixTranspose(const Array<Array<float>>& A, Array<Array<float>>& result);
    void AddMatrices(const Array<Array<float>>& A, const Array<Array<float>>& B, Array<Array<float>>& result);
    void AddMatrixToItself(Array<Array<float>>& A, const Array<Array<float>>& B);
    void MatrixGetRowOrColumn(const Array<Array<float>>& A, bool row, int index, Array<Array<float>>& result);
    void returnMainDiagonal(const Array<Array<float>>& A, Array<float> &result);
    void cholesky(const Array<Array<float>>& A, Array<Array<float>> &result);
    void MatrixMultiplyWithConstant(const Array<Array<float>>& A, float c, Array<Array<float> > &result);
    void MatrixReplaceColumn(Array<Array<float> > &A, const Array<Array<float>>& B, int colindex);
    void eye(int dim, Array<Array<float> > &result);
    void dotSlashOperation(const Array<Array<float>>& A, float c, Array<Array<float> > &result);
    double getMatrixDeterminant(const Array<Array<float> > &A);
    double getMatrixDeterminant(const Array<Array<long double>>& A);
    double getMatrixDeterminant(const Array<Array<double> > &A);
    double getMatrixLogDeterminant(const Array<Array<float> > &A);
    void getLikelihoodPerNeuron(int FoundNeuronIndex, float tmp);
    void collectSamplesForSpikeObject(int electrodeIndex, float trigSample);
    void addNewSampleAndLikelihoodsToCurrentSpikeObject(float sample, MidiBuffer &eventBuffer, int chan);
    int findNeuronID();
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

    Array<float> lthr;
    Array<Array<float>> lon;

    void GetMinor(long double **src, long double **dest, int row, int col, int order);
    void GetMinor(double **src, double **dest, int row, int col, int order);
    double CalcDeterminant(long double **in_matrix, int n);
    double CalcDeterminant(double **in_matrix, int n);
    double CalcDeterminant(float **in_matrix, int n);
    void MatrixInversion(long double **A, int order, long double **Y);
    void MatrixInversion(double **A, int order, double **Y);

    bool buffersArePlush;

    bool spikeDetected;
    bool samplesBeingCollected;

    int number;
    double startTime;
    double stopTime;

    float likelihoodNoSpike;
    float maximumLikelihoodPerNeuron;
    Array <float> likelihoodPerNeuron;

    uint8_t* spikeBuffer;
    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(SpikeSorter);

};


#endif  // __SPIKESORTER_H_4BC9FE8__
