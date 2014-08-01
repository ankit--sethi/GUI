/*
  ==============================================================================

    SpikeSorter.cpp
    Created: 27 May 2014 1:34:38pm
    Author:  sethisan

  ==============================================================================
*/


//#include "SpikeDetector.h"
//#include "ParameterEstimation.h"
#include "Channel.h"
#include "SpikeSorter.h"
#include <math.h>
#include <limits.h>
#include <fstream>


using namespace std;


#define PI 3.1415926535897932384626433832795


// calculate the cofactor of element (row,col)

int CircularQueue::getsize()
{
    return size;
}

void CircularQueue::setsize(int siz)
{
    front = back = -1;
    size = siz - 1;
    array = new float[size];
    count = 0;

    for(int i = 0; i <= size; i++)
    {
     array[i] = 0.0;
    }
    //std::cout << size << std::endl;
}

void CircularQueue::enqueue(float item)
{
    if (front == 0 && back == size || front == back + 1) {
        //cout << "Queue is full\n";
    }
    else if (front == -1 && back == -1)
    {
        front = 0;
        back = 0;
        array[front] = item;
        count++;
    }
    else if (back == size)
    {
        back = 0;
        array[back] = item;
        count++;
    }
    else
    {
        back++;
        array[back] = item;
        count++;
    }
}

void CircularQueue::dequeue()
{
    if (front == -1 && back == -1)
    {
        std::cout << "Queue is empty\n";
    }
    else {
        if (front == back)
        {
            array[front] = 0;
            front = -1;
            back = -1;
            count--;
        }
        else if (front == size)
        {
            array[front] = 0;
            front = 0;
            count--;
        }
        else
        {
            array[front] = 0;
            front++;
            count--;
        }
    }
}

juce::Array<float> CircularQueue::show(int index)
{
    juce::Array <float> xwithPsamples;
    if (index <= size + 1)
    {
        for(int i = 0; i < index; i++)
            xwithPsamples.add(array[i]);

        for (int i = 0; i < xwithPsamples.size(); i++)
        {
            //std::cout << xwithPsamples[i];
        }
    }
    return xwithPsamples;

}
bool CircularQueue::isBufferPlush(int minsize)
{
    if (count >= minsize)
        return true;
    else
        return false;
}

SpikeSorter::SpikeSorter()
    : GenericProcessor("Spike Sorter"), nullbuffer(2,100), dataBuffer(nullbuffer), stopTime(0), startTime(0) //, threshold(200.0), state(true)
{

    // THIS CONSTRUCTOR TAKES CARE OF EVERYTHING EXCEPT ParameterEstimator* node and Array<CircularBuffer> circbuffer
    // THEY ARE HANDLED ONCE SpikeSorter::process() STARTS AND CALLS SpikeSorter::setParameter()

    alpha = 1e-1;
    kappa0 = 0.01;
    nu0 = 0.1;
    buffersArePlush = false;


    K = 3; // number of SVD components to use; need to make this a control
    phi0 = (Eigen::MatrixXd::Identity(K,K).array()*0.1).matrix();
    //check phi0!!!
    number = 0;
    apii = 1;
    bpii = 1e7;
    beta = 1/(30*(40000));
    //std::cout << " Sample RATE is " << getSampleRate() << "//";
    tau = 10;

    P = 40;

    checkIfAllParametersEstimated = false;
    paramsCopied = false;


    Cmax = 50;  //maximum possible number of neurons present
    curndx = 0; //used to index current location in buffer
    lookahead = 500; //functionally, it is the size of the buffer
    range = 20; // defines basically the width of a spike
    // setting sigma
    sigma = Eigen::MatrixXd::Zero(P,P);

    pii=apii/bpii;

    Eigen::MatrixXd phi0nu0forlamclus = ((phi0.inverse().array()*nu0)).matrix();
    Eigen::MatrixXd phi0nu0forlamclusinv = phi0nu0forlamclus;
    Eigen::MatrixXd phi0nu0forR = ((phi0.inverse().array()*nu0)*0.2).matrix();
    Eigen::MatrixXd phi0nu0forRinv = phi0nu0forR;
    nu = Eigen::RowVectorXd::Constant(1,Cmax, nu0);
    kappa = Eigen::RowVectorXd::Constant(1,Cmax, kappa0);
    ltheta = Eigen::RowVectorXd::Zero(1,Cmax);
    ngamma = Eigen::RowVectorXd::Zero(1,Cmax);
    tlastspike = Eigen::RowVectorXd::Zero(1,Cmax);
    likelihoodPerNeuron = Eigen::RowVectorXd::Zero(1,Cmax);

    for (int i = 0; i < Cmax; i++)
    {
        phi.add(phi0);
        lamclus.add(phi0nu0forlamclus);
        lamclusinv.add(phi0nu0forlamclusinv);
        R.add(phi0nu0forR);
        Rinv.add(phi0nu0forRinv);
    }
    muu = Eigen::MatrixXd::Constant(K, Cmax, 0);
    muu0 = Eigen::MatrixXd::Constant(K, Cmax, 0);
    Q = Eigen::MatrixXd::Constant(P, P, 0);
    Qinv = Eigen::MatrixXd::Constant(P, P, 0);
    Qhat = Eigen::MatrixXd::Constant(K, K, 0);
    Qmat = Eigen::MatrixXd::Constant(K, K, 0);
    mhat = Eigen::VectorXd::Constant(K, 1, 0);
    yhat = Eigen::VectorXd::Constant(K, 1, 0);
    ReducedDictionary = Eigen::MatrixXd::Constant(P, K, 0);
    ReducedDictionaryTranspose = Eigen::MatrixXd::Constant(K, P, 0);

    neuronCount = 0; // This is 'C' in the MATLAB code.
    //nz = 0; // Not sure what this does. Possible the same thing.
    setParameter(1,0);
    xwind = Eigen::VectorXd::Zero(P);

    threshold = std::log(pii/(1-pii));
    std::cout<< "Threshold is " << threshold;

    suppresslikelihood = 1e5;
    suppress = false;
    spikeDetected = false;
    samplesBeingCollected = false;
    likelihoodThreshold = 0;
    idx = 0;
    Hadj = 0;
    setLambda = true;
    masterSampleIndex = 0;
    totalSpikesFound = 0;
    spikeBuffer = new uint8_t[MAX_SORTED_SPIKE_BUFFER_LEN];
}

void SpikeSorter::handleEvent(int eventType, MidiMessage& event, int sampleNum)
{

    if (eventType == TIMESTAMP)
    {
        const uint8* dataptr = event.getRawData();

        memcpy(&timestamp, dataptr + 4, 8); // remember to skip first four bytes
    }


}

bool SpikeSorter::enable()
{

    SpikeSorterEditor* editor = (SpikeSorterEditor*) getEditor();
     editor->enable();

    return true;
}

bool SpikeSorter::disable()
{

}

AudioProcessorEditor* SpikeSorter::createEditor()
{
    editor = new SpikeSorterEditor(this, true);
    return editor;
}

void SpikeSorter::updateSettings()
{

}

bool SpikeSorter::checkIfLogIsInf(float  num, float den)
{
    if (num == 0)
        return true;
    if (den == 0)
        return true;

    return false;
}

void SpikeSorter::setParameter(int parameterIndex, float newValue)
{

    if (parameterIndex == 1 )
    {

        for (int i = 0; i < ltheta.size(); i++)
        {
            if (!checkIfLogIsInf(ngamma[i],(alpha+totalSpikesFound)))
            {
                ltheta(i) = (log(ngamma[i]/(alpha+totalSpikesFound)));

            }
            else
            {
                ltheta(i) = -1*PI;
            }

        }
        ltheta(neuronCount) = log(alpha/(alpha+totalSpikesFound));

    }

}

float SpikeSorter::getNextSample(int& chan)
{
        if (sampleIndex < dataBuffer.getNumSamples())
            return *dataBuffer.getSampleData(chan, sampleIndex);
        else
            return 0;
}

int SpikeSorter::findNeuronID()
{
        float min = 0;
        idx = 0;
        min = lthr(0);
        for (int i= 0; i < lthr.size(); i++)
        {
            if (min < lthr(i))
            {
                min = lthr(i);
                idx = i;
            }
        }
        float max = 0, idx1 = 0;
        max = lon(idx, 0);
        for (int i= 0; i < lon.row(idx).size(); i++)
        {
            if (max < lon(idx, i))
            {
                max = lon(idx, i);
                idx1 = i;
            }
        }

        if (idx1 == neuronCount)
        {
            neuronCount = idx1;
        }
        return idx1;
}

void SpikeSorter::PackageCurrentSortedSpikeIntoBuffer(MidiBuffer& eventBuffer1)
{
    //std::cout << "ohh noes";
    int numBytes = packSortedSpike(currentSpike, spikeBuffer, MAX_SPIKE_BUFFER_LEN);

    //std::cout << "reached heres" << std::endl;

    if (numBytes > 0)
        eventBuffer1.addEvent(spikeBuffer, numBytes, int(currentSpike.timestamp));

}


void SpikeSorter::collectSamplesForSpikeObject(int electrodeIndex, float trigSample)
{
    //cout<<"reached here!!!!11111vvvvvvvvv";
    samplesBeingCollected = true;
    currentIndex = 0;
    currentSpike.eventType = SORTEDSPIKE;
    //std::cout<< "The EVENTTYPE IS" << int(currentSpike.eventType) << std::endl;

    currentSpike.nChannels = 1;//node->electrodes[electrodeIndex]->numChannels;
    currentSpike.nSamples = P;
    currentSpike.source = electrodeIndex;
    int chan = *(node->electrodes[electrodeIndex]->channels+ 0); // HARDCODING ONE CHANNEL FOR NOW
    for (int i = 0; i < MAX_NUMBER_OF_SPIKE_CHANNELS; i++)
    {
    currentSpike.threshold[i] = threshold;
    currentSpike.gain[i] = (int)(1.0f / channels[chan]->bitVolts)*1000;
    }

    currentSpike.timestamp = masterSampleIndex + sampleIndex;
    //std::cout << "The timestamp is = " << timestamp;
    lthr = Eigen::VectorXd::Zero(range + 1);
    lon = Eigen::MatrixXd::Zero(range + 1, neuronCount + 1);
    xwindLonger = Eigen::VectorXd::Zero(xwind.size() + range);
    xwindLonger.head(xwind.size()) = xwind;
    //cout<<"reached here!!!!11111vvvvvvvvv";
    lthr(currentIndex) = likelihoodThreshold;
    for (int i = 0; i <= neuronCount; i++)
    {
        lon(currentIndex, i) = likelihoodPerNeuron(i);
    }
    currentIndex++;
}

void SpikeSorter::addNewSampleAndLikelihoodsToCurrentSpikeObject(float sample, MidiBuffer& eventBuffer, int chan)
{
    if (currentIndex <= range)
    {
        lthr(currentIndex) = likelihoodThreshold;
        //cout<< "reached till just before lon";
        for (int i = 0; i <= neuronCount; i++)
        {
            lon(currentIndex, i) = likelihoodPerNeuron(i);
        }
        //cout << xwind.size() << " is the xwindsize " << "and the currentIndex is" << currentIndex << "// and the xwindlonger size is " << xwindLonger.size() << "//";
        xwindLonger(xwind.size() + currentIndex - 1) = sample;
        currentIndex++;
    }
    else
    {
        //cout<< "reached start of else";
        spikeDetected = false;
        samplesBeingCollected = false;
        double minthr = lthr.minCoeff(&idx);
        int Cnew;
        double maxlon = lon.row(idx).maxCoeff(&Cnew);
        if (Cnew > neuronCount)
        {
            neuronCount = Cnew;
        }
        currentSpike.neuronID = Cnew;
        currentIndex = 0;
        //cout << "The INDEX IS //" << idx << "//";
        for (int i = 0; i < P; i ++)
        {
        currentSpike.data[i] = uint16(xwindLonger(idx + i)/ channels[chan]->bitVolts + 32768);
        }
        PackageCurrentSortedSpikeIntoBuffer(eventBuffer);

        //std:: cout << int(currentSpike.eventType) << "//" << currentSpike.gain[1] << "//" << currentSpike.nChannels << "//" << currentSpike.neuronID << "//" << currentSpike.nSamples << "//" <<
        //currentSpike.source << "//" << currentSpike.threshold[1] << "//" << currentSpike.timestamp << std::endl;

        int numBytes = packSortedSpike(currentSpike, spikeBuffer, MAX_SPIKE_BUFFER_LEN);
        if (numBytes > 0)
        eventBuffer.addEvent(spikeBuffer, numBytes, int(currentSpike.timestamp));

        //clear everything that needs to be cleared now! we look for new spikes
        cout<< "reached till just before sorting parameters";
        updateAllSortingParameters();
        std::cout << "One Spike Handled!" << std::endl;
        totalSpikesFound += 1;

    }
}


void SpikeSorter::updateAllSortingParameters()
{
    int neuronID = currentSpike.neuronID;
    Qmat = ReducedDictionaryTranspose*lambda*ReducedDictionary + lamclus[neuronID];
    yhat = Qmat.inverse()*(ReducedDictionaryTranspose*lambda*(xwindLonger.segment(idx, P)) + lamclus[neuronID]*muu.col(neuronID));
    ngamma(neuronID) = ngamma(neuronID) + 1;
    deltaT = masterSampleIndex + sampleIndex -P - range + idx - tlastspike(neuronID);
    tlastspike(neuronID) = masterSampleIndex + sampleIndex - P - range + idx;
    double ebet = std::exp(-1*double(beta)*double(deltaT));
    mhat = (muu0.col(neuronID).array()*(1-ebet) + muu.col(neuronID).array()*ebet).matrix();
    muu0.col(neuronID) = ((kappa(neuronID)*muu0.col(neuronID).array() + yhat.array())/(kappa(neuronID)+1)).matrix();
    Qhat = ((1/tau)*Eigen::MatrixXd::Identity(K,K).array()*(1-ebet*ebet) + R[neuronID].inverse().array()*(ebet*ebet)).matrix();
    R.setUnchecked(neuronID, Qhat.inverse() + lamclus[neuronID]);
    muu.col(neuronID) = R[neuronID].inverse()*(Qhat.inverse()*(mhat) + lamclus[neuronID]*yhat);
    Eigen::MatrixXd temp = ReducedDictionaryTranspose*lambda*xwindLonger.segment(idx,P) + lamclus[neuronID]*muu.col(neuronID);
    yhat = Qmat.inverse()*(temp);
    double constant = (kappa(neuronID)/(kappa(neuronID) + 1));
    Eigen::MatrixXd trans = (yhat - muu.col(neuronID))*((yhat - muu.col(neuronID)).transpose());
    Eigen::MatrixXd newphi = (constant*trans.array()).matrix();
    Eigen::MatrixXd temp1 = newphi + phi[neuronID] + Qmat.inverse();
    phi.setUnchecked(neuronID, temp1);
    kappa(neuronID) = kappa(neuronID) + 1;
    nu(neuronID) = nu(neuronID) + 1;
    lamclus.setUnchecked(neuronID, (phi[neuronID].inverse().array()*nu(neuronID)).matrix());
}


void SpikeSorter::process(AudioSampleBuffer& buffer,
                          MidiBuffer& events,
                          int& nSamples)
{

    int activeChannels = 0;
    dataBuffer = buffer;


    if (setLambda)
    {

        ProcessorGraph* gr = getProcessorGraph();
        juce::Array<GenericProcessor*> p = gr->getListOfProcessors();

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
        }

        CircularQueue* queue = new CircularQueue[node->electrodes.size()];
        circbuffer.clear();

        for ( int i = 0; i < node->electrodes.size(); i++ )
        {
            circbuffer.add((queue+i));
            circbuffer[i]->setsize(P);
        }
        setLambda = false;
        std::cout<<"Parameter Estimator contacted and Circular Buffers setup and the buffer size and electrode number is " << circbuffer[0]->getsize() << "  and   " << circbuffer.size() << std::endl;

    }



    //checkForEvents(events);

    if(node->allParametersEstimated)
    {
        // spike sort waveform using Vogelstein code
        Logger *log1 = Logger::getCurrentLogger();
        startTime = Time::getMillisecondCounterHiRes();

        if(!paramsCopied)
        {
            std::cout<<"Entered Params copied if condition" << std::endl;


            //Array<float>* test = sigma.getRawDataPointer();
            Eigen::ArrayXd powers = Eigen::VectorXd::Zero(P).array();
            double start = 1;

            for (int i = 0; i < P; i++)
            {
                powers(i) = start;
                start *= node->acfLag1;
            }

            float powerIndex;
            for (int i = 0; i < P; i++)
            {
                powerIndex = -1*i;
                for (int j = 0; j < P; j++)
                {

                    if ( i == j )
                        sigma(i,j) = node->getElectrodeNoiseVar(0);
                    else
                        sigma(i,j) = node->getElectrodeNoiseVar(0)*powers(std::abs(powerIndex));
                    powerIndex++;
                }
            }

            std::cout<< " ACF LAG in Spike Sorter is " << node->acfLag1 << " and the Var is " << node->getElectrodeNoiseVar(0) << "  "  << std::endl;
            double startTime1 = Time::getMillisecondCounterHiRes();
            lambda = sigma.inverse();
            double stopTime1 = Time::getMillisecondCounterHiRes();
            log1->writeToLog("Time for 40 x 40 inverse is = " + String(stopTime1 - startTime1));
            Eigen::HouseholderQR<Eigen::MatrixXd> lambdaQR(lambda);
            logDeterminantOfLambda = lambdaQR.logAbsDeterminant();
            std::cout<<"reached line 3 and det is = " << logDeterminantOfLambda << " and " << ReducedDictionary.col(0).size() << "x" << ReducedDictionary.row(0).size() << "//" << std::endl;

            for (int i = 0; i < P; i++)
            {
                for (int j = 0; j < K; j++)
                {
                    ReducedDictionary(i,j) = node->dictionary.get(i,j);
                }
            }
            ReducedDictionaryTranspose = ReducedDictionary.transpose();
            paramsCopied = true;
            std::cout<<"Params copied and nSamples is" << nSamples << std::endl;
        }

        for (int i = 0; i < (*node).electrodes.size(); i++)
        {
            sampleIndex = 0;
            activeChannels = 0;

            for (int chan = 0; chan < 1; chan++) // right now, ONLY ONE CHANNEL
            {
                if ( *(node->electrodes[i]->isActive +chan) )
                {
                    activeChannels++;
                    while (sampleIndex < nSamples)
                    {
                        int currentChannel = *(node->electrodes[i]->channels+chan);

                        float sample = getNextSample(currentChannel);
                        sampleIndex++;

                        circbuffer[i]->enqueue(sample);

                        pii = (apii + totalSpikesFound)/(bpii + masterSampleIndex);
                        threshold = std::log(pii/(1-pii));


                        if(circbuffer[i]->isBufferPlush(P))
                        {

                            // find likelihood of no spike
                            Array <float> temp = circbuffer[i]->show(P);
                            for (int i = 0; i < temp.size(); i++)
                            {
                                    xwind(i) = temp[i];
                            }

                            xwindloop = xwind;
                            for (int i = 0; i < ltheta.size(); i++)
                            {
                                if (!checkIfLogIsInf(ngamma(i),(alpha+totalSpikesFound)))
                                {
                                    ltheta(i) = (std::log(ngamma(i)/(alpha+totalSpikesFound)));

                                }
                                else
                                {
                                    ltheta(i) = -1*PI;
                                }

                            }
                            ltheta(neuronCount) = std::log(alpha/(alpha+totalSpikesFound));

                            double Quad = (xwind.dot(lambda*xwind));
                            likelihoodNoSpike = (-1*P/2)*std::log(2*PI) + 0.5*logDeterminantOfLambda - 0.5*Quad;

                            for (int j = 0; j < neuronCount; j++)
                            {
                                Q = sigma + ReducedDictionary*(R[neuronCount].inverse() + lamclus[neuronCount].inverse())*(ReducedDictionaryTranspose);
                                xwindloop = xwind - ReducedDictionary*muu.col(j);
                                Eigen::HouseholderQR<Eigen::MatrixXd> QQR(Q);
                                double sum = QQR.logAbsDeterminant();
                                cout << xwindloop*(Q.householderQr().solve(xwindloop));
                                if( ( (masterSampleIndex + sampleIndex) - tlastspike(j)) < 40000*50/10000 )  // THIS NEEDS TO BE INVESTIGATED
                                    likelihoodPerNeuron(j) = -(0.5*P)*std::log(2*PI) - sum - suppresslikelihood - 0.5*(xwindloop.dot(Q.inverse()*xwindloop));
                                else
                                    likelihoodPerNeuron(j) = -(0.5*P)*std::log(2*PI) - sum - 0.5*(xwindloop.dot(Q.inverse()*xwindloop));
                            }

                            Q = sigma + ReducedDictionary*(R[neuronCount].inverse() + lamclus[neuronCount].inverse())*(ReducedDictionaryTranspose);
                            Eigen::HouseholderQR<Eigen::MatrixXd> QQR(Q);
                            double sum = QQR.logAbsDeterminant();
                            likelihoodPerNeuron(neuronCount) = -(0.5*P)*std::log(2*PI) - sum - 0.5*(xwind.dot(Q.inverse()*xwind));
                            //cout<< R[neuronCount];
                            //cout<< "BA" <<  lamclus[neuronCount];
                            for (int i = 0; i <= neuronCount; i++ )
                            {
                                likelihoodPerNeuron(i) =  likelihoodPerNeuron(i) + ltheta(i);
                                //std::cout << "And " << likelihoodPerNeuron[i] << " is the Likelihood per Neuron for " << i << " and " << " ltheta is " << ltheta[i].getVar() << std::endl;
                            }

                            for (int i = 0; i <= neuronCount; i++ )
                            {
                                if (i == 0)
                                    {maximumLikelihoodPerNeuron = likelihoodPerNeuron(i);}
                                else if(maximumLikelihoodPerNeuron < likelihoodPerNeuron(i))
                                    maximumLikelihoodPerNeuron = likelihoodPerNeuron(i);
                            }
                            for (int i = 0; i <= neuronCount; i++ )
                            {
                                likelihoodPerNeuron(i) = likelihoodPerNeuron(i) - maximumLikelihoodPerNeuron;
                            }
                            sum = 0;
                            for (int i = 0; i <= neuronCount; i++ )
                            {
                                sum += std::exp(likelihoodPerNeuron(i));
                            }

                            Hadj = std::log(sum);

                            likelihoodThreshold = likelihoodNoSpike - maximumLikelihoodPerNeuron - Hadj;
                            //cout<< Hadj;
                            //std::cout<< "Likelihood Threshold is = " << likelihoodThreshold << " and the threshold is " << threshold << std::endl;
                            if (likelihoodThreshold < threshold)
                                spikeDetected = true;

                            if(spikeDetected)
                            {
                                if (samplesBeingCollected)
                                {
                                    //std::cout<< "came here to samplesBeingCollected " << std::endl;
                                    addNewSampleAndLikelihoodsToCurrentSpikeObject(sample, events, chan);
                                }
                                else
                                {
                                    //std::cout<< "threshold entered" << std::endl;
                                    collectSamplesForSpikeObject(i, sample);
                                }
                            }

                            circbuffer[i]->dequeue(); // After everything is done, remove most past sample, ie, xwind(0)
                        }

                    }
                }

            }
        }
        masterSampleIndex += sampleIndex;

        stopTime = Time::getMillisecondCounterHiRes();
        //log1->writeToLog("Duration taken was " + String(stopTime - startTime) + "ms for nSamples = " + String(nSamples));

    }


}
