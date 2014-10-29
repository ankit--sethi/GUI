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
#define LOG2PIBY2 0.39908993417


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

void CircularQueue::show(int index, Eigen::Matrix<float, 30, 1> &returnedX) // caution: works right only if size of circular buffer = P
{
    int count  = 0, i = front;
    while (count < index)
    {
        returnedX(count) = array[i];
        i++;
        if (i > size)
        {
            i = 0;
        }
        count++;
    }
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
    cout << "Reached here";
    n = 0;
    timetaken = 0;
    timeflag = false;

    alpha = 1e-10;
    kappa0 = 1;
    nu0 = 0.1;
    buffersArePlush = false;


    K = 3; // number of SVD components to use; need to make this a control
    phi0 = (Eigen::Matrix3f::Identity(K,K).array()*0.1).matrix();
    //check phi0!!!
    number = 0;
    apii = 1;
    bpii = 1e7;
    beta = 1/(30*(10000));
    std::cout << " Sample RATE is " << getSampleRate() << "//";
    tau = 10;

    P = 30;

    checkIfAllParametersEstimated = false;

    Cmax = 50;  //maximum possible number of neurons present
    curndx = 0; //used to index current location in buffer
    lookahead = 500; //functionally, it is the size of the buffer
    range = 40; // defines basically the width of a spike

    pii=apii/bpii;

    Eigen::Matrix3f phi0nu0forlamclus = ((phi0.inverse().array()*nu0)).matrix();
    Eigen::Matrix3f phi0nu0forlamclusinv = phi0nu0forlamclus.inverse();
    Eigen::Matrix3f phi0nu0forR = ((phi0.inverse().array()*nu0)*0.2).matrix();
    Eigen::Matrix3f phi0nu0forRinv = phi0nu0forR.inverse();

    for (int i = 0; i < MAX_ELECTRODES; i++)
    {
        se[i].thingsHaveChanged = 1;
        se[i].nu = Eigen::VectorXf::Constant(Cmax, nu0);
        se[i].kappa = Eigen::VectorXf::Constant(Cmax, kappa0);
        se[i].ltheta = Eigen::VectorXf::Zero(Cmax);
        se[i].ngamma = Eigen::VectorXf::Zero(Cmax);
        se[i].tlastspike = Eigen::VectorXf::Zero(Cmax);
        se[i].likelihoodPerNeuron = Eigen::VectorXf::Zero(Cmax);
        se[i].Qhat = Eigen::MatrixXf::Constant(K, K, 0);
        se[i].Qmat = Eigen::Matrix3f::Constant(K, K, 0);
        se[i].mhat = Eigen::VectorXf::Constant(K, 1, 0);
        se[i].yhat = Eigen::VectorXf::Constant(K, 1, 0);
        se[i].ReducedDictionary.fill(0);
        se[i].ReducedDictionaryTranspose.fill(0);
        se[i].neuronCount = 0;
        se[i].threshold = std::log(pii/(1-pii));
        std::cout<< "Threshold is " << se[i].threshold;
        se[i].spikeDetected = false;
        se[i].samplesBeingCollected = false;
        se[i].likelihoodThreshold = 0;
        se[i].idx = 0;
        se[i].Hadj = 0;
        se[i].setLambda = true;
        se[i].totalSpikesFound = 0;
        se[i].spikeBuffer = new uint8_t[MAX_SORTED_SPIKE_BUFFER_LEN];
        se[i].sigma = Eigen::MatrixXf::Zero(P,P);

        for (int j = 0; j < MAX_CHANNELS; j++)
        {
            se[i].sec[j].cLL = Eigen::VectorXf::Zero(Cmax);
            se[i].sec[j].muu = Eigen::MatrixXf::Constant(K, Cmax, 0);
            se[i].sec[j].muu0 = Eigen::MatrixXf::Constant(K, Cmax, 0);
            se[i].sec[j].xwind.fill(0);

            for (int k = 0; k < Cmax; k++)
            {
                se[i].sec[j].phi.add(phi0);
                se[i].sec[j].lamclus.add(phi0nu0forlamclus);
                se[i].sec[j].lamclusinv.add(phi0nu0forlamclusinv);
                se[i].sec[j].R.add(phi0nu0forR);
                se[i].sec[j].Rinv.add(phi0nu0forRinv);
            }

        }
    }

    setParameter(1,0);
    suppresslikelihood = 1e5;
    masterSampleIndex = 0;
    allParametersEstimated = false;
    downsamplingFactor = 4e5/1e5;
}

SpikeSorter::~SpikeSorter()
{

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
    int electrodeIndex = int(newValue);
    if (parameterIndex == 1 )
    {
        for(int k = 0; k < MAX_ELECTRODES; k++)
        {
            for (int i = 0; i < se[k].ltheta.size(); i++)
            {
                if (!checkIfLogIsInf(se[k].ngamma[i],(alpha+se[k].totalSpikesFound)))
                {
                    se[k].ltheta(i) = (std::log(se[k].ngamma[i]/(alpha+se[k].totalSpikesFound)));

                }
                else
                {
                    se[k].ltheta(i) = -1*PI;
                }

            }
            se[k].ltheta(se[k].neuronCount) = std::log(alpha/(alpha+se[k].totalSpikesFound));
        }
    }

    if (parameterIndex == 2 )
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
            std::cout << "Could not find a the Parameter Estimator." << std::endl;
        }

        se[electrodeIndex].ReducedDictionary = node->dictionaryThread.dictionary[electrodeIndex];
        se[electrodeIndex].ReducedDictionaryTranspose = se[electrodeIndex].ReducedDictionary.transpose();
        std::ofstream myfile;
        myfile.open("dictionary.csv", std::fstream::app);
        myfile << se[electrodeIndex].ReducedDictionary;
        //cout << "ElectrodeIndex is " << node->dictionaryThread.dictionary[electrodeIndex] << endl;

        se[electrodeIndex].RTLR = se[electrodeIndex].ReducedDictionaryTranspose*se[electrodeIndex].lambda*se[electrodeIndex].ReducedDictionary;
        se[electrodeIndex].RTL = se[electrodeIndex].ReducedDictionaryTranspose*se[electrodeIndex].lambda;
        se[electrodeIndex].LR = se[electrodeIndex].lambda*se[electrodeIndex].ReducedDictionary;
        se[electrodeIndex].detSigma = se[electrodeIndex].sigma.determinant();

        for (int j = 0; j < MAX_CHANNELS; j++)
        {
            for (int i = 0; i < Cmax; i++)
            {
                se[electrodeIndex].C.noalias() = se[electrodeIndex].sec[j].Rinv[i] + se[electrodeIndex].sec[j].lamclusinv[i];
                se[electrodeIndex].bracket.noalias() = se[electrodeIndex].C.inverse() + se[electrodeIndex].RTLR;

                se[electrodeIndex].sec[j].Q.add(se[electrodeIndex].lambda - se[electrodeIndex].LR*se[electrodeIndex].bracket.inverse()*se[electrodeIndex].RTL);
                se[electrodeIndex].sec[j].detQ.add(se[electrodeIndex].bracket.determinant()*se[electrodeIndex].C.determinant()*se[electrodeIndex].detSigma);
                se[electrodeIndex].sec[j].cLL(i) = -(float(P)*LOG2PIBY2) - 0.5*(std::log(std::abs(se[electrodeIndex].sec[j].detQ[i])));
            }
        }
    }

    if (parameterIndex == 3)
    {
        for (int k = 0; k < MAX_CHANNELS; k++)
        {
            se[electrodeIndex].sec[k].circbuffer.setsize(P);
        }
    }
}

float SpikeSorter::getNextSample(int& chan)
{
        if (sampleIndex < dataBuffer.getNumSamples())
            return *dataBuffer.getSampleData(chan, sampleIndex);
        else
            return 0;
}

void SpikeSorter::PackageCurrentSortedSpikeIntoBuffer(MidiBuffer& eventBuffer1, int electrodeIndex)
{
    int numBytes = packSortedSpike(se[electrodeIndex].currentSpike, se[electrodeIndex].spikeBuffer, MAX_SPIKE_BUFFER_LEN);

    if (numBytes > 0)
        eventBuffer1.addEvent(se[electrodeIndex].spikeBuffer, numBytes, int(se[electrodeIndex].currentSpike.timestamp));
}


void SpikeSorter::collectSamplesForSpikeObject(int electrodeIndex, int channelCount)
{
    se[electrodeIndex].samplesBeingCollected = true;
    se[electrodeIndex].currentIndex = 0;
    se[electrodeIndex].currentSpike.eventType = SORTEDSPIKE;

    se[electrodeIndex].currentSpike.nChannels = 4;//node->electrodes[electrodeIndex]->numChannels;
    se[electrodeIndex].currentSpike.nSamples = P;
    se[electrodeIndex].currentSpike.source = electrodeIndex;
     // HARDCODING ONE CHANNEL FOR NOW

    for (int i = 0; i < MAX_NUMBER_OF_SPIKE_CHANNELS; i++)
    {
    int chan = *(node->electrodes[electrodeIndex]->channels+ i);
    se[electrodeIndex].currentSpike.threshold[i] = se[electrodeIndex].threshold;
    se[electrodeIndex].currentSpike.gain[i] = (int)(1.0f / channels[chan]->bitVolts)*1000;
    }

    se[electrodeIndex].currentSpike.timestamp = masterSampleIndex + sampleIndex;
    se[electrodeIndex].lthr = Eigen::VectorXf::Zero(range + 1);
    se[electrodeIndex].lon = Eigen::MatrixXf::Zero(range + 1, se[electrodeIndex].neuronCount + 1);

    for(int k = 0; k < channelCount; k++)
    {
    se[electrodeIndex].sec[k].xwindLonger = Eigen::VectorXf::Zero(se[electrodeIndex].sec[k].xwind.size() + range);
    se[electrodeIndex].sec[k].xwindLonger.head(se[electrodeIndex].sec[k].xwind.size()) = se[electrodeIndex].sec[k].xwind;
    }

    se[electrodeIndex].lthr(se[electrodeIndex].currentIndex) = se[electrodeIndex].likelihoodThreshold;
    for (int i = 0; i <= se[electrodeIndex].neuronCount; i++)
    {
        se[electrodeIndex].lon(se[electrodeIndex].currentIndex, i) = se[electrodeIndex].likelihoodPerNeuron(i);
    }
    se[electrodeIndex].currentIndex++;
}

void SpikeSorter::addNewSampleAndLikelihoodsToCurrentSpikeObject(Eigen::Matrix<float, 4, 1> sample, MidiBuffer& eventBuffer, int electrodeIndex, int channelCount)
{
    if (se[electrodeIndex].currentIndex <= range)
    {
        se[electrodeIndex].lthr(se[electrodeIndex].currentIndex) = se[electrodeIndex].likelihoodThreshold;
        for (int i = 0; i <= se[electrodeIndex].neuronCount; i++)
        {
            se[electrodeIndex].lon(se[electrodeIndex].currentIndex, i) = se[electrodeIndex].likelihoodPerNeuron(i);
        }
        for (int k = 0; k < MAX_CHANNELS; k++)
        {
        se[electrodeIndex].sec[k].xwindLonger(se[electrodeIndex].sec[k].xwind.size() + se[electrodeIndex].currentIndex - 1) = sample(k);
        }
        se[electrodeIndex].currentIndex++;
    }
    else
    {
        se[electrodeIndex].spikeDetected = false;
        se[electrodeIndex].samplesBeingCollected = false;
        double minthr = se[electrodeIndex].lthr.segment(0,range).minCoeff(&(se[electrodeIndex].idx));
        //cout << se[electrodeIndex].idx << "// ";

                //std::ofstream myfile;
        //myfile.open("xwindlonger.csv", fstream::app);
        //myfile << se[electrodeIndex].sec[0].xwindLonger.segment(se[electrodeIndex].idx, P).transpose() << endl;

        int Cnew;
        double maxlon = se[electrodeIndex].lon.row(se[electrodeIndex].idx).maxCoeff(&Cnew);
        if (Cnew == se[electrodeIndex].neuronCount)
        {
            cout<< "reached here so apparently new neuron" << endl;
            se[electrodeIndex].neuronCount++;
        }
        se[electrodeIndex].currentSpike.neuronID = Cnew;
        se[electrodeIndex].currentSpike.nDictionary = K;
        se[electrodeIndex].currentIndex = 0;

        for (int k = 0; k < channelCount; k++)
        {
            for (int i = 0; i < P; i ++)
            {
                se[electrodeIndex].currentSpike.data[P*k+i] = uint16(se[electrodeIndex].sec[k].xwindLonger(se[electrodeIndex].idx + i)/ channels[k]->bitVolts + 32768);
            }
        }

        // --- Update just enough to get yhat and pack spike
        se[electrodeIndex].deltaT = masterSampleIndex + sampleIndex -P - range + se[electrodeIndex].idx - se[electrodeIndex].tlastspike(se[electrodeIndex].currentSpike.neuronID);
        se[electrodeIndex].tlastspike(se[electrodeIndex].currentSpike.neuronID) = masterSampleIndex + sampleIndex - P - range + se[electrodeIndex].idx;
        // ----

        for (int k = 0; k < MAX_CHANNELS; k++)
        updateAllSortingParameters(electrodeIndex, k);

        PackageCurrentSortedSpikeIntoBuffer(eventBuffer, electrodeIndex);

        int numBytes = packSortedSpike(se[electrodeIndex].currentSpike, se[electrodeIndex].spikeBuffer, MAX_SPIKE_BUFFER_LEN);
        if (numBytes > 0)
        eventBuffer.addEvent(se[electrodeIndex].spikeBuffer, numBytes, int(se[electrodeIndex].currentSpike.timestamp));

        se[electrodeIndex].totalSpikesFound += 1;

    }
}


void SpikeSorter::updateAllSortingParameters(int electrodeIndex, int chan)
{
    //Logger *log1 = Logger::getCurrentLogger();
    //double startTimet = Time::getMillisecondCounterHiRes();
    int neuronID = se[electrodeIndex].currentSpike.neuronID;
    se[electrodeIndex].Qmat.noalias() = se[electrodeIndex].RTLR + se[electrodeIndex].sec[chan].lamclus[neuronID];
    se[electrodeIndex].Qmatinv = se[electrodeIndex].Qmat.inverse();

    se[electrodeIndex].yhat.noalias() = se[electrodeIndex].Qmatinv*(se[electrodeIndex].RTL*(se[electrodeIndex].sec[chan].xwindLonger.segment(se[electrodeIndex].idx, P))
                                        + se[electrodeIndex].sec[chan].lamclus[neuronID]*se[electrodeIndex].sec[chan].muu.col(neuronID));
    std::ofstream myfile;
    myfile.open("pcomponents.csv", fstream::app);
    if (chan == 0)
    {
        for (int i = 0; i < se[electrodeIndex].yhat.size(); i++) // right now only the last will get saved. use average?
        {
            se[electrodeIndex].currentSpike.principalComponent[i] = uint16(se[electrodeIndex].yhat(i)/ channels[0]->bitVolts + 32768);
            myfile << se[electrodeIndex].yhat(i);
            if ( i != se[electrodeIndex].yhat.size() -1 )
            {
                myfile << ",";
            }
        }
        myfile << endl;
    }
    se[electrodeIndex].ngamma(neuronID) = se[electrodeIndex].ngamma(neuronID) + 1;
    double ebet = std::exp(-1*double(beta)*double(se[electrodeIndex].deltaT));
    se[electrodeIndex].mhat.noalias() = (se[electrodeIndex].sec[chan].muu0.col(neuronID).array()*(1-ebet) + se[electrodeIndex].sec[chan].muu.col(neuronID).array()*ebet).matrix();
    se[electrodeIndex].sec[0].muu0.col(neuronID) = ((se[electrodeIndex].kappa(neuronID)*se[electrodeIndex].sec[chan].muu0.col(neuronID).array() + se[electrodeIndex].yhat.array())/(se[electrodeIndex].kappa(neuronID)+1)).matrix();
    se[electrodeIndex].Qhat.noalias() = ((1/tau)*Eigen::MatrixXf::Identity(K,K).array()*(1-ebet*ebet) + se[electrodeIndex].sec[chan].R[neuronID].inverse().array()*(ebet*ebet)).matrix();
    se[electrodeIndex].sec[chan].R.setUnchecked(neuronID, se[electrodeIndex].Qhat.inverse() + se[electrodeIndex].sec[chan].lamclus[neuronID]);
    se[electrodeIndex].sec[chan].Rinv.setUnchecked(neuronID, se[electrodeIndex].sec[chan].R[neuronID].inverse());
    se[electrodeIndex].sec[chan].muu.col(neuronID) = se[electrodeIndex].sec[chan].Rinv[neuronID]*(se[electrodeIndex].Qhat.inverse()*(se[electrodeIndex].mhat) + se[electrodeIndex].sec[chan].lamclus[neuronID]*se[electrodeIndex].yhat);
    Eigen::MatrixXf temp = se[electrodeIndex].RTL*se[electrodeIndex].sec[chan].xwindLonger.segment(se[electrodeIndex].idx,P) + se[electrodeIndex].sec[chan].lamclus[neuronID]*se[electrodeIndex].sec[chan].muu.col(neuronID);
    se[electrodeIndex].yhat = se[electrodeIndex].Qmatinv*(temp);






    double constant = (se[electrodeIndex].kappa(neuronID)/(se[electrodeIndex].kappa(neuronID) + 1));
    Eigen::Matrix3f trans = (se[electrodeIndex].yhat - se[electrodeIndex].sec[chan].muu.col(neuronID))*((se[electrodeIndex].yhat - se[electrodeIndex].sec[chan].muu.col(neuronID)).transpose());
    Eigen::Matrix3f newphi = (constant*trans.array()).matrix();
    Eigen::Matrix3f temp1 = newphi + se[electrodeIndex].sec[chan].phi[neuronID] + se[electrodeIndex].Qmatinv;
    se[electrodeIndex].sec[chan].phi.setUnchecked(neuronID, temp1);
    se[electrodeIndex].kappa(neuronID) = se[electrodeIndex].kappa(neuronID) + 1;
    se[electrodeIndex].nu(neuronID) = se[electrodeIndex].nu(neuronID) + 1;
    se[electrodeIndex].sec[chan].lamclus.setUnchecked(neuronID, (se[electrodeIndex].sec[chan].phi[neuronID].inverse().array()*se[electrodeIndex].nu(neuronID)).matrix());
    se[electrodeIndex].sec[chan].lamclusinv.setUnchecked(neuronID, se[electrodeIndex].sec[chan].lamclus[neuronID].inverse());
    se[electrodeIndex].thingsHaveChanged = 1;
    se[electrodeIndex].C.noalias() = se[electrodeIndex].sec[chan].Rinv[neuronID] + se[electrodeIndex].sec[chan].lamclusinv[neuronID];
    se[electrodeIndex].bracket.noalias() = se[electrodeIndex].C.inverse() + se[electrodeIndex].RTLR;
    se[electrodeIndex].sec[chan].Q.setUnchecked(neuronID, se[electrodeIndex].lambda - se[electrodeIndex].LR*se[electrodeIndex].bracket.inverse()*se[electrodeIndex].RTL);
    se[electrodeIndex].sec[chan].detQ.setUnchecked(neuronID, se[electrodeIndex].bracket.determinant()*se[electrodeIndex].C.determinant()*se[electrodeIndex].detSigma);
    se[electrodeIndex].sec[chan].cLL(neuronID) = -(float(P)*LOG2PIBY2) - 0.5*std::log(std::abs(se[electrodeIndex].sec[chan].detQ[neuronID]));
    //double stopTimet = Time::getMillisecondCounterHiRes();
    //log1->writeToLog("UPDATE DURATION is " + String(stopTimet - startTimet));
}


void SpikeSorter::process(AudioSampleBuffer& buffer,
                          MidiBuffer& events,
                          int& nSamples)
{
    dataBuffer = buffer;

    //checkForEvents(events);

    if(allParametersEstimated)
    {
        // spike sort waveform using Vogelstein code

        Logger *log1 = Logger::getCurrentLogger();
        startTime = Time::getMillisecondCounterHiRes();

        for (int i = 0; i < (*node).electrodes.size(); i++)
        {
            sampleIndex = 0;

            int channelCount = (*node).electrodes[i]->numChannels;

            //for (int chan = 0; chan < 1; chan++) // right now, ONLY ONE CHANNEL
            {
                //if ( *(node->electrodes[i]->isActive +chan) )
                {
                    while (sampleIndex < nSamples)
                    {
                        for (int k = 0 ; k < channelCount; k++)
                        {

                            se[i].sample(k) = *dataBuffer.getSampleData(k, sampleIndex);
                            se[i].sec[k].circbuffer.enqueue(se[i].sample(k));
                        }

                        pii = (apii + se[i].totalSpikesFound)/(bpii + masterSampleIndex);
                        se[i].threshold = std::log(pii/(1-pii));

                        if(se[i].sec[0].circbuffer.isBufferPlush(P))
                        {

                            for (int k = 0 ; k < channelCount; k++)
                            {
                                se[i].sec[k].circbuffer.show(P, se[i].sec[k].xwind);
                            }

                            if(se[i].thingsHaveChanged)
                            {

                                for (int i = 0; i < se[i].neuronCount; i++)
                                {
                                    {
                                        se[i].ltheta(i) = std::log(se[i].ngamma(i)) - std::log(alpha+se[i].totalSpikesFound);

                                    }
                                }

                                se[i].ltheta(se[i].neuronCount) = std::log(alpha) - std::log(alpha+se[i].totalSpikesFound);
                                se[i].thingsHaveChanged = 0;
                            }

                            se[i].likelihoodNoSpike = se[i].logPlusDetTermForNoiseLL -0.5*(se[i].sec[0].xwind.transpose()*(se[i].lambda*(se[i].sec[0].xwind)))(0);


                            for (int k = 1 ; k < channelCount; k++)
                            {
                            se[i].likelihoodNoSpike += se[i].logPlusDetTermForNoiseLL -0.5*(se[i].sec[k].xwind.transpose()*(se[i].lambda*(se[i].sec[k].xwind)))(0);
                            }

                            se[i].likelihoodNoSpike /= 4;

                            for (int j = 0; j < se[i].neuronCount; j++)
                            {
                                se[i].sec[0].xwindloop = se[i].sec[0].xwind;
                                se[i].sec[0].xwindloop.noalias() -=  se[i].ReducedDictionary*(se[i].sec[0].muu.col(j)); // <-- think about this later

                                if( ( (masterSampleIndex + float(sampleIndex)) - se[i].tlastspike(j)) < 10000*50/10000 )// THIS NEEDS TO BE INVESTIGATED
                                {
                                    se[i].likelihoodPerNeuron(j) = se[i].sec[0].cLL(j) - suppresslikelihood;
                                    se[i].likelihoodPerNeuron(j) -= 0.5*((se[i].sec[0].xwindloop.transpose())*(se[i].sec[0].Q[j]*(se[i].sec[0].xwindloop)))(0);
                                }
                                else
                                {
                                    se[i].likelihoodPerNeuron(j) = se[i].sec[0].cLL(j);
                                    se[i].likelihoodPerNeuron(j) -= 0.5*((se[i].sec[0].xwindloop.transpose())*(se[i].sec[0].Q[j]*(se[i].sec[0].xwindloop)))(0);
                                }
                            }

                            for (int j = 0; j < se[i].neuronCount; j++)
                            {
                                for (int k = 1 ; k < channelCount; k++)
                                {
                                    se[i].sec[k].xwindloop = se[i].sec[k].xwind;
                                    se[i].sec[k].xwindloop.noalias() -=  se[i].ReducedDictionary*(se[i].sec[k].muu.col(j)); // <-- think about this later

                                    if( ( (masterSampleIndex + float(sampleIndex)) - se[i].tlastspike(j)) < 10000*50/10000 )// THIS NEEDS TO BE INVESTIGATED
                                    {
                                        se[i].likelihoodPerNeuron(j) += se[i].sec[k].cLL(j) - suppresslikelihood;
                                        se[i].likelihoodPerNeuron(j) -= 0.5*((se[i].sec[k].xwindloop.transpose())*(se[i].sec[k].Q[j]*(se[i].sec[k].xwindloop)))(0);
                                    }
                                    else
                                    {
                                        se[i].likelihoodPerNeuron(j) += se[i].sec[k].cLL(j);
                                        se[i].likelihoodPerNeuron(j) -= 0.5*((se[i].sec[k].xwindloop.transpose())*(se[i].sec[k].Q[j]*(se[i].sec[k].xwindloop)))(0);
                                    }
                                }
                                se[i].likelihoodPerNeuron(j) /= 4;
                            }

                            se[i].likelihoodPerNeuron(se[i].neuronCount) = se[i].sec[0].cLL(se[i].neuronCount);
                            se[i].likelihoodPerNeuron(se[i].neuronCount) -= 0.5*(se[i].sec[0].xwind.transpose()*(se[i].sec[0].Q[se[i].neuronCount]*(se[i].sec[0].xwind)))(0);

                            for (int k = 1 ; k < channelCount; k++)
                            {
                                se[i].likelihoodPerNeuron(se[i].neuronCount) += se[i].sec[k].cLL(se[i].neuronCount);
                                se[i].likelihoodPerNeuron(se[i].neuronCount) -= 0.5*(se[i].sec[k].xwind.transpose()*(se[i].sec[k].Q[se[i].neuronCount]*(se[i].sec[k].xwind)))(0);
                            }

                            se[i].likelihoodPerNeuron(se[i].neuronCount) /= 4;

                            se[i].likelihoodPerNeuron.head(1+se[i].neuronCount) += se[i].ltheta.head(1+se[i].neuronCount);
                            se[i].maximumLikelihoodPerNeuron = se[i].likelihoodPerNeuron.head(1+se[i].neuronCount).maxCoeff();
                            se[i].likelihoodPerNeuron.head(1+se[i].neuronCount).array() -= se[i].maximumLikelihoodPerNeuron;
                            se[i].Hadjsum = 0;

                            for (int k = 0; k <= se[i].neuronCount; k++ )
                            {
                                se[i].Hadjsum += std::exp(se[i].likelihoodPerNeuron(k));
                            }

                            se[i].Hadj = std::log(se[i].Hadjsum);

                            se[i].likelihoodThreshold = se[i].likelihoodNoSpike - se[i].maximumLikelihoodPerNeuron - se[i].Hadj;

                            if (se[i].likelihoodThreshold < se[i].threshold)
                                se[i].spikeDetected = true;

                            if(se[i].spikeDetected)
                            {
                                if (se[i].samplesBeingCollected)
                                {
                                    addNewSampleAndLikelihoodsToCurrentSpikeObject(se[i].sample, events, i, channelCount);
                                }
                                else
                                {
                                    collectSamplesForSpikeObject(i, channelCount);
                                }
                            }

                            for (int k = 0 ; k < channelCount; k++)
                            {
                                se[i].sec[k].circbuffer.dequeue();
                            }
                            // After everything is done, remove most past sample, ie, xwind(0)
                        }
                    sampleIndex += 4;
                    }
                }
            }
        }
        masterSampleIndex += sampleIndex;
        stopTime = Time::getMillisecondCounterHiRes();

        timetaken += stopTime - startTime;
        n++;
        if (timeflag == false && n > 500)
        {
            log1->writeToLog("Average Time is " + String(timetaken/n) + " ms.");
            timeflag = true;
        }
    }
}


      //likelihoodPerNeuron(neuronCount) -= 0.5*((xwind.transpose())*(QLLT[neuronCount].solve(xwind)))(0);

                            /*if (flag == false)
                            {
                                //cout << " The quotient is " << Eigen::Matrix<float,1,40>::Ones(1,40)*(Q[neuronCount]*(Eigen::Matrix<float,40,1>::Ones(40,1))) << endl;
                                //cout << "The solver is " << Eigen::Matrix<float,1,40>::Ones(1,40)*(QLLT[neuronCount].solve(Eigen::Matrix<float,40,1>::Ones(40,1))) << endl;
                                //cout << xwind.transpose();
                                likelihoodPerNeuron(neuronCount) = cLL(neuronCount);
                                likelihoodPerNeuron(neuronCount) -= 0.5*(xwind.transpose()*(Q[neuronCount]*(xwind)))(0);
                                likelihoodPerNeuron.head(1+neuronCount) += ltheta.head(1+neuronCount);
                                maximumLikelihoodPerNeuron = likelihoodPerNeuron.head(1+neuronCount).maxCoeff();
                                likelihoodPerNeuron.head(1+neuronCount).array() -= maximumLikelihoodPerNeuron;
                                Hadjsum = 0;
                                for (int i = 0; i <= neuronCount; i++ )
                                {
                                    Hadjsum += std::exp(likelihoodPerNeuron(i));
                                }

                                Hadj = std::log(Hadjsum);
                                likelihoodThreshold = likelihoodNoSpike - maximumLikelihoodPerNeuron - Hadj;
                                cout << "The inverse is "<< likelihoodThreshold << " :: "<< likelihoodNoSpike << " :: "<< maximumLikelihoodPerNeuron << " :: ";
                                likelihoodPerNeuron(neuronCount) = cLL(neuronCount);
                                likelihoodPerNeuron(neuronCount) -= 0.5*((xwind.transpose())*(QLLT[neuronCount].solve(xwind)))(0);
                                likelihoodPerNeuron.head(1+neuronCount) += ltheta.head(1+neuronCount);
                                maximumLikelihoodPerNeuron = likelihoodPerNeuron.head(1+neuronCount).maxCoeff();
                                likelihoodPerNeuron.head(1+neuronCount).array() -= maximumLikelihoodPerNeuron;
                                Hadjsum = 0;
                                for (int i = 0; i <= neuronCount; i++ )
                                {
                                    Hadjsum += std::exp(likelihoodPerNeuron(i));
                                }

                                Hadj = std::log(Hadjsum);

                                likelihoodThreshold = likelihoodNoSpike - maximumLikelihoodPerNeuron - Hadj;
                                likelihoodPerNeuron(neuronCount) -= 0.5*((xwind.transpose())*(QLLT[neuronCount].solve(xwind)))(0);

                                cout << "The solver is"<< likelihoodThreshold << " :: "<< likelihoodNoSpike << " :: "<< maximumLikelihoodPerNeuron << " :: ";
                                //flag = true;
                            }

                            likelihoodPerNeuron(neuronCount) = cLL(neuronCount);
                            likelihoodPerNeuron(neuronCount) -= 0.5*((xwind.transpose())*(QLLT[neuronCount].solve(xwind)))(0);*/

                            /*for (int i = 0; i <= neuronCount; i++ )
                            {
                                likelihoodPerNeuron(i) =  likelihoodPerNeuron(i) + ltheta(i);
                            }*/
                            /*for (int i = 0; i <= neuronCount; i++ )
                            {
                                    if(maximumLikelihoodPerNeuron < likelihoodPerNeuron(i))
                                    maximumLikelihoodPerNeuron = likelihoodPerNeuron(i);
                            }*/
                                                        /*for (int i = 0; i <= neuronCount; i++ )
                            {
                                likelihoodPerNeuron(i) = likelihoodPerNeuron(i) - maximumLikelihoodPerNeuron;
                            }*/
