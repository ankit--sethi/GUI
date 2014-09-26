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

void CircularQueue::show(int index, Eigen::VectorXf &returnedX) // caution: works right only if size of circular buffer = P
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

UpdateThread::UpdateThread(SpikeSorter *spikeSorter) : Thread("UpdateThread")
{
    ss = spikeSorter;
}

void UpdateThread::run()
{
    /*int neuronID = ss->currentSpike.neuronID;
    ss->ngamma(neuronID) = ss->ngamma(neuronID) + 1;
    ss->deltaT = ss->masterSampleIndex + ss->sampleIndex - ss->P - ss->range + ss->idx - ss->tlastspike(neuronID);
    ss->tlastspike(neuronID) = ss->masterSampleIndex + ss->sampleIndex - ss->P - ss->range + ss->idx;
    double ebet = std::exp(-1*double(ss->beta)*double(ss->deltaT));
    ss->mhat = (ss->muu0.col(neuronID).array()*(1-ebet) + ss->muu.col(neuronID).array()*ebet).matrix();
    ss->muu0.col(neuronID) = ((ss->kappa(neuronID)*ss->muu0.col(neuronID).array() + ss->yhat.array())/(ss->kappa(neuronID)+1)).matrix();
    ss->Qhat = ((1/ss->tau)*Eigen::MatrixXf::Identity(K,K).array()*(1-ebet*ebet) + ss->R[neuronID].inverse().array()*(ebet*ebet)).matrix();
    ss->R.setUnchecked(neuronID, ss->Qhat.inverse() + ss->lamclus[neuronID]);
    ss->Rinv.setUnchecked(neuronID, ss->R[neuronID].inverse());
    ss->muu.col(neuronID) = ss->R[neuronID].inverse()*(ss->Qhat.inverse()*(ss->mhat) + ss->lamclus[neuronID]*ss->yhat);
    Eigen::MatrixXf temp = ss->ReducedDictionaryTranspose*ss->lambda*ss->xwindLonger.segment(ss->idx,ss->P) + ss->lamclus[neuronID]*ss->muu.col(neuronID);
    ss->yhat = ss->Qmat.inverse()*(temp);
    double constant = (ss->kappa(neuronID)/(ss->kappa(neuronID) + 1));
    Eigen::MatrixXf trans = (ss->yhat - ss->muu.col(neuronID))*((ss->yhat - ss->muu.col(neuronID)).transpose());
    Eigen::MatrixXf newphi = (constant*trans.array()).matrix();
    Eigen::MatrixXf temp1 = newphi + ss->phi[neuronID] + ss->Qmat.inverse();
    ss->phi.setUnchecked(neuronID, temp1);
    ss->kappa(neuronID) = ss->kappa(neuronID) + 1;
    ss->nu(neuronID) = ss->nu(neuronID) + 1;
    ss->lamclus.setUnchecked(neuronID, (ss->phi[neuronID].inverse().array()*ss->nu(neuronID)).matrix());
    ss->lamclusinv.setUnchecked(neuronID, ss->lamclus[neuronID].inverse());
    thingsHaveChanged = 1;
    ss->QQR[neuronID].compute(ss->sigma + ss->ReducedDictionary*(ss->Rinv[neuronID] + ss->lamclusinv[neuronID])*ss->ReducedDictionaryTranspose);
    ss->QLogAbsDeterminant(neuronID) = ss->QQR[neuronID].logAbsDeterminant();*/
    ss->updateAllSortingParameters();
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
    phi0 = (Eigen::MatrixXf::Identity(K,K).array()*0.1).matrix();
    //check phi0!!!
    number = 0;
    apii = 1;
    bpii = 1e7;
    beta = 1/(30*(10000));
    //std::cout << " Sample RATE is " << getSampleRate() << "//";
    tau = 10;

    P = 40;

    checkIfAllParametersEstimated = false;
    paramsCopied = false;


    Cmax = 50;  //maximum possible number of neurons present
    curndx = 0; //used to index current location in buffer
    lookahead = 500; //functionally, it is the size of the buffer
    range = 15; // defines basically the width of a spike
    // setting sigma
    sigma = Eigen::MatrixXf::Zero(P,P);

    pii=apii/bpii;

    thingsHaveChanged = 1;

    Eigen::MatrixXf phi0nu0forlamclus = ((phi0.inverse().array()*nu0)).matrix();
    Eigen::MatrixXf phi0nu0forlamclusinv = phi0nu0forlamclus.inverse();
    Eigen::MatrixXf phi0nu0forR = ((phi0.inverse().array()*nu0)*0.2).matrix();
    Eigen::MatrixXf phi0nu0forRinv = phi0nu0forR.inverse();
    nu = Eigen::VectorXf::Constant(Cmax, nu0);
    kappa = Eigen::VectorXf::Constant(Cmax, kappa0);
    ltheta = Eigen::VectorXf::Zero(Cmax);
    ngamma = Eigen::VectorXf::Zero(Cmax);
    tlastspike = Eigen::VectorXf::Zero(Cmax);
    likelihoodPerNeuron = Eigen::VectorXf::Zero(Cmax);
    cLL = Eigen::VectorXf::Zero(Cmax);
    QLogAbsDeterminant = Eigen::VectorXf::Zero(Cmax);

    for (int i = 0; i < Cmax; i++)
    {
        phi.add(phi0);
        lamclus.add(phi0nu0forlamclus);
        lamclusinv.add(phi0nu0forlamclusinv);
        R.add(phi0nu0forR);
        Rinv.add(phi0nu0forRinv);
    }
    muu = Eigen::MatrixXf::Constant(K, Cmax, 0);
    muu0 = Eigen::MatrixXf::Constant(K, Cmax, 0);
    //Q = Eigen::MatrixXf::Random(P, P);
    Qinv = Eigen::MatrixXf::Constant(P, P, 0);
    Qhat = Eigen::MatrixXf::Constant(K, K, 0);
    Qmat = Eigen::MatrixXf::Constant(K, K, 0);
    mhat = Eigen::VectorXf::Constant(K, 1, 0);
    yhat = Eigen::VectorXf::Constant(K, 1, 0);
    ReducedDictionary = Eigen::MatrixXf::Constant(P, K, 0);
    ReducedDictionaryTranspose = Eigen::MatrixXf::Constant(K, P, 0);


    neuronCount = 0; // This is 'C' in the MATLAB code.
    //nz = 0; // Not sure what this does. Possible the same thing.
    setParameter(1,0);
    xwind = Eigen::VectorXf::Zero(P);
    xRDmu = Eigen::MatrixXf::Constant(Cmax, P, 0);

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
    allParametersEstimated = false;
}

SpikeSorter::~SpikeSorter()
{
    for (int i = 0; i < circbuffer.size(); i++)
    {
       // delete circbuffer[i];
    }
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
                ltheta(i) = (std::log(ngamma[i]/(alpha+totalSpikesFound)));

            }
            else
            {
                ltheta(i) = -1*PI;
            }

        }
        ltheta(neuronCount) = std::log(alpha/(alpha+totalSpikesFound));

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

        ReducedDictionary = node->dictionaryThread.dictionary;
        ReducedDictionaryTranspose = ReducedDictionary.transpose();

        Eigen::HouseholderQR<Eigen::MatrixXf> t;
        Eigen::LLT<Eigen::MatrixXf> tllt;
        for (int i = 0; i < Cmax; i++)
        {
            t.compute(sigma + ReducedDictionary*(Rinv[i] + lamclusinv[i])*ReducedDictionaryTranspose);
            tllt.compute(sigma + ReducedDictionary*(Rinv[i] + lamclusinv[i])*ReducedDictionaryTranspose);
            QQR.add(t);
            QLLT.add(tllt);
            QLogAbsDeterminant(i) = QQR[i].logAbsDeterminant();
            cLL(i) = -(float(P)*LOG2PIBY2) - 0.5*QLogAbsDeterminant(i);
            xRDmu.row(i) = (ReducedDictionary*muu.col(i)).transpose();
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

void SpikeSorter::PackageCurrentSortedSpikeIntoBuffer(MidiBuffer& eventBuffer1)
{
    int numBytes = packSortedSpike(currentSpike, spikeBuffer, MAX_SPIKE_BUFFER_LEN);

    if (numBytes > 0)
        eventBuffer1.addEvent(spikeBuffer, numBytes, int(currentSpike.timestamp));
}


void SpikeSorter::collectSamplesForSpikeObject(int electrodeIndex, float trigSample)
{
    samplesBeingCollected = true;
    currentIndex = 0;
    currentSpike.eventType = SORTEDSPIKE;

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
    lthr = Eigen::VectorXf::Zero(range + 1);
    lon = Eigen::MatrixXf::Zero(range + 1, neuronCount + 1);
    xwindLonger = Eigen::VectorXf::Zero(xwind.size() + range);
    xwindLonger.head(xwind.size()) = xwind;

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
        for (int i = 0; i <= neuronCount; i++)
        {
            lon(currentIndex, i) = likelihoodPerNeuron(i);
        }
        xwindLonger(xwind.size() + currentIndex - 1) = sample;
        currentIndex++;
    }
    else
    {
        spikeDetected = false;
        samplesBeingCollected = false;
        double minthr = lthr.minCoeff(&idx);

        int Cnew;
        double maxlon = lon.row(idx).maxCoeff(&Cnew);
        if (Cnew == neuronCount)
        {
            cout<< "reached here so apparently new neuron" << endl;
            neuronCount++;
        }
        currentSpike.neuronID = Cnew;
        currentSpike.nDictionary = K;
        currentIndex = 0;
        for (int i = 0; i < P; i ++)
        {
        currentSpike.data[i] = uint16(xwindLonger(idx + i)/ channels[chan]->bitVolts + 32768);

        }

        // --- Update just enough to get yhat and pack spike
        deltaT = masterSampleIndex + sampleIndex -P - range + idx - tlastspike(currentSpike.neuronID);
        tlastspike(currentSpike.neuronID) = masterSampleIndex + sampleIndex - P - range + idx;
        // ----
        updateAllSortingParameters();

        PackageCurrentSortedSpikeIntoBuffer(eventBuffer);

        int numBytes = packSortedSpike(currentSpike, spikeBuffer, MAX_SPIKE_BUFFER_LEN);
        if (numBytes > 0)
        eventBuffer.addEvent(spikeBuffer, numBytes, int(currentSpike.timestamp));

        totalSpikesFound += 1;

    }
}


void SpikeSorter::updateAllSortingParameters()
{
    Logger *log1 = Logger::getCurrentLogger();
    double startTimet = Time::getMillisecondCounterHiRes();
    int neuronID = currentSpike.neuronID;
    Qmat = ReducedDictionaryTranspose*lambda*ReducedDictionary + lamclus[neuronID];
    yhat = Qmat.inverse()*(ReducedDictionaryTranspose*lambda*(xwindLonger.segment(idx, P)) + lamclus[neuronID]*muu.col(neuronID));
    for (int i = 0; i < yhat.size(); i++)
        currentSpike.principalComponent[i] = uint16(yhat(i)/ channels[0]->bitVolts + 32768);
    ngamma(neuronID) = ngamma(neuronID) + 1;
    double ebet = std::exp(-1*double(beta)*double(deltaT));
    mhat = (muu0.col(neuronID).array()*(1-ebet) + muu.col(neuronID).array()*ebet).matrix();
    muu0.col(neuronID) = ((kappa(neuronID)*muu0.col(neuronID).array() + yhat.array())/(kappa(neuronID)+1)).matrix();
    Qhat = ((1/tau)*Eigen::MatrixXf::Identity(K,K).array()*(1-ebet*ebet) + R[neuronID].inverse().array()*(ebet*ebet)).matrix();
    R.setUnchecked(neuronID, Qhat.inverse() + lamclus[neuronID]);
    Rinv.setUnchecked(neuronID, R[neuronID].inverse());
    muu.col(neuronID) = R[neuronID].inverse()*(Qhat.inverse()*(mhat) + lamclus[neuronID]*yhat);
    Eigen::MatrixXf temp = ReducedDictionaryTranspose*lambda*xwindLonger.segment(idx,P) + lamclus[neuronID]*muu.col(neuronID);
    yhat = Qmat.inverse()*(temp);
    double constant = (kappa(neuronID)/(kappa(neuronID) + 1));
    Eigen::MatrixXf trans = (yhat - muu.col(neuronID))*((yhat - muu.col(neuronID)).transpose());
    Eigen::MatrixXf newphi = (constant*trans.array()).matrix();
    Eigen::MatrixXf temp1 = newphi + phi[neuronID] + Qmat.inverse();
    phi.setUnchecked(neuronID, temp1);
    kappa(neuronID) = kappa(neuronID) + 1;
    nu(neuronID) = nu(neuronID) + 1;
    lamclus.setUnchecked(neuronID, (phi[neuronID].inverse().array()*nu(neuronID)).matrix());
    lamclusinv.setUnchecked(neuronID, lamclus[neuronID].inverse());
    thingsHaveChanged = 1;
    trans = sigma + ReducedDictionary*(Rinv[neuronID] + lamclusinv[neuronID])*ReducedDictionaryTranspose;
    QQR[neuronID].compute(trans);
    QLLT[neuronID].compute(trans);
    QLogAbsDeterminant(neuronID) = QQR[neuronID].logAbsDeterminant();
    cLL(neuronID) = -(float(P)*LOG2PIBY2) - 0.5*QLogAbsDeterminant[neuronID];
    //xRDmu.row(neuronID) = (ReducedDictionary*muu.col(neuronID)).transpose();
    double stopTimet = Time::getMillisecondCounterHiRes();
    //log1->writeToLog("UPDATE DURATION is " + String(stopTimet - startTimet));
}


void SpikeSorter::process(AudioSampleBuffer& buffer,
                          MidiBuffer& events,
                          int& nSamples)
{
    dataBuffer = buffer;
    double stopTime0, stopTime2, stopTime1;
    double startTime0, startTime2, startTime1;

    //checkForEvents(events);

    if(allParametersEstimated)
    {
        // spike sort waveform using Vogelstein code

        Logger *log1 = Logger::getCurrentLogger();
        startTime = Time::getMillisecondCounterHiRes();

        for (int i = 0; i < (*node).electrodes.size(); i++)
        {
            sampleIndex = 0;

            for (int chan = 0; chan < 1; chan++) // right now, ONLY ONE CHANNEL
            {
                if ( *(node->electrodes[i]->isActive +chan) )
                {
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
                            startTime2 = Time::getMillisecondCounterHiRes();
                            circbuffer[i]->show(P, xwind);
                            //xwindloop = xwind;

                            if(thingsHaveChanged)
                            {

                                for (int i = 0; i < neuronCount; i++)
                                {
                                    {
                                        ltheta(i) = std::log(ngamma(i)) - std::log(alpha+totalSpikesFound);

                                    }
                                }

                                ltheta(neuronCount) = std::log(alpha) - std::log(alpha+totalSpikesFound);

                            }

                            likelihoodNoSpike = logPlusDetTermForNoiseLL -0.5*float(xwind.transpose()*(lambda*(xwind)));


                            stopTime2 = Time::getMillisecondCounterHiRes();
                            //log1->writeToLog("First Duration is " + String(stopTime2 - startTime2) + "ms for nSamples = " + String(nSamples));
                            //log1->writeToLog("Step as percentage is " + String((frac10-frac1)*100/(stopTime2 - startTime2)) + "%");
                            startTime1 = Time::getMillisecondCounterHiRes();
                            //#pragma omp parallel for
                            for (int j = 0; j < neuronCount; j++)
                            {
                                xwindloop = xwind - ReducedDictionary*muu.col(j); // <-- think about this later

                                if( ( (masterSampleIndex + float(sampleIndex)) - tlastspike(j)) < 10000*50/10000 )  // THIS NEEDS TO BE INVESTIGATED
                                    likelihoodPerNeuron(j) = cLL(j) - suppresslikelihood - 0.5*double(xwindloop.transpose()*(QLLT[j].solve(xwindloop)));
                                else
                                    likelihoodPerNeuron(j) = cLL(j) - 0.5*double(xwindloop.transpose()*(QLLT[j].solve(xwindloop)));
                            }
                            if (thingsHaveChanged)
                            {
                                //logPlusDetTermForNewNeuronLL = -(float(P)*LOG2PIBY2) - 0.5*QLogAbsDeterminant[neuronCount];
                                thingsHaveChanged = 0;
                            }
                            likelihoodPerNeuron(neuronCount) = cLL(neuronCount) - 0.5*double((xwind.transpose()*(QLLT[neuronCount].solve(xwind))));
                            stopTime1 = Time::getMillisecondCounterHiRes();
                            //log1->writeToLog("Second Duration is " + String(stopTime1 - startTime1) + "ms for nSamples = " + String(nSamples));

                            startTime0 = Time::getMillisecondCounterHiRes();
                            /*for (int i = 0; i <= neuronCount; i++ )
                            {
                                likelihoodPerNeuron(i) =  likelihoodPerNeuron(i) + ltheta(i);
                            }*/
                            likelihoodPerNeuron.head(1+neuronCount) += ltheta.head(1+neuronCount);

                            maximumLikelihoodPerNeuron = likelihoodPerNeuron.head(1+neuronCount).maxCoeff();

                            /*for (int i = 0; i <= neuronCount; i++ )
                            {
                                    if(maximumLikelihoodPerNeuron < likelihoodPerNeuron(i))
                                    maximumLikelihoodPerNeuron = likelihoodPerNeuron(i);
                            }*/
                            likelihoodPerNeuron.head(1+neuronCount).array() -= maximumLikelihoodPerNeuron;

                            /*for (int i = 0; i <= neuronCount; i++ )
                            {
                                likelihoodPerNeuron(i) = likelihoodPerNeuron(i) - maximumLikelihoodPerNeuron;
                            }*/
                            Hadjsum = 0;
                            for (int i = 0; i <= neuronCount; i++ )
                            {
                                Hadjsum += std::exp(likelihoodPerNeuron(i));
                            }

                            Hadj = std::log(Hadjsum);

                            likelihoodThreshold = likelihoodNoSpike - maximumLikelihoodPerNeuron - Hadj;
                            if (likelihoodThreshold < threshold)
                                spikeDetected = true;

                            if(spikeDetected)
                            {
                                if (samplesBeingCollected)
                                {
                                    addNewSampleAndLikelihoodsToCurrentSpikeObject(sample, events, chan);
                                }
                                else
                                {
                                    collectSamplesForSpikeObject(i, sample);
                                }
                            }

                            circbuffer[i]->dequeue(); // After everything is done, remove most past sample, ie, xwind(0)
                            double stopTime0 = Time::getMillisecondCounterHiRes();
                            //log1->writeToLog("Third Duration is " + String(stopTime0 - startTime0) + "ms for nSamples = " + String(nSamples));
                            //likelihoodNoSpike = logPlusDetTermForNoiseLL;
                        }

                    }
                }
                // here
            }

        }
        masterSampleIndex += sampleIndex;
         stopTime = Time::getMillisecondCounterHiRes();
        //log1->writeToLog("Duration taken was " + String(stopTime - startTime) + "ms for nSamples = " + String(nSamples));
        //log1->writeToLog("The three durations as percentage are: " + String((stopTime2 - startTime2)*100*464/(stopTime - startTime)) + " AND " +
                         //String((stopTime1 - startTime1)*100*464/(stopTime - startTime)) + " AND " + String((stopTime0 - startTime0)*100*464/(stopTime - startTime)) + " % ");
    }
}
