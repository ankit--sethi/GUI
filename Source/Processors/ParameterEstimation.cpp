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

    ParameterEstimation.cpp
    Created: 15 May 2014 9:29:25am
    Author:  sethisan

  ==============================================================================
*/

#include <stdio.h>
#include <algorithm>
#include "ParameterEstimation.h"
#include "SpikeSorter.h"
#include <fstream>

#include "Channel.h"
#define LOG2PIBY2 0.39908993417
ParameterEstimator::ParameterEstimator()
    : GenericProcessor("Parameter Estimator"),overflowBuffer(2,1), dataBuffer(overflowBuffer), overflowBufferSize(100), currentElectrode(-1), SVDCols(3)
{

        electrodeTypes.clear();
        electrodeCounter.clear();
        electrodeTypes.add("single electrode");
        electrodeTypes.add("stereotrode");
        electrodeTypes.add("tetrode");

    for (int i = 0; i < electrodeTypes.size()+1; i++)
        {
        electrodeCounter.add(0);
        }
        SVDmethod = true;
        cssp = 0;
        delta = 0;
        sampleRate = 10000;
        for (int i = 0; i < 16; i++)
        allParametersEstimated.add(false);
        tempSpike.fill(0);

}

ParameterEstimator::~ParameterEstimator()
{

}

bool ParameterEstimator::disable()
{

    for (int n = 0; n < electrodes.size(); n++)
    {
        resetElectrode(electrodes[n]);
    }

    return true;
}

void ParameterEstimator::updateSettings()
{

    if (getNumInputs() > 0)
        overflowBuffer.setSize(getNumInputs(), overflowBufferSize);

    for (int i = 0; i < electrodes.size(); i++)
    {

        Channel* ch = new Channel(this, i);
        ch->isEventChannel = true;
        ch->eventType = SPIKE_BASE_CODE + electrodes[i]->numChannels;
        ch->name = electrodes[i]->name;

        eventChannels.add(ch);
    }

}

bool ParameterEstimator::enable()
{

    useOverflowBuffer = false;
    ParameterEstimatorEditor* editor = (ParameterEstimatorEditor*) getEditor();
     editor->enable();

    return true;
}

bool ParameterEstimator::addElectrode(int nChans)
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
        firstChan = 0; // make sure we don't overflow available channels
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
    newElectrode->prePeakSamples = 15;
    newElectrode->postPeakSamples = 15;
    newElectrode->paramAveragingCount = 0;
    newElectrode->mu = 0;
    newElectrode->sigma = 0;
    newElectrode->musqrd = 0;
    newElectrode->sumOfSquaresOfDifferences = 0;
    newElectrode->thresholds = new double[nChans];
    newElectrode->isActive = new bool[nChans];
    newElectrode->channels = new int[nChans];
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

float ParameterEstimator::getDefaultThreshold()
{
    return 50.0f;
}

AudioProcessorEditor* ParameterEstimator::createEditor()
{
    editor = new ParameterEstimatorEditor(this, true);
    return editor;
}

void ParameterEstimator::resetElectrode(Electrode* e)
{
    e->lastBufferIndex = 0;
}

bool ParameterEstimator::removeElectrode(int index)
{

    // std::cout << "Spike detector removing electrode" << std::endl;

    if (index > electrodes.size() || index < 0)
        return false;

    electrodes.remove(index);
    return true;
}

void ParameterEstimator::setElectrodeName(int index, String newName)
{
    electrodes[index-1]->name = newName;
}

void ParameterEstimator::setSVDColumnsToUse(int dim)
{
    SVDCols = dim;
}

void ParameterEstimator::setChannel(int electrodeIndex, int channelNum, int newChannel)
{

    std::cout << "Setting electrode " << electrodeIndex << " channel " << channelNum <<
              " to " << newChannel << std::endl;

    *(electrodes[electrodeIndex]->channels+channelNum) = newChannel;
}

int ParameterEstimator::getNumChannels(int index)
{

    if (index < electrodes.size())
        return electrodes[index]->numChannels;
    else
        return 0;
}

int ParameterEstimator::getChannel(int index, int i)
{
    return *(electrodes[index]->channels+i);
}


void ParameterEstimator::setChannelActive(int electrodeIndex, int subChannel, bool active)
{


    currentElectrode = electrodeIndex;
    currentChannelIndex = subChannel;

    std::cout << "Setting channel active to " << active << std::endl;

    if (active)
        setParameter(98, 1);
    else
        setParameter(98, 0);

}


void ParameterEstimator::setChannelThreshold(int electrodeNum, int channelNum, float thresh)
{
    currentElectrode = electrodeNum;
    currentChannelIndex = channelNum;
    std::cout << "Setting electrode " << electrodeNum << " channel threshold " << channelNum << " to " << thresh << std::endl;
    setParameter(99, thresh);
}

double ParameterEstimator::getChannelThreshold(int electrodeNum, int channelNum)
{
    return *(electrodes[electrodeNum]->thresholds+channelNum);
}

void ParameterEstimator::setParameter(int parameterIndex, float newValue)
{
    //editor->updateParameterButtons(parameterIndex);

    if (parameterIndex == 99 && currentElectrode > -1)
    {
        *(electrodes[currentElectrode]->thresholds+currentChannelIndex) = newValue;
    }
    else if (parameterIndex == 98 && currentElectrode > -1)
    {
        if (newValue == 0.0f)
            *(electrodes[currentElectrode]->isActive+currentChannelIndex) = false;
        else
            *(electrodes[currentElectrode]->isActive+currentChannelIndex) = true;
    }
}
StringArray ParameterEstimator::getElectrodeNames()
{
    StringArray names;

    for (int i = 0; i < electrodes.size(); i++)
    {
        names.add(electrodes[i]->name);
    }

    return names;
}

SVDcomputingThread::SVDcomputingThread() : Thread("SVD")
{
    J.reportDone = false;
    for (int i = 0; i < 16; i++)
    dictionary.add(Eigen::MatrixXf::Constant(30, 3, 0));
}

SVDjob::~SVDjob()
{

}

bool ParameterEstimator::samplesAvailable(int& nSamples)
{

    if (sampleIndex >= nSamples - overflowBufferSize/2)
    {
        return false;
    }
    else
    {
        return true;
    }


}

float ParameterEstimator::getNextSample(int& chan)
{


    if (sampleIndex < 0)
    {
        int ind = overflowBufferSize + sampleIndex;

        if (ind < overflowBuffer.getNumSamples())
            return *overflowBuffer.getSampleData(chan, ind);
        else
            return 0;

    }
    else
    {
        if (sampleIndex < dataBuffer.getNumSamples())
            return *dataBuffer.getSampleData(chan, sampleIndex);
        else
            return 0;
    }

}
float ParameterEstimator::getCurrentSample(int& chan)
{

    if (sampleIndex < 1)
    {
        return *overflowBuffer.getSampleData(chan, overflowBufferSize + sampleIndex - 1);
    }
    else
    {
        return *dataBuffer.getSampleData(chan, sampleIndex - 1);
    }

}
bool ParameterEstimator::isChannelActive(int electrodeIndex, int i)
{
    return *(electrodes[electrodeIndex]->isActive+i);
}
void ParameterEstimator::addWaveformToSpikeObject(SpikeObject* s,
                                             int& peakIndex,
                                             int& electrodeNumber,
                                             int& currentChannel)
{
    int spikeLength = electrodes[electrodeNumber]->prePeakSamples +
                      + electrodes[electrodeNumber]->postPeakSamples;

    s->timestamp = timestamp + peakIndex;

    s->nSamples = spikeLength;

    s->eventType = SPIKE_EVENT_CODE;

    int chan = *(electrodes[electrodeNumber]->channels+currentChannel);

    s->gain[currentChannel] = (int)(1.0f / channels[chan]->bitVolts)*1000;
    s->threshold[currentChannel] = (int) *(electrodes[electrodeNumber]->thresholds+currentChannel); // / channels[chan]->bitVolts * 1000;

    // cycle through buffer

    std::ofstream myfile;
    myfile.open("detectedspikes.csv", std::fstream::app);

    if (isChannelActive(electrodeNumber, currentChannel))
    {
        for (int sample = 0; sample < spikeLength; sample++)
        {

            // warning -- be careful of bitvolts conversion
            s->data[currentIndex] = uint16(getNextSample(*(electrodes[electrodeNumber]->channels+currentChannel)) / channels[chan]->bitVolts + 32768);
            myfile << getNextSample(*(electrodes[electrodeNumber]->channels+currentChannel));
            if ( sample != spikeLength - 1)
            {
                myfile << ",";
            }
            currentIndex++;
            sampleIndex++;

            //std::cout << currentIndex << std::endl;

        }
        myfile << std::endl;
    }
    else
    {
        for (int sample = 0; sample < spikeLength; sample++)
        {

            // insert a blank spike if the
            s->data[currentIndex] = 0;
            currentIndex++;
            sampleIndex++;

            //std::cout << currentIndex << std::endl;

        }
    }

    sampleIndex -= spikeLength; // reset sample index


}


#define SIGN(a,b) ((b) > 0.0 ? fabs(a) : - fabs(a))

static double maxarg1,maxarg2;
#define FMAX(a,b) (maxarg1 = (a),maxarg2 = (b),(maxarg1) > (maxarg2) ? (maxarg1) : (maxarg2))

static int iminarg1,iminarg2;
#define IMIN(a,b) (iminarg1 = (a),iminarg2 = (b),(iminarg1 < (iminarg2) ? (iminarg1) : iminarg2))

static double sqrarg;
#define SQR(a) ((sqrarg = (a)) == 0.0 ? 0.0 : sqrarg * sqrarg)

void SVDjob::SVDsetdim(Array<SpikeObject> _spikes, bool _reportDone, int electrodeIndex)
{
    dim = _spikes[0].nChannels*_spikes[0].nSamples;
    spikes = _spikes;
    reportDone = _reportDone;
    SVDjob::electrodeIndex = electrodeIndex;
}

float SVDjob::pythag(float a, float b) {
  float absa,absb;

  absa = fabs(a);
  absb = fabs(b);

  if(absa > absb)
    return(absa * sqrt(1.0 + SQR(absb/absa)));
  else
    return(absb == 0.0 ? 0.0 : absb * sqrt(1.0 + SQR(absa / absb)));
}
int SVDjob::svdcmp(float **a, int nRows, int nCols, float *w, float **v) {


  int flag,i,its,j,jj,k,l,nm;
  float anorm,c,f,g,h,s,scale,x,y,z,*rv1;

  rv1 = new float[nCols];
  if(rv1 == NULL) {
    printf("svdcmp(): Unable to allocate vector\n");
    return(-1);
  }

  g = scale = anorm = 0.0;
  for(i=0;i<nCols;i++) {
    l = i+1;
    rv1[i] = scale*g;
    g = s = scale = 0.0;
    if(i < nRows) {
      for(k=i;k<nRows;k++) scale += fabs(a[k][i]);
      if(scale) {
    for(k=i;k<nRows;k++) {
      a[k][i] /= scale;
      s += a[k][i] * a[k][i];
    }
    f = a[i][i];
    g = -SIGN(sqrt(s),f);
    h = f * g - s;
    a[i][i] = f - g;
    for(j=l;j<nCols;j++) {
      for(s=0.0,k=i;k<nRows;k++) s += a[k][i] * a[k][j];
      f = s / h;
      for(k=i;k<nRows;k++) a[k][j] += f * a[k][i];
    }
    for(k=i;k<nRows;k++) a[k][i] *= scale;
      }
    }
    w[i] = scale * g;
    g = s = scale = 0.0;
    if(i < nRows && i != nCols-1) {
      for(k=l;k<nCols;k++) scale += fabs(a[i][k]);
      if(scale)  {
    for(k=l;k<nCols;k++) {
      a[i][k] /= scale;
      s += a[i][k] * a[i][k];
    }
    f = a[i][l];
    g = - SIGN(sqrt(s),f);
    h = f * g - s;
    a[i][l] = f - g;
    for(k=l;k<nCols;k++) rv1[k] = a[i][k] / h;
    for(j=l;j<nRows;j++) {
      for(s=0.0,k=l;k<nCols;k++) s += a[j][k] * a[i][k];
      for(k=l;k<nCols;k++) a[j][k] += s * rv1[k];
    }
    for(k=l;k<nCols;k++) a[i][k] *= scale;
      }
    }
    anorm = FMAX(anorm, (fabs(w[i]) + fabs(rv1[i])));


  }

  for(i=nCols-1;i>=0;i--) {
    if(i < nCols-1) {
      if(g) {
    for(j=l;j<nCols;j++)
      v[j][i] = (a[i][j] / a[i][l]) / g;
    for(j=l;j<nCols;j++) {
      for(s=0.0,k=l;k<nCols;k++) s += a[i][k] * v[k][j];
      for(k=l;k<nCols;k++) v[k][j] += s * v[k][i];
    }
      }
      for(j=l;j<nCols;j++) v[i][j] = v[j][i] = 0.0;
    }
    v[i][i] = 1.0;
    g = rv1[i];
    l = i;
  }

  for(i=IMIN(nRows,nCols) - 1;i >= 0;i--) {
    l = i + 1;
    g = w[i];
    for(j=l;j<nCols;j++) a[i][j] = 0.0;
    if(g) {
      g = 1.0 / g;
      for(j=l;j<nCols;j++) {
    for(s=0.0,k=l;k<nRows;k++) s += a[k][i] * a[k][j];
    f = (s / a[i][i]) * g;
    for(k=i;k<nRows;k++) a[k][j] += f * a[k][i];
      }
      for(j=i;j<nRows;j++) a[j][i] *= g;
    }
    else
      for(j=i;j<nRows;j++) a[j][i] = 0.0;
    ++a[i][i];
  }

  for(k=nCols-1;k>=0;k--) {
    for(its=0;its<30;its++) {
      flag = 1;
      for(l=k;l>=0;l--) {
    nm = l-1;
    if((fabs(rv1[l]) + anorm) == anorm) {
      flag =  0;
      break;
    }
    if((fabs(w[nm]) + anorm) == anorm) break;
      }
      if(flag) {
    c = 0.0;
    s = 1.0;
    for(i=l;i<=k;i++) {
      f = s * rv1[i];
      rv1[i] = c * rv1[i];
      if((fabs(f) + anorm) == anorm) break;
      g = w[i];
      h = pythag(f,g);
      w[i] = h;
      h = 1.0 / h;
      c = g * h;
      s = -f * h;
      for(j=0;j<nRows;j++) {
        y = a[j][nm];
        z = a[j][i];
        a[j][nm] = y * c + z * s;
        a[j][i] = z * c - y * s;
      }
    }
      }
      z = w[k];
      if(l == k) {
    if(z < 0.0) {
      w[k] = -z;
      for(j=0;j<nCols;j++) v[j][k] = -v[j][k];
    }
    break;
      }
      //if(its == 29) printf("no convergence in 30 svdcmp iterations\n");
      x = w[l];
      nm = k-1;
      y = w[nm];
      g = rv1[nm];
      h = rv1[k];
      f = ((y - z) * (y + z) + (g - h) * (g + h)) / (2.0 * h * y);
      g = pythag(f,1.0);
      f = ((x - z) * (x + z) + h * ((y / (f + SIGN(g,f))) - h)) / x;
      c = s = 1.0;
      for(j=l;j<=nm;j++) {
    i = j+1;
    g = rv1[i];
    y = w[i];
    h = s * g;
    g = c * g;
    z = pythag(f,h);
    rv1[j] = z;
    c = f/z;
    s = h/z;
    f = x * c + g * s;
    g = g * c - x * s;
    h = y * s;
    y *= c;
    for(jj=0;jj<nCols;jj++) {
      x = v[jj][j];
      z = v[jj][i];
      v[jj][j] = x * c + z * s;
      v[jj][i] = z * c - x * s;
    }
    z = pythag(f,h);
    w[j] = z;
    if(z) {
      z = 1.0 / z;
      c = f * z;
      s = h * z;
    }
    f = c * g + s * y;
    x = c * y - s * g;
    for(jj=0;jj < nRows;jj++) {
      y = a[jj][j];
      z = a[jj][i];
      a[jj][j] = y * c + z * s;
      a[jj][i] = z * c - y * s;
    }
      }
      rv1[l] = 0.0;
      rv1[k] = f;
      w[k] = x;
    }
  }

  delete rv1;

  return(0);
}
void SVDcomputingThread::run()
{
    while (jobs.size() > 0)
    {

        J = jobs.front();
        jobs.pop();


        float **eigvec, *sigvalues, **U;
        sigvalues = new float[J.dim];
        eigvec = new float*[J.dim];
        U = new float*[J.spikes.size()];



        for (int k=0;k<J.dim;k++) {
            eigvec[k] = new float[J.dim];
            for (int j=0;j<J.dim;j++)
            {
                eigvec[k][j] = 0;
            }
        }


        for (int k = 0; k < J.spikes.size(); k++)
        {
            U[k] = new float[J.dim];
            for (int j = 0; j < J.dim; j++)
            {
                SpikeObject spike = J.spikes[k];
                float v = spikeDataIndexToMicrovolts(&spike, j) ;

                U[k][j] = v;

            }
        }

        J.svdcmp(U,J.spikes.size(),J.spikes[1].nSamples,sigvalues,eigvec);

        Eigen::Matrix<float, 30, 3> tempmatrix;
        tempmatrix.fill(0);

        std::ofstream myfile;
        myfile.open("dictionary.csv", std::fstream::app);
        for (int k = 0; k < dictionary[J.electrodeIndex].rows(); k++)
        {
            for (int j = 0; j < dictionary[J.electrodeIndex].cols(); j++)
            {
                myfile << eigvec[k][j];
                tempmatrix(k,j) = eigvec[k][j];

                if (j != J.dim - 1)
                {
                    myfile << ",";
                }
            }
            myfile << std::endl;
        }

        dictionary.setUnchecked(J.electrodeIndex, tempmatrix);
        //std::cout << "DICTIONARY BEING MADE" << dictionary[J.electrodeIndex];
        //clear memory
        for (int k=0;k<J.dim;k++)
        {
            delete eigvec[k];
        }
        for (int k=0;k<J.spikes.size();k++)
        {
            delete U[k];
        }
        delete U;
        delete sigvalues;
        delete eigvec;
        (J.reportDone) = true;

    }
}
void SVDcomputingThread::addSVDjob(SVDjob job)
{
    jobs.push(job);
    if (!isThreadRunning())
    {
        startThread();
        //std::cout<<"Thread has started"<<std::endl;
    }
}

void ParameterEstimator::handleEvent(int eventType, MidiMessage& event, int sampleNum)
{

    if (eventType == TIMESTAMP)
    {
        const uint8* dataptr = event.getRawData();

        memcpy(&timestamp, dataptr + 4, 8); // remember to skip first four bytes
    }


}
float ParameterEstimator::getElectrodeNoiseVar(int index)
{
    return ((electrodes[index]->sigma)*(electrodes[index]->sigma));
}

void ParameterEstimator::process(AudioSampleBuffer& buffer,
                               MidiBuffer& events,
                               int& nSamples)
{

    Electrode* electrode;
    dataBuffer = buffer;
    // TBD if one needs only one SVD over all electrodes, or one per electrode



    checkForEvents(events);

    // this marks beginning of parameter detection --------------------------
    for(int i = 0; i < electrodes.size(); i++)
    {
        electrode = electrodes[i];
        sampleIndex = electrode->lastBufferIndex;

        while (samplesAvailable(nSamples))
        {
            electrode->paramAveragingCount++;

            for(int chan=0; chan < 1 ; chan++) //detecting on 1 channel only currently. to be expanded. CUSUM -> multichannel/distributed CUSUM
            {
                if (*(electrode->isActive+chan))
                {

                    int currentChannel = *(electrode->channels+chan);

                    double var = getNextSample(currentChannel);


                    //std::cout<< var << std::endl;
                    if (electrode->paramAveragingCount < sampleRate*30) //averaging for 30 s to find noise mu and sigma (needs to be changed to that each electrode gets its own noise measurement)
                    {
                        //electrode->paramAveragingCount++;
                        delta = (var - electrode->mu);
                        electrode->mu = electrode->mu + (delta/electrode->paramAveragingCount);
                        electrode->sumOfSquaresOfDifferences = electrode->sumOfSquaresOfDifferences + delta*(var - electrode->mu);


                    }
                    if (electrode->paramAveragingCount == int(sampleRate*30))
                    {
                        electrode->musqrd = (electrode->mu)*(electrode->mu);
                        electrode->sigma = std::sqrt(electrode->sumOfSquaresOfDifferences/electrode->paramAveragingCount); //DANGER?
                        std::cout<<"End of first 30 seconds" << std::endl;
                        std::cout << "The mu and sigma for ripple detection are - " << electrode->mu << " and " << electrode->sigma << std::endl;
                    }
                    if (electrode->paramAveragingCount >= sampleRate*60 && electrode->paramAveragingCount < sampleRate*120) // START PRELIM SPIKE DETECTION
                    {
                        if (electrode->paramAveragingCount == (int)(sampleRate*60))
                        {
                            std::cout<< "Starting spike detection." << std::endl;
                        }
                        if (var > 5*electrode->sigma) // trigger spike
                        {

                            //std::cout << "Spike detected on electrode " << i << std::endl;
                            // find the peak
                            int peakIndex = sampleIndex;
                            for (int i = 0 ; i < 5; i++)
                            {
                                tempSpike(i) = getCurrentSample(currentChannel);
                                sampleIndex++;
                            }
                            int offset;
                            int non = tempSpike.maxCoeff(&offset);
                            peakIndex += offset - 1;
                            //sampleIndex -= 15;

                            //while (-getCurrentSample(currentChannel) < -getNextSample(currentChannel))//&& sampleIndex < peakIndex + electrode->postPeakSamples)
                            {
                                //sampleIndex++;
                            }
                            //--sampleIndex;

                            sampleIndex = peakIndex - (electrode->prePeakSamples-1);

                            SpikeObject newSpike;
                            newSpike.timestamp = peakIndex;
                            newSpike.source = i;
                            newSpike.nChannels = 1;  // HARDCODING TO 1

                            currentIndex = 0;

                            // package spikes;
                            for (int channel = 0; channel < 1; channel++)
                            {
                                addWaveformToSpikeObject(&newSpike, peakIndex, i, channel);
                            }


                            //add spikes to an array here
                            //if (SVDmethod == true)
                            {
                                detectedSpikesAllElectrodes.add(newSpike);
                            }
                            //else
                            {
                                //detectedSpikesPerElectrode[i].add(newSpike);
                            }

                            // advance the sample index
                            sampleIndex = peakIndex + electrode->postPeakSamples + 10;

                            break; // quit spike "for" loop
                        } // end spike trigger
                        //std::cout << "Spike detected on electrode " << i << std::endl;
                    }

                    if (electrode->paramAveragingCount < sampleRate*60 && electrode->paramAveragingCount >= sampleRate*30)
                    {
                        sampleIndex += 4;
                        double var0 = getNextSample(currentChannel);
                        sampleIndex -= 4;
                        cssp = cssp + (var0*var - (electrode->mu)*(var0 + var) + electrode->musqrd);
                    }
                    if (electrode->paramAveragingCount == int(sampleRate*60))
                    {
                        acfLag1 = cssp/(electrode->sigma*electrode->sigma);
                        acfLag1 = acfLag1/(sampleRate*30);
                        std::cout << " The ACF LAG1 is - " << acfLag1 << std::endl;
                    }
                }
            }
            sampleIndex +=4;
        }
        //this marks the ending of parameter detection -----------------------------

        if (electrode->paramAveragingCount >= (sampleRate*120) && !allParametersEstimated[i])
        {

            if (SVDmethod) // starting the SVD Thread
            {
                std::cout << "A total of " << detectedSpikesAllElectrodes.size() << "spikes have been detected and stored." << std::endl;
                std::cout<<detectedSpikesAllElectrodes.size()<<" spikes objects stored"<<std::endl;
                SVDjob Tjob;
                Tjob.SVDsetdim(detectedSpikesAllElectrodes, false, i);
                job.add(Tjob);
                std::cout<<"job was created" << std::endl;
                dictionaryThread.addSVDjob(job[i]);
                SVDmethod = false;
            }
            if(!allParametersEstimated[i])
            {

                if(dictionaryThread.J.reportDone && i == dictionaryThread.J.electrodeIndex)
                {
                    std::cout << "done with svd";
                    ProcessorGraph* gr = getProcessorGraph();
                    juce::Array<GenericProcessor*> p = gr->getListOfProcessors();

                    bool flag = false;
                    for (int k=0;k<p.size();k++)
                    {
                        if (p[k]->getName() == "Spike Sorter")
                        {
                            node = (SpikeSorter*)p[k];
                            flag = true;
                        }

                    }
                    if (!flag)
                    {
                        std::cout << "Could not find a the Spike Sorter." << std::endl;
                    }
                    node->setParameter(3,i);


                    std::cout<<"Parameter Estimator contacted and Circular Buffers set. Buffer size and electrode number is " << node->se[i].sec[0].circbuffer.getsize() << std::endl;

                    Eigen::ArrayXf powers = Eigen::VectorXf::Zero(node->P).array();
                    double start = 1;

                    for (int i = 0; i < node->P; i++)
                    {
                        powers(i) = start;
                        start *= acfLag1;
                    }

                    float powerIndex;
                    for (int k = 0; k < node->P; k++)
                    {
                        powerIndex = -1*i;
                        for (int j = 0; j < node->P; j++)
                        {

                            if ( k == j )
                                node->se[i].sigma(k,j) = getElectrodeNoiseVar(0);
                            else
                                node->se[i].sigma(k,j) = getElectrodeNoiseVar(0)*powers(std::abs(powerIndex));

                            powerIndex++;
                        }
                    }

                    std::cout<< " ACF LAG in Spike Sorter is " << acfLag1 << " and the Var is " << getElectrodeNoiseVar(0) << "  "  << std::endl;
                    double startTime1 = Time::getMillisecondCounterHiRes();
                    node->se[i].lambda = node->se[i].sigma.inverse();
                    double stopTime1 = Time::getMillisecondCounterHiRes();
                    Logger *log1 = Logger::getCurrentLogger();
                    log1->writeToLog("Time for 30 x 30 inverse is = " + String(stopTime1 - startTime1));
                    node->se[i].lambdaQR.compute(node->se[i].lambda);
                    node->se[i].logPlusDetTermForNoiseLL = (-1*node->P)*LOG2PIBY2 + 0.5*node->se[i].lambdaQR.logAbsDeterminant();
                    std::cout<<"Determinant is = " << node->se[i].lambdaQR.logAbsDeterminant() << " and " << node->se[i].ReducedDictionary.col(0).size() << "x" << node->se[i].ReducedDictionary.row(0).size() << "//" << std::endl;
                    node->setParameter(2,i);
                    allParametersEstimated.setUnchecked(i,true);
                    std::cout << " All parameters estimated for electrode " << i << "." << std::endl;
                    bool test = true;
                    for ( int k = 0; k < electrodes.size(); k++)
                    {
                        test = test && allParametersEstimated[k];
                    }
                    if (test)
                    node->allParametersEstimated = true;

                }
            }
        }
        else
        {

        }
        electrode->lastBufferIndex = sampleIndex - nSamples;



    }
    if (nSamples > overflowBufferSize)
    {

        for (int i = 0; i < overflowBuffer.getNumChannels(); i++)
        {

            overflowBuffer.copyFrom(i, 0,
                                    buffer, i,
                                    nSamples-overflowBufferSize,
                                    overflowBufferSize);

            useOverflowBuffer = true;
        }

    }
    else
    {
        useOverflowBuffer = false;
    }
}

void ParameterEstimator::loadCustomParametersFromXml()
{

}

void ParameterEstimator::saveCustomParametersToXml(XmlElement *electrodeNode)
{


}
