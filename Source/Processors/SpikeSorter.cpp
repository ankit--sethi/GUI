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


using namespace Eigen;


#define PI 3.1415926535897932384626433832795


// calculate the cofactor of element (row,col)
void SpikeSorter::GetMinor(long double **src, long double **dest, int row, int col, int order)
{
    // indicate which col and row is being copied to dest
    int colCount=0,rowCount=0;

    for(int i = 0; i < order; i++ )
    {
        if( i != row )
        {
            colCount = 0;
            for(int j = 0; j < order; j++ )
            {
                // when j is not the element
                if( j != col )
                {
                    dest[rowCount][colCount] = src[i][j];
                    colCount++;
                }
            }
            rowCount++;
        }
    }
}
void SpikeSorter::GetMinor(double **src, double **dest, int row, int col, int order)
{
    // indicate which col and row is being copied to dest
    int colCount=0,rowCount=0;

    for(int i = 0; i < order; i++ )
    {
        if( i != row )
        {
            colCount = 0;
            for(int j = 0; j < order; j++ )
            {
                // when j is not the element
                if( j != col )
                {
                    dest[rowCount][colCount] = src[i][j];
                    colCount++;
                }
            }
            rowCount++;
        }
    }
}


// Calculate the determinant recursively.
double SpikeSorter::CalcDeterminant(long double **in_matrix, int n)
{
   int i, j, k;
  long double **matrix;
  double det = 1;

  matrix = new long double*[n];

  for ( i = 0; i < n; i++ )
    matrix[i] = new long double[n];

  for ( i = 0; i < n; i++ ) {
    for ( j = 0; j < n; j++ )
      matrix[i][j] = in_matrix[i][j];
  }

  for ( k = 0; k < n; k++ ) {
    if ( matrix[k][k] == 0 ) {
      bool ok = false;

      for ( j = k; j < n; j++ ) {
        if ( matrix[j][k] != 0 )
          ok = true;
      }

      if ( !ok )
        return 0;

      for ( i = k; i < n; i++ )
        std::swap ( matrix[i][j], matrix[i][k] );

      det = -det;
    }

    det *= matrix[k][k];

    if ( k + 1 < n ) {
      for ( i = k + 1; i < n; i++ ) {
        for ( j = k + 1; j < n; j++ )
          matrix[i][j] = matrix[i][j] - matrix[i][k] *
          matrix[k][j] / matrix[k][k];
      }
    }
  }

  for ( i = 0; i < n; i++ )
    delete [] matrix[i];

  delete [] matrix;

  return det;
}

double SpikeSorter::CalcDeterminant(double **in_matrix, int n)
{
   int i,j,j1,j2;
   double det = 0;
   double **m = NULL;

   if (n < 1) { /* Error */

   } else if (n == 1) { /* Shouldn't get used */
      det = in_matrix[0][0];
   } else if (n == 2) {
      det = in_matrix[0][0] * in_matrix[1][1] - in_matrix[1][0] * in_matrix[0][1];
   } else {
      det = 0;
      for (j1=0;j1<n;j1++) {
         m = (double**)malloc((n-1)*sizeof(double *));
         for (i=0;i<n-1;i++)
            m[i] = (double*)malloc((n-1)*sizeof(double));
         for (i=1;i<n;i++) {
            j2 = 0;
            for (j=0;j<n;j++) {
               if (j == j1)
                  continue;
               m[i-1][j2] = in_matrix[i][j];
               j2++;
            }
         }
         det += pow(-1.0,1.0+j1+1.0) * in_matrix[0][j1] * CalcDeterminant(m,n-1);
         for (i=0;i<n-1;i++)
            free(m[i]);
         free(m);
      }
   }
   return(det);
}
double SpikeSorter::CalcDeterminant(float **in_matrix, int n)
{
   int i,j,j1,j2;
   double det = 0;
   float **m = NULL;

   if (n < 1) { /* Error */

   } else if (n == 1) { /* Shouldn't get used */
      det = in_matrix[0][0];
   } else if (n == 2) {
      det = in_matrix[0][0] * in_matrix[1][1] - in_matrix[1][0] * in_matrix[0][1];
   } else {
      det = 0;
      for (j1=0;j1<n;j1++) {
         m = (float**)malloc((n-1)*sizeof(float *));
         for (i=0;i<n-1;i++)
            m[i] = (float*)malloc((n-1)*sizeof(float));
         for (i=1;i<n;i++) {
            j2 = 0;
            for (j=0;j<n;j++) {
               if (j == j1)
                  continue;
               m[i-1][j2] = in_matrix[i][j];
               j2++;
            }
         }
         det += pow(-1.0,1.0+j1+1.0) * in_matrix[0][j1] * CalcDeterminant(m,n-1);
         for (i=0;i<n-1;i++)
            free(m[i]);
         free(m);
      }
   }
   return(det);
}
void SpikeSorter::MatrixInversion(long double **A, int order, long double **Y)
{
    /*bool flag = false;
    int count = 0;

    for (int i = 0; i < order; i++)
    {
        for(int j = 0; j < order; j++)
        {
            if (A[i][j] != A[i][j])
            {
                flag = true;
                count++;
                //test = false;
            }
        }
    }
    if (flag)
    {
    //std::cout << "Matrix arrived NAN with " << count << " elements." << std::endl;
    return;
    }*/

    long double tempest = CalcDeterminant(A,order);
    long double det = 1.0/tempest;

    long double *temp = new long double[(order-1)*(order-1)];
    long double **minor = new long double*[order-1];
    for(int i=0;i<order-1;i++)
        minor[i] = temp+(i*(order-1));
    for(int j=0;j<order;j++)
    {
        for(int i=0;i<order;i++)
        {
            // get the co-factor (matrix) of A(j,i)
            GetMinor(A,minor,j,i,order);
            Y[i][j] = det*CalcDeterminant(minor,order-1);
            if (Y[i][j] != Y[i][j])
            {
                std::cout << "possible nan";
            }
            if( (i+j)%2 == 1) // in case of zero, values become negative zero
                Y[i][j] = -Y[i][j];
        }
    }

    delete [] temp;
    delete [] minor;
}

void SpikeSorter::MatrixInversion(double **A, int order, double **Y)
{

    /*bool flag = false;
    int count = 0;

    for (int i = 0; i < order; i++)
    {
        for(int j = 0; j < order; j++)
        {
            if (A[i][j] != A[i][j])
            {
                flag = true;
                count++;
            }
        }
    }
    if (flag)
    {
    //std::cout << "Matrix arrived NAN with " << count << " elements." << std::endl;
    }*/


    double tempest = CalcDeterminant(A,order);
    double det = 1.0/tempest;

    double *temp = new double[(order-1)*(order-1)];
    double **minor = new double*[order-1];
    for(int i=0;i<order-1;i++)
        minor[i] = temp+(i*(order-1));

    for(int j=0;j<order;j++)
    {
        for(int i=0;i<order;i++)
        {
            // get the co-factor (matrix) of A(j,i)

            GetMinor(A,minor,j,i,order);
            Y[i][j] = det*CalcDeterminant(minor,order-1);
            if (Y[i][j] != Y[i][j])
            {
                std::cout << "possible nan";
            }
            if( (i+j)%2 == 1) // in case of zero, values become negative zero
                Y[i][j] = -Y[i][j];
        }
    }

    delete [] temp;
    delete [] minor;
}
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

VectorXf CircularQueue::show(int index)
{
    Array <float> xwithPsamples;
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
    phi0 = MatrixXd::Zero(K,K);
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
    range = 40; // defines basically the width of a spike
    // setting sigma
    sigma = MatrixXd::Zero(P,P);

    pii=apii/bpii;

    MatrixXd phi0invnu0forR = ((phi0.inverse().array() + nu0)*2).matrix();
    nu = VectorXd::Constant(1,Cmax, nu0);
    kappa = VectorXd::Constant(1,Cmax, kappa0);
    ltheta = VectorXd::Constant(1,Cmax, 0);
    ngamma = VectorXd::Constant(1,Cmax, 0);
    tlastspike = VectorXd::Constant(1,Cmax, 0);
    likelihoodPerNeuron = VectorXd::Constant(1,Cmax, 0);

    for (int i = 0; i < Cmax; i++)
    {
        phi.add(phi0);
        lamclus.add(phi0invnu0);
        lamclusinv.add(phi0invnu0);
        R.add(phi0invnu0forR);
        Rinv.add(phi0invnu0forR);
    }
    muu = MatrixXd::Constant(K, Cmax, 0);
    muu0 = MatrixXd::Constant(K, Cmax, 0);
    Q = MatrixXd::Constant(P, P, 0);
    Qinv = MatrixXd::Constant(P, P, 0);
    Qhat = MatrixXd::Constant(K, K, 0);
    Qmat = MatrixXd::Constant(K, K, 0);
    mhat = VectorXd::Constant(K, 1, 0);
    yhat = VectorXd::Constant(K, 1, 0);
    ReducedDictionary = Matrix::Constant(P, K, 0);
    ReducedDictionaryTranspose = ReducedDictionary.transposeInPlace();

    neuronCount = 0; // This is 'C' in the MATLAB code.
    //nz = 0; // Not sure what this does. Possible the same thing.
    setParameter(1,0);

    threshold = log(pii/(1-pii));
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

    if (parameterIndex == 1 ) //
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
    if (parameterIndex == 2 )
    {
        for (int i = 0; i < P; i++)
        {
            float power = -1*i;
            for (int j = 0; j <P; j++)
            {

                sigma[i].set(j,node->getElectrodeNoiseVar(1)*pow(node->acfLag1, power++)); // NEEDS TO BE MULTIPED WITH COV(X);

            }
        }
        for (int i = 0; i < P; i++)
        {
            for (int j = 0; j <P; j++)
            {
                if (i==j)
                    sigma[i].set(j,node->getElectrodeNoiseVar(1));
            }
        }

    }

    if (parameterIndex == 3 )
    {

    }

    if (parameterIndex == 4 ) // this exists to get the SVD
    {
        float temp; Array<float> temparray;
        for (int i = 0; i < P; i++)
        {
            for (int j = 0; j < K; j++)
            {
                temp = node->dictionary.get(i,j);
                temparray.add(temp);
            }
        }
         MatrixTranspose(ReducedDictionary, ReducedDictionaryTranspose);

    }

    if (parameterIndex == 5)
    {
        ProcessorGraph* gr = getProcessorGraph();
        Array<GenericProcessor*> p = gr->getListOfProcessors();

        bool flag = false;
        for (int k=0;k<p.size();k++)
        {
            if (p[k]->getName() == "Parameter Estimator")
            {
                node = (ParameterEstimator*)p[k];
                checkIfAllParametersEstimated = node->allParametersEstimated;
                //P = node->dictionary.getcoldim();
                flag = true;
            }

        }
        if (!flag)
        {
            std::cout << "Could not find a Parameter Estimator." << std::endl;
        }
    }
    if (parameterIndex == 6)
    {


    }

}

float SpikeSorter::getNextSample(int& chan)
{
        if (sampleIndex < dataBuffer.getNumSamples())
            return *dataBuffer.getSampleData(chan, sampleIndex);
        else
            return 0;
}



void SpikeSorter::getLikelihoodPerNeuron(int FoundNeuronIndex, float tmp)
{
    //std::cout << "The inner product here is " << getaBaInnerProductSum(xwind, Qinv);

}

int SpikeSorter::findNeuronID()
{
        float min = 0;
        idx = 0;
        min = lthr[0];
        for (int i= 0; i < lthr.size(); i++)
        {
            if (min < lthr[i])
            {
                min = lthr[i];
                idx = i;
            }
        }
        float max = 0, idx1 = 0;
        max = lon[idx][0];
        for (int i= 0; i < lon[idx].size(); i++)
        {
            if (max < lon[idx][i])
            {
                max = lon[idx][i];
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
    lthr = VectorXd::Zero(range);
    lon = MatrixXd::Zero(range, likelihoodPerNeuron.size());
    xwindLonger = VectorXd::Zero(xwind.size() + range);
    xwindLonger.head(xwind.size()) = xwind;

    currentSpike.data[currentIndex++] = uint16(trigSample/ channels[chan]->bitVolts + 32768);
}

void SpikeSorter::addNewSampleAndLikelihoodsToCurrentSpikeObject(float sample, MidiBuffer& eventBuffer, int chan)
{
    if (currentIndex < range)
    {

        currentSpike.data[currentIndex++] = uint16(sample/ channels[chan]->bitVolts + 32768);
        lthr(currentIndex) = likelihoodThreshold;


        for (int i = 0; i< likelihoodPerNeuron.size(); i++)
        {
            lon(currentIndex, i) = likelihoodPerNeuron[i];
        }
        xwindLonger(xwind.size() + currentIndex) = sample;

    }
    else
    {

        spikeDetected = false;
        samplesBeingCollected = false;
        currentIndex = 0;
        currentSpike.neuronID = findNeuronID();

        PackageCurrentSortedSpikeIntoBuffer(eventBuffer);
        //std:: cout << int(currentSpike.eventType) << "//" << currentSpike.gain[1] << "//" << currentSpike.nChannels << "//" << currentSpike.neuronID << "//" << currentSpike.nSamples << "//" <<
        //currentSpike.source << "//" << currentSpike.threshold[1] << "//" << currentSpike.timestamp << std::endl;

        int numBytes = packSortedSpike(currentSpike, spikeBuffer, MAX_SPIKE_BUFFER_LEN);
        if (numBytes > 0)
        eventBuffer.addEvent(spikeBuffer, numBytes, int(currentSpike.timestamp));

        //clear everything that needs to be cleared now! we look for new spikes

        updateAllSortingParameters();
        std::cout << "One Spike Handled!" << std::endl;
        totalSpikesFound += 1;

    }
}


void SpikeSorter::updateAllSortingParameters()
{
    int neuronID = currentSpike.neuronID;
    lambda = sigma.inverse();
    Qmat = ReducedDictionaryTranspose*lambda*ReducedDictionary + lamclus(neuronID);

    HouseholderQR<MatrixXd> QmatQR(Qmat);
    yhat = QmatQR.solve(ReducedDictionaryTranspose*lambda*xwindLonger.segment(idx,idx + P - 1)) + lamclus(neuronID)*muu.row(neuronID);

    ngamma(neuronID) = ngamma(neuronID) + 1;
    deltaT = masterSampleIndex + sampleIndex -P - range + idx - tlastspike(neuronID);
    tlastspike(neuronID) = masterSampleIndex + sampleIndex - P - range + idx;
    double ebet = exp(-1*double(beta)*double(deltaT));
    //std::cout << "Ebet and delta T is" << ebet << "//" << deltaT << "//";

    //MatrixGetRowOrColumn(muu,false, neuronID,t);
    mhat = (muu0.row(neuronID).array()*(1-ebet) + muu.row(neuronID).array()*ebet).matrix();
    muu0.row(neuronID) = (kappa(neuronID)*muu0(neuronID).array() + yhat.array())/(kappa(neuronID)+1).matrix();
    Qhat = (itau*MatrixXd::Identity(K).array()*(1-ebet*ebet) + R(neuronID).inverse().array()*(ebet*ebet)).matrix();
    R(neuronID) = Qhat.inverse() + lamclus(neuronID);
    muu.row(neuronID) = R(neuronID).HouseholderQR().solve(Qhat.HouseholderQR().solve(mhat) + lamclus(neuronID)*yhat);
    yhat = Qmat.HouseholderQR().solve(ReducedDictionaryTranspose*lambda*xwindLonger.segment(idx,idx+ P-1) + lamclus(neuronID)*muu.row(neuronID));
    phi(neuronID) = phi(neuronID) + ((kappa(neuronID)/(kappa(neuronID)+1))*((yhat - muu(neuronID))*(yhat - muu(neuronID)).transpose()).array()).matrix() + Qmat.inverse;
    kappa(neuronID) = kappa(neuronID) + 1;
    nu(neuronID) = nu(neuronID) + 1;
    lamclus(neuronID) = phi.HouseholderQR().solve(nu(neuronID));
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
        Array<GenericProcessor*> p = gr->getListOfProcessors();

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
            Array<double> powers;
            double start = 1;

            for (int i = 0; i < P; i++)
            {
                powers.add(start);
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
                        sigma(i,j) = node->getElectrodeNoiseVar(0)*powers.getUnchecked(abs(powerIndex));
                    powerIndex++;
                }
            }


            std::cout<< " ACF LAG in Spike Sorter is " << node->acfLag1 << " and the Var is " << node->getElectrodeNoiseVar(0) << "  "  << std::endl;
            //Array<Array<long double>> proxy;


            //double startTime1 = Time::getMillisecondCounterHiRes();
            invertSquareTwoDimMatrix(sigma, proxy);
            //double stopTime1 = Time::getMillisecondCounterHiRes();
            //log1->writeToLog("Time for 40 x 40 inverse is = " + String(stopTime1 - startTime1));

            lambda = sigma.inverse();

            logDeterminantOfLambda = lambda.determinant();

            std::cout<<"reached line 3 and det is" << logDeterminantOfLambda << ReducedDictionary.size() << "x" << ReducedDictionary[0].size() << "//" << std::endl;

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
                        threshold = log(pii/(1-pii));


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
                                    ltheta(i) = (log(ngamma(i)/(alpha+totalSpikesFound)));

                                }
                                else
                                {
                                    ltheta(i) = -1*PI;
                                }

                            }
                            ltheta(neuronCount) = log(alpha/(alpha+totalSpikesFound));

                            HouseholderQR<MatrixXd> sigmaQR(sigma);

                            double Quad = (xwind*sigmaQR.householderQr().solve(xwind)).squaredNorm();

                            likelihoodNoSpike = (-1*P/2)*log(2*PI) + 0.5*logDeterminantOfLambda - 0.5*Quad;

                            for (int j = 0; j < neuronCount; j++)
                            {
                                HouseholderQR<MatrixXd> lamclusQR(lamclus[j]);
                                HouseholderQR<MatrixXd> RQR(R[j]);
                                Q = sigma + ReducedDictionary*(RQR.solve(ReducedDictionaryTranspose) + lamclusQR.solve(ReducedDictionaryTranspose));
                                xwindloop = xwind - ReducedDictionary*muu.col(j);
                                double sum = Q.householderQr().logAbsDeterminant();


                                HouseholderQR<MatrixXd> QQR(Q);

                                if( ( (masterSampleIndex + sampleIndex) - tlastspike(j)) < 40000*50/10000 )  // THIS NEEDS TO BE INVESTIGATED
                                    likelihoodPerNeuron(FoundNeuronIndex) = -(0.5*P)*log(2*PI) - tmp - suppresslikelihood - 0.5*(xwindloop*QQR.solve(xwindloop)).squaredNorm();
                                else
                                    likelihoodPerNeuron(FoundNeuronIndex) = -(0.5*P)*log(2*PI) - tmp - 0.5*(xwindloop*QQR.solve(xwindloop)).squaredNorm();
                            }
                            HouseholderQR<MatrixXd> lamclusQR(lamclus[neuronCount]);
                            HouseholderQR<MatrixXd> RQR(R[neuronCount]);
                            Q = sigma + ReducedDictionary*(RQR.solve(ReducedDictionaryTranspose) + lamclusQR.solve(ReducedDictionaryTranspose));
                            xwindloop = xwind - ReducedDictionary*muu.col(j);
                            double sum = Q.householderQr().logAbsDeterminant();
                            HouseholderQR<MatrixXd> QQR(Q);
                            likelihoodPerNeuron(neuronCount) = -(0.5*P)*log(2*PI) - sum - 0.5*(xwind*QQR.solve(xwind)).squaredNorm();


                            for (int i = 0; i <= neuronCount; i++ )
                            {
                                likelihoodPerNeuron(i) =  likelihoodPerNeuron(i) + ltheta(i);
                                //std::cout << "And " << likelihoodPerNeuron[i] << " is the Likelihood per Neuron for " << i << " and " << " ltheta is " << ltheta[i].getVar() << std::endl;
                            }
                            for (int i = 0; i <= neuronCount; i++ )
                            {
                                if (i == 0)
                                    maximumLikelihoodPerNeuron = likelihoodPerNeuron(i);
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
                                sum += exp(likelihoodPerNeuron(i);
                            }

                            Hadj = log(sum);

                            likelihoodThreshold = likelihoodNoSpike - maximumLikelihoodPerNeuron - Hadj;

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

                            circbuffer[i]->dequeue(); // After everything is done, remove most past sample
                        }

                    }
                }

            }
        }
        masterSampleIndex += sampleIndex;

        stopTime = Time::getMillisecondCounterHiRes();
        log1->writeToLog("Duration taken was " + String(stopTime - startTime) + "ms for nSamples = " + String(nSamples));

    }


}
