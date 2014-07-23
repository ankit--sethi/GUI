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

Array<float> CircularQueue::show(int index)
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
    : GenericProcessor("Spike Sorter"), nullbuffer(2,100), dataBuffer(nullbuffer) //, threshold(200.0), state(true)

{

    // THIS CONSTRUCTOR TAKES CARE OF EVERYTHING EXCEPT ParameterEstimator* node and Array<CircularBuffer> circbuffer
    // THEY ARE HANDLED ONCE SpikeSorter::process() STARTS AND CALLS SpikeSorter::setParameter()

    alpha = 1e-1;
    kappa0 = 0.01;
    nu0 = 0.1;
    buffersArePlush = false;

    K = 3; // number of SVD components to use; need to make this a control
    //phi0 = new float*[K];

    Array<float> temp;
    for (int i = 0; i < K; i++)
    {
       temp.add(0);
    }
    for (int j = 0; j< K; j++)
    {
        phi0.add(temp);
    }
    Array<float>* t = phi0.getRawDataPointer();
    for (int i = 0; i < K; i++)
    {
        for (int j = 0; j< K; j++)
        {
            if (i == j)
                t->setUnchecked(j,1*0.1);
            else
                t->setUnchecked(j,0);
        }
        t++;
    }

    //Phi_0=.1*eye(K);
    number = 0;
    apii = 1;
    bpii = 1e7;
    beta = 1/(30*(40000));
    std::cout << " Sample RATE is " << getSampleRate() << "//";
    tau = 10;

    P = 40;

    checkIfAllParametersEstimated = false;
    paramsCopied = false;


    Cmax = 50;  //maximum possible number of neurons present
    curndx = 0; //used to index current location in buffer
    lookahead = 500; //functionally, it is the size of the buffer
    range = 40; // defines basically the width of a spike
    // setting sigma
    temp.clear();
    for (int i = 0; i < P; i++)
    {
        temp.add(0);
    }
    for (int j = 0; j< P; j++)
    {
        sigma.add(temp);
    }
    // end of setting sigma

    pii=apii/bpii;



    Array<Array<float>> phi0inv;
    invertSquareTwoDimMatrix(phi0, phi0inv);

    Array<Array<float>> phi0invnu0;
    MatrixMultiplyWithConstant(phi0inv, nu0, phi0invnu0);

    Array<Array<float>> phi0invnu0forR;
    MatrixMultiplyWithConstant(phi0invnu0,0.2, phi0invnu0forR);

    for (int i = 0; i < Cmax; i++)
    {
        expandedfloat zerovar;
        zerovar.setNonInf();
        zerovar.setVar(0.0);
        nu.add(nu0);
        phi.add(phi0);
        kappa.add(kappa0);
        lamclus.add(phi0invnu0);
        lamclusinv.add(phi0invnu0);
        R.add(phi0invnu0forR);
        Rinv.add(phi0invnu0forR);
        ltheta.add(zerovar);
        ngamma.add(0);
        tlastspike.add(0);
        likelihoodPerNeuron.add(0);
    }
    temp.clear();


    for (int j = 0; j < Cmax; j++)
    {
            temp.add(0); // careful: from this point on, Phi0 is actually 0.1*Phi0
    }
    for (int i = 0; i < K; i++)
    {
        muu.add(temp);
        muu0.add(temp);
    }

    temp.clear();
    for (int j = 0; j < P; j++)
    {
            temp.add(0);
    }
    for (int i = 0; i < P; i++)
    {
        Q.add(temp);
        Qinv.add(temp);
    }
    temp.clear();
    for (int j = 0; j< K; j++)
    {
            temp.add(0);
    }
    for (int i = 0; i < K; i++)
    {
        Qhat.add(temp);
        Qmat.add(temp);
    }
    temp.clear();
    for (int j = 0; j< 1; j++)
    {
            temp.add(0);
    }
    for (int i = 0; i < K; i++)
    {
        mhat.add(temp);
        yhat.add(temp);
    }
    temp.clear();
    for (int i = 0; i < K; i++)
    {
        temp.add(0);
    }
    for (int j = 0; j < P; j++)
    {
        ReducedDictionary.add(temp);
    }
    MatrixTranspose(ReducedDictionary,ReducedDictionaryTranspose);

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

    //muuS possibly does not need initialization. More on this as the story unfolds.
    //lamclusS possibly does not need initialization. More on this as the story unfolds.


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

void SpikeSorter::invertSquareTwoDimMatrix(const Array<Array<float>>& A, Array<Array<long double>>& result )//float* calcDeterminant)
{

    result.clear();

    int size = A.size();

    Array<long double> temp;
    for (int i = 0; i < A[0].size(); i++)
    {
        temp.add(0);
    }
    for (int i = 0; i < A.size(); i++)
    {
        result.add(temp);
    }

    long double** matrixForInversion = new long double*[size];
    long double** matrixForInverse = new long double*[size];

    for (int i = 0; i < A.size(); i++)
    {
        matrixForInversion[i] = new long double[size];
        matrixForInverse[i] = new long double[size];
        for (int j = 0; j < A[0].size(); j++)
        {
            matrixForInversion[i][j] = (long double)(A[i][j]);
            matrixForInverse[i][j] = 0;
        }
    }

    MatrixInversion(matrixForInversion,size,matrixForInverse);

    Array<long double>* test = result.getRawDataPointer();

    for (int i = 0; i < A.size(); i++)
    {

        for (int j = 0; j < A[0].size(); j++)
        {
            test->setUnchecked(j,matrixForInverse[i][j]);
        }
        test++;
    }
    for (int i = 0; i < A.size(); i++)
    {
        delete [] matrixForInversion[i];
        delete [] matrixForInverse[i];
    }
     delete [] matrixForInversion;
     delete [] matrixForInverse;
}

void SpikeSorter::invertSquareTwoDimMatrix(const Array<Array<float>>& A, Array<Array<float>>& result)
{

    result.clear();

    int size = A.size();

    Array<float> temp;
    for (int i = 0; i < A[0].size(); i++)
    {
        temp.add(0);
    }
    for (int i = 0; i < A.size(); i++)
    {
        result.add(temp);
    }

    double** matrixForInversion = new double*[size];
    double** matrixForInverse = new double*[size];

    for (int i = 0; i < A.size(); i++)
    {
        matrixForInversion[i] = new double[size];
        matrixForInverse[i] = new double[size];
        for (int j = 0; j < A[0].size(); j++)
        {
            matrixForInversion[i][j] = double(A[i][j]);
            matrixForInverse[i][j] = 0;
        }
    }

    MatrixInversion(matrixForInversion,size,matrixForInverse);

    Array<float>* test = result.getRawDataPointer();

    for (int i = 0; i < A.size(); i++)
    {

        for (int j = 0; j < A[0].size(); j++)
        {
            test->setUnchecked(j,float(matrixForInverse[i][j]));
        }
        test++;
    }

    for (int i = 0; i < A.size(); i++)
    {
        delete [] matrixForInversion[i];
        delete [] matrixForInverse[i];
    }
     delete [] matrixForInversion;
     delete [] matrixForInverse;
}

double SpikeSorter::getMatrixDeterminant(const Array<Array<double> > &A)
{
    int size = A.size();
    double** matrixForInversion = new double*[size];

    for (int i = 0; i < A.size(); i++)
    {
        matrixForInversion[i] = new double[size];
        for (int j = 0; j < A[0].size(); j++)
        {
            matrixForInversion[i][j] = A[i][j];
        }
    }

    double calcDeterminant = CalcDeterminant(matrixForInversion, size);
    for (int i = 0; i < A.size(); i++)
    {
        delete [] matrixForInversion[i];
    }
    delete [] matrixForInversion;
    return calcDeterminant;
}
double SpikeSorter::getMatrixDeterminant(const Array<Array<float>>& A)
{
    int size = A.size();
    float** matrixForInversion = new float*[size];

    for (int i = 0; i < A.size(); i++)
    {
        matrixForInversion[i] = new float[size];
        for (int j = 0; j < A[0].size(); j++)
        {
            matrixForInversion[i][j] = A[i][j];
        }
    }

    double calcDeterminant = CalcDeterminant(matrixForInversion, size);
    for (int i = 0; i < A.size(); i++)
    {
        delete [] matrixForInversion[i];
    }
    delete [] matrixForInversion;
    return calcDeterminant;
}

double SpikeSorter::getMatrixDeterminant(const Array<Array<long double>>& A)
{
    int size = A.size();
    long double** matrixForInversion = new long double*[size];

    for (int i = 0; i < A.size(); i++)
    {
        matrixForInversion[i] = new long double[size];
        for (int j = 0; j < A[0].size(); j++)
        {
            matrixForInversion[i][j] = A[i][j];
        }
    }

    double calcDeterminant = CalcDeterminant(matrixForInversion, size);
    for (int i = 0; i < A.size(); i++)
    {
        delete [] matrixForInversion[i];
    }
    delete [] matrixForInversion;
    return calcDeterminant;
}

double SpikeSorter::getMatrixLogDeterminant(const Array<Array<float>>& A)
{
    Array<Array<float>> chA;
    cholesky(A, chA);
    Array<float> Diag;
    returnMainDiagonal(chA, Diag);
    float sum = 0;

    for (int i = 0; i < Diag.size(); i++)
    {
        sum += log(Diag[i]);
    }
    return 2*sum;
}

void SpikeSorter::setParameter(int parameterIndex, float newValue)
{

    if (parameterIndex == 1 ) //
    {

        for (int i = 0; i < ltheta.size(); i++)
        {
            if (!checkIfLogIsInf(ngamma[i],(alpha+totalSpikesFound)))
            {
                ltheta[i].setVar(log(ngamma[i]/(alpha+totalSpikesFound)));

            }
            else
            {
                ltheta[i].setInf();
            }

        }
        ltheta[neuronCount].setVar(log(alpha/(alpha+totalSpikesFound)));

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

float SpikeSorter::getaBaInnerProductSum(Array<float>& xwindm, Array<Array<float>>& InnerMatrix)
{
    //Array<float> xwind;
    //xwind = circbuffer[electrodeIndex].show(P);
    //std::cout << "Inner Matrix size = " <<InnerMatrix[0].size() << " and " << InnerMatrix.size() << " and xwind size is " << xwind.size() << std::endl;

    float sum = 0, totalsum = 0;
    if (InnerMatrix[0].size() == xwind.size()) // This should ALWAYS be true
    {
        for (int i = 0; i < InnerMatrix.size(); i++)
        {
            sum = 0;
            for (int j = 0; j < InnerMatrix[0].size(); j++)
            {
                sum = sum + InnerMatrix[i][j]*xwind[j];
                //std::cout << InnerMatrix[i][j] << "is the innerm and " << xwind[j] << std::endl;
            }

            sum = sum*xwind[i];
            totalsum = totalsum + sum;
        }
        return totalsum;
    }
    else
    {
        std::cout<<"aBa product failed"<<std::endl;
        return -1;
    }

}

void SpikeSorter::MatrixMultiplication (const Array<Array<float>>& A, const Array<Array<float>>& B, Array<Array<float>>& result)
{
    result.clear();
    float temp;  Array<float> temp1;
     for (unsigned int i = 0; i < B[0].size(); i++)
    {
       temp1.add(0);
    }
    for (unsigned int j = 0; j < A.size(); j++)
    {
       result.add(temp1);
    }
    Array <float>* t = result.getRawDataPointer();

    for (unsigned int i = 0; i < A.size(); i++)
    {
        for (unsigned int j = 0; j < B[0].size(); j++)
        {
            temp = A[i][0]*B[0][j];
            for (unsigned int k = 1; k < B.size(); k++)
            {
                temp += A[i][k]*B[k][j];
            }
            t->setUnchecked(j,temp);
        }
        t++;
    }
}



void SpikeSorter::MatrixTranspose(const Array<Array<float>>& A, Array<Array<float>>& result)
{
    result.clear();
    Array<float> temp;
    for (int i = 0; i < A[0].size(); i++)
    {
        temp.clear();
        for(int j = 0; j < A.size(); j++)
        {
            temp.add(A[j][i]);
        }
        result.add(temp);
    }
}

void SpikeSorter::MatrixGetRowOrColumn(const Array<Array<float>>& A, bool row, int index, Array<Array<float>>& result)
{
    Array<float> temp;
    result.clear();
    if (row)
    {
        for (int i = 0; i < A[0].size(); i++)
        {
            temp.add(A[index][i]);
        }
        result.add(temp);
    }
    else
    {
        for (int i = 0; i < A.size(); i++)
        {
            temp.clear();
            temp.add(A[i][index]);
            result.add(temp);
        }
    }

}


//Cholesky decomposition of matrix A
void SpikeSorter::cholesky(const Array<Array<float>>& A, Array<Array<float>>& result)
{
     Array <float> tmp;
    result.clear();

    for (int j = 0; j<A[0].size(); j++)
    {
        tmp.add(0);
    }
    for (int i = 0; i < A.size(); i++)
    {

        result.add(tmp);
    }
    int n = A.size();

    Array<float>* temp = result.getRawDataPointer();
     for (int i = 0; i < n; i++)
        for (int j = 0; j < (i+1); j++)
        {
            double s = 0;
            for (int k = 0; k < j; k++)
                s += result[i][k] * result[j][k];

            (temp+i)->setUnchecked(j,(i == j) ?
                        sqrt(A[i][i] - s) :
                (1.0 / result[j][j] * (A[i][j] - s)));
        }

}

void SpikeSorter::AddMatrices(const Array<Array<float>>& A, const Array<Array<float>>& B, Array<Array<float>>& result)
{
    result.clear();
    Array<float> flimsy;
    for (int i = 0; i < A[0].size(); i++)
    {
        flimsy.add(0);
    }
    for (int j = 0; j < A.size(); j++)
    {
        result.add(flimsy);
    }
    Array<float>* temp = result.getRawDataPointer();

    for (int j = 0; j < A.size(); j++)
    {
        for (int k = 0; k < A[0].size(); k++)
        {
            (temp+j)->setUnchecked(k, A[j][k] + B[j][k]);
        }
    }
}

void SpikeSorter::AddMatrixToItself(Array<Array<float>>& A, const Array<Array<float>>& B) // Adds B to A and saves in A
{
    Array<float>* temp = A.getRawDataPointer();

    for (int j = 0; j < A.size(); j++)
    {
        for (int k = 0; k < A[0].size(); k++)
        {
            (temp+j)->setUnchecked(k, A[j][k] + B[j][k]);
        }
    }
}

void SpikeSorter::getLikelihoodPerNeuron(int FoundNeuronIndex, float tmp)
{
    //std::cout << "The inner product here is " << getaBaInnerProductSum(xwind, Qinv);
    if (suppress)
    likelihoodPerNeuron.set(FoundNeuronIndex, -(0.5*P)*log(2*PI) - tmp - suppresslikelihood - 0.5*getaBaInnerProductSum(xwind, Qinv));
    else
    likelihoodPerNeuron.set(FoundNeuronIndex, -(0.5*P)*log(2*PI) - tmp - 0.5*getaBaInnerProductSum(xwind, Qinv));
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
    lthr.clear();
    lon.clear();
    xwindLonger.clear();

    currentSpike.data[currentIndex++] = uint16(trigSample/ channels[chan]->bitVolts + 32768);

    xwindLonger = xwind;

}

void SpikeSorter::addNewSampleAndLikelihoodsToCurrentSpikeObject(float sample, MidiBuffer& eventBuffer, int chan)
{
    if (currentIndex < range)
    {

        currentSpike.data[currentIndex++] = uint16(sample/ channels[chan]->bitVolts + 32768);
        lthr.add(likelihoodThreshold);
        Array<float> temp;
        for (int i = 0; i<likelihoodPerNeuron.size(); i++)
        {
            temp.add(likelihoodPerNeuron[i]);
        }
        lon.add(temp);
        xwindLonger.add(sample);

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
void SpikeSorter::MatrixMultiplyWithConstant (const Array<Array<float>>& A, float c, Array<Array<float>>& result)
{
    result.clear();
    Array<float> flimsy;
    for (int i = 0; i < A[0].size(); i++)
    {
        flimsy.add(0);
    }
    for (int j = 0; j < A.size(); j++)
    {
        result.add(flimsy);
    }

    Array<float>* temp = result.getRawDataPointer();

    for (int i = 0; i < A.size(); i++)
    {
        for (int j = 0; j < A[0].size(); j++)
        {
            (temp+i)->setUnchecked(j, A[i][j]*c);
        }
    }
}

void SpikeSorter::MatrixReplaceColumn (Array<Array<float>>& A, const Array<Array<float>>& B, int colindex)
{

    Array<float>* temp = A.getRawDataPointer();
    //B here is a column vector, i.e., an Nx1 matrix
    if (B[0].size() > 1)
    {
        std::cout<< "Error in MatrixReplaceColumn: B should be an Nx1 matrix"<<std::endl;
    }
    for (int i = 0; i < A.size(); i++)
    {
        (temp+i)->setUnchecked(colindex, B[i][0]);
    }
}

void SpikeSorter::eye (int dim, Array<Array<float>>& result)
{
    Array<float> temp;
    result.clear();
    for (int i = 0; i < dim; i++)
    {
        temp.add(0);
    }
    for (int i = 0; i < dim; i++)
    {
        result.add(temp);
    }
    Array<float>* t = result.getRawDataPointer();

    for (int i = 0; i < dim; i++)
    {
        for (int j = 0; j < dim; j++)
        {
            if (i == j)
                t->setUnchecked(j, 1);
        }
        t++;
    }

}

void SpikeSorter::dotSlashOperation(const Array<Array<float>>& A, float c, Array<Array<float>>& result)
{
    result.clear();
    Array<float> flimsy;
    for (int i = 0; i < A[0].size(); i++)
    {
        flimsy.add(0);
    }
    for (int j = 0; j < A.size(); j++)
    {
        result.add(flimsy);
    }
    Array<float>* temp = result.getRawDataPointer();
                    for (int i = 0; i < A.size(); i++)
{
    for (int j = 0; j < A[0].size(); j++)
    {
        temp->setUnchecked(j, c/A[i][j]);
    }
    temp++;
}

}

void SpikeSorter::returnMainDiagonal(const Array<Array<float>> &A, Array<float> &result)
{
    result.clear();
    for (int i = 0; i< A.size(); i ++)
    {
        for (int j = 0; j < A[0].size(); j++)
        {
            if (i == j)
            {
                result.add(A[i][j]);
            }
        }
    }

}

void SpikeSorter::updateAllSortingParameters()
{
    int neuronID = currentSpike.neuronID;
    Array<Array<float>> temp, temp2;

    MatrixMultiplication(ReducedDictionaryTranspose,lambda, temp);
    MatrixMultiplication(temp,ReducedDictionary, Qmat);

    AddMatrices(Qmat, lamclus[neuronID], temp2);
    Qmat = temp2;

    Array<Array<float>> Qmatinv;
    invertSquareTwoDimMatrix(Qmat, Qmatinv);

    Array<float> den;
    Array<Array<float>> x;
    for (int i = 0; i < P; i++)
    {
        den.clear();
        den.add(xwindLonger[idx + i]);
        x.add(den);
    }



    Array<Array<float>> t;
    MatrixGetRowOrColumn(muu, false, neuronID, t);    
    Array<Array<float>> t1;
    MatrixMultiplication(lamclus[neuronID],t, t1);
    Array<Array<float>> t2;
    MatrixMultiplication(temp,x,t2);
    Array<Array<float>> t3;
    AddMatrices(t1,t2,t3);

    MatrixMultiplication(Qmatinv,t3,yhat);


    ngamma.set(neuronID, ngamma[neuronID] + 1);
    deltaT = masterSampleIndex + sampleIndex -P - range + idx - tlastspike[neuronID];
    tlastspike.set(neuronID, masterSampleIndex + sampleIndex - P - range + idx);
    double ebet = exp(-1*double(beta)*double(deltaT));
    //std::cout << "Ebet and delta T is" << ebet << "//" << deltaT << "//";

    //MatrixGetRowOrColumn(muu,false, neuronID,t);

    MatrixMultiplyWithConstant(t, (ebet), t1);

    MatrixGetRowOrColumn(muu0,false, neuronID,t);

    MatrixMultiplyWithConstant(t, (1-ebet),t3);

    AddMatrices(t3 ,t1, mhat);




    //MatrixGetRowOrColumn(muu0,false,neuronID,t);
    //std::cout << "KAPPA is" << kappa[neuronID] << "//";
    MatrixMultiplyWithConstant(t,kappa[neuronID],t1);
    AddMatrices(t1,yhat,t2);
    MatrixMultiplyWithConstant(t2,1/(kappa[neuronID]+1),t3);
    MatrixReplaceColumn(muu0,t3,neuronID);

    invertSquareTwoDimMatrix(R[neuronID],t);



    MatrixMultiplyWithConstant(t,pow(ebet,2),t1);
    Array<Array<float>> eyeK;
    eye(K,eyeK);
    MatrixMultiplyWithConstant(eyeK,(1/tau)*(1-pow(ebet,2)),t2);
    AddMatrices(t2,t1,Qhat);

    invertSquareTwoDimMatrix(Qhat,t);
    Array<Array<float>> Qhatinv;
    Qhatinv = t;

    t2 = R[neuronID];
    AddMatrices(Qhatinv, lamclus[neuronID], t2);
    R.set(neuronID,t2);
    /*std::ofstream myfile3 ("R[neuronID].txt");

    for (int i = 0; i < R[neuronID].size(); i++)
    {
        for (int j = 0; j < R[neuronID][0].size(); j++)
        {
            myfile3 << R[neuronID][i][j];
            if ( j != R[neuronID][0].size() -1 )
            {
                myfile3 << ",";
            }
        }
        myfile3 << std::endl;
    }
    myfile3.close();*/

    MatrixMultiplication(lamclus[neuronID],yhat,t);

    //t1 = Qhatinv;
    MatrixMultiplication(Qhatinv,mhat, t2);
    AddMatrices(t2,t, t3);
    Array<Array<float>> t4;
    invertSquareTwoDimMatrix(R[neuronID],t4);

    Array<Array<float>> t5;
    MatrixMultiplication(t4,t3,t5);
    MatrixReplaceColumn(muu,t5,neuronID);
    //t = Qmatinv;
    //MatrixMultiplication(ReducedDictionaryTranspose,lambda, t1);
    MatrixMultiplication(temp,x, t2);
    std::cout << "The Neuron ID is //" << neuronID << "//";
    MatrixGetRowOrColumn(muu,false,neuronID, t3);
    MatrixMultiplication(lamclus[neuronID],t3, t4);
    AddMatrices(t2,t4,t5);
    MatrixMultiplication(Qmatinv,t5,yhat);

    //std::cout << "The projections of the spike on the dictionary are: " << std::endl;
    /*for (int i = 0; i < yhat.size(); i++)
    {
        std::cout << yhat[i][0] << "//";
    }
    std::cout << std::endl;*/
    // t is already Qmatinv
    t1 = phi[neuronID];
    AddMatrices(phi[neuronID], Qmatinv, t1); // need to replace this with AddMatrixWithItself
    phi.set(neuronID, t1);
    temp.clear();
    MatrixMultiplyWithConstant(yhat,-1, t);
    MatrixGetRowOrColumn(muu, false, neuronID, t1);
    AddMatrices(t, t1, temp);
    Array<Array<float>> temptr;
    MatrixTranspose(temp, temptr);
    MatrixMultiplication(temp, temptr, t);
    MatrixMultiplyWithConstant(t,(kappa[neuronID]+1), t1);
    dotSlashOperation(t1,kappa[neuronID], t2);
    Array<Array<float>> tempt = phi[neuronID];
    AddMatrices(tempt,t2, t3);
    phi.set(neuronID, t3);
    kappa.set(neuronID, kappa[neuronID] + 1);
    nu.set(neuronID, nu[neuronID] + 1);
    tempt = phi[neuronID];
    invertSquareTwoDimMatrix(tempt, t);
    //t1 = lamclus[neuronID];
    MatrixMultiplyWithConstant(t,nu[neuronID], t1);
    lamclus.set(neuronID, t1);
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

        if(!paramsCopied)
        {
            std::cout<<"Entered Params copied if condition" << std::endl;

            Array<float>* test = sigma.getRawDataPointer();

            for (int i = 0; i < P; i++)
            {
                float power = -1*i;
                for (int j = 0; j < P; j++)
                {

                    (test+i)->setUnchecked(j,node->getElectrodeNoiseVar(0)*pow(node->acfLag1, abs(power)));

                    if ( i == j )
                    (test+i)->setUnchecked(j,node->getElectrodeNoiseVar(0));

                    power++;
                }
            }


            std::cout<< " ACF LAG in Spike Sorter is " << node->acfLag1 << " and the Var is " << node->getElectrodeNoiseVar(0) << "  "  << std::endl;
            Array<Array<long double>> proxy;
            invertSquareTwoDimMatrix(sigma, proxy);

            Array<float> tempz;
            for (int i = 0; i < sigma[0].size(); i++)
            {
                tempz.add(0);
            }
            for (int i = 0; i < sigma.size(); i++)
            {
                lambda.add(tempz);
            }
            Array<float>* ptr = lambda.getRawDataPointer();

            for (int i = 0; i < proxy.size(); i++)
            {
                for(int j = 0; j < proxy[0].size(); j++)
                {

                    (ptr+i)->setUnchecked(j,float(proxy[i][j]));
                }

            }
            logDeterminantOfLambda = (getMatrixLogDeterminant(lambda));

            std::cout<<"reached line 3 and det is" << logDeterminantOfLambda << ReducedDictionary.size() << "x" << ReducedDictionary[0].size() << "//" << std::endl;

            Array<float>* ptr1 = ReducedDictionary.getRawDataPointer();
            for (int i = 0; i < P; i++)
            {
                for (int j = 0; j < K; j++)
                {
                    (ptr1+i)->setUnchecked(j,node->dictionary.get(i,j));
                }
            }
            MatrixTranspose(ReducedDictionary, ReducedDictionaryTranspose);
            paramsCopied = true;
            std::cout<<"Params copied and nSamples is" << nSamples << std::endl;
        }

        Array<Array<float>> t;
        Array<Array<float>> t1;
        Array<Array<float>> t2;
        Array<Array<float>> t3;

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
                        number++;
                        int currentChannel = *(node->electrodes[i]->channels+chan);

                        float sample = getNextSample(currentChannel);
                        sampleIndex++;

                        circbuffer[i]->enqueue(sample);

                        pii = (apii + totalSpikesFound)/(bpii + masterSampleIndex);
                        //std::cout << "Mastersampleindex is //" << masterSampleIndex << "//";
                        threshold = log(pii/(1-pii));


                        if(circbuffer[i]->isBufferPlush(P))
                        {

                            // find likelihood of no spike
                            xwind = circbuffer[i]->show(P);
                            xwindloop = xwind;
                            //std::cout <<" XWIND - //" << xwind.size() << "  " << lambda.size() << "   " << lambda[0].size() << std::endl;
                            //std::cout<<"WOOOOOOO" << logDeterminantOfLambda << " XYZ " << getaBaInnerProductSum(xwind, lambda) << std::endl;



                            for (int i = 0; i < ltheta.size(); i++)
                            {
                                if (!checkIfLogIsInf(ngamma[i],(alpha+totalSpikesFound)))
                                {
                                    ltheta[i].setVar(log(ngamma[i]/(alpha+totalSpikesFound)));

                                }
                                else
                                {
                                    ltheta[i].setInf();
                                }

                            }
                            ltheta[neuronCount].setVar(log(alpha/(alpha+totalSpikesFound)));


                            likelihoodNoSpike = (-1*P/2)*log(2*PI) + 0.5*logDeterminantOfLambda - 0.5*getaBaInnerProductSum(xwind, lambda);
                            for (int j = 0; j < neuronCount; j++)
                            {
                                //std::cout << "inverse 11" << std::endl;
                                invertSquareTwoDimMatrix(lamclus[j], t);
                                lamclusinv.setUnchecked(j,t);
                                //std::cout << "inverse 12" << std::endl;
                                invertSquareTwoDimMatrix(R[j], t);
                                Rinv.setUnchecked(j,t);

                                AddMatrices(lamclusinv[j],Rinv[j], Qt);
                                //ReducedDictionaryTranspose =  MatrixTranspose(ReducedDictionary);
                                MatrixMultiplication(Qt,ReducedDictionaryTranspose, t1);
                                MatrixMultiplication(ReducedDictionary,t1, t2);
                                //std::cout<<"reached addition";
                                AddMatrices(sigma,t2, Q);

                                Array<Array<long double>> Qinvld;
                                invertSquareTwoDimMatrix(Q, Qinvld);
                                Array<float>* ptr = Qinv.getRawDataPointer();

                                for (int i = 0; i < Qinvld.size(); i++)
                                {
                                    for(int j = 0; j < Qinvld[0].size(); j++)
                                    {

                                        (ptr+i)->setUnchecked(j,float(Qinvld[i][j]));
                                    }

                                }
                                //take projection of one neuron and then subtract it from
                                MatrixGetRowOrColumn(muu, false, j, t3);
                                Array<Array<float>> proj;
                                MatrixMultiplication(ReducedDictionary, t3, proj);


                                //
                                for (int k = 0; k < proj.size(); k++)
                                {
                                    xwindloop.set(k, xwind[k] - proj[k][i]);
                                }

                                //subtracted
                                //std::cout<<"reached line 1280" << std::endl;
                                // now take Cholesky factorization of Q
                                Array<Array<float>> chQ;
                                cholesky(Q, chQ);
                                Array<float> Diag;
                                returnMainDiagonal(chQ, Diag);
                                float sum = 0;

                                for (int i = 0; i < Diag.size(); i++)
                                {
                                    sum += log(Diag[i]);
                                }

                                if( ( (masterSampleIndex + sampleIndex) - tlastspike[j]) < 40000*50/10000 )  // THIS NEEDS TO BE INVESTIGATED
                                    suppress = true;
                                else
                                    suppress = false;

                                getLikelihoodPerNeuron(j, sum);
                            }
                            //std::cout<<"reached line 1300" << std::endl;
                            //std::cout << "inverse 13" << std::endl;

                            invertSquareTwoDimMatrix(lamclus[neuronCount], t);
                            lamclusinv.setUnchecked(neuronCount, t);
                            //std::cout << "inverse 14" << std::endl;
                            invertSquareTwoDimMatrix(R[neuronCount], t1);
                            //std::cout << "after inverse 14" << std::endl;
                            Rinv.setUnchecked(neuronCount, t1);

                            AddMatrices(lamclusinv[neuronCount] , Rinv[neuronCount], Qt);
                            //std::cout<<"reached line 1308" << std::endl;
                            MatrixTranspose(ReducedDictionary, ReducedDictionaryTranspose);
                            Array<Array<float>> temp1;
                            MatrixMultiplication(Qt, ReducedDictionaryTranspose, temp1);

                            MatrixMultiplication(ReducedDictionary,temp1, Q);


                            AddMatrixToItself(Q,sigma);
                            //std::cout << "the size of Q is" << Q.size() << " and  "  << Q[0].size() <<std::endl;
                            test = true;
                            Array<Array<long double>> Qinvld;
                            invertSquareTwoDimMatrix(Q, Qinvld);
                            Array<float>* ptr = Qinv.getRawDataPointer();

                            for (int i = 0; i < Qinvld.size(); i++)
                            {
                                for(int j = 0; j < Qinvld[0].size(); j++)
                                {

                                    (ptr+i)->setUnchecked(j,float(Qinvld[i][j]));
                                }

                            }
                            if (!test)
                            {
                            std::cout<<"boo";
                            }

                            Array<Array<float>> t2;
                            cholesky(Q, t2);

                            Array<float> Diag;
                            returnMainDiagonal(t2, Diag);

                            float sum = 0;                            
                            for (int i = 0; i < Diag.size(); i++)
                            {
                                sum += log(Diag[i]);
                            }

                            suppress = false;
                            getLikelihoodPerNeuron(neuronCount,sum);
                            for (int i = 0; i <= neuronCount; i++ )
                            {
                                likelihoodPerNeuron.set(i, likelihoodPerNeuron[i] + ltheta[i].getVar());
                                //std::cout << "And " << likelihoodPerNeuron[i] << " is the Likelihood per Neuron for " << i << " and " << " ltheta is " << ltheta[i].getVar() << std::endl;
                            }
                            for (int i = 0; i <= neuronCount; i++ )
                            {
                                if (i == 0)
                                    maximumLikelihoodPerNeuron = likelihoodPerNeuron[i];
                                else if(maximumLikelihoodPerNeuron < likelihoodPerNeuron[i])
                                    maximumLikelihoodPerNeuron = likelihoodPerNeuron[i];
                            }
                            for (int i = 0; i <= neuronCount; i++ )
                            {
                                likelihoodPerNeuron.set(i, likelihoodPerNeuron[i] - maximumLikelihoodPerNeuron);
                            }
                            //std::cout<<"reached line 1342" << std::endl;
                            sum = 0;
                            for (int i = 0; i <= neuronCount; i++ )
                            {
                                sum += exp(likelihoodPerNeuron[i]);
                            }

                            Hadj = log(sum);

                            likelihoodThreshold = likelihoodNoSpike - maximumLikelihoodPerNeuron - Hadj;

                            //std::cout<< "Likelihood Threshold is = " << likelihoodThreshold << " and the threshold is " << threshold << std::endl;

                            std::ofstream myfile;
                            myfile.open ("test.txt", std::ofstream::out | std::ofstream::app);
                            myfile << likelihoodThreshold << "," << number << std::endl;
                            myfile.close();

                            if (likelihoodThreshold < 2)
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
        //std::cout << "UPDATED MASTER SAMPLE INDEX";
        masterSampleIndex += sampleIndex;
        //std::cout << "UPDATED MASTER SAMPLE INDEX";

    }

}
