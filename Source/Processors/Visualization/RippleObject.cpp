/*
  ==============================================================================

    RippleObject.cpp
    Created: 5 Jun 2013 2:45:06pm
    Author:  sethi-san

  ==============================================================================
*/

#include "RippleObject.h"
#include "SpikeObject.h"

#include "memory.h"
#include <stdlib.h>
#include "time.h"

// Simple method for serializing a RippleObject into a string of bytes
int packRipple(RippleObject* s, uint8_t* buffer, int bufferSize)
{

    //int reqBytes = 1 + 4 + 2 + 2 + 2 + 2 * s->nChannels * s->nSamples + 2 * s->nChannels * 2;
    //std::cout<<"reached to packripple beginning"<<std::endl;
    int idx = 0;

    s->eventType = RIPPLE_EVENT_CODE;

    //std::cout<<"the buffer address is "<<std::endl;
   // std::cout<<buffer<<std::endl;
    //std::cout<<"reached to just before memwriting"<<std::endl;
    memcpy(buffer+idx, &(s->eventType), 1);
    idx += 1;
    memcpy(buffer+idx, &(s->eventType), 1);
    idx += 1;
    //std::cout<<"reached to just before packing timestamp "<<s->timestamp<<std::endl;
    memcpy(buffer+idx, &(s->eventId), 1);
    idx += 1;
    //std::cout<<"reached to just before source"<<std::endl;
    memcpy(buffer+idx, &(s->eventChannel), 1);
    idx +=1;


    if (idx >= MAX_RIPPLE_BUFFER_LEN)
    {
        std::cout<<"Spike is larger than it should be. Size was:"<<idx<<" Max Size is:"<<MAX_RIPPLE_BUFFER_LEN<<std::endl;

    }
  //  makeBufferValid(buffer, idx);

    return idx;

}

// Simple method for deserializing a string of bytes into a Spike object
bool unpackRipple(RippleObject* s, const uint8_t* buffer, int bufferSize)
{
    //if (!isBufferValid(buffer, bufferSize))
    //    return false;

    int idx = 0;
    //std::cout<<"reached to just before memwriting"<<(int)*(buffer+idx)<<std::endl;
    memcpy(&(s->eventType), buffer+idx, 1);
    idx += 1;

    memcpy(&(s->eventId), buffer+idx, 1);
    //std::cout<<"reached to just after timestamp"<<s->timestamp<<std::endl;
    idx += 1;

    memcpy(&(s->eventChannel), buffer+idx, 1);
    //std::cout<<"reached to just after source"<<(s->source)<<std::endl;
    idx += 1;

    return true;


}

// Checks the validity of the buffer, this should be run before unpacking and after packing the buffer
bool isBufferValid(uint8_t* buffer, int bufferSize)
{

    if (! CHECK_BUFFER_VALIDITY)
        return true;

    uint16_t runningSum = 0;
    uint16_t value = 0;

    int idx;

    for (idx = 0; idx < bufferSize-2; idx += 2)
    {
        memcpy(buffer + idx, &value, 2);
        runningSum += value;
    }

    uint16_t integrityCheck = 0;
    memcpy(buffer + idx, &integrityCheck, 2);

    std::cout<<integrityCheck<<" == "<< runningSum <<std::endl;

    return (integrityCheck == runningSum);
}
/*
void makeBufferValid(uint8_t* buffer, int bufferSize)
{
    if (! CHECK_BUFFER_VALIDITY)
        return;

    uint16_t runningSum = 0;
    uint16_t value = 0;

    int idx;

    for (idx = 0; idx < bufferSize-2; idx += 2)
    {
        memcpy(buffer + idx, &value, 2);
        runningSum += value;
    }

    memcpy(&runningSum, buffer + idx, 2);

}*/

/*void generateSimulatedRipple(RippleObject* s, uint64_t timestamp, int noise)
{
    //std::cout<<"generateSimulatedSpike()"<<std::endl;

    uint16_t trace[][32] =
    {
        {
            880,	900,	940,	1040,	1290,	1790,	2475,	2995, 	3110, 	2890,
            2505,	2090,	1720,	1410, 	1155,  	945,	775,	635,	520, 	420,
            340,	265,	205,	155,	115,	80,		50,		34,		10, 	34,  	50,		80
        },
        {
            1040,   1090,   1190,	1350,	1600,	1960,   2380,   2790,   3080,	3140,
            2910,	2430,	1810,	1180,	680,	380,	270,	320,	460, 	630,
            770,	870,	940,	970,	990,	1000,	1000,	1000,	1000,	1000,  1000,	1000
        },
        {
            1000,	1000,	1000,	1000,	1000,	1040,	1140,	1440,	2040,	2240,
            2400,	2340,	2280,	1880,	1640,	920,	520,	300,	140,	040,
            20,		20,		40,		100,	260,	500,	740,	900,	960,	1000,	1000,	1000
        }
    };


    // We don't want to shift the waveform but scale it, and we don't want to scale
    // the baseline, just the peak of the waveform
    float scale[32] =
    {
        1.0, 1.0, 1.0, 1.0, 1.1, 1.2, 1.3, 1.5, 1.7, 2.0, 2.1, 2.2, 2.1, 2.0, 1.7, 1.5,
        1.3, 1.2, 1.1, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0
    };

    uint16_t gain = 2000;

    s->eventType = SPIKE_EVENT_CODE;
    s->timestamp = timestamp;
    s->source = 0;
    s->nChannels = 4;
    s->nSamples = 32;
    int idx=0;

    int waveType = rand()%2; // Pick one of the three predefined waveshapes to generate
    int shift = 1000 + 32768;

    for (int i=0; i<4; i++)
    {
        s->gain[i] = gain;
        s->threshold[i] = 4000;
        double scaleExponent = (double)(rand()%26+2) / 10.0f;  // Scale the wave between 50% and 150%

        for (int j=0; j<32; j++)
        {

            int n = 0;
            if (noise>0)
            {
                n = rand() % noise - noise/2;
            }

           // s->data[idx] = (trace[waveType][j] + n)  * pow(double(scale[j]),scaleExponent) + shift;
            idx = idx+1;
        }
    }

}
void generateEmptyRipple(RippleObject* s, int nChannels)
{

    s->eventType = RIPPLE_EVENT_CODE;
    s->timestamp = 0;
    s->source = 0;
    s->nChannels = 4;
    s->nSamples = 32;

    int idx = 0;
    for (int i=0; i<4; i++)
    {
        s->gain[i] = 0;
        s->threshold[i] = 0;
        for (int j=0; j<32; j++)
        {
           // s->data[idx] = 0;
            idx = idx+1;
        }
    }
}

void printRipple(RippleObject* s)
{

    std::cout<< " SpikeObject:\n";
    std::cout<< "\tTimestamp:" << s->timestamp;
    std::cout<< "\tSource:" << s->source;
    std::cout<< "\tnChannels:" <<s->nChannels;
    std::cout<<"\tnSamples" << s->nSamples;
    std::cout<<"\n\t 8 Data Samples:";
    for (int i=0; i<8; i++)
        //std::cout<<s->data+i<<" ";
    std::cout<<std::endl;
}
*/
