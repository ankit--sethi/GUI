/*
  ==============================================================================

    RippleDisplayCanvas.cpp
    Created: 20 Jun 2013 3:32:29pm
    Author:  sethi-san

  ==============================================================================
*/
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
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

*/

#include "RippleDisplayCanvas.h"

RippleDisplayCanvas::RippleDisplayCanvas(RippleDisplayNode* n) :
    processor(n), newRipple(false)
{

    RippleBuffer = processor->getRippleBufferAddress();

    viewport = new Viewport();
    rippleDisplay = new RippleDisplay(this, viewport);

    viewport->setViewedComponent(rippleDisplay, false);
    viewport->setScrollBarsShown(true, false);

    scrollBarThickness = viewport->getScrollBarThickness();

    clearButton = new UtilityButton("Clear plots", Font("Small Text", 13, Font::plain));
    clearButton->setRadius(3.0f);
    clearButton->addListener(this);
    addAndMakeVisible(clearButton);

    addAndMakeVisible(viewport);

    setWantsKeyboardFocus(true);

    update();

}

RippleDisplayCanvas::~RippleDisplayCanvas()
{

}

void RippleDisplayCanvas::beginAnimation()
{
    std::cout << "RippleDisplayCanvas beginning animation." << std::endl;

    startCallbacks();
}

void RippleDisplayCanvas::endAnimation()
{
    std::cout << "RippleDisplayCanvas ending animation." << std::endl;

    stopCallbacks();
}

void RippleDisplayCanvas::update()
{

    std::cout << "Updating RippleDisplayCanvas" << std::endl;

    int nPlots = processor->getNumElectrodes();
    rippleDisplay->removePlots();

    for (int i = 0; i < nPlots; i++)
    {
        rippleDisplay->addRipplePlot(processor->getNumberOfChannelsForElectrode(i), i);
    }

    rippleDisplay->resized();
    rippleDisplay->repaint();
}


void RippleDisplayCanvas::refreshState()
{
    // called when the component's tab becomes visible again
    resized();
}

void RippleDisplayCanvas::resized()
{
    viewport->setBounds(0,0,getWidth(),getHeight()-90);

    rippleDisplay->setBounds(0,0,getWidth()-scrollBarThickness, rippleDisplay->getTotalHeight());

    clearButton->setBounds(10, getHeight()-40, 100,20);

}

void RippleDisplayCanvas::paint(Graphics& g)
{

    g.fillAll(Colours::darkgrey);

}

void RippleDisplayCanvas::refresh()
{
    processRippleEvents();

    repaint();
}

void RippleDisplayCanvas::processRippleEvents()
{
   // std::cout<<"Entered processRipple Events"<<std::endl;
  //  std::cout<<RippleBuffer->getNumEvents()<<std::endl;

   /* if (RippleBuffer->getNumEvents() > 0)
    {

        MidiBuffer::Iterator i(*RippleBuffer);
        MidiMessage message(0xf4);

        int samplePosition = -40;

        i.setNextSamplePosition(samplePosition);
        std::cout<<i.getNextEvent(message, samplePosition)<<std::endl;
        while (i.getNextEvent(message, samplePosition))
        {

            const uint8_t* dataptr = message.getRawData();
            int bufferSize = message.getRawDataSize();

            RippleObject newRipple;
          //  std::cout<<"befpre unpackRipple "<<"bufferSize is :"<<bufferSize<<std::endl;
            bool isValid = unpackRipple(&newRipple, dataptr, bufferSize);
        //    std::cout<<"after unpackRipple"<<std::endl;
            int electrodeNum = newRipple.source;

            if (isValid)
                rippleDisplay->plotRipple(newRipple, electrodeNum);

        }

    }

    RippleBuffer->clear();
  //  std::cout<<"Exited processRipple Events"<<std::endl;*/

}

bool RippleDisplayCanvas::keyPressed(const KeyPress& key)
{

    KeyPress c = KeyPress::createFromDescription("c");

    if (key.isKeyCode(c.getKeyCode())) // C
    {
        rippleDisplay->clear();
        
        std::cout << "Clearing display" << std::endl;
        return true;
    }

    return false;

}

void RippleDisplayCanvas::buttonClicked(Button* button)
{

    if (button == clearButton)
    {
        rippleDisplay->clear();
    }
}



// ----------------------------------------------------------------

RippleDisplay::RippleDisplay(RippleDisplayCanvas* sdc, Viewport* v) :
    canvas(sdc), viewport(v)
{

    totalHeight = 1000;

}

RippleDisplay::~RippleDisplay()
{

}

void RippleDisplay::clear()
{
   if (ripplePlots.size() > 0)
   {
        for (int i = 0; i < ripplePlots.size(); i++)
        {
            ripplePlots[i]->clear();
        }
   }
        
}

void RippleDisplay::removePlots()
{
   ripplePlots.clear();
        
}

void RippleDisplay::addRipplePlot(int numChannels, int electrodeNum)
{

    std::cout << "Adding new Ripple plot." << std::endl;

    RipplePlot* ripplePlot = new RipplePlot(canvas, electrodeNum, 1000 + numChannels);
    ripplePlots.add(ripplePlot);
    addAndMakeVisible(ripplePlot);
}

void RippleDisplay::paint(Graphics& g)
{

    g.fillAll(Colours::grey);

}

void RippleDisplay::resized()
{
    // this is kind of a mess -- is there any way to optimize it?

    if (ripplePlots.size() > 0)
    {

        int w = getWidth();

        int numColumns = 1;
        int column, row;

        int stereotrodeStart = 0;
        int tetrodeStart = 0;

        int singlePlotIndex = -1;
        int stereotrodePlotIndex = -1;
        int tetrodePlotIndex = -1;
        int index;

        float width, height;


        float maxHeight = 0;

        for (int i = 0; i < ripplePlots.size(); i++)
        {

            if (ripplePlots[i]->nChannels == 1)
            {
                index = ++singlePlotIndex;
                numColumns = (int) jmax(w / ripplePlots[i]->minWidth, 1.0f);
                width = jmin((float) w / (float) numColumns, (float) getWidth());
                height = width * ripplePlots[i]->aspectRatio;

            }
            else if (ripplePlots[i]->nChannels == 2)
            {
                index = ++stereotrodePlotIndex;
                numColumns = (int) jmax(w / ripplePlots[i]->minWidth, 1.0f);
                width = jmin((float) w / (float) numColumns, (float) getWidth());
                height = width * ripplePlots[i]->aspectRatio;

            }
            else if (ripplePlots[i]->nChannels == 4)
            {
                index = ++tetrodePlotIndex;
                numColumns = (int) jmax(w / ripplePlots[i]->minWidth, 1.0f);
                width = jmin((float) w / (float) numColumns, (float) getWidth());
                height = width * ripplePlots[i]->aspectRatio;
            }

            column = index % numColumns;

            row = index / numColumns;

            ripplePlots[i]->setBounds(width*column, row*height, width, height);

            maxHeight = jmax(maxHeight, row*height + height);

            if (ripplePlots[i]->nChannels == 1)
            {
                stereotrodeStart = (int)(height*(float(row)+1));
            }
            else if (ripplePlots[i]->nChannels == 2)
            {
                tetrodeStart = (int)(height*(float(row)+1));
            }

        }


        for (int i = 0; i < ripplePlots.size(); i++)
        {

            int x = ripplePlots[i]->getX();
            int y = ripplePlots[i]->getY();
            int w2 = ripplePlots[i]->getWidth();
            int h2 = ripplePlots[i]->getHeight();

            if (ripplePlots[i]->nChannels == 2)
            {
                ripplePlots[i]->setBounds(x, y+stereotrodeStart, w2, h2);
                maxHeight = jmax(maxHeight, (float) y+stereotrodeStart+h2);

            }
            else if (ripplePlots[i]->nChannels == 4)
            {
                ripplePlots[i]->setBounds(x, y+stereotrodeStart+tetrodeStart, w2, h2);
                maxHeight = jmax(maxHeight, (float) y+stereotrodeStart+tetrodeStart+h2);
            }


        }

        totalHeight = (int) maxHeight + 50; 

       // std::cout << "New height = " << totalHeight << std::endl;

        setBounds(0, 0, getWidth(), totalHeight);
    }

}

void RippleDisplay::mouseDown(const MouseEvent& event)
{

}

void RippleDisplay::plotRipple(const RippleObject& Ripple, int electrodeNum)
{
    ripplePlots[electrodeNum]->processRippleObject(Ripple);
}


// ----------------------------------------------------------------

RipplePlot::RipplePlot(RippleDisplayCanvas* sdc, int elecNum, int p) :
     canvas(sdc), isSelected(false), electrodeNumber(elecNum),  plotType(p),
    limitsChanged(true)

{

    font = Font("Default", 15, Font::plain);

    switch (p)
    {
        case SINGLE_PLOT:
            // std::cout<<"RipplePlot as SINGLE_PLOT"<<std::endl;
            nWaveAx = 1;
            nProjAx = 0;
            nChannels = 1;
            minWidth = 200;
            aspectRatio = 1.0f;
            break;
        case STEREO_PLOT:
            //  std::cout<<"RipplePlot as STEREO_PLOT"<<std::endl;
            nWaveAx = 2;
            nProjAx = 1;
            nChannels = 2;
            minWidth = 300;
            aspectRatio = 0.5f;
            break;
        case TETRODE_PLOT:
            // std::cout<<"RipplePlot as TETRODE_PLOT"<<std::endl;
            nWaveAx = 4;
            nProjAx = 6;
            nChannels = 4;
            minWidth = 500;
            aspectRatio = 0.5f;
            break;
            //        case HIST_PLOT:
            //            nWaveAx = 1;
            //            nProjAx = 0;
            //            nHistAx = 1;
            //            break;
        default: // unsupported number of axes provided
            std::cout << "RipplePlot as UNKNOWN, defaulting to SINGLE_PLOT" << std::endl;
            nWaveAx = 1;
            nProjAx = 0;
            plotType = SINGLE_PLOT;
            nChannels = 1;
    }

    initAxes();

    for (int i = 0; i < nChannels; i++)
    {
        UtilityButton* rangeButton = new UtilityButton("250", Font("Small Text", 10, Font::plain));
        rangeButton->setRadius(3.0f);
        rangeButton->addListener(this);
        addAndMakeVisible(rangeButton);

        rangeButtons.add(rangeButton);
    }

}

RipplePlot::~RipplePlot()
{

}

void RipplePlot::paint(Graphics& g)
{

    g.setColour(Colours::white);
    g.drawRect(0,0,getWidth(),getHeight());

    g.setFont(font);

    g.drawText(String(electrodeNumber+1),10,0,50,20,Justification::left,false);

}

void RipplePlot::processRippleObject(const RippleObject& s)
{
    //std::cout<<"ElectrodePlot::processRippleObject()"<<std::endl;

    for (int i = 0; i < nWaveAx; i++)
        wAxes[i]->updateRippleData(s);

    for (int i = 0; i < nProjAx; i++)
        pAxes[i]->updateRippleData(s);
}

void RipplePlot::select()
{
    isSelected = true;
}

void RipplePlot::deselect()
{
    isSelected = false;
}

void RipplePlot::initAxes()
{
    initLimits();

    for (int i = 0; i < nWaveAx; i++)
    {
        RippleWaveAxes* wAx = new RippleWaveAxes(WAVE1 + i);
        wAxes.add(wAx);
        addAndMakeVisible(wAx);
        ranges.add(250.0f); // default range is 250 microvolts
    }

    for (int i = 0; i < nProjAx; i++)
    {
        RippleProjectionAxes* pAx = new RippleProjectionAxes(PROJ1x2 + i);
        pAxes.add(pAx);
        addAndMakeVisible(pAx);
    }

    setLimitsOnAxes(); // initialize the ranges
}

void RipplePlot::resized()
{

    float width = getWidth()-10;
    float height = getHeight()-25;

    float axesWidth, axesHeight;

    // to compute the axes positions we need to know how many columns of proj and wave axes should exist
    // using these two values we can calculate the positions of all of the sub axes
    int nProjCols, nWaveCols;

    switch (plotType)
    {
        case SINGLE_PLOT:
            nProjCols = 0;
            nWaveCols = 1;
            axesWidth = width;
            axesHeight = height;
            break;

        case STEREO_PLOT:
            nProjCols = 1;
            nWaveCols = 2;
            axesWidth = width/2;
            axesHeight = height;
            break;
        case TETRODE_PLOT:
            nProjCols = 3;
            nWaveCols = 2;
            axesWidth = width/4;
            axesHeight = height/2;
            break;
    }

    for (int i = 0; i < nWaveAx; i++)
    {
        wAxes[i]->setBounds(5 + (i % nWaveCols) * axesWidth/nWaveCols, 20 + (i/nWaveCols) * axesHeight, axesWidth/nWaveCols, axesHeight);
        rangeButtons[i]->setBounds(8 + (i % nWaveCols) * axesWidth/nWaveCols, 
                                   20 + (i/nWaveCols) * axesHeight + axesHeight - 18, 
                                   25, 15);
    }

    for (int i = 0; i < nProjAx; i++)
        pAxes[i]->setBounds(5 + (1 + i%nProjCols) * axesWidth, 20 + (i/nProjCols) * axesHeight, axesWidth, axesHeight);


}

void RipplePlot::buttonClicked(Button* button)
{
    UtilityButton* buttonThatWasClicked = (UtilityButton*) button;

    int index = rangeButtons.indexOf(buttonThatWasClicked);
    String label;

    if (ranges[index] == 250.0f)
    {
        ranges.set(index, 500.0f);
        label = "500";
    } else if (ranges[index] == 500.0f)
    {
        ranges.set(index, 100.0f);
        label = "100";
    } else if  (ranges[index] == 100.0f)
    {
        ranges.set(index, 250.0f);
        label = "250";
    } 

    buttonThatWasClicked->setLabel(label);

    setLimitsOnAxes();

}

void RipplePlot::setLimitsOnAxes()
{
    //std::cout<<"RipplePlot::setLimitsOnAxes()"<<std::endl;

    for (int i = 0; i < nWaveAx; i++)
         wAxes[i]->setRange(ranges[i]);

    // Each projection sets its limits using the limits of the two waveform dims it represents.
    // Convert projection number to indices, and then set the limits using those indices
    int j1, j2;
    for (int i = 0; i < nProjAx; i++)
    {
        pAxes[i]->n2ProjIdx(pAxes[i]->getType(), &j1, &j2);
        pAxes[i]->setRange(ranges[j1], ranges[j2]);
    }
}

void RipplePlot::initLimits()
{
    for (int i = 0; i < nChannels; i++)
    {
        limits[i][0] = 1209;//-1*pow(2,11);
        limits[i][1] = 11059;//pow(2,14)*1.6;
    }

}

void RipplePlot::getBestDimensions(int* w, int* h)
{
    switch (plotType)
    {
        case TETRODE_PLOT:
            *w = 4;
            *h = 2;
            break;
        case STEREO_PLOT:
            *w = 2;
            *h = 1;
            break;
        case SINGLE_PLOT:
            *w = 1;
            *h = 1;
            break;
        default:
            *w = 1;
            *h = 1;
            break;
    }
}

void RipplePlot::clear()
{
    std::cout << "RipplePlot::clear()" << std::endl;

    for (int i = 0; i < nWaveAx; i++)
        wAxes[i]->clear();
    for (int i = 0; i < nProjAx; i++)
        pAxes[i]->clear();
}



// --------------------------------------------------


RippleWaveAxes::RippleWaveAxes(int channel) : RippleGenericAxes(channel), drawGrid(true),
    bufferSize(10), RippleIndex(0), thresholdLevel(0.5f), range(250.0f),
    isOverThresholdSlider(false), isDraggingThresholdSlider(false), x(0.00)
{
/*
    addMouseListener(this, true);

    thresholdColour = Colours::red;

    font = Font("Small Text",10,Font::plain);

    for (int n = 0; n < bufferSize; n++)
    {
        RippleObject so;
        generateEmptyRipple(&so, 4);
        
        RippleBuffer.add(so);
    }*/
}

void RippleWaveAxes::setRange(float r)
{

    //std::cout << "Setting range to " << r << std::endl;

    range = r;

    repaint();
}

void RippleWaveAxes::paint(Graphics& g)
{
    g.setColour(Colours::black);
    g.fillRect(0,0,getWidth(), getHeight());

    int chan = 0;

    // draw the grid lines for the waveforms

    if (drawGrid)
        drawWaveformGrid(g);

    // draw the threshold line and labels
    drawThresholdSlider(g);

    // if no Ripples have been received then don't plot anything
    if (!gotFirstRipple)
    {
        return;
    }
   

    for (int RippleNum = 0; RippleNum < bufferSize; RippleNum++)
    {

       // if (RippleNum != RippleIndex)
        {
            g.setColour(Colours::white);
            plotRipple(RippleBuffer[RippleNum], g);
         }

    }

  //  g.setColour(Colours::white);
  //  plotRipple(RippleBuffer[RippleIndex], g);
    

    

}

void RippleWaveAxes::plotRipple(const RippleObject& s, Graphics& g)
{
/*
    float h = getHeight();

    //compute the spatial width for each waveform sample
  //  float dx = getWidth()/float(RippleBuffer[0].nSamples);



    float dx = getWidth()/float(9250); //(180ms/6ms)*(232+40) with some rounding off
    std::cout<<"the width by default and dx is "<<getWidth()<<" "<<dx<<std::endl;
    
    // type corresponds to channel so we need to calculate the starting
    // sample based upon which channel is getting plotted
    int sampIdx = 40*type; //RippleBuffer[0].nSamples * type; //

    int dSamples = 1;

    if (s.seqNumber <= 1)
    {
         x = 0.0f;
    }



     for (int i = 0; i < s.nSamples-1; i++)
    {
        //std::cout << s.data[sampIdx] << std::endl;

        if (*s.gain != 0)
        {
          //  float s1 = h/2 + float(s.data[sampIdx]-32768)/float(*s.gain)*1000.0f / range * h;
          //  float s2 =  h/2 + float(s.data[sampIdx+1]-32768)/float(*s.gain)*1000.0f / range * h;
            //temporarily removing ripple data
             g.drawLine(x,
                 s1, 
                 x+dx, 
                 s2);
        }

        

        sampIdx += dSamples;
        x += dx;
    }
*/
}

void RippleWaveAxes::drawThresholdSlider(Graphics& g)
{

    float h = getHeight()*thresholdLevel;

    g.setColour(thresholdColour);
    g.drawLine(0, h, getWidth(), h);

}

void RippleWaveAxes::drawWaveformGrid(Graphics& g)
{

    float h = getHeight();
    float w = getWidth();

    g.setColour(Colours::darkgrey);

    for (float y = -range/2; y < range/2; y += 25.0f)
    {
        g.drawLine(0,h/2 + y/range*h, w, h/2+ y/range*h);
    }
   
}

void RippleWaveAxes::updateRippleData(const RippleObject& s)
{
    if (!gotFirstRipple)
    {
        gotFirstRipple = true;
    }

    RippleObject newRipple = s;

    RippleIndex++;
    RippleIndex %= bufferSize;

    RippleBuffer.set(RippleIndex, newRipple);

    

}

void RippleWaveAxes::clear()
{
/*
    RippleBuffer.clear();
    RippleIndex = 0;

    for (int n = 0; n < bufferSize; n++)
    {
        RippleObject so;
        generateEmptyRipple(&so, 4);
        
        RippleBuffer.add(so);
    }*/
}

void RippleWaveAxes::mouseMove(const MouseEvent& event)
{

   // Point<int> pos = event.getPosition();

    float y = event.y;

    float h = getHeight()*thresholdLevel;

   // std::cout << y << " " << h << std::endl;

    if (y > h - 10.0f && y < h + 10.0f && !isOverThresholdSlider)
    {
        thresholdColour = Colours::yellow; 

      //  std::cout << "Yes." << std::endl;
        
        repaint();

        isOverThresholdSlider = true;

       // cursorType = MouseCursor::DraggingHandCursor;

    } else if ((y < h - 10.0f || y > h + 10.0f) && isOverThresholdSlider){

        thresholdColour = Colours::red;
        repaint();

        isOverThresholdSlider = false;

     //   cursorType = MouseCursor::NormalCursor;
        
    }


}

void RippleWaveAxes::mouseDown(const MouseEvent& event)
{
    // if (isOverThresholdSlider)
    // {
    //     cursorType = MouseCursor::DraggingHandCursor;
    // }
}

void RippleWaveAxes::mouseDrag(const MouseEvent& event)
{
    if (isOverThresholdSlider)
    {
        thresholdLevel = float(event.y) / float(getHeight());
        repaint();
    }
}

// MouseCursor RippleWaveAxes::getMouseCursor()
// {
//     MouseCursor c = MouseCursor(cursorType);

//     return c;
// }

void RippleWaveAxes::mouseExit(const MouseEvent& event)
{
    if (isOverThresholdSlider)
     {
        isOverThresholdSlider = false;
        thresholdColour = Colours::red; 
        repaint();
    }
}

// --------------------------------------------------

RippleProjectionAxes::RippleProjectionAxes(int projectionNum) : RippleGenericAxes(projectionNum), imageDim(500),
                                                    rangeX(250), rangeY(250)
{
    projectionImage = Image(Image::RGB, imageDim, imageDim, true);

    clear();
    //Graphics g(projectionImage);
    //g.setColour(Colours::red);
    //g.fillEllipse(20, 20, 300, 200);
    
    n2ProjIdx(projectionNum, &ampDim1, &ampDim2);


}

void RippleProjectionAxes::setRange(float x, float y)
{
    rangeX = (int) x;
    rangeY = (int) y;

    //std::cout << "Setting range to " << x << " " << y << std::endl;

    repaint();
}

void RippleProjectionAxes::paint(Graphics& g)
{
    //g.setColour(Colours::orange);
    //g.fillRect(5,5,getWidth()-5, getHeight()-5);

    g.drawImage(projectionImage,
                0, 0, getWidth(), getHeight(),
                0, imageDim-rangeY, rangeX, rangeY);
}

void RippleProjectionAxes::updateRippleData(const RippleObject& s)
{
    if (!gotFirstRipple)
    {
        gotFirstRipple = true;
    }

    int idx1, idx2;
    calcWaveformPeakIdx(s, ampDim1, ampDim2, &idx1, &idx2);

    // add peaks to image

    //updateProjectionImage(s.data[idx1], s.data[idx2], *s.gain);

}

void RippleProjectionAxes::updateProjectionImage(uint16_t x, uint16_t y, uint16_t gain)
{
    Graphics g(projectionImage);

   // h/2 + float(s.data[sampIdx]-32768)/float(*s.gain)*1000.0f / range * h;

    if (gain != 0)
    {
        float xf = float(x-32768)/float(gain)*1000.0f; // in microvolts
        float yf = float(imageDim) - float(y-32768)/float(gain)*1000.0f; // in microvolts

        g.setColour(Colours::white);
        g.fillEllipse(xf,yf,2.0f,2.0f); 
    }

}

void RippleProjectionAxes::calcWaveformPeakIdx(const RippleObject& s, int d1, int d2, int* idx1, int* idx2)
{

    int max1 = -1*pow(2.0,15);
    int max2 = max1;

  /*  for (int i = 0; i < s.nSamples; i++)
    {
        if (s.data[d1*s.nSamples + i] > max1)
        {
            *idx1 = d1*s.nSamples+i;
            max1 = s.data[*idx1];
        }
        if (s.data[d2*s.nSamples+i] > max2)
        {
            *idx2 = d2*s.nSamples+i;
            max2 = s.data[*idx2];
        }
    }*/
}



void RippleProjectionAxes::clear()
{
    projectionImage.clear(Rectangle<int>(0, 0, projectionImage.getWidth(), projectionImage.getHeight()),
                         Colours::black);
}

void RippleProjectionAxes::n2ProjIdx(int proj, int* p1, int* p2)
{
    int d1, d2;
    if (proj==PROJ1x2)
    {
        d1 = 0;
        d2 = 1;
    }
    else if (proj==PROJ1x3)
    {
        d1 = 0;
        d2 = 2;
    }
    else if (proj==PROJ1x4)
    {
        d1 = 0;
        d2 = 3;
    }
    else if (proj==PROJ2x3)
    {
        d1 = 1;
        d2 = 2;
    }
    else if (proj==PROJ2x4)
    {
        d1 = 1;
        d2 = 3;
    }
    else if (proj==PROJ3x4)
    {
        d1 = 2;
        d2 = 3;
    }
    else
    {
        std::cout<<"Invalid projection:"<<proj<<"! Cannot determine d1 and d2"<<std::endl;
        *p1 = -1;
        *p2 = -1;
        return;
    }
    *p1 = d1;
    *p2 = d2;
}

// --------------------------------------------------

RippleGenericAxes::RippleGenericAxes(int t)
    : gotFirstRipple(false), type(t)
{
    ylims[0] = 0;
    ylims[1] = 1;

    xlims[0] = 0;
    xlims[1] = 1;

    font = Font("Default", 12, Font::plain);

}

RippleGenericAxes::~RippleGenericAxes()
{

}

void RippleGenericAxes::updateRippleData(const RippleObject& newRipple)
{
    if (!gotFirstRipple)
    {
        gotFirstRipple = true;
    }

    s = newRipple;
}

void RippleGenericAxes::setYLims(double ymin, double ymax)
{

    //std::cout << "setting y limits to " << ymin << " " << ymax << std::endl;
    ylims[0] = ymin;
    ylims[1] = ymax;
}
void RippleGenericAxes::getYLims(double* min, double* max)
{
    *min = ylims[0];
    *max = ylims[1];
}
void RippleGenericAxes::setXLims(double xmin, double xmax)
{
    xlims[0] = xmin;
    xlims[1] = xmax;
}
void RippleGenericAxes::getXLims(double* min, double* max)
{
    *min = xlims[0];
    *max = xlims[1];
}


void RippleGenericAxes::setType(int t)
{
    if (t < WAVE1 || t > PROJ3x4)
    {
        std::cout<<"Invalid Axes type specified";
        return;
    }
    type = t;
}

int RippleGenericAxes::getType()
{
    return type;
}

int RippleGenericAxes::roundUp(int numToRound, int multiple)
{
    if (multiple == 0)
    {
        return numToRound;
    }

    int remainder = numToRound % multiple;
    if (remainder == 0)
        return numToRound;
    return numToRound + multiple - remainder;
}


void RippleGenericAxes::makeLabel(int val, int gain, bool convert, char* s)
{
    if (convert)
    {
        double volt = ad16ToUv(val, gain)/1000.;
        if (abs(val)>1e6)
        {
            //val = val/(1e6);
            sprintf(s, "%.2fV", volt);
        }
        else if (abs(val)>1e3)
        {
            //val = val/(1e3);
            sprintf(s, "%.2fmV", volt);
        }
        else
            sprintf(s, "%.2fuV", volt);
    }
    else
    {
        sprintf(s,"%d", (int)val);
    }
}

double RippleGenericAxes::ad16ToUv(int x, int gain)
{
    int result = (double)(x * 20e6) / (double)(gain * pow(2.0,16));
    return result;
}
