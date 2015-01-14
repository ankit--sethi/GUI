/*
  ==============================================================================

    SortedSpikeCanvas.cpp
    Created: 18 Aug 2014 1:35:47pm
    Author:  sethisan

  ==============================================================================
*/

#include "SortedSpikeCanvas.h"
#include "../ParameterEstimation.h"

SortedSpikeCanvas::SortedSpikeCanvas(SpikeSortingDisplay* n) :
    processor(n), newSpike(false)
{
    viewport = new Viewport();
    pcDisplay = new PCDisplay(this, viewport);
    viewport->setScrollBarsShown(true, false);
    viewport->setViewedComponent(pcDisplay, false);
    scrollBarThickness = viewport->getScrollBarThickness();

    clearButton = new UtilityButton("Clear plots", Font("Small Text", 13, Font::plain));
    clearButton->setRadius(3.0f);
    clearButton->addListener(this);
    addAndMakeVisible(clearButton);

    lockThresholdsButton = new UtilityButton("Lock Thresholds", Font("Small Text", 13, Font::plain));
    lockThresholdsButton->setRadius(3.0f);
    lockThresholdsButton->addListener(this);
    lockThresholdsButton->setClickingTogglesState(true);
    addAndMakeVisible(lockThresholdsButton);

    addAndMakeVisible(viewport);

    setWantsKeyboardFocus(true);

    update();

}

SortedSpikeCanvas::~SortedSpikeCanvas()
{

}

void SortedSpikeCanvas::update()
{

    //std::cout << "Updating SpikeDisplayCanvas" << std::endl;

    int nPlots = processor->getNumElectrodes();
    pcDisplay->removePlots();
    processor->removePCPlots();

    pcDisplay->svdCols = processor->svdCols;

    for (int i = 0; i < nPlots; i++)
    {
        PCPlot* sp = pcDisplay->addPCPlot(processor->getNumberOfChannelsForElectrode(i), i,processor->getNameForElectrode(i));
        processor->addPCPlotForElectrode(sp, i);
    }
    pcDisplay->resized();
    pcDisplay->repaint();
}

void SortedSpikeCanvas::paint(Graphics& g)
{

    g.fillAll(Colours::grey);

}

void SortedSpikeCanvas::refresh()
{
    processSpikeEvents();

    repaint();
}

void SortedSpikeCanvas::processSpikeEvents()
{

    processor->setParameter(2, 0.0f); // request redraw

}


void SortedSpikeCanvas::beginAnimation()
{
    std::cout << "SortedSpikeCanvas beginning animation." << std::endl;

    startCallbacks();
}

void SortedSpikeCanvas::endAnimation()
{
    std::cout << "SpikeDisplayCanvas ending animation." << std::endl;

    stopCallbacks();
}

void SortedSpikeCanvas::refreshState()
{
    // called when the component's tab becomes visible again
    resized();
}

void SortedSpikeCanvas::resized()
{
    viewport->setBounds(0,0,getWidth(),getHeight()-90);

    pcDisplay->setBounds(0,0,getWidth()-scrollBarThickness, pcDisplay->getTotalHeight());

    clearButton->setBounds(10, getHeight()-20, 100,20);

    lockThresholdsButton->setBounds(130, getHeight()-40, 130,20);

}

bool SortedSpikeCanvas::keyPressed(const KeyPress& key)
{

    KeyPress c = KeyPress::createFromDescription("c");

    if (key.isKeyCode(c.getKeyCode())) // C
    {

        std::cout << "Clearing display" << std::endl;
        return true;
    }

    return false;

}

void SortedSpikeCanvas::buttonClicked(Button* button)
{

    if (button == clearButton)
    {

    }
    else if (button == lockThresholdsButton)
    {

    }
}

PCDisplay::PCDisplay(SortedSpikeCanvas* sdc, Viewport* v) :
    canvas(sdc), viewport(v)
{

    totalHeight = 1000;

}

PCDisplay::~PCDisplay()
{

}

void PCDisplay::clear()
{
    if (pcPlots.size() > 0)
    {
        for (int i = 0; i < pcPlots.size(); i++)
        {
            pcPlots[i]->clear();
        }
    }

}


void PCDisplay::removePlots()
{
    pcPlots.clear();

}

PCPlot* PCDisplay::addPCPlot(int numChannels, int electrodeNum, String name_)
{

    std::cout << "Adding new spike plot and svdCols is." << svdCols << std::endl;

    PCPlot* pcPlot = new PCPlot(canvas, electrodeNum, svdCols, name_);
    //std::cout << "Done making Spike plot" << std:: endl;
    pcPlots.add(pcPlot);
    addAndMakeVisible(pcPlot);
    return pcPlot;
}

void PCDisplay::paint(Graphics& g)
{

    g.fillAll(Colours::yellow);

}

void PCDisplay::resized()
{
    // this is kind of a mess -- is there any way to optimize it?

    if (pcPlots.size() > 0)
    {

        int w = getWidth();

        int numColumns = 1;
        int column, row;

        int PlotIndex = -1;
        int index = -1;

        float width, height;


        float maxHeight = 0;
        for (int i = 0; i < pcPlots.size(); i++)
        {
            index = ++PlotIndex;
            numColumns = (int) jmax(w / pcPlots[i]->minWidth, 1.0f);
            width = jmin((float) w / (float) numColumns, (float) getWidth());
            height = width * pcPlots[i]->aspectRatio + 150;

            column = index % numColumns;

            row = index / numColumns;

            pcPlots[i]->setBounds(width*column, row*height, width, height);

            //pcPlots[i]->setBounds(5, 5, width, height);

            maxHeight = jmax(maxHeight, row*height + height);
        }

        totalHeight = (int) maxHeight + 50;

        //std::cout << "New height = " << totalHeight << std::endl;

        setBounds(0, 0, getWidth(), totalHeight);
        //setBounds(0, 0, 10, 10);
    }

}

void PCDisplay::mouseDown(const MouseEvent& event)
{

}

void PCDisplay::plotSpike(const SortedSpikeObject& spike, int electrodeNum)
{
    pcPlots[electrodeNum]->processSortedSpikeObject(spike);
}


PCPlot::PCPlot(SortedSpikeCanvas *sdc, int elecNum, int dim, String name_) :
    canvas(sdc), isSelected(false), electrodeNumber(elecNum), svdCol(dim), limitsChanged(true), name(name_)
{

    font = Font("Default", 15, Font::plain);
    switch (svdCol + 7000)
    {
    case K_IS_TWO:
        nPCAx = 1;
        minWidth = 200;
        aspectRatio = 1.0f;
        break;
    case K_IS_THREE:
        nPCAx = 3;
        minWidth = 600;
        aspectRatio = 0.334f;
        break;
    case K_IS_FOUR:
        nPCAx = 6;
        minWidth = 600;
        aspectRatio = 0.667f;
        break;
    case K_IS_FIVE:
        nPCAx = 10;
        minWidth = 1000;
        aspectRatio = 0.4f;
        break;

    default: // unsupported number of axes provided
        std::cout << "SpikePlot as UNKNOWN, defaulting to SINGLE_PLOT" << std::endl;
        nPCAx = 3;
        svdCol = K_IS_THREE;
    }

    initAxes();
    raster = new RasterPlot;
    raster->setBounds(5,0.1*getHeight() ,10,10);
    addAndMakeVisible(raster);

}

PCPlot::~PCPlot()
{

}

void PCPlot::paint(Graphics& g)
{

    g.fillAll(Colours::green);

    g.setFont(font);

    g.drawText(name,10,0,200,20,Justification::left,false);

    g.fillAll(Colours::green);

}

void PCPlot::resized()
{

    float width = getWidth()-10;
    float height = getHeight() - 150 -25;

    float axesWidth, axesHeight;

    // to compute the axes positions we need to know how many columns of proj and wave axes should exist
    // using these two values we can calculate the positions of all of the sub axes
    int nPCCols;
    std::cout << "ATTENTION " << svdCol;
    switch (svdCol + 7000)
    {
    case K_IS_TWO:
        nPCCols = 1;
        axesWidth = width;
        axesHeight = height;
        break;
    case K_IS_THREE:
        nPCCols = 3;
        axesWidth = width/3;
        axesHeight = height;
        break;
    case K_IS_FOUR:
        nPCCols = 3;
        axesWidth = width/3;
        axesHeight = height/2;
        break;
    case K_IS_FIVE:
        nPCCols = 5;
        axesWidth = width/5;
        axesHeight = height/2;
        break;
    }

    for (int i = 0; i < nPCAx; i++)
    {
        std::cout << nPCAx;
        pAxes[i]->setBounds(5 + (i%nPCCols) * axesWidth, 20 + (i/nPCCols) * axesHeight, axesWidth, axesHeight);
        pAxes[i]->projectionImage = pAxes[i]->projectionImage.rescaled(pAxes[i]->getWidth(),pAxes[i]->getHeight());
    }
    raster->setBounds(5,axesHeight + 27 ,nPCCols*axesWidth,50);
    raster->ProjectionImage = raster->ProjectionImage.rescaled(nPCCols*axesWidth, 50);
    //raster->startTimeStampPixel = 0.025*raster->getWidth();
    //raster->endTimeStampPixel = 0.975*raster->getWidth();
}


void PCPlot::initAxes()
{
    initLimits();
    //std::cout << "nPCAx" << nPCAx << "//";

    for (int i = 0; i < nPCAx; i++)
    {
        PCAxes* pAx = new PCAxes(PC1x2 + i, svdCol);
        pAxes.add(pAx);
        addAndMakeVisible(pAx);
    }

    setLimitsOnAxes(); // initialize the ranges
}

void PCPlot::initLimits()
{
    /*for (int i = 0; i < nChannels; i++)
    {
        limits[i][0] = 1209;//-1*pow(2,11);
        limits[i][1] = 11059;//pow(2,14)*1.6;
    }*/

}

void PCPlot::setLimitsOnAxes()
{
    for (int i = 0; i < nPCAx; i++)
    {
        pAxes[i]->setRange(20, 20);
    }
}

void PCPlot::select()
{
    isSelected = true;
}

void PCPlot::deselect()
{
    isSelected = false;
}

void PCPlot::processSortedSpikeObject(const SortedSpikeObject& s)
{
    //std::cout << "processing";
    for (int i = 0; i < nPCAx; i++)
        pAxes[i]->updateSpikeData(s, i);

}

void PCPlot::processRasterPlot(const RasterArray& r, int electrodeNum)
{
    raster->updateRasterData(r, electrodeNum);
    //raster->repaint();
}

void PCPlot::clear()
{
    std::cout << "SpikePlot::clear()" << std::endl;

    for (int i = 0; i < nPCAx; i++)
        pAxes[i]->clear();
}

void PCPlot::buttonClicked(Button* button)
{
    UtilityButton* buttonThatWasClicked = (UtilityButton*) button;

    int index = rangeButtons.indexOf(buttonThatWasClicked);
    String label;

    if (ranges[index] == 0.1f)
    {
        ranges.set(index, 0.1f);
        label = "0.1";
    }
    else if (ranges[index] == 0.5f)
    {
        ranges.set(index, 0.5f);
        label = "0.5";
    }
    else if (ranges[index] == 1.0f)
    {
        ranges.set(index, 1.0f);
        label = "1.0";
    }

    buttonThatWasClicked->setLabel(label);

    setLimitsOnAxes();

}
RasterPlot::RasterPlot()
{
    ProjectionImage = Image(Image::RGB, 50, 50, true);

    clear();
    repaint();

}

void RasterPlot::situateRasterData(const RasterData& s, unsigned long startTimestamp, unsigned long endTimeStamp)
{

    Graphics g(ProjectionImage);
    g.setColour(Colours::white);
    g.fillAll();
    g.setColour(Colours::darkblue);

    switch(s.neuronID)
    {
    case 1:
        g.setColour(Colours::red);
        break;
    case 2:
        g.setColour(Colours::blue);
        break;
    case 3:
        g.setColour(Colours::green);
        break;
    case 4:
        g.setColour(Colours::yellow);
        break;
    case 5:
        g.setColour(Colours::violet);
        break;
    case 6:
        g.setColour(Colours::gold);
        break;
    case 7:
        g.setColour(Colours::darkorchid);
        break;
    case 8:
        g.setColour(Colours::mintcream);
        break;
    case 9:
        g.setColour(Colours::oldlace);
        break;
    case 10:
        g.setColour(Colours::orange);
        break;
    case 11:
        g.setColour(Colours::palegreen);
        break;
    case 12:
        g.setColour(Colours::peachpuff);
        break;
    case 13:
        g.setColour(Colours::powderblue);
        break;
    case 14:
        g.setColour(Colours::saddlebrown);
        break;
    case 15:
        g.setColour(Colours::seashell);
        break;
    case 16:
        g.setColour(Colours::slateblue);
        break;
    case 17:
        g.setColour(Colours::steelblue);
        break;
    case 18:
        g.setColour(Colours::tomato);
        break;
    case 19:
        g.setColour(Colours::crimson);
        break;
    }

    if (s.timestamp < startTimestamp )
    {
    std::cout << "Error: Neuron detected before Ripple Start!" << std::endl;
    }
    std::cout << "start and end is " << startTimestamp << "   " << endTimeStamp << std::endl;
    float length = float(endTimeStamp - startTimestamp);

    float ratio = float(s.timestamp - startTimestamp)/length;

    g.drawLine(ratio*getWidth(), 0.9*getHeight(), ratio*getWidth(), 0.1*getHeight(), 1.0f);

}

bool RasterPlot::updateRasterData(const RasterArray& s, int electrodeNum)
{
    unsigned long int startTimestamp, stopTimestamp;

    startTimestamp = s.startRipple;
    stopTimestamp = s.stopRipple;

    int count = 0;

    for(int i = 0; i < s.accruedRasterMarks.size(); i++)
    {
        if(s.accruedRasterMarks[i].neuronID != -1)
        {
            count++;
            situateRasterData(s.accruedRasterMarks[i], startTimestamp, stopTimestamp);
        }
    }
    std::cout << "In one plot, I entered " << count << " times." << std::endl;
}

void RasterPlot::clear()
{
    ProjectionImage.clear(Rectangle<int>(0, 0, ProjectionImage.getWidth(), ProjectionImage.getHeight()),
                          Colours::white);

    repaint();
}
void RasterPlot::paint(Graphics& g)
{

   g.setColour(Colours::pink);
   g.fillRect(1,1,getWidth(), getHeight());

   g.drawImage(ProjectionImage, 1, 1, getWidth() - 1, getHeight() - 1, 0, 0, getWidth(), getHeight());
   g.setColour(Colours::black);
   //g.drawLine(startTimeStampPixel, 0.9*getHeight(), endTimeStampPixel, 0.9*getHeight(), 2.0f);

   // draw x axes
}

PCAxes::PCAxes(int projectionNum, int dictlength) : GenericAxes(projectionNum), rangeX(30.0), rangeY(30.0), spikesReceivedSinceLastRedraw(0), svdCols(dictlength)
{
    projectionImage = Image(Image::RGB, 50, 50, true);

    clear();
    repaint();
}

bool PCAxes::updateSpikeData(const SortedSpikeObject& s)
{
    if (!gotFirstSpike)
    {
        gotFirstSpike = true;
    }

    return true;
}

bool PCAxes::updateSpikeData(const SortedSpikeObject& s, int index)
{
    if (!gotFirstSpike)
    {
        gotFirstSpike = true;
    }

    updatePrincipalComponents(s, index);

    return true;
}

void  PCAxes::updatePrincipalComponents(const SortedSpikeObject& s, int index)
{
    float h = getHeight();
    float w = getWidth();

    int choose2[2] = {0,1};
    int choose3[6] = {0,1,0,2,1,2};
    int choose4[12] = {0,1,0,2,0,3,1,2,1,3,2,3};
    int choose5[20] = {0,1,0,2,0,3,0,4,1,2,1,3,1,4,2,3,2,4,3,4};
    float xf, yf;
    switch(svdCols)
    {
    case 2:
        xf = float(s.principalComponent[choose2[2*(index)]]-32768)/float(*s.gain)*1000.0f;
        yf = float(s.principalComponent[choose2[2*(index) + 1]]-32768)/float(*s.gain)*1000.0f;
        break;
    case 3:
        xf = float(s.principalComponent[choose3[2*(index)]]-32768)/float(*s.gain)*1000.0f;
        yf = float(s.principalComponent[choose3[2*(index) + 1]]-32768)/float(*s.gain)*1000.0f;
        break;
    case 4:
        xf = float(s.principalComponent[choose4[2*(index)]]-32768)/float(*s.gain)*1000.0f;
        yf = float(s.principalComponent[choose4[2*(index) + 1]]-32768)/float(*s.gain)*1000.0f;
        break;
    case 5:
        xf = float(s.principalComponent[choose5[2*(index)]]-32768)/float(*s.gain)*1000.0f;
        yf = float(s.principalComponent[choose5[2*(index) + 1]]-32768)/float(*s.gain)*1000.0f;
        break;
    }




    xf = 0.5*w + 0.5*(w*xf/rangeX);
    yf = 0.5*h + 0.5*(w*yf/rangeY);

    Graphics g(projectionImage);

    g.setColour(Colours::white);
    g.fillAll();
    //g.fillEllipse(xf,yf,4.0f,4.0f);

    g.setColour(Colours::orange);

    switch(s.neuronID)
    {
    case 1:
        g.setColour(Colours::green);
        break;
    case 2:
        g.setColour(Colours::red);
        break;
    case 3:
        g.setColour(Colours::gold);
        break;
    case 4:
        g.setColour(Colours::yellow);
        break;
    case 5:
        g.setColour(Colours::violet);
        break;
    case 6:
        g.setColour(Colours::blue);
        break;
    case 7:
        g.setColour(Colours::darkorchid);
        break;
    case 8:
        g.setColour(Colours::mintcream);
        break;
    case 9:
        g.setColour(Colours::oldlace);
        break;
    case 10:
        g.setColour(Colours::darkblue);
        break;
    case 11:
        g.setColour(Colours::palegreen);
        break;
    case 12:
        g.setColour(Colours::peachpuff);
        break;
    case 13:
        g.setColour(Colours::powderblue);
        break;
    case 14:
        g.setColour(Colours::saddlebrown);
        break;
    case 15:
        g.setColour(Colours::seashell);
        break;
    case 16:
        g.setColour(Colours::slateblue);
        break;
    case 17:
        g.setColour(Colours::steelblue);
        break;
    case 18:
        g.setColour(Colours::tomato);
        break;
    case 19:
        g.setColour(Colours::crimson);
        break;
    }
        //std::cout << "XF / width IS" << xf  << "/" << w << " And YF / height is" << yf << "/" << h << std::endl;
        g.fillEllipse(xf,yf,6.0f,6.0f);
        //g.fillEllipse(3,4,50.0f,50.0f);

}
void PCAxes::drawGrid(Graphics& g)
{

    float h = getHeight();
    float w = getWidth();

    g.setColour(Colours::blue);

    g.drawLine(0, h/2, w, h/2, 2.0f);
    g.drawLine(w/2, 0, w/2, h, 2.0f);


}
void PCAxes::paint(Graphics& g)
{

    g.setColour(Colours::white);
    g.fillRect(1,1,getWidth() - 1 , getHeight() -1);
    //std::cout << "drawing again";
    g.drawImage(projectionImage, 1, 1, getWidth() - 1, getHeight() - 1, 0, 0, getWidth(), getHeight());
    drawGrid(g);
}

void PCAxes::clear()
{
    projectionImage.clear(Rectangle<int>(0, 0, projectionImage.getWidth(), projectionImage.getHeight()),
                          Colours::white);

    repaint();
}

void PCAxes::setRange(float x, float y)
{
    rangeX = (int) x;
    rangeY = (int) y;

    //std::cout << "Setting range to " << x << " " << y << std::endl;

    repaint();
}
