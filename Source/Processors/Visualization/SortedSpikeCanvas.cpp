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

    g.fillAll(Colours::black);

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
            height = width * pcPlots[i]->aspectRatio;

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
}

PCPlot::~PCPlot()
{

}

void PCPlot::paint(Graphics& g)
{

    g.fillAll(Colours::white);

    g.setFont(font);

    g.drawText(name,10,0,200,20,Justification::left,false);

}

void PCPlot::resized()
{

    float width = getWidth()-10;
    float height = getHeight()-25;

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
        pAxes[i]->setRange(2, 2);
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
        pAxes[i]->updateSpikeData(s);

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

PCAxes::PCAxes(int projectionNum, int dictlength) : GenericAxes(projectionNum), rangeX(2.0), rangeY(2.0), spikesReceivedSinceLastRedraw(0), svdCols(dictlength)
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

    updatePrincipalComponents(s);

    return true;
}

void  PCAxes::updatePrincipalComponents(const SortedSpikeObject& s)
{
    float h = getHeight();
    float w = getWidth();
    Graphics g(projectionImage);

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
    }

    int choose2[2] = {1,2};
    int choose3[6] = {1,2,1,3,2,3};
    int choose4[12] = {1,2,1,3,1,4,2,3,2,4,3,4};
    int choose5[20] = {1,2,1,3,1,4,1,5,2,3,2,4,2,5,3,4,3,5,4,5};
    float xf, yf;
    switch(svdCols)
    {
    case 2:
        xf = s.principalComponent[choose2[2*(type - 1)]];
        yf = s.principalComponent[choose2[2*type - 1]];
        break;
    case 3:
        xf = s.principalComponent[choose3[2*(type - 1)]];
        yf = s.principalComponent[choose3[2*type - 1]];
        break;
    case 4:
        xf = s.principalComponent[choose4[2*(type - 1)]];
        yf = s.principalComponent[choose4[2*type - 1]];
        break;
    case 5:
        xf = s.principalComponent[choose5[2*(type - 1)]];
        yf = s.principalComponent[choose5[2*type - 1]];
        break;
    }

        xf = 0.5*w + (w*xf/rangeX);
        yf = 0.5*h + (w*yf/rangeY);
        g.fillEllipse(xf,yf,2.0f,2.0f);
        g.fillEllipse(3,4,50.0f,50.0f);

}
void PCAxes::drawGrid(Graphics& g)
{

    float h = getHeight();
    float w = getWidth();

    g.setColour(Colours::pink);

    g.drawLine(0, h/2, w, h/2, 2.0f);
    g.drawLine(w/2, 0, w/2, h, 2.0f);


}
void PCAxes::paint(Graphics& g)
{

    g.setColour(Colours::black);
    g.fillRect(1,1,getWidth() - 1 , getHeight() -1);

    g.drawImage(projectionImage, 1, 1, getWidth() - 1, getHeight() - 1, 0, 0, getWidth(), getHeight());
    drawGrid(g);
}

void PCAxes::clear()
{
    projectionImage.clear(Rectangle<int>(0, 0, projectionImage.getWidth(), projectionImage.getHeight()),
                          Colours::black);

    repaint();
}

void PCAxes::setRange(float x, float y)
{
    rangeX = (int) x;
    rangeY = (int) y;

    //std::cout << "Setting range to " << x << " " << y << std::endl;

    repaint();
}
