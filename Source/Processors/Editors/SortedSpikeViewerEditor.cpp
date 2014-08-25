/*
  ==============================================================================

    SortedSpikeViewerEditor.cpp
    Created: 18 Aug 2014 1:22:17pm
    Author:  sethisan

  ==============================================================================
*/

#include "SortedSpikeViewerEditor.h"
#include <string>

SortedSpikeDisplayEditor::SortedSpikeDisplayEditor(GenericProcessor* parentNode)
    : VisualizerEditor(parentNode,200)

{
    tabText = "Spike Sorting";

}

SortedSpikeDisplayEditor::~SortedSpikeDisplayEditor()
{
    deleteAllChildren();
}


void SortedSpikeDisplayEditor::buttonCallback(Button* button)
{

}

Visualizer* SortedSpikeDisplayEditor::createNewCanvas()
{

    SpikeSortingDisplay* processor = (SpikeSortingDisplay*) getProcessor();

    return new SortedSpikeCanvas(processor);
}
