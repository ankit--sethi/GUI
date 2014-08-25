/*
  ==============================================================================

    SortedSpikeViewerEditor.h
    Created: 18 Aug 2014 1:22:17pm
    Author:  sethisan

  ==============================================================================
*/

#ifndef __SORTEDSPIKEVIEWEREDITOR_H_D53F42F6__
#define __SORTEDSPIKEVIEWEREDITOR_H_D53F42F6__

#include "../../../JuceLibraryCode/JuceHeader.h"
#include "GenericEditor.h"
#include "../../UI/UIComponent.h"
#include "../../UI/DataViewport.h"
#include "../Visualization/DataWindow.h"
#include "../SortedSpikeViewer.h"
#include "../Visualization/SortedSpikeCanvas.h"
#include "VisualizerEditor.h"

class Visualizer;
class UtilityButton;

/**

  User interface for the Sorted Spike Viewer sink.

  @see SpikeDisplayNode, SpikeDisplayCanvas

*/

class SortedSpikeDisplayEditor : public VisualizerEditor
{
    public:
    SortedSpikeDisplayEditor(GenericProcessor*);
    ~SortedSpikeDisplayEditor();

    void buttonCallback(Button* button);
    Visualizer* createNewCanvas();

    private:

    UtilityButton* panUpBtn;
    UtilityButton* panDownBtn;
    UtilityButton* zoomInBtn;
    UtilityButton* zoomOutBtn;
    UtilityButton* clearBtn;
    UtilityButton* saveImgBtn;

    Label* panLabel;
    Label* zoomLabel;

    UtilityButton* allSubChansBtn;

    int nSubChannels;
    Label* subChanLabel;

    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(SortedSpikeDisplayEditor);
};

#endif  // __SORTEDSPIKEVIEWEREDITOR_H_D53F42F6__
