/*
  ==============================================================================

    ParameterEstimatorEditor.h
    Created: 21 May 2014 9:33:15am
    Author:  sethisan

  ==============================================================================
*/

#ifndef __PARAMETERESTIMATOREDITOR_H_2CD92CAA__
#define __PARAMETERESTIMATOREDITOR_H_2CD92CAA__


#include "../../../JuceLibraryCode/JuceHeader.h"
#include "GenericEditor.h"
#include "ElectrodeButtons.h"
#include "SpikeDetectorEditor.h"

class TriangleButton;
class UtilityButton;



/**

  Used to change the spike detection threshold.

  @see ParameterEstimatorEditor

*/

class ThresholdSlider; /*: public Slider
{
public:
    ThresholdSlider(Font f);
    ~ThresholdSlider() {}

    void setActive(bool);

    void setValues(Array<double>);

private:
    void paint(Graphics& g);

    Path makeRotaryPath(double, double, double);

    Font font;

    bool isActive;

    Array<double> valueArray;

};*/

/**

  User interface for the SpikeDetector processor.

  Allows the user to add single electrodes, stereotrodes, or tetrodes.

  Parameters of individual channels, such as channel mapping, threshold,
  and enabled state, can be edited.

  @see SpikeDetector

*/

class ParameterEstimatorEditor : public GenericEditor,
    public Label::Listener,
    public ComboBox::Listener

{
public:
    ParameterEstimatorEditor(GenericProcessor* parentNode, bool useDefaultParameterEditors);
    virtual ~ParameterEstimatorEditor();
    void buttonEvent(Button* button);
    void labelTextChanged(Label* label);
    void comboBoxChanged(ComboBox* comboBox);
    void sliderEvent(Slider* slider);

    void channelChanged(int chan);

    bool addElectrode(int nChans);
    void removeElectrode(int index);

    void checkSettings();
    void refreshElectrodeList();

private:

    void drawElectrodeButtons(int);

    ComboBox* electrodeTypes;
    ComboBox* electrodeList;
    ComboBox* svdCol;
    UtilityButton* setSVD;
    Label* numElectrodes;
    TriangleButton* upButton;
    TriangleButton* downButton;
    UtilityButton* plusButton;

    OwnedArray<ElectrodeButton> electrodeButtons;
    Array<ElectrodeEditorButton*> electrodeEditorButtons;
    
    void editElectrode(int index, int chan, int newChan);

    int lastId;
    bool isPlural;

    Font font;

    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(ParameterEstimatorEditor);

};




#endif  // __PARAMETERESTIMATOREDITOR_H_2CD92CAA__
