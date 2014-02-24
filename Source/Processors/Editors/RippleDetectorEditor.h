/*
  ==============================================================================

    RippleDetectorEditor.h
    Created: 30 May 2013 12:59:19pm
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
#ifndef __RIPPLEDETECTOREDITOR_H_A4E5A8A7__
#define __RIPPLEDETECTOREDITOR_H_A4E5A8A7__



#include "../../../JuceLibraryCode/JuceHeader.h"
#include "SpikeDetectorEditor.h"
#include "GenericEditor.h"
#include "ElectrodeButtons.h"


class TriangleButton;
class UtilityButton;
class ThresholdSlider;


class RippleDetectorEditor : public GenericEditor,
    public Label::Listener,
    public ComboBox::Listener

{
public:
    RippleDetectorEditor(GenericProcessor* parentNode, bool useDefaultParameterEditors);
    virtual ~RippleDetectorEditor();
    void buttonEvent(Button* button);
    void labelTextChanged(Label* label);
    void comboBoxChanged(ComboBox* comboBox);
    void sliderEvent(Slider* slider);

    void channelChanged(int chan);

    bool addElectrode(int nChans);


    void checkSettings();
    void refreshElectrodeList();

private:

    void drawElectrodeButtons(int);



    ComboBox* electrodeTypes;
    ComboBox* electrodeList;
    Label* numElectrodes;
    Label* thresholdLabel;
    TriangleButton* upButton;
    TriangleButton* downButton;
    UtilityButton* plusButton;

    ThresholdSlider* thresholdSlider;

    OwnedArray<ElectrodeButton> electrodeButtons;
    Array<ElectrodeEditorButton*> electrodeEditorButtons;

    
    void removeElectrode(int index);
    void editElectrode(int index, int chan, int newChan);

    int lastId;
    bool isPlural;

    Font font;

    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(RippleDetectorEditor);

};




#endif  // __RIPPLEDETECTOREDITOR_H_A4E5A8A7__
