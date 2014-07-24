/*
  ==============================================================================

    SpikeSorterEditor.h
    Created: 23 Jun 2014 6:17:31pm
    Author:  sethisan

  ==============================================================================
*/

#ifndef __SPIKESORTEREDITOR_H_48FCE7AC__
#define __SPIKESORTEREDITOR_H_48FCE7AC__

#include "../../../JuceLibraryCode/JuceHeader.h"
#include "GenericEditor.h"
#include "ElectrodeButtons.h"

class SpikeSorterEditor : public GenericEditor,
    public Label::Listener,
    public ComboBox::Listener

{
    public:
    SpikeSorterEditor(GenericProcessor* parentNode, bool useDefaultParameterEditors);
    virtual ~SpikeSorterEditor();
    void labelTextChanged(Label* label);
    void comboBoxChanged(ComboBox* comboBox);


    private:
    ComboBox* electrodeTypes;

    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(SpikeSorterEditor);

};

#endif  // __SPIKESORTEREDITOR_H_48FCE7AC__
