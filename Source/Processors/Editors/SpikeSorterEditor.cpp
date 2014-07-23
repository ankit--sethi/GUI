/*
  ==============================================================================

    SpikeSorterEditor.cpp
    Created: 23 Jun 2014 6:17:31pm
    Author:  sethisan

  ==============================================================================
*/

#include "../SpikeSorter.h"
#include "SpikeSorterEditor.h"
#include "ChannelSelector.h"
#include "../../UI/EditorViewport.h"
#include <stdio.h>

SpikeSorterEditor::SpikeSorterEditor(GenericProcessor* parentNode, bool useDefaultParameterEditors=true)
    : GenericEditor(parentNode, useDefaultParameterEditors)

{



}

SpikeSorterEditor::~SpikeSorterEditor()
{


}
void SpikeSorterEditor::labelTextChanged(Label* label)
{

    getEditorViewport()->makeEditorVisible(this, false, true);

}

void SpikeSorterEditor::comboBoxChanged(ComboBox* comboBox)
{


}
