/*
  ==============================================================================

    This file contains the basic framework code for a JUCE plugin editor.

  ==============================================================================
*/

#include "PluginProcessor.h"
#include "PluginEditor.h"

//==============================================================================
PitchCorrectionAudioProcessorEditor::PitchCorrectionAudioProcessorEditor(PitchCorrectionAudioProcessor& p)
    : AudioProcessorEditor(&p), audioProcessor(p), genericParameterEditor(audioProcessor.getParameterManager())
{
    // Make sure that before the constructor has finished, you've set the
    // editor's size to whatever you need it to be.
    int height = static_cast<int>(audioProcessor.getParameterManager().getParameters().size())
        * genericParameterEditor.parameterWidgetHeight;
    setSize(300, height);
    addAndMakeVisible(genericParameterEditor);
}

PitchCorrectionAudioProcessorEditor::~PitchCorrectionAudioProcessorEditor()
{
}

//==============================================================================
void PitchCorrectionAudioProcessorEditor::paint(juce::Graphics& g)
{
    // (Our component is opaque, so we must completely fill the background with a solid colour)
    g.fillAll(getLookAndFeel().findColour(juce::ResizableWindow::backgroundColourId));

    //g.setColour (juce::Colours::white);
    //g.setFont (15.0f);
    //g.drawFittedText ("Hello World!", getLocalBounds(), juce::Justification::centred, 1);
}

void PitchCorrectionAudioProcessorEditor::resized()
{
    // This is generally where you'll want to lay out the positions of any
    // subcomponents in your editor..
    genericParameterEditor.setBounds(getLocalBounds());
}
