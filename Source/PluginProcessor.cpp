/*
  ==============================================================================

    This file contains the basic framework code for a JUCE plugin processor.

  ==============================================================================
*/

#include "PluginProcessor.h"
#include "PluginEditor.h"
#include <algorithm>
#define M_PI 3.14159265358979323846264338327950288

static const std::vector<mrta::ParameterInfo> ParameterInfos
{
    // "units", default, min, max, step, skew;
    { Param::ID::Bypass,   Param::Name::Bypass,    "OFF", "ON", true },
    { Param::ID::Semitones,Param::Name::Semitones, "", 0.f, -12.f, 12.f, 0.01f, 1.f },
    { Param::ID::Autotune, Param::Name::Autotune,  "Autotune Off", "Autotune On", false },
    { Param::ID::Noob,  Param::Name::Noob,   "Noob Off", "Noob On", false },
    { Param::ID::Blend,    Param::Name::Blend,       "", 1.f, 0.f, 1.f, 0.01f, 1.f },
    { Param::ID::Volume,   Param::Name::Volume,    "dB", 0.0f, -60.f, 12.f, 0.1f, 3.8018f } // 0 at exact midpoint 3.8018f
};


//==============================================================================
PitchCorrectionAudioProcessor::PitchCorrectionAudioProcessor() : parameterManager(*this, ProjectInfo::projectName, ParameterInfos),
                                                fft(fftOrder)

{

    //parameterManager.registerParameterCallback(Param::ID::Gain,
    //    [this](float value, bool /*forced*/)
    //    {
    //        DBG(Param::Name::Gain + ": " + juce::String{ value });
    //        TanhDist1.setGain(value);
    //    });

    parameterManager.registerParameterCallback(Param::ID::Bypass,
        [this](float newValue, bool force)
        {
            enableRamp.setCurrentAndTargetValue(std::fmin(std::fmax(newValue, 0.f), 1.f));
        });

    parameterManager.registerParameterCallback(Param::ID::Semitones,
        [this](float value, bool /*forced*/)
        {
            DBG(Param::Name::Semitones + ": " + juce::String{ value });
            //shift = powf(2.f, (float)value / 12.f);
            shiftRamp.setCurrentAndTargetValue(powf(2.f, (float)value / 12.f));
        });
    parameterManager.registerParameterCallback(Param::ID::Autotune,
        [this](float newValue, bool force)
        {
            autoTuneRamp.setCurrentAndTargetValue(std::fmin(std::fmax(newValue, 0.f), 1.f));
        });

    parameterManager.registerParameterCallback(Param::ID::Noob,
        [this](float newValue, bool force)
        {
            noobRamp.setCurrentAndTargetValue(std::fmin(std::fmax(newValue, 0.f), 1.f));
        });

    parameterManager.registerParameterCallback(Param::ID::Blend,
        [this](float value, bool /*forced*/)
        {
            //DBG(Param::Name::Semitones + ": " + juce::String{ value });
            //shift = powf(2.f, (float)value / 12.f);
            blendRamp.setCurrentAndTargetValue(value);
        });

    parameterManager.registerParameterCallback(Param::ID::Volume,
        [this](float value, bool forced)
        {
            DBG(Param::Name::Volume + ": " + juce::String{ value });
            //TanhDist1.setVolume(value);

            float dbValue{ 0.f }; // mute
            if (value > -60.f)
                dbValue = std::pow(10.f, value * 0.05f);

            if (forced)
                outputGain.setCurrentAndTargetValue(dbValue);
            else
                outputGain.setTargetValue(dbValue);

        });

}

PitchCorrectionAudioProcessor::~PitchCorrectionAudioProcessor()
{
}

//==============================================================================
const juce::String PitchCorrectionAudioProcessor::getName() const
{
    return JucePlugin_Name;
}

bool PitchCorrectionAudioProcessor::acceptsMidi() const
{
#if JucePlugin_WantsMidiInput
    return true;
#else
    return false;
#endif
}

bool PitchCorrectionAudioProcessor::producesMidi() const
{
#if JucePlugin_ProducesMidiOutput
    return true;
#else
    return false;
#endif
}

bool PitchCorrectionAudioProcessor::isMidiEffect() const
{
#if JucePlugin_IsMidiEffect
    return true;
#else
    return false;
#endif
}

double PitchCorrectionAudioProcessor::getTailLengthSeconds() const
{
    return 0.0;
}

int PitchCorrectionAudioProcessor::getNumPrograms()
{
    return 1;   // NB: some hosts don't cope very well if you tell them there are 0 programs,
                // so this should be at least 1, even if you're not really implementing programs.
}

int PitchCorrectionAudioProcessor::getCurrentProgram()
{
    return 0;
}

void PitchCorrectionAudioProcessor::setCurrentProgram(int index)
{
}

const juce::String PitchCorrectionAudioProcessor::getProgramName(int index)
{
    return {};
}

void PitchCorrectionAudioProcessor::changeProgramName(int index, const juce::String& newName)
{
}

//==============================================================================
void PitchCorrectionAudioProcessor::prepareToPlay(double sampleRate, int samplesPerBlock)
{
    // Use this method as the place to do any pre-playback
    // initialisation that you need..
    //TanhDist1.prepare(sampleRate, samplesPerBlock);
    
    f0_instance.prepare(sampleRate, fftSize);

    outputGain.reset(sampleRate, 0.01f);
    enableRamp.reset(sampleRate, 0.001f);
    autoTuneRamp.reset(sampleRate, 0.001f);
    autoShift.reset(sampleRate, 0.00009f);
    shiftRamp.reset(sampleRate, 0.001f);
    blendRamp.reset(sampleRate, 0.001f);
    noobRamp.reset(sampleRate, 0.001f);

    makeHannWindow(anaWindow, fftSize);
    computeOLA_coeff();
    inputBufferLength = fftSize;
    inputBufferWritePosition = 0;
    inputBuffer.clear();
    inputBuffer.setSize(getTotalNumInputChannels(), inputBufferLength);
    float maxRatio = powf(2.f, -12.f / 12.f);
    outputBufferLength = (int)floorf((float)fftSize / maxRatio);
    outputBufferWritePosition = 0;
    outputBufferReadPosition = 0;
    outputBuffer.clear();
    outputBuffer.setSize(getTotalNumInputChannels(), outputBufferLength);



    fftDataTime.realloc(fftSize);
    fftDataTime.clear(fftSize);

    fftDataFreq.realloc(fftSize);
    fftDataFreq.clear(fftSize);

    samplesSinceLastFFT = 0;


    omega.realloc(fftSize);
    for (int index = 0; index < fftSize; ++index)
        omega[index] = 2.0f * M_PI * index / (float)fftSize;

    inputPhase.clear();
    inputPhase.setSize(getTotalNumInputChannels(), outputBufferLength);

    outputPhase.clear();
    outputPhase.setSize(getTotalNumInputChannels(), outputBufferLength);

}

void PitchCorrectionAudioProcessor::releaseResources()
{
    // When playback stops, you can use this as an opportunity to free up any
    // spare memory, etc.
}

#ifndef JucePlugin_PreferredChannelConfigurations
bool PitchCorrectionAudioProcessor::isBusesLayoutSupported(const BusesLayout& layouts) const
{
#if JucePlugin_IsMidiEffect
    juce::ignoreUnused(layouts);
    return true;
#else
    // This is the place where you check if the layout is supported.
    // In this template code we only support mono or stereo.
    // Some plugin hosts, such as certain GarageBand versions, will only
    // load plugins that support stereo bus layouts.
    if (layouts.getMainOutputChannelSet() != juce::AudioChannelSet::mono()
        && layouts.getMainOutputChannelSet() != juce::AudioChannelSet::stereo())
        return false;

    // This checks if the input layout matches the output layout
#if ! JucePlugin_IsSynth
    if (layouts.getMainOutputChannelSet() != layouts.getMainInputChannelSet())
        return false;
#endif

        return true;
#endif
}
#endif


float PitchCorrectionAudioProcessor::princArg(const float phase)
{
    //return fmod(phase + M_PI, -2.0f * M_PI) + M_PI;
    if (phase >= 0.0f)
        return fmod(phase + M_PI, 2.0f * M_PI) - M_PI;
    else
        return fmod(phase + M_PI, -2.0f * M_PI) + M_PI;
}

void PitchCorrectionAudioProcessor::makeHannWindow(float* window, const int winSize)
{
    for (int sample = 0; sample < winSize; ++sample)
        window[sample] = 0.5f - 0.5f * std::cosf(2.0f * M_PI * (float)sample / (float)(winSize - 1));
}
void PitchCorrectionAudioProcessor::computeOLA_coeff()
{
    float windowSum = 0.0f;
    for (int sample = 0; sample < fftSize; ++sample)
        windowSum += anaWindow[sample];

    ola_coeff = 0.0f;
    if (windowSum != 0.0f)
        ola_coeff = 1.0f / 8.f / windowSum * (float)fftSize;
}
const float* PitchCorrectionAudioProcessor::FindClosest(const float* v, size_t N, float value)
{
    // Assume that v is sorted in descending order
    if (value >= v[0]) {
        return &v[0];
    }
    else if (value <= v[N - 1]) {
        return &v[N - 1];
    }
    else {
        const float* it = std::lower_bound(v, v + N, value, std::greater<float>());
        return std::abs(value - *it) < std::abs(value - *(it - 1)) ? it : it - 1;
    }
}
//TanhDist1.process(channelData, buffer.getNumSamples());
void PitchCorrectionAudioProcessor::processBlock(juce::AudioBuffer<float>& buffer, juce::MidiBuffer& midiMessages)
{
    juce::ScopedNoDenormals noDenormals;
    
    parameterManager.updateParameters();
    auto totalNumInputChannels = getTotalNumInputChannels();
    auto totalNumOutputChannels = getTotalNumOutputChannels();
    auto numSamples = buffer.getNumSamples();


    int currentInputBufferWritePosition;
    int currentOutputBufferWritePosition;
    int currentOutputBufferReadPosition;
    int currentSamplesSinceLastFFT;

    float shift = shiftRamp.getNextValue();
    if (autoTuneRamp.getTargetValue() == 1.f)
    {
        ratio = roundf(autoShift.getNextValue() * (float)hopSize) / (float)hopSize;
    }
    else
    {
        ratio = roundf(shift * (float)hopSize) / (float)hopSize;
    }
        
    int resampledLength = floorf((float)fftSize / ratio);
    juce::HeapBlock<float> resampledOutput(resampledLength, true);
    juce::HeapBlock<float> synthesisWindow(resampledLength, true);
    makeHannWindow(synthesisWindow, resampledLength);
    

    for (int channel = 0; channel < totalNumInputChannels; ++channel) 
    {
        float* channelData = buffer.getWritePointer(channel);

        currentInputBufferWritePosition = inputBufferWritePosition;
        currentOutputBufferWritePosition = outputBufferWritePosition;
        currentOutputBufferReadPosition = outputBufferReadPosition;
        currentSamplesSinceLastFFT = samplesSinceLastFFT;

        for (int sample = 0; sample < numSamples; ++sample) 
        {

            //======================================

            const float in = channelData[sample];
            float mix = blendRamp.getNextValue();
            if (enableRamp.getTargetValue() == 1.f)
                channelData[sample] = (outputBuffer.getSample(channel, currentOutputBufferReadPosition) * mix) + in*(1.f - mix);

            //======================================

            outputBuffer.setSample(channel, currentOutputBufferReadPosition, 0.0f);
            if (++currentOutputBufferReadPosition >= outputBufferLength)
                currentOutputBufferReadPosition = 0;

            //======================================

            inputBuffer.setSample(channel, currentInputBufferWritePosition, in);
            if (++currentInputBufferWritePosition >= inputBufferLength)
                currentInputBufferWritePosition = 0;

            //======================================

            if (++currentSamplesSinceLastFFT >= hopSize) {
                currentSamplesSinceLastFFT = 0;

                //======================================

                int inputBufferIndex = currentInputBufferWritePosition;
                for (int index = 0; index < fftSize; ++index) 
                {
                    fftDataTime[index].real((anaWindow[index]) * inputBuffer.getSample(channel, inputBufferIndex));
                    fftDataTime[index].imag(0.0f);
                    dataTimeToEstim[index] = fftDataTime[index].real();

                    if (++inputBufferIndex >= inputBufferLength)
                        inputBufferIndex = 0;
                }

                //======================================

                f0_instance.process(dataTimeToEstim, fftSize);

                currentF0 = f0_instance.get_f0();
                correctF0 = FindClosest(lookupChromatic, chromaticSize, currentF0);

                
                diffF0 = *correctF0 - currentF0;
                if ((std::abs(diffF0) < 6.5f) && !noobRamp.getTargetValue())
                    diffF0 = 0.f;

                autoShift.setTargetValue( powf(2.f, (smtHzRatio * diffF0) / 12.f) );// convert hz to semitones, then semitones to shift ratio

                fft.perform(fftDataTime, fftDataFreq, false);


                if (shiftRamp.isSmoothing())
                    needToResetPhases = true;
                if (shift == shiftRamp.getTargetValue() && needToResetPhases) 
                {
                    inputPhase.clear();
                    outputPhase.clear();
                    needToResetPhases = false;
                }

                for (unsigned int index = 0; index < fftSize; ++index)
                {
                    float magnitude = std::abs(fftDataFreq[index]);

                    float phase = std::arg(fftDataFreq[index]);

                    float phaseDeviation = phase - inputPhase.getSample(channel, index) - omega[index] * (float)hopSize;
                    float deltaPhi = omega[index] * hopSize + princArg(phaseDeviation);
                    float newPhase = princArg(outputPhase.getSample(channel, index) + deltaPhi * ratio);

                    inputPhase.setSample(channel, index, phase);
                    outputPhase.setSample(channel, index, newPhase);
                    fftDataFreq[index] = std::polar(magnitude, newPhase);
                }


                fft.perform(fftDataFreq, fftDataTime, true);

                for (int index = 0; index < resampledLength; ++index) 
                {
                    float x = (float)index * (float)fftSize / (float)resampledLength;
                    int ix = (int)floorf(x);
                    float dx = x - (float)ix;

                    float sample1 = fftDataTime[ix].real();
                    float sample2 = fftDataTime[(ix + 1) % fftSize].real();
                    resampledOutput[index] = sample1 + dx * (sample2 - sample1);
                    resampledOutput[index] *=(synthesisWindow[index]);
                }

                for (int index = 0; index < fftSize; ++index)
                {
                    fftDataTime[index].real((anaWindow[index]) * inputBuffer.getSample(channel, inputBufferIndex));
                    fftDataTime[index].imag(0.0f);
                    dataTimeToEstim[index] = fftDataTime[index].real();

                    if (++inputBufferIndex >= inputBufferLength)
                        inputBufferIndex = 0;
                }


                //======================================

                int outputBufferIndex = currentOutputBufferWritePosition;
                for (int index = 0; index < resampledLength; ++index) 
                {
                    float out = outputBuffer.getSample(channel, outputBufferIndex);
                    out += resampledOutput[index] / 2.99f; // ((float)hopSize / (float)fftSize); // 
                    //out += resampledOutput[index] * ola_coeff;
                    outputBuffer.setSample(channel, outputBufferIndex, out);

                    if (++outputBufferIndex >= outputBufferLength)
                        outputBufferIndex = 0;
                }

                //======================================

                currentOutputBufferWritePosition += hopSize;
                if (currentOutputBufferWritePosition >= outputBufferLength)
                    currentOutputBufferWritePosition = 0;
            }

            //======================================
        }
    }

    inputBufferWritePosition = currentInputBufferWritePosition;
    outputBufferWritePosition = currentOutputBufferWritePosition;
    outputBufferReadPosition = currentOutputBufferReadPosition;
    samplesSinceLastFFT = currentSamplesSinceLastFFT;

    //======================================
    outputGain.applyGain(buffer, buffer.getNumSamples());
    //enableRamp.applyGain(buffer, buffer.getNumSamples());
    // Apply the enableRamp as a switch between processed output and unprocessed input

    for (int channel = totalNumInputChannels; channel < totalNumOutputChannels; ++channel)
        buffer.clear(channel, 0, numSamples);
}

//==============================================================================
bool PitchCorrectionAudioProcessor::hasEditor() const
{
    return true; // (change this to false if you choose to not supply an editor)
}

juce::AudioProcessorEditor* PitchCorrectionAudioProcessor::createEditor()
{
    return new PitchCorrectionAudioProcessorEditor(*this);
}

//==============================================================================
void PitchCorrectionAudioProcessor::getStateInformation(juce::MemoryBlock& destData)
{
    // You should use this method to store your parameters in the memory block.
    // You could do that either as raw data, or use the XML or ValueTree classes
    // as intermediaries to make it easy to save and load complex data.
}

void PitchCorrectionAudioProcessor::setStateInformation(const void* data, int sizeInBytes)
{
    // You should use this method to restore your parameters from this memory block,
    // whose contents will have been created by the getStateInformation() call.
}

//==============================================================================
// This creates new instances of the plugin..
juce::AudioProcessor* JUCE_CALLTYPE createPluginFilter()
{
    return new PitchCorrectionAudioProcessor();
}
