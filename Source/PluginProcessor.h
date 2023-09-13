/*
  ==============================================================================

    This file contains the basic framework code for a JUCE plugin processor.

  ==============================================================================
*/

#pragma once

#include <JuceHeader.h>
#include "TanhDist.h"
#include "Ramp.h"
#include "F0Estim.h"
namespace Param
{
    namespace ID
    {
        static const juce::String Bypass{ "bypass" };
        static const juce::String Semitones{ "semitones" };
        static const juce::String Autotune{ "autotune" };
        static const juce::String Noob{ "noob" };
        static const juce::String Blend{ "blend" };
        static const juce::String Volume{ "volume" };
    }
    namespace Name
    {
        static const juce::String Bypass{ "Bypass" };
        static const juce::String Semitones{ "Semitones" };
        static const juce::String Autotune{ "Autotune" };
        static const juce::String Noob{ "Noob" };
        static const juce::String Blend{ "Blend" };
        static const juce::String Volume{ "Volume" };
    }
}

//==============================================================================
/**
*/
class PitchCorrectionAudioProcessor : public juce::AudioProcessor
    #if JucePlugin_Enable_ARA
        , public juce::AudioProcessorARAExtension
    #endif
{
public:
    //==============================================================================
    PitchCorrectionAudioProcessor();
    ~PitchCorrectionAudioProcessor() override;

    //==============================================================================
    void prepareToPlay(double sampleRate, int samplesPerBlock) override;
    void releaseResources() override;

#ifndef JucePlugin_PreferredChannelConfigurations
    bool isBusesLayoutSupported(const BusesLayout& layouts) const override;
#endif

    void processBlock(juce::AudioBuffer<float>&, juce::MidiBuffer&) override;

    mrta::ParameterManager& getParameterManager() { return parameterManager; }
    //==============================================================================
    juce::AudioProcessorEditor* createEditor() override;
    bool hasEditor() const override;

    //==============================================================================
    const juce::String getName() const override;

    bool acceptsMidi() const override;
    bool producesMidi() const override;
    bool isMidiEffect() const override;
    double getTailLengthSeconds() const override;

    //==============================================================================
    int getNumPrograms() override;
    int getCurrentProgram() override;
    void setCurrentProgram(int index) override;
    const juce::String getProgramName(int index) override;
    void changeProgramName(int index, const juce::String& newName) override;

    //==============================================================================
    void getStateInformation(juce::MemoryBlock& destData) override;
    void setStateInformation(const void* data, int sizeInBytes) override;
    
    int fftOrder = 10;
    int fftSize = 1 << fftOrder; // 2 ^ 9 = 1024
    int hopSize = 1024 * 0.125f;
    

private:


    mrta::ParameterManager parameterManager;
    //TanhDist TanhDist1;
    
    
    void computeOLA_coeff();
    float ola_coeff; // window scaling factor for size 1024 with hop 128

    float analysisMagnitude[1024];
    float analysisFrequency[1024];
    float synthesisMagnitude[1024];
    float synthesisFrequency[1024];
    float lastInputPhase[1024];
    float lastOutputPhase[1024];

    float princArg(const float phase);          // principal argument function
    void makeHannWindow(float* window, const int winSize);
    

    juce::HeapBlock<juce::dsp::Complex<float>> fftDataTime, fftDataTime_y;            // time domain
    juce::HeapBlock<juce::dsp::Complex<float>> fftDataFreq;            // frequency domain
    //juce::HeapBlock<juce::dsp::Complex<float>> cepstrum_x,cepstrum_X, cepstrum_y, cepstrum_Y;
    unsigned int numLPFCepstrum = 3;
    float dataTimeToEstim[1024];
    float currentF0 = 0.f;
    const float* correctF0;
    float diffF0;

    juce::dsp::FFT fft;     // fft object
    float anaWindow[1024]; // gonna use Hann

    int inputBufferLength;
    int inputBufferWritePosition;
    juce::AudioSampleBuffer inputBuffer;

    int outputBufferLength;
    int outputBufferWritePosition;
    int outputBufferReadPosition;
    juce::AudioSampleBuffer outputBuffer;

    int samplesSinceLastFFT;


    //======================================

    juce::HeapBlock<float> omega;
    juce::AudioSampleBuffer inputPhase;
    juce::AudioSampleBuffer outputPhase;

    bool needToResetPhases = true;

    //float shift = 1.f;
    juce::SmoothedValue<float> outputGain {1.f};
    juce::SmoothedValue<float> shiftRamp{ 1.f };
    juce::SmoothedValue<float> enableRamp {1.f};
    juce::SmoothedValue<float> autoTuneRamp{ 0.f };
    juce::SmoothedValue<float> autoShift{ 1.f };
    juce::SmoothedValue<float> blendRamp{ 1.f };
    juce::SmoothedValue<float> noobRamp{ 0.f };
    
    F0Estim f0_instance;

    
    const size_t chromaticSize = 26;
    const float* FindClosest(const float* v, size_t N, float value);
    float lookupChromatic[26] =
    {
    293.66f, 277.18f, 261.63f, 246.94f, 233.08f, 220.0f, 207.65f, 196.0f,
    185.0f, 174.61f, 164.81f, 155.56f, 146.83f, 138.59f, 130.81f, 123.47f,
    116.54f, 110.0f, 103.83f, 98.0f, 92.5f, 87.31f, 82.41f, 77.78f, 73.42f ,0.f,
    };
    float smtHzRatio = 1.0594630943592953f; // Hz per 1 semitone
    //float autoShift = 0.f;
    float ratio = 1.f; // actual ratio for pitch shifting
    //==============================================================================
    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(PitchCorrectionAudioProcessor)
};
