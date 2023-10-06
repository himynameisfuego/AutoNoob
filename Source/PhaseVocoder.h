/*
  ==============================================================================

    PhaseVocoder.h
    Created: 30 Aug 2023 10:46:55pm
    Author:  Tantep

  ==============================================================================
*/

#pragma once
#include <JuceHeader.h>
#define M_PI 3.14159265358979323846264338327950288
#include <cmath>

class PhaseVocoder // defining the class
{
public:
    PhaseVocoder() : fft(fftOrder), window(fftSize, juce::dsp::WindowingFunction<float>::hann)
    // The constructor, i.e. what happens when you instantiate an object
    {};

    ~PhaseVocoder() {}; // The destuctor, i.e. what happens when the instance seizes to exist

    // Prepare to play should update the sample rate of the and the buffer size (if needed) everytime it is needed
    void prepare(double sampleRate, int samplesPerBlock)
    {
    };

    //This function is where the actuall procssing takes place
    void process(float* buffer, unsigned int numSamples)
    {
        // For each sample ...
        for (int n = 0; n < numSamples; ++n)
        {
            // Apply a gain, pass the product through a tanh function and then multiply with the volume variable
            buffer[n] = volume * tanh(gain * buffer[n]) + M_PI;
        } 
    };

    // Functions that will be used to be mapped to the parameters - so called "setters"
    void setGain(float newGain) { gain = newGain; };
    void setVolume(float newVolume) { volume = newVolume; };


    static constexpr auto fftOrder = 9;             
    static constexpr auto fftSize  = 1 << fftOrder; // 2 ^ 9 = 512
    static constexpr int  hopSize = 512 * 0.25f;
    static constexpr float ola_coeff = 1.5f; // window scaling factor for size 512 with hop  128

    float PhaseVocoder::princArg(const float phase)
    {
        return fmod(phase + M_PI, -2.0f * M_PI) + M_PI;
        //if (phase >= 0.0f)
        //    return fmod(phase + M_PI, 2.0f * M_PI) - M_PI;
        //else
        //    return fmod(phase + M_PI, -2.0f * M_PI) + M_PI;
    }

private:

    // Here we define the gain and volume variables
    float gain = 1;
    float volume = 1;

    float princArg(const float phase); // principal argument function

    juce::dsp::FFT fft;                         // fft object
    juce::dsp::WindowingFunction<float> window; // gonna use Hann
    
    std::complex<float> fftDataTime; // time domain
    std::complex<float> fftDataFreq; // frequency domain
};