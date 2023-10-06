/*
  ==============================================================================

    TanhDist.h
    Created: 30 Aug 2023 5:04:29pm
    Author:  Tantep

  ==============================================================================
*/

#pragma once
// We need the math library to use the tanh function
#include <math.h>
class TanhDist // defining the class
{
public:
    TanhDist(){}; // The constructor, i.e. what happens when you instantiate an object
    ~TanhDist(){}; // The destuctor, i.e. what happens when the instance seizes to exist
    // Prepare to play should update the sample rate of the and the buffer size (if needed) everytime it is needed
    void prepare(double sampleRate, int samplesPerBlock)
    {
    };
    //This function is where the actuall procssing takes place
    void process(float* buffer, unsigned int numSamples)
    {
        // For each sample ...
        for (int n = 0 ; n < numSamples; ++n)
        {
            // Apply a gain, pass the product through a tanh function and then multiply with the volume variable
            buffer[n] = volume*tanh(gain*buffer[n]);
        }
    };
// Functions that will be used to be mapped to the parameters - so called "setters"
void setGain(float newGain){ gain = newGain;};
void setVolume(float newVolume){ volume = newVolume;};

private:
    // Here we define the gain and volume variables
    float gain = 1;
    float volume = 1;
};