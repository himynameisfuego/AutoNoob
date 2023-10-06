/*
  ==============================================================================

    f0Estim.h
    Created: 31 Aug 2023 10:29:43pm
    Author:  Tantep

  ==============================================================================
*/


#pragma once

#include <cmath>
#include <algorithm>

#define HUMAN_F0_LOWER (70.f)
#define HUMAN_F0_HIGHER (200.f)

class F0Estim // defining the class
{
public:

    F0Estim() {}; // The constructor, i.e. what happens when you instantiate an object

    ~F0Estim() {}; // The destuctor, i.e. what happens when the instance seizes to exist


    // Prepare to play should update the sample rate of the and the buffer size (if needed)
    // everytime it is needed
    void prepare(double sampleRate, unsigned int samplesPerBlock) {
        //round, set these for process
        float min_conveyable_f0 = sampleRate / samplesPerBlock;
        f0_min = std::max(min_conveyable_f0, HUMAN_F0_LOWER);
        f0_max = std::max(min_conveyable_f0, HUMAN_F0_HIGHER);
        min_lag = (unsigned int)lroundf(sampleRate / f0_max);
        max_lag = (unsigned int)lroundf(sampleRate / f0_min);
        this->sampleRate = sampleRate;
    };

    //processing for one channel
    void process(float* buffer, unsigned int numSamples) {
        int best_lag = 0;
        float highest_corr = -1.f;
        float buff_mean = mean(buffer, numSamples);

        //look in known range for maximum autocor val;
        for (int t = min_lag; t < max_lag; ++t) {
            float cov = 0; // Numerator
            float var = 0; // Denominator
            unsigned int seq_len = max_lag - min_lag;
            for (unsigned int i = 0; i < seq_len; ++i) {
                float xim = buffer[i] - buff_mean;
                cov += xim * (buffer[(i + t) % numSamples] - buff_mean);
                var += xim * xim;
            }
            float cur_corr = cov / var;
            if (cur_corr > highest_corr) {
                highest_corr = cur_corr;
                best_lag = t;
            }
        }
        //safety measure
        if (best_lag) {
            double new_val = sampleRate / best_lag;
            f0 = f0 ? (float)((1.f - smoothing_factor) * new_val + (smoothing_factor * f0)) :
                (float)new_val;
            return;
        }
        f0 = 0.f;
    };

    // Functions that will be used to be mapped to the parameters - so called "setters"

    void set_f0_min(unsigned int new_f0) { f0_min = new_f0; }
    void set_smoothing_factor(double new_factor)
    {
        double intpart;
        smoothing_factor = fabs(modf(new_factor, &intpart));
    }

    void set_f0_max(unsigned int new_f0) { f0_min = new_f0; }

    float get_f0() { return f0; }

private:
    double sampleRate = 48000; //can be changed later
    float f0 = 0.f;
    double f0_min = HUMAN_F0_LOWER; //lower limit of human voice
    double f0_max = HUMAN_F0_HIGHER; //higher limit of human voice
    unsigned int max_lag;
    unsigned int min_lag;
    double smoothing_factor = 0.5;

    float mean(float* buffer, unsigned int numSamples) {
        float sum = 0.f;
        for (unsigned int n = 0; n < numSamples; ++n) {
            sum += buffer[n];
        }
        return sum / numSamples;
    }
};