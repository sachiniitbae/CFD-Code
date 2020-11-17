#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>
#include "print_vector.h"
#include "linspace.h"
using namespace std;

// Function to create Gaussian filter
void FilterCreation(int imax, int jmax, vector<float> x, vector<float> y, vector<vector<float>> &GKernel)
{
    // Initializing standard deviation to 1.0
    float sigma = 1.0;
    float r, s = 2.0 * sigma * sigma;

    // sum is for normalization
    float sum = 0.0;

    // generating 5x5 kernel
    for (int i = 0; i < imax; i++)
    {
        for (int j = 0; j < jmax; j++)
        {
            r = sqrt(x[i] * x[i] + y[j] * y[j]);
            GKernel[i][j] = (exp(-(r * r) / s)) / (M_PI * s);
            sum += GKernel[i][j];
        }
    }

    // Normalizing the Kernel
    for (int i = 0; i < imax; ++i)
        for (int j = 0; j < jmax; ++j)
            GKernel[i][j] /= sum;
}

