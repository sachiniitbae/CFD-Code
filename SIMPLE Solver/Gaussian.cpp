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

// Driver program to test above function
int main()
{
    int imax = 33;
    int jmax = 33;
    vector<vector<float>> GKernel(imax,vector<float>(jmax,0));
    float dx = 1.0/(imax-1);
    float dy = 1.0/(jmax-1);
    cout << dx << endl;
    vector<float> x = linspace_interval(0.0, 1.0, dx);;
    vector<float> y = linspace_interval(0.0, 1.0, dy);

    print_vector_1d(x);
    FilterCreation(imax,jmax,x,y,GKernel);
    print_vector_2d(GKernel);
}
