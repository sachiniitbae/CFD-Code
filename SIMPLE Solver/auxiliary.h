#include <iostream>
#include <vector>
using namespace std;

float A(float F, float D)
{
    return fmax(0, pow((1-0.1 * abs(F/D)),5) );
}

