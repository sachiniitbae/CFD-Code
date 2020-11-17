#include <iostream>
#include <vector>
#include <cmath>
#include <utility>
#include <tuple>
#include <cstdlib>
#include "get_diagonals.h"
#include "PDMA.h"
using namespace std;

vector<vector<float>> pressure_correct(int imax, int jmax, vector<float> &rhsp, vector<vector<float>> &Ap, vector<vector<float>> p, vector<vector<float>> &p_prime, float alpha)
{
    int i,j,z = 0;
    vector<vector<float>> pressure = p;
    vector <float> p_prime_interior;
    vector<vector <float>> diagonals(5, vector<float> (imax*jmax));
    get_diagonals(imax,jmax,Ap,diagonals);

    p_prime_interior = PDMA(Ap,rhsp,diagonals);
    //convert pressure correction in to a matrix
    //update pressure values

    for(j=0;j<jmax;j++)
    {
        for(i=0;i<imax;i++)
        {
            p_prime[i][j] = p_prime_interior[z];
            z = z + 1;
            pressure[i][j] = p[i][j] + alpha*p_prime[i][j];
        }
    }
    pressure[0][0] = 0;

    return pressure;

}
