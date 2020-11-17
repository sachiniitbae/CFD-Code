#include <iostream>
#include <vector>

using namespace std;

void check_divergence(int imax,int jmax, float dx, float dy, vector<vector<float>> u, vector<vector<float>> v, vector<vector<float>> &div)
{
    int i,j;
    for(i=0;i<imax;i++)
    {
        for(j=0;j<jmax;j++)
        {
            div[i][j] = (1/dx)*(u[i][j]-u[i+1][j]) + (1/dy)*(v[i][j]-v[i][j+1]);
        }
    }

}
