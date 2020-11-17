#include <iostream>
#include <vector>
#include <cmath>
#include <utility>
#include <tuple>
#include <cstdlib>
using namespace std;

void get_rhs(int imax, int jmax, float dx, float dy, float rho, vector<vector<float>> u_star, vector<vector<float>> v_star, vector<float> &rhsp)
{
    int i,j;
    int stride = jmax;
    int position;


    // RHS is the same for all nodes except the p_prime(1,1)
    // Because p(1,1) is set to be zero, it has no pressure correction

    for(j=0;j<jmax;j++)
    {
        for(i=0;i<imax;i++)
        {
            position = i+(j-0)*stride;
            rhsp[position] = rho * (u_star[i][j]*dy - u_star[i+1][j]*dy  + v_star[i][j]*dx - v_star[i][j+1]*dx);
        }
    }

    // modify for p_prime(1,1)
    rhsp[0] = 0.0;
}
