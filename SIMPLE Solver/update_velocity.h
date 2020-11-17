#include <iostream>
#include <vector>
#include <cmath>
#include <utility>
#include <tuple>
#include <cstdlib>
using namespace std;

void update_velocity(int imax, int jmax, vector<vector<float>> &u, vector<vector<float>> u_star, vector<vector<float>> &v, vector<vector<float>> v_star, vector<vector<float>> p_prime, vector<vector<float>> d_u, vector<vector<float>> d_v, float velocity)
{
    int i,j;
    //update interior nodes of u and v
    for(i=1;i<imax;i++)
    {
        for(j=1;j<jmax-1;j++)
        {
            u[i][j] = u_star[i][j] + d_u[i][j]*(p_prime[i-1][j] - p_prime[i][j]);
        }
    }

    for(i=1;i<imax-1;i++)
    {
        for(j=1;j<jmax;j++)
        {
            v[i][j] = v_star[i][j] + d_v[i][j]*(p_prime[i][j-1] - p_prime[i][j]);
        }
    }

    //update BCs
    for(j=0;j<jmax+1;j++)
    {
        v[0][j] = 0.0; //left wall
        v[imax-1][j]= 0.0; //right wall
    }

    for(i=0;i<imax;i++)
    {
        v[i][0]  = -v[i][1]; //bottom wall
        v[i][jmax]  = -v[i][jmax-1] ; //top wall
    }

    for(j=0;j<jmax;j++)
    {
        u[0][j] = -u[1][j]; //left wall
        u[imax][j] = -u[imax-1][j]; //right wall
    }

    for(i=0;i<imax+1;i++)
    {
        u[i][0] = 0.0; //bottom wall
        u[i][jmax-1] = velocity; //top wall
    }
}


