#include <iostream>
#include <vector>
#include <cmath>
#include <utility>
#include <tuple>
#include <cstdlib>
using namespace std;

float AT(float F, float Diffusion)
{
    return fmax(0, pow((1 - 0.1 * abs(F/Diffusion)),5));
}

void energy_equation(int imax,int jmax,float dx,float dy, vector<vector<float>> u, vector<vector<float>> v, vector<vector<float>> &T, vector<vector<float>> Told, float D, float alpha,float T_top, float T_bottom)
{
    int i,j;
    float Fe,Fw,Fn,Fs,aE,aW,aN,aS,aP,pressure_term;
    float De,Dw,Dn,Ds,U,V;

    De  = D*dy / dx;  //convective coefficients
    Dw  = D*dy / dx;
    Dn  = D*dx / dy;
    Ds  = D*dx / dy;

    for(i=1;i<imax-1;i++)
    {
        for(j=1;j<jmax-1;j++)
        {
            Fe  = dy*(u[i+1][j]);
            Fw  = dy*(u[i][j]);
            Fn  = dx*(v[i][j+1]);
            Fs  = dx*(v[i][j]);

            U = 0.5*(u[i+1][j] + u[i][j]);
            V = 0.5*(v[i][j+1] + v[i][j]);

            aE = De * AT(Fe,De) + fmax(-Fe,0);
            aW = Dw * AT(Fw,Dw) + fmax(Fw,0);
            aN = Dn * AT(Fn,Dn) + fmax(-Fn,0);
            aS = Ds * AT(Fs,Ds) + fmax(Fs,0);
            aP = aE + aW + aN + aS + (Fe-Fw) + (Fn-Fs);

            T[i][j] = (alpha/aP) * ( (aE*Told[i+1][j] + aW*Told[i-1][j] + aN*Told[i][j+1] + aS*Told[i][j-1])) + (1-alpha)*Told[i][j];

        }
    }

    //Apply BCs
    for(j=0;j<jmax;j++)
    {
        T[0][j] = T[1][j];           //left wall
        T[imax-1][j] = T[imax-2][j];   //right wall
    }
    for(i=0;i<imax;i++)
    {
        T[i][0] = T_bottom;               //bottom wall
        T[i][jmax-1] = T_top;            //top wall
    }

}
