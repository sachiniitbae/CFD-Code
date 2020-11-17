#include <iostream>
#include <vector>
#include <cmath>
#include <utility>
#include <tuple>
#include <cstdlib>
using namespace std;

float AU(float F, float D)
{
    return fmax(0, pow((1 - 0.1 * abs(F/D)),5));
}

void u_momentum_new(int imax,int jmax,float dx,float dy,float rho, float mu, vector<vector<float>> u, vector<vector<float>> &u_star, vector<vector<float>> &d_u, vector<vector<float>> v, vector<vector<float>> p, float velocity, float alpha)
{
    int i,j;
    float Fe,Fw,Fn,Fs,aE,aW,aN,aS,aP,pressure_term;
    float De,Dw,Dn,Ds;

    De  = mu*dy / dx;  //convective coefficients
    Dw  = mu*dy / dx;
    Dn  = mu*dx / dy;
    Ds  = mu*dx / dy;

    for(i=1;i<imax;i++)
    {
        for(j=1;j<jmax-1;j++)
        {
            Fe  = 0.5*rho*dy*(u[i+1][j]+u[i][j]);
            Fw  = 0.5*rho*dy*(u[i-1][j]+u[i][j]);
            Fn  = 0.5*rho*dx*(v[i][j+1]+v[i-1][j+1]);
            Fs  = 0.5*rho*dx*(v[i][j]+v[i-1][j]);

            aE = De * AU(Fe,De) + fmax(-Fe,0);
            aW = Dw * AU(Fw,Dw) + fmax(Fw,0);
            aN = Dn * AU(Fn,Dn) + fmax(-Fn,0);
            aS = Ds * AU(Fs,Ds) + fmax(Fs,0);
            aP = aE + aW + aN + aS + (Fe-Fw) + (Fn-Fs);

            pressure_term = (p[i-1][j] - p[i][j]) * dy;

            u_star[i][j] = (alpha/aP) * ( (aE*u[i+1][j] + aW*u[i-1][j] + aN*u[i][j+1] + aS*u[i][j-1]) + pressure_term ) + (1-alpha)*u[i][j];

            d_u[i][j] = alpha * dy / aP;   //refer to Versteeg CFD book
        }
    }

    //set d_u for top and bottom BCs
    //they will be later used by the pressure correction equation
    //they should not be zero, or BCs of pressure correction will get messed up

    j = 0; //bottom
    for(i=1;i<imax;i++)
    {
        Fe  = 0.5*rho*dy*(u[i+1][j]+u[i][j]);
        Fw  = 0.5*rho*dy*(u[i-1][j]+u[i][j]);
        Fn  = 0.5*rho*dx*(v[i][j+1]+v[i-1][j+1]);
        Fs  = 0;

        aE = De * AU(Fe,De) + fmax(-Fe,0);
        aW = Dw * AU(Fw,Dw) + fmax(Fw,0);
        aN = Dn * AU(Fn,Dn) + fmax(-Fn,0);
        aS = 0.0;
        aP = aE + aW + aN + aS + (Fe-Fw) + (Fn-Fs);
        d_u[i][j] = alpha * dy / aP;
    }

    j = jmax-1; //top
    for(i=1;i<imax;i++)
    {
        Fe  = 0.5*rho*dy*(u[i+1][j]+u[i][j]);
        Fw  = 0.5*rho*dy*(u[i-1][j]+u[i][j]);
        Fn  = 0;
        Fs  = 0.5*rho*dx*(v[i][j]+v[i-1][j]);

        aE = De * AU(Fe,De) + fmax(-Fe,0);
        aW = Dw * AU(Fw,Dw) + fmax(Fw,0);
        aN = 0.0;
        aS = Ds * AU(Fs,Ds) + fmax(Fs,0);
        aP = aE + aW + aN + aS + (Fe-Fw) + (Fn-Fs);
        d_u[i][j] = alpha * dy / aP;
    }

    //Apply BCs
    for(j=0;j<jmax;j++)
    {
        u_star[0][j] = -u_star[1][j];           //left wall
        u_star[imax][j] = -u_star[imax-1][j];   //right wall
    }
    for(i=0;i<imax+1;i++)
    {
        u_star[i][0] = 0.0;                      //bottom wall
        u_star[i][jmax-1] = velocity;            //top wall
    }

}
