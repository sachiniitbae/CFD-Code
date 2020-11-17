#include <iostream>
#include <vector>
#include <cmath>
#include <utility>
#include <tuple>
#include <cstdlib>
using namespace std;

float AV(float F, float D)
{
    return fmax(0, pow((1 - 0.1 * abs(F/D)),5) );
}

void v_momentum_new(int imax,int jmax,float dx,float dy,float rho, float mu, vector<vector<float>> u, vector<vector<float>> v, vector<vector<float>> &v_star, vector<vector<float>> &d_v, vector<vector<float>> p, vector<vector<float>> T, float alpha)
{
    int i,j;
    float Fe,Fw,Fn,Fs,aE,aW,aN,aS,aP,pressure_term;
    float Boussinesq_source = 0.0;
    float De  = mu*dy / dx;  //convective coefficients
    float Dw  = mu*dy / dx;
    float Dn  = mu*dx / dy;
    float Ds  = mu*dx / dy;

    for(i=1;i<imax-1;i++)
    {
        for(j=1;j<jmax;j++)
        {
            Fe  = 0.5*rho*dy*(u[i+1][j]+u[i+1][j-1]);
            Fw  = 0.5*rho*dy*(u[i][j]+u[i][j-1]);
            Fn  = 0.5*rho*dx*(v[i][j]+v[i][j+1]);
            Fs  = 0.5*rho*dx*(v[i][j-1]+v[i][j]);

            aE = De * AV(Fe,De) + fmax(-Fe,0);
            aW = Dw * AV(Fw,Dw) + fmax(Fw,0);
            aN = Dn * AV(Fn,Dn) + fmax(-Fn,0);
            aS = Ds * AV(Fs,Ds) + fmax(Fs,0);
            aP = aE + aW + aN + aS + (Fe-Fw) + (Fn-Fs);

            pressure_term = (p[i][j-1] - p[i][j]) * dx;
            Boussinesq_source = 0.5*(T[i][j-1]+T[i][j]) * dx*dy ;

            v_star[i][j] = alpha/aP * ( (aE*v[i+1][j] + aW*v[i-1][j] + aN*v[i][j+1] + aS*v[i][j-1]) + pressure_term + Boussinesq_source) + (1-alpha)*v[i][j];

            d_v[i][j] = alpha * dx / aP;   //refer to Versteeg CFD book
        }
    }

    //set d_v for left and right BCs
    //they will be later used by the pressure correction equation
    //they should not be zero, or BCs of pressure correction will get messed up
    //Apply BCs

    i = 0; //Left
    for(j=1;j<jmax;j++)
    {
        Fe  = 0.5*rho*dy*(u[i+1][j]+u[i+1][j-1]);
        Fw  = 0.0;
        Fn  = 0.5*rho*dx*(v[i][j]+v[i][j+1]);
        Fs  = 0.5*rho*dx*(v[i][j-1]+v[i][j]);

        aE = De * AV(Fe,De) + fmax(-Fe,0);
        aW = 0.0;
        aN = Dn * AV(Fn,Dn) + fmax(-Fn,0);
        aS = Ds * AV(Fs,Ds) + fmax(Fs,0);
        aP = aE + aW + aN + aS + (Fe-Fw) + (Fn-Fs);
        d_v[i][j] = alpha * dx / aP;
    }

    i = imax-1; //top
    for(j=1;j<imax;j++)
    {
        Fe  = 0.0;
        Fw  = 0.5*rho*dy*(u[i][j]+u[i][j-1]);
        Fn  = 0.5*rho*dx*(v[i][j]+v[i][j+1]);
        Fs  = 0.5*rho*dx*(v[i][j-1]+v[i][j]);

        aE = 0.0;
        aW = Dw * AV(Fw,Dw) + fmax(Fw,0);
        aN = Dn * AV(Fn,Dn) + fmax(-Fn,0);
        aS = Ds * AV(Fs,Ds) + fmax(Fs,0);
        aP = aE + aW + aN + aS + (Fe-Fw) + (Fn-Fs);
        d_v[i][j] = alpha * dx / aP;
    }

    //Apply BCs
    for(j=0;j<jmax+1;j++)
    {
        v_star[0][j] = 0.0;      //left wall
        v_star[imax-1][j] = 0.0;  //right wall
    }
    for(i=0;i<imax;i++)
    {
        v_star[i][0] = - v_star[i][1];                     //bottom wall
        v_star[i][jmax] = - v_star[i][jmax-1] ;          //top wall
    }

}
