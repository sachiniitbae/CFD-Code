#include <cstdio>
#include <iostream>
#include <algorithm>
#include <array>
using namespace std;

void Scheme(float domain[], float u[], float un[], float grid, int nt, float sigma)
{
    //Discretize the second-order derivative with a
    //Central Difference scheme: a combination of
    //Forward Difference and Backward Difference of the first derivative
    FILE *fp;
    fp = fopen("1D_Diffusion.csv","w");
    int i,j;
    for (j = 0; j < nt; j++)
    {
       for(i=0;i<grid;i++)
        {
            un[i]= u[i];
        }

        for ( i = 1; i<grid-1; i++)
        {
            u[i] = un[i] + sigma * (un[i+1] - 2*un[i] + un[i-1]);
        }
    }
    for (i=0;i<grid;i++)
    {
        fprintf(fp,"%f,%f\n",domain[i],u[i]);
    }
}


int main()
{
    int grid = 41,nt,i;
    float dx,dt,nu,sigma;
    float domain[grid],u[grid], un[grid];
    dx = 2.0/(grid-1);
    nu = 0.3;
    nt = 20;
    sigma = 0.2;
    dt = (sigma * dx*dx )/ nu;

    for (i=0; i < grid; i++)
    {
        domain[i] = i*dx;
    }

    fill_n(u,grid,1.0);
    fill_n(un,grid,1.0);
    for (i= int(0.5/dx) ; i <= int(1.0/dx+1.0) ; i++ )
    {
        u[i] = 2.0;
    }

    Scheme(domain,u,un,grid,nt,sigma);
}
