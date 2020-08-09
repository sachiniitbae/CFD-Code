#include <stdio.h>
#include <math.h>
#define M_PI 3.14


float funct(float t, float x, float nu)
{
    float value,phi,phiprime;
    phi = exp(-pow((-4*t + x - 2*M_PI),2)/(4*nu*(t + 1))) + exp(-pow((-4*t + x),2)/(4*nu*(t + 1)));
    phiprime = -(-8*t + 2*x)*exp(-pow((-4*t + x),2)/(4*nu*(t + 1)))/(4*nu*(t + 1)) - (-8*t + 2*x - 4*M_PI)*exp(-pow((-4*t + x - 2*M_PI),2)/(4*nu*(t + 1)))/(4*nu*(t + 1));
    value = -2 * nu * (phiprime / phi) + 4 ;
    return value;
}

void Scheme(float domain[], float u[], float un[],float uanalytical[], int grid, int nt,float nu,float sigma)
{
    FILE *fp;
    fp = fopen("1D_Burgers.csv","a");
    int i,j;
    for (j = 0; j < nt; j++)
    {
       for(i=0;i<grid;i++)
        {
            un[i]= u[i];
        }

        for ( i = 1; i<grid-1; i++)
        {
            u[i] = un[i] - un[i] * nu *(un[i] - un[i-1]) + sigma *(un[i+1] - 2 * un[i] + un[i-1]) ;
            u[0] = un[0] - un[0] * nu* (un[0] - un[-2]) + sigma*(un[1] - 2 * un[0] + un[-2]);
            u[-1] = u[0];
        }
    }
    fprintf(fp,"At time step %d\n",nt);
    for (i=0;i<grid;i++)
    {
        fprintf(fp,"%f,%f,%f\n",domain[i],u[i],uanalytical[i]);
    }
}


int main()
{
    FILE *fp;
    fp = fopen("1D_Burgers.csv","w");
    int grid = 101,nt,i;
    float dx,dt,nu,t=0.0,sigma;
    float domain[grid],u[grid], un[grid];
    float uanalytical[grid];
    dx = (2.0* M_PI)/(grid-1);
    nu = 0.07;
    nt = 100;
    dt = dx*nu;
    sigma = (dt *nu)/(dx*dx);


    for (i=0; i < grid; i++)
    {
        domain[i] = i*dx;
    }


    for (i = 0; i<grid; i++)
    {
        u[i] = funct(t,domain[i],nu);
        uanalytical[i] = funct((nt*dt),domain[i],nu);
        un[i] = 1.0;
    }
    fprintf(fp,"Initial Condition\n");

    for (i=0;i<grid;i++)
    {
        fprintf(fp,"%f,%f\n",domain[i],u[i]);
    }
    Scheme(domain,u,un,uanalytical,grid,nt,nu,sigma);
}
