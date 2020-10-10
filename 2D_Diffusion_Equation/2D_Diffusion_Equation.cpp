#include <iostream>
#include <vector>
#include <math.h>
using namespace std;
#define nx 81
#define ny 81
#define nt 100

void Diffusion(float u[][ny], float un[][ny],float nu, float dt, float dx, float dy)
{
    FILE *fp;
    fp = fopen("2D_Diffusion_final_Condition.dat","w");
    int i,j,k;
    float y[ny],x[nx];
    for (j = 0; j < ny;j++)
    {
        y[j] = j*dy;
        x[j] = j*dx;
    }

    for (k = 0;k < nt; k++)
    {
       for(i=0;i<nx;i++)
        {
            for(j=0;j<ny;j++)
            {
                un[i][j]= u[i][j];
            }
        }

        fprintf(fp,"VARIABLES = \"x\"\t\"y\"\t\"u\"\n");
        fprintf(fp,"ZONE T=\"%d\", I=%d, J=%d, F = point\n",k,nx,ny);

        for ( i = 1; i < nx-1; i++)
        {
            for (j = 1; j < ny-1; j++)
            {

                u[i][j] = (un[i][j] + nu * dt / pow(dx,2) * (un[i+1][j] - 2 * un[i][j] + un[i-1][j]) + nu * dt / pow(dy,2) * (un[i][j+1] - 2 * un[i][j] + un[i][j-1])) ;
                u[0][j] = 1.0;
                u[nx][j] = 1.0;
                u[i][0] = 1.0;
                u[i][ny] = 1.0;
            }
        }

        for (i=0 ; i < nx ; i++)
            {
                for (j=0 ; j < ny ; j++)
                {
                    fprintf(fp,"%f\t %f\t %f\n",x[i],y[j],u[i][j]);
                }
            }
            fprintf(fp,"\n\n");
    }
/*
    fprintf(fp,"VARIABLES = \"x\"\t,\"y\"\t,\"u\"\n");
    for (i=0 ; i < nx ; i++)
    {
        for (j=0 ; j < ny ; j++)
        {
            fprintf(fp,"%0.2f\t %0.2f\t %f\n",x[i],y[j],u[i][j]);
        }
    }*/
}

int main()
{
    int i,j;
    float dx = 2.0 / (nx - 1.0);
    float dy = 2.0 / (ny - 1.0);
    float sigma = 0.25;
    float nu = 0.5;
    float dt = sigma * dx * dy / nu;
    FILE *fp;
    fp = fopen("2D_Diffusion_Initial_Condition.dat","w");
    float x[nx],y[ny],u[nx][ny], un[nx][ny];

    for (i = 0; i < nx ;i++)
    {
        x[i] = i*dx;
        y[i] = i*dy;
    }

    for (i = 0; i<nx ;i++)
    {
        for(j=0;j<ny;j++)
        {
            u[i][j] = 1.0;
            un[i][j] = 1.0;
        }
    }

    for (i = int(0.5/dx) ; i <= int(1.0/dx+1.0) ; i++ )
    {
        for(j = int(0.5/dy) ; j <= int(1.0/dy+1.0) ; j++)
        {
            u[i][j] = 2.0;
        }
    }

    fprintf(fp,"VARIABLES = \"x\"\t,\"y\"\t,\"u\"\n");
    for (i = 0; i < nx; i++)
    {
      for (j = 0; j < ny; j++)
      {
          fprintf(fp,"%0.2f\t %0.2f\t %f\n",x[i],y[j],u[i][j]);
      }
   }
   Diffusion(u,un,nu,dt,dx,dy);
}
