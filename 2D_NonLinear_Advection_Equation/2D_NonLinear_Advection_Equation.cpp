#include <iostream>
#include <vector>
using namespace std;
#define nx 81
#define ny 81
#define nt 100

void Scheme(float u[][ny], float un[][ny],float v[][ny], float vn[][ny], float dt, float dx, float dy, float c)
{
    FILE *fp;
    fp = fopen("2D_NonLinear_Advection_final_Condition.dat","w");
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
                vn[i][j]= v[i][j];
            }
        }

        fprintf(fp,"VARIABLES = \"x\"\t\"y\"\t\"u\"\t\"v\"\n");
        fprintf(fp,"ZONE T=\"%d\", I=%d, J=%d, F = point\n",k,nx,ny);

        for ( i = 1; i < nx-1; i++)
        {
            for (j = 1; j < ny-1; j++)
            {
                u[i][j] = (un[i][j] - (un[i][j] * dt / dx * (un[i][j] - un[i-1][j])) - (vn[i][j]* dt / dy * (un[i][j] - un[j-1][i])));
                u[0][j] = 1.0;
                u[nx][j] = 1.0;
                u[i][0] = 1.0;
                u[i][ny] = 1.0;

                v[i][j] = (vn[i][j] - (un[i][j] * dt / dx * (vn[i][j] - vn[i-1][j])) - (vn[i][j]* dt / dy * (vn[i][j] - vn[j-1][i])));
                v[0][j] = 1.0;
                v[nx][j] = 1.0;
                v[i][0] = 1.0;
                v[i][ny] = 1.0;
            }
        }

        for (i=0 ; i < nx ; i++)
            {
                for (j=0 ; j < ny ; j++)
                {
                    fprintf(fp,"%f\t %f\t %f\t %f\n",x[i],y[j],u[i][j],v[i][j]);
                }
            }
            fprintf(fp,"\n\n");
    }
/*
    fprintf(fp,"VARIABLES = \"x\"\t,\"y\"\t,\"u\"\t,\"v\"\n");
    for (i=0 ; i < nx ; i++)
    {
        for (j=0 ; j < ny ; j++)
        {
            fprintf(fp,"%0.2f\t %0.2f\t %f\t %f\n",x[i],y[j],u[i][j],v[i][j]);
        }
    }*/
}

int main()
{
    int i,j;
    float c = 1.0;
    float dx = 2.0 / (nx - 1.0);
    float dy = 2.0 / (ny - 1.0);
    float sigma = 0.2;
    float dt = sigma * dx;
    FILE *fp;
    fp = fopen("2D_NonLinear_Advection_Initial_Condition.dat","w");
    float x[nx],y[ny],u[nx][ny], un[nx][ny],v[nx][ny], vn[nx][ny];

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
            v[i][j] = 1.0;
            un[i][j] = 1.0;
            vn[i][j] = 1.0;
        }
    }

    for (i = int(0.5/dx) ; i <= int(1.0/dx+1.0) ; i++ )
    {
        for(j = int(0.5/dy) ; j <= int(1.0/dy+1.0) ; j++)
        {
            u[i][j] = 2.0;
            v[i][j] = 2.0;
        }
    }

    fprintf(fp,"VARIABLES = \"x\"\t,\"y\"\t,\"u\"\t,\"v\"\n");
    for (i = 0; i < nx; i++)
    {
      for (j = 0; j < ny; j++)
      {
          fprintf(fp,"%0.2f\t %0.2f\t %f\t %f\n",x[i],y[j],u[i][j],v[i][j]);
      }
   }
   Scheme(u,un,v,vn,dt,dx,dy,c);
}
