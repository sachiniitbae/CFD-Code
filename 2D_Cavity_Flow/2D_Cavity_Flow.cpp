#include <iostream>
#include <vector>
#include <math.h>
using namespace std;
#define nx 41
#define ny 41
#define nt 100
#define nit 50

float build_up_b(float b[][ny],float rho,float dt,float u[][ny],float v[][ny], float dx, float dy)
{
    int i,j;
    b[i][j] = (rho * (1 / dt *
                ((u[i+1][j] - u[i-1][j]) /
                (2 * dx) + (v[i][j+1] - v[i][j-1]) / (2 * dy)) -
                pow(((u[i+1][j] - u[i-1][j]) / (2 * dx)),2) -
                2 * ((u[i][j+1] - u[i][j-1]) / (2 * dy) *
                (v[i+1][j] - v[i-1][j]) / (2 * dx))-
                (pow(((v[i][j+1] - v[i][j-1]) / (2 * dy)),2)))) ;

    return b[i][j];

}

float pressure_poisson(float p[][ny],float dx, float dy, float b[][ny])
{
    int i,j,k;
    static float pn[nx][ny];
    for (i = 0; i<nx ;i++)
    {
        for(j=0;j<ny;j++)
        {
            pn[i][j] = 0.0;
        }
    }

    for(k=0;k<nit ; k++)
    {
        for(i=0;i<nx;i++)
        {
            for(j=0;j<ny;j++)
            {
                pn[i][j]= p[i][j];
            }
        }

        p[i][j] = (((pn[i+1][j] + pn[i-1][j]) * dy*dy +
                          (pn[i][j+1] + pn[i][j-1]) * dx*dx) /
                          (2 * (dx*dx + dy*dy)) -
                          dx*dx * dy*dy / (2 * (dx*dx + dy*dy)) *
                          b[i][j]);

        p[nx][j] = p[nx-1][j]; // dp/dx = 0 at x = 2
        p[i][0] = p[i][1];   // dp/dy = 0 at y = 0
        p[0][j] = p[1][j];   // dp/dx = 0 at x = 0
        p[i][ny] = 0;    // p = 0 at y = 2
    }

    return p[i][j];


}


void Cavity_flow(float u[][ny],float v[][ny],float p[][ny],float rho, float nu, float dt, float dx, float dy)
{
    FILE *fp;
    fp = fopen("2D_Cavity_Flow_final_Solution.dat","w");
    int i,j,k;
    static float y[ny],x[nx],b[nx][ny];
    static float un[nx][ny],vn[nx][ny];


    for (j = 0; j < ny;j++)
    {
        y[j] = j*dy;
        x[j] = j*dx;
    }

    for (i = 0; i<nx ;i++)
    {
        for(j=0;j<ny;j++)
        {
            b[i][j] = 0.0;
            un[i][j] = 0.0;
            vn[i][j] = 0.0;
        }
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
        b[i][j] = build_up_b(b,rho,dt,u,v,dx,dy);
        p[i][j] = pressure_poisson(p,dx,dy,b);

        for ( i = 1; i < nx; i++)
        {
            for (j = 1; j < ny; j++)
            {
                //b[i][j] = build_up_b(b,rho,dt,u,v,dx,dy);
                //p[i][j] = pressure_poisson(p,pn,dx,dy,b);
                u[i][j] = (un[i][j]-
                         un[i][j] * dt / dx *
                        (un[i][j] - un[i-1][j]) -
                         vn[i][j] * dt / dy *
                        (un[i][j] - un[i][j-1]) -
                         dt / (2 * rho * dx) * (p[i+1][j] - p[i-1][j]) +
                         nu * (dt / dx*dx *
                        (un[i+1][j] - 2 * un[i][j]+ un[i-1][j]) +
                         dt / dy*dy *
                        (un[i][j+1] - 2 * un[i][j] + un[i][j-1]))) ;


                v[i][j] = (vn[i][j]-
                         un[i][j] * dt / dx *
                        (vn[i][j] - vn[i-1][j]) -
                         vn[i][j] * dt / dy *
                        (vn[i][j] - vn[i][j-1]) -
                         dt / (2 * rho * dx) * (p[i+1][j] - p[i-1][j]) +
                         nu * (dt / dx*dx *
                        (vn[i+1][j] - 2 * vn[i][j]+ vn[i-1][j]) +
                         dt / dy*dy *
                        (vn[i][j+1] - 2 * vn[i][j] + vn[i][j-1]))) ;

                u[0][j] = 0.0;
                u[nx][j] = 0.0;
                u[i][0] = 0.0;
                u[i][ny] = 1.0;
                v[0][j] = 0.0;
                v[nx][j] = 0.0;
                v[i][0] = 0.0;
                v[i][ny] = 0.0;
            }
        }
    }

    fprintf(fp,"VARIABLES = \"x\"\t,\"y\"\t,\"u\"\t,\"v\"\t,\"p\"\n");
    for (i=0 ; i < nx ; i++)
    {
        for (j=0 ; j < ny ; j++)
        {
            fprintf(fp,"%0.2f\t %0.2f\t %f\t %f\t %f\n",x[i],y[j],u[i][j],v[i][j],p[i][j]);
        }
    }
}

int main()
{
    int i,j;
    float c = 1.0;
    float dx = 2.0 / (nx - 1.0);
    float dy = 2.0 / (ny - 1.0);
    float rho = 1.0;
    float nu = 0.5;
    float dt = .001;
    FILE *fp;
    fp = fopen("2D_Cavity_Initial_Condition.dat","w");
    static float x[nx],y[ny],u[nx][ny], v[nx][ny], p[nx][ny], b[nx][ny];
    for (i = 0; i < nx ;i++)
    {
        x[i] = i*dx;
        y[i] = i*dy;
    }

    for (i = 0; i<nx ;i++)
    {
        for(j=0;j<ny;j++)
        {
            u[i][j] = 0.0;
            v[i][j] = 0.0;
            p[i][j] = 0.0;
            b[i][j] = 0.0;
        }
    }
/*
    for (i = int(0.5/dx) ; i <= int(1.0/dx+1.0) ; i++ )
    {
        for(j = int(0.5/dy) ; j <= int(1.0/dy+1.0) ; j++)
        {
            u[i][j] = 2.0;
        }
    }*/

    fprintf(fp,"VARIABLES = \"x\"\t,\"y\"\t,\"u\"\t,\"v\"\t,\"p\"\n");
    for (i = 0; i < nx; i++)
    {
      for (j = 0; j < ny; j++)
      {
          fprintf(fp,"%0.2f\t %0.2f\t %f\t %f\t %f\n",x[i],y[j],u[i][j],v[i][j],p[i][j]);
      }
   }
   Cavity_flow(u,v,p,rho,nu,dt,dx,dy);
}
