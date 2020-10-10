#include <iostream>
#include <math.h>
using namespace std;
#define nx 243
#define ny 81
#define nt 400
#define nit 50

void Cavity_flow(float u[nx][ny],float v[nx][ny],float p[nx][ny],float u_doublet[nx][ny],float v_doublet[nx][ny],float rho, float nu, float dt, float dx, float dy)
{
    FILE *fp; FILE *fpv; FILE *fpu; FILE *fpp;
    fp = fopen("2D_Cavity_Flow_Doublet_Uniform_final_Solution.dat","w");
    fpv = fopen("V_solution.csv","w");
    fpu = fopen("U_solution.csv","w");
    fpp = fopen("P_solution.csv","w");
    int i,j,k,l;
    float y[ny],x[nx];
    float un[nx][ny],vn[nx][ny];
    float pn[nx][ny];

//Domain Creation
    for (j = 0; j < ny; j++)
    {
        y[j] = j*dy;
    }

    for(i=0; i<nx ;i++)
    {
        x[i] = i*dx;
    }


//Creating a dummy variables to store old values of U ,V and P
    for (i = 0; i<nx ;i++)
    {
        for(j=0;j<ny;j++)
        {
            un[i][j] = 0.0;
            vn[i][j] = 0.0;
            pn[i][j] = 0.0;
        }
    }


    for (k = 0;k < nt; k++)
    {
       for(i=0;i<nx;i++)
        {
            for(j=0;j<ny;j++)
            {
                un[i][j]= u[i][j] ;
                vn[i][j]= v[i][j] ;
            }
        }

        fprintf(fp,"VARIABLES = \"x\"\t\"y\"\t\"u\"\t\"v\"\t\"p\"\n");
        fprintf(fp,"ZONE T=\"%d\", I=%d, J=%d, F = point\n",k,nx,ny);

// Calculation of pressure
        for(l=0; l<nit ; l++)
            {
                for(i=0;i<nx;i++)
                {
                    for(j=0;j<ny;j++)
                    {
                        pn[i][j]= p[i][j];
                    }
                }
                for(i=1;i<nx-1;i++)
                    {
                        for(j=1;j<ny-1;j++)
                        {
                            p[i][j] = ((((pn[i+1][j] + pn[i-1][j]) * (dy*dy) +
                                    (pn[i][j+1] + pn[i][j-1]) * (dx*dx)) /(2 * (dx*dx + dy*dy))) -
                                    ((dx*dx * dy*dy )/
                                    (2 * (dx*dx + dy*dy)))*((rho * ((1 / dt) *((u[i+1][j] - u[i-1][j]) /(2 * dx) + (v[i][j+1] - v[i][j-1]) / (2 * dy))-
                              pow(((u[i+1][j] - u[i-1][j]) / (2 * dx)),2) -
                              2 * ((u[i][j+1] - u[i][j-1]) / (2 * dy) *(v[i+1][j] - v[i-1][j]) / (2 * dx))-
                              (pow(((v[i][j+1] - v[i][j-1]) / (2 * dy)),2))))));

                        }
                    }
// Boundary Conditions
                    for(i=0; i<nx ; i++)
                        {
                            p[i][0] = 0.0;
                            //p[i][0] = p[i][1]; //dp/dy = 0 at y = 0
                            p[i][ny-1] = 0.0;  //p = 0 at y = 2
                            //p[i][ny-1] = p[i][ny-2];
                        }

                    for(j=0; j<ny; j++)
                        {
                            //p[nx-1][j] = 0.0;
                            p[nx-1][j] = p[nx-2][j]; //dp/dx = 0 at x = 2
                            //p[0][j] = p[1][j];   //dp/dx = 0 at x = 0
                            p[0][j] = 0.0;
                        }

            }
// Calculation of x-momentum

        for ( i = 1; i < nx-1; i++)
        {
            for (j = 1; j < ny-1; j++)
            {
                u[i][j] = (un[i][j]- un[i][j] * (dt / dx) *
                        (un[i][j] - un[i-1][j]) -
                         vn[i][j] * (dt / dy) *
                        (un[i][j] - un[i][j-1]) -
                         dt / (2 * rho * dx) * (p[i+1][j] - p[i-1][j]) +
                         nu * ((dt / (dx*dx)) *
                        (un[i+1][j] - 2 * un[i][j]+ un[i-1][j]) +
                         (dt / (dy*dy)) *
                        (un[i][j+1] - 2 * un[i][j] + un[i][j-1]))) ;
            }
        }

// Calculation of y-momentum

        for(i=1;i<nx-1;i++)
        {
            for(j=1;j<ny-1;j++)
            {
                v[i][j] = (vn[i][j]- un[i][j] * (dt / dx )*
                        (vn[i][j] - vn[i-1][j]) -
                         vn[i][j] * (dt / dy) *
                        (vn[i][j] - vn[i][j-1]) -
                         dt / (2 * rho * dy) * (p[i][j+1] - p[i][j-1]) +
                         nu * ((dt / (dx*dx)) *
                        (vn[i+1][j] - 2 * vn[i][j]+ vn[i-1][j]) +
                         (dt / (dy*dy)) *
                        (vn[i][j+1] - 2 * vn[i][j] + vn[i][j-1]))) ;
            }
        }

// Boundary Conditions
        for(j=1;j<ny-1;j++)
            {
                u[0][j] = 1.0;
                v[0][j] = 0.0;
                u[nx-1][j] = 0.0;
                v[nx-1][j] = 0.0;
            }

        for(i=0;i<nx;i++)
            {
                u[i][0] = 1.0;
                v[i][0] = 0.0;
                u[i][ny-1] = 1.0;  // Setting Velocity at lid to  1.0
                v[i][ny-1] = 0.0;
            }

         for (i=0 ; i < nx ; i++)
            {
                for (j=0 ; j < ny ; j++)
                {
                    fprintf(fp,"%f\t %f\t %f\t %f\t %f\n",x[i],y[j],u[i][j] + u_doublet[i][j],v[i][j] + v_doublet[i][j],p[i][j]);
                }
            }
            fprintf(fp,"\n\n");


    }
/*
// Tecplot format
    fprintf(fp,"VARIABLES = \"x\"\t,\"y\"\t,\"u\"\t,\"v\"\t,\"p\"\n");
    for (i=0 ; i < nx ; i++)
    {
        for (j=0 ; j < ny ; j++)
        {
            fprintf(fp,"%f\t %f\t %f\t %f\t %f\n",x[i],y[j],u[i][j],v[i][j],p[i][j]);
        }
    }
*/
    for(j=0;j<ny;j++)
   {
       for(i=0;i<nx;i++)
       {
           fprintf(fpu,"%f,",u[i][j]); fprintf(fpv,"%f,",v[i][j]);
           fprintf(fpp,"%f,",p[i][j]);
       }
       fprintf(fpu,"\n"); fprintf(fpv,"\n"); fprintf(fpp,"\n");
   }
}

int main()
{
    int i,j;
    float dx = 6.0 / (nx - 1.0);
    float dy = 2.0 / (ny - 1.0);
    float rho = 1.0;
    float nu = 0.1;
    float dt = 0.001;
    float kappa = 0.035;
    float x_doublet , y_doublet;
    int x_position,y_position;
    float u[nx][ny], v[nx][ny], p[nx][ny];
    float u_doublet[nx][ny] , v_doublet[nx][ny];
    float x[nx],y[ny];


    for (j = 0; j < ny;j++)
    {
        y[j] = j*dy;
    }

    for(i=0;i<nx;i++)
    {
        x[i] = i*dx;
    }

    for(i=0;i<nx;i++)
    {
        for(j=0;j<ny;j++)
        {
            u_doublet[i][j] = 0.0;
            v_doublet[i][j] = 0.0;
        }
    }

    x_position = 40  ;
    y_position = 40 ;
    x_doublet = x[x_position];
    y_doublet = y[y_position];

    for (i = 0; i<nx ; i++)
    {
        for(j=0; j<ny ; j++)
        {
            if(i == x_position && j == y_position)
            {
                u_doublet[i][j] = 0.0;
                v_doublet[i][j] = 0.0;
            }
            else
            {
                u_doublet[i][j] = (- kappa / (2 * M_PI) * (pow((x[i] - x_doublet),2) - pow((y[j] - y_doublet),2))) / pow((pow((x[i] - x_doublet),2) + pow((y[j] - y_doublet),2)),2);
                v_doublet[i][j] = (- kappa/ (2 * M_PI) * 2 * (x[i] - x_doublet) * (y[j] - y_doublet)) / pow((pow((x[i] - x_doublet),2) + pow((y[j] - y_doublet),2)),2);
            }

        }
    }

//Initializing values of U, V and P
    for (i = 0; i<nx ; i++)
    {
        for(j=0; j<ny ; j++)
        {
            u[i][j] = 0.0;
            v[i][j] = 0.0;
            p[i][j] = 0.0;
        }
    }
   Cavity_flow(u,v,p,u_doublet,v_doublet,rho,nu,dt,dx,dy);
}
