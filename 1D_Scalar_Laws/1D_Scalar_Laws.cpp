#include <iostream>
#include <cmath>
#include <vector>
#include <cstring>
#define nx 51      // Number of grid points
#define t_end1 2.0  // End time for 1st equation
#define t_end2 0.5  // End time for 2nd equation
#define t_end3 0.5  // End time for 3rd equation
#define c1 0.8      // dt/dx for 1st equation
#define c2 0.5      // dt/dx for 2nd equation
#define c3 0.5      // dt/dx for 3rd equation
using namespace std;


vector <float> Domain()    // Returns a vector of x (domain)
{
    vector <float> x;
    float dx = 0.0;
    float x_start = -1.0;
    float x_end = 1.0;
    dx = (x_end - x_start) / (nx-1);

    for (int i = 0; i < nx ; i++)
    {
        x.push_back(x_start + i*dx);
    }
    return x;
}

vector <float> U_IC()     // Initialization of U a.k.a initial condition and returns a vector
{
    vector <float> a,x;
    x = Domain();
    for(int i=0; i<nx; i++)
    {
        if(-0.33 <= x[i] && x[i] <= 0.33 )
        {
            a.push_back(1.5);
        }
        else
        {
            a.push_back(0.5);
        }
    }
    return a;
}

float sigma(float u1,float u2, float f1, float f2)  // Signal function for 1st Order
{
    float sig = 0.0;
    if(((f1-f2) > 0 && (u1-u2) > 0) || ((f1-f2) < 0 && (u1-u2) < 0))
        sig = 1.0;
    else if(((f1-f2) > 0 && (u1-u2) < 0) || ((f1-f2) < 0 && (u1-u2) > 0))
        sig = -1.0;
    else
        sig = 0.0;
    return sig;
}

float NU(float u1, float u2, float f1, float f2, float dt, float dx)  // Signal function for 2nd Order
{
    float nu = 0.0;
    if((f1-f2) == (u1-u2))
        nu = 1.0 * (dt/dx);
    else
        nu = ((f1-f2)/(u1-u2))*(dt/dx);
    return nu;
}

float fluc(float f1, float f2)  // Returns fluctuation
{
    return f1-f2;
}

vector <float> First_Order(vector <float> &u_1, vector <float> &un_1, vector <float> &fn_1, float &dt, float &dx)
{
    float sigmap = 0.0,sigman = 0.0, flp_1 = 0.0, fln_1 = 0.0;
    for(int i =0;i<nx;i++)
    {
        if(i == 0)   //Periodic Boundary condition at i=0
            {
                sigmap = sigma(un_1[i+1],un_1[i],fn_1[i+1],fn_1[i]);
                sigman = sigma(un_1[i],un_1[nx-1],fn_1[i],fn_1[nx-1]);
                flp_1 = fluc(fn_1[i+1],fn_1[i]);
                fln_1 = fluc(fn_1[i],fn_1[nx-1]);
                u_1[i] = un_1[i] - (dt/dx)*(0.5*(1+sigman)*fln_1 + 0.5*(1-sigmap)*flp_1);
            }

        else if(i == nx-1)   //Periodic Boundary condition at i = last_grid_point
            {
                sigmap = sigma(un_1[0],un_1[i],fn_1[0],fn_1[i]);
                sigman = sigma(un_1[i],un_1[i-1],fn_1[i],fn_1[i-1]);
                flp_1 = fluc(fn_1[0],fn_1[i]);
                fln_1 = fluc(fn_1[i],fn_1[i-1]);
                u_1[i] = un_1[i] - (dt/dx)*(0.5*(1+sigman)*fln_1 + 0.5*(1-sigmap)*flp_1);
            }

        else
            {
                sigmap = sigma(un_1[i+1],un_1[i],fn_1[i+1],fn_1[i]);
                sigman = sigma(un_1[i],un_1[i-1],fn_1[i],fn_1[i-1]);
                flp_1 = fluc(fn_1[i+1],fn_1[i]);
                fln_1 = fluc(fn_1[i],fn_1[i-1]);
                u_1[i] = un_1[i] - (dt/dx)*(0.5*(1+sigman)*fln_1 + 0.5*(1-sigmap)*flp_1);
            }
    }
    return u_1;
}

vector <float> Second_Order(vector <float> &u_2, vector <float> &un_2, vector <float> &fn_2, float &dt, float &dx)
{
    float nup = 0.0, nun = 0.0, flp_2 = 0.0, fln_2 = 0.0;
    for(int i = 0 ; i<nx ; i++)
    {
        if(i == 0)    //Periodic Boundary condition at i=0
            {
                nup = NU(un_2[i+1],un_2[i],fn_2[i+1],fn_2[i],dt,dx);
                nun = NU(un_2[i],un_2[nx-1],fn_2[i],fn_2[nx-1],dt,dx);
                flp_2 = fluc(fn_2[i+1],fn_2[i]);
                fln_2 = fluc(fn_2[i],fn_2[nx-1]);
                u_2[i] = un_2[i] - (dt/dx)*(0.5*(1+nun)*fln_2 + 0.5*(1-nup)*flp_2);
            }

        else if(i == nx-1)   //Periodic Boundary condition at i = last_grid_point
            {
                nup = NU(un_2[0],un_2[i],fn_2[0],fn_2[i],dt,dx);
                nun = NU(un_2[i],un_2[i-1],fn_2[i],fn_2[i-1],dt,dx);
                flp_2 = fluc(fn_2[0],fn_2[i]);
                fln_2 = fluc(fn_2[i],fn_2[i-1]);
                u_2[i] = un_2[i] - (dt/dx)*(0.5*(1+nun)*fln_2 + 0.5*(1-nup)*flp_2);
            }

        else
            {
                nup = NU(un_2[i+1],un_2[i],fn_2[i+1],fn_2[i],dt,dx);
                nun = NU(un_2[i],un_2[i-1],fn_2[i],fn_2[i-1],dt,dx);
                flp_2 = fluc(fn_2[i+1],fn_2[i]);
                fln_2 = fluc(fn_2[i],fn_2[i-1]);
                u_2[i] = un_2[i] - (dt/dx)*(0.5*(1+nun)*fln_2 + 0.5*(1-nup)*flp_2);
            }
    }
    return u_2;
}


void linear(vector <float> &x, vector <float> &u_1, vector <float> &u_2, vector <float> &un_1, vector <float> &un_2,float &dt, float &dx)
{
    FILE *fpf;
    fpf = fopen("1st_equation_final.csv","w");
    float sigmap = 0.0,sigman = 0.0, flp_1 = 0.0, fln_1 = 0.0, flp_2 = 0.0, fln_2 = 0.0, nup = 0.0, nun= 0.0;
    vector <float> fn_1,fn_2;

    fn_1 = u_1;  //Initialization of flux term
    fn_2 = u_2;  //Initialization of flux term
    for (float t = dt; t <= t_end1 ; t+=dt)   //Time loop
    {
        un_1 = u_1;  //Saving the old values of u in another same kind of vector for 1st Order
        fn_1 = u_1;  //Saving the old values of flux in another same kind of vector for 1st Order
        un_2 = u_2;  //Saving the old values of u in another same kind of vector for 2nd Order
        fn_2 = u_2;  //Saving the old values of flux in another same kind of vector for 2nd Order

        u_1 = First_Order(u_1,un_1,fn_1,dt,dx);  // Calling 1st order function defined earlier to calculate
        u_2 = Second_Order(u_2,un_2,fn_2,dt,dx); // Calling 2nd order function  defined earlier to calculate

    }

    fprintf(fpf,"x,u_1st_order,u_2nd_order\n");

    for (int i=0;i<nx;i++)
    {
        fprintf(fpf,"%f,%f,%f\n",x[i],u_1[i],u_2[i]);  //Printing in CSV format
    }
}

void burgers(vector <float> &x, vector <float> &u_1, vector <float> &u_2, vector <float> &un_1, vector <float> &un_2,float &dt, float &dx)
{
    FILE *fpf;
    fpf = fopen("2nd_equation_final.csv","w");
    float sigmap = 0.0,sigman = 0.0, flp_1 = 0.0, fln_1 = 0.0, flp_2 = 0.0, fln_2 = 0.0, nup = 0.0,nun= 0.0;
    vector <float> fn_1,fn_2;

    for(int i=0 ;i<nx ;i++)
    {
        fn_1.push_back((u_1.at(i) * u_1.at(i))/2);   //Initialization of flux term
        fn_2.push_back((u_2.at(i) * u_2.at(i))/2);   //Initialization of flux term
    }

    for (float t = dt; t <= t_end2 ; t+=dt)  //Time loop
    {

        for(int i=0;i<nx;i++)
        {
            fn_1[i] = (u_1.at(i) * u_1.at(i))/2 ;  //Saving the old values of flux in another same kind of vector for 1st Order
            fn_2[i] = (u_2.at(i) * u_2.at(i))/2 ;  //Saving the old values of flux in another same kind of vector for 2nd Order
        }
        un_1 = u_1;   //Saving the old values of u in another same kind of vector for 1st Order
        un_2 = u_2;   //Saving the old values of u in another same kind of vector for 2nd Order

        u_1 = First_Order(u_1,un_1,fn_1,dt,dx);
        u_2 = Second_Order(u_2,un_2,fn_2,dt,dx);
    }

    fprintf(fpf,"x,u_1st_order,u_2nd_order\n");
    for (int i=0;i<nx;i++)
    {
        fprintf(fpf,"%f,%f,%f\n",x[i],u_1[i],u_2[i]);
    }
}

void concave(vector <float> &x, vector <float> &u_1, vector <float> &u_2, vector <float> &un_1, vector <float> &un_2,float &dt, float &dx)
{
    FILE *fpf;
    fpf = fopen("3rd_equation_final.csv","w");
    float sigmap = 0.0,sigman = 0.0, flp_1 = 0.0, fln_1 = 0.0, flp_2 = 0.0, fln_2 = 0.0, nup = 0.0,nun= 0.0;
    vector <float> fn_1,fn_2;

    for(int i=0 ;i<nx ;i++)
    {
        fn_1.push_back(u_1.at(i) * (1 - u_1.at(i)));
        fn_2.push_back(u_2.at(i) * (1 - u_2.at(i)));
    }

    for (float t = dt; t <= t_end3 ; t+=dt)  //Time loop
    {

        for(int i=0;i<nx;i++)
        {
            fn_1[i] = (u_1.at(i) * (1 - u_1.at(i)));
            fn_2[i] = (u_2.at(i) * (1 - u_2.at(i))) ;
        }
        un_1 = u_1;
        un_2 = u_2;

        u_1 = First_Order(u_1,un_1,fn_1,dt,dx);
        u_2 = Second_Order(u_2,un_2,fn_2,dt,dx);
    }

    fprintf(fpf,"x,u_1st_order,u_2nd_order\n");

    for (int i=0;i<nx;i++)
    {
        fprintf(fpf,"%f,%f,%f\n",x[i],u_1[i],u_2[i]);
    }
}

void burgers_exact(vector <float> &x, vector <float> &u_exa)
{
    FILE *fpf;
    fpf = fopen("2nd_equation_exact.csv","w");
    for(int i=0;i<nx;i++)
    {
        u_exa[i] = 0.5;
        if(x[i]<=0.83 && x[i]>=0.42)
            u_exa[i] = 1.5;
        else if(x[i]>-0.08 && x[i]<0.42)
            u_exa[i] = 0.5 + ((x[i]+0.08)/0.5);
        fprintf(fpf,"%f,%f\n",x[i],u_exa[i]);
    }
}

void concave_exact(vector <float> &x, vector <float> &u_exa)
{
    FILE *fpf;
    fpf = fopen("3rd_equation_exact.csv","w");
    for(int i=0;i<nx;i++)
    {
        u_exa[i] = 0.5;
        if(x[i] >= -0.83 && x[i] <= -0.67)
            u_exa[i] = 1.5;
        else if(x[i] < 0.33 && x[i] > -0.67)
            u_exa[i] = ((0.32 - x[i])) + 0.5;
        fprintf(fpf,"%f,%f\n",x[i],u_exa[i]);
    }
}

int main()
{
    vector <float> x,u_1,u_2,un_1,un_2,u_exa;
    float dt1 = 0.0 ,dt2 = 0.0, dt3 = 0.0 ;
    float dx = 0.0;
    float x_start = -1.0;
    float x_end = 1.0;
    dx = (x_end - x_start) / (nx-1);
    FILE *fpi;
    fpi = fopen("Initial_Condition.csv","w");
    dt1 = c1*dx;
    dt2 = c2*dx;
    dt3 = c3*dx;

    x = Domain();
    u_1 = U_IC();
    u_2 = U_IC();
    un_1 = U_IC();
    un_2 = U_IC();
    linear(x,u_1,u_2,un_1,un_2,dt1,dx);

    x = Domain();
    u_1 = U_IC();
    u_2 = U_IC();
    un_1 = U_IC();
    un_2 = U_IC();
    burgers(x,u_1,u_2,un_1,un_2,dt2,dx);

    x = Domain();
    u_1 = U_IC();
    u_2 = U_IC();
    un_1 = U_IC();
    un_2 = U_IC();
    concave(x,u_1,u_2,un_1,un_2,dt3,dx);

    x = Domain();
    u_exa = U_IC();
    burgers_exact(x,u_exa);

    x = Domain();
    u_exa = U_IC();
    concave_exact(x,u_exa);

    for (int i=0;i<nx;i++)
    {
        fprintf(fpi,"%f,%f\n",x[i],u_1[i]);
    }
}
