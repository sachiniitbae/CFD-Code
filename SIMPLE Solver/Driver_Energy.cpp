#include <iostream>
#include <cmath>
#include <vector>
#include <cstdlib>
#include <algorithm>
#include "print_vector.h"
#include "linspace.h"
#include "u_momentum.h"
#include "v_momentum.h"
#include "energy_equation.h"
#include "get_rhs.h"
#include "coefficient_matrix.h"
#include "pressure_correct.h"
#include "update_velocity.h"
#include "check_divergence.h"
using namespace std;
#define x_start        -0.0
#define x_end           1.0
#define y_start        -0.0
#define y_end           1.0

int main()
{
    FILE *fp;
    fp = fopen("Natural_Convection_Ra1e4_33_t0.5_Cavity.plt","w+t");
    FILE *fpr;
    fpr = fopen("Residual_Ra1e5_33_t0.5.csv","w");
    int imax = 11;
    int jmax = 11;
    int max_iteration = 2;
    float maxRes = 1000;
    int iteration = 6000;
    float Pr = 0.71;
    float Ra = 1e4;
    float mu = sqrt(Pr/Ra);                         //viscosity
    float D = 1 / sqrt(Ra*Pr);
    float rho = 1.0;                                //density
    float velocity = 0.0;                           //velocity = lid velocity
    float alphaP = 0.3;                             //pressure under-relaxation
    float alphaU = 0.9;                             //velocity under-relaxation
    float alphaV = 0.9;
    float alphaT = 1.0;
    float tol = 1e-5;
    float T_top = -0.5;
    float T_bottom = 0.5;
    unsigned int i,j;
    float dx = (x_end - x_start) / (imax-1);
    float dy = (y_end - y_start) / (jmax-1);

    vector<float> x = linspace_interval(x_start , x_end , dx,imax);
    //cout << x.size() << endl;
    vector<float> y = linspace_interval(y_start, y_end, dy,jmax);
    //cout << y.size() << endl;
    //print_vector_1d(x);
    //print_vector_1d(y);

    //Variable declaration
    vector<vector<float>> p(imax, vector<float> (jmax,0));                //   p = Pressure
    vector<vector<float>> p_star(imax, vector<float> (jmax,0));           //  Guessed Pressure
    vector<vector<float>> p_prime(imax, vector<float> (jmax,0));          //   pressure correction
    vector<vector<float>> divergence(imax, vector<float> (jmax,0));
    vector<float> rhsp(imax + (jmax-1)*jmax, 0);             //   Right hand side vector of pressure correction equation
    vector<vector<float>> Ap(imax*jmax, vector<float> (imax*jmax,0));
    vector<vector <float>> diagonals(5, vector<float> (imax*jmax));
    vector<vector <float>> div(imax, vector<float> (jmax,0));

    vector < vector <float>> T(imax,vector<float>(jmax,0));
    vector < vector <float>> Told(imax,vector<float>(jmax,0));
    vector < vector <float>> TRes(imax,vector<float>(jmax,0));

    //Vertical velocity
    vector<vector<float>> v_star(imax, vector<float> (jmax+1,0));
    vector<vector<float>> vold (imax, vector<float> (jmax+1,0));
    vector<vector<float>> vRes(imax, vector<float> (jmax+1,0));
    vector<vector<float>> v (imax, vector<float> (jmax+1,0));
    vector<vector<float>> d_v(imax, vector<float> (jmax+1,0));    //velocity correction coefficient

    // Horizontal Velocity
    vector<vector<float>> u_star(imax+1, vector<float> (jmax,0));
    vector<vector<float>> uold(imax+1, vector<float> (jmax,0));
    vector<vector<float>> uRes(imax+1, vector<float> (jmax,0));
    vector<vector<float>> u(imax+1, vector<float> (jmax,0));
    vector<vector<float>> d_u(imax+1, vector<float> (jmax,0));  //velocity correction coefficient


    //Boundary condition
    //Lid velocity (Top wall is moving with 1m/s)
    for(i=0;i<u_star.size();i++)
    {
        u_star[i][jmax-1] = velocity;
        u[i][jmax-1] = velocity;
    }

    for(i=0;i<T.size();i++)
    {
        T[i][jmax-1] = T_top;     //top temperature
        Told[i][jmax-1] = T_top;
        T[i][0] = T_bottom;       //bottom temperature
        Told[i][0] = T_bottom;
    }

    //print_vector_2d(T);


    // ---------- iterations -------------------//
    while ( (iteration <= max_iteration) && (maxRes >= tol) )
    {
        //cout << "Confirmation" << endl;
        iteration = iteration + 1;
        //print_vector_2d(u);
        //print_vector_2d(v);
        //print_vector_2d(u_star);
        u_momentum_new(imax,jmax,dx,dy,rho,mu,u,u_star,d_u,v,p_star,velocity,alphaU);
        //print_vector_2d(u);
        //print_vector_2d(u_star);
        v_momentum_new(imax,jmax,dx,dy,rho,mu,u,v,v_star,d_v,p_star,T,alphaV);
        //print_vector_2d(v);
        //print_vector_2d(v_star);
        uold = u;
        //print_vector_2d(uold);
        vold = v;
        //print_vector_2d(vold);
        get_rhs(imax,jmax,dx,dy,rho,u_star,v_star,rhsp);
        //print_vector_1d(rhsp);
        coefficient_matrix(imax,jmax,dx,dy,rho,d_u,d_v,Ap);
        //print_vector_2d(Ap);
        p = pressure_correct(imax,jmax,rhsp,Ap,p_star,p_prime,alphaP);
        //print_vector_2d(p);
        //print_vector_2d(p_prime);
        update_velocity(imax,jmax,u,u_star,v,v_star,p_prime,d_u,d_v,velocity);
        //print_vector_2d(u);
        //print_vector_2d(v);
        check_divergence(imax,jmax,dx,dy,u,v,div);
        //print_vector_2d(div);
        energy_equation(imax,jmax,dx,dy,u,v,T,Told,D,alphaT,T_top,T_bottom);
        //print_vector_2d(T);


        //find maximum residual in the domain
        for(i=0;i<v.size();i++)
        {
            for(j=0;j<v[0].size();j++)
            {
                vRes[i][j] = abs(v[i][j] - vold[i][j]);
            }
        }

        for(i=0;i<u.size();i++)
        {
            for(j=0;j<u[0].size();j++)
            {
                uRes[i][j] = abs(u[i][j] - uold[i][j]);
            }
        }
        for(i=0;i<T.size();i++)
        {
            for(j=0;j<T[0].size();j++)
            {
                TRes[i][j] = abs(T[i][j] - Told[i][j]);
            }
        }

        vector<float> oneDimU,oneDimV,oneDimT;
        for(int i = 0; i < imax; i++)
        {
            for(int j = 0; j < jmax; j++)
            {
                oneDimU.push_back(uRes[i][j]);
                oneDimV.push_back(vRes[i][j]);
                oneDimT.push_back(TRes[i][j]);
            }
        }
        float maxRes_u = *max_element(oneDimU.begin(),oneDimU.end());
        //cout << maxRes_u << endl;
        float maxRes_v = *max_element(oneDimV.begin(),oneDimV.end());
        //cout << maxRes_v << endl;
        float maxRes_T = *max_element(oneDimT.begin(),oneDimT.end());
        //cout << maxRes_T << endl;
        vector<float> vec = {maxRes_u, maxRes_v, maxRes_T};
        maxRes = *max_element(vec.begin(),vec.end());
        //maxRes = fmax(maxRes1, maxRes_T);

        cout << "It: "<<iteration<< " " <<"Res: "<< maxRes <<endl;

        if (maxRes > 2)
        {
            cout << "Solution is not converging!" << endl;
            break;
        }

        p_star = p;
        Told = T;

/*
        fprintf(fp,"VARIABLES = \"x\"\t\"y\"\t\"u\"\t\"v\"\t\"vres\"\t\"T\"\t\"p\"\n");
        fprintf(fp,"ZONE T=\"%d\", I=%d, J=%d, F = point\n",iteration,imax,jmax);
        for (i=0 ; i < imax ; i++)
        {
            for (j=0 ; j < jmax ; j++)
            {
                fprintf(fp,"%f\t %f\t %f\t %f\t %f\t %f\t %f\n",x[i],y[j],u[i][j],v[i][j],sqrt(pow(u[i][j],2)+ pow(v[i][j],2)), T[i][j], p[i][j]);
            }
        }
        fprintf(fp,"\n\n");
*/
        //fprintf(fpr,"Iteration,Residual");
        fprintf(fpr,"%d,%f\n",iteration,maxRes);

    }


    fprintf(fp,"VARIABLES = \"x\"\t\"y\"\t\"u\"\t\"v\"\t\"vres\"\t\"T\"\t\"p\"\n");
    fprintf(fp,"ZONE T=\"%d\", I=%d, J=%d, F = point\n",iteration,imax,jmax);
    for (i=0 ; i < imax ; i++)
        {
            for (j=0 ; j < jmax ; j++)
            {
                fprintf(fp,"%f\t %f\t %f\t %f\t %f\t %f\t %f\n",x[i],y[j],u[i][j],v[i][j],sqrt(pow(u[i][j],2)+ pow(v[i][j],2)), T[i][j], p[i][j]);
            }
        }


    return 0;

}
