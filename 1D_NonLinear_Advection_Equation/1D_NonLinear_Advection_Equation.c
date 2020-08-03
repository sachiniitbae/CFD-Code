#include<stdio.h>
#include<math.h>

void FTBS(float domain[], float u[], float un[], float grid, int nt, float dt, float dx)
{
    FILE *fp;
    fp = fopen("FTBS.csv","w");
    int i,j;
    for (j = 0; j < nt; j++)
    {
       for(i=0;i<grid;i++)
        {
            un[i]= u[i];
        }
        for ( i = 1; i<grid; i++)
        {
            u[i] = un[i] - un[i] * dt / dx * (un[i] - un[i-1]);
        }
    }
    for (i=0;i<grid;i++)
    {
        fprintf(fp,"%f,%f\n",domain[i],u[i]);
    }
}

void FTCS(float domain[], float u[], float un[], float grid, int nt, float dt, float dx)
{
    FILE *fp;
    fp = fopen("FTCS.csv","w");
    int i,j;
    for (j = 0; j < nt; j++)
    {
       for(i=0;i<grid;i++)
        {
            un[i]= u[i];
        }
        for ( i = 1; i<grid; i++)
        {
            u[i] = un[i] - un[i] * dt / 2 * dx * (un[i+1] - un[i-1]);
        }
    }
    for (i=0;i<grid;i++)
    {
        fprintf(fp,"%f,%f\n",domain[i],u[i]);
    }
}

void FTFS(float domain[], float u[], float un[], float grid, int nt, float dt, float dx)
{
    FILE *fp;
    fp = fopen("FTFS.csv","w");
    int i,j;
    for (j = 0; j < nt; j++)
    {
       for(i=0;i<grid;i++)
        {
            un[i]= u[i];
        }

        for ( i = 1; i<grid; i++)
        {
            u[i] = un[i] - un[i]* dt / dx * (un[i+1] - un[i]);
        }
    }
    for (i=0;i<grid;i++)
    {
        fprintf(fp,"%f,%f\n",domain[i],u[i]);
    }
}

int main()
{
    int grid = 41,i;
    float dx,nt,dt;
    float domain[grid],u[grid], un[grid];
    dx = 2.0/(grid-1);
    dt = 0.025;
    nt = 25;

    for (i=0; i < grid; i++)
    {
        domain[i] = i*dx;
    }

    for (i = 0; i < grid  ; i++)
    {
        if(0.5 <= domain[i] && domain[i] <= 1.0)
        {
            u[i] = 2.0;
        }
        else
        {
            u[i] = 1.0;
        }
        un[i] = 1.0;
    }
    FTBS(domain,u,un,grid,nt,dt,dx);
    FTCS(domain,u,un,grid,nt,dt,dx);
    FTFS(domain,u,un,grid,nt,dt,dx);
}
