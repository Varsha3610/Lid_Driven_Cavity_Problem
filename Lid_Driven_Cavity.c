#include<stdio.h>
#include<math.h>
#include<stdlib.h>
int main()
{
    int i,j,k=0,Re,m,n;                      //m and n are no. of grids
    printf("Enter the value of Re\tm\tn\t");
    scanf("%d%d%d",&Re,&m,&n);
    double dx=1.0/(m-1);
    double dy=1.0/(n-1);
    double bet;
    bet=dx/dy;
    double psi[m][n],om[m][n],psi_prev[m][n],om_prev[m][n],u[m][n],v[m][n],x,y,err_psi,err_om,err1_psi=0.0,err1_om=0.0;

//Initializing the interior points
for(j= 1; j < n-1; j++)
{
   for (i = 1; i < m-1; i++)
   {
       psi[i][j] = 0;
       u[i][j] = 0;
       v[i][j] = 0;
       om[i][j] = 0;
   }
}
//Boundary conditions
for(j= 0; j < n; j++)
{
   psi[0][j] =0;
   psi[m-1][j] =0;
   u[0][j] = 0;
   u[m-1][j] = 0;
   v[0][j] =0;
   v[m-1][j] =0;
   om[0][j] =0;
   om[m-1][j] =0;
}
for(i = 0; i < m; i++)
{
    psi[i][0] = 0;
    psi[i][n-1] = 0;
    u[i][0] = 0;
    u[i][n-1] = 1;
    v[i][0] = 0;
    v[i][n-1] = 0;
    om[i][0] = 0;
    om[i][n-1] = -2/dy;
}
//Gauss-seidel Method
    FILE *g;
    g=fopen("result.dat","w");
      do
      {
         for(j=0; j<n; j++)
             {
                 for(i=0; i<m; i++)
                 {
                   psi_prev[i][j]=psi[i][j];
                   om_prev[i][j]=om[i][j];
                 }
             }

//Solving for stream function
       for (j=1; j<(n-1); j++)
       {
         for(i=1; i<(m-1); i++)
         {
            psi[i][j] = (psi[i+1][j]+psi[i-1][j]+bet*bet*(psi[i][j+1]+psi[i][j-1])+dx*dx*om[i][j])/(2*(1+pow(bet,2.0)));
         }
       }

//Solving for vorticity function

    for (j=1; j<(n-1); j++)
    {
       for(i=1; i<(m-1); i++)
       {
          om[i][j] = (0.5/(1.0+pow(bet,2)))*((1.0-(psi[i][j+1]-psi[i][j-1])*((bet*Re)/4.0))*om[i+1][j]+ (1.0+(psi[i][j+1]-psi[i][j-1])*((bet*Re)/4.0))*om[i-1][j]
                     + (1.0+(psi[i+1][j]-psi[i-1][j])*(Re/(4.0*bet)))*(pow(bet,2)*om[i][j+1])+ (1.0-(psi[i+1][j]-psi[i-1][j])*(Re/(4.0*bet)))*(pow(bet,2)*om[i][j-1]));
       }
    }
//Update vorticity at boundaries

    for(j=0; j<n; j++)
       {
           om[0][j]=(-2*(psi[1][j]-psi[0][j]))/(dx*dx);
           om[m-1][j]=(-2*(psi[m-2][j]-psi[m-1][j]))/(dx*dx);
       }
    for(i=0; i<m; i++)
       {
            om[i][0]=(-2*(psi[i][1]-psi[i][0]))/(dy*dy);
            om[i][n-1]=(-2*(psi[i][n-2]-psi[i][n-1]+dy))/(dy*dy);
       }

//Error calculation
err_psi=0.0;
err_om=0.0;
  for (j=1; j<(n-1); j++)
  {
     for(i=1; i<(m-1); i++)
     {
         err_psi+=pow((psi[i][j]-psi_prev[i][j]),2.0);
         err_om+=pow((om[i][j]-om_prev[i][j]),2.0);
     }
  }
      err1_psi=sqrt(err_psi)/((m-2)*(n-2));
      err1_om=sqrt(err_om)/((m-2)*(n-2));

      printf("iteration=%d\t",k);
      printf("err1_psi=%.9lf\terr1_om=%.9lf\n",err1_psi,err1_om);
      k=k+1;
    }while(err1_psi>pow(10,-6) || err1_om>pow(10,-6));

//Update velocities u and v
for (j=1; j<(n-1); j++)
{
    for(i=1; i<(m-1); i++)
    {
        u[i][j]=(psi[i][j+1]-psi[i][j-1])/(2*dy);
        v[i][j]=-(psi[i+1][j]-psi[i-1][j])/(2*dx);
    }
}

fprintf(g,"ZONE I=%d, J=%d\n",m,n);
for(j = 0; j < n; j++)
{
   y=j*dy;
for(i = 0; i < m; i++)
{
   x=i*dx;
fprintf(g,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",x,y,u[i][j],v[i][j],psi[i][j],om[i][j]);
}
}
fclose(g);
FILE *uv,*vv;
uv=fopen("U_velocity100.dat","w");
 fprintf(uv,"\nU-velocity\t\tY\n",n);
 for(j = 0; j < n; j++)
 {
 fprintf(uv,"%lf \t %lf \n",u[50][j],j*dy); /* saving of u-velocity values into file along with grids */
 }
 fclose(uv);
 vv=fopen("V_velocity100.dat","w");
 fprintf(vv,"\nX\t\tV-velocity\n",m);
 for(i = 0; i < m; i++)
 {
 fprintf(vv,"%lf \t %lf \n",i*dx,v[i][50]); /* saving of v-velocity values into file along with grids */
 }
 fclose(vv);
return 0;
}




