#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>

#define PI 3.14159
#define K  8.988*10E9

int main()
{

  clock_t start, finish;//defines the clock variables
  double duration;
  start = clock();//starts the clock to determine how long the code will run for

  int N=0;//dimension of potential matrix
  printf("Dimension of Matrix (even number greater than or equal to 10): \n");
  scanf("%d", &N);
  int n=N/10;//scaling factor

  double boundary=0;//setting the outer boundary potential to zero
  double check=0.00000001;//convergence criterion
  double omega=1.50;//2/(1+(PI/N));//if running Gauss-Seidel instead of SOR then omega=1
  int counter=0;//sets the iteration counter to initialize at zero
  double dx=1;
  double dy=1;
  double dz=1;
  double q=1;

  double a[N][N];//defines the 2 dimensional array used for potential

  double a_d[N][N][N];//defines the 3 dimensional array used for potential

  //defines the theoretical potential and electric field matrices
  double V[N][N];
  double E_x_t[N][N];
  double E_y_t[N][N];


  //defines the 2 dimensional arrays used for electric field in both x and y
  double E_x[N][N];
  double E_y[N][N];

  //defines the 3 dimensional arrays used for the electric field in both x,y and z
  double E_x_d[N][N][N];
  double E_y_d[N][N][N];
  double E_z_d[N][N][N];


  double r, dmax, E_x_temp, E_y_temp, E_z_temp, E_x_t_temp, E_y_t_temp, temp_dmax, temp_a, temp_V, potential_barrier1, potential_barrier2;//defining some variables...

//**************************Defining all necessary values for the 4 geometries***********************************//

//******Hollow Metallic Plate******//
  int width_hmp=4*n;
  double potential_hmp=10*n;

//****Parallel-Plate Capacitor****//
  int height_pp=5*n;
  int x_pp=3*n;
  int y_pp=3*n;
  double potential_ppl=100;
  double potential_ppr=-100;

//*******Plate and Conductor*****//
  int height_mp=6*n;
  int x_mp=2*n;
  int y_mp=2*n;
  int width_mic=2*n;
  int x_mic=5*n;
  int y_mic=4*n;
  double potential_mp=-100;
  double potential_mic=100;

//*********Lightning Rod*********//
  int height_ned=7*n;
  int x_ned=5*n;
  int y_ned=1*n;
  int len_c=8*n;
  int x_c=1*n;
  int y_c=9*n;
  double potential_ned=-100;
  double potential_c=100;

//*********Point Charge*********//
  int x_pc=5*n;
  int y_pc=5*n;
  int z_pc=5*n;
  double potential_pc=100;

//****************************Determine which Geometry and Method to Use****************************************//
  int geometry=0;
  int method=0;
  printf("Which geometry would you like to implement? (Hollow Metallic Prism (1)\n");
  printf("Parallel Plate Capacitor (2), Metal plate and metallic conductor (3),\n");
  printf("Lightning Rod Model(4), or a Point Charge in 3 dimensions (5))\n");
  scanf("%d", &geometry);

  printf("Which method would you like to use? (Jacobi (1), Gauss-Seidel (2),\n");
  printf("or Successive Over-Relaxation (3))\n");
  scanf("%d", &method);

  if(method==2)
  {
    omega=1;//defining omega=1 for Gauss-Seidel Method
  }

  if(geometry==1)
  {
    potential_barrier1=potential_hmp;
    potential_barrier2=potential_hmp;
  }

  if(geometry==2)
  {
    potential_barrier1=potential_ppr;
    potential_barrier2=potential_ppl;
  }
  if(geometry==3)
  {
    potential_barrier1=potential_mp;
    potential_barrier2=potential_mic;
  }
  if(geometry==4)
  {
    potential_barrier1=potential_ned;
    potential_barrier2=potential_c;
  }
  if(geometry==5)
  {
    potential_barrier1=potential_pc;
    potential_barrier2=potential_pc;
  }


//**************************************************************************************************************//

  FILE *Ep;
  FILE *Vp;
  FILE *fp;
  FILE *ep;
  Ep=fopen("theoreticEp.txt","w");//output calculated electric field
  Vp=fopen("theoretic.txt", "w");//output calculated potential for point charge (x and y comps)
  ep=fopen("efield.txt", "w");//output electric field will be stored in this file
  fp=fopen("laplace.txt","w"); //output potential will be stored in this file


  if(geometry!=5)
  {
    for(int i=0;i<N;i++)//assigning boundary values begins (2d)
    {
      for(int j=0;j<N;j++)
      {
        a[i][j]=boundary;
      }
    }//assigning boundary values ends (2d)
  }


  if(geometry==5)
  {
    for(int i=0;i<N;i++)//assigning boundary values begins (3d)
    {
      for(int j=0;j<N;j++)
      {
        for(int k=0;k<N;k++)
        {
          a_d[i][j][k]=boundary;
        }
      }
    }//ends (3d)
  }

//*************************************************GEOMETRIES**********************************************************//

//********************Hollow metallic prism cross-section********************//
  if(geometry==1)
  {
    for(int i=0;i<=width_hmp/2;i++)
    {
      for(int j=0;j<=width_hmp/2;j++)
      {
        a[width_hmp+i][width_hmp+j]=potential_hmp;
      }//cross section of a hollow metallic prism with the center located on the origin
    }
  }

//****************************Parallel-Plate Capacitor***********************//
  if(geometry==2)
  {
    for(int i=0;i<height_pp;i++)
    {
      a[y_pp+i][x_pp]=potential_ppl;//left
      a[y_pp+i][x_pp+x_pp]=potential_ppr;//right
    }//cross section of a parallel-plate capacitor with the the center located on the origin
  }

//*******Cross-Section of a Metal Plate and a metallic inner conductor*****//
  if(geometry==3)
  {
    for(int i=0;i<height_mp;i++)
    {
      a[y_mp+i][x_mp]=potential_mp;
    }//cross section of a metal plate
    for(int i=0;i<width_mic/2;i++)
    {
      for(int j=0;j<width_mic/2;j++)
      {
        a[y_mic+i][x_mic+j]=potential_mic;
      }
    }//cross section of a metallic inner conductor
  }

//**************************A Model for a Lightning Rod*******************//
  if(geometry==4)
  {
    for(int i=0;i<height_ned;i++)
    {
      a[y_ned+i][x_ned]=potential_ned;
    }//model for the needle
    for(int i=0;i<len_c;i++)
    {
      a[y_c][x_c+i]=potential_c;
    }//model for the cloud
  }


  if(geometry==5)
  {
    a_d[x_pc][y_pc][z_pc]=potential_pc;
  }

//*********************************************************************************************************************//
//************************************Jacobi, Gauss-Seidel, and SOR Methods Implemented*******************************//
  while(1<2)
  {
    dmax=0;
    counter=counter+1;//counts the number of iterations for convergence


//*************************Gauss-Seidel/SOR Method***********************************//
    if(method==2 || method==3  && geometry!=5)
    {
      for(int i=1;i<=N-2;i++)
      {
        for(int j=1;j<=N-2;j++)
        {
          if(a[i][j]==potential_barrier1 || a[i][j]==potential_barrier2)//change the potential value to correspond to the geometry you are working with
          {
            continue;
          }
          temp_a=(a[i+1][j]+a[i-1][j]+a[i][j+1]+a[i][j-1])*0.25;
          r=omega*(temp_a-a[i][j]);//change omega to 1 for Gauss-Seidel
          temp_a=a[i][j]+r;
          temp_dmax=temp_a-a[i][j];
          a[i][j]=temp_a;
          E_x_temp=-(a[i+1][j]-a[i-1][j])/(2*dx);
          E_y_temp=-(a[i][j+1]-a[i][j-1])/(2*dy);
          E_x[i][j]=E_x_temp;
          E_y[i][j]=E_y_temp;

//********Convergence Criterion Check***********//
          if(temp_dmax<0)
          {
            temp_dmax=-1*temp_dmax;
          }

          if(temp_dmax>dmax)
          {
            dmax=temp_dmax;
          }

        }
      }
    }

//******************************Jacobi Method***************************************//
    if(method==1 && geometry!=5)
    {
      double a_old[N][N];
      for(int i=0;i<N;i++)
      {
        for(int j=0;j<N;j++)
        {
          a_old[i][j]=a[i][j];
        }
      }
      for(int i=1;i<=N-2;i++)
      {
        for(int j=1;j<=N-2;j++)
        {
          if(a[i][j]==potential_barrier1 || a[i][j]==potential_barrier2)//change the potential value to correspond to the geometry you are working with
          {
            continue;
          }
          temp_a=(a_old[i+1][j]+a_old[i-1][j]+a_old[i][j+1]+a_old[i][j-1])*0.25;
          temp_dmax=temp_a-a_old[i][j];
          a[i][j]=temp_a;
          E_x_temp=-(a[i+1][j]-a[i-1][j])/(2*dx);
          E_y_temp=-(a[i][j+1]-a[i][j-1])/(2*dy);
          E_x[i][j]=E_x_temp;
          E_y[i][j]=E_y_temp;

//********Convergence Criterion Check***********//
          if(temp_dmax<0)
          {
            temp_dmax=-1*temp_dmax;
          }

          if(temp_dmax>dmax)
          {
            dmax=temp_dmax;
          }

        }
      }
    }

//***********************************Point Charge SOR Method***************************//
    if(method==2 || method==3  && geometry==5)
    {
      for(int i=1;i<=N-2;i++)
      {
        for(int j=1;j<=N-2;j++)
        {
          for(int k=1;k<=N-2;k++)
          {
            if(a_d[i][j][k]==potential_barrier1 || a_d[i][j][k]==potential_barrier2)//change the potential value to correspond to the geo$
            {
              continue;
            }
            temp_a=(a_d[i+1][j][k]+a_d[i-1][j][k]+a_d[i][j+1][k]+a_d[i][j-1][k]+a_d[i][j][k+1]+a_d[i][j][k-1])*0.1666666;
            r=omega*(temp_a-a_d[i][j][k]);//change omega to 1 for Gauss-Seidel
            temp_a=a_d[i][j][k]+r;
            temp_dmax=temp_a-a_d[i][j][k];
            a_d[i][j][k]=temp_a;
            E_x_temp=-(a_d[i+1][j][k]-a_d[i-1][j][k])/(2*dx);
            E_y_temp=-(a_d[i][j+1][k]-a_d[i][j-1][k])/(2*dy);
            E_z_temp=-(a_d[i][j][k+1]-a_d[i][j][k-1])/(2*dz);
            E_x_d[i][j][k]=E_x_temp;
            E_y_d[i][j][k]=E_y_temp;
            E_z_d[i][j][k]=E_z_temp;

//********Convergence Criterion Check***********//
            if(temp_dmax<0)
            {
              temp_dmax=-1*temp_dmax;
            }

            if(temp_dmax>dmax)
            {
              dmax=temp_dmax;
            }

          }
        }
      }
    }

        if(dmax<=check)//exits code if dmax is less than or equal to the desired convergence value 'check'
        {
          for(int i=0;i<N;i++)
          {
            for(int j=0;j<N;j++)
            {
              if(geometry!=5)
              {
                fprintf(ep," %d   %d   %.10lf   %.10lf\n", j, i, E_y[i][j], E_x[i][j]);
              }
              if(geometry==5)
              {
                for(int k=0;k<N;k++)
                {
                  if(k==(n*7))
                  {
                    fprintf(ep,"%d  %d  %.10lf  %.10lf\n",  j, i, E_y_d[i][j][k], E_x_d[i][j][k]);
                  }
                }
              }
            }
          }

          for(int i=0;i<N;i++)
          {
            for(int j=0;j<N;j++)
            {
              if(geometry!=5)
              {
                fprintf(fp,"%.10lf \n",a[i][j]);
              }
              if(geometry==5)
              {
                for(int k=0;k<N;k++)
                {
                if(k==(n*7))
                  {
                    fprintf(fp,"%.10lf \n", a_d[i][j][k]);
                  }
                }
              }
            }
          }
          break;
        }
  }

  for(int i=0;i<N;i++)
  {
    for(int j=0;j<N;j++)
    {
      temp_V=K*q*pow((N/2-i)*(N/2-i)+(N/2-j)*(N/2-j),0.5);
      V[i][j]=temp_V;
    }
  }
  for(int i=0;i<N;i++)
  {
    for(int j=0;j<N;j++)
    {
      fprintf(Vp, "%.10lf\n", V[i][j]);
    }
  }

  for(int i=0;i<N;i++)
  {
    for(int j=0;j<N;j++)
    {
      E_x_t_temp=q/((N/2-i)*(N/2-i));
      E_y_t_temp=q/((N/2-j)*(N/2-j));
      E_x_t[i][j]=E_x_t_temp;
      E_y_t[i][j]=E_y_t_temp;
      fprintf(Ep, "%d  %d %.10lf  %.10lf\n", j, i, E_y_t[i][j], E_x_t[i][j]);
    }
  }

  fclose(Ep);
  fclose(Vp);
  fclose(ep);
  fclose(fp);
  printf("\n");
  printf("Number of iterations required for convergence: %d \n", counter);
  printf("\nData stored in laplace.txt\n");
  printf("\n");

  finish = clock();//stops the clock
  duration = (double)(finish-start)/ CLOCKS_PER_SEC;//gives the final value of the clock to determine how long our code ran for
  printf("Total time taken by CPU: %18.20f seconds \n", duration);//prints the amount of cycles the code ran for


  return 0;
}

