#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void TDMA(double*,double*,double*,double*,double*,int);

main()
{
	printf("Starting....\n");
	//Input Nx,Ny,Nt,Lx,Ly,alpha
	double alpha,Lx,Ly,T1,T2,T3,T4;
	Lx=1.0;
	Ly=1.0;
	alpha=1;

	//Wall Temperatures:
	T1=800;
	T2=200;
	T3=400;
	T4=600;

	//Initial temperature:
	T0=T2;

	int i,j,k,Nx,Ny,Nt;
	double lx,ly,***T,*ax,*ay,*bx,*by,*cx,*cy,*dx,*dy,*tempX,*tempY;

	Nx=5;
	Ny=5;
	Nt=100000;

	double dt,delx,dely;

	dt=0.01;
	delx=Lx/Nx;
	dely=Ly/Ny;

	lx=alpha*dt/(delx*delx);
	ly=alpha*dt/(dely*dely);

	printf("lamdas: %.14f %.14f\n",lx,ly);

	printf("Temperature array in creation\n");
	T=(double***)malloc(Nt*sizeof(double**));
	for (k = 0; k < Nt; ++k)
	{
		T[k]=(double**)malloc(Nx*sizeof(double*));
		for (i = 0; i < Nx; ++i)
		{
			T[k][i]=(double*)malloc(Ny*sizeof(double));
		}
	}
	printf("Temperature array created\n");
	//initial condition:
	for (i = 0; i < Nx; ++i)
	{
		for (j = 0; j < Ny; ++j)
		{
			T[0][i][j]=T0;
		}
	}
	//boundary condition:
	for (k = 0; k < Nt; ++k)
	{
		for (i = 0; i < Nx; ++i)
		{
			T[k][i][0]=T2;
		}
		for (i = 0; i < Nx; ++i)
		{
			T[k][i][Ny-1]=T4;
		}
		for (j = 0; j < Ny; ++j)
		{
			T[k][0][j]=T1;
		}
		for (j = 0; j < Ny; ++j)
		{
			T[k][Nx-1][j]=T3;
		}
	}
	printf("Boundary conditions assigned\n");

	ax=(double*)malloc(Nx*sizeof(double));
	bx=(double*)malloc(Nx*sizeof(double));
	cx=(double*)malloc(Nx*sizeof(double));
	dx=(double*)malloc(Nx*sizeof(double));
	tempX=(double*)malloc(Nx*sizeof(double));

	ay=(double*)malloc(Ny*sizeof(double));
	by=(double*)malloc(Ny*sizeof(double));
	cy=(double*)malloc(Ny*sizeof(double));
	dy=(double*)malloc(Ny*sizeof(double));
	tempY=(double*)malloc(Ny*sizeof(double));

	ax[0]=1;
	bx[0]=0;
	cx[0]=0;
	for (i = 1; i <  Nx-1; ++i)
	{
		ax[i]=1+lx;
		bx[i]=-0.5*lx;
		cx[i]=-0.5*lx;
	}
	ax[Nx-1]=1;
	bx[Nx-1]=0;
	cx[Nx-1]=0;

	ay[0]=1;
	by[0]=0;
	cy[0]=0;
	for (i = 1; i <  Ny-1; ++i)
	{
		ay[i]=1+ly;
		by[i]=-0.5*ly;
		cy[i]=-0.5*ly;
	}
	ay[Ny-1]=1;
	by[Ny-1]=0;
	cy[Ny-1]=0;	

	printf("Coefficients of TDMA calculated\n");

	for (k = 0; k < Nt-1; ++k)
	{
		//printf("%d\n",k );
		if(k%2==0)
		{
			// for (i = 0; i < Nx; ++i)
			// {
			// 	T[k+1][i][0]=T2;
			// }
			for (j = 1; j < Ny-1; ++j)
			{
				dx[0]=T1;
				for (i = 1; i < Nx-1; ++i)
				{
					dx[i]=(T[k][i][j]+0.5*lx*(T[k][i+1][j]-2*T[k][i][j]+T[k][i-1][j])+ly*(T[k][i][j+1]-2*T[k][i][j]+T[k][i][j-1]));
					//printf("\t%.14f\n",dx[i]);
				}
				dx[Nx-1]=T3;

				TDMA(tempX,ax,cx,bx,dx,Nx);

				for (i = 0; i < Nx; ++i)
				{
					T[k+1][i][j]=tempX[i];
					//printf("\t\t%.14f\n",tempX[i]);
				}
			}
			// for (i = 0; i < Nx; ++i)
			// {
			// 	T[k+1][i][Ny-1]=T4;
			// }
		}
		else
		{
			// for (j = 0; j < Ny; ++j)
			// {
			// 	T[k+1][0][j]=T1;
			// }
			for (i = 1; i < Nx-1; ++i)
			{
				dy[0]=T2;
				for (j = 1; j < Ny-1; ++j)
				{
					dy[j]=(T[k][i][j]+lx*(T[k][i+1][j]-2*T[k][i][j]+T[k][i-1][j])+0.5*ly*(T[k][i][j+1]-2*T[k][i][j]+T[k][i][j-1]));
					//printf("%f\n",dy[j]);
				}
				dy[Ny-1]=T4;

				TDMA(tempY,ay,cy,by,dy,Ny);

				for (j = 0; j < Ny; ++j)
				{
					T[k+1][i][j]=tempY[j];
				}
			}
			// for (j = 0; j < Ny; ++j)
			// {
			// 	T[k+1][Nx-1][j]=T3;
			// }
		}
	}

	FILE*f_out;

	f_out=fopen("Temperature_T_x_y.txt","w");
	for(k=0;k<Nt;++k)
	{
		for(i=0;i<Nx;i++)
		{
			for(j=0;j<Ny;j++)
			{
				fprintf(f_out,"%.14f\n",T[k][i][j]);
			}
		}
	}
	fclose(f_out);

	for(i=0;i<Nx;i++)
	{
		for(j=0;j<Ny;j++)
		{
			printf("%.14f\t",T[Nt-1][i][j]);
		}
		printf("\n");
	}
}

void TDMA(double*x,double*a,double*b,double*c,double*d,int n)
{
    int i;
    double*b_,*d_;
    b_=(double*)malloc(n*sizeof(double));
    d_=(double*)malloc(n*sizeof(double));
    b_[0]=b[0]/a[0];
    d_[0]=d[0]/a[0];
    //printf("%f\n",a[0]);
    //printf("%f\t%f\n",b_[0],d_[0]);
    for(i=1;i<n-1;i++)
    {
        b_[i]=b[i]/(a[i]-c[i]*b_[i-1]);
        d_[i]=(d[i]-c[i]*d_[i-1])/(a[i]-c[i]*b_[i-1]);
        //printf("%f\t%f\n",b_[i],d_[i]);
    }
    d_[n-1]=(d[n-1]-c[n-1]*d_[n-2])/(a[n-1]-c[n-1]*b_[n-2]);
    //printf("\t\t%f\n",d_[n-1]);
    x[n-1]=d_[n-1];
    for(i=n-2;i>=0;i--)
    {
        x[i]=d_[i]-b_[i]*x[i+1];
        //printf("\t%f\n",x[i]);
    }
    return;
}