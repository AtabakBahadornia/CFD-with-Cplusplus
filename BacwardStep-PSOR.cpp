////////////Project 1 - Cavity////////////

#include <iostream>
#include <iomanip>
#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <fstream>
using namespace std;
int main()
{
/////////////////////////////////Variables////////////////////////////////////////////
double L,H,Delta_x,Delta_y,U_0,u_x_0,v_x_0,u_x_L,v_x_L,u_y_0,v_y_0,u_y_H,v_y_H,psi_x_0,psi_x_L,psi_y_0,psi_y_H,beta,Ws,Wv,nu,ERR,XLength,YHeight,h;
int imax,jmax;
/////////////////////////////////User inputs//////////////////////////////////////////
L=1;
h=L;
H=L;
XLength=4*L;YHeight=h+H;imax=201;jmax=101;
U_0=1;
Delta_x=XLength/(imax-1);
Delta_y=YHeight/(jmax-1);
beta=Delta_x/Delta_y;
Ws=1.71;
Wv=1.71;
nu=1;
ERR=1e-5;
double psi_old[imax+2][jmax+2],omega_old[imax+2][jmax+2],psi[imax+2][jmax+2],omega[imax+2][jmax+2],u[imax+2][jmax+2],v[imax+2][jmax+2];
double CXP,CXO,CXM,CYP,CYO,CYM,Cell_Rex,Cell_Rey,error=1,counter=0;

//////////////////////////Initialize psi and omega values////////////////////////////
for(int i=0;i<=imax+1;i++)
	{
		for(int j=0;j<=jmax+1;j++)
		{
			u[i][j]=0;
			v[i][j]=0;
			psi_old[i][j]=0;
			psi[i][j]=0;
			omega[i][j]=0;
			omega_old[i][j]=0;
		}
		}
/////////////////////////////////Psi Boundary conditions/////////////////////////////////
for(int i=1;i<=imax;i++)
	{
			psi_old[i][jmax]=U_0*H;
			psi[i][jmax]=U_0*H;
		}	
for(int i=1;i<=L/Delta_y+1;i++)
	{
			int k=L/Delta_y+1;
			psi_old[i][k]=0;
			psi[i][k]=0;
		}
for(int i=L/Delta_y+1;i<=imax;i++)
	{
			psi_old[i][1]=0;
			psi[i][1]=0;
		}
for(int j=1;j<=L/Delta_y+1;j++)
	{
			int k=L/Delta_x+1;
			psi_old[k][j]=0;
			psi[k][j]=0;
		}
for(int j=L/Delta_y+1;j<=jmax;j++)
	{
			double y=Delta_y*(j-(L/Delta_y+1));
			psi_old[1][j]=3.0/L*U_0*pow(y,2)-2/pow(L,2)*U_0*pow(y,3);
			psi[1][j]=3.0/L*U_0*pow(y,2)-2/pow(L,2)*U_0*pow(y,3);
			u[1][j]=6.0/L*U_0*y-6/pow(L,2)*U_0*pow(y,2);

		}
for(int j=1;j<=jmax;j++)
	{
			//double y=Delta_y*(j-1);
			//psi_old[imax][j]=3.0/4.0/L*U_0*pow(y,2)-1.0/4.0/pow(L,2)*U_0*pow(y,3);
			//psi[imax][j]=3.0/4.0/L*U_0*pow(y,2)-1.0/4.0/pow(L,2)*U_0*pow(y,3);
			//u[imax][j]=3.0/2.0/L*U_0*y-3.0/4.0/pow(L,2)*U_0*pow(y,2);

	}
////////////////////////////////////////////////////////////////////////////////////

for(int iteration=1;iteration<=100000 && error>ERR ;iteration++)
{

	counter=counter+1;
	for(int j=1;j<=jmax;j++)
	{
			//double y=Delta_y*(j-1);
			psi_old[imax][j]=1.0/2.0*(psi_old[imax-3][j]-4*psi_old[imax-2][j]+5*psi_old[imax-1][j]);
			psi[imax][j]=1.0/2.0*(psi[imax-3][j]-4*psi[imax-2][j]+5*psi[imax-1][j]);
			//if(j>=2 && j<=(jmax-1))
			u[imax][j]=1.0/3.0*(-u[imax-2][j]+4*u[imax-1][j]);
	}
//////////////////////////Updating psi Values////////////////////////////
	for(int i=2;i<=(imax-1);i++)
	{
		for(int j=2;j<=(jmax-1);j++)
		{
			if(i<=(L/Delta_x+1) && j<=(L/Delta_y+1))
			{
				continue;
			}
			psi[i][j]=(1-Ws)*psi_old[i][j]+Ws/(2*(1+pow(beta,2)))*(psi[i-1][j]+psi[i+1][j]+pow(beta,2)*(psi[i][j-1]+psi[i][j+1])+pow(Delta_x,2)*omega_old[i][j]);
		}
	}
				for(int j=jmax;j>=1;j--)
	{
		for(int i=1;i<=imax;i++)
		{
			//cout<<setw(8)<<fixed<<setprecision(2)<< psi_old[i][j];
		}
		//cout<<endl;
	}
	//cout<<psi_old[imax][jmax]<<endl;
//////////////////////////Updating Omega boundary values////////////////////////////
for(int i=1;i<=imax;i++)
	{
			omega[i][jmax]=2.0/pow(Delta_y,2)*(U_0*H-psi[i][jmax-1]);
			omega_old[i][jmax]=2.0/pow(Delta_y,2)*(U_0*H-psi[i][jmax-1]);
		}	
for(int i=1;i<=L/Delta_y+1;i++)
	{
			int k=L/Delta_y+1;
			omega[i][k]=-2.0/pow(Delta_y,2)*(psi[i][k+1]);
			omega_old[i][k]=-2.0/pow(Delta_y,2)*(psi[i][k+1]);
		}
for(int i=L/Delta_y+1;i<=imax;i++)
	{

			omega[i][1]=-2.0/pow(Delta_y,2)*psi[i][2];
			omega_old[i][1]=-2.0/pow(Delta_y,2)*psi[i][2];
		}
for(int j=1;j<=L/Delta_y+1;j++)
	{
			int k=L/Delta_x+1;
			omega[k][j]=-2.0/pow(Delta_x,2)*(psi[k+1][j]);
			omega_old[k][j]=-2.0/pow(Delta_x,2)*(psi[k+1][j]);
		}
for(int j=L/Delta_y+1;j<=jmax;j++)
	{
		if(j=1)
		{
		//omega_old[1][j]=-1.0/(2.0*Delta_y)*(-3*u[1][j]+4*u[1][j+1]-u[1][j+2]);
		//omega[1][j]=-1.0/(2.0*Delta_y)*(-3*u[1][j]+4*u[1][j+1]-u[1][j+2]);
		omega[1][j]=2/pow(Delta_x,2)*(psi[1][j]-psi[2][j])-1/pow(Delta_y,2)*(psi[1][j+2]-2*psi[1][j+1]+psi[1][j]);
		omega_old[1][j]=2/pow(Delta_x,2)*(psi[1][j]-psi[2][j])-1/pow(Delta_y,2)*(psi[1][j+2]-2*psi[1][j+1]+psi[1][j]);
		}
		if(j=jmax)
		{
		//omega_old[1][j]=-1.0/(2.0*Delta_y)*(3*u[1][j]-4*u[1][j-1]+u[1][j-2]);
		//omega[1][j]=-1.0/(2.0*Delta_y)*(3*u[1][j]-4*u[1][j-1]+u[1][j-2]);
		omega[1][j]=2/pow(Delta_x,2)*(psi[1][j]-psi[2][j])-1/pow(Delta_y,2)*(psi[1][j-2]-2*psi[1][j-1]+psi[1][j]);
		omega_old[1][j]=2/pow(Delta_x,2)*(psi[1][j]-psi[2][j])-1/pow(Delta_y,2)*(psi[1][j-2]-2*psi[1][j-1]+psi[1][j]);
		}
		if(j>1 && j<jmax)
		{
		//omega_old[1][j]=-1.0/(2.0*Delta_y)*(u[1][j+1]-u[1][j-1]);
		//omega[1][j]=-1.0/(2.0*Delta_y)*(u[1][j+1]-u[1][j-1]);
		omega[1][j]=2/pow(Delta_x,2)*(psi[1][j]-psi[2][j])-1/pow(Delta_y,2)*(psi[1][j+1]-2*psi[1][j]+psi[1][j-1]);
		omega_old[1][j]=2/pow(Delta_x,2)*(psi[1][j]-psi[2][j])-1/pow(Delta_y,2)*(psi[1][j+1]-2*psi[1][j]+psi[1][j-1]);
		
		}

		}
for(int j=1;j<=jmax;j++)
	{
		if(j=1)
		{
		omega_old[imax][j]=2/pow(Delta_x,2)*(psi[imax][j]-psi[imax-1][j])-1/pow(Delta_y,2)*(psi[imax][j+2]-2*psi[imax][j+1]+psi[imax][j]);
		omega[imax][j]=2/pow(Delta_x,2)*(psi[imax][j]-psi[imax-1][j])-1/pow(Delta_y,2)*(psi[imax][j+2]-2*psi[imax][j+1]+psi[imax][j]);
		}
		if(j=jmax)
		{
		omega_old[imax][j]=2/pow(Delta_x,2)*(psi[imax][j]-psi[imax-1][j])-1/pow(Delta_y,2)*(psi[imax][j-2]-2*psi[imax][j-1]+psi[imax][j]);
		omega[imax][j]=2/pow(Delta_x,2)*(psi[imax][j]-psi[imax-1][j])-1/pow(Delta_y,2)*(psi[imax][j-2]-2*psi[imax][j-1]+psi[imax][j]);
		}
		if(j>1 && j<jmax)
		{
		//omega_old[imax][j]=-1.0/(2.0*Delta_y)*(u[imax][j+1]-u[imax][j-1]);
		//omega[imax][j]=-1.0/(2.0*Delta_y)*(u[imax][j+1]-u[imax][j-1]);
		omega_old[imax][j]=2/pow(Delta_x,2)*(psi[imax][j]-psi[imax-1][j])-1/pow(Delta_y,2)*(psi[imax][j+1]-2*psi[imax][j]+psi[imax][j-1]);
		omega[imax][j]=2/pow(Delta_x,2)*(psi[imax][j]-psi[imax-1][j])-1/pow(Delta_y,2)*(psi[imax][j+1]-2*psi[imax][j]+psi[imax][j-1]);
		}


	}
//////////////////////////Updating Omega interior values////////////////////////////	
		for(int i=2;i<=(imax-1);i++)
		{
		for(int j=2;j<=(jmax-1);j++)
		{
			if(i<=(L/Delta_x+1) && j<=(L/Delta_y+1))
			{
				continue;
			}
			u[i][j]=(psi[i][j+1]-psi[i][j-1])/(2*Delta_y);
			v[i][j]=-(psi[i+1][j]-psi[i-1][j])/(2*Delta_x);
			Cell_Rex=abs(u[i][j]*Delta_x/nu);
			Cell_Rey=abs(v[i][j]*Delta_y/nu);
			if (Cell_Rex<=2)
			{
				CXP=1.0/2;
				CXO=0;
				CXM=-1.0/2;
			}
			if (Cell_Rey<=2)
			{
				CYP=1.0/2;
				CYO=0;
				CYM=-1.0/2;					
			}
			if (Cell_Rex>2)
			{
				if (u[i][j]>0) { 
				CXP=0;
				CXO=1;
				CXM=-1.0;
				}
				else
				{
				CXP=1.0;
				CXO=-1;
				CXM=0;					
				}
			}
			if (Cell_Rey>2)
			{
				if (v[i][j]>0) { 
				CYP=0;
				CYO=1;
				CYM=-1.0;
				}
				else
				{
				CYP=1.0;
				CYO=-1;
				CYM=0;					
				}
			}
			omega[i][j]=(1-Wv)*omega_old[i][j]
			+Wv/(2*(1+pow(beta,2))+u[i][j]*Delta_x/(nu)*CXO+beta*v[i][j]*Delta_x/(nu)*CYO)*(
			(1-u[i][j]*Delta_x/(nu)*CXM)*omega[i-1][j]+
			(1-u[i][j]*Delta_x/(nu)*CXP)*omega[i+1][j]+
			(pow(beta,2)-beta*v[i][j]*Delta_x/(nu)*CYM)*omega[i][j-1]+
			(pow(beta,2)-beta*v[i][j]*Delta_x/(nu)*CYP)*omega[i][j+1]);
		}
		}
//////////////////////////Calculate error////////////////////////////	

			//cout<<endl;
			 error=0;
		for(int j=2;j<=(jmax-1);j++)
		{
		for(int i=2;i<=(imax-1);i++)
		{
		if (abs(omega_old[i][j]-omega[i][j])>error)
		{
			error=abs(omega_old[i][j]-omega[i][j]);
		}
		if (abs(psi_old[i][j]-psi[i][j])>error)
		{
			error=abs(psi_old[i][j]-psi[i][j]);
		}
		}
		}
		cout<<"error equals to:   "<<error<<endl;
	 
			for(int j=0;j<=(jmax+1);j++)
		{
		for(int i=0;i<=(imax+1);i++)
		{
			psi_old[i][j]=psi[i][j];
			omega_old[i][j]=omega[i][j];
		}
		}

}
for(int j=jmax;j<=L/Delta_y;j--)
{
	for(int i=1;i<=L/Delta_x;i++)
		{
			u[i][j]=0;
			v[i][j]=0;
			psi_old[i][j]=0;
			psi[i][j]=0;
		}
	}
	if (error<10)
	cout<<"Converged at    "<<counter<<"     iterations.";
	if (error>10)
	cout<<"Diverged";


//////////////////////////Save psi-omega-velocity////////////////////////////				
  ofstream myfile;
  myfile.open ("omega.txt");
  		for(int j=0;j<=(jmax+1);j++)
	{
		for(int i=0;i<=(imax+1);i++)
		{
			myfile << setw(15)<<fixed<<setprecision(10)<< omega[i][j];
		}
		myfile <<endl;
	}
  myfile.close();
	
  myfile.open ("psi.txt");
  		for(int j=0;j<=(jmax+1);j++)
	{
		for(int i=0;i<=(imax+1);i++)
		{
			myfile << setw(15)<<fixed<<setprecision(10)<< psi[i][j];
		}
		myfile <<endl;
	}
  myfile.close();
  myfile.open ("u.txt");
  		for(int j=0;j<=(jmax+1);j++)
	{
		for(int i=0;i<=(imax+1);i++)
		{
			myfile << setw(15)<<fixed<<setprecision(10)<< u[i][j];
		}
		myfile <<endl;
	}
  myfile.close();
	
}	

