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
double L,H,Delta_x,Delta_y,U_0,u_x_0,v_x_0,u_x_L,v_x_L,u_y_0,v_y_0,u_y_H,v_y_H,psi_x_0,psi_x_L,psi_y_0,psi_y_H,beta,Ws,Wv,nu,ERR;
int imax,jmax;
/////////////////////////////////User inputs//////////////////////////////////////////
L=1;H=1;imax=57;jmax=57;
U_0=1;
Delta_x=L/(imax-1);
Delta_y=H/(jmax-1);
beta=Delta_x/Delta_y;
Ws=1.1;
Wv=1.1;
nu=0.01;
ERR=1e-5;
/////////////////////////////////Boundary conditions/////////////////////////////////
u_x_0 =0;v_x_0 =0;
u_x_L =0;v_x_L =0;
u_y_0 =0;v_y_0 =0;
u_y_H = U_0;v_y_H =0;
psi_x_0 =0;
psi_x_L =0;	
psi_y_0 =0;	
psi_y_H =0;	
//////////////////////////Initialize psi and omega values////////////////////////////

double psi_old[imax+2][jmax+2],omega_old[imax+2][jmax+2],psi[imax+2][jmax+2],omega[imax+2][jmax+2],u[imax+2][jmax+2],v[imax+2][jmax+2],p[imax+2][jmax+2];
double CXP,CXO,CXM,CYP,CYO,CYM,Cell_Rex,Cell_Rey,error=1,counter=0;
for(int i=0;i<=(imax+1);i++)
	{
		for(int j=0;j<=(jmax+1);j++)
		{
			u[i][j]=0;
			v[i][j]=0;
			p[i][j]=0;
			psi_old[i][j]=0;
			psi[i][j]=0;
			omega_old[i][j]=-double(j)/jmax*(2*U_0/Delta_y);
			omega[i][j]=-double(j)/jmax*(2*U_0/Delta_y);
		}
		u[i][jmax]=U_0;
	}
for(int iteration=1;iteration<=60000 && error>ERR && error<1e6 ;iteration++)
{
	counter=counter+1;

//////////////////////////Updating psi Values////////////////////////////
	for(int i=2;i<=(imax-1);i++)
	{
		for(int j=2;j<=(jmax-1);j++)
		{
			psi[i][j]=(1-Ws)*psi_old[i][j]+Ws/(2*(1+pow(beta,2)))*(psi[i-1][j]+psi[i+1][j]+pow(beta,2)*(psi[i][j-1]+psi[i][j+1])+pow(Delta_x,2)*omega_old[i][j]);
		}
	}
//////////////////////////Updating Omega boundary values////////////////////////////
	for(int i=2;i<=(imax-1);i++)	
		{
			omega[i][1]=-2*psi[i][2]/(pow(Delta_y,2));
			omega[i][jmax]=-2*(psi[i][jmax-1]+U_0*Delta_y)/(pow(Delta_y,2));
		} 
for(int j=2;j<=(jmax-1);j++)
		{
			omega[1][j]=-2*psi[2][j]/pow(Delta_x,2);
			omega[imax][j]=-2*psi[imax-1][j]/pow(Delta_x,2);
		}
//////////////////////////Updating Omega interior values////////////////////////////	
		for(int i=2;i<=(imax-1);i++)
		{
		for(int j=2;j<=(jmax-1);j++)
		{
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

			cout<<endl;
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
		cout<<"psi-omega calculation error:    "<<error<<endl;
	 
			for(int j=0;j<=(jmax+1);j++)
		{
		for(int i=0;i<=(imax+1);i++)
		{
			psi_old[i][j]=psi[i][j];
			omega_old[i][j]=omega[i][j];
		}
		}

}
	if (error<10)
	cout<<"Converged at    "<<counter<<"     iterations."<<endl;
	if (error>10)
	cout<<"Diverged"<<endl;
	
	
//////////////////////////Calculate Pressure in flow field////////////////////////////	
	double rho=1;
	double mu=rho*nu;	
	p[1][1]=1;
	p[2][1]=p[1][1]-mu*beta/2.0*(-3*omega[2][1]+4*omega[2][2]-omega[2][3]);
	p[1][2]=p[1][1]+mu/beta/2.0*(-3*omega[1][2]+4*omega[2][2]-omega[3][2]);
	for(int j=2;j<=(jmax-1);j++)
	{
		p[1][j+1]=p[1][j-1]+mu/beta*(-3*omega[1][j]+4*omega[2][j]-omega[3][j]);
	}
	for(int i=2;i<=(imax-1);i++)
	{
		p[i+1][1]=p[i-1][1]-mu*beta*(-3*omega[i][1]+4*omega[i][2]-omega[i][3]);
	}	
	p[2][jmax]=p[1][jmax]-mu*beta/2.0*(omega[2][jmax-2]-4*omega[2][jmax-1]+3*omega[2][jmax]);
	p[imax][2]=p[imax][1]+mu/beta/2.0*(omega[imax-2][2]-4*omega[imax-1][2]+3*omega[imax][2]);
	for(int j=2;j<=(jmax-1);j++)
	{
		p[imax][j+1]=p[imax][j-1]+mu/beta*(omega[imax-2][j]-4*omega[imax-1][j]+3*omega[imax][j]);
	}
	for(int i=2;i<=(imax-1);i++)
	{
		p[i+1][jmax]=p[i-1][jmax]-mu*beta*(omega[i][jmax-2]-4*omega[i][jmax-1]+3*omega[i][jmax]);
	}
	error=1;
			for(int kk=1;kk<=10000 && error>=1e-5;kk++)
			{
	for (int j=2;j<=(jmax-1);j++)
	{
		for (int i=2;i<=(imax-1);i++)
		{
			double p_new;
			p_new=1.0/4.0*(p[i-1][j]+p[i][j-1]+p[i+1][j]+p[i][j+1])-2*rho*pow(Delta_x,2)*(1/pow(Delta_x,2)*(psi[i-1][j]-2*psi[i][j]+psi[i+1][j])
			*1/pow(Delta_y,2)*(psi[i][j-1]-2*psi[i][j]+psi[i][j+1])-pow((1/(2*Delta_y)*1/(2*Delta_x)*(psi[i+1][j+1]-psi[i-1][j+1]-psi[i+1][j-1]+psi[i-1][j-1])),2));
			if (abs(p_new-p[i][j])>= ERR);
			error=abs(p_new-p[i][j]);
			p[i][j]=p_new;
		}
	}
	cout<<"pressure calculation error:    "<<error<<endl;
	}
//////////////////////////Save psi-omega-velocity-pressure////////////////////////////				
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
  myfile.open ("pressure.txt");
  		for(int j=1;j<=(jmax);j++)
	{
		for(int i=1;i<=(imax);i++)
		{
			myfile << setw(15)<<fixed<<setprecision(10)<< p[i][j];
		}
		myfile <<endl;
	}
  myfile.close();	
	
	
}

