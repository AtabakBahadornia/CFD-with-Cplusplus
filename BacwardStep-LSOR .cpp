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
XLength=8*L;YHeight=h+H;imax=201;jmax=51;
U_0=1;
Delta_x=XLength/(imax-1);
Delta_y=YHeight/(jmax-1);
beta=Delta_x/Delta_y;
Ws=1;
Wv=1;
nu=0.002;
ERR=1e-5;
double psi_old[imax+2][jmax+2],omega_old[imax+2][jmax+2],psi[imax+2][jmax+2],omega[imax+2][jmax+2],u[imax+2][jmax+2],v[imax+2][jmax+2],p[imax+2][jmax+2]
,Coeff[imax+2][jmax+2],R[imax+2],temp[imax+2],Coeff_psi[imax+2][jmax+2],R_psi[imax];
double CXP,CXO,CXM,CYP,CYO,CYM,Cell_Rex,Cell_Rey,error=1,counter=0;
//////////////////////////Initialize psi and omega values////////////////////////////
for(int i=0;i<=imax+1;i++)
	{
		for(int j=0;j<=jmax+1;j++)
		{
			u[i][j]=0;
			v[i][j]=0;
			p[i][j]=0;
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
			double y=Delta_y*(j-1);
			psi_old[imax][j]=3.0/4.0/L*U_0*pow(y,2)-1.0/4.0/pow(L,2)*U_0*pow(y,3);
			psi[imax][j]=3.0/4.0/L*U_0*pow(y,2)-1.0/4.0/pow(L,2)*U_0*pow(y,3);

	}
////////////////////////////////////////////////////////////////////////////////////

for(int iteration=1;iteration<=20000 && error>ERR ;iteration++)
{
//////////////////////////Updating psi Values////////////////////////////
	counter=counter+1;
	for(int j=2;j<=(jmax-1);j++)
	{
		if(j>(L/Delta_x+1))	
	{
			for(int i=2;i<=(imax-1);i++)
		{
			Coeff_psi[i-1][i]=Ws;
			Coeff_psi[i][i]=-2*(1+pow(beta,2));
			Coeff_psi[i+1][i]=Ws;
			R_psi[i]=-(1-Ws)*2*(1+pow(beta,2))*psi_old[i][j]-Ws*pow(beta,2)*(psi_old[i][j+1]+psi[i][j-1])-pow(Delta_x,2)*omega_old[i][j];
			
		}	
				R_psi[2]=R_psi[2]-Ws*psi[1][j];
				R_psi[imax-1]=R_psi[imax-1]-Ws*psi[imax][j];
///THOMAS algorithm//
			for(int k=2;k<=(imax-1);k++)
			{
				Coeff_psi[k+1][k]=Coeff_psi[k+1][k]/Coeff_psi[k][k];
				R_psi[k]=R_psi[k]/Coeff_psi[k][k];
				Coeff_psi[k][k]=Coeff_psi[k][k]/Coeff_psi[k][k];
				Coeff_psi[k+1][k+1]=Coeff_psi[k+1][k+1]-Coeff_psi[k+1][k]*Coeff_psi[k][k+1];
				R_psi[k+1]=R_psi[k+1]-R_psi[k]*Coeff_psi[k][k+1];
				Coeff_psi[k][k+1]=Coeff_psi[k][k+1]-Coeff_psi[k][k]*Coeff_psi[k][k+1];
			}
			psi[imax-1][j]=R_psi[imax-1];
			for(int k=imax-2;k>=2;k--)
			{
				psi[k][j]=R_psi[k]-Coeff_psi[k+1][k]*psi[k+1][j];
			}
		}
	if(j<=(L/Delta_x+1))
	{
			for(int i=L/Delta_x+2;i<=(imax-1);i++)
		{
			Coeff_psi[i-1][i]=Ws;
			Coeff_psi[i][i]=-2*(1+pow(beta,2));
			Coeff_psi[i+1][i]=Ws;
			R_psi[i]=-(1-Ws)*2*(1+pow(beta,2))*psi_old[i][j]-Ws*pow(beta,2)*(psi_old[i][j+1]+psi[i][j-1])-pow(Delta_x,2)*omega_old[i][j];
			
		}	
		int xx=L/Delta_x+2;
				R_psi[xx]=R_psi[xx]-Ws*psi[xx-1][j];
				R_psi[imax-1]=R_psi[imax-1]-Ws*psi[imax][j];
///THOMAS algorithm//
			for(int k=L/Delta_x+2;k<=(imax-1);k++)
			{
				Coeff_psi[k+1][k]=Coeff_psi[k+1][k]/Coeff_psi[k][k];
				R_psi[k]=R_psi[k]/Coeff_psi[k][k];
				Coeff_psi[k][k]=Coeff_psi[k][k]/Coeff_psi[k][k];
				Coeff_psi[k+1][k+1]=Coeff_psi[k+1][k+1]-Coeff_psi[k+1][k]*Coeff_psi[k][k+1];
				R_psi[k+1]=R_psi[k+1]-R_psi[k]*Coeff_psi[k][k+1];
				Coeff_psi[k][k+1]=Coeff_psi[k][k+1]-Coeff_psi[k][k]*Coeff_psi[k][k+1];
			}
			psi[imax-1][j]=R_psi[imax-1];
			for(int k=imax-2;k>=L/Delta_x+2;k--)
			{
				psi[k][j]=R_psi[k]-Coeff_psi[k+1][k]*psi[k+1][j];
			}
	}
	
	
	}
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
		omega[1][j]=2/pow(Delta_x,2)*(psi[1][j]-psi[2][j])-1/pow(Delta_y,2)*(psi[1][j+2]-2*psi[1][j+1]+psi[1][j]);
		omega_old[1][j]=2/pow(Delta_x,2)*(psi[1][j]-psi[2][j])-1/pow(Delta_y,2)*(psi[1][j+2]-2*psi[1][j+1]+psi[1][j]);
		}
		if(j=jmax)
		{
		omega[1][j]=2/pow(Delta_x,2)*(psi[1][j]-psi[2][j])-1/pow(Delta_y,2)*(psi[1][j-2]-2*psi[1][j-1]+psi[1][j]);
		omega_old[1][j]=2/pow(Delta_x,2)*(psi[1][j]-psi[2][j])-1/pow(Delta_y,2)*(psi[1][j-2]-2*psi[1][j-1]+psi[1][j]);
		}
		if(j>1 && j<jmax)
		{
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
		omega_old[imax][j]=2/pow(Delta_x,2)*(psi[imax][j]-psi[imax-1][j])-1/pow(Delta_y,2)*(psi[imax][j+1]-2*psi[imax][j]+psi[imax][j-1]);
		omega[imax][j]=2/pow(Delta_x,2)*(psi[imax][j]-psi[imax-1][j])-1/pow(Delta_y,2)*(psi[imax][j+1]-2*psi[imax][j]+psi[imax][j-1]);
		}
	}
//////////////////////////Updating Omega interior values////////////////////////////			
		for(int j=(L/Delta_x+2);j<=(jmax-1);j++)
	{
			for(int i=2;i<=(imax-1);i++)
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
			Coeff[i-1][i]=Wv*(1-u[i][j]*Delta_x/(nu)*CXM);
			Coeff[i][i]=-(2*(1+pow(beta,2))+u[i][j]*Delta_x/(nu)*CXO+beta*v[i][j]*Delta_x/(nu)*CYO);
			Coeff[i+1][i]=Wv*(1-u[i][j]*Delta_x/(nu)*CXP);
			R[i]=-(1-Wv)*(2*(1+pow(beta,2))+u[i][j]*Delta_x/(nu)*CXO+beta*v[i][j]*Delta_x/(nu)*CYO)*omega[i][j]-
			Wv*(pow(beta,2)-beta*v[i][j]*Delta_x/(nu)*CYM)*(omega[i][j-1])-
			Wv*(pow(beta,2)-beta*v[i][j]*Delta_x/(nu)*CYP)*(omega[i][j+1]);			
		}	
		R[2]=R[2]-Wv*(1-u[2][j]*Delta_x/(nu)*CXM)*omega[1][j];
		R[imax-1]=R[imax-1]-Wv*(1-u[2][j]*Delta_x/(nu)*CXP)*omega[imax][j];
///THOMAS algorithm//
			for(int k=2;k<=(imax-1);k++)
			{
				Coeff[k+1][k]=Coeff[k+1][k]/Coeff[k][k];
				R[k]=R[k]/Coeff[k][k];
				Coeff[k][k]=Coeff[k][k]/Coeff[k][k];
				Coeff[k+1][k+1]=Coeff[k+1][k+1]-Coeff[k+1][k]*Coeff[k][k+1];
				R[k+1]=R[k+1]-R[k]*Coeff[k][k+1];
				Coeff[k][k+1]=Coeff[k][k+1]-Coeff[k][k]*Coeff[k][k+1];
			}
			omega[imax-1][j]=R[imax-1];
			for(int k=imax-2;k>=2;k--)
			{
				omega[k][j]=R[k]-Coeff[k+1][k]*omega[k+1][j];
			}
	}
	
	
	
	
	
	
for(int j=2;j<=(L/Delta_x+1);j++)
	{
			for(int i=(L/Delta_x+2);i<=(imax-1);i++)
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
			Coeff[i-1][i]=Wv*(1-u[i][j]*Delta_x/(nu)*CXM);
			Coeff[i][i]=-(2*(1+pow(beta,2))+u[i][j]*Delta_x/(nu)*CXO+beta*v[i][j]*Delta_x/(nu)*CYO);
			Coeff[i+1][i]=Wv*(1-u[i][j]*Delta_x/(nu)*CXP);
			R[i]=-(1-Wv)*(2*(1+pow(beta,2))+u[i][j]*Delta_x/(nu)*CXO+beta*v[i][j]*Delta_x/(nu)*CYO)*omega[i][j]-
			Wv*(pow(beta,2)-beta*v[i][j]*Delta_x/(nu)*CYM)*(omega[i][j-1])-
			Wv*(pow(beta,2)-beta*v[i][j]*Delta_x/(nu)*CYP)*(omega[i][j+1]);			
		}	
		int xx=L/Delta_x+2;
		R[xx]=R[xx]-Wv*(1-u[xx][j]*Delta_x/(nu)*CXM)*omega[xx-1][j];
		R[imax-1]=R[imax-1]-Wv*(1-u[2][j]*Delta_x/(nu)*CXP)*omega[imax][j];
///THOMAS algorithm//
			for(int k=xx;k<=(imax-1);k++)
			{
				Coeff[k+1][k]=Coeff[k+1][k]/Coeff[k][k];
				R[k]=R[k]/Coeff[k][k];
				Coeff[k][k]=Coeff[k][k]/Coeff[k][k];
				Coeff[k+1][k+1]=Coeff[k+1][k+1]-Coeff[k+1][k]*Coeff[k][k+1];
				R[k+1]=R[k+1]-R[k]*Coeff[k][k+1];
				Coeff[k][k+1]=Coeff[k][k+1]-Coeff[k][k]*Coeff[k][k+1];
			}
			omega[imax-1][j]=R[imax-1];
			for(int k=imax-2;k>=xx;k--)
			{
				omega[k][j]=R[k]-Coeff[k+1][k]*omega[k+1][j];
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
//////////////////////////Calculate Pressure in flow field////////////////////////////	
	double rho=1;
	double mu=rho*nu;
	p[imax][1]=1;
	p[imax-1][1]=p[imax][1]+mu*beta/2.0*(-3*omega[imax][1]+4*omega[imax][2]-omega[imax][3]);
	int vv=L/Delta_x+1;
	int hh=L/Delta_y+1;
		for(int i=imax-1;i>=(L/Delta_x+2);i--)
	{
		p[i-1][1]=p[i+1][1]+mu*beta*(-3*omega[i][1]+4*omega[i][2]-omega[i][3]);
	}
	p[vv][2]=p[vv][1]+mu/beta/2.0*(-3*omega[vv][2]+4*omega[vv+1][2]-omega[vv+2][2]);
		for(int j=2;j<=(L/Delta_x);j++)
	{
		p[vv][j+1]=p[vv][j-1]+mu/beta*(-3*omega[vv][j]+4*omega[vv+1][j]-omega[vv+2][j]);
	}
	p[vv-1][hh]=p[vv][vv]-mu*beta/2.0*(-3*omega[vv][hh]+4*omega[vv][hh+1]-omega[vv][hh+2]);
		for(int i=L/Delta_x;i>=1;i--)
	{
		p[i-1][hh]=p[i+1][hh]+mu*beta*(-3*omega[i][hh]+4*omega[i][hh+1]-omega[i][hh+2]);
	for(int j=hh;j<=jmax;j++)
	{
		p[1][j]=p[1][hh];
	}
	p[2][jmax]=p[1][jmax]-mu*beta/2.0*(omega[2][jmax-2]-4*omega[2][jmax-1]+3*omega[2][jmax]);
	}
		for(int i=2;i<=imax-1;i++)
	{
		p[i+1][jmax]=p[i-1][jmax]-mu*beta*(omega[i][jmax-2]-4*omega[i][jmax-1]+3*omega[i][jmax]);
	}
		for(int j=1;j<=jmax;j++)
	{
		p[imax][j]=p[imax][1];
	}
	//p[imax-1][1]=p[imax][1]-mu*beta/2.0*(-3*omega[imax][1]+4*omega[imax][2]-omega[imax][3]);
	//p[1][2]=p[1][1]+mu/beta/2.0*(-3*omega[1][2]+4*omega[2][2]-omega[3][2]);
					for(int j=jmax;j>=1;j--)
	{
		for(int i=1;i<=imax;i++)
		{
			//cout<<setw(8)<<fixed<<setprecision(2)<< p[i][j];
		}
		//cout<<endl;
	}
	error=1;
			for(int kk=1;kk<=10000 && error>=1e-5;kk++)
			{
	for (int j=2;j<=(jmax-1);j++)
	{
		for (int i=2;i<=(imax-1);i++)
		{
			if(i<=(L/Delta_x+1) && j<=(L/Delta_y+1))
			{
				continue;
			}
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
		if (error<10)
	cout<<"Converged at    "<<counter<<"     iterations.";
	if (error>10)
	cout<<"Diverged";
	cout<<endl;
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
	 myfile.open ("v.txt");
  		for(int j=0;j<=(jmax+1);j++)
	{
		for(int i=0;i<=(imax+1);i++)
		{
			myfile << setw(15)<<fixed<<setprecision(10)<< v[i][j];
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
	


