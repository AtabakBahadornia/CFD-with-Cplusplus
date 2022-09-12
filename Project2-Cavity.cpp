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
double L,H,Delta_x,Delta_y,U_0,u_x_0,v_x_0,u_x_L,v_x_L,u_y_0,v_y_0,u_y_H,v_y_H,psi_x_0,psi_x_L,psi_y_0,psi_y_H,beta,nu,ContinuityERROR,XMOM_ERROR,YMOM_ERROR,Re;
double FeU,FwU,FsU,FnU,DeU,DwU,DnU,DsU,FeV,FwV,FsV,FnV,DeV,DwV,DnV,DsV,Wu,Wv,Wp;
double Due,Dvn,Dvp,Dup;
int imax,jmax;
/////////////////////////////////User inputs//////////////////////////////////////////
L=1;H=1;imax=57;jmax=57;
U_0=1;
Delta_x=L/(imax-1);
Delta_y=H/(jmax-1);
Re=100;
Wu=0.5;
Wv=0.5;
Wp=0.1;
ContinuityERROR=1e-6;
XMOM_ERROR=1e-6;
YMOM_ERROR=1e-6;
/////////////////////////////////Boundary conditions/////////////////////////////////
u_x_0 =0;v_x_0 =0;
u_x_L =0;v_x_L =0;
u_y_0 =0;v_y_0 =0;
u_y_H = U_0;v_y_H =0;
//////////////////////////Initialize u-velocity v-velocity and pressure values////////////////////////////
double psi_old[imax+2][jmax+2],omega_old[imax+2][jmax+2],psi[imax+2][jmax+2],omega[imax+2][jmax+2],u[imax+2][jmax+2],P[imax+2][jmax+2]
,v[imax+2][jmax+2],Coeff_v[imax+2][jmax+2],R_v[imax+2],temp[imax+2],Coeff_u[imax+2][jmax+2],R_u[imax+2],
Ae[imax+2][jmax+2],Aw[imax+2][jmax+2],An[imax+2][jmax+2],As[imax+2][jmax+2],Ap[imax+2][jmax+2],Coeff_Pc[imax+2][jmax+2],R_Pc[imax+2],
AeU[imax+2][jmax+2],AwU[imax+2][jmax+2],AnU[imax+2][jmax+2],AsU[imax+2][jmax+2],ApU[imax+2][jmax+2],PPU[imax+2][jmax+2],
AeV[imax+2][jmax+2],AwV[imax+2][jmax+2],AnV[imax+2][jmax+2],AsV[imax+2][jmax+2],ApV[imax+2][jmax+2],PPV[imax+2][jmax+2],b[imax+2][jmax+2]
,Pc[imax+2][jmax+2],d_UP[imax+2][jmax+2],d_VP[imax+2][jmax+2],d_UE[imax+2][jmax+2],d_VN[imax+2][jmax+2],ustar[imax+2][jmax+2],vstar[imax+2][jmax+2],Pstar[imax+2][jmax+2];
double CXP,CXO,CXM,CYP,CYO,CYM,Cell_Rex,Cell_Rey,error=1,counter=0;
for(int i=0;i<=(imax+1);i++)
	{
		for(int j=0;j<=(jmax+1);j++)
		{
			u[i][j]=0;
			ustar[i][j]=0;
			v[i][j]=0;
			vstar[i][j]=0;
			P[i][j]=0;
			Pstar[i][j]=0;
			Pc[i][j]=0;
			Coeff_v[i][j]=0;
			Coeff_u[i][j]=0;
			R_u[i]=0;
			R_v[i]=0;
		}
		u[i][jmax]=U_0;
		ustar[i][jmax]=U_0;
	}
	
//////////////////////////Main loop iteration////////////////////////////

double source=1;
double errorXmomentum=1;
double errorYmomentum=1;
for(int counter=1;counter<=6000 && (source>ContinuityERROR || errorXmomentum>XMOM_ERROR || errorYmomentum>YMOM_ERROR);counter++)
{

	source=0;
	errorXmomentum=0;
	errorYmomentum=0;
//////////////////////////Updating U Values////////////////////////////

for(int iteration=1;iteration<=1;iteration++)
{

	for(int j=2;j<=(jmax-1);j++)
	{
			
	
			for(int i=3;i<=(imax-1);i++)
		{
			FeU=0.5*(u[i][j]+u[i+1][j])*Delta_y;
			FwU=0.5*(u[i][j]+u[i-1][j])*Delta_y;
			FnU=0.5*(v[i][j+1]+v[i-1][j+1])*Delta_x;
			FsU=0.5*(v[i][j]+v[i-1][j])*Delta_x;
			DeU=Delta_y/(Re*Delta_x);
			DwU=Delta_y/(Re*Delta_x);
			DnU=Delta_x/(Re*Delta_y);
			DsU=Delta_x/(Re*Delta_y);
			if (j==2)
			{
				DsU=2.0*DsU;
			}
			if (j==jmax-1)
			{
				DnU=2.0*DnU;
			}
		//Pe=F/D
		AeU[i][j]=DeU*max(0.0,(1.0-0.5*abs(FeU/DeU)))+max(-FeU,0.0);
		AwU[i][j]=DwU*max(0.0,(1.0-0.5*abs(FwU/DwU)))+max(FwU,0.0);
		AnU[i][j]=DnU*max(0.0,(1.0-0.5*abs(FnU/DnU)))+max(-FnU,0.0);
		AsU[i][j]=DsU*max(0.0,(1.0-0.5*abs(FsU/DsU)))+max(FsU,0.0);
		ApU[i][j]=AeU[i][j]+AwU[i][j]+AnU[i][j]+AsU[i][j]+(FeU+FnU-FwU-FsU);
		PPU[i][j]=(Pstar[i-1][j]-Pstar[i][j])*Delta_y;
		d_UP[i][j]=Delta_y/ApU[i][j];
		d_UE[i][j]=Delta_y/AeU[i][j];
		errorXmomentum=max(errorXmomentum,abs(ustar[i][j]-u[i][j]));		

			Coeff_u[i-2][i-1]=-AwU[i][j];
			Coeff_u[i-1][i-1]=ApU[i][j];
			Coeff_u[i][i-1]=-AeU[i][j];
			R_u[i-1]=PPU[i][j]+AsU[i][j]*u[i][j-1]+AnU[i][j]*u[i][j+1];	
			}
			R_u[2]=R_u[2]+AwU[3][j]*u[2][j];
			R_u[imax-2]=R_u[imax-2]+AeU[imax][j]*u[imax][j];

///THOMAS algorithm//
			for(int k=2;k<=(imax-2);k++)
			{
				Coeff_u[k+1][k]=Coeff_u[k+1][k]/Coeff_u[k][k];
				R_u[k]=R_u[k]/Coeff_u[k][k];
				Coeff_u[k][k]=Coeff_u[k][k]/Coeff_u[k][k];
				Coeff_u[k+1][k+1]=Coeff_u[k+1][k+1]-Coeff_u[k+1][k]*Coeff_u[k][k+1];
				R_u[k+1]=R_u[k+1]-R_u[k]*Coeff_u[k][k+1];
				Coeff_u[k][k+1]=Coeff_u[k][k+1]-Coeff_u[k][k]*Coeff_u[k][k+1];
			}
			error=0;
			ustar[imax-1][j]=R_u[imax-2];
			for(int k=imax-2;k>=3;k--)
			{
				error=max(error,(u[k][j]-R_u[k-1]+Coeff_u[k][k-1]*u[k+1][j]));
				ustar[k][j]=R_u[k-1]-Coeff_u[k][k-1]*ustar[k+1][j];
			}
		}
	}
	

	
//////////////////////////Updating V values////////////////////////////	
	for(int iteration=1;iteration<=1;iteration++)
{
	for(int i=2;i<=(imax-1);i++)
	{
			for(int j=3;j<=(jmax-1);j++)
		{
			FeV=0.5*(u[i+1][j]+u[i+1][j-1])*Delta_y;
			FwV=0.5*(u[i][j]+u[i][j-1])*Delta_y;
			FnV=0.5*(v[i][j]+v[i][j+1])*Delta_x;
			FsV=0.5*(v[i][j]+v[i][j-1])*Delta_x;
			DeV=Delta_y/(Re*Delta_x);
			DwV=Delta_y/(Re*Delta_x);
			DnV=Delta_x/(Re*Delta_y);
			DsV=Delta_x/(Re*Delta_y);
			if (i==2)
			{
				DwV=2.0*DwV;
			}
			if (i==jmax-1)
			{
				DeV=2.0*DeV;
			}
		/////Pe=F/D////
		AeV[i][j]=DeV*max(0.0,(1.0-0.5*abs(FeV/DeV)))+max(-FeV,0.0);
		AwV[i][j]=DwV*max(0.0,(1.0-0.5*abs(FwV/DwV)))+max(FwV,0.0);
		AnV[i][j]=DnV*max(0.0,(1.0-0.5*abs(FnV/DnV)))+max(-FnV,0.0);
		AsV[i][j]=DsV*max(0.0,(1.0-0.5*abs(FsV/DsV)))+max(FsV,0.0);
		ApV[i][j]=AeV[i][j]+AwV[i][j]+AnV[i][j]+AsV[i][j]+(FeV+FnV-FwV-FsV);
		PPV[i][j]=(Pstar[i][j-1]-Pstar[i][j])*Delta_x;
		d_VP[i][j]=Delta_x/ApV[i][j];
		d_VN[i][j]=Delta_y/AnV[i][j];
		errorYmomentum=max(errorYmomentum,abs(vstar[i][j]-v[i][j]));
			Coeff_v[j-2][j-1]=-AsV[i][j];
			Coeff_v[j-1][j-1]=ApV[i][j];
			Coeff_v[j][j-1]=-AnV[i][j];
			R_v[j-1]=PPV[i][j]+AwV[i][j]*v[i-1][j]+AeV[i][j]*v[i+1][j];	
			}
			R_v[2]=R_v[2]+AsV[i][3]*v[i][2];
			R_v[imax-2]=R_v[imax-2]+AnV[i][jmax]*v[i][jmax];

///THOMAS algorithm//
			for(int k=2;k<=(imax-2);k++)
			{
				Coeff_v[k+1][k]=Coeff_v[k+1][k]/Coeff_v[k][k];
				R_v[k]=R_v[k]/Coeff_v[k][k];
				Coeff_v[k][k]=Coeff_v[k][k]/Coeff_v[k][k];
				Coeff_v[k+1][k+1]=Coeff_v[k+1][k+1]-Coeff_v[k+1][k]*Coeff_v[k][k+1];
				R_v[k+1]=R_v[k+1]-R_v[k]*Coeff_v[k][k+1];
				Coeff_v[k][k+1]=Coeff_v[k][k+1]-Coeff_v[k][k]*Coeff_v[k][k+1];
			}
			error=0;
			vstar[i][jmax-1]=R_v[jmax-2];
			for(int k=jmax-2;k>=3;k--)
			{
				error=max(error,(v[i][k]-R_v[k-1]+Coeff_v[k][k-1]*vstar[i][k+1]));
				vstar[i][k]=R_v[k-1]-Coeff_v[k][k-1]*vstar[i][k+1];
			}
		}
	}


source=0;
//////////////////////////Updating Pc values////////////////////////////
	for(int j=2;j<=(jmax-1);j++)
	{
			for(int i=2;i<=(imax-1);i++)
		{

		Ae[i][j]=d_UP[i+1][j]*Delta_y;
		Aw[i][j]=d_UP[i][j]*Delta_y;
		An[i][j]=d_VP[i][j+1]*Delta_x;
		As[i][j]=d_VP[i][j]*Delta_x;
		Ap[i][j]=Ae[i][j]+Aw[i][j]+An[i][j]+As[i][j];
		b[i][j]=(ustar[i][j]-ustar[i+1][j])*Delta_y+(vstar[i][j]-vstar[i][j+1])*Delta_x;
		source=max(source,abs(b[i][j]));
}
}
double errorPc=1;	
	for(int iteration=1;iteration<=1000 && errorPc>=1e-6;iteration++)
{
 errorPc=0;
	for(int j=2;j<=(jmax-1);j++)
	{
			for(int i=2;i<=(imax-1);i++)
		{
			Coeff_Pc[i-1][i]=-Aw[i][j];
			Coeff_Pc[i][i]=Ap[i][j];
			Coeff_Pc[i+1][i]=-Ae[i][j];
			R_Pc[i]=An[i][j]*Pc[i][j+1]+As[i][j]*Pc[i][j-1]+b[i][j];
		}
				R_Pc[2]=R_Pc[2]+Aw[2][j]*Pc[1][j];
				R_Pc[imax-1]=R_Pc[imax-1]+Ae[imax-1][j]*Pc[imax][j];
			for(int k=2;k<=(imax-1);k++)
			{
				Coeff_Pc[k+1][k]=Coeff_Pc[k+1][k]/Coeff_Pc[k][k];
				R_Pc[k]=R_Pc[k]/Coeff_Pc[k][k];
				Coeff_Pc[k][k]=Coeff_Pc[k][k]/Coeff_Pc[k][k];
				Coeff_Pc[k+1][k+1]=Coeff_Pc[k+1][k+1]-Coeff_Pc[k+1][k]*Coeff_Pc[k][k+1];
				R_Pc[k+1]=R_Pc[k+1]-R_Pc[k]*Coeff_Pc[k][k+1];
				Coeff_Pc[k][k+1]=Coeff_Pc[k][k+1]-Coeff_Pc[k][k]*Coeff_Pc[k][k+1];
			}
			errorPc=abs(Pc[imax-1][j]-R_Pc[imax-1]);
			Pc[imax-1][j]=R_Pc[imax-1];
			for(int k=imax-2;k>=2;k--)
			{
				errorPc=max(errorPc,(Pc[k][j]-R_Pc[k]+Coeff_Pc[k+1][k]*Pc[k+1][j]));
				Pc[k][j]=R_Pc[k]-Coeff_Pc[k+1][k]*Pc[k+1][j];
			}
		}	
	}
//////////////////////////Correct Pressure and velocity field values ////////////////////////////

	for(int j=2;j<=(jmax-1);j++)
	{
			
	
			for(int i=2;i<=(imax-1);i++)
		{
			P[i][j]=Pstar[i][j]+Wp*Pc[i][j];
			Pstar[i][j]=P[i][j];
		}
	}
		for(int j=2;j<=(jmax-1);j++)
	{
			
	
			for(int i=3;i<=(imax-1);i++)
		{
			u[i][j]=ustar[i][j]+d_UP[i][j]*(Pc[i-1][j]-Pc[i][j]);
		}
	}
		for(int i=2;i<=(imax-1);i++)
	{
			for(int j=3;j<=(jmax-1);j++)
		{
			v[i][j]=vstar[i][j]+d_VP[i][j]*(Pc[i][j-1]-Pc[i][j]);
		}
	}
if (counter%10==0 || (errorXmomentum<XMOM_ERROR*10 && errorYmomentum<YMOM_ERROR*10 && source<ContinuityERROR*10))
{
cout<<endl<<"..................     iteration number:   "<<counter<<"     ................."<<endl;
cout<<"X-momentum error=   "<<errorXmomentum<<endl;
cout<<"Y-momentum error=   "<<errorYmomentum<<endl;
cout<<"Continuity error=   "<<source<<endl;
}
		for(int p=1;p<=(jmax);p++)	
		{
		for(int t=1;t<=(imax);t++)
		{
			Pc[t][p]=0;
		
		}
		}

}
//////////////////////////Save psi-omega-velocity-pressure////////////////////////////				
  ofstream myfile;
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
			myfile << setw(15)<<fixed<<setprecision(10)<< P[i][j];
		}
		myfile <<endl;
	}
  myfile.close();	
}

