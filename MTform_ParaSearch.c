/*============================================================================
 Name        : 	MTformation_ParaSearch.c
 Author      : 	Chieh-Ting(jimmy) Hsu: MCGILL UNIVERSITY: April 2019
 Description : 	Calculation of growing experiment to get comparable Glat, Glong
		for GMPCPP depolymerization rate from experiment
		Lateral energy only calculated once for each bond strength
		Use to test the principle of David Odde's simple model with 
		paramater compression
		This is s a new code 2016/11 for smaller compile outputs
		This is a new code 2017/09 for different onrate
		update 2019/02/02 for running of 7 parameters
		This is a new updated code to do parameter search to match 
		experimental data
 ============================================================================*/
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include "random.c"
#include "ran2.c"
#define MIN(A,B)  ((A)<(B) ? (A) : (B))
#define MAX(A,B)  ((A)>(B) ? (A) : (B))
#define pi 3.14159265358979323846264338327950
#define RAND00 genrand_real1()		//[0,1]
#define RAND01 ran2(&idum)		//[0,1)

//Physical constants
#define nPF 13				//Number of protofilament
//Simulation constants
#define N0 2				//Original # of total dimers/PF
#define N 10000				//# of total dimers/PF
#define MAXtime 600.2			//Simulation time [s]
#define prn 7				//Number of parameters in the model
#define kruns 1E2			//Number of runs to build probability distribution

/*Function for calculating Glat depending on # of dimers in a ring for Alpha-Alpha*/
double GlatCalAA(int NumbOfDim ,double AAbond)
{
	double GlatAA=0;
	GlatAA = AAbond;
	return GlatAA;
}
/*Function for calculating Glat depending on # of dimers in a ring for Beta-Beta*/
double GlatCalBB(int NumbOfDim, double BBbond)
{
	double GlatBB=0;
	GlatBB = BBbond;
	return GlatBB;
}
/*Function for calulating OnRate due too different number of neighbors*/
double OnRateCal(int (**lattice), int row, int col,double Kplus,double Ks)
{
	int Nneighbour=0;//Number of Neighbour
	if(col==0)	//Get left and right position
	{
		if (lattice[1][row]!=0)//has left neighbor
		{
			if (lattice[(nPF-1)][row-1]!=0&&(row-1>=0)&&lattice[(nPF-1)][row-2]!=0&&(row-2>=0)) 
				Nneighbour = 2;		//Same as 2 neighbors
			else if (lattice[(nPF-1)][row-1]==0&&(row-1>=0)&&lattice[(nPF-1)][row-2]!=0&&(row-2>=0))
				Nneighbour = 2; 		//Like 1.5 neighbors
			else if (lattice[(nPF-1)][row-1]==0&&(row-1>=0)&&lattice[(nPF-1)][row-2]==0&&(row-2>=0))
				Nneighbour = 1;		//Same as 1 neighbors
		}	
		else
		{
			if (lattice[(nPF-1)][row-1]!=0&&(row-1>=0)&&lattice[(nPF-1)][row-2]!=0&&(row-2>=0)) 
				Nneighbour = 1;		//Same as 1 neighbors
			else if (lattice[(nPF-1)][row-1]==0&&(row-1>=0)&&lattice[(nPF-1)][row-2]!=0&&(row-2>=0))
				Nneighbour = 1; 		//Same as 0.5 neighbors
			else if (lattice[(nPF-1)][row-1]==0&&(row-1>=0)&&lattice[(nPF-1)][row-2]==0&&(row-2>=0))
				Nneighbour = 0;		//Same as 0 neighbors
		}
	}	
	else if (col==(nPF-1))
	{
		if (lattice[col-1][row]!=0)//has right neighbor
		{
			if (lattice[0][row+2]!=0&&lattice[0][row+1]!=0) 
				Nneighbour = 2;		//Same as 2 neighbors
			else if (lattice[0][row+2]==0&&lattice[0][row+1]!=0)
				Nneighbour = 2; 		//Like 1.5 neighbors
			else if (lattice[0][row+2]==0&&lattice[0][row+1]==0)
				Nneighbour = 1;		//Same as 1 neighbors
		}	
		else
		{
			if (lattice[0][row+2]!=0&&lattice[0][row+1]!=0) 
				Nneighbour = 1;		//Same as 1 neighbors
			else if (lattice[0][row+2]==0&&lattice[0][row+1]!=0)
				Nneighbour = 1; 		//Same as 0.5 neighbors
			else if (lattice[0][row+2]==0&&lattice[0][row+1]==0)
				Nneighbour = 0;		//Same as 0 neighbors
		}
	}
	else
	{
		if(lattice[col-1][row]!= 0&&lattice[col+1][row] != 0) //2 neighbors			
			Nneighbour = 2;			//Same as 2 neighbors			
		else if(lattice[col-1][row]!= 0||lattice[col+1][row] != 0) //1 Neighbor
			Nneighbour = 1;			//Same as 1 neighbors
		else
			Nneighbour = 0;			//Same as 0 neighbors
	}
	//printf("Kplus %lf\n",Kplus);
	if (Nneighbour==0) //Return diffrate
		return Kplus;
	else
		return (Kplus/(Ks*Nneighbour));
}
/*Function for calulating Offrate due too different number of neighbors*/
double OffRateCal(int (**lattice), int row, int col,double ABbond, double Glong,double Eglat[nPF],double Kplus, double Ks,double ddG)
{
	double Energy=0;
	double Nneighbour=0;
	if (lattice[col][row-1]!=0&&(row-1>=0)) 
		Energy += Glong;
	//double ddG = 3.15;
	if(lattice[col][row]==1)//current dimer is a GTP
	{
		//printf("Energy0: %lf\n",Energy);
		if(col==0)						//At Seam
		{
			if (lattice[(int)(nPF-1)][row-1]==1&&(row-1>=0)) 	//GTP
			{ 
				Energy +=  ABbond;
				Nneighbour += 0.5;
			}
			if (lattice[(int)(nPF-1)][row-2]==1&&(row-2>=0)) 	//GTP
			{
				Energy +=  ABbond;
				Nneighbour += 0.5;
			}
			if (lattice[(int)(nPF-1)][row-1]==2&&(row-1>=0)) 	//GDP 
			{
				Energy +=  ABbond+(ddG/2);
				Nneighbour += 0.5;
			}
			if (lattice[(int)(nPF-1)][row-2]==2&&(row-2>=0)) 	//GDP
			{
				Energy +=  ABbond+(ddG/2);
				Nneighbour += 0.5;
			}
			//Lateral Bond	
			if (lattice[col+1][row]==1)			//Check only left
			{
				Energy += Eglat[1]+Eglat[1];
				Nneighbour += 1;
			}
			else if (lattice[col+1][row]==2)
			{
				Energy += Eglat[1]+Eglat[1]+ddG;
				Nneighbour += 1;
			}
		}	
		else if(col==(nPF-1))					//At Seam
		{
			if (lattice[0][row+1]==1)			//GTP
			{
				Energy +=  ABbond;	
				Nneighbour += 0.5;
			}		
			if (lattice[0][row+2]==1)			//GTP
			{
				Energy +=  ABbond;
				Nneighbour += 0.5;	
			}		
			if (lattice[0][row+1]==2) 			//GDP 
			{
				Energy +=  ABbond+(ddG/2);
				Nneighbour += 0.5;
			}
			if (lattice[0][row+2]==2) 			//GDP
			{
				Energy +=  ABbond+(ddG/2);
				Nneighbour += 0.5;
			}
			//Lateral Bond
			if (lattice[col-1][row]==1)			//Check only right
			{
				Energy += Eglat[1]+Eglat[1];
				Nneighbour += 1;
			}
			else if (lattice[col-1][row]==2)
			{
				Energy += Eglat[1]+Eglat[1]+ddG;
				Nneighbour += 1;
			}
		}	
		else
		{
			if (lattice[col+1][row]==1)			//Check left GTP
			{
				Energy +=  Eglat[1]+Eglat[1];	
				Nneighbour += 1;
			}
			else if (lattice[col+1][row]==2)			//Check left GDP
			{
				Energy +=  Eglat[1]+Eglat[1]+ddG; 
				Nneighbour += 1;
			}
			if (lattice[col-1][row]==1)			//Check right GTP
			{
				Energy += Eglat[1]+Eglat[1];	
				Nneighbour += 1;
			}
			else if (lattice[col-1][row]==2)			//Check left GDP
			{
				Energy +=  Eglat[1]+Eglat[1]+ddG;
				Nneighbour += 1;
			} 
			//printf("%d\n",tempring);	
		}
	}
	else if(lattice[col][row]==2)//GDP lateral bond is low
	{
		//printf("Energy0: %lf\n",Energy);
		if(col==0)						//At Seam
		{
			if (lattice[(int)(nPF-1)][row-1]!=0&&(row-1>=0))//Check on right
			{	
				Energy +=  ABbond+(ddG/2);
				Nneighbour += 0.5;
			}
			if (lattice[(int)(nPF-1)][row-2]!=0&&(row-2>=0)) 
			{	
				Energy +=  ABbond+(ddG/2);
				Nneighbour += 0.5;
			}
			if (lattice[col+1][row]!=0)			//Check only left
			{
				Energy += Eglat[1]+Eglat[1]+ddG;
				Nneighbour += 1;
			}
		}	
		else if(col==(nPF-1))					//At Seam
		{
			if (lattice[0][row+1]!=0)			//GTP
			{
				Energy +=  ABbond+(ddG/2);	
				Nneighbour += 0.5;	
			}	
			if (lattice[0][row+2]!=0)			//GTP
			{		 
				Energy +=  ABbond+(ddG/2);
				Nneighbour += 0.5;
			}
			//Lateral Bond
			if (lattice[col-1][row]!=0)			//Check only right
			{
				Energy += Eglat[1]+Eglat[1]+ddG;
				Nneighbour += 1;
			}
		}	
		else
		{
			if (lattice[col+1][row]!=0)			//Check left 	
			{		
				Energy +=  Eglat[1]+Eglat[1]+ddG; 
				Nneighbour += 1;
			}
			if (lattice[col-1][row]!=0)			//Check right	
			{
				Energy +=  Eglat[1]+Eglat[1]+ddG; 
				Nneighbour += 1;
			}
			//printf("%d\n",tempring);	
		}
	}
	if (Nneighbour==0)
		return (Kplus/exp(-Energy));	//Offrate
	else if (Nneighbour>0&&Nneighbour<=1)
		return ((Kplus/Ks)/exp(-Energy));	//Offrate
	else 
		return ((Kplus/(Ks*2))/exp(-Energy));	//Offrate
}
void MAINRUN(double prval[prn],double mfactor[prn],double C)
{
	//INITIALIZATION OF RANDOM NUMBER GENERATOR	
	double seed;
	long int idum=-time(NULL);
	int kk;
	for(kk=0;kk<18573;kk++)
		ran2(&idum);
	seed = floor(4294967296.0*ran2(&idum));
	//printf("seed = %lf\n",seed);
	init_genrand(seed);
	for(kk=0;kk<18573;kk++)// warm up RNG
		genrand_int32();
	//VARIABLES FOR SIMULATION
	double Glong;			//Lower limit Longitudinal bond energy/dimer [KbT] Odde.
	double Glat;			//Lower limit Lateral bond enery/dimer [KbT] Odde.
	double Hrate;			//Hydrolysis rate
	double Ks;			//On rate factor
	double ddG;			//Change in hydrolysis effect
	double Kplus;			//On rate
	int** lattice;   			//MT dimers lattice
	int PF[(int)nPF];			//Number of dimers per filament
	int GTP[(int)nPF];			//Number of GTP perfilament
	int row, col;			//indices
	int totaldimer=0;	//Keep track of total dimers in the configuration
	int totalGTP=0;		//Keep track of total number of GTP
	int potdimer=0;		//Keep track of the number of dimers not at tip
	int max,min;		//Max and min number of dimers
	lattice = (int**) malloc((int)(nPF) *sizeof(int*));  
	FILE *output;		//write to file
	double t=0;		//Keep track of time
	int temprow,tempcol;	//Save position and event
	double alpha,tau,R2,Psum;
	double Parray[(int)(2*nPF+1)];		//probability array
	double Eglat[nPF];			//for Lateral energy
	int i;
	char filename[50];// long enough to hold filename
    	for(col=0;col<nPF;col++)  
       		lattice[col] = (int*) malloc((N) *sizeof(int)); //2*N: monomer number	
	double AAbond;		//Bond between alpha tubulin
	double ABbond;		//Bond between alpha/beta tubulin [KbT] at Seam David Sept.
	double avg2,avg,var;
	int count;
	int flag=1;			//Flag for output
	double tout,delt,delt2;	//For outputing purpose	
	Glong = prval[0];
	Glat = prval[1];
	ABbond = prval[2];
	Hrate = prval[3];
	Ks = prval[4];
	Kplus = prval[5];
	ddG = prval[6];
	int tpoint = 0;
	//int zerocount =0;

	/*START ACTUAL SIMULATION*/
	AAbond=(Glat/2.0);		//Bond between alpha tubulin
	ABbond = ABbond/2.0;		//Bond at seam
	//printf("2 %lf %lf %lf %lf %lf %lf %lf\n",Glong,Glat,ABbond,Hrate*1000,Ks*1000,Kplus*1E-4,ddG);
	for(col=0;col<nPF;col++)  
	{	
		if (col==0)
			Eglat[col] = 0;
		else
			Eglat[col] =  GlatCalAA(col,AAbond);//Half of bond strength
		// printf("Glat: %lf\n",Eglat[col]);
	}
	for (kk=0;kk<kruns;kk++)
	{
		tout = 0.1;		//First output
		delt = 0.1; 		//Delta time for output
		delt2 = 0.1;		
		sprintf(filename,"MT_%1.0lf_%1.0lf_%1.0lf_%1.0lf_%1.0lf_%1.0lf_%1.0lf_%1.0lf_%d",-Glong*mfactor[0],-Glat*mfactor[1],-ABbond*2*mfactor[2],Hrate*mfactor[3],Ks*mfactor[4],Kplus*mfactor[5],ddG*mfactor[6],C*1E7,kk+1);//File name for the particular simulation
		output = fopen(filename,"a+");//Open the file
		fprintf(output,"%1.6lf %1.6lf %1.6lf\n",0.0,0.0,0.0);//check point
		for(col=0;col<nPF;col++)  
			PF[col]=GTP[col]=0;
		totaldimer=totalGTP=0;
		for(row=0;row<N;row++)
		{
			for(col=0;col<nPF;col++)
			{
				if(row<N0)//x,y position of dimers
				{
					lattice[col][row] = 1;
					PF[col]++;
					GTP[col]++;
					totaldimer++;
					totalGTP++;
				}
				else
					lattice[col][row] = 0;
			}
		}
		t=0;
		count=0;
		tpoint = 0;
		flag=1;
		//while (count<=countmax) //Gillespie Algorithm
		while (t<MAXtime&&totaldimer<=26000) //Gillespie Algorithm
		{
			alpha=Psum=0;			//For Gillespie
			for(col=0;col<nPF;col++)	
			{
				Parray[col] = OnRateCal(lattice,PF[col],col,Kplus,Ks)*C;
				alpha += Parray[col];
			}
			for(col=0;col<nPF;col++)
			{
				if (PF[col]>N0)		//Dissociation above seed
				{
					Parray[(int)(col+nPF)] = OffRateCal(lattice,(PF[col]-1),col,ABbond,Glong,Eglat,Kplus,Ks,ddG);	//Offrate				
					alpha += Parray[(int)(col+nPF)];
				}
				else
					Parray[(int)(col+nPF)] = 0;		
			}
			potdimer = 0;
			for(col=0;col<nPF;col++)
			{
				if ((GTP[col]-N0)>1)	
					potdimer += (GTP[col]-N0)-1;	//seed and tip doesn't count
			}
			//if(kk==0)
			//	printf("potdimer: %d time: %lf\n",potdimer,t);
			if(potdimer>0)			//Hydrolysis event condition meet
			{
				Parray[(int)(2*nPF)] = Hrate*potdimer;
				alpha += Parray[(int)(2*nPF)];
			}
			else
				Parray[(int)(2*nPF)] = 0;
			//for(kk=0;kk<(2*nPF+1);kk++)
			//	printf("Parray: %lf\n",Parray[kk]);
			//printf("\n");
			tau = (1/alpha)*log(1.0/RAND01);	//time step
			R2 = RAND01;				//Second random number
			for (i=0;i<(2*nPF+1);i++)			//Find event
			{
				if (i==0)
				{
					if (R2<Parray[i]/alpha)
						break;
					else 
						Psum += Parray[i];
				}
				else
				{
					if(Psum/alpha<=R2&&R2<(Psum+Parray[i])/alpha)
						break;
					else
						Psum += Parray[i];
				}
			}	
			if (i<nPF)		//Association
			{
				tempcol = i;
				temprow = PF[tempcol];
				lattice[tempcol][temprow]=1;
				PF[tempcol]++;
				GTP[tempcol]++;
				totaldimer++;
				totalGTP++;
			}
			else if (nPF<=i&&i<2*nPF)//Dissociation
			{
				tempcol = (i-nPF);	//corresponding column
				temprow = PF[tempcol]-1;
				if (lattice[tempcol][temprow]==1)
				{
					GTP[tempcol]--;
					totalGTP--;
				}
				lattice[tempcol][temprow]=0;
				PF[tempcol]--;
				totaldimer--;
			}
			else if (i==(nPF*2))		//Coupled random Hydrolysis
			{	
				while (1)
				{
					tempcol = floor((double)nPF*RAND01);
					if(PF[tempcol]>3 && GTP[tempcol]>N0)//seed and terminal dimers don't count
					{
						temprow = floor((PF[tempcol]-3)*RAND01)+N0;	//get row within the PF above the seed
						if (lattice[tempcol][temprow]==1)
						{
							lattice[tempcol][temprow]=2;
							GTP[tempcol]--;
							totalGTP--;
							break;
						}
					}
				}
			}
			t = t + tau;	//Increment time
			count = count+1;
			if(t>tout)
			{
				flag = 0;
			}
			if(flag==0)
			{
				if (tout<10)
					tout += delt; //update to 1 second
				else
					tout += delt2;//update to 10 second
				flag = 1;	//Update flag									
				avg2=avg=var=0;
				for (col=0;col<nPF;col++)
				{
					avg += (PF[col]-2)*8;
					avg2 += (PF[col]-2)*8*(PF[col]-2)*8;
					if(col==0)				//Max and min calulation
					{
						max = PF[col]-2;
						min = PF[col]-2;
						//maxcol = col;			//save where the max # of PF occur
						//mincol = col;			//Save the least # of PF
					}
					else
					{
						if (PF[col]-2>max)
						{
							max = PF[col]-2;
							//maxcol = col;
						}
						if (PF[col]-2<min)
					 	{
							min = PF[col]-2;
							//mincol = col;
						}
					}	
				}
				avg2 = avg2/13.0;
				avg = avg/13.0;
				var = avg2-avg*avg;
				tpoint = tpoint + 1; 
				//sprintf(filename,"MTlength_%1.0lf_%1.0lf_%1.0lf_%d",-Glong*mfactor[0],-Glat*mfactor[1],C*1E6,kk+1);//partitio
				fprintf(output,"%1.6lf %1.6lf %1.6lf\n",t,(double)(totaldimer-26)*8.0/13.0/1000,(double)min*8/1000);//check point
				/*if (min>=15)//To calculate decay constant
				{
					sprintf(filename,"Decay_%1.0lf_%1.0lf_%1.0lf_%1.0lf_%1.0lf_%1.0lf_%1.0lf_%d",-Glong*mfactor[0],-Glat*mfactor[1],-ABbond*2*mfactor[2],Hrate*mfactor[3],Ks*mfactor[4],Kplus*mfactor[5],ddG*mfactor[6],tpoint);//partitio
					output = fopen(filename,"a+");
					zerocount=0;	//Don't output zeros
					if((min+N0-1)>=N0)
					{
						for(row=min+N0-1;row>=N0;row--)//go throug one row a time
						{
							ringGTP = 0;
							for(col=0;col<nPF;col++)
							{
								if (lattice[col][row]==1)
									ringGTP = ringGTP+1; 
							}
							fprintf(output,"%d ",ringGTP);//check point
							if(ringGTP==0)
								zerocount=zerocount+1;
							if(zerocount>15)
								break;
						}
					}
					else
						fprintf(output,"%d %d ",0,0);//check point
					fprintf(output,"\n");//check point
					fclose(output);
				}*/
				flag = 1;	//Update flag
				if(tout>MAXtime)
					break;
			}//End flag
		}//End simulation loop
		fclose(output);//Closing output file
	}//End kk loop
free(lattice);
}

int main(int argc,char** argv)
{
	// Simulation setup
	double c=0;
	double CON=8E-6;//Concentration
	double dCON = 1E-6;//Change in Concentration
	double GLONG,GLAT,GSEAM,HRATE,KS,KPLUS,DDG;//All possible parameters	
	//double dGLONG,dGLAT,dGSEAM,dHRATE,dKS,dKPLUS,dDDG;//All possible difference in parameters
	double mfactor[prn] = {100,100,100,1000,1000,1E-4,100};	//For multiplying factors for output
	GLONG=-7.0;
	GLAT=-5.7;
	GSEAM=-5.13;
	HRATE=0.35;
	KS=3.0;
	KPLUS=5E6;
	DDG=3.15;
	double prval[prn]={GLONG,GLAT,GSEAM,HRATE,KS,KPLUS,DDG};	//Set up array for value of parameter
	// Start of simulation
	for (c=(CON-dCON);c<=CON+dCON;c+=(dCON/5))//Loop through Concentration
	{
		printf("Con: %1.2lfuM, Para: %1.2lf %1.2lf %1.2lf %1.2lf %1.2lf %1.2lf %1.2lf\n",c*1E6,prval[0],prval[1],prval[2],prval[3],prval[4],prval[5]*1E-6,prval[6]);
		MAINRUN(prval,mfactor,c);
	}		
	//double tempval[prn];				//Temp variable
	//prdel = [dGLONG,dGLAT,dGSEAM,dHrate,dKs];
	//double prdel[prn]={-0.05,-0.05,-0.05,0.02,0.15,2.5E5,0.05};//Difference for derivative		
	//Run the other  distributions for FIM
	/*for(i=0;i<prn;i++)
	{
		for(k=0;k<prn;k++)			//Set values of tempval 
			tempval[k] = prval[k];
		for(j=0;j<2;j++)
		{
			if (j==0)	//Run prval+prdel
				tempval[i] = prval[i] + prdel[i];
			else if (j==1) 	//Run Rrval-prdel
				tempval[i] = prval[i] - prdel[i];
			printf("1 %lf %lf %lf %lf %lf %lf %lf\n",tempval[0],tempval[1],tempval[2],tempval[3]*1000,tempval[4]*1000,tempval[5]*1E-4,tempval[6]*100);
			MAINRUN(tempval,mfactor);	//Run simulation	
		}
	}*/
return 0;
}
