/*
AUTHORS: 
*SIHAM TABIK
*LUIS FELIPE ROMERO GOMEZ
*ANTONIO MANUEL RODRIGUEZ CERVILLA 
DEPT: ARQUITECTURA DE LOS COMPUTADORES
UNIVERSIDAD DE MALAGA
 
*/

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <math.h>
#include <omp.h>


#include "Defs.h"
#include "Sector.h"


#define _USE_MATH_DEFINES
using namespace std;


int utm_n=UTMN,utm_e=UTME,step=10;
int bw;
int xdim=2000,ydim=2000;
double *heights;	

string pathfsort= CACHED;		//Directory order file 

void precomputing()
{
FILE *fs;
char filename[100];
string nfile;     
bw=((xdim+ydim)/2)+1;
sectorgeometry G;
G.dimx=xdim;
G.dimy=ydim;
G.step=step;
G.obs_h=0;
G.bw=bw;
Sector *A=new Sector( G);
A->precomputed=false;
A->nocompute=true;
A->nostore=true;
A->fullstore=false;
A->init_storage();
string command="mkdir "+pathfsort;
system(command.c_str());
for (int s = 0; s < 180; s++)
	{
	sprintf(filename,"%s/%s%03d%s",pathfsort.c_str(),"nn_list_",s,".bin");
	printf("Precomputing %s \n",filename);
	fs=fopen(filename,"wb");
	A->listfileW = fs;
	A->change(s);
	A->loop();
	fclose(fs);
	}
delete A;
}



int main()
{	
string nmde,nmap;
precomputing();
}
