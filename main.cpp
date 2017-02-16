/*
AUTHORS: 
*SIHAM TABIK
*LUIS FELIPE ROMERO GOMEZ
*ANTONIO MANUEL RODRIGUEZ CERVILLA 
DEPT: ARQUITECTURA DE LOS COMPUTADORES
UNIVERSIDAD DE MALAGA
 
*/

/*
 Things to do:
 * Check step!=10
 * Test cell-sequence kernel
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <math.h>
#include <omp.h>
#include "Defs.h"
#include "Sector.h"
#include "AuxFunc.h"
#include "Earth.h"


using namespace std;


double *heights;	
double obsheight =OBS_H;
int utm_n=UTMN,utm_e=UTME,step=STEP;
int bw;
int dim,dimx,dimy;
extern char filename[100];
extern char base_name[100];
extern char maskname[100];
extern char listname[100];
extern char volcname[100];
extern char horzname[100];
extern char cvrgfname[100];
extern char surname[100], volname[100];
extern char surlname[100], vollname[100];
extern bool towernotgiven;


int tower=1760986;       // Single point viewshed
int *towers;             // Multiple viewshed
int ntowers;
int cnttower;
bool *mask;              // Masked viewshed 
bool *coverage;          // Viewshed for towers
bool *cumcoverage;          // Viewshed for towers
float *cumvolcoverage;	 // Vol viewshed fortowers
float *isolatedMatrix;
unsigned char *horizonsMatrix;
nodoscolores *suelocolores;

float *sur;	         //Total viewshed and volumetric viewshed
float *vol;




//This are default values. To configure them, use configure() carefully.
bool quienmeve=false;
bool nocompute=false;    // For timing the structure manteinance
bool fullstore=false;    // Save all viewshad topologies
bool nostore=false;	 // For timing algorithms
bool volume=false;	 // Volumetric viewshed
bool horizon=false;	 // Horizon kernel
bool inputmask=false;	 // Cell tower kernels
bool sequential=false;	 // Cell tower kernels. Cover in a tower is mask for next one 
bool coveroutput=false;  // Single or multiple point viewshed
bool isolated=false;
bool horizons=false;
bool showtime=false;
bool trackmode=false;   //Code mod to calculate many towers coverage 
                         //for paper at ICCS with random places, equally distrib, etc...
//bool harmonicmean=true;  // Isolation index (instead of max))


int identifier=0;
extern char idname[100];

//==================================================================================
//	Set working mode
//==================================================================================
int mode=0;
bool verbose=false;
bool silent=false;
sectorgeometry G;
int Angle;

void configure(int argc, char *argv[])
{

checkargs(argc,argv);    


G.dimx=dimx;
G.dimy=dimy;
G.bw=bw;
G.step=step;
G.obs_h=obsheight;

    


switch(mode)
{
    case 0://Total viewshed
        break;
    case 1://Total volume visibility
        volume=true;
        break;
    case 2://Single tower viewshed
        fullstore=false;//override command line arg
        coveroutput=true;
		nostore=true; //nostore total data
        break;
    case 3://Several towers viewshed
        fullstore=false;//override command line arg
        coveroutput=true;
        nostore=true; //nostore total data
        trackmode=true;
        sequential=true; // Opcional
        break;
    case 4://Single tower volume visibility
        //inputmask=false; //optional
        fullstore=false;//override command line arg
        coveroutput=true;
	    nostore=true; //nostore total data
        volume=true;
        break;
    case 5://Several towers volume visibility
        //inputmask=true;
        fullstore=false;//override command line arg
        coveroutput=true;
		nostore=true; //nostore total data
        trackmode=true;
        volume=true;
        break;
    case 6://Several towers viewshed sequential mode
        ntowers=5;
        fullstore=false;//override command line arg
        coveroutput=true;
	    nostore=true; //nostore total data
        sequential=true; //Optional. For cell towers algoritm
        trackmode=false;
        break;
    case 7://Isolated areas
        isolated=true;
        break;
    case 8://Horizon computation
        horizons=true;
        break;
    case 9://Timing
        fullstore=false;//override command line arg
        coveroutput=false;
        nocompute=true;
        nostore=true;   //set true for timing
        break;
    case 10://Quien me ve
        quienmeve=true;//override command line arg
        nostore=true;   //set true for timing
        
        break;
}
if(coveroutput&&!trackmode)setTowers();
if(coveroutput&&trackmode)setTracks();
if(quienmeve)setTracks();

if(!setfilenames(fullstore)){printf("Error: An input folder with the bil DEM is required\n");exit(0);}

}



//==================================================================================
//      OPENMP Collect data in natural ordering
//==================================================================================
//openmp shared
bool first=false;
bool first2=false;
bool firstisol=false;


void collectData(int idt, float*threadSur, float* threadVol, bool *threadSurf1pt, float* threadVol1pt, float* threadIsolated, unsigned char * &threadHorizons)
{
#ifdef _OPENMP
#pragma omp critical
	{
#endif
	//if(verbose)printf("Thread %02d - T=%f\n",idt, (timer2-timer0));
if(!nostore){
	if(first==false){
		for(int i=0;i<dim;i++)sur[i]=threadSur[N2T(i)];
		if(volume)
			for(int i=0;i<dim;i++)vol[i]=threadVol[N2T(i)];
		first=true;
		}
	else
		{
		for(int i=0;i<dim;i++)sur[i]+=threadSur[N2T(i)];
		if(volume)
			for(int i=0;i<dim;i++)vol[i]+=threadVol[N2T(i)];
		}
	}


if(coveroutput){
	if(first2==false){
		for(int i=0;i<dim;i++)coverage[i]=threadSurf1pt[i];
                if(volume)for(int i=0;i<dim;i++)cumvolcoverage[i]=threadVol1pt[i];
		first2=true;
		}
	else
		{
		for(int i=0;i<dim;i++)coverage[i]=coverage[i]||threadSurf1pt[i];
                if(volume)for(int i=0;i<dim;i++)cumvolcoverage[i]=max(cumvolcoverage[i],threadVol1pt[i]);
		}
	}

if(isolated){
	if(firstisol==false){memcpy(isolatedMatrix,threadIsolated,dim*sizeof(float));firstisol=true;}
	else {
        //0 Max
        //1 aritmetica
        //2 Harmonica
        //3 Geométrica
            if(isolationindex==0)     
                    for(int i=0;i<dim;i++)isolatedMatrix[i]=max(isolatedMatrix[i],threadIsolated[i]);
            if(isolationindex==1)     
                    for(int i=0;i<dim;i++)isolatedMatrix[i]=isolatedMatrix[i]+threadIsolated[i];
            if(isolationindex==2)     
                    for(int i=0;i<dim;i++)isolatedMatrix[i]=isolatedMatrix[i]+threadIsolated[i];
            if(isolationindex==3)     
                    for(int i=0;i<dim;i++)isolatedMatrix[i]=isolatedMatrix[i]+threadIsolated[i];

        }
}

//Falta por hacer horizons

#ifdef _OPENMP
	}
#endif

}




//==================================================================================
//      Allocate & deallocate shared data
//==================================================================================
bool floatsDEM=false,endianOK=true;

void allocate()
{
    
if(!nocompute)heights=readHeights(filename,floatsDEM,endianOK);
  
if(coveroutput){
	coverage=new bool[dim];
        if(trackmode){
            cumcoverage=new bool[dim];
            for(int j=0;j<dim;j++) cumcoverage[j]=false;
            }
	if(volume)
		{
			cumvolcoverage=new float[dim];
			for(int j=0;j<dim;j++) cumvolcoverage[j]=0;
		}
        //adds tower coverage to input mask (already covered area))

	}

/*
cumvolcoverage=readSurf(volcfname);
for(int j=0;j<dim;j++) cumvolcoverage[j]=step*(cumvolcoverage[j]+heights[2000*(j%2000)+j/2000]);
saveFloats(cumvolcoverage,"prueba.flt");
exit(0);
*/


// Allocate data
if(!nostore){
	sur=new float[dim];
	if(volume) vol=new float[dim];
	}
if(quienmeve)
{
    suelocolores=new nodoscolores[NUMNC*dim];
    nodoscolores init={-1,-1};
    for(int j=0;j<dim*NUMNC;j++) suelocolores[j]=init;
}

if(inputmask)
	{
	if(verbose)printf("Reading coverage as input mask %s\n",maskname);
	mask=new bool[dim];
	inputmask=readCover(mask,maskname);
        //exit(0);
	if(sequential&&inputmask)inputmask=readCover(coverage,maskname); //Add mask as previous coverage
        if(verbose&&inputmask)printf("Reading coverage as input mask %s success\n",maskname);
        if(!inputmask)delete mask; //failed
	}
if(isolated)isolatedMatrix=new float[dim];
if(horizons)horizonsMatrix=new unsigned char[360*dim];
}



void deallocate()
{
if(!nostore){
      	if(volume){
	        delete vol;
		}
		delete sur;
	}
if(coveroutput)
{
	if(volume)
	{
            //delete volcoverage;
	    delete cumvolcoverage;
	}
	delete coverage;
        if(trackmode)delete cumcoverage;
}
if(inputmask)delete mask;
if(!nocompute) 	delete heights;
if(isolated)delete isolatedMatrix;
if(horizons)delete horizonsMatrix;
if(quienmeve) delete suelocolores;
}

//==================================================================================
//      Allocate & Deallocate private memory 
//==================================================================================


void privateAllocate(float * &threadSur, float * &threadVol,bool * &threadSurf1pt, float * &threadVol1pt, float * &threadIsolated, unsigned char * &threadHorizons)
{

	// Per thread storage
	if(!nostore){
			threadSur=new float[dim];
			for (int i = 0; i < dim; i++) threadSur[i] = 0;
		if(volume){
			threadVol=new float[dim];
        	for (int i = 0; i < dim; i++) threadVol[i] = 0;
		}
		}
if(coveroutput){        
        threadSurf1pt=new bool[dim];
        for (int i = 0; i < dim; i++) threadSurf1pt[i] = false;
        threadVol1pt=new float[dim];
        for (int i = 0; i < dim; i++) threadVol1pt[i] = 0;
        }
        
if(isolated)   {
    threadIsolated=new float[dim];     
    for (int i = 0; i < dim; i++) threadIsolated[i] = 0;}
if(horizons)   {
    threadHorizons=new unsigned char[2*dim];
                }
        
    
}


void privateDeallocate(float *&threadSur, float *&threadVol,bool *&threadSurf1pt, float *&threadVol1pt, float * &threadIsolated, unsigned char * &threadHorizons)
{
	if(!nostore)delete threadSur;
	if(!nostore&&volume)delete threadVol;
	if(coveroutput)delete threadSurf1pt;
	if(coveroutput)delete threadVol1pt;
        if(isolated)delete threadIsolated;
        if(horizons)delete threadHorizons;
}


//==================================================================================
//      Save data
//==================================================================================


void saveData(bool last)
{
if(verbose)printf("storing results...\n");
int test_pt;

if(coveroutput){


                if(volume){        
		int cntx=saveCover(coverage,cvrgfname,false);  //append if manytowers
                
                writeheader(utm_n,utm_e,step,cvrgfname,true);
		if(!silent)printf("2D Viewshed %f ha\n",cntx*1.0/(step*step));
                            for(int i=0;i<dim;i++){
                                //cumvolcoverage[i]=max(cumvolcoverage[i],volcoverage[i]);
                                cumvolcoverage[i]*=step;
                                //volcoverage[i]=(cumvolcoverage[i]+heights[i])*step;
                            }
                    if(trackmode){
                            float totalize=0;
                            for(int i=0;i<dim;i++)totalize+=cumvolcoverage[i];
                            if(!silent)printf("3D Viewshed: %f Hm3\n",totalize/10000);
                            saveFloats(cumvolcoverage,volcname,true); //recovers natural ordering
                            writeheader(utm_n,utm_e,step,volcname,false);
                            }
                    else
                    {
//        int imxv=0;
//        float mxv=0;
//        for(int i=0;i<dim;i++)if(cumvolcoverage[i]>mxv){imxv=i;mxv=cumvolcoverage[i];}
//        printf("%d %f\n",imxv,mxv);
                            saveFloats(cumvolcoverage,volcname,true);
                            writeheader(utm_n,utm_e,step,volcname,false);
                            float totalize=0;
                            for(int i=0;i<dim;i++)totalize+=cumvolcoverage[i];
                            //for(int i=0;i<dim;i++)
                            //    volcoverage[i]=cumvolcoverage[i]+heights[i]*step;
                            if(!silent)printf("3D Viewshed: %f Hm3\n",totalize/10000);
                    }
                }
                else{        
                    if(!trackmode){
                        int cntx=saveCover(coverage,cvrgfname,false);  
                        writeheader(utm_n,utm_e,step,cvrgfname,true);
                    }
                    int cnt3=0;
                    if(trackmode){
                            for(int i=0;i<dim;i++){
                                cumcoverage[i]=cumcoverage[i]||coverage[i];
                                cnt3+=(cumcoverage[i])?1:0;
                                }
                            printf("Updating cover %d\n",cnt3);
                        if(last)
                        {
                            printf("Saving cover\n");
                            saveCover(cumcoverage,cvrgfname,false);  //append in tracks, when desired                            
                            writeheader(utm_n,utm_e,step,cvrgfname,true);
                        }
                    }
                }
                


    }


if(!nostore&&!coveroutput){
	//string path="./output/";
       	
       	if(volume){
            saveFloats(vol,volname,false);
       		//fs=fopen(volfname,"wb");
       		//fwrite(vol,dim,4,fs);
       		//fclose(fs);
            if(!strcmp(base_name,"4070000_0310000_010")){
                test_pt=get_point_id(4062820,324430);
                printf("Surface at %d %f\n",test_pt,sur[T2N(test_pt)]);
                printf("Volume  at %d %f\n",test_pt,vol[T2N(test_pt)]);
                test_pt=get_point_id(4059940,313470);
                printf("Surface at %d %f\n",test_pt,sur[T2N(test_pt)]);
                printf("Volume  at %d %f\n",test_pt,vol[T2N(test_pt)]);
                test_pt=get_point_id(4063880,324460);
                printf("Surface at %d %f\n",test_pt,sur[T2N(test_pt)]);
                printf("Volume  at %d %f\n",test_pt,vol[T2N(test_pt)]);
            }
            for(int i=0;i<dim;i++)vol[i]=(vol[i]<0.01)?-2:log10(vol[i]);
            saveFloats(vol,vollname,false);
       		//fs=fopen(vollfname,"wb");
       		//fwrite(vol,dim,4,fs);
       		//fclose(fs);
	        writeheader(utm_n,utm_e,step,volname,false);
	        writeheader(utm_n,utm_e,step,vollname,false);
		}
        saveFloats(sur,surname,false);
       	//fs=fopen(surfname,"wb");
       	//fwrite(sur,dim,4,fs);
       	//fclose(fs);
        
                
        for(int i=0;i<dim;i++)sur[i]=(sur[i]<0.01)?-2:log10(sur[i]);
        saveFloats(sur,surlname,false);
       	//fs=fopen(surlfname,"wb");
       	//fwrite(sur,dim,4,fs);
       	//fclose(fs);
	writeheader(utm_n,utm_e,step,surname,false);
	writeheader(utm_n,utm_e,step,surlname,false);
	}

if(isolated){
        if(isolationindex==0){
            strcat(horzname,"MX");
        }
        if(isolationindex==1){
            for(int i=0;i<dim;i++)isolatedMatrix[i]=isolatedMatrix[i]/360.0; //aritemetica
            strcat(horzname,"MA");
        }
        if(isolationindex==2){
         for(int i=0;i<dim;i++)isolatedMatrix[i]=360.0/isolatedMatrix[i]; //harmonica   
            strcat(horzname,"MH");
        }
        if(isolationindex==3){
            for(int i=0;i<dim;i++)
            {
                float x=isolatedMatrix[i]/360.0;
                if(x<0)printf("%d %e\n",i,x);
                if(x>10)printf("%d %e\n",i,x);
                isolatedMatrix[i]=exp(x);
                if(x<0)printf("%d %e\n",i,isolatedMatrix[i]);
                if(x>10)printf("%d %e\n",i,isolatedMatrix[i]);
            }
            strcat(horzname,"MG");
        }
        saveFloats(isolatedMatrix,horzname,true);
	writeheader(utm_n,utm_e,step,horzname,false);
            }
if(horizons){
        if(verbose)printf("Filename: %s\n",horzname); 
        saveHorizons(horizonsMatrix,horzname);
	writeheader(utm_n,utm_e,step,horzname,false);
            }


if(quienmeve)
{
    if(verbose)printf("Almacenando coloreado del territorio\n");
    saveColores(suelocolores);
}

}



//==================================================================================
//	Computationally intensive section - Parallelization kernel
//==================================================================================

void compute(){
// Generate files with ring sector of all points  


#ifdef _OPENMP
int ncores=omp_get_num_procs()-1;
//ncores=1;
while((NSECTOR%ncores))ncores--;
if(verbose)printf("Using OPENMP with %d cores\n",ncores); 
omp_set_num_threads(ncores);
#pragma omp parallel private(Angle)
#endif
	{
//==================================================================================
//	Parallel job
//==================================================================================
	int nthreads=1;
	int idt=0;
#ifdef _OPENMP
	nthreads =omp_get_num_threads();	
	idt=omp_get_thread_num();
#endif
        //openmp private
        float *cumsur;	//Per thread job
        float *cumvol;
        bool *threadSur1pt=NULL; 
        float *threadIsolated=NULL; 
        unsigned char *threadHorizons=NULL; 
        float *threadVol1pt=NULL;
        double timer0,timer1,timer2,timer3,timer4;
        //if(verbose)printf("Thread %d of %d\n",idt,nthreads);
	Sector *A;
        //G.obs_h=obsheight;
        
        /*
	if(!coveroutput) A=new Sector(G,nostore, fullstore, volume);
	else             A=new Sector(G,tower,volume);
	*/
//==================================================================================
        A=new Sector(G,nostore, fullstore, volume,coveroutput);

	A->verbose=verbose;
        A->idt=idt;
        A->inputmask=inputmask;
	A->mask=mask;
	A->nocompute=nocompute;
	A->isolated=isolated;
	A->horizons=horizons;
        A->quienmeve=quienmeve;
        if(coveroutput||quienmeve){
            A->ntowers=sequential?1:ntowers;
            A->towersloc=sequential?(new int[1]):towers;
            A->towerloc=tower;
            if(sequential)A->towersloc[0]=tower;
            }
	A->setType();
	A->init_storage(); 
        if(!nocompute)A->setHeights(heights);
        privateAllocate(cumsur,cumvol,threadSur1pt,threadVol1pt,threadIsolated,threadHorizons);
	if(coveroutput)A->coverage=threadSur1pt; 
	if(coveroutput)A->volcoverage=threadVol1pt; 
	if(isolated)A->isolatedMatrix=threadIsolated;
	if(horizons)A->horizonsMatrix=threadHorizons;

//==================================================================================
//      Repeat loop for sequential viewshed 
//==================================================================================
        bool repeat=false;
        int iterate=0;
        do{
#pragma omp barrier
            if(sequential){
                 if(idt==0)printf("New tower, towers[%d]=%d\n",cnttower-1,tower);
                 resetCover(coverage);
                 resetCover(A->coverage);
                 if(verbose&&(idt==0))printf("Iteration: %d\n",iterate++);
            }
            repeat=false;
            timer0=dtime();
            int s=0; int ss=0, nsect=NSECTOR;
#ifdef _OPENMP
            nsect=NSECTOR/nthreads;
            ss=omp_get_thread_num()*nsect;
#endif
        
            //==================================================================================
            //      Sector loop
            //==================================================================================
            //if(verbose)printf("Thread %d starts sector loop\n",idt,nthreads);
            for(s=ss;s<ss+nsect;s++){ //sector loop
        	FILE *fs;
                Angle=s;
	       	char nfile[100];
                sprintf(nfile,"%s/%d/nn_list%d/nn_list_%03d.bin",CACHED,NSECTOR,bw,s);
		//nfile= CACHED+"/nn_list" + converInt(BW) + "/nn_list_" + converInt(s) + ".bin";
                //if(verbose)printf("%d in sector loop with %s\n",idt,nfile);
		bool precalc=false;
                if(precalc=!isCached(bw,s)){
			fs=fopen(nfile,"wb");
			A->listfileW = fs;         
		}
		else{
			fs=fopen(nfile,"rb");
			A->listfileR = fs;         
		}


		A->precomputed=!precalc;
		A->change(s);
		timer1=dtime();
		A->loop();
		timer2=dtime();
		if(showtime)printf("Thread %02d Sector %d - T=%f\n",idt,s, (timer2-timer1));
		fclose(fs);                
	        
		if(!nostore) {
		        char fnin[100];
                        sprintf(fnin,"%s/full/rsector_%03d.bin",OUTPUT,s);
			if (fullstore) A->recordsectorRS(fnin);
        		for (int i = 0; i < dim; i++) cumsur[i]+=A->surfaceF[i]+A->surfaceB[i]; 
        		if(volume) 
				for (int i = 0; i < dim; i++) cumvol[i] +=A->volumeF[i]+A->volumeB[i];
			}
	        //if(verbose)printf("Thread %02d Sector %d - T=%f\n",idt,s, (timer2-timer1));
                if(horizons)
                        for(int i=0;i<dim;i++)
                                {
                                horizonsMatrix[i*360+((Angle+90)%360)]=threadHorizons[T2N(i)*2+0];
                                horizonsMatrix[i*360+180+((Angle+90)%360)]=threadHorizons[T2N(i)*2+1];
                            }
		} //end sector loop

            //==================================================================================
            // Data Collection from threads
            //==================================================================================
            //printf("%02d\n",idt);
            collectData(idt,cumsur,cumvol,threadSur1pt,threadVol1pt,threadIsolated,threadHorizons);        

            if(sequential){
                if(trackmode&&(ntowers-cnttower)){
                    repeat=true;
                    A->towersloc[0]=tower=towers[cnttower];
                    }
                if(!trackmode&&(iterate<1)){ //No implementado. Debe alternar mode tower con mode no tower
                    repeat=true;
                    if(idt==0)setTowers();
                    }
                }
#pragma omp barrier
		timer3=dtime();
	        if(showtime)printf("Thread %02d Sector %d - Tfinal=%f\n",idt,s, (timer3-timer0));
            if(idt==0){
                cnttower++;
                saveData(!repeat);        
            }
        }while(repeat);
        
        
        //==================================================================================
        // End Parallel section
        //==================================================================================
        privateDeallocate(cumsur,cumvol,threadSur1pt,threadVol1pt,threadIsolated,threadHorizons);
        delete A; //Cada thread maneja su propia estructura
        } //end parallel
}


//==================================================================================
//	MAIN
//==================================================================================



int main(int argc,char*argv[])
{	
//------------------------------------------------------------
//	Kernel selection (usually, viewshed)
//------------------------------------------------------------

configure(argc,argv);    
if(verbose)printf("Input filename: %s\n",filename);
allocate();
if(verbose)printf("Start computing\n");
compute();
if(verbose)printf("End computing\n");
deallocate();
if(verbose)printf("End write results\n");

if(!nocompute&&!coveroutput)savegoogleearth(volume);

#ifdef WIN32
//getchar();
#endif
exit(0);
}
