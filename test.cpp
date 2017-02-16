/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <math.h>
#include <omp.h>
#include "Defs.h"
#include "AuxFunc.h"
#include "Sector.h"
#include "Earth.h"
#include "lodepng.h"

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

double *heights;	
double obsheight =OBS_H;
int utm_n=UTMN,utm_e=UTME,step=STEP;
int bw;
int dim,dimx,dimy;
extern char idname[100];
int mode=0;
bool verbose=false;
bool silent=false;
sectorgeometry G;
int Angle;

int tower=1760986;       // Single point viewshed
int *towers;             // Multiple viewshed
int ntowers;
int cnttower;
bool floatsDEM=false,endianOK=true;
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

int identifier=0;



void configure(int argc, char *argv[])
{
checkargs(argc,argv);    


G.dimx=dimx;
G.dimy=dimy;
G.bw=bw;
G.step=step;
G.obs_h=obsheight/G.step;

switch(mode)
{
    case 0://Total viewshed
        break;
}
}


double **profile;
void get_prof(int n,double a,int length)
{
    a=PI*a/180;
    int col;
for(int row=0;row<length;row++){
    float x= tan(a)*row;
    int ix=floor(x);
    float dx=x-ix;
    float cdx=1-dx;
    int col=n+ix;
    profile[n][row]=cdx*heights[dimx*row+col]+dx*heights[dimx*row+col+1];
    //profile[0][row]=heights[dimx*row+start];
    //if(row==length-1)printf("%e %e %e \n",a,cdx,dx);
    }
    
//for(int i=0;i<length;i+=100)printf("%d %d\n",(int)(10*profile[0][i]),(int)(10*profile[1][i]));
}

void allocate()
{
    heights=readHeights(filename);
    profile=(double **)malloc(dimx*sizeof(double *));
    //printf(" %d %d %d %d   \n",(int)heights[0],(int)heights[dimx-1],(int)heights[dim-dimx],(int)heights[dim-1]);
    for(int i=0;i<dimx;i++)profile[i]=(double *)malloc(sqrt(dim)*sizeof(double));
}
void deallocate()
{
    free(heights);
    for(int i=0;i<dimx;i++)free(profile[i]);
    free(profile);
}
float **resultado;
void alojador()
{
    resultado=(float **)malloc(5000*sizeof(float *));
    for(int i=0;i<2500;i++){
        resultado[i]=(float *)malloc(2500*sizeof(float));
        resultado[2500+i]=(float *)malloc(2500*sizeof(float));
        for(int j=0;j<2500;j++){
            resultado[2500+i][j]=heights[dimx*i+j];
            resultado[i][j]=0;
        }
        }
}
void desalojador()
{
    for(int i=0;i<2500;i++){
        free(resultado[i]);
        free(resultado[2500+i]);
    }
    free(resultado);
}










bool visible=false;
double aw=3;
double range;

//determina si el punto del final es visible o no, y devuelve el conjunto de puntos intermedios visibles
int get_vs(int prof, int length, double a)
{
    	float d = 0;
	float h = profile[prof][0] + G.obs_h;
        float delta_d=0;
        float delta_h=0;
        float max_angle=-2000;
        int cnt=0;
        double icos=1.0/cos(a);
        //icos=1;
        float open_delta_d=0;
        float open_delta_h=0;
        visible=false;
        
        for(int i=1;i<length;i++){
            //delta_d=i*icos-d;
            delta_d=i*icos;
            delta_h=profile[prof][i]-h;
            float angle = delta_h/delta_d;
            bool above = (angle > max_angle);
            bool opening = above && (!visible);
            bool closing = (!above) && (visible);
            //if(i<100)printf("%d %e %e \n",above,angle, 180*atan(angle)/PI);
            //if(opening)printf("%d",i);
            //if(closing)printf("-%d ",i);
            visible = above;
            max_angle = max(angle, max_angle);
            if(visible)cnt++;
        }
        //printf("Viewshed: %d",cnt);
        return cnt;
}

//determina si el punto del final es visible o no, y devuelve el vs por rs
double get_vs2(int prof, int length, double a)
{
    	float d = 0;
	float h = profile[prof][0] + G.obs_h;
        float delta_d=0;
        float delta_h=0;
        float max_angle=-2000;
        double cnt=0;
        double icos=10.0/cos(a);
        double icos2=icos/2;
        float open_delta_d=0;
        visible=false;
        float cte2=2*aw*3.1415926/360.0;
        range=0;
        for(int i=1;i<length;i++){
            //delta_d=i*icos-d;
            delta_d=i*icos;
            delta_h=profile[prof][i]-h;
            float angle = delta_h/delta_d;
            bool above = (angle > max_angle);
            bool opening = above && (!visible);
            bool closing = (!above) && (visible);
            visible = above;
            max_angle = max(angle, max_angle);
            if (opening)open_delta_d = delta_d-icos2;
            if (closing){
                range=(delta_d-icos2);
                cnt += ( range*range - open_delta_d * open_delta_d);
                }
        }
        if(visible)range=(delta_d-icos2);
        //printf("Viewshed: %d",cnt);
        return cnt*cte2;
}


int central(int i, double angle)
{
        get_prof(i,angle,2500);
        double x=get_vs2(i,2500,angle*PI/180);
        return (int)(x);
}

int total(int point, double angle)
{
    int cnt=0;
    int cnt2=0;
    for(int i=1;i<2500;i++) //todas las filas
    for(int j=0;j<2500;j++) //todas las columnas
    {
        double a=atan2(j-point,i)*180/PI;
        if( (a>(angle-aw))&&(a < (angle+aw))){
            cnt2++;
            int farthest=i;
            get_prof(point,a,farthest);
            int x=get_vs(point,farthest,a*PI/180);
            if(visible)cnt++;
    }
    }
    return cnt;
}


void calculador()
{
    
    double angle=22;
    //angle*=(PI/180.0);
    int cs[2000];
    int ts[2000];
    int rs[2000];
    for(int i=0;i<2000;i+=1){
        cs[i]=central(i,angle);
        ts[i]=total(i,angle);
        rs[i]=(int)range;
        printf("%d %d Range: %d\n",cs[i],100*ts[i],(int)range);
    }
    FILE *f=fopen("stats22q.bin","wt");
    for(int i=0;i<2000;i+=1)fprintf(f,"%d,%d,%d\n",cs[i],100*ts[i],rs[i]);
    fclose(f);
    
    
}

void borrador()
{
    for(int i=0;i<2500;i++)
        for(int j=0;j<2500;j++)
        {
            resultado[i][j]=0;
            resultado[i+2500][j]=0;
        }
    
}


bool compactar=false;
bool borrar=true;


void girador1(int a)
{
    if(borrar)borrador();
    if(a<=45){
    for(int i=0;i<2500;i++){
        for(int j=0;j<2500;j++){
            int k=tan(a*PI/180.0)*j;
            int l=tan(a*PI/180.0)*i;
            if(i-k>=0){  //el punto cae en la parte inferior y no hay que compactar a izqda la matriz
//            resultado[2500+i][j]=0;
            resultado[2500+i-k][j]=heights[dimx*i+j];
            }
            else
            {
            resultado[2500+i-k][j]=heights[dimx*i+j];
            }
        }
        }
    if(compactar)
    for(int i=0;i<2500;i++){
        int ji=0;
        for(ji=0;ji<2500;ji++)if(resultado[i][ji])break;
        if(ji<2500)for(int j=0;j<2500;j++)resultado[i][j]=((ji+j)<2500)?resultado[i][ji+j]:0;
            
    }
    }
    else
    {
    for(int i=0;i<2500;i++){
        for(int j=0;j<2500;j++){
            int b=90-a;
            int k=tan(b*PI/180.0)*j;
            int l=tan(b*PI/180.0)*i;
            int jt=2500-j;
            if(j-l>=0){
              resultado[jt+l][i]=heights[dimx*i+j];
            }
            else
            {
            //resultado[2500+i][j]=0;
            //resultado[2500+i-k][l]=heights[dimx*i+j];
              resultado[jt+l][i]=heights[dimx*i+j];
            }
        }
        }
    
    if(compactar)
       for(int i=0;i<2500;i++){

        int ji=0;
        //buscar la primera fila vacía
        for(ji=0;ji<2500;ji++)if(resultado[2500+i][ji])break;
        if(ji<2500)
            for(int j=0;j<2500;j++)
                resultado[2500+i][j]=((ji+j)<2500)?resultado[2500+i][ji+j]:0;
       }           
 
    
    }
    
}

double timer0,timer1,timer2,timer3,timer4;

void girador3(int a)
{
    int dimx=(2500/256)*256;
    int dimy=(2500/256)*256;
    int nsbx=2; //Num of superblocks(core)
    int nsby=2;
    int sblockx=dimx/nsbx;//tamaño superbloque x
    int sblocky=dimy/nsby;
    int nbx=sblockx/128;
    int nby=sblocky/128;
    int blockx=sblockx/nbx;
    int blocky=sblocky/nby;
    int bandsize=blockx+blocky; //crea las franjas de blockx+blocky de alto
    int sbandsize=sblockx+sblocky; //crea las franjas de blockx+blocky de alto
    int sfilabase;
    printf("%dx%d superbloques de %dx%d, %dx%d bloques de %dx%d\n",nsby,nsbx,sblocky,sblockx,nby,nbx,blocky,blockx);
    if(borrar)borrador();
    float ta=tan(a*PI/180.0);
    for(int isn=0;isn<nsby;isn++){
        for(int jsn=0;jsn<nsbx;jsn++){
            int srow=isn*sblocky;
            int scol=jsn*sblockx;
            if(a<=45)
            {
               sfilabase=isn*sbandsize; 
            }
            else
            {
               sfilabase=jsn*sbandsize; 
            }
            for(int in=0;in<nby;in++){
                for(int jn=0;jn<nbx;jn++)
                {
                    if(a<=45)
                        for(int ii=0;ii<blocky;ii++)
                            for(int jj=0;jj<blockx;jj++){
                            int i=ii+in*blocky+srow;    
                            int j=jj+jn*blockx+scol;    
                            int k=ta*jj;
                            int filabase=sfilabase+in*bandsize+blockx; //de 0 a 45 ocupa la mitad baja de la franja y avanza blockx a la alta
                            int fila=filabase+ii;
                            resultado[fila-k][j]=heights[2500*i+j];
                            }
                    if(a>45)
                        for(int ii=0;ii<blocky;ii++)
                            for(int jj=0;jj<blockx;jj++){
                            int i=ii+in*blocky+srow;    
                            int j=jj+jn*blockx+scol;    
                            int l=ii/ta;

                            int jt=dimx-j-1; //de derecha a izqda
                           
                            int filabase=(jt/sblockx)*sbandsize+((jt%sblockx)/blockx)*bandsize; //de 46 a 90 ocupa la mitad alta de la franja y retrae desde la baja
                            int fila=filabase+jt%blocky;
                            if(fila+l>=5000)printf("%d %d %d %d\n",jt, sbandsize,bandsize,filabase);
                            resultado[fila+l][i]=heights[2500*i+j];
                            }
                }}}}
    
}



void girador2(int a)
{
    int blockx=10;
    int blocky=100;
    int nybandas=2500/blocky;
    int nxbandas=2500/blockx;
    int banda=blockx+blocky; //crea las franjas de blockx+blocky de alto
    if(borrar)borrador();
    if(a<=45){
    for(int i=0;i<2500;i++){
        for(int j=0;j<2500;j++){
            int k=tan(a*PI/180.0)*(j%blockx);
            //resultado[2500+i-k][j]=heights[dimx*i+j];
            int filabase=(i/blocky)*banda+blockx; //de 0 a 45 ocupa la mitad baja de la franja y avanza blockx a la alta
            int fila=filabase+i%blocky;
            resultado[fila-k][j]=heights[dimx*i+j];
        }
        }
    }
    else
    {
    for(int i=0;i<2500;i++){
        for(int j=0;j<2500;j++){
            int b=90-a;
            int l=tan(b*PI/180.0)*(i%blockx);
            int jt=2500-j-1; //de derecha a izqda
            int filabase=(jt/blocky)*banda; //de 46 a 90 ocupa la mitad alta de la franja y retrae desde la baja
            int fila=filabase+jt%blocky;
            resultado[fila+l][i]=heights[dimx*i+j];
        }
        }
    
 
    
    }
    
    if(compactar)
        
        for(int ii=0;ii<nybandas;ii++)
            for(int jj=0;jj<nxbandas;jj++)
            {
            int first=0;
        //buscar la primera fila vacía
            for(int i=0;i<banda;i++){
                for(first=0;first<blockx;first++)if(resultado[ii*banda+i][jj*blockx+first])break;
                if(first<blockx)
                    for(int j=0;j<blockx;j++)
                        resultado[ii*banda+i][jj*blockx+j]=((first+j)<blockx)?resultado[ii*banda+i][jj*blockx+first+j]:0;
            }
                
            }

        
        for(int i=0;i<2500;i++){

       }           

    
}


extern color *my_palette;
std::vector<unsigned char>  prep_image(float **data, float maxval)
{
my_palette=build_palette();
std::vector<unsigned char> image;
image.resize(dimx * 2* dimy * 4);
if(verbose)printf("Preparando imagen \n");
for(int i = 0; i < 2*dimy; i++)
for(int j = 0; j < dimx; j++)
  {
    int k=dimx*i+j;
    int v;
    v=1024.0*(data[i][j]/maxval);
    if(v<0)v=0;
    if(v>1023)v=1023;
    //int x=(i%dimx)*dimy+i/dimx;
    image[4 * k + 0] = my_palette[v].R;
    image[4 * k + 1] = my_palette[v].G;
    image[4 * k + 2] = my_palette[v].B;
    image[4 * k + 3] = 255;
  }
return image;
}

bool modotiming=false;
bool mododibujo=true;




int main(int argc,char*argv[])
{	
int folder=7;
configure(argc,argv);   
if(verbose)printf("Input filename: %s\n",filename);
allocate();
alojador();
timer0=dtime();
if(true)
for(int i=0;i<=90;i++){
    if(verbose)printf("Girando en modo %d %02d \n",mode,i);
    if(mode==1)girador1(i);
    if(mode==2)girador2(i);
    if(mode==3)girador3(i);
        
    if(!modotiming){
        std::vector<unsigned char> image=prep_image(resultado,200);
        if(verbose)printf("Creando imagen \n");
        char fn[100];
        sprintf(fn,"carpeta%d/im%02d.png",folder,i);
        lodepng::encode(fn, image, dimx, 2*dimy);
        }
}
for(int i=91;i<=90;i++){
    if(!modotiming){
    if(verbose)printf("Girando en modo %d %02d \n",mode,i);
    if(mode==2)
    girador2(i);
    else
    girador1(i);
        std::vector<unsigned char> image=prep_image(resultado,200);
        if(verbose)printf("Creando imagen \n");
        char fn[100];
        sprintf(fn,"carpeta%d/im%02d.png",folder,i);
        lodepng::encode(fn, image, dimx, 2*dimy);
        }
}
timer1=dtime();
printf("Tiempo: %f\n",timer1-timer0);
if(mododibujo)desalojador();
if(mododibujo)deallocate();
if(mododibujo)exit(0);
calculador();
desalojador();
deallocate();

#ifdef WIN32
//getchar();
#endif
exit(0);
}