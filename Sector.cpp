/*
AUTHORS: 
*SIHAM TABIK
*LUIS FELIPE ROMERO GOMEZ
*ANTONIO MANUEL RODRIGUEZ CERVILLA 
DEPT: ARQUITECTURA DE LOS COMPUTADORES
UNIVERSIDAD DE MALAGA
 
*/

#include <stdio.h>
#include <cmath>
#include <omp.h>
#include "Sector.h"

using namespace std;
	 
extern nodoscolores* suelocolores;
extern bool alineados(int a, int b, int c);

void Sector::default_parms(sectorgeometry G)
{
surscale=PI_F/(360*G.step*G.step);  //results in hectometers ^2 ^3
volscale=surscale/(3*G.step);
dimx=G.dimx;
dimy=G.dimy;
bw=G.bw;
hbw=(G.bw-1)/2;
N=G.dimx*G.dimy;
origin = SHIFT; 
obs_h=G.obs_h/G.step;
step=G.step;
//printf("step=%e\n",step);
}

Sector::Sector(sectorgeometry G)
{ 
	default_parms(G);
        precomputed = true;
        nocompute=true;
        nostore=true;
        fullstore=true;
        volume=false;
}


Sector::Sector(sectorgeometry G,int i, bool lvolume)
{
	default_parms(G);
        precomputed = true;
        nocompute=true;
        nostore=true;
        fullstore=true;
        volume=lvolume;
        coverstore=true;
        towerloc=i;
        //printf("%d %d\n",i,towerloc);
}

// Total surface, Total Volume and (future) Horizon kernels
Sector::Sector(sectorgeometry G,bool nostore_i, bool fullstore_i, bool volume_i, bool cover_i)
{ 
	default_parms(G);
        vectorize=false;
        //vectorize=false;
        precomputed = true;
        nocompute=cover_i;
        coverstore=cover_i;
        inputmask=false;
        nostore=nostore_i;
        fullstore=fullstore_i;
        volume=volume_i;
 }
	// ============================================================================================
	// ============================================================================================
  // Initial tasks:
	// ============================================================================================
	// ============================================================================================

void Sector::init_storage()
{
nodes= new node[N];
MyLinkedList = LinkedList(bw); 

isin = new double[N]; 
icos = new double[N]; 
icot = new double[N]; 
itan = new double[N]; 
orderiAT = new int[N];
tmp1 = new int[N]; 
tmp2 = new int[N];
if(!nostore){
		surfaceF=new float[N];
		surfaceB=new float[N];
		if(volume)volumeF=new float[N];
		if(volume)volumeB=new float[N];
        	if(fullstore){
        		rsectorF = new int*[N];
        		rsectorB = new int*[N];
        		size_dsF=new unsigned short[N];
        		size_dsB=new unsigned short[N];
        		}
		}
if(coverstore){
		//coverage_out=new bool[dim];
    rsectorF = new int*[1];
    rsectorB = new int*[1];
    size_dsF=new unsigned short[1];
    size_dsB=new unsigned short[1];
	  }
if(vectorize)init_vec(bw);

#ifdef USE_CLOSEST_POINTS
                    cls_flag = new bool[2*dimx*dimy]; //forw and backw
                    cls_angle = new float[2*dimx*dimy];
#endif

}

void Sector::setHeights(double *lheights)
{
heights=lheights;
for (int i = 0; i < N; i++) {
	nodes[i].idx = i;
       	nodes[i].h = (float)lheights[i]; //decameters
     	}
}
 
	// ============================================================================================
	// ============================================================================================
  // Finishing tasks:
	// ============================================================================================
	// ============================================================================================

Sector::~Sector()
{
delete isin;
delete icos;
delete icot;
delete itan;

delete orderiAT;
delete tmp1;
delete tmp2;
#ifdef USE_CLOSEST_POINTS
                    delete cls_flag; //forw and backw
                    delete cls_angle;
#endif
if(!nostore){
        delete surfaceB;
        delete surfaceF;
        if(volume)delete volumeB;
        if(volume)delete volumeF;
        if(fullstore||coverstore){
         	delete rsectorF;
         	delete rsectorB;
         	delete size_dsF;
         	delete size_dsB;
	 	}
//if(coverstore)delete coverage_out;
	//if(coverinput)delete coverage_in;
	}
MyLinkedList.~LinkedList();
if(vectorize)dispose_vec();
}



	// ============================================================================================
	// ============================================================================================
  // per sector specific functions:
	// ============================================================================================
	// ============================================================================================

//Task to do when a new direction is taken


void Sector::change(int n)
{
sector = n;
sector_step=(180.0/NSECTOR);
sector_angle=n*sector_step;


sect_angle_g=(origin+sector_angle+(sector_step/2));
sect_langle_g=(origin+sector_angle);
sect_hangle_g=(origin+sector_angle+sector_step);
quad = (sector_angle >= 90) ? 1 : 0;
if (quad == 1) sector_angle -= 90; 
pre_computed_sin();
sort();
if (quad == 1) sector_angle += 90; 
setDistances(sector_angle);
MyLinkedList.Clear();
MyLinkedList.FirstNode(nodes[orderiAT[0]]);
}

//precomputed values for trigonometry
void Sector::pre_computed_sin()
{
double ac, al, ar;
double s, c, ct, tn;

ac = (origin + sector_angle + (sector_step/2)) * torad;
al = (origin + sector_angle) * torad;
ar = (origin + sector_angle + sector_step) * torad;
sect_angle = ac;
s = sin(sect_angle); 
c = cos(sect_angle); 
tn = tan(sect_angle); 
ct = 1 / tn;

for (int i = 0; i < N; i++)
{
    icos[i] = i * c;
    isin[i] = i * s;
    icot[i] = i * ct;
    itan[i] = i * tn;
 }
}

        

//Distance to segment over point 0 (NW)
void Sector::setDistances(int sector)
{
double val;
for (int j = 0; j < dimx; j++)
	for (int i = 0; i < dimy; i++)
		{
		val= icos[j] + isin[i];
            	if (quad == 1) val = icos[i] - isin[j];
		nodes[j * dimy + i].d = val;
    }
}   



void Sector::presort()
{ 
double ct,tn;
tmp1[0] = 0;
tmp2[0] = 0;
	
tn = tan(sect_angle); 
ct = 1 / tn;
for (int j = 1; j < dimx; j++)
	tmp1[j] = tmp1[j - 1] + (int)min(dimy, (int)floor(ct * j));
for (int i = 1; i < dimy; i++)
        tmp2[i] = tmp2[i - 1] + (int)min(dimx, (int)floor(tn * i));
}
      

void Sector::sort()
{
presort();
double lx = dimx - 1;
double ly = dimy - 1;
double x,y;
int ind,xx,yy,p,np;
for (int j = 1; j <= dimx; j++)
{
x = (j - 1);
    for (int i = 1; i <= dimy; i++)
        {
        y = (i - 1);
        ind = i * j;
        ind += ((ly - y) < (icot[j - 1])) ? ((dimy - i) * j - tmp2[dimy - i] - (dimy - i)) : tmp1[j - 1];
        ind += ((lx - x) < (itan[i - 1])) ? ((dimx - j) * i - tmp1[dimx - j] - (dimx - j)) : tmp2[i - 1];
        xx = j - 1;
        yy = i - 1;
        p = xx * dimy + yy;
        np = (dimx - 1 - yy) * dimy + xx;
        if (quad == 0)
        {
            nodes[p].oa = ind - 1;
            orderiAT[ind - 1] = np;
        }else
        {
            nodes[np].oa = ind - 1;
            orderiAT[N - ind] = p;
        }
    }
    
}
       
}

	// ============================================================================================
	// ============================================================================================
	// Data structure manteinance
	// ============================================================================================
	// ============================================================================================


void Sector::post_loop(int i, bool opening, bool closing)
{
	atright =false;
	atright_trans =false;
	NEW_position = -1;
	NEW_position_trans = -1;
	node tmp = newnode_trans;
	
	if(!closing)
	{
		newnode = nodes[orderiAT[opening ? 2*i+1:hbw+i+1]];
		atright = newnode.oa > MyLinkedList.LL[POV].Value.oa;
	}
	if (opening)
	{
		newnode_trans = nodes[orderiAT [2*i+2]];
		atright_trans = newnode_trans.oa > MyLinkedList.LL[POV].Value.oa;
	}	
	
  if(precomputed){
            if (!closing)
            {
              fread(&NEW_position,4,1,listfileR);
              insert_node(newnode, NEW_position, !opening);
            }
            if(opening)
            {
             fread(&NEW_position_trans,4,1,listfileR);
             insert_node(newnode_trans,NEW_position_trans,false);
            }
		}
	else
            {
            if (!opening&&!closing)
                          {
			                newnode = nodes[orderiAT[opening ? 2 * i + 1 : hbw + i + 1]];
                    	NEW_position = -1;
                    	atright = newnode.oa > MyLinkedList.LL[POV].Value.oa;
                    	calcule_pos_newnode(true);
              }
            if (opening){
                        newnode = nodes[orderiAT[2 * i + 1]];
                        NEW_position = -1;
                        atright = newnode.oa > MyLinkedList.LL[POV].Value.oa;
                        calcule_pos_newnode(false);
                        newnode = nodes[orderiAT[2 * i + 2]];
                        atright = newnode.oa > MyLinkedList.LL[POV].Value.oa;
                        calcule_pos_newnode(false);
                        }
		}
	if (closing) MyLinkedList.Remove_two();
 }

void Sector::calcule_pos_newnode(bool remove)
{
            int sweep = 0;
            if (newnode.oa < MyLinkedList.LL[MyLinkedList.First].Value.oa)
            {
                NEW_position = -2; // Before first
                MyLinkedList.AddFirst(newnode, remove);//MyLL
                fwrite(&NEW_position,4,1,listfileW);
                sweep = MyLinkedList.First;
            }
            else if (newnode.oa > MyLinkedList.LL[MyLinkedList.Last].Value.oa)
            {
                NEW_position = -1;
                MyLinkedList.AddLast(newnode, remove);
                fwrite(&NEW_position,4,1,listfileW);
            }

            else
            {
                sweep = MyLinkedList.LL[MyLinkedList.First].next;
                bool go_on = true;
                do
                {
                    if (newnode.oa < MyLinkedList.LL[sweep].Value.oa)
                    {
                        NEW_position = MyLinkedList.LL[sweep].prev;
                        MyLinkedList.Add(newnode, NEW_position, remove);//MyLL
                        fwrite(&NEW_position,4,1,listfileW);
                        sweep = MyLinkedList.LL[NEW_position].next;
                        go_on = false;
                    }
                    else
                        sweep = MyLinkedList.LL[sweep].next;
                } while (go_on);
            }
}


void Sector::insert_node(node newnode, int position, bool remove)
{
	if(position > -1)
	{
		 MyLinkedList.Add(newnode,position,remove);
	}else
	{
		if(position==-1)
		{
			MyLinkedList.AddLast(newnode,remove);
		}else
		{
			MyLinkedList.AddFirst(newnode,remove);
		}
	}
	 
}



	// ============================================================================================
	// ============================================================================================
	// Main loop accross point
	// ============================================================================================
	// ============================================================================================

void Sector::setType()
{
type=1;
if(vectorize)type=3;
if(volume)type=4;
if(!volume&&inputmask)type=5;
if(volume&&inputmask)type=6;
if(coverstore)type=2;
if(nocompute&&!coverstore)type=0;
if(isolated)type=7;
if(horizons)type=8;
if(quienmeve)type=10;
//if((idt==0)&&verbose)printf("Sector kernel Type %d\n",type);
}


void Sector::loop()
{

#ifdef USE_CLOSEST_POINTS
                    closest_parameters(sect_angle_g);
                    //cls_flag = new bool[2*dimx*dimy]; //forw and backw
                    //cls_angle = new float[2*dimx*dimy];
                    cls_compute(sect_angle_g,heights,obs_h, dimx,dimy);
#endif



bool opening,closing;
//int diffCount=0;
//MyLinkedList.saveLL(MyLinkedList.savedLLnew);
for(int i=0; i<N-1; i++) {
	POV=i%bw;  
	kernel_select(i, type);
	opening= (i<hbw);
	closing= (i>= (N-hbw-1));
	//MyLinkedList.saveLL(MyLinkedList.savedLLold);
	post_loop(i,opening,closing);
	//MyLinkedList.saveLL(MyLinkedList.savedLLnew);
	//diffCount+=MyLinkedList.compareLL(MyLinkedList.savedLLold,MyLinkedList.savedLLnew);
	}
}




	// ============================================================================================
	// ============================================================================================
  // End of direction specific functions: Store results
	// ============================================================================================
	// ============================================================================================


//Fulldata to store ring sector topology 
//
void Sector::storers(int i, bool volume, bool fulldata)
{
surfaceF[i]=max((float)0.0,sur_dataF * surscale) ;
surfaceB[i]=max((float)0.0,sur_dataB * surscale) ;

if(volume)
{
volumeF[i]=max((float)0.0,vol_dataF * volscale) ;
volumeB[i]=max((float)0.0,vol_dataB * volscale) ;
//printf("%f\n",volumeF[i]);
}
if(fulldata){
	rsectorF[i] = new int[nrsF * 2];
	rsectorB[i] = new int[nrsB * 2];
	size_dsF[i]=2*nrsF;
	size_dsB[i]=2*nrsB;
	for (int j = 0; j < nrsF; j++)
	{
		rsectorF[i][2 * j + 0] = rsF[j][0];
		rsectorF[i][2 * j + 1] = rsF[j][1];
	}
	for (int j = 0; j < nrsB; j++)
	{
		rsectorB[i][2 * j + 0] = rsB[j][0];
		rsectorB[i][2 * j + 1] = rsB[j][1];
	}
}
}


void Sector::storers()
{
	rsectorF[0] = new int[nrsF * 2];
	rsectorB[0] = new int[nrsB * 2];
	size_dsF[0]=2*nrsF;
	size_dsB[0]=2*nrsB;
	for (int j = 0; j < nrsF; j++)
	{
		rsectorF[0][2 * j + 0] = rsF[j][0];
		rsectorF[0][2 * j + 1] = rsF[j][1];
	}
	for (int j = 0; j < nrsB; j++)
	{
		rsectorB[0][2 * j + 0] = rsB[j][0];
		rsectorB[0][2 * j + 1] = rsB[j][1];
	}
	//printf("%03d %03d\n",2*nrsF,2*nrsB);
}



// When fullstore active, viewshed geometry is saved
//
void Sector::recordsectorRS(char *  name)
{
	FILE *fs;
	//const char *c = name.c_str();
        rsectorF[orderiAT[N - 1]] = new int[0];
        rsectorF[orderiAT[N - 1]] = new int[0];
        size_dsB[orderiAT[N - 1]]=0;
        size_dsB[orderiAT[N - 1]]=0;
	fs=fopen(name,"wb");
	
	for (int i = 0;i<N;i++)
	{
		int j = size_dsF[i];
		fwrite(&j,4,1,fs);
		for (int k = 0;k <j;k++)
		{
			fwrite(&rsectorF[i][k],4,1,fs);
		}
		
		j = size_dsB[i];
		fwrite(&j,4,1,fs);
		for (int k = 0;k<j;k++)
		{
			fwrite(&rsectorB[i][k],4,1,fs);
		}
		
	}
	
	fclose(fs);
	
	for(int i=0; i<N; i++)
	{
		delete rsectorF[i];
		delete rsectorB[i];
	}
		
}











	// ============================================================================================
	// ==========================  Kernels         =======================================
	// ============================================================================================

  // Sector oriented sweep algorithm can be used with different kernels
  // Wrapping sweep procedure, kernel procedure and additional functions:


void Sector::kernel_select(int &i, int type)
{
switch(type)
{
	case 0:
	break;

	case 1:
	sweepS(i);
	if(!nostore)storers(orderiAT[i],false,fullstore);
	break;

	case 2:
            for(int j=0;j<ntowers;j++)
        	if(orderiAT[i]==towersloc[j]){
		//printf("Towerfound\n");
		sweepS(i);
                towerloc=towersloc[j];
		save_cover();
		//i=dim;
		}
	break;

	case 3:
	sweepS_vec(i);
	break;

	case 4:
	sweepV(i);
	if(!nostore)storers(orderiAT[i],true,fullstore);
	break;

	case 5:
	sweepSmask(i);
	if(!nostore)storers(orderiAT[i],false,fullstore);
	break;

      	case 6:
	sweepVmask(i);
	if(!nostore)storers(orderiAT[i],true,fullstore);
	break;
        
	case 7:
        horizonF=horizonB=0;
	sweepS(i);
        //0 Max
        //1 aritmetica
        //2 Harmonica
        //3 Geométrica
        if(isolationindex==0)isolatedMatrix[orderiAT[i]]=max(isolatedMatrix[orderiAT[i]],max(horizonF,horizonB));
        if(isolationindex==1)isolatedMatrix[orderiAT[i]]=isolatedMatrix[orderiAT[i]]+(horizonF)+(horizonB);
        if(isolationindex==2)isolatedMatrix[orderiAT[i]]=isolatedMatrix[orderiAT[i]]+(1.0/horizonF)+(1.0/horizonB);
        if(isolationindex==3)isolatedMatrix[orderiAT[i]]=isolatedMatrix[orderiAT[i]]+log(max(1.0F,horizonF))+log(max(1.0F,horizonB));
        //if(orderiAT[i]==3502550)printf("%d\n",(int)exp(isolatedMatrix[orderiAT[i]]/2));
	break;
        
	case 8:
        horizonhF=horizonhB=0;
	sweepS(i);
        horizonsMatrix[2*orderiAT[i]+0]=max((unsigned char)0,horizonhF);
        horizonsMatrix[2*orderiAT[i]+1]=max((unsigned char)0,horizonhB);
        //if(orderiAT[i]==3502550)printf("%d\n",(int)isolatedMatrix[orderiAT[i]]);
	break;
        
        case 10:
            for(int j=0;j<ntowers;j++)//
            {
                //int j=50;
                int a=towersloc[j];
                int b=towersloc[j-1];
                int c=orderiAT[i];
        	//if(alineados(a,b,c)){
                towerloc=towersloc[j];
                idxtowerloc=j;
		if(towerloc==c){
                    printf("Tower %d found at %d\n",idxtowerloc,towerloc);
                    sweepColores(i);
                    }
		//save_cover();
		//i=dim;
		//}
            }
        break;

}
}


	// ============================================================================================
	// ==========================  Kernel 1: Simple viewshed
	// ============================================================================================



void Sector::sweepS(int i)
{
	float d = MyLinkedList.LL[POV].Value.d;
	float h = MyLinkedList.LL[POV].Value.h + obs_h;


	//printf("i:  %d  sweeppt:  %d\n",i,sweeppt);
	int sweep=presweepF();

#ifdef USE_CLOSEST_POINTS
	float mindist;
if (POV != MyLinkedList.Last){
	get_cls_params(MyLinkedList.LL[POV].Value.idx, 0, visible, max_angle, mindist);
	if(!visible)rsF[nrsF][1]=rsF[nrsF][0];
	while((sweep!=-1)&&( MyLinkedList.LL[sweep].Value.d - d<mindist))sweep = MyLinkedList.LL[sweep].next;
}
#endif

	if ((POV != MyLinkedList.Last)&&(sweep!=-1)) do{



                    delta_d = MyLinkedList.LL[sweep].Value.d - d; 
                    delta_h = MyLinkedList.LL[sweep].Value.h - h;
                    sweeppt = MyLinkedList.LL[sweep].Value.idx;
                    kernelS(rsF,nrsF,sur_dataF);
			//if(espt==0)printf("%f\t%f\t%f\t%d\n",delta_d,delta_h,h,visible);
                    sweep = MyLinkedList.LL[sweep].next;
	}while(sweep!=-1);

	closeprof(visible,true,false,false);
        horizonF=(nrsF>0)?step*(nodes[rsF[nrsF-1][1]].d-d):1;
        horizonhF=(unsigned char)(atan(max_angle)*tograd);
        //horizonF=(nrsF>0)?distance(MyLinkedList.LL[POV].Value.idx,rsF[nrsF-1][1]):0;
        //horizonF=step*sqrt(horizonF);
	if(sweeppt==-1){
	printf("i:  %d  sweeppt:  %d\n",i,sweeppt);
	printf("%f %f\n",d,h);
	printf("%d %f \n",nrsF ,sur_dataF);
	}


	sweep=presweepB();

#ifdef USE_CLOSEST_POINTS
if (POV != MyLinkedList.First){
	get_cls_params(MyLinkedList.LL[POV].Value.idx, 1, visible, max_angle, mindist);
	if(!visible)rsB[nrsB][1]=rsB[nrsB][0];
	while((sweep!=-2)&&(d- MyLinkedList.LL[sweep].Value.d<mindist))sweep = MyLinkedList.LL[sweep].prev;
	}
#endif

	if ((POV != MyLinkedList.First)&&(sweep!=-2)) do{
                    delta_d = d-MyLinkedList.LL[sweep].Value.d; 
                    delta_h = MyLinkedList.LL[sweep].Value.h - h;
                    sweeppt = MyLinkedList.LL[sweep].Value.idx;
                    kernelS(rsB,nrsB, sur_dataB);
                    sweep = MyLinkedList.LL[sweep].prev;
	}while(sweep!=-2);
	//printf("i:  %d  sweeppt:  %d\n",i,sweeppt);

	closeprof(visible,false,false,false);
        horizonB=(nrsB>0)?step*(d-nodes[rsB[nrsB-1][1]].d):1;
        horizonhB=(unsigned char)(atan(max_angle)*tograd);
        //horizonB=(nrsB>0)?distance(MyLinkedList.LL[POV].Value.idx,rsB[nrsB-1][1]):0;
        //horizonB=step*sqrt(horizonB);
        
}



	// ============================================================================================
	// ==========================  Kernel 1: Simple viewshed
	// ============================================================================================



void Sector::sweepColores(int i)
{
	float d = MyLinkedList.LL[POV].Value.d;
	float h = MyLinkedList.LL[POV].Value.h + obs_h;
        forward=true;
	int sweep=presweepF();
	if ((POV != MyLinkedList.Last)&&(sweep!=-1)) do{
                    delta_d = MyLinkedList.LL[sweep].Value.d - d; 
                    delta_h = MyLinkedList.LL[sweep].Value.h - h;
                    sweeppt = MyLinkedList.LL[sweep].Value.idx;
                    kernelColores();
			//if(espt==0)printf("%f\t%f\t%f\t%d\n",delta_d,delta_h,h,visible);
                    sweep = MyLinkedList.LL[sweep].next;
	}while(sweep!=-1);

	closeprof(visible,true,false,false);

        forward=false;
	sweep=presweepB();
	if ((POV != MyLinkedList.First)&&(sweep!=-2)) do{
                    delta_d = d-MyLinkedList.LL[sweep].Value.d; 
                    delta_h = MyLinkedList.LL[sweep].Value.h - h;
                    sweeppt = MyLinkedList.LL[sweep].Value.idx;
                    kernelColores();
                    sweep = MyLinkedList.LL[sweep].prev;
	}while(sweep!=-2);
	//printf("i:  %d  sweeppt:  %d\n",i,sweeppt);

	closeprof(visible,false,false,false);
}








        inline void Sector::kernelS(int rs[][2], int &nrs, float &sur_data)
        {
        float angle = delta_h/delta_d;
        bool above = (angle > max_angle);
        bool opening = above && (!visible);
        bool closing = (!above) && (visible);
        visible = above;
        max_angle = max(angle, max_angle);
        if (opening)
        {
                open_delta_d = delta_d;
                nrs++;
                rs[nrs][0] = sweeppt;
        }
        if (closing)
        {
                sur_data += (delta_d * delta_d - open_delta_d * open_delta_d);
                rs[nrs][1] = sweeppt;
        }
        }

    
        float Sector::intoringsector(int point, float mind, float maxd)
                {
                int tx=towerloc/dimy;
                int ty=towerloc%dimy;
                int x=point/dimy;
                int y=point%dimy;
                float angle =atan2((float)(y - ty),(float) (x - tx));
                if(angle<0)angle+=2*PI_F;
                angle/=torad;
                bool insector=(forward)?belongtoFsector(angle):belongtoBsector(angle);
                if(!insector)return -1;
                float d = sqrt(distance(tx, x, ty, y));
                return ((d<mind)||(d>maxd))?-1:d;
                }

        
        
        void Sector::anotarquienmeve(float distancia, int nodo)
        {
            //no olvidar inicializar con menos uno
            //Buscar si la distancia es inferior a aalguna de las 5 guardadas, o existe algún -1
            //crear sección crítica
            //reemplazar
            //destruir seccion critica
            float maxd=-100;
            int idx=-1;
            int i;
            int cte=nodo*NUMNC;
            for(i=0;i<NUMNC;i++)
                if(suelocolores[cte+i].idx==-1){ //Todavía no hay 5 posiciones del dron que ven este lugar
                    suelocolores[cte+i]={idxtowerloc,distancia};
                    return;}
                else
                    if(suelocolores[cte+i].d>maxd){idx=i;maxd=suelocolores[cte+i].d;} //busco el indice más alejado
            if(maxd>distancia)
                suelocolores[cte+idx]={idxtowerloc,distancia};
            
            
            
        }
        
        void Sector::anotarquienmeve2(float d1,float d2)
        {
            float d=d1;
            for(int j=0;j<N;j++) 
                if((d=intoringsector(j,d1,d2))>=0)
                {
                float maxd=-100;
                int idx=-1;
                int i;
#pragma omp critical
                for(i=0;i<NUMNC;i++)
                    if(suelocolores[j*NUMNC+i].idx==-1){ //Todavía no hay 5 posiciones del dron que ven este lugar
                        suelocolores[j*NUMNC+i]={idxtowerloc,d};
                        break;}
                    else
                        if(suelocolores[j*NUMNC+i].d>maxd){idx=i;maxd=suelocolores[j*NUMNC+i].d;} //busco el indice más alejado
                if((maxd>d)&&(i==NUMNC))
                    suelocolores[j*NUMNC+idx]={idxtowerloc,d};                    
                }
        }
    
    
    
        inline void Sector::kernelColores()
        {
        float angle = delta_h/delta_d;
        bool above = (angle > max_angle);
        bool opening = above && (!visible);
        bool closing = (!above) && (visible);
        visible = above;
        max_angle = max(angle, max_angle);
        if (opening)
                open_delta_d = delta_d;
        if (closing)
            anotarquienmeve2(open_delta_d,delta_d);
        }
        
        
        

	// ============================================================================================
	// ==========================  Kernel 2: Simple viewshed - vectorized 
	// ============================================================================================

  //Very expensive. Vectorization requires additional storage and
  //produces cache faults

void Sector::sweepS_vec(int i)
{
	float d = MyLinkedList.LL[POV].Value.d;
	float h = MyLinkedList.LL[POV].Value.h + obs_h;
	cnt_vec=0;
	int pos;
	int sweep=presweepF();
	if (POV != MyLinkedList.Last) do{
		    //sweeppt_vec[cnt_vec]=MyLinkedList.LL[sweep].Value.idx;
		    delta_d_vec[cnt_vec]=MyLinkedList.LL[sweep].Value.d;
		    height_vec[cnt_vec++]=MyLinkedList.LL[sweep].Value.h;
                    sweep = MyLinkedList.LL[sweep].next;
	}while(sweep!=-1);
	pos=cnt_vec;
	sweep=presweepB();
	if (POV != MyLinkedList.First) do{
		    //sweeppt_vec[cnt_vec]=MyLinkedList.LL[sweep].Value.idx;
		    delta_d_vec[cnt_vec]=MyLinkedList.LL[sweep].Value.d;
		    height_vec[cnt_vec++]=MyLinkedList.LL[sweep].Value.h;
                    sweep = MyLinkedList.LL[sweep].prev;
	}while(sweep!=-2);

//	for(int j=0;j<BW-1;j++)

	return;
	visible=true;
	max_angle=-2000;
	open_delta_d=0;
	delta_d=0;
#pragma ivdep
	for(int j=0;j<pos;j++){
		delta_d_vec[j]-=d;
		delta_d_vec[j]/=(height_vec[j]-h);
		/*
                bool above = (angle > max_angle);
                bool opening = above && (!visible);
                bool closing = (!above) && (visible);
                visible = above;
        	if (opening)
        	{
                	open_delta_d = delta_d;
                	nrsF++;
                	rsF[nrsF][0] = MyLinkedList.LL[sweeppt_vec[j]].Value.idx;
        	}
        	if (closing)
        	{
                	sur_dataF += (delta_d * delta_d - open_delta_d * open_delta_d);
                	rsF[nrsF][1] =  MyLinkedList.LL[sweeppt_vec[j]].Value.idx;
        	}*/
		}
	closeprof(visible,true,false,false);
	visible=true;
	max_angle=-2000;
	open_delta_d=0;

	delta_d=0;
#pragma ivdep
	for(int j=pos;j<cnt_vec;j++){
		delta_d_vec[j]-=d;
		delta_d_vec[j]/=(height_vec[j]-h);
		//float delta_d=(d-MyLinkedList.LL[sweeppt_vec[j]].Value.d);
		//height_vec[j]=(MyLinkedList.LL[sweeppt_vec[j]].Value.h-h)/delta_d;
		/*
                bool above = (angle > max_angle);
                bool opening = above && (!visible);
                bool closing = (!above) && (visible);
                visible = above;
 		if (opening)
        	{
                	open_delta_d = delta_d;
                	nrsB++;
                	rsB[nrsB][0] = MyLinkedList.LL[sweeppt_vec[j]].Value.idx;
        	}
        	if (closing)
        	{
                	sur_dataB += (delta_d * delta_d - open_delta_d * open_delta_d);
                	rsB[nrsB][1] =  MyLinkedList.LL[sweeppt_vec[j]].Value.idx;
        	}*/
		}
	closeprof(visible,false,false,false);



/*
*/
}


  void Sector::kernelS_vec(int rs[][2], int &nrs, float &sur_data)
  {
  bool opening = above_vec[cnt_vec] && (!visible);
  bool closing = (!above_vec[cnt_vec]) && (visible);
  visible = above_vec[cnt_vec];
  if (opening)
  {
          open_delta_d = delta_d_vec[cnt_vec];
          nrs++;
          rs[nrs][0] = sweeppt_vec[cnt_vec];
  }
  if (closing)
  {
          sur_data += (delta_d_vec[cnt_vec] * delta_d_vec[cnt_vec] - open_delta_d * open_delta_d);
          rs[nrs][1] = sweeppt_vec[cnt_vec];
  }
  }
  
  
	void Sector::init_vec(int bw)
	{
	above_vec= new bool[bw];
	delta_d_vec= new float[bw];
	height_vec= new float[bw];
	sweeppt_vec= new int[bw];
	}
	void Sector::dispose_vec()
	{
	delete above_vec;
	delete delta_d_vec;
	delete sweeppt_vec;
	delete height_vec;
	}


	// ============================================================================================
	// ==========================  Kernel 3: Volumetric viewshed
	// ============================================================================================


void Sector::sweepV(int i)
{
	float d = MyLinkedList.LL[POV].Value.d;
	float h = MyLinkedList.LL[POV].Value.h + obs_h;
	int sweep;

	visible=true;
	nrsF=0;
	sur_dataF=0;
	vol_dataF=0;
	delta_d=0;
	delta_h=0;
	open_delta_d=0;
	open_delta_h=0;
        
	sweep = MyLinkedList.LL[POV].next;
	sweeppt= MyLinkedList.LL[POV].Value.idx;
	rsF[0][0]=sweeppt;
	max_angle=-2000; 

#ifdef USE_CLOSEST_POINTS
	float mindist;
if (POV != MyLinkedList.Last){
	get_cls_params(MyLinkedList.LL[POV].Value.idx, 0, visible, max_angle, mindist);
	if(!visible)rsF[nrsF][1]=rsF[nrsF][0];
	while((sweep!=-1)&&( MyLinkedList.LL[sweep].Value.d - d<mindist))
		sweep = MyLinkedList.LL[sweep].next;
}
#endif        
        
	if ((POV != MyLinkedList.Last) &&(sweep!=-1))
		
		do{
                    delta_d = MyLinkedList.LL[sweep].Value.d - d; 
                    delta_h = MyLinkedList.LL[sweep].Value.h - h;
                    sweeppt = MyLinkedList.LL[sweep].Value.idx;
                    kernelV(rsF, nrsF, sur_dataF, vol_dataF);
                    sweep = MyLinkedList.LL[sweep].next;
	}while(sweep!=-1);
	//if((rsF[0][0]==2906611)&&(vol_dataF!=0)&&(sector==105))printf("%e\n",vol_dataF*volscale);

	closeprof(visible,true,true,false);

	visible=true;
	nrsB=0;
	sur_dataB=0;
	vol_dataB=0;
	delta_d=0;
	delta_h=0;
	open_delta_d=0;
	open_delta_h=0;
        
	sweep = MyLinkedList.LL[POV].prev;
	sweeppt= MyLinkedList.LL[POV].Value.idx;
	rsB[0][0]=sweeppt;
	max_angle=-2000;
        
#ifdef USE_CLOSEST_POINTS
if (POV != MyLinkedList.First){
	get_cls_params(MyLinkedList.LL[POV].Value.idx, 1, visible, max_angle, mindist);
	if(!visible)rsB[nrsB][1]=rsB[nrsB][0];
	while((sweep!=-2)&&(d- MyLinkedList.LL[sweep].Value.d<mindist))sweep = MyLinkedList.LL[sweep].prev;
	}
#endif
        
	if ((POV != MyLinkedList.First) &&(sweep!=-2))
		
		
		do{
                    delta_d = d-MyLinkedList.LL[sweep].Value.d; 
                    delta_h = MyLinkedList.LL[sweep].Value.h - h;
                    sweeppt = MyLinkedList.LL[sweep].Value.idx;
                    kernelV(rsB, nrsB, sur_dataB, vol_dataB);
                    sweep = MyLinkedList.LL[sweep].prev;
	}while(sweep!=-2);

	closeprof(visible,false,true,false);
}


  void Sector::kernelV(int rs[][2], int &nrs, float &sur_data, float &vol_data)
  {
  float idd=1/delta_d;
  float angle = delta_h*idd;
  bool above = (angle > max_angle);
  bool opening = above && (!visible);
  bool closing = (!above) && (visible);
  visible = above;
  max_angle = max(angle, max_angle);
  if (opening)
  {
          open_delta_d = delta_d;
          open_delta_h = delta_h;
          nrs++;
          rs[nrs][0] = sweeppt;
  }
  if (closing)
  {
          sur_data += (delta_d * delta_d - open_delta_d * open_delta_d);
          float mean = (delta_d + open_delta_d);
          vol_data += mean * fabs(open_delta_h * delta_d - open_delta_d * delta_h);
          rs[nrs][1] = sweeppt;
  }
  }


	// ============================================================================================
	// ==========================  Kernel 4: Viewshed with "already covered areas" for coverage
	// ============================================================================================

void Sector::sweepSmask(int i)
{
	float d = MyLinkedList.LL[POV].Value.d;
	float h = MyLinkedList.LL[POV].Value.h + obs_h;
	delta_d=0;


	//printf("i:  %d  sweeppt:  %d\n",i,sweeppt);
	int sweep=presweepF();
	if (POV != MyLinkedList.Last) do{
                    sweeppt = MyLinkedList.LL[sweep].Value.idx;
		    bool alreadycovered=mask[sweeppt];
		    if(alreadycovered&&visible)
                    sur_dataF+=delta_d*delta_d;
                    delta_d = MyLinkedList.LL[sweep].Value.d - d; 
		    if(alreadycovered&&visible)
                    sur_dataF-=delta_d*delta_d;
                    delta_h = MyLinkedList.LL[sweep].Value.h - h;
                    kernelS(rsF,nrsF,sur_dataF);
			//if(espt==0)printf("%f\t%f\t%f\t%d\n",delta_d,delta_h,h,visible);
                    sweep = MyLinkedList.LL[sweep].next;
	}while(sweep!=-1);

	closeprof(visible,true,false,false);

	if(sweeppt==-1){
	printf("i:  %d  sweeppt:  %d\n",i,sweeppt);
	printf("%f %f\n",d,h);
	printf("%d %f \n",nrsF ,sur_dataF);
	}

	delta_d=0;
	sweep=presweepB();
	if (POV != MyLinkedList.First) do{
                    sweeppt = MyLinkedList.LL[sweep].Value.idx;
		    bool alreadycovered=mask[sweeppt];
		    if(alreadycovered&&visible)
                    sur_dataB+=delta_d*delta_d;
                    delta_d = d-MyLinkedList.LL[sweep].Value.d; 
		    if(alreadycovered&&visible)
                    sur_dataB-=delta_d*delta_d;
                    delta_h = MyLinkedList.LL[sweep].Value.h - h;
                    kernelS(rsB,nrsB, sur_dataB);
                    sweep = MyLinkedList.LL[sweep].prev;
	}while(sweep!=-2);
	//printf("i:  %d  sweeppt:  %d\n",i,sweeppt);

	closeprof(visible,false,false,false);
}




	// ============================================================================================
	// ==========================  Kernel 4: Viewshed with "already covered areas" for coverage
	// ============================================================================================


void Sector::sweepVmask(int i)
{
	float d = MyLinkedList.LL[POV].Value.d;
	float h = MyLinkedList.LL[POV].Value.h + obs_h;
        
	visible=true;
	nrsF=0;
	sur_dataF=0;
	vol_dataF=0;
	delta_d=0;
	delta_h=0;
	open_delta_d=0;
	open_delta_h=0;


	//printf("i:  %d  sweeppt:  %d\n",i,sweeppt);
	int sweep=presweepF();
        float prev_delta_d=0;
        float prev_delta_h=0;
	if (POV != MyLinkedList.Last) do{
                    sweeppt = MyLinkedList.LL[sweep].Value.idx;
		    bool alreadycovered=mask[sweeppt];
		    if(alreadycovered&&visible){
                        prev_delta_d=delta_d;
                        prev_delta_h=delta_h;
                        sur_dataF+=delta_d*delta_d;
                    }
                    delta_d = MyLinkedList.LL[sweep].Value.d - d; 
                    delta_h = MyLinkedList.LL[sweep].Value.h - h;
		    if(alreadycovered&&visible){
                        sur_dataF-=delta_d*delta_d;
                        vol_dataF-=(delta_d+prev_delta_d)*fabs(delta_d*prev_delta_h-delta_h*prev_delta_d);
                    }
                    kernelV(rsB, nrsB, sur_dataF, vol_dataF);
			//if(espt==0)printf("%f\t%f\t%f\t%d\n",delta_d,delta_h,h,visible);
                    sweep = MyLinkedList.LL[sweep].next;
	}while(sweep!=-1);

	closeprof(visible,true,true,false);

	if(sweeppt==-1){
	printf("i:  %d  sweeppt:  %d\n",i,sweeppt);
	printf("%f %f\n",d,h);
	printf("%d %f \n",nrsF ,sur_dataF);
	}

	visible=true;
	nrsB=0;
	sur_dataB=0;
	vol_dataB=0;
	delta_d=0;
	delta_h=0;
	open_delta_d=0;
	open_delta_h=0;
        
	sweep=presweepB();
	if (POV != MyLinkedList.First) do{
                    sweeppt = MyLinkedList.LL[sweep].Value.idx;
		    bool alreadycovered=mask[sweeppt];
		    if(alreadycovered&&visible){
                        prev_delta_d=delta_d;
                        prev_delta_h=delta_h;
                        sur_dataB+=delta_d*delta_d;
                    }
                    delta_d = d-MyLinkedList.LL[sweep].Value.d; 
                    delta_h = MyLinkedList.LL[sweep].Value.h - h;
		    if(alreadycovered&&visible){
                        sur_dataB-=delta_d*delta_d;
                        vol_dataB-=(delta_d+prev_delta_d)*fabs(delta_d*prev_delta_h-delta_h*prev_delta_d);
                    }
                    kernelV(rsB, nrsB, sur_dataB, vol_dataB);
                    sweep = MyLinkedList.LL[sweep].prev;
	}while(sweep!=-2);
	//printf("i:  %d  sweeppt:  %d\n",i,sweeppt);

	closeprof(visible,false,true,false);
}











//Distance between two points
float Sector::distance(int a, int b)
{
            float ay = a % dimy;
            float ax = a / dimy;
            float by = b % dimy;
            float bx = b / dimy;
            return (float)((ax - bx) * (ax - bx) + (ay - by) * (ay - by));
}

float Sector::distance(int x1, int x2, int y1, int y2)
{
return (float)((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1));
}

bool Sector::belongtoFsector(float angle)
{
return ((angle>=sect_langle_g)&&(angle<sect_hangle_g)) ;
}

bool Sector::belongtoBsector(float angle)
{
if(sector!=NSECTOR-1)return ((angle>=sect_langle_g+180)&&(angle<sect_hangle_g+180));
else 
{
	if(angle < origin)return true;
	if(angle >= (360-sector_step)+origin )return true; 
	return false;
}
}


//Determine if two points are visible
//
bool Sector::pointeval(int i, int cpx, int cpy)
{
//just for malaga test
//if(i%2000 < 900) return  true;
int x=i/dimy;
int y=i%dimy;
float angle =atan2((float)(y - cpy),(float) (x - cpx));
if(angle<0)angle+=2*PI_F;
angle/=torad;
if(angle<0)printf("%f ",angle);
if(angle>=360)printf("%f ",angle);

bool FS=belongtoFsector(angle);
bool BS=belongtoBsector(angle);
if(!FS&&!BS)return false;

float d = distance(cpy, y, cpx, x);
if(FS)
for (int k1 = 0; k1 < nrsF; k1++)
      if ((d >= lsF[k1*2]) && (d <= lsF[k1*2+1])) {cntcoverage++;return true; }

if(BS)
for (int k1 = 0; k1 < nrsB; k1++)
      if ((d >= lsB[k1*2]) && (d <= lsB[k1*2+1]))  {cntcoverage++;return true; }
return false;
}



float Sector::pointvolume(int i, int cpx, int cpy,bool *coverage)
{
//just for malaga test
//if(i%2000 < 900) return  true;
int x=i/dimy;
int y=i%dimy;
float angle =atan2((float)(y - cpy),(float) (x - cpx));
if(angle<0)angle+=2*PI_F;
angle/=torad;
if(angle<0)printf("%f ",angle);
if(angle>=360)printf("%f ",angle);
bool FS=belongtoFsector(angle);
bool BS=belongtoBsector(angle);
if(!FS&&!BS)
return 0;

double d = distance(cpx, x, cpy, y);
double sd,techoB,techoF,suelo;
float H=0;
if(FS){
for (int k1 = 0; k1 < nrsF; k1++){
      if ((d >= lsF[k1*2]) && (d <= lsF[k1*2+1])) { //squared values
                  sd=sqrt(d);
                  techoF=towerheight+sd*ceil_angle_F;
                  suelo=heights[i];
		  cntcoverage++;
		  coverage[i]=true;
		  H=techoF-suelo;
		  return (H>0)?H:0; 
	  }
      if(k1>0)if ((d > lsF[k1*2-1]) && (d < lsF[k1*2])){//Only air mass is visible
                  sd=sqrt(d);
                  techoF=towerheight+sd*ceil_angle_F;
                  suelo=heights[i];
		  double floor_angle_F=(heights[rsF[k1][0]]-towerheight)/sqrt(lsF[k1*2]);
                  suelo=max(suelo,towerheight+sd*floor_angle_F); //Max, because sometimes, BOS doesn't represents sector
		  suelo=towerheight+sd*floor_angle_F;
		  H=techoF-suelo;
		  return (H>0)?H:0; 
	  }
}
}
if(BS){
for (int k1 = 0; k1 < nrsB; k1++){
      if ((d >= lsB[k1*2]) && (d <= lsB[k1*2+1]))  {
                  sd=sqrt(d);
                  techoB=towerheight+sd*ceil_angle_B;
		  suelo=heights[i];
		  cntcoverage++;
		  coverage[i]=true;
                  H=techoB-suelo;
		  return (H>0)?H:0; 
	  }
      if(k1>0)if ((d > lsB[k1*2-1]) && (d < lsB[k1*2])){//Only air mass is visible
                  sd=sqrt(d);
                  techoB=towerheight+sd*ceil_angle_B;
                  suelo=heights[i];
		  double floor_angle_B=(heights[rsB[k1][0]]-towerheight)/sqrt(lsB[k1*2]);
                  suelo=max(suelo,towerheight+sd*floor_angle_B); //Max, because sometimes, BOS doesn't represents sector
                  suelo=towerheight+sd*floor_angle_B;
		  H=techoB-suelo;
		  return (H>0)?H:0; 
	  }
}
}
return 0;
}


void Sector::save_cover()
{
		storers();
		cntcoverage=0;
		int cpx=towerloc/dimy;
		int cpy=towerloc%dimy;
		lsF=new float[2*nrsF]; //squared distance with step as unit
		lsB=new float[2*nrsB];
		for(int j=0;j<nrsF;j++){lsF[2*j]=distance(towerloc,rsF[j][0]);lsF[2*j+1]=distance(towerloc,rsF[j][1]);}
		for(int j=0;j<nrsB;j++){lsB[2*j]=distance(towerloc,rsB[j][0]);lsB[2*j+1]=distance(towerloc,rsB[j][1]);}
		if(volume)
		{
			towerheight=obs_h+heights[towerloc];
			if(nrsF>0){
			int horizonF=rsF[nrsF-1][1];
			double farhF=heights[horizonF];
			double distF=sqrt(lsF[2*nrsF-1]);
			ceil_angle_F=(farhF-towerheight)/distF;
			}
			if(nrsB>0){
			int horizonB=rsB[nrsB-1][1];
			double farhB=heights[horizonB];
			double distB=sqrt(lsB[2*nrsB-1]);
			ceil_angle_B=(farhB-towerheight)/distB;
			}
			for(int j=0;j<N;j++) 
				if((!inputmask)||(mask[j]==0))
					{
					//volcoverage[j]=max(volcoverage[j],pointvolume(j,cpx,cpy,coverage));
					//if(sector==169&&j==2684662)
					volcoverage[j]=max(volcoverage[j],pointvolume(j,cpx,cpy,coverage));
                                        }
		}
		else
			for(int j=0;j<N;j++) if(pointeval(j,cpx,cpy))coverage[j]=true;
		delete lsF;
		delete lsB;
		//printf("%d \n",cntcoverage);
}




	// ============================================================================================
	// ==========================  Kernels: Common functions:
	// ============================================================================================




int Sector::presweepF()
{
visible=true;
nrsF=0;
sur_dataF=0;
delta_d=0;
delta_h=0;
open_delta_d=0;
open_delta_h=0;
int sweep = MyLinkedList.LL[POV].next;
sweeppt= MyLinkedList.LL[POV].Value.idx;
rsF[0][0]=sweeppt;
max_angle=-2000; 
return sweep;
}

int Sector::presweepB()
{
visible=true;
nrsB=0;
sur_dataB=0;
delta_d=0;
delta_h=0;
open_delta_d=0;
open_delta_h=0;
int sweep = MyLinkedList.LL[POV].prev;
sweeppt= MyLinkedList.LL[POV].Value.idx;
rsB[0][0]=sweeppt;
max_angle=-2000; 
return sweep;
}



void Sector::closeprof(bool visible, bool fwd, bool vol, bool full)
{
if(fwd) {
	nrsF++;
	if(visible)
		{
		rsF[nrsF-1][1]=sweeppt;
		sur_dataF+=delta_d*delta_d-open_delta_d*open_delta_d;
		if(vol)vol_dataF+=(delta_d+open_delta_d)*fabs(delta_d*open_delta_h-delta_h*open_delta_d);
		}
	if((nrsF==1)&&(delta_d<1.5F))nrsF=0;
	}
else    {
	nrsB++;
	if(visible)
		{
		rsB[nrsB-1][1]=sweeppt;
		sur_dataB+=delta_d*delta_d-open_delta_d*open_delta_d;
		if(vol)vol_dataB+=(delta_d+open_delta_d)*fabs(delta_d*open_delta_h-delta_h*open_delta_d);
		}
	if((nrsB==1)&&(delta_d<1.5F))nrsB=0; // safe data with small vs
	}
}











