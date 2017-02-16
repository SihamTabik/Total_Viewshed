/*
AUTHORS: 
*SIHAM TABIK
*LUIS FELIPE ROMERO GOMEZ
*ANTONIO MANUEL RODRIGUEZ CERVILLA 
DEPT: ARQUITECTURA DE LOS COMPUTADORES
UNIVERSIDAD DE MALAGA
 
*/

#include "Defs.h"
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <iostream>
#include <sstream>
#include <sys/stat.h>
#include "AuxFunc.h"
#include "lodepng.h"

#define IS_BIG_ENDIAN (*(uint16_t *)"\0\xff" < 0x100)

using namespace std;
#ifdef WIN32
void create_directory(const char *name)
{
	CreateDirectoryA(name,NULL);
}
#else
void create_directory(const char *name)
{
	mkdir(name,S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
}
#endif



//==================================================================================
//	Filenames
//==================================================================================


//string pathfrsect = "./cached/";	//Directory ring sector file (fullstore only)
char filename[100]; // Input filename
char base_name[100]; // Input filename without path and extension
char maskname[100]; // Mask file with boolean data
char listname[100]; // List of tower locations

char surname[100]; // Total viewshed result
char surlname[100]; // Total viewshed result
char surcname[100]; // Volumetric visibility results single point
char volname[100]; // Volumetric visibility results
char vollname[100]; // Volumetric visibility results 
char volcname[100]; // Volumetric visibility results single point
char horzname[100]; // Total viewshed result

char pngfname[100]; // Graphical output filename
char pngname[100];  // Same, without path


char pngvfname[100]; // Graphical output filename - volume
char pngvname[100];  // Same, without path - volume
char cvrgfname[100]; // Coverage output
char idname[100]; // Track or tower identifier



extern int dim,step,utm_n,utm_e,dimx,dimy;
extern int mode,identifier;
extern int bw;
extern bool verbose,silent,inputmask,volume,fullstore,showtime;
extern int tower,*towers,ntowers,cnttower;
extern double obsheight;
extern bool floatsDEM, endianOK;
bool towernotgiven=true;



		

//------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------
//				Useful functions
//------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------


bool exists(char * filename) {
    
    struct stat fileinfo;
    return !stat(filename, &fileinfo);
}

long filesize(char *filename)
{
	long size=0;
    struct stat fileinfo;
	if(!stat(filename, &fileinfo))size=fileinfo.st_size;
	return size;
}

string converInt(int number){ 
string Result;			
ostringstream convert; 
	convert<<number;		
	Result=convert.str();
	if (number<10){
		Result= "00"+ Result;
	}else if (number<100){
		Result= "0"+ Result;
	}
return Result; 
}


#ifdef WIN32
double dtime() {
double tseconds = 0.0;
FILETIME  ft;
GetSystemTimeAsFileTime(&ft);
    unsigned long long tt = ft.dwHighDateTime;
    tt <<=32;
    tt |= ft.dwLowDateTime;
    double  x =tt ;
	x/=10000000;
	return x;
}
#else
double dtime() {
double tseconds = 0.0;
struct timeval mytime;
	gettimeofday(&mytime,(struct timezone *) 0);
	tseconds = (double) (mytime.tv_sec + (double)mytime.tv_usec * 1.0e-6);
return (tseconds) ;
}
#endif


//------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------
//                           Graphical output section
//------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------
color *my_palette;

std::vector<unsigned char>  data2vectColores(nodoscolores *nc)
{
std::vector<unsigned char> image;
image.resize(dimx * dimy * 4);


for(int i = 0; i < dimx*dimy; i++)
  {
    int ii=N2T(i);
    float best=-1;
    int ibest=ii*NUMNC;
    for(int j=0;j<NUMNC;j++)
        if(nc[ii*NUMNC+j].d>best)
        {best=nc[ii*NUMNC+j].d;ibest=ii*NUMNC+j;}
    if(best==-1)
    {
    image[4 * i + 0] = 0;
    image[4 * i + 1] = 0;
    image[4 * i + 2] = 0;
    image[4 * i + 3] = 0;
    }
    else
    {
    int v=1024.0*(nc[ibest].idx/256.0);
    if(v>1023)v=1023;
    //int x=(i%dimx)*dimy+i/dimx;
    image[4 * i + 0] = my_palette[v].R;
    image[4 * i + 1] = my_palette[v].G;
    image[4 * i + 2] = my_palette[v].B;
    image[4 * i + 3] = 255;
    }
  }
for(int i=0;i<ntowers;i++){
    
}
return image;
}

//image tiene el orden natural
std::vector<unsigned char>  data2vect(float * data, float maxval)
{
std::vector<unsigned char> image;
image.resize(dimx * dimy * 4);
for(int i = 0; i < dimx*dimy; i++)
  {
    int v=1024.0*(data[i]/maxval);
    if(v<0)v=0;
    if(v>1023)v=1023;
    //int x=(i%dimx)*dimy+i/dimx;
    image[4 * i + 0] = my_palette[v].R;
    image[4 * i + 1] = my_palette[v].G;
    image[4 * i + 2] = my_palette[v].B;
    image[4 * i + 3] = 255;
  }
return image;
}

void saveFileColores(nodoscolores *nc, char *fname)
{
char ffname[100];    
sprintf(ffname,"%s/%s.flt",OUTPUT,fname);
FILE * f;
int n=0;
f=fopen(ffname,"wb");
nodoscolores x1,x2;
for(int i=0;i<dim;i++){
        int dist=-1;
        x2=nc[NUMNC*N2T(i)];
        for(int j=0;j<NUMNC;j++){
                x1=nc[NUMNC*N2T(i)+j];
                if(x1.d>dist)
                    {x2=x1;
                    dist=nc[NUMNC*N2T(i)+j].d;
                    }
        }
        fwrite(&x2,sizeof(nodoscolores),1,f);
}
fclose(f);
}




void saveColores(nodoscolores *nc)
{
my_palette=build_palette();
saveFileColores(nc,"datosColoreado.dat");
//saveFileColores(nc,NUMNC,"datosColoreado5.dat");
std::vector<unsigned char> image=data2vectColores(nc);
printf("Almacenando archivo png de colores del territorio\n");
lodepng::encode("output/prueba.png", image, dimx, dimy);
}

void data2image(const char* filename, float *data ,float maxval)
{
my_palette=build_palette();
//printf("aqui\n");
std::vector<unsigned char> image=data2vect(data,maxval);
lodepng::encode(filename, image, dimx, dimy);
delete my_palette;
}

float interp( float val, float y0, float x0, float y1, float x1 ) {
    return (val-x0)*(y1-y0)/(x1-x0) + y0;
}

bool alineados(int a, int b, int c)
{
    if(c==a)return true;
    if(c==b)return true;
    int ay=a%dimy;
    int ax=a/dimx;
    int by=b%dimy;
    int bx=b/dimx;
    int cy=c%dimy;
    int cx=c/dimx;
    if((ax-cx)*(cx-bx)<0)return false;
    if((ay-cy)*(cy-by)<0)return false;
    int ac=atan2(ax-cx,ay-cy);
    int cb=atan2(cx-bx,cy-by);
    if(abs(ac-cb)<0.01745*2)return true;
    return false;
}

float base( float val ) {
    if ( val <= -0.75 ) return 0;
    else if ( val <= -0.25 ) return interp( val, 0.0, -0.75, 1.0, -0.25 );
    else if ( val <= 0.25 ) return 1.0;
    else if ( val <= 0.75 ) return interp( val, 1.0, 0.25, 0.0, 0.75 );
    else return 0.0;
}

float red( float gray ) {
    return 256*base( gray - 0.5 );
}
float green( float gray ) {
    return 256*base( gray );
}
float blue( float gray ) {
    return 256*base( gray + 0.5 );
}
color get_color_jet(int index)
{
color c ;
c.R=(int)red(index/1024.0);
c.G=(int)green(index/1024.0);
c.B=(int)blue(index/1024.0);
return c;
}

color get_color_jet2(int index)
{
   float v=index;
   float dv=1024;
   color c={255,255,255};
   if (v < (0.25 * dv)) {
      c.R = 0;
      c.G = 256*(4 * v / dv);
   } else if (v < (0.5 * dv)) {
      c.R = 0;
      c.B = 256*(1 + 4 * (0.25 * dv - v) / dv);
   } else if (v < (0.75 * dv)) {
      c.R = 256*(4 * (v - 0.5 * dv) / dv);
      c.B = 0;
   } else {
      c.G = 256*(1 + 4 * (0.75 * dv - v) / dv);
      c.B = 0;
   }
   return(c);
}

color * build_palette()
{
color* cm = new color[1024];
for (int i = 0; i < 1024; i++)
	cm[i] = get_color_jet2(i);
return cm;
}



//------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------
// 			Mask file read/write (Cell tower location kernel)
//------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------


//This function saves a boolean array
int saveCover(bool * c, char *fname,bool append)
{
char ffname[100];    
sprintf(ffname,"%s/%s.bil",OUTPUT,fname);    
FILE * f;
int n=0;
if(append)f=fopen(ffname,"ab");
else f=fopen(ffname,"wb");
for(int i=0;i<dim;i++){
	unsigned char a=0;
	if(c[N2T(i)]){a=255;n++;}
	fwrite(&a,1,1,f);
	}
fclose(f);
return n;
}

void saveFloats(float * c, char *fname,bool transpose)
{
char ffname[100];    
sprintf(ffname,"%s/%s.flt",OUTPUT,fname);
FILE * f;
int n=0;
f=fopen(ffname,"wb");
if(transpose)
    for(int i=0;i<dim;i++){
            float x=c[N2T(i)];
            fwrite(&x,1,4,f);
            }
else            fwrite(c,4, dim,f);

fclose(f);
}

void saveHorizons(unsigned char * c, char *fname)
{
char ffname[100];    
sprintf(ffname,"%s/%s.bin",OUTPUT,fname);
FILE * f;
int n=0;
f=fopen(ffname,"wb");
fwrite(c,1, 360*dim,f);
fclose(f);
}


void resetCover(bool* c)
{
for(int i=0;i<dim;i++)c[i]=false;
}

void resetVolCover(float* c)
{
for(int i=0;i<dim;i++)c[i]=0;
}




bool readCover(bool * c, char *fname)
{
try{
    FILE * f;
    f=fopen(fname,"rb");
    //Data stored as bytes with natural ordering (west to east inner loop))
    for(int i=0;i<dim;i++)
            fread(&c[N2T(i)],1,1,f);
    fclose(f);
    return true;
}catch(int E)
    {
    return false;
    }
}


//------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------
//			Height file 
//------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------

double *readHeights(char * nmde, bool floats, bool sameendian)
{//Read heights of the points on the map and positions according to order A
        FILE * f;
        double *h=new double[dim*sizeof(double)];
        f=fopen(nmde,"rb");
        if(f==NULL)
        {
		printf("Error opening %s\n",nmde);
        }else
        {
			if(sameendian&&!floats)
                for(int i=0;i<dim;i++)
                {
				        short num;
                        fread(&num,2,1,f);
                        h[N2T(i)]=(1.0*num)/step; //internal representation from top to bottom (inner loop)
                }
			if(sameendian&&floats)
                for(int i=0;i<dim;i++)
                {
						float num;
                        fread(&num,4,1,f);
                        h[N2T(i)]=num/step; //internal representation from top to bottom (inner loop)
                }
        }
        return h;
}

double *readHeights(char * nmde)
{//Read heights of the points on the map and positions according to order A
        FILE * f;
        double *h=new double[dim*sizeof(double)];
        f=fopen(nmde,"rb");
        if(f==NULL)
        {
		printf("Error opening %s\n",nmde);
        }else
        {
                for(int i=0;i<dim;i++)
                {
		short num;
                        fread(&num,2,1,f);
                        h[i]=(1.0*num)/step; //internal representation from top to bottom (inner loop)
                }
        }
        return h;
}


//------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------
//			Headr for flt arcgis
//------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------
void writeheader(int n, int e, int step, char *fname,bool isbool)
{
char ffname[100];    
sprintf(ffname,"%s/%s.hdr",OUTPUT,fname);
    
char *endian="LSBFIRST";
if(IS_BIG_ENDIAN)strcpy(endian,"MSBFIRST");
FILE *f=fopen(ffname,"wt");
fprintf(f,"NROWS        %d\n", dimx);
fprintf(f,"NCOLS        %d\n", dimy);
fprintf(f,"XLLCORNER    %d\n", e);
fprintf(f,"YLLCORNER    %d\n", n-step*dimy);
fprintf(f,"CELLSIZE     %d\n", step);
if(!isbool)
fprintf(f,"NODATA_VALUE -1.0\n");
else 
fprintf(f,"NODATA_VALUE 0\n");
if(!isbool)
fprintf(f,"BYTEORDER    %s\n",endian);
else 
fprintf(f,"NBITS        8\n");
fclose(f);
}

//------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------
//			Get point with max coverage, for cell tower kernel
//------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------

int maxSurf(  char * nmds)
{
char fullfilename[100];
sprintf(fullfilename,"%s/%s.flt",OUTPUT,nmds);

FILE * f;
        float num;
        float maxnum=-10000;
        int idx=0;
        f=fopen(nmds,"rb");
        if(f==NULL)
        {
		if(verbose)printf("Error opening existing %s/%s.flt\n",OUTPUT,nmds);
                return 0;
        }else
        {
                for(int i=0;i<dim;i++)
                {
                        fread(&num,4,1,f);
                        if(num<0)num=0;
                        if(num>maxnum){ //input with natural ordering
                                maxnum=num;idx=N2T(i);
                                }
                }
                if(verbose)printf("Max %d value: %f\n",idx,maxnum);
        }
        return idx;
}

float * readFloats(  char  * nmds)
{
float *f=new float[dim];
FILE *F=fopen(nmds,"rb");
fread(f,dim,4,F);
fclose(F);
return f;
}




void readTowers()
{
sprintf(listname,"%s/list.dat",INPUTDIR);
if(verbose)printf("Reading towers from %s\n",listname);
FILE *F=fopen(listname,"rb");
//int n;
if(F)
fread(&ntowers,1,4,F);
if(verbose)printf("%d towers readed\n",ntowers);
towers=new int[ntowers];
fread(towers,ntowers,4,F);
fclose(F);
if(verbose)printf("%d towers readed. First, at %d\n",ntowers,towers[0]);
}



int get_point_id(int n, int e)
{
    return (utm_n-n)/step + ((e-utm_e)/step)*dimy;
    
}
int get_point_n(int pt)
{
    return  utm_n-(pt%dimy)*step;
}
int get_point_e(int pt)
{
    return utm_e+step*(pt/dimy);
    
}


char *strupr(char *s) { 
            unsigned c; 
            unsigned char *p = (unsigned char *)s; 
            while (c = *p) *p++ = toupper(c); 
			return s;
        } 



/* TO DO
while (fgets(buffer, sizeof buffer, file) != NULL) {
		name = strupr(strtok (buffer," ,.-=:"));
		strcpy(value, strtok(NULL, "\r\n"));
		if (strcmp(name, "NROWS") == 0){
		}
		if (strcmp(name, "NCOLS") == 0){
		}
		if (strcmp(name, "NBANDS") == 0){
		}
		if (strcmp(name, "NBITS") == 0){
		}
		if (strcmp(name, "BYTEORDER") == 0){
		}
		if (strcmp(name, "LAYOUT") == 0){
		}
		if (strcmp(name, "SKIPBYTES") == 0){
		}
		if (strcmp(name, "ULXMAP") == 0){
		}
		if (strcmp(name, "ULYMAP") == 0){
		}
		if (strcmp(name, "XDIM") == 0){
		}
		if (strcmp(name, "YDIM") == 0){
		}
		if (strcmp(name, "BANDROWBYTES") == 0){
		}
		if (strcmp(name, "TOTALROWBYTES") == 0){
		}
		if (strcmp(name, "BANDGAPBYTES") == 0){
		}
		if (strcmp(name, "NODATA") == 0){
		}
		if (strcmp(name, "XLLCORNER") == 0){
		}
		if (strcmp(name, "YLLCORNER") == 0){
		}
		if (strcmp(name, "CELLSIZE") == 0){
		}
		if (strcmp(name, "NO_DATAVALUE") == 0){
		}
		}
return pp;
}

*/



bool isCached(int bw, int sector)
{
    char file[100];
    sprintf(file,"%s/%d/nn_list%d/nn_list_%03d.bin",CACHED,NSECTOR,bw,sector);
    return exists(file);
}


//==================================================================================
//	Create Google Earth kml files
//==================================================================================

void savegoogleearth(bool volume){
char filename1[100];
char filename2[100];
        sprintf(filename1,"%s/%s.flt",OUTPUT,surname);
	float * surf=readFloats(filename1);
        sprintf(filename1,"%s/%s.png",OUTPUT,surname);
	data2image(filename1,surf, 10000);
	delete surf;
        sprintf(filename1,"%s/%s.kml",OUTPUT,surname);
        sprintf(filename2,"%s.png",surname);
	earth(filename1,filename2, utm_n,utm_e,dimx,step,ZONE);
	if(volume)
	{
        sprintf(filename1,"%s/%s.flt",OUTPUT,volname);
	float * volf=readFloats(filename1);
        sprintf(filename1,"%s/%s.png",OUTPUT,volname);
	data2image(filename1,volf, 10000);
	delete volf;
        sprintf(filename1,"%s/%s.kml",OUTPUT,volname);
        sprintf(filename2,"%s.png",volname);
	earth(filename1,filename2, utm_n,utm_e,dimx,step,ZONE);
	}
}



//==================================================================================
//	Read tracks and towers (single point or multiple viewshed)
//==================================================================================

//Custom function: 
//Input: Track identifier
//Output: Global array of ints and number of points
void setCustomTrack(int identifier)
{
int ntracks[17];
int *tracks[17];             // Multiple viewshed for 15 tracks in Sierra de las Nieves

	char a[100];
	int b,c,d;
	int longi[16]={0,193,315,561,711,860,1171,1493,1719,2034,2303,2617,2907,3074,3283,3350};
	int *listas[17];
	for(int i=0;i<15;i++)tracks[i]=(int *)malloc(4*(ntracks[i]=longi[i+1]-longi[i]));
        ntracks[15]=ntracks[1]+ntracks[2]+ntracks[3];
        ntracks[16]=ntracks[11]+ntracks[12];
        
        tracks[15]=(int *)malloc(4*ntracks[15]);
        tracks[16]=(int *)malloc(4*ntracks[16]);
	FILE *f1=fopen("input/tracks.txt","r");
	for(int i=0;i<15;i++){
		for(int j=0;j<ntracks[i];j++){fscanf(f1,"%d;%d;%d\n",&b,&c,&d);
		if(i!=b)
			break;
		else
			tracks[i][j]=((d-utm_e)/step)*dimy+(utm_n-c)/step;}
	}
        if(verbose)printf("Reading tracks\n");
        int cnt=0;
        for(int j=0;j<ntracks[1];j++)tracks[15][cnt++]=tracks[1][j];
        for(int j=0;j<ntracks[2];j++)tracks[15][cnt++]=tracks[2][j];
        for(int j=0;j<ntracks[3];j++)tracks[15][cnt++]=tracks[3][j];
        cnt=0;
        for(int j=0;j<ntracks[11];j++)tracks[16][cnt++]=tracks[11][j];
        for(int j=0;j<ntracks[12];j++)tracks[16][cnt++]=tracks[12][j];
        if(verbose)printf("Successfully read tracks\n");
                        towers=tracks[identifier];
                        ntowers=ntracks[identifier];

}


int t1[149]={681990,697988,715975,733963,751951,767939,785927,803915,821903,837891,855879,873867,891855,907843,925831,943819,961807,977795,995783,
        1013771,1031759,1047747,1065735,1083723,1101711,1117698,1135686,1153674,1171662,1187650,1205638,1223626,1241614,1257602,1275590,
        1293578,1311566,1327554,1345542,1363530,1381518,1397506,1415494,1433482,1451470,1467458,1485446,1503434,1521421,1537409,1555397,
        1573385,1591373,1607361,1625349,1643337,1661325,1677313,1695301,1713289,1731277,1747265,1765253,1783241,1801229,1817217,1835205,
        1853193,1871181,1841182,1813183,1783185,1755186,1725187,1697189,1667190,1639192,1611193,1581194,1553196,1523197,1495198,1465200,
        1437201,1407203,1379204,1351205,1321207,1293208,1263209,1235211,1205212,1177213,1147215,1119216,1091218,1083232,1075247,1067261,
        1059276,1053290,1045304,1037319,1029333,1021348,1015362,1007377,999391,991406,983420,977435,969449,961464,953478,945493,939507,
        931522,923536,915551,907565,901580,893594,885609,877623,869638,863652,855667,847681,839695,831710,825724,817739,809753,801768,
        793782,787797,779811,771826,763840,755855,749869,741884,733898,725913,717927,711942,703956,695971,687985};

int t2[162]={681990,695987,711975,727962,743950,759937,775925,791912,807900,821888,837875,853863,869850,885838,901825,917813,933800,947788,
963775,979763,995751,1011738,1027726,1043713,1059701,1073688,1089676,1105663,1121651,1137639,1153626,1169614,1185601,1199589,1215576,1231564,
1247551,1263539,1279526,1295514,1311502,1311487,1311472,1313457,1313442,1315427,1315412,1315397,1317382,1317367,1319353,1343345,1369337,1393329,
1419321,1443313,1469306,1493298,1519290,1543282,1569274,1595267,1619259,1645251,1669243,1695235,1719227,1745220,1769212,1795204,1819196,1845188,
1871181,1841182,1813183,1783185,1755186,1725187,1697189,1667190,1639192,1611193,1581194,1553196,1523197,1495198,1465200,1437201,1407203,1379204,
1351205,1321207,1293208,1263209,1235211,1205212,1177213,1147215,1119216,1091218,1065225,1041232,1015239,991247,965254,941261,915268,891275,
865283,841290,815297,791304,767312,753325,739338,725351,711365,697378,685391,671404,657417,643431,629444,617457,603470,589484,575497,561510,
549523,535537,521550,507563,493576,481590,487604,493618,501632,507646,515660,521674,529689,535703,543717,549731,555745,563759,569773,577788,
583802,591816,597830,605844,611858,617872,625886,631901,639915,645929,653943,659957,667971,673985};

//algoritmo montes de malaga
int t3[235]={
681990,707995,735990,763985,791981,819976,847971,873967,901962,929957,
957952,985948,1013943,1039938,1067933,1095929,1123924,1151919,1179915,1193903,
1207892,1221881,1235869,1249858,1263847,1277836,1249830,1223824,1197818,1169813,
1143807,1117801,1091796,1063790,1037784,1011778,983773,957767,931761,905756,
915742,927729,937715,949702,961689,969676,979663,989650,997637,1007624,
1017611,1039602,1061594,1083585,1107577,1129569,1151560,1175552,1197544,1219535,
1243527,1265518,1287510,1311502,1311487,1311472,1313457,1313442,1315427,1315412,
1315397,1317382,1317367,1319353,1341345,1363336,1385328,1409320,1431312,1453304,
1475296,1499288,1521280,1543272,1565264,1589256,1613249,1639242,1665235,1691228,
1717221,1741215,1767208,1793201,1819194,1845187,1871181,1865167,1859152,1855138,
1849124,1845110,1839096,1835082,1829068,1825054,1819040,1815026,1809012,1804998,
1776993,1750989,1722985,1696981,1668977,1642973,1614969,1588965,1576951,1564937,
1552924,1540910,1530896,1518883,1506869,1494855,1484842,1460851,1438859,1414868,
1392877,1370886,1346895,1324904,1302913,1278922,1256931,1234940,1210949,1188958,
1166967,1142976,1120985,1098994,1075003,1053012,1031021,1007029,985037,963046,
939054,917063,895071,871080,849088,827096,803105,781113,759122,783129,
809136,835144,861151,885159,911166,937173,963181,987188,1013195,1039203,
1065210,1091218,1071229,1053240,1033251,1015263,995274,977285,957296,939308,
909308,881309,853310,823310,795311,767312,753325,739338,725351,711365,
697378,685391,671404,657417,643431,629444,617457,603470,589484,575497,
561510,549523,535537,521550,507563,493576,481590,487604,493618,501632,
507646,515660,521674,529689,535703,543717,549731,555745,563759,569773,
577788,583802,591816,597830,605844,611858,617872,625886,631901,639915,
645929,653943,659957,667971,673985};

void setCustomTrack_a(int identifier)
{
    printf("estoy aqui\n");
    switch(identifier)
    {
        case 0:
            towers=t1;
            ntowers=149;
        return;
            break;
        case 1:
            towers=t2;
            ntowers=162;
        return;
            break;
        case 2:
            towers=t3;
            ntowers=235;
            
        return;
    }

}




int setCustomTower(int id)
{
    switch(id)
    {
        case 0:
        //tower=get_point_id(4063790,324410); //Caina
                        return get_point_id(4062820,324430); //Carnicer
                        break;
                        //return;
        case 1:
                        
                        return get_point_id(4059940,313470); //conejeras
                        break;
                        //return;
        case 2:
                        
                        //return get_point_id(4063890,324530); //Caina
                        return get_point_id(4063880,324460); //Caina
                        break;
		default:
						return 0;
						break;
    }
}


//==================================================================================
//	Get towers (if coverstore)
//==================================================================================
void setTowers()
{
    char stmp[100];
    if(towernotgiven)
	try{
		if(verbose)printf("Computing max value from existing output: %s/%s.flt\n",OUTPUT,surname);
                sprintf(stmp,"%s/%s.flt",OUTPUT,surname);
		tower=maxSurf(stmp);				
		if(verbose&&(tower!=0))printf("Tower suggested at %d (%d,%d)\n",tower,get_point_n(tower),get_point_e(tower));
		if(verbose&&(tower==0)){printf("No tower given and no candidate found. Using a tower at point 0\n");}
	}catch(int E) { }
        ntowers=1;
        towers=(int *)malloc(4);
        towers[0]=tower;
        if(verbose)printf("Computing Viewshed for tower %s at point %d\n",idname,tower);
}
void setTracks()
{
		//if(!volume)readTowers();
		//if(volume)setCustomTrack(identifier);
                setCustomTrack_a(identifier);
                        
                cnttower=0;
                tower=towers[cnttower++];
                if(verbose)printf("Viewshed for track %s with %d pt. First tower at %d\n",idname,ntowers, towers[0]);
}




//==================================================================================
//	Intro
//==================================================================================

header parse_hdr(char *filename){
	header hdr;
char tmpfn[100];
char hdrfile[100];
hdr.nitems=0;

strcpy(tmpfn,filename);
char *prefix = strtok(tmpfn,".");
strcpy(hdrfile,prefix);
strcat(hdrfile,".hdr");
if(!exists(hdrfile)){
    printf("Missing header file %s. Using default cell size %d and UTM coordinates .\n",hdrfile,STEP);
    return hdr;
    }
FILE *file=fopen(hdrfile,"r");
char buffer[200];
char *name,*value;
char names[200];
//char value[200];
while (fgets(buffer, sizeof buffer, file) != NULL) {
		name = strtok (buffer," ,.-=:");
                strupr(name);
		value = strtok(NULL, "\r\n");
                
                char *tmp;
                tmp=(char *)malloc(strlen(name)+1);
                strcpy(tmp,name);
                hdr.keys[hdr.nitems]=tmp;
                tmp=(char *)malloc(strlen(value)+1);
                strcpy(tmp,value);
                hdr.values[hdr.nitems]=tmp;
                hdr.nitems++;
                //printf("%s %s\n",name,value);
                }
return hdr;
}


bool setfilenames(bool fullstore)
{
    char tmp[100];
    sprintf(tmp,"%s/%d/nn_list%d",CACHED,NSECTOR,bw);
    if(!exists(tmp))create_directory(tmp);
    sprintf(tmp,"%s/full",OUTPUT);
    if(fullstore)if(!exists(tmp))create_directory(tmp);
    
    sprintf(cvrgfname,"%s_cvrg_%s",base_name,idname);
    sprintf(surcname,"%s_surc_%s",base_name,idname);
    sprintf(volcname,"%s_volc_%s",base_name,idname);
    sprintf(horzname,"%s_horz",base_name);
    
    return true;
}






bool checkfile(char *argv)
{
header hdr;
char cadena[100];
char *extension;
char *basen;
long size;
strcpy(base_name,argv);
strtok(base_name,".");
//strcpy(cadena,INPUTDIR);
if(!exists(INPUTDIR))create_directory(INPUTDIR);
if(!exists(OUTPUT))create_directory(OUTPUT);
if(!exists(CACHED))create_directory(CACHED);
sprintf(cadena,"%s/%d",CACHED,NSECTOR);
if(!exists(cadena))create_directory(cadena);
sprintf(surname,"%s_sur",base_name);
sprintf(surlname,"%s_surlog",base_name);
sprintf(volname,"%s_vol",base_name);
sprintf(vollname,"%s_vollog",base_name);


if(!silent)printf(filename,"%s/%s\n",INPUTDIR,argv);
sprintf(filename,"%s/%s",INPUTDIR,argv);
hdr.nitems=0;
int i=sscanf(argv,"%d_%d_%d.bil",&utm_n,&utm_e,&step);
int hdimx=0,hdimy=0,hstep=10,hutmn=0,hutme=0;
if(i!=3){ //Este archivo no es nuestro
    hdr=parse_hdr(filename);
    if(hdr.nitems==0){
        for(int i=0;i<hdr.nitems;i++){
            if(!strcmp(hdr.keys[i],"NROWS")){sscanf(hdr.values[i],"%d",&hdimy);}
            if(!strcmp(hdr.keys[i],"NCOLS")){sscanf(hdr.values[i],"%d",&hdimx);}
            if(!strcmp(hdr.keys[i],"CELLSIZE")){sscanf(hdr.values[i],"%d",&hstep);}
            if(!strcmp(hdr.keys[i],"XLLCORNER")){sscanf(hdr.values[i],"%d",&hutme);}
            if(!strcmp(hdr.keys[i],"YLLCORNER")){sscanf(hdr.values[i],"%d",&hutmn);}
        }
        utm_e=hutme;
        utm_n=hutmn+hdimy*hstep;
        step=hstep;
        printf("%d %d %d\n",utm_e,utm_n,step);
        }
    
}
if((extension = strrchr(filename,'.')) != NULL )floatsDEM=!strcmp(extension,".flt");
//strcpy(filename,cadena);
size=filesize(filename);
if(size>0){
        dimx=dimy=(int)sqrt((double)(size/(floatsDEM?4:2)));
}else return false;
if((hdr.nitems!=0)&&((dimx!=hdimx)||(dimy!=hdimy)))return false;
dim=dimx*dimy;
bw=1+((dimx+dimy)/4);

return true;
}


void checkargs(int argc, char **argv){
 char cadena[100];
 bool ayuda=false;
 towernotgiven=true;
 if(IS_BIG_ENDIAN){
     printf("Big Endian architectures not supported now\n");
     exit(0); //Comment this line if you can test with your own big endian models
     }
 endianOK=true;
 
 
 
 
 if(argc==1)ayuda=true;
 else ayuda=!checkfile(argv[1]);
 if(argc==2)return;
 int j;
 for(int i=2;i<argc;i++)
 {
     if(argv[i][0]!='-'){ayuda=true;break;}
     switch(argv[i][1])
     {
         case 'H':
             j=sscanf(argv[i],"-H%lf",&obsheight);
             if(j==0)obsheight=OBS_H;
             if(verbose)printf("Height: %f\n",(float)obsheight);
             break;
         case 'v':
             verbose=true;
             break;
         case 'q':
             silent=true;
             break;
         case 'F':
             fullstore=true;
             break;
         case 'm':
             inputmask=true;
             j=sscanf(argv[i],"-m%s",maskname);
             if(strlen(maskname)==0)sprintf(maskname,"%s/mask.dat",INPUTDIR);
             break;
         case 't':
             j=sscanf(argv[i],"-t%d",&mode);
             printf("Modo: %d\n",mode);
             if(j==0)mode=0;
             break;
         case '?': 
             break;
         case 'n': 
             sscanf(argv[i],"-n%s",cadena);
             int utmn,utme;
             try{
             j=sscanf(cadena,"%d,%d,%s",&utmn,&utme,idname);
             if(j==3){
                 tower=get_point_id(utmn,utme);
                 if(tower<dim)towernotgiven=false;}
             else{
                 j=sscanf(cadena,"%d,%s",&identifier,idname);    
                 if(j==2){
                    tower=setCustomTower(identifier);
                    if(tower<dim)towernotgiven=false;}
                 else{
                    j=sscanf(cadena,"%d",&identifier);;
                    if(j==1){
                        strcpy(idname,cadena);
                        tower=setCustomTower(identifier);
                        if(tower<dim)towernotgiven=false;
                        //printf("%d %d\n",tower,dim);
                        }
                    }
                }
             }catch(int E){
                towernotgiven=true;
             }
             break;
         default:
             ayuda=true;
             break;
     }
 }
    if(ayuda){printf("Usage: %s model [-v] [-q] [-ttype][-Hheight][-mmaskfilename] [-F] [-nitem]\n\n"
            "-v Verbose mode\n"
            "-q Silent mode\n"
            "-tN Execution mode:\n"
            "\t N=0 Total viewshed (default)\n"
            "\t N=1 Total 3D-Viewshed\n"
            "\t N=2 Single tower viewshed\n"
            "\t N=3 Track Viewshed\n"
            "\t N=4 Single tower 3D-viewshed\n"
            "\t N=5 Track 3D-Viewshed\n"
            "\t N=6 Find sequential towers (cell coverage algorithm)\n"
            "\t N=7 Isolation index\n"
            "\t N=8 Horizons\n"
            "\t N=10 Modo quien me ve\n"
            "\t N=9 Timing mode\n"
            "-Hheight Observer's height in meters (default, 1.5m)\n"
            "-m Use default mask file (input/mask.dat for area of interest\n"
            "-mfile Use file for area of interest\n"
			"-F Fullstore flag: Save the geometry of total viewshed (types 1 and 2)"
            "-nitem In single-tower and several-tower (tracks) viewshed execution modes, specifies the track:\n"
            "\t item:utmn,utme,name\n"
            "\t item:index,name\n"
            "\t item:index\n"
            "\t Use custom function to assign indexes to towers\n\n"
            "",argv[0]);exit(0);}
    //exit(0);
}
    
