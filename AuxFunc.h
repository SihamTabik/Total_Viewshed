/*
AUTHORS: 
*SIHAM TABIK
*LUIS FELIPE ROMERO GOMEZ
*ANTONIO MANUEL RODRIGUEZ CERVILLA 
DEPT: ARQUITECTURA DE LOS COMPUTADORES
UNIVERSIDAD DE MALAGA
 
*/

// auxiliary functions


#include <string>
#include <time.h>
#ifdef WIN32
#include <winsock2.h>
#include <windows.h>
#else
#include <sys/time.h>
//#include <boost/filesystem.hpp>
#include <sys/types.h>
#include <unistd.h>
#endif



#ifndef AUXFUN_H
#define AUXFUN_H

using namespace std;


struct color
{
     unsigned char R;
     unsigned char G;
     unsigned char B;
};



string converInt(int number); // Convert integer to string
double dtime();
bool readCover(bool *data,   char  *fname);
void resetCover(bool* c);
int saveCover(bool *data,   char  *fname, bool append);
void resetVolCover(float* c);
double *readHeights(char *x, bool floats, bool sameendian);
double *readHeights(char *x);
int maxSurf(char * x);
color *build_palette();
void data2image(const char* filename, float *data ,float maxval=10000);
float * readFloats(char *nmds);
void saveFloats( float * d, char  * filename, bool transpose);
void saveHorizons( unsigned char * d, char  * filename);
void readTowers();
void writeheader(int n, int e, int step, char *filename,bool isbool);
int get_point_id(int n, int e);
bool isCached(int bw, int sector);
void earth( char *KFN,  char *JFN, int utm_n, int utm_e, int Dim, int step,  const char* zone );
void savegoogleearth(bool volume);
bool alineados(int a, int b, int c);
bool setfilenames(bool fullstore);
void checkargs(int argc, char **argv);


void setCustomTrack(int i);
void setTowers();
void setTracks();
void saveColores(nodoscolores * sc);

#endif
