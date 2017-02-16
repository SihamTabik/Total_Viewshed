/*
AUTHORS: 
*SIHAM TABIK
*LUIS FELIPE ROMERO GOMEZ
*ANTONIO MANUEL RODRIGUEZ CERVILLA 
DEPT: ARQUITECTURA DE LOS COMPUTADORES
UNIVERSIDAD DE MALAGA
 
*/


// Definition of structures and variables used.
#define _USE_MATH_DEFINES
#include <math.h>
#include <string>
#ifndef DEFS_H
#define DEFS_H

//Actually overrided by inputfilename
//#define SCALE 10.0
#define STEP 10
#define UTMN 4080000 
#define UTME 360000
#define CACHED "../cached"
#define OUTPUT "./output"
#define INPUTDIR "./input"
#define _USE_CLOSEST_POINTS


// UTM zone
#define ZONE "30S"


// obs_h=10;       observer's height in meters
#define OBS_H 1.5


// Numero maximo de puntos del vuelo que ven un punto del territorio
#define NUMNC 5

  // openmp default number of threads
#ifndef WIN32
#define NTHREAD 12
#else
#define NTHREAD 4
#endif



  // -Single Viewshed- Default location of a tower in north-to-south first ordering 
//#define TOWERLOC 2013016



	// Initial sector shift in degrees, to avoid point alignment
	// Bandwith and half bandwith of data structure
	// Order of sqrt(dim)
	
#define SHIFT 0.001
//052345245534643563456

#define NSECTOR 180
//#define BW 1251

//T fila es i%dimy
//T colu es i/dimy
#define T2N(i) ((dimy*((i)%dimy))+((i)/dimy))

//N fila es i/dimx
//N colu es i%dimx
#define N2T(i) ((dimx*((i)%dimx))+((i)/dimx))

using namespace std;
		
const double PI = M_PI;				
const double mpi= M_PI/2;			
const double tmpi= 3*M_PI/2;
const double tograd = 360/(2*M_PI);	// conver const of radians to degree	
const double torad =(2*M_PI)/360;		// conver const of degree to radians
const float  PI_F=3.14159265358979f;
const int isolationindex=2; //0 max 1 m_arit 2 m_harm 3 m_geom

	//Struct node
	struct node {	//	This structure has the information of one point
        int idx;	//	identifier of poit.
        int oa;		//	orden given the point in our algorithm.
        float h;	//	height at that id it.
        float d;	//	distance the point from axis x.
    	};

      	struct nodoscolores {	//
        int idx;	//	identifier of poit.
        float d;	//	distance the point from axis x.
    	};

	
	// Node used in the Linkded List
	struct LinkedListNode{ //It is a node whit next an previous identifier
		public:		
		node Value;
        short next;// Identifier the next node in Linked List.
        short prev;// Identifier the Previous node in Linked List.
	};

        
	struct sectorgeometry {	//
        int dimx;	//
        int dimy;		//
        int bw;
        double obs_h;	//
        double step;	//
    	};        
	
	typedef  struct header{
		int nitems;
		char * keys[50];
		char * values[50];
	}header;

#ifdef USE_CLOSEST_POINTS

	const int cls_shift[]={
			0,1,1,1,0,2,1,2, 
			0,1,1,1,1,2,2,2, 
			1,1,1,0,2,2,2,1, 
			1,1,1,0,2,1,2,0, 
			1,0,1,-1,2,0,2,-1, 
			1,0,1,-1,2,-1,2,-2, 
			1,-1,0,-1,2,-2,1,-2,
			1,-1,0,-1,1,-2,0,-2};

        const double cls_mindistance[]=
        { sqrt(5.0), sqrt(8.0), sqrt(8.0),sqrt(5.0),sqrt(5.0), sqrt(8.0), sqrt(8.0),sqrt(5.0) }; //in steps

#endif


#endif
