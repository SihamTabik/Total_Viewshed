/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   sector2.h
 * Author: felipe
 *
 * Created on 19 de abril de 2016, 8:49
 */



#ifndef SECTOR2_H
#define SECTOR2_H

class Sector2 {
public:
    Sector2();
    Sector2(const Sector2& orig);
    virtual ~Sector2();
    
    
    
	double *heights;
        float open_delta_d;
        float delta_d;
        float open_delta_h;
        float delta_h;
        bool visible;
      	float max_angle;
        
    
    
private:
    inline void Sector2::kernelS(int rs[][2], int &nrs, float &sur_data);

};

#endif /* SECTOR2_H */

