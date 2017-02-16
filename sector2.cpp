/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   sector2.cpp
 * Author: felipe
 * 
 * Created on 19 de abril de 2016, 8:49
 */
#include <stdio.h>
#include <cmath>
#include <omp.h>
#include "Sector2.h"

Sector2::Sector2() {
}

Sector2::Sector2(const Sector2& orig) {
}

Sector2::~Sector2() {
}

inline void Sector2::kernelS(int rs[][2], int &nrs, float &sur_data)
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
