#define MAX(x, y) (((x) > (y)) ? (x) : (y))


        int cls_range;

        //Index of closest points, for each angular range, for each distance level, for left and right delimiters, for row and col

        //
		/*
        int cls_shift[8][2][2][2]= {
        {{ { 0, 1 }, { 1, 1 } } , { { 0, 2 }, { 1, 2 } } } ,
        {{ { 0, 1 }, { 1, 1 } } , { { 1, 2 }, { 2, 2 } } } ,
        {{ { 1, 1 }, { 1, 0 } } , { { 2, 2 }, { 2, 1 } } } ,
        {{ { 1, 1 }, { 1, 0 } } , { { 2, 1 }, { 2, 0 } } } ,
        {{ { 1, 0 }, { 1,-1 } } , { { 2, 0 }, { 2,-1 } } } ,
        {{ { 1, 0 }, { 1,-1 } } , { { 2,-1 }, { 2,-2 } } } ,
        {{ { 1,-1 }, { 0,-1 } } , { { 2,-2 }, { 1,-2 } } } ,
        {{ { 1,-1 }, { 0,-1 } } , { { 1,-2 }, { 0,-2 } } }
        };
		*/


        //Coeficients for left point. For each range, and each level
        double cls_coefs[2];
        double cls_dist;
        float cls_mindist;


        void closest_parameters(double angulo)
        {
            double torad = PI / 180;
            double tn = tan(angulo*torad);
            double co = cos(angulo * torad);
            double se = sin(angulo * torad);
            double ct = 1 / tn;
               cls_range = 7;
            if (angulo < 153.44) cls_range = 6;
            if (angulo < 135.00) cls_range = 5;
            if (angulo < 116.56) cls_range = 4;
            if (angulo < 90.00)  cls_range = 3;
            if (angulo < 63.44)  cls_range = 2;
            if (angulo < 45.00)  cls_range = 1;
            if (angulo < 26.56)  cls_range = 0;

            if ((cls_range == 0) || (cls_range == 1)) cls_dist = 1 / co;
            if ((cls_range == 2) || (cls_range == 3)) cls_dist = 1 / se;
            if ((cls_range == 4) || (cls_range == 5)) cls_dist = 1 / se;
            if ((cls_range == 6) || (cls_range == 7)) cls_dist = -1 / co;

			//In step units
            //cls_dist *= step;
            cls_mindist = (float) (cls_mindistance[cls_range] );

             if (cls_range == 0){cls_coefs[0]=tn;cls_coefs[1]=2 * tn ;};
             if (cls_range == 1){cls_coefs[0]=tn;cls_coefs[1]=2 * tn - 1; };
             if (cls_range == 2){cls_coefs[0]=1-ct;cls_coefs[1]=2 - 2* ct; };
             if (cls_range == 3){cls_coefs[0]=1-ct;cls_coefs[1]=1 -2 * ct ;};
             if (cls_range == 4){cls_coefs[0]=-ct;cls_coefs[1]=- 2 * ct ;};
             if (cls_range == 5){cls_coefs[0]=-ct;cls_coefs[1]=- 2 * ct -1; };
             if (cls_range == 6){cls_coefs[0]=1+tn;cls_coefs[1]= 2 + 2* tn ;};
             if (cls_range == 7) { cls_coefs[0] = 1 + tn; cls_coefs[1] = 1 + 2 * tn; };



        }

bool *cls_flag;
float *cls_angle;
float cls_2angles[2];
bool  cls_2flags[2];

void get_closest(double a, int ydim, double* h, double obs_h,int row, int col)
        {
            double cls_l0, cls_l1, cls_l0o, cls_l1o;
            double k0 = cls_coefs[0];
            double k1=  cls_coefs[1];
            cls_l0 =
                h[(row + cls_shift[cls_range*8+0*4+0*2+0]) + ydim * (col + cls_shift[cls_range*8+0*4+0*2+1])] * (1 - k0) +
                h[(row + cls_shift[cls_range*8+0*4+1*2+0]) + ydim * (col + cls_shift[cls_range*8+0*4+1*2+1])] * (k0);
            cls_l1 =
                h[(row + cls_shift[cls_range*8+1*4+0*2+0]) + ydim * (col + cls_shift[cls_range*8+1*4+0*2+1])] * (1 - k1) +
                h[(row + cls_shift[cls_range*8+1*4+1*2+0]) + ydim * (col + cls_shift[cls_range*8+1*4+1*2+1])] * (k1);
            cls_l0o =
                h[(row - cls_shift[cls_range*8+0*4+0*2+0]) + ydim * (col - cls_shift[cls_range*8+0*4+0*2+1])] * (1 - k0) +
                h[(row - cls_shift[cls_range*8+0*4+1*2+0]) + ydim * (col - cls_shift[cls_range*8+0*4+1*2+1])] * (k0);
            cls_l1o =
                h[(row - cls_shift[cls_range*8+1*4+0*2+0]) + ydim * (col - cls_shift[cls_range*8+1*4+0*2+1])] * (1 - k1) +
                h[(row - cls_shift[cls_range*8+1*4+1*2+0]) + ydim * (col - cls_shift[cls_range*8+1*4+1*2+1])] * (k1);
            double d = cls_dist;
            double povh=obs_h+h[row + ydim* col];
            double f0 = (cls_l0 -povh ) / d;
            double f1 = (cls_l1 - povh) / (2*d);
            double b0 = (cls_l0o - povh) / d;
            double b1 = (cls_l1o - povh) / (2 * d);
            cls_2flags[0] = f1 > f0;
            cls_2flags[1] = b1 > b0;
	    cls_2angles[0]=(float)MAX(f0, f1);
            cls_2angles[1]=(float)MAX(b0, b1);
        }


        void cls_borders(int xdim, int ydim)
        {
            for (int j = 0; j < xdim; j++)
                for (int i = 0; i < 2; i++)
                {
                    cls_angle[2*(j * ydim + i)+0] = -500;
                    cls_angle[2*(j * ydim + i)+1] = -500;
                     cls_flag[2*(j * ydim + i)+0] = true;
                     cls_flag[2*(j * ydim + i)+1] = true;
                    cls_angle[2*(j * ydim + ydim - 2 + i)+0] = -500;
                    cls_angle[2*(j * ydim + ydim - 2 + i)+1] = -500;
                     cls_flag[2*(j * ydim + ydim - 2 + i)+0] = true;
                     cls_flag[2*(j * ydim + ydim - 2 + i)+1] = true;
                }
            for (int j = 0; j < 2; j++)
                for (int i = 2; i < ydim-2; i++)
                {
                    cls_angle[2*(j * ydim + i)+0] = -500;
                    cls_angle[2*(j * ydim + i)+1] = -500;
                     cls_flag[2*(j * ydim + i)+0] =true;
                     cls_flag[2*(j * ydim + i)+1] =true;
                    cls_angle[2*((j + xdim - 2) * ydim + i)+0] = -500;
                    cls_angle[2*((j + xdim - 2) * ydim + i)+1] = -500;
                     cls_flag[2*((j + xdim - 2) * ydim + i)+0] = true;
                     cls_flag[2*((j + xdim - 2) * ydim + i)+1] = true;
                }
        }


		void cls_compute(double ang, double* height,double obs_h, int xdim, int ydim)
        {
            cls_borders(xdim, ydim);
            for(int j=2;j < xdim-2; j++)
                for (int i = 2; i < ydim - 2; i++)
                {
		   get_closest(ang, xdim, height,obs_h, i, j);
                   cls_angle[2*(j * ydim + i)] = cls_2angles[0];
                   cls_angle[2*(j * ydim + i)+1] = cls_2angles[1];
                   cls_flag[2*(j * ydim + i)] = cls_2flags[0];
                   cls_flag[2*(j * ydim + i)+1] =cls_2flags[1];
                }
        }


		void get_cls_params(int pt, int B, bool &visible, float &mangle, float &mindist)
        {
            visible = cls_flag[2*pt+B];
            mangle = cls_angle[2*pt+B];
            mindist = (float)cls_mindist;  //Only points beyond this distance
        }