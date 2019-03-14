#include <stdio.h>
#include <iostream>
#include <cstdlib>
#include <cmath>

using namespace std;

class RAM_Access{

    public:

    int x;
    int y;
    int address;

    RAM_Access(){}
    RAM_Access(int a, int b){
        x = a;
        y = b;
    }

    RAM_Access(const RAM_Access &old){
        x = old.x;
        y = old.y;
    }
//    bool operator == (RAM_Access a){
//
//        if(x == a.x && y == a.y){
//            return true;
//        }
//
//        return false;
//    }
};

int RAM_SIZE, mode;

long long SRAM_access = 0, DRAM_access = 0;

long long ram_load = 0;
long ten2eight = 100000000;
long DRAM_access_m = 0;
/*
 *  VR Projection
 */
double PI = 3.14159265358979323846;

int w ,h;

double toRadian(double a){
    return a / 180.0 * PI;
}

int nearestNeighbor(double num){

  int res = (int)(num + 0.5);

  return res;
}

void spherical2cartesian(double the, double phi, double result[3]){

    double x = sin(phi) * cos(the);
    double y = sin(phi) * sin(the);
    double z = cos(phi);


    result[0] = x;
    result[1] = y;
    result[2] = z;

}

void spherical2coordinates(double the, double phi, double result[2]){

    double i,j;

    if(the > PI){
        i = (the - PI) / 2.0 / PI * w;
    }
    else{
        i = (PI + the) / 2.0 / PI * w;
    }

    j = phi /  PI * h;

    result [0] = i;
    result [1] = j;
}

void cartesian2coordinates(double x, double y, double z, double result[2]){

    double the;

    if(x != 0) {
        the = atan2(y, x);

    } else {
        the = toRadian(90.0);
    }

    double phi = acos(z);
    spherical2coordinates(the, phi, result);
}

void matrixMultiplication(double* vector, double matrix[3][3], double res[3]) {

    for (int i=0; i<3; i++) {
        for (int j=0; j<3; j++) {

            res[i] += matrix[i][j] * vector[j];

        }
    }

}

void loadSRAM(RAM_Access ***dram, RAM_Access **sram, int i, int mode){

    int m = i, n = 0;

    switch (mode) {

        case 0:
            for (int s = 0; s < RAM_SIZE; s++) {
                sram[s] = dram[m][n];
                n++;
//                printf("sram i: %d, j:%d\n", sram[s]->x, sram[s]->y);
                if (m > w - 1)
                    break;

                if (n == h) {
                    m++;
                    n = 0;
                }

                DRAM_access++;
                if(DRAM_access == ten2eight){
                    DRAM_access_m++;
                    DRAM_access = 0;
                }
            }
            break;
    }

}

void initDRAM(RAM_Access*** dram){

    dram = new RAM_Access **[h];

    for(int i = 0; i < w; i++){

        dram[i] = new RAM_Access*[h];

        for(int j = 0; j < h; j++){

            dram[i][j] = new RAM_Access(i, j);
//            printf("i: %d, j:%d\n", dram[i][j]->x, dram[i][j]->y);
        }
    }
}

void initSRAM(RAM_Access** sram){

    sram = new RAM_Access *[RAM_SIZE];

    for(int i = 0; i < RAM_SIZE; i++){

        sram[i] = new RAM_Access(-1, -1);
        //printf("i: %d, j:%d\n", sram[i]->x, sram[i]->y);
    }
}

bool sram_contains(int i, int j, RAM_Access** sram){

    for(int s = 0; s < RAM_SIZE; s++){

        if(sram[s]->x == i && sram[s]->y == j){

            return  true;
        }
    }

    return false;
}

bool check_sram(int i, int j, RAM_Access** sram){

    return sram_contains(i, j, sram);
}

int main(int argc, char** argv) {

    // number of pixels sram could store
    // 627 is in KB
    RAM_SIZE = 627 * 1024 / 3;

    // get width and height
    w = 3840;
    h = 2160;
    int fw = 1174, fh = 1080;
    int fovX = 110, fovY = 90;
    double hp = 45,ht = 45;

    RAM_Access **SRAM = (RAM_Access**)calloc(RAM_SIZE, sizeof(RAM_Access));

    SRAM = new RAM_Access *[RAM_SIZE];

    for(int i = 0; i < RAM_SIZE; i++){

        SRAM[i] = new RAM_Access(-1, -1);
        //printf("i: %d, j:%d\n", sram[i]->x, sram[i]->y);
    }

    RAM_Access ***DRAM = (RAM_Access***)calloc(w, sizeof(RAM_Access**));

    DRAM = new RAM_Access **[w];

    for(int i = 0; i < w; i++){

        DRAM[i] = new RAM_Access*[h];

        for(int j = 0; j < h; j++){

            DRAM[i][j] = new RAM_Access(i, j);
//            printf("i: %d, j:%d\n", dram[i][j]->x, dram[i][j]->y);
        }
    }
//
//    initDRAM(DRAM);
//    initSRAM(SRAM);

    // convert to radian
    double htr = toRadian(ht);
    double hpr = toRadian(hp);

    // rotation matrices
    double rot_y [3][3] = {
            {cos(hpr), 0, -sin(hpr)},
            {0, 1, 0},
            {sin(hpr), 0, cos(hpr)}
    };

    double rot_z [3][3] = {
            {cos(htr), sin(htr), 0},
            {-sin(htr), cos(htr), 0},
            {0, 0, 1}
    };

    int a = 0, b = 0;

    int pixel_count = 0;
    // default head orientation is 0,90
    for (double i = 90  - fovY/2.0; i < 90 + fovY/2.0; i+= fovY*1.0/fh, b++) {
        for (double j = -fovX/2.0; j < fovX/2.0; j+= fovX*1.0/fw, a++) {

            double p1[] = {0.0, 0.0, 0.0};
            spherical2cartesian(toRadian((j < 0) ? j + 360 : j), toRadian((i < 0) ? (i + 180) : i), p1);
            // rotation along y axis
            double p2[] = {0.0, 0.0, 0.0};
            matrixMultiplication(p1, rot_y, p2);

            // rotate along z axis
            double p3[] = {0.0, 0.0, 0.0};
            matrixMultiplication(p2, rot_z, p3);

            double res[] = {0.0,0.0};

            // convert 3D catesian to 2D coordinates
            cartesian2coordinates(p3[0], p3[1], p3[2], res);

            if (b >= fh || a >= fw){
                break;
            }

            int temp_x = nearestNeighbor (res[0]);
            int temp_y = nearestNeighbor (res[1]);

            if(!check_sram(temp_x, temp_y, SRAM)){

                loadSRAM(DRAM, SRAM, temp_x, mode);
                ram_load++;
            }
            SRAM_access++;

//            printf("Pixel: %d\n", pixel_count++);
        }
        a = 0;
    }

    printf("DRAM Access: %dE+08 %d\n", DRAM_access_m, DRAM_access);
    printf("SRAM Access: %d\n", SRAM_access);
    printf("SRAM Load: %d\n", ram_load);

    return 0;
}
