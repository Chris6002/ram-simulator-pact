#include <stdio.h>
#include <iostream>
#include <cstdlib>
#include <cmath>

using namespace std;

class Address{

    public:

    int x;
    int y;
    int address;

    Address(){}
    Address(int a, int b) : x(a), y(b) {}
    bool operator == (Address a){

        if(x == a.x && y == a.y){
            return true;
        }

        return false;
    }
};

double PI = 3.14159265358979323846;

int w ,h;

int RAM_SIZE;

int SRAM_access = 0, DRAM_access = 0;

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

void loadSRAM(Address *sram, int size, Address src){



}

void initDRAM(Address*** dram){

    dram = new Address **[h];

    for(int i = 0; i < w; i++){

        dram[i] = new Address*[h];
        
        for(int j = 0; j < h; j++){

            dram[i][j] = new Address(i, j);
        }
    }
}

int main(int argc, char** argv) {

    // get width and height
    w = 3840;
    h = 2160;
    RAM_SIZE = 627 * 1024 / 3;
    int fw = 1174, fh = 1080;
    int fovX = 110, fovY = 90;
    double hp = 45,ht = 45;

    Address SRAM[RAM_SIZE];
    Address ***DRAM;

    initDRAM(DRAM);
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
            // assign the pixel value
           // fov.at<Vec3b>(b,a) = image.at<Vec3b>(nearestNeighbor (res[1]), nearestNeighbor (res[0]));
            
        }
        a = 0;
    }

    return 0;
}
