#include <stdio.h>
#include <iostream>
#include <cstdlib>
#include <sys/types.h>
#include <errno.h>
#include <unistd.h>
#include <sys/fcntl.h>
#include <time.h>
#include <opencv2/opencv.hpp>
#include <thread>
#include <atomic>

using namespace cv;

double PI = 3.14159265358979323846;

int w ,h;

double tileSize;
double tileSizeX,tileSizeY;


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


void findPixel(int index, double x,double y, double result[2]) {

    int vertical;

    if (index > 2) {
        vertical = 1;
    }
    else {
        vertical = 0;
    }

    double n = (tileSize * (index % 3))  + y * tileSize;
    double m = (tileSize * vertical) + x * tileSize;

    result [0] = n;
    result [1] = m;

}


//from wikipedia: https://en.wikipedia.org/wiki/Cube_mapping
void convert_xyz_to_cube_uv(double x, double y, double z, double result[2]) {

    double maxAxis, uc, vc;
    double u, v;
    int index;

    double absX = fabs(x);
    double absY = fabs(y);
    double absZ = fabs(z);

    int isXPositive = x > 0 ? 1 : 0;
    int isYPositive = y > 0 ? 1 : 0;
    int isZPositive = z > 0 ? 1 : 0;


    // POSITIVE X
    // Front
    if (isXPositive && absX >= absY && absX >= absZ) {
        // u (0 to 1) goes from +z to -z
        // v (0 to 1) goes from -y to +y
        index = 4;
        maxAxis = absX;
        uc = -z;
        vc = y;

    }

    // NEGATIVE X
    // Back
    else if (!isXPositive && absX >= absY && absX >= absZ) {
        // u (0 to 1) goes from -z to +z
        // v (0 to 1) goes from -y to +y
        index = 5;
        maxAxis = absX;
        uc = -z;
        vc = -y;

    }

    // POSITIVE Y
    // Left
    else if (isYPositive && absY >= absX && absY >= absZ) {
        // u (0 to 1) goes from -x to +x
        // v (0 to 1) goes from +z to -z
        index = 1;
        maxAxis = absY;
        uc = -z;
        vc = -x;

    }

    // NEGATIVE Y
    // Right
    else  if (!isYPositive && absY >= absX && absY >= absZ) {
        // u (0 to 1) goes from -x to +x
        // v (0 to 1) goes from -z to +z
        index = 0;
        maxAxis = absY;
        uc = -z;
        vc = x;

    }

    // POSITIVE Z
    // Up
    else  if (isZPositive && absZ >= absX && absZ >= absY) {
        // u (0 to 1) goes from -x to +x
        // v (0 to 1) goes from -y to +y
        index = 2;
        maxAxis = absZ;
        uc = x;
        vc = y;

    }

    // NEGATIVE Z
    // Down
    else  if (!isZPositive && absZ >= absX && absZ >= absY) {
        // u (0 to 1) goes from +x to -x
        // v (0 to 1) goes from -y to +y

        index = 3;
        maxAxis = absZ;
        uc = -x;
        vc = y;

    }

    // Convert range from -1 to 1 to 0 to 1
    u = 0.5f * (uc / maxAxis + 1.0f);
    v = 0.5f * (vc / maxAxis + 1.0f);
  
    findPixel(index, u, (1 - v), result);
}

void findPixel_EAC(int index, double u, double v, double result[2]){

    double n,m;

    // Left Front Right
    if(index <= 2){

        n = (tileSizeX * (index % 3))  + v * tileSizeX;
      m = u * tileSizeY;

    }
    // Down Back Up
    else{

        switch(index){
            case 3:
                n = u * tileSizeX;
                m = tileSizeY + (1 - v) * tileSizeY;
                break;

            case 4:
                n = tileSizeX + (1 - u) * tileSizeX + 1;
                m = tileSizeY + v * tileSizeY + 1;
                break;

            case 5:
                n = tileSizeX * 2.0f  +  u * tileSizeX;
                m = tileSizeY  + (1 - v) * tileSizeY;
                break;
        }
    }
    result[0] = n - 1;
    result[1] = m - 1;

}

void convert_EAC(double x, double y, double z, double result[2]){

    double maxAxis, uc, vc;
    double u, v;
    int index;

    double absX = fabs(x);
    double absY = fabs(y);
    double absZ = fabs(z);

    int isXPositive = x > 0 ? 1 : 0;
    int isYPositive = y > 0 ? 1 : 0;
    int isZPositive = z > 0 ? 1 : 0;

    // Front
    if (isXPositive && absX >= absY && absX >= absZ) {

        index = 1;
        maxAxis = absX;
        uc = -z;
        vc = y;

    }
    // Back
    else if (!isXPositive && absX >= absY && absX >= absZ) {

        index = 4;
        maxAxis = absX;
        uc = -z;
        vc = -y;

    }
    // Left
    else if (isYPositive && absY >= absX && absY >= absZ) {

        index = 2;
        maxAxis = absY;
        uc = -z;
        vc = -x;

    }
    // Right
    else if (!isYPositive && absY >= absX && absY >= absZ) {

        index = 0;
        maxAxis = absY;
        uc = -z;
        vc = x;

    }
    // Up
    else if (isZPositive && absZ >= absX && absZ >= absY) {

        index = 5;
        maxAxis = absZ;
        uc = x;
        vc = y;

    }
    // Down
    else if (!isZPositive && absZ >= absX && absZ >= absY) {

        index = 3;
        maxAxis = absZ;
        uc = -x;
        vc = y;

    }


    u = 2.0f * atan((uc / maxAxis))/PI + 0.5f;
    v = 2.0f * atan((vc / maxAxis))/PI + 0.5f;

    findPixel_EAC(index,u,v, result);

}

void get_power(bool *done) {
    int fd = open("/sys/devices/3160000.i2c/i2c-0/0-0040/iio_device/in_power1_input", O_RDONLY | O_NONBLOCK);
    if (fd < 0) {
      perror("open()");
      exit(1);
    }
	int cnt = 0;

 	double sum = 0;

 	while (!(*done)) {
      char buf[32];
      lseek(fd, 0, 0);
      int n = read(fd, buf, 32);
      if (n > 0) {
      	buf[n] = 0;
      	char *o = NULL;
      	sum += strtod(buf, &o);
      }
	  cnt++;
	}
	fprintf(stdout, "avg power: %lf (%d)\n", sum/cnt, cnt);
}


int main(int argc, char** argv) {

    int option = argv[6][0] - '0';
    // load image
    Mat image = imread(argv[1], CV_LOAD_IMAGE_COLOR);
    //Mat pat = imread(argv[1], CV_LOAD_IMAGE_COLOR);

	//not lock-free, but whatever...
	bool done = false;
	std::thread first(get_power, &done);
    struct timespec tv;
 	double start = 0;
    clock_gettime(CLOCK_MONOTONIC_RAW, &tv);
    start = tv.tv_sec + tv.tv_nsec * 1e-9;

    // get width and height
    w = image.cols;
    h = image.rows;
    tileSize = w/3.0;
    tileSizeX = w/3.0;
    tileSizeY = h/2.0;

    // parameters for FoV
    //int fovX = 90,fovY = 90,fw = w * (fovX / 360.0) + 1,fh = h * (fovY / 360.0) + 1;
    int fw = atof(argv[2]), fh = atof(argv[3]);
    int fovX = atof(argv[4]), fovY = atof(argv[5]);
    // ht is theta (horizontal), goes toward left first
    // hp is phi (vertical), goes toward up first
    // both are relative rotation angles
    double hp = atof(argv[8]),ht = atof(argv[7]);

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

    // initialize fov image
    Mat fov(fh, fw, CV_8UC3);

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
            if (option == 0) {

                // convert 3D catesian to 2D coordinates
                cartesian2coordinates(p3[0], p3[1], p3[2], res);

            }
            else if(option == 1){

                // convert 3D catesian to cube UVs
                convert_xyz_to_cube_uv(p3[0], p3[1], p3[2], res);

            }
            else {
                // convert 3D catesian to EAC UVs
                convert_EAC(p3[0], p3[1], p3[2], res);

            }

            if (b >= fh || a >= fw){
                break;
            }
            // assign the pixel value
            //fprintf(stdout, "indices: %d, %d\n", nearestNeighbor(res[0]),nearestNeighbor(res[1]));
            fov.at<Vec3b>(b,a) = image.at<Vec3b>(nearestNeighbor (res[1]), nearestNeighbor (res[0]));
            //pat.at<Vec3b>(nearestNeighbor (res[1]), nearestNeighbor (res[0]))[0] = 255;
            //pat.at<Vec3b>(nearestNeighbor (res[1]), nearestNeighbor (res[0]))[1] = 255;
            //pat.at<Vec3b>(nearestNeighbor (res[1]), nearestNeighbor (res[0]))[2] = 255;
        }
        a = 0;
    }

    clock_gettime(CLOCK_MONOTONIC_RAW, &tv);
    double end = tv.tv_sec + tv.tv_nsec * 1e-9;
    fprintf(stderr, "Latency: %.3f ms\n", (end - start) * 1000);
	done = true;;
	first.join();

    // save the fov image
    imwrite("output.jpg", fov);
    //imwrite("input.jpg",pat);

    return 0;
}
