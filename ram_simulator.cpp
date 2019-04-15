#include <iostream>
#include <cmath>
#include <queue>
#include <cstdbool>
#include <climits>
#include <cstddef>
#include <cstdint>
#include <iomanip>
#include <fstream>
#include <string>
#include <chrono>
#include <opencv2/opencv.hpp>

using namespace std;


const auto start = std::chrono::system_clock::now();

class RAM_Access{

    public:

    uintptr_t address;
    long long date;
    // Image Test
    cv::Vec3b pixelVal;

    RAM_Access(){}

    RAM_Access(const RAM_Access &old){

        address = old.address;
        date = old.date;

        // Image Test
        // pixelVal = old.pixelVal;
    }

    bool operator != (const RAM_Access &another) const{

        return address != another.address || date != another.date;
    }

};

// ostream& operator <<(ostream &strm, const RAM_Access &ra) {

//     strm << "[Addr: 0x"<< uppercase << hex << setfill('0') << setw(12) << reinterpret_cast<uintptr_t>(ra.address) << ", Freq: " << ra.date << "]" << endl;

// }


template<
    class T,
    class Container = vector<T>,
    class Compare = less<typename Container::value_type>
> class MyQueue : public priority_queue<T, Container, Compare> {

public:

    bool isFull(int RAM_SIZE){

        return (this -> size()) >= RAM_SIZE;
    }

    bool contains(const uintptr_t val){

        auto first = this->c.begin();
        auto last = this->c.end();

        while (first != last) {

            if ((*first).address == val)
                return true;

            ++first;
        }
        return false;
    }

    // Image Test
    // cv::Vec3b get_pixel(const uintptr_t val){

    //     auto first = this->c.begin();
    //     auto last = this->c.end();

    //     while (first != last) {

    //         if ((*first).address == val)
    //             return (*first).pixelVal;

    //         ++first;
    //     }

    //     return 0;
    // }
    //
};


// template<typename A> void print_queue(A& pq)
// {
//     while (!pq.empty())
//     {
//         cout << pq.top() << endl;
//         pq.pop();
//     }
// }


class CompareFreq{
public:

    bool operator()(const RAM_Access& lhs, const RAM_Access& rhs){

        return  lhs.date - rhs.date > 0;

    }
};

// ***********************Global Variable**************************
int RAM_SIZE, mode;

long long SRAM_access = 0, DRAM_access = 0;

long long ram_load = 0;
/*
 *  VR Projection
 */
double PI = 3.14159265358979323846;

int w, h;
int fw, fh;
int fovX, fovY;
double hp, ht;

static MyQueue<RAM_Access, vector<RAM_Access>, CompareFreq> SRAM;
static int **DRAM;
static int **DRAM_output;

// ****************Image Test*********************
// cv::Mat image = cv::imread("/home/rhein/Desktop/ram-simulator/1080p.jpg");
// *********************************************************

long long order;
//ofstream trace;

//****************************************************************

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

//New approach
void coordinates2spherical(double i, double j, double result[2]){
    double theta, phi;

    if (i >= w / 2.0){

        theta = (2 * i * PI / w) - PI;

    }
    else{

        theta = (2 * i * PI / w) + PI;

    }
    phi = j * PI / h;
    //printf("theta: %lf , phi: %lf\n", theta,phi);
    result[0] = theta;
    result[1] = phi;
}

void cartesian2coordinates_inverse(double x, double y, double z, double result[2]){

    double the;
    // pay atentions to atan2 vs atan
    if (x != 0){

        the = atan2(y, x);

    }
    else{

        the = toRadian(90.0);

    }

    double phi = acos(z);

    the = the / PI * 180.0;
    phi = phi / PI * 180.0;

    if(the >= -fovX/2.0 && the <= fovX/2 && phi >= 90 -fovY/2.0 && phi <= 90 +fovY/2.0){

        result[0] = (the + fovX/2.0)* fw /fovX;
        result[1] = (phi -  90  + fovY/2.0) * fh/ fovY;
    }
    else{

        result[0] = 0;
        result[1] = 0;
    }
}

void coordinates2cartesian(double i, double j, double result[3]){

    double angles [] = {0.0, 0.0};

    coordinates2spherical(i, j, angles);

    spherical2cartesian(angles[0],angles[1],result);

}

void matrixMultiplication(double* vector, double matrix[3][3], double res[3]) {

    for (int i=0; i<3; i++) {
        for (int j=0; j<3; j++) {

            res[i] += matrix[i][j] * vector[j];
        }
    }

}
/*
 *  End of VR Algorithm
 */
bool check_SRAM(uintptr_t addr){

    return SRAM.contains(addr);
}

void loadSequentially(int i, int j){

    int row = j, col = 0;

    for (int m = 0; m < w * 2; m++) {

        RAM_Access temp;

        if(col == w - 1){
            col = 0;
            row++;
        }
        // Image test
        // temp.pixelVal = image.at<cv::Vec3b>(row, col);
        //

        temp.address = (uintptr_t)&DRAM[row][col];

        DRAM_access++;

        auto elapsed = chrono::high_resolution_clock::now() - start;
        temp.date = chrono::duration_cast<std::chrono::microseconds>(elapsed).count();


        SRAM.push(temp);

        col++;
    }
}

void loadSRAM(int i, int j, int mode) {

    if (SRAM.isFull(RAM_SIZE)) {
        for (int m = 0; m < w * 2; m++)
            SRAM.pop();
    }

    if(mode == 0){

        RAM_Access temp;
        temp.address = (uintptr_t)&DRAM[j][i];
        auto elapsed = chrono::high_resolution_clock::now() - start;
        temp.date = chrono::duration_cast<std::chrono::microseconds>(elapsed).count();
        SRAM.push(temp);
    }

    else if(mode == 1){

        loadSequentially(i ,j);
    }
    ram_load++;

    

}


int main(int argc, char** argv) {

//    trace.open("trace.trc");
    // parameter for simulator
    // number of pixels sram could store
    // 4.75 MB
    RAM_SIZE = 4.75 * 1024 * 1024 / 3;
    mode = stoi(argv[1]);

    int resolution = stoi(argv[4]);

    switch(resolution){

        case 0:
            // 480
            // input image size
            w = 720;
            h = 480;

            // parameters for FOV
            fw = 220;
            fh = 240;
            break;

        case 1:
            // 720p
            // input image size
            w = 1280;
            h = 720;

            // parameters for FOV
            fw = 392;
            fh = 360;
            break;
        case 2:
            // 1080p
            // input image size
            w = 1920;
            h = 1080;

            // parameters for FOV
            fw = 587;
            fh = 540;
            break;
        case 3:
            // 2K
            // input image size
            w = 2048;
            h = 1080;

            // parameters for FOV
            fw = 626;
            fh = 540;
            break;
        case 4:
            // 4K
            // input image size
            w = 3840;
            h = 2160;

            // parameters for FOV
            fw = 1174;
            fh = 1080;
            break;
    }

    fovX = 110;
    fovY = 90;
    ht = stoi(argv[2]);
    hp = stoi(argv[3]);

    // ****************Image Test*********************
    // cv::Mat fov(fh, fw, CV_8UC3);
// *********************************************************
    RAM_SIZE = (RAM_SIZE / w) * w;

    DRAM = new int*[h];

    for(int i = 0; i < h; i++){

        DRAM[i] = new int[w];

        for(int j = 0; j < w; j++){

            DRAM[i][j] = 0;
        }
    }

    DRAM_output = new int*[fh];

    for(int i = 0; i < fh; i++){

        DRAM_output[i] = new int[fw];

        for(int j = 0; j < fw; j++){

            DRAM_output[i][j] = 0;
        }
    }

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

    // reverse
    if(mode == 0) {

        double rot_y_inverse[3][3] = {
                {cos(hpr),  0, sin(hpr)},
                {0,         1, 0},
                {-sin(hpr), 0, cos(hpr)}};

        double rot_z_inverse[3][3] = {
                {cos(htr), -sin(htr), 0},
                {sin(htr), cos(htr),  0},
                {0,        0,         1}};

        //border on the input frame that map to output pixels

        double maxX = -INFINITY;
        double minX = INFINITY;
        double maxY = -INFINITY;
        double minY = INFINITY;

        double jT = -fovX / 2.0;
        double jR = fovX / 2.0;
        double jB = -fovX / 2.0;
        double jL = 360 - fovX / 2.0;
        double iT = 90 - fovY / 2.0;
        double iR = 90 - fovY / 2.0;
        double iB = 90 + fovY / 2.0;
        double iL = 90 - fovY / 2.0;

        double i = 0.0;
        double j = 0.0;

        for(int k = 0; k < 2*(fh+fw); k++){

            //Top
            if (k < fw){
                i = iT;
                j = jT;
                jT +=  fovX * 1.0 / fw;
            }

            //Right
            if ((k >= fw) && (k < fw+fh)){
                i = iR;
                j = jR;
                iR += fovY * 1.0 / fh;
            }

            //Bottom
            if ((k >= fw+fh) && (k < 2*fw+fh)){
                i = iB;
                j = jB;
                jB +=  fovX * 1.0 / fw;
            }

            //Left
            if ((k >= 2*fw+fh) && (k < 2*(fw+fh))){
                i = iL;
                j = jL;
                iL += fovY * 1.0 / fh;
            }

            // rotation along y axis
            double p1[] = {0.0, 0.0, 0.0};
            spherical2cartesian(toRadian((j < 0)? (j + 360): j), toRadian((i < 0) ? (i + 180) : i), p1);

            double p2[] = {0.0, 0.0, 0.0};
            matrixMultiplication(p1, rot_y, p2);

            // rotate along z axis
            double p3[] = {0.0, 0.0, 0.0};
            matrixMultiplication(p2, rot_z, p3);

            double res[] = {0.0, 0.0};

            // convert 3D catesian to 2D coordinates
            cartesian2coordinates(p3[0], p3[1], p3[2], res);

            if (b >= fh) break;

            if (minX > res[0]) minX = res[0];
            if (maxX < res[0]) maxX = res[0];
            if (minY > res[1]) minY = res[1];
            if (maxY < res[1]) maxY = res[1];

        }

        if (hp <= -45 || hp >= 315){

            maxY = h;
            maxX = w;
            minX = 0.0;

        }
        if(hp >= 45) {

            minY = 0.0;
            maxX = w;
            minX = 0.0;

        }

        //for input pixel in the output range, calculate the outpout cordinnates
        int x , y;

        for (y = 0; y < h; y++){
            for (x = 0; x < w; x++){

                //if pixel map to output get input index
                if (x <= maxX && x >= minX && y <= maxY && y >= minY){

                    // ****************Image Test*********************
                    double cartesian []  ={0.0, 0.0, 0.0};
                    coordinates2cartesian(x, y, cartesian);

                    double p1[] = {0.0, 0.0, 0.0};
                    matrixMultiplication(cartesian, rot_z_inverse , p1);

                    // rotate along z axis
                    double p2[] = {0.0, 0.0, 0.0};
                    matrixMultiplication( p1, rot_y_inverse, p2);


                    double res[] = {0.0,0.0};
                    cartesian2coordinates_inverse(p2[0], p2[1], p2[2], res);

                    int temp_x = nearestNeighbor(res[0]);
                    int temp_y = nearestNeighbor(res[1]);

                    //fov.at<Vec3b>(nearestNeighbor(res[1]), nearestNeighbor(res[0])) = image.at<Vec3b>(y,x);
                    // ************************************************

                    uintptr_t addr = (uintptr_t)&DRAM[y][x];

                    // if(!check_SRAM(addr)){

                    loadSRAM(x, y, mode);

                    //     SRAM.set_freq(addr);

                    // }
                    
                    // cout << "0x"<< uppercase << hex << setfill('0') << setw(12) << reinterpret_cast<uintptr_t>(SRAM.get_addr(addr)) << " P_MEM_RD " << dec << order++ << endl;

                    // Image Test
                    // fov.at<cv::Vec3b>(temp_y, temp_x) = SRAM.get_pixel(addr);

                    // DRAM_output[temp_y][temp_x] = DRAM[y][x]; 
                    cout << "0x" << uppercase << hex <<  setfill('0') << setw(8) << reinterpret_cast<uintptr_t>(addr) << " P_MEM_RD ";
                    cout << dec << chrono::duration_cast<std::chrono::milliseconds>(chrono::high_resolution_clock::now() - start).count() << " ms" << endl;
                    SRAM_access++;
                }
            }
        }

        for(int i = 0; i < fh; i++){

            for(int j = 0; j < fw; j++){

                uintptr_t addr_out = (uintptr_t)&DRAM_output[i][j];
                
                cout <<"0x" << uppercase << hex <<  setfill('0') << setw(8) << reinterpret_cast<uintptr_t>(addr_out) << " P_MEM_WR ";
                cout << dec << chrono::duration_cast<std::chrono::milliseconds>(chrono::high_resolution_clock::now() - start).count() << " ms" << endl;
            }
        }
//        cout << "Max X, Y: " << maxX << ", " << maxY << " Min X, Y: " << minX << ", " << minY << endl;
    }
//
    else {

        //    int pixel_count = 0;
        // default head orientation is 0,90
        for (double i = 90 - fovY / 2.0; i < 90 + fovY / 2.0; i += fovY * 1.0 / fh, b++) {
            for (double j = -fovX / 2.0; j < fovX / 2.0; j += fovX * 1.0 / fw, a++) {

                double p1[] = {0.0, 0.0, 0.0};
                spherical2cartesian(toRadian((j < 0) ? j + 360 : j), toRadian((i < 0) ? (i + 180) : i), p1);
                // rotation along y axis
                double p2[] = {0.0, 0.0, 0.0};
                matrixMultiplication(p1, rot_y, p2);

                // rotate along z axis
                double p3[] = {0.0, 0.0, 0.0};
                matrixMultiplication(p2, rot_z, p3);

                double res[] = {0.0, 0.0};

                // convert 3D catesian to 2D coordinates
                cartesian2coordinates(p3[0], p3[1], p3[2], res);

                if (b >= fh || a >= fw) {
                    break;
                }

                int temp_x = nearestNeighbor(res[0]);
                int temp_y = nearestNeighbor(res[1]);

                uintptr_t addr = (uintptr_t)&DRAM[temp_y][temp_x];

                if(!check_SRAM(addr)){

                    loadSRAM(temp_x, temp_y, mode);
                }

                //   auto elapsed = chrono::high_resolution_clock::now() - start;
                //   temp.date = chrono::duration_cast<std::chrono::microseconds>(elapsed).count();    
               
                cout << "0x" << uppercase << hex <<  setfill('0') << setw(8) << reinterpret_cast<uintptr_t>(addr) << " P_MEM_RD ";
                cout << dec << chrono::duration_cast<std::chrono::milliseconds>(chrono::high_resolution_clock::now() - start).count() << " ms" << endl;

                // Image Test
                // fov.at<cv::Vec3b>(b, a) = SRAM.get_pixel(addr);
                // DRAM_output[b][a] = DRAM[temp_y][temp_x];
                SRAM_access++;

                uintptr_t addr_out = (uintptr_t)&DRAM_output[b][a];

                cout <<"0x" << uppercase << hex <<  setfill('0') << setw(8) << reinterpret_cast<uintptr_t>(addr_out) << " P_MEM_WR ";
                cout << dec << chrono::duration_cast<std::chrono::milliseconds>(chrono::high_resolution_clock::now() - start).count() << " ms" << endl;
            }
            a = 0;
        }
    }


//
   // printf("DRAM Access: %d\n", DRAM_access);
   // printf("SRAM Access: %d\n", SRAM_access);
   // printf("SRAM Load: %d\n", ram_load);

    // imwrite("output.jpg", fov);
    return 0;
}
