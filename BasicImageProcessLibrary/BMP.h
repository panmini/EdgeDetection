#include<iostream>
#include<string>
#include<fstream>
#include <cstdint>
#include <cstring>
#include <math.h>
#include <algorithm>
#include "FileInfo.h"
#include "ImageInfo.h"

typedef unsigned char BYTE;
typedef unsigned short WORD;
typedef uint32_t DWORD;

class BMP
{
    public:
        BMP();

        FileInfo * fileinfo;
        ImageInfo * imageinfo;

        void readImage(string filename); //Reads Img from BMP File
        void saveGrayScale(string filename); //Convert bmp to intensity
        void saveGrayScale(BYTE * ,string filename);//Saves Intensity to same path as bmp.h (filename is the name of the file) 
        void saveImage(string filename);//Save the BMP file
        
        BYTE* zoom(int x1, int y1, int x2, int y2, int w, int h, BYTE* buffer1);
        BYTE* Conv(BYTE* raw, int width, int height, float* convMatrix, int cW);
        BYTE* SobelFiltering(BYTE* img, int w,int h,float* edge_direction);
        void NonMaxSuppression(float* edge_direction, float* edge_magnitude,BYTE* value,int width, int height);
        void Hysteresis(BYTE lowThreshold, BYTE highThreshold, int w, int h, BYTE* value);
        void houghTransform(BYTE* value, int w, int h,float* edge_direction);
        //int Transform(BYTE* value, int w, int h);
        void Line( float x1, float y1, float x2, float y2);
        void finalisation(); //Final function includes other func.
      
    private:


        BYTE * data; //Points RGB Image
        BYTE * grayData;//Points Intensity Image


        //BYTE * tempMatrix;

        int padding; //Width or Height must be dividible by 4 while operating img file conversions.
};