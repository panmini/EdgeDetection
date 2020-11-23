#include "BMP.h"

int main()
{
    BMP bmp;
    bmp.readImage("saha.bmp");
    bmp.saveGrayScale("grayScale");
    bmp.finalisation();

    return 0;
}
