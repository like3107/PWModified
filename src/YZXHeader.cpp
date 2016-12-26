#include "../include/YZXHeader.h"
#include "../PINK/mcimage.h"
#include "../PINK/mccodimage.h"

using namespace cv;

void xvImage2CVMat(xvimage* input, cv::Mat& output)
{
    //注意，他的xvImage的row和col是和cv的反过来的（就是他的是一般图像默认的那种row和col排序的）。
    int nRows = input->col_size;
    int nCols = input->row_size;
    output.create(nRows, nCols, CV_8UC1);
    
    for (int x=0; x<nRows; x++) {
        for (int y=0; y<nCols; y++) {
            output.at<uchar>(x,y) = pixel(input,y,x);
        }
    }

}

void Data2CVMat(unsigned int* data, int nRows, int nCols, cv::Mat& output)
{
    uint* dataPtr = data;
    output.create(nRows, nCols, CV_8UC1);
    
    for (int x=0; x<nRows; x++) {
        for (int y=0; y<nCols; y++) {
            output.at<uchar>(x,y) = saturate_cast<int>( *dataPtr++ );
        }
    }

}





























