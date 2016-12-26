#pragma once
#include <iostream>
#include <opencv2/core/core.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/highgui/highgui.hpp>

struct xvimage;

void xvImage2CVMat(xvimage* input, cv::Mat& output);

void Data2CVMat(unsigned int* data, int nRows, int nCols, cv::Mat& output);