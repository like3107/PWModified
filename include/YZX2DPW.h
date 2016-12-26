#pragma once
#include <opencv2/core/core.hpp>
using namespace cv;

namespace YZX {
    
    void YZX2DPowerWatershed(Mat& input, vector<Point2i> fgSeeds, vector<Point2i> bgSeeds, Mat& output);
    
}