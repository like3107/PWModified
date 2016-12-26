#include "../include/YZX2DPW.h"
#include "../include/YZXPWSegmentation.h"
#include <glog/logging.h>

namespace YZX {
    
    void YZX2DPowerWatershed(Mat& input, vector<Point2i> _fgSeeds, vector<Point2i> _bgSeeds, Mat& output)
    {
        typedef itk::Image<double, 3> itk_3d_double;
        itk_3d_double::Pointer inputItk = itk_3d_double::New();
        
        itk_3d_double::IndexType coord;
        itk_3d_double::SizeType size;
        coord[0] = coord[1] = coord[2] = 0;
        size[0] = input.cols; size[1] = input.rows; size[2] = 1;
        
        itk_3d_double::RegionType rg;
        rg.SetIndex(coord);
        rg.SetSize(size);
        inputItk->SetRegions(rg);
        inputItk->Allocate();
        inputItk->FillBuffer(0);
        
        for (int x=0; x<input.rows; x++) {
            coord[1] = x;
            for (int y=0; y<input.cols; y++) {
                coord[0] = y;
                inputItk->SetPixel(coord, static_cast<double>(input.at<uchar>(x,y)));
            }
        }
        
        vector<Point3i> fgSeeds, bgSeeds;
        for (int i=0; i<_fgSeeds.size(); i++) {
            fgSeeds.push_back(Point3i(_fgSeeds[i].x,_fgSeeds[i].y,0));
        }
        for (int i=0; i<_bgSeeds.size(); i++) {
            bgSeeds.push_back(Point3i(_bgSeeds[i].x,_bgSeeds[i].y,0));
        }

        PWSegmentation pwSeg(inputItk, fgSeeds, bgSeeds);
        pwSeg.Solve();
        
        //output
        output.create(input.rows, input.cols, CV_8UC1);
        size = pwSeg.resultImg->GetLargestPossibleRegion().GetSize();
        CHECK_EQ(input.rows, size[1]);
        CHECK_EQ(input.cols, size[0]);
        for (int x=0; x<input.rows; x++) {
            for (int y=0; y<input.cols; y++) {
                coord[0] = y;
                coord[1] = x;
                output.at<uchar>(x,y) = saturate_cast<uchar>(pwSeg.resultImg->GetPixel(coord));
            }
        }

    }
    
}