/*Copyright ESIEE (2009)
 
 Author :
 Camille Couprie (c.couprie@esiee.fr)
 
 Contributors :
 Hugues Talbot (h.talbot@esiee.fr)
 Leo Grady (leo.grady@siemens.com)
 Laurent Najman (l.najman@esiee.fr)
 
 This software contains some image processing algorithms whose purpose is to be
 used primarily for research.
 
 This software is governed by the CeCILL license under French law and
 abiding by the rules of distribution of free software.  You can  use,
 modify and/ or redistribute the software under the terms of the CeCILL
 license as circulated by CEA, CNRS and INRIA at the following URL
 "http://www.cecill.info".
 
 As a counterpart to the access to the source code and  rights to copy,
 modify and redistribute granted by the license, users are provided only
 with a limited warranty  and the software's author,  the holder of the
 economic rights,  and the successive licensors  have only  limited
 liability.
 
 In this respect, the user's attention is drawn to the risks associated
 with loading,  using,  modifying and/or developing or reproducing the
 software by the user in light of its specific status of free software,
 that may mean  that it is complicated to manipulate,  and  that  also
 therefore means  that it is reserved for developers  and  experienced
 professionals having in-depth computer knowledge. Users are therefore
 encouraged to load and test the software's suitability as regards their
 requirements in conditions enabling the security of their systems and/or
 data to be ensured and,  more generally, to use and operate it in the
 same conditions as regards security.
 
 The fact that you are presently reading this means that you have had
 knowledge of the CeCILL license and that you accept its terms.
 */


#include <stdio.h>
#include <sys/types.h>
#include <sys/time.h>
#include <time.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include "../PINK/mccodimage.h"
#include "../PINK/mcimage.h"
#include "../include/cccodimage.h"
#include "../include/lMSF.h"
#include "../include/MSF_RW.h"
#include "../include/powerwatsegm.h"
#include "../include/image_toolbox.h"
#include <unistd.h>
//#include "../argv/argv.h"
#include "../include/YZXHeader.h"
#include "../include/YZX2DPW.h"

//#include <fstream>
#include <iostream>
#include <vector>

#include <glog/logging.h>

#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>

using namespace std;
using namespace cv;


/* =============================================================== */
int main(int argc, char **argv)
/* =============================================================== */
{
    google::InitGoogleLogging(argv[0]);
    google::LogToStderr();
    
    //start my own PW program.
#if 1
    string fileName = "/Users/yanzixu/Downloads/TestImage.bmp";
    Mat input = cv::imread(fileName, CV_LOAD_IMAGE_GRAYSCALE);
    
    Mat outResult;
    
    vector<Point2i> fgSeed;
    fgSeed.push_back(Point2i(150,150));
    vector<Point2i> bgSeed;
    bgSeed.push_back(Point2i(1,1));
    
    YZX::YZX2DPowerWatershed(input, fgSeed, bgSeed, outResult);
    imshow("wnd", outResult*254);
    waitKey();
#endif
    
    //generate image and seed.
//    {
//        string fileName = "/Users/yanzixu/Downloads/TestImage.bmp";
//        Mat input = cv::imread(fileName, CV_LOAD_IMAGE_COLOR);
//        cv::imwrite("TestImage.ppm", input);
//    }
//    
//    string fileName = "/Users/yanzixu/Downloads/TestImage.bmp";
//    Mat input = cv::imread(fileName, CV_LOAD_IMAGE_GRAYSCALE);
//    input.setTo(Scalar(0));
//    input.at<uchar>(150,150) = 255;
//    input.at<uchar>(1,1) = 100;
//    cv::imwrite("TestSeed.ppm", input);
    
    int algo =-1;
    char * image_name = NULL;
    char * seeds_name =NULL;
    char * multi_seeds_name =NULL;
    char * overlay_name =NULL;
    char * output_name =NULL;
    bool mult = false;
    bool geod = false;
//    argv_t args[] = {
//        { 'a', 0L, ARGV_INT, &algo, (char*)" 1|2|3", (char*)"algo: Kruskal(1) PW q=2(2) Prim(3)" },
//        { 'i', 0L, ARGV_CHAR_P, &image_name, (char*)" *.pgm |*.ppm",  (char*)"image_name" },
//        { 's', 0L, ARGV_CHAR_P, &seeds_name, (char*)" *.pgm",  (char*)"seed_name (see README)" },
//        { 'm', 0L, ARGV_CHAR_P, &multi_seeds_name, (char*)" *.pgm",  (char*)"multi_seed_name (see README)" },
//        { 'g', 0L, ARGV_BOOL, &geod, (char*)" 0|1", (char*)"geod: reconstruction of the weights(1)" },
//        { 'o', 0L, ARGV_CHAR_P, &output_name, (char*)" *.pgm",  (char*)"output mask name" },
//        { 'v', 0L, ARGV_CHAR_P, &overlay_name, (char*)" *.ppm",  (char*)"output overlay name" },
//        { ARGV_LAST, 0L, 0, 0L, 0L, 0L}
//    };
//    
//    argv_process(args, argc, argv);
    
    
//    cout<<"image_name : "<<(string)image_name<<endl;
//    
//    if(NULL != seeds_name)
//    cout<<"seed_name : "<<(string)seeds_name<<endl;
//    
//    if(NULL != multi_seeds_name)
//    cout<<"multi_seeds_name : "<<(string)multi_seeds_name<<endl;
//    
//    if(NULL != overlay_name)
//    cout<<"overlay_name : "<<(string)overlay_name<<endl;
//    
//    if(NULL != output_name)
//    cout<<"output_name : "<<(string)output_name<<endl;
//    
//    cout<<"geod : "<<geod<<endl;
//    cout<<"mult : "<<mult<<endl;
    
    //default: mult = false, geodesic = false;
    //input of image name and multi_seeds name.
//    image_name = "/Users/yanzixu/Downloads/PW_1.0.1/BUILD/Release/images/2D/241004.ppm";
//    seeds_name = "/Users/yanzixu/Downloads/PW_1.0.1/BUILD/Release/images/2D/seeds_241004_MULT.pgm"; //this seed image contains 15 different labels from 1 to 15
    
    image_name = "/Users/yanzixu/Downloads/PW_1.0.1/BUILD/Release/TestImage.ppm";
    seeds_name = "/Users/yanzixu/Downloads/PW_1.0.1/BUILD/Release/TestSeed.ppm";
    
    algo = 2;//PW
    
    if(seeds_name ==NULL)
    {
        seeds_name= multi_seeds_name;
        mult = true;
    }
    
//    Mat img = cv::imread(multi_seeds_name, CV_LOAD_IMAGE_GRAYSCALE);
//    double min,max;
//    min = max = -1;
//    cv::minMaxLoc(img, &min, &max);
//    cout<<"min = "<<min<<" max = "<<max<<endl;
    //cv::imshow("output", img * 17);
    //waitKey(0);
    
    
//    if ((algo==-1) ||(image_name ==NULL  )||(seeds_name ==NULL))
//    {
//        fprintf(stderr, "usage: %s -a algo -i image(.pgm or .ppm) <-s seeds_2_labels.pgm | -m seeds_m_labels.pgm > \n", argv[0]);
//        fprintf(stderr, "options : [-g geod] [-o output mask name] [-v image overlay name]\n");
//        fprintf(stderr, "algo : 1 : Kruskal \n");
//        fprintf(stderr, "       2 : PW q=2 \n");
//        fprintf(stderr, "       3 : Prim \n");
//        fprintf(stderr, "seeds[MULT].pgm  (see the README for more precision) \n");
//        fprintf(stderr, "geod : 1 : geodesic reconstruction \n");
//        fprintf(stderr, "       0 : no reconstruction \n");
//        exit(1);
//    }
    
    int t1=clock();
    bool quicksort = false;
    int32_t nblabels, i,j;
    struct xvimage * image_r;
    struct xvimage * image_v;
    struct xvimage * image_b;
    
    struct xvimage * output = NULL;
    struct xvimage * seeds;
    
    unsigned char * s;
    int rs, cs, ds, N, M;
    
    // Reading the seed image
    
    seeds = readimage(seeds_name);
    if (seeds == NULL) { fprintf(stderr, "msf_rw: readimage failed\n"); exit(1); }
    s = UCHARDATA(seeds);
    
    rs = rowsize(seeds);
    cs = colsize(seeds);
    ds = colsize(seeds);
    bool color=false;
    
    //这边读进来的图像变成了3个通道的。
    int size = strlen(image_name);
    ds=1; //default
    if (strcmp(&image_name[size-3], "pgm")==0) // grey levels image
    {
        image_r = readimage(image_name);
        ds = depth_pw(image_r);
    }
    else
    {
        color = true;
        readrgbimage(image_name,  &image_r, &image_v, &image_b); //green in french is vert.
        Mat rImg, bImg, gImg;
        xvImage2CVMat(image_r, rImg);
        xvImage2CVMat(image_v, gImg);
        xvImage2CVMat(image_b, bImg);
        
        //imshow("rImg", rImg);
        //imshow("gImg", gImg);
        //imshow("bImg", bImg);
        //waitKey();
    }
    
    
    
    
    
    //rs : rows. cs : cols. ds: depth.
    N = rs * cs * ds;
    M = ds*rs*(cs-1)+ds*(rs-1)*cs+(ds-1)*cs*rs;  //number of edges 6连通的边
    int ** edges;
    
    int * index_seeds = (int*)malloc(N*sizeof(int));
    uint8_t * index_labels = (uint8_t*)malloc(N*sizeof(uint8_t));
    j=0;
    
    //s是seed iamge的内部data,一个 uchar*的数组。
    
    //multilabel seed image
    if (mult == true)
    {
        nblabels = 0;
        for (i=0;i<rs*cs*ds;i++)
            if(s[i]>0)
            {
                index_seeds[j]=i;
                index_labels[j]=s[i];
                j++;
                if(s[i]>nblabels) nblabels = s[i];
            }
    }
    else
    {
        nblabels=2;
        for (i=0;i<rs*cs*ds;i++)
        {
            if(s[i]>155)
            {
                index_seeds[j] = i;
                index_labels[j]=1;
                j++;
            }
            else if(s[i] == 100)
            {
                index_seeds[j] = i;
                index_labels[j]=2;
                j++;
            }
        }
    }
    int size_seeds = j;
    freeimage(seeds);
    
    cout<<"nlabel = "<<nblabels<<" seed size = "<<j<<endl;
    
    
    //我不太懂，这边他分配的edges的id的时候，两个方向的id是不一样的，为什么要这样？
    //我知道了，这个edge[0][i]和edge[1][i]分别代表第i个edge的两个顶点的idx.
    //顶点的idx用的是一般的那种方法编号的，行x列x深度依次序编号。
    edges =  (int**)malloc(2*sizeof(int*));
    for(i=0;i<2;i++) edges[i] = (int*)malloc(M*sizeof(int));
    compute_edges(edges,rs, cs, ds);
    
    cout<<"edge number = "<<M<<endl;
    
    if (algo == 2) // Kruskal & RW on plateaus multiseeds linear time
    {
        //这个就是PW的代码位置。
        //M是边的个数。N是顶点的个数
        struct xvimage * img_proba;
        uint32_t * weights = (uint32_t *)malloc(sizeof(uint32_t)*M);
        uint32_t * normal_weights ;
        uint32_t max_weight = 255;
        
        //这边似乎是算一个edge weight，大致似乎是将seed点附近的weight (和梯度有关)扩散到全图？
        //我觉得需要做一个实验来确认一下是不是这样。
        if (color == true)
        {
            normal_weights = color_standard_weights_PW( image_name , weights, edges, index_seeds, size_seeds, &max_weight, quicksort);
        }
        else
        {
            normal_weights = grey_weights_PW(image_name, edges,index_seeds, size_seeds, weights, quicksort);
        }
#ifndef SPEEDY
        img_proba = allocimage(NULL, rs, cs, ds, VFF_TYP_1_BYTE);
#endif
        //test weight.
        unsigned char* img_r = UCHARDATA(image_r);
        for(int i=0; i<M; i++)
        {
            int e1 = edges[0][i];
            int e2 = edges[1][i];
            
            //cout<<"e1 = "<<e1<<" e2 = "<<e2<<endl;
            double diff = fabs(img_r[e1] - img_r[e2]);
            double cweight = exp(-1.0 * diff * diff * (1.0/ (255.0*255.0)))*100000;
            weights[i] = cweight;
            normal_weights[i] = cweight;
            max_weight = 100000;
            
            //cout<<"weight = "<<weights[i]<<endl;
        }
        
        for(int i=0; i<size_seeds; i++)
        {
            cout<<"seed "<<i<<" idx = "<<index_seeds[i]<<" label = "<<(int)index_labels[i]<<endl;
        }
        
        /*
         edges =  (int**)malloc(2*sizeof(int*));
         for(i=0;i<2;i++) edges[i] = (int*)malloc(M*sizeof(int)); //edges是两个通道，2M，正反应该是一样的（无向图）。
         
         */
        
        //main functional.
        if (geod ==true)
            output = PowerWatershed_q2(edges, weights, weights, max_weight,index_seeds, index_labels, size_seeds,rs, cs, ds, nblabels, quicksort, img_proba);
        else
            output = PowerWatershed_q2(edges, weights, normal_weights,max_weight,index_seeds, index_labels, size_seeds,rs, cs, ds, nblabels, quicksort, img_proba);
#ifndef SPEEDY
        writeimage(img_proba, (char*)"proba.pgm");
        freeimage(img_proba);
#endif
        free(weights);
        free(normal_weights);
        
        cout<<"segmentation complete,"<<endl;
        
    }else{
        printf("%s\n","error, only PW is allowed.");
        return 1;
    }
    
    int t2=clock();
    assert(output != NULL);
    
    free(index_seeds);
    free(index_labels);
    if (output_name == NULL) 
        output_name =(char*)"mask.pgm"; 
    writeimage(output, output_name);
    
    // overlay for 2D images only
    if (ds==1){
        if (overlay_name == NULL) 
            overlay_name =(char*)"overlay.ppm"; 
        overlay(algo, image_r, image_v, image_b, output, color, overlay_name);
    }
    if (color) 
    {
        freeimage(image_v);
        freeimage(image_b);
    }
    
    freeimage(image_r);
    freeimage(output);
    
    printf("Computation time : %.6lf seconds elapsed\n", ((double)t2-t1)/CLOCKS_PER_SEC);
    return 0;
} 


