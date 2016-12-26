#pragma once
#include <itkImage.h>
#include <opencv2/core/core.hpp>
#include <vector>
#include "../PINK/mclifo.h"

using std::vector;

typedef unsigned int WEIGHT_TYPE;



namespace YZXV3 {
    
    class PWSegmentation{
        
        
        
    public:
        
        typedef itk::Image<double, 3> itk_3d_double;//if necessary, this can be removed.
        
        PWSegmentation(itk_3d_double::Pointer& inputImage,
                       vector<cv::Point3i> fgSeeds,
                       vector<cv::Point3i> bgSeeds);
        
        ~PWSegmentation();
        
        //build graph including two virtual nodes.
        //我们的排序是这样的：Idx为：0-totalVertex为Dummy vertex，然后后面的作为普通的边。
        // 2*totalVoxelNum - XXX 后面的乱七八糟的六链接，顺序为：
        // 先x方向，再y方向，最后z方向。
        void BuildVertexesAndEdges();
        
        void BuildSeed();//build seed.
        
        void NeighborhoodEdges();
        
        void DeepCopy(itk_3d_double::Pointer& input, itk_3d_double::Pointer& output);
        
        int  Idx(int x, int y, int z);
        void Coordinate(int Idx, int& x, int& y, int& z);
        
        void GetNeighborhood(int vertIdx, int forbidEdgeIdx, vector<int>& neighborhoods);
        
        void GenerateNeighborhood(int edgeIdx, vector<int>& neighborhood);
        
        void GenerateWeight();
        
        double UnarySource(int vertIdx);
        
        double UnarySink(int vertIdx);
        
        double GetPixelValue(int vertIdx);
        
        
        
        //////////////////Functionals from MSF_PW///////////////////
        
        void Solve();
        
        void memory_allocation_PW(Lifo ** LIFO,       /* stack for plateau labeling */
                                  Lifo ** LCP,        /* list of the edges belonging to a plateau */
                                  bool ** indic_E,    /* indicator for edges */
                                  bool ** indic_P,    /* indicator for edges */
                                  int ** indic_VP,    /* indicator for nodes */
                                  int ** Rnk,         /* array needed for union-find efficiency */
                                  int ** Fth,         /* array for storing roots of merged nodes trees */
                                  int ** local_seeds, /* array for storing the index of seeded nodes of a plateau */
                                  int ** LCVP ,       /* list of vertices of a plateau */
                                  int ** Es,          /* array of sorted edges according their reconstructed weight */
                                  int ** NEs,         /* array of sorted edges according their weight */
                                  int N,              /* number of vertices */
                                  int M);             /* number of edges*/

        
        void merge_node (int e1,            /* index of node 1 */
                         int e2,            /* index of node 2 */
                         int * Rnk,         /* array needed for union-find efficiency */
                         int *Fth,          /* array for storing roots of merged nodes trees */
                         float ** proba, /* array for storing the result x */
                         int nb_labels);     /* nb of labels */

        
        void PowerWatershed_q2(int ** edges,              /*array of node indexes composing edges */
                                           uint32_t * weights,        /* reconstructed weights */
                                           uint32_t * normal_weights, /* original weights */
                                           int max_weight,            /* maximum weight value */
                                           int *seeds,                /* array of seeded nodes indexes */
                                           uint8_t * labels,          /* label values on the seeded nodes */
                                           int size_seeds,            /* nb of seeded nodes */
                                           int nb_labels,             /* number of different labels */
                                           bool quicksort,            /* true : bucket sort used; false : stochastic sort o(n log n) */
                               itk_3d_double::Pointer& img_proba); /* output image of potential/proba map x minimizing Epq*/

        
        void PowerWatershed_q2(int ** edges,
                               uint32_t * weights,
                               uint32_t * normal_weights,
                               int maxweight,
                               int *seeds,
                               uint8_t * labels,
                               int size_seeds,
                               int rs,  
                               int cs,  
                               int ds,  
                               int nb_labels,
                               bool quicksort,
                               itk_3d_double::Pointer& img_proba);

        
        
    //protected:
        itk_3d_double::Pointer   inputImage;
        vector<cv::Point3i>      fgSeeds;
        vector<cv::Point3i>      bgSeeds;
        int                      nbLabels;//=2
        
        int                      imgDim[3];
        int                      totalVoxelNum;
        
        vector<vector<int> >     neighborhoodEdgeTable;
        
        //graph-related functional.
        int**                    edges;//edges[0][i]-edges[1][i] represent a valid edge in graph.
        int                      edgeNum;//=M
        int                      edgeXOffset;
        int                      edgeYOffset;
        WEIGHT_TYPE*             weights;//input weights.
        int                      max_weight;
        int*                     seed;//seed idx array.
        int                      seed_size;//size of the seed.
        unsigned char*           labels;//labels of current seeds starting from 1. size of label equals to size of seed.
        
        itk_3d_double::Pointer   resultImg;
        
    };
    
}