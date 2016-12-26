#include "../include/YZXPWSegmentation.h"
#include "../include/MSF_RW_YZX.h"

#include <itkImageRegionIterator.h>
#include <itkImageConstIterator.h>

#include <glog/logging.h>
#include <string>

#include "../PINK/mccodimage.h"
#include "../PINK/mcimage.h"
#include "../PINK/mcutil.h"
#include "../include/cccodimage.h"
#include "../include/MSF_RW.h"
#include "../include/union_find.h"
#include "../include/ccsort.h"
#include "../include/powerwatsegm.h"
#include "../include/random_walker.h"

namespace YZX {
    
    void PWSegmentation::Solve()
    {
        resultImg = itk_3d_double::New();
        this->PowerWatershed_q2(edges, weights, weights, max_weight, seed, labels, seed_size, nbLabels, false, resultImg);
    }
    
    ///////////////////Functionals from MSF_RW ////////////////////////////
    void PWSegmentation::memory_allocation_PW(Lifo ** LIFO,       /* stack for plateau labeling */
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
                              int M)             /* number of edges*/
    {
        std::string F_NAME = "memory_allocation_PW";
        
        *LIFO = CreeLifoVide(M);
        if (*LIFO == NULL) { fprintf(stderr, "%s : CreeLifoVide failed\n", F_NAME.c_str()); exit(0); }
        
        *LCP = CreeLifoVide(M);
        if (*LCP == NULL) { fprintf(stderr, "%s : CreeLifoVide failed\n", F_NAME.c_str()); exit(0); }
        
        *indic_E = (bool*)calloc(M ,sizeof(bool));
        if (*indic_E == NULL) { fprintf(stderr, "%s: calloc indic_E failed\n", F_NAME.c_str()); exit(0);}
        *indic_P = (bool*)calloc(M ,sizeof(bool));
        if (*indic_P == NULL) { fprintf(stderr, "%s: calloc indic_P failed\n", F_NAME.c_str()); exit(0);}
        *indic_VP = (int*)calloc(N ,sizeof(int));
        if (*indic_VP == NULL) { fprintf(stderr, "%s: calloc indic_VP failed\n", F_NAME.c_str()); exit(0);}
        *Rnk = (int*)calloc(N, sizeof(int));
        if (*Rnk == NULL) { fprintf(stderr, "%s: malloc Rnk failed\n", F_NAME.c_str()); exit(0);}
        *Fth = (int*)malloc(N*sizeof(int));
        if (*Fth == NULL) { fprintf(stderr, "%s: malloc Fth failed\n", F_NAME.c_str()); exit(0);}
        
        *local_seeds = (int*)malloc(N*sizeof(int));
        if (*local_seeds == NULL) { fprintf(stderr, "%s: malloc local_seeds failed\n", F_NAME.c_str()); exit(0);}
        
        *LCVP = (int*)malloc(N*sizeof(int)); // vertices of a plateau.
        if (*LCVP == NULL) { fprintf(stderr, "%s: malloc LCVP failed\n", F_NAME.c_str()); exit(0);}
        
        *Es = (int*)malloc(M*sizeof(int));
        if (*Es == NULL) { fprintf(stderr, "%s: malloc Es failed\n", F_NAME.c_str()); exit(0);}
        
        *NEs = (int*)malloc(M*sizeof(int));
        if (*NEs == NULL) { fprintf(stderr, "%s: malloc NEs failed\n", F_NAME.c_str()); exit(0);}
    }
    
    void PWSegmentation::merge_node (int e1,            /* index of node 1 */
                                     int e2,            /* index of node 2 */
                                     int * Rnk,         /* array needed for union-find efficiency */
                                     int *Fth,          /* array for storing roots of merged nodes trees */
                                     float ** proba, /* array for storing the result x */
                                     int nb_labels)     /* nb of labels */
    {
        int k,re1, re2;
        re1 = element_find(e1, Fth );
        re2 = element_find(e2, Fth );
        
        if ((re1 != re2) && (!(proba[0][re1]>=0 && proba[0][re2]>=0)))
        {
            element_link(re1,re2, Rnk, Fth);
            if (proba[0][re2]>=0 && proba[0][re1]<0)
                for(k=0;k<nb_labels-1;k++)
                    proba[k][re1]= proba[k][re2];
                    else if (proba[0][re1]>=0 && proba[0][re2]<0)
                        for(k=0;k<nb_labels-1;k++)
                            proba[k][re2]= proba[k][re1];
        }

    }

    void PWSegmentation::PowerWatershed_q2(int ** edges,              /*array of node indexes composing edges */
                                                uint32_t * weights,        /* reconstructed weights */
                                                uint32_t * normal_weights, /* original weights */
                                                int max_weight,            /* maximum weight value */
                                                int *seeds,                /* array of seeded nodes indexes */
                                                uint8_t * labels,          /* label values on the seeded nodes */
                                                int size_seeds,            /* nb of seeded nodes */
                                                int nb_labels,             /* number of different labels */
                                                bool quicksort,            /* true : bucket sort used; false : stochastic sort o(n log n) */
                                                itk_3d_double::Pointer& img_proba) /* output image of potential/proba map x minimizing Epq*/
    {
        std::string F_NAME = "PowerWatershed_q2";
        int i, j, k, x, y, e1, e2, re1,re2, p, xr;
        int N = this->totalVoxelNum + 2;                /* number of vertices */
        int M = this->edgeNum;    /* number of edges*/
        double val;
        int argmax;
        int nb_vertices, e_max, Ne_max, nb_edges, Nnb_edges;
        bool success, different_seeds;
        uint32_t wmax;
        Lifo * LIFO;
        Lifo * LCP;
        bool * indic_E;
        bool * indic_P;
        int * indic_VP;
        int * Rnk;
        int * Fth;
        int * local_seeds;
        int * LCVP ;
        int * Es;
        int * NEs ;
        
        memory_allocation_PW( &LIFO, &LCP, &indic_E,  &indic_P,  &indic_VP, &Rnk,  &Fth, &local_seeds, &LCVP, &Es, &NEs, N, M);
        
        float **proba = (float **)malloc((nb_labels-1) *sizeof(float*));
        for (i=0; i<nb_labels-1; i++)
        {
            proba[i]= (float *)malloc(N *sizeof(float));
            for (j=0; j<N; j++) proba[i][j]=-1;
        }
        int ** edgesLCP = (int**)malloc(2*sizeof(int*));
        if (edgesLCP == NULL) { fprintf(stderr, "%s: malloc edgesLCP failed\n", F_NAME.c_str()); exit(0);}
        for(k=0;k<2;k++)
        {
            edgesLCP[k] = (int*)malloc(M*sizeof(int));
            if ( edgesLCP[k]== NULL)
            { fprintf(stderr, "%s: malloc edgesLCP failed\n", F_NAME.c_str()); exit(0);}
        }
        for (i=0;i<size_seeds;i++)
            for (j=0;j<nb_labels-1;j++)
            {
                if (labels[i]==j+1)
                    proba[j][seeds[i]] = 1;
                else proba[j][seeds[i]] = 0;
            }
        
        for(k=0;k<N;k++) Fth[k]=k;
        
        float ** local_labels = (float**)malloc((nb_labels-1)*sizeof(float*));
        if (local_labels == NULL) { fprintf(stderr, "%s: malloc local_labels failed\n", F_NAME.c_str()); exit(0);}
        for (i=0; i<nb_labels-1; i++)
        {
            local_labels[i]= (float *)malloc(N *sizeof(float));
            if ( local_labels[i]== NULL)
            { fprintf(stderr, "%s: malloc local_labels failed\n", F_NAME.c_str()); exit(0);}
        }
        
        uint32_t * sorted_weights = (uint32_t*)malloc(M*sizeof(uint32_t));
        if (sorted_weights == NULL) { fprintf(stderr, "%s: malloc sorted_weights failed\n", F_NAME.c_str()); exit(0);}
        
        for(k=0;k<M;k++)
        {
            sorted_weights[k]=weights[k];
            Es[k]=k;
        }
        if (quicksort == true)
            BucketSort(sorted_weights, Es, M, max_weight+1);
        else
            TriRapideStochastique_dec(sorted_weights,Es, 0,M-1);
        int cpt_aretes = 0;
        int Ncpt_aretes = 0;
        
        /* beginning of main loop */
        while (cpt_aretes < M)
        {
            do
            {
                e_max=Es[cpt_aretes];
                cpt_aretes=cpt_aretes+1;
                if(cpt_aretes==M) break;
            }while(indic_E[e_max]==true);
            
            if(cpt_aretes==M) break;
            
            //1. Computing the edges of the plateau LCP linked to the edge e_max
            LifoPush(LIFO, e_max);
            indic_P[e_max]=true;
            indic_E[e_max]=true;
            LifoPush(LCP, e_max);
            nb_vertices=0;
            nb_edges = 0;
            wmax = weights[e_max];
            
            // 2. putting the edges and vertices of the plateau into arrays
            while (! LifoVide(LIFO))
            {
                x = LifoPop(LIFO);
                e1 = edges[0][x]; e2 = edges[1][x];
                re1 = element_find(e1, Fth );
                re2 = element_find(e2, Fth );
                if (proba[0][re1]<0 || proba[0][re2]<0)
                {
                    if (indic_VP[e1]==0)
                    {
                        LCVP[nb_vertices]=e1;
                        nb_vertices++;
                        indic_VP[e1]=1;
                    }
                    if (indic_VP[e2]==0)
                    {
                        LCVP[nb_vertices]=e2;
                        nb_vertices++;
                        indic_VP[e2]=1;
                    }
                    edgesLCP[0][ nb_edges] = e1;
                    edgesLCP[1][ nb_edges] = e2;
                    NEs[nb_edges]=x;
                    
                    nb_edges ++;
                }
                
                vector<int> curNeighbors;
                GenerateNeighborhood(x, curNeighbors);
                
                for(int k=0; k<curNeighbors.size(); k++)
                {
                    y = curNeighbors[k];
                    if ((indic_P[y]==false) && (weights[y] == wmax))
                    {
                        indic_P[y]=true;
                        LifoPush(LIFO, y);
                        LifoPush(LCP, y);
                        indic_E[y]= true;
                    }
                    
                }
//                for (k = 1; k <= nb_neighbor_edges; k++)
//                {
//                    if (ds>1)
//                        y = neighbor_edge_3D(e1, e2, x, k, rs, cs, ds);
//                    else y = neighbor_edge(x, k, rs, cs, ds);
//                    if (y != -1)
//                        if ((indic_P[y]==false) && (weights[y] == wmax))
//                        {
//                            indic_P[y]=true;
//                            LifoPush(LIFO, y);
//                            LifoPush(LCP, y);
//                            indic_E[y]= true;
//                        }
//                }
            }
            for (j=0;j<nb_vertices;j++)
                indic_VP[LCVP[j]]=0;
            for (j=0;j<LCP->Sp;j++)
                indic_P[LCP->Pts[j]]=false;
            
            // 3. If e_max belongs to a plateau
            if (nb_edges > 0)
            {
                // 4. Evaluate if there are differents seeds on the plateau
                p=0;
                different_seeds = false;
                
                for (i=0;i<nb_labels-1;i++)
                {
                    val = -0.5;
                    for (j=0;j<nb_vertices;j++)
                    {
                        
                        x = LCVP[j];
                        xr = element_find(x, Fth);
                        if(fabs(proba[i][xr]-val)>epsilon_WP && proba[i][xr]>=0 )
                        {
                            p++; val = proba[i][xr];
                        }
                    }
                    if (p>=2)
                    {
                        different_seeds = true;
                        break;
                    }
                    else p=0;
                }
                
                if (different_seeds == true)
                {
                    // 5. Sort the edges of the plateau according to their normal weight
                    for(k=0;k<nb_edges;k++)
                        sorted_weights[k]=normal_weights[NEs[k]];
                    
                    
                    
                    if (quicksort == true)
                        BucketSort(sorted_weights, NEs, nb_edges , max_weight+1);
                    else
                        TriRapideStochastique_dec(sorted_weights,NEs, 0,nb_edges-1);
                    
                    
                    // Merge nodes for edges of real max weight
                    nb_vertices=0;
                    Nnb_edges = 0;
                    for(Ncpt_aretes = 0; Ncpt_aretes< nb_edges; Ncpt_aretes++)
                    {
                        Ne_max=NEs[Ncpt_aretes];
                        e1 = edges[0][ Ne_max];
                        e2 = edges[1][ Ne_max];
                        if (normal_weights[Ne_max] != wmax)
                            merge_node (e1, e2,  Rnk, Fth, proba, nb_labels);
                        else
                        {
                            re1 = element_find(e1, Fth );
                            re2 = element_find(e2, Fth );
                            if ((re1 !=re2)&&((proba[0][re1]<0 || proba[0][re2]<0)))
                            {
                                if (indic_VP[re1]==0)
                                {
                                    LCVP[nb_vertices]=re1;
                                    nb_vertices++;
                                    indic_VP[re1]=1;
                                }
                                if (indic_VP[re2]==0)
                                {
                                    LCVP[nb_vertices]=re2;
                                    nb_vertices++;
                                    indic_VP[re2]=1;
                                }
                                edgesLCP[0][ Nnb_edges] = re1;
                                edgesLCP[1][ Nnb_edges] = re2;
                                Nnb_edges ++;
                            }
                        }
                    }
                    for (i=0;i<nb_labels-1;i++)
                    {
                        k=0;
                        for (j=0;j<nb_vertices;j++)
                        {
                            xr = LCVP[j];
                            if (proba[i][xr]>=0)
                            {
                                local_labels[i][k] = proba[i][xr];
                                local_seeds[k] = xr;
                                k++;
                            }
                        }
                    }
                    
                    // 6. Execute Random Walker on plateaus
                    
                    if(nb_vertices<SIZE_MAX_PLATEAU)
                        success = RandomWalker(edgesLCP, Nnb_edges, LCVP, indic_VP, nb_vertices, local_seeds, local_labels, k, nb_labels, proba);
                    if ((nb_vertices>=SIZE_MAX_PLATEAU)||(success==false))
                    {
                        printf("Plateau of a big size (%d vertices,%d edges) the RW is not performed on it\n", nb_vertices, Nnb_edges);
                        for (j=0;j<Nnb_edges;j++)
                        {
                            e1 = edgesLCP[0][j];
                            e2 = edgesLCP[1][j];
                            merge_node (e1, e2,  Rnk, Fth, proba, nb_labels);
                        }
                    }
                    
                    for (j=0;j<nb_vertices;j++)
                        indic_VP[LCVP[j]]=0;
                }
                else // if different seeds = false
                    // 7. Merge nodes for edges of max weight
                {
                    for (j=0;j<nb_edges;j++)
                    {
                        e1 = edgesLCP[0][j];
                        e2 = edgesLCP[1][j];
                        merge_node (e1, e2,  Rnk, Fth, proba, nb_labels);
                    }
                }
            }
            LifoFlush(LCP);
        } // end main loop
        
        //building the final proba map (find the root vertex of each tree)
        for (i=0; i<N; i++)
        {
            j=i;
            xr = i;
            while(Fth[i] != i)
            {
                i = xr;
                xr = Fth[i];
            }
            for(k=0; k< nb_labels-1;k++) proba[k][j] =proba[k][i];
            i=j;
        }
        
        //writing results
        DeepCopy(this->inputImage, img_proba);
        img_proba->FillBuffer(0);
        double maxi;
        itk_3d_double::IndexType coord;
        for (j = 2; j < N; j++)//start from img idxes.
        {
            maxi=0; argmax = 0; val =1;
            for(k=0; k< nb_labels-1;k++)
            {
                if(proba[k][j]> maxi)
                {
                    maxi = proba[k][j] ;
                    argmax = k;
                }
                val = val - proba[k][j];
            }
            if (val>maxi) argmax = k;
            double value = argmax / (nbLabels - 1);
            
            int x,y,z;
            Coordinate(j, x, y, z);
            coord[0] = x;
            coord[1] = y;
            coord[2] = z;
            img_proba->SetPixel(coord, value);
        }
        
        
        
//        struct xvimage * temp;
//        int rs = this->imgDim[0];
//        int cs = this->imgDim[1];
//        int ds = this->imgDim[2];
//        temp = allocimage(NULL, rs, cs, ds, VFF_TYP_1_BYTE);
//        
//        unsigned char *Temp = UCHARDATA(temp);
//#ifndef SPEEDY
//        unsigned char *Temp2 = UCHARDATA(img_proba);
//#endif     
//        double maxi;
//        for (j = 0; j < N; j++)
//        {
//            maxi=0; argmax = 0; val =1;
//            for(k=0; k< nb_labels-1;k++)
//            {
//                if(proba[k][j]> maxi) 
//                { 
//                    maxi = proba[k][j] ;
//                    argmax = k;
//                }
//                val = val - proba[k][j];
//                
//            }  
//            if (val>maxi) argmax = k;
//            Temp[j] = ((argmax)*255)/(nb_labels-1);
//        } 
//        
//#ifndef SPEEDY
//        for (j = 0; j < N; j++)
//            Temp2[j] = (unsigned char)(255-255*proba[0][j]); 
//#endif     
        
        // free memory 
        LifoTermine(LCP);
        LifoTermine(LIFO);
        
        for (i=0;i<2;i++) 
            free(edges[i]); free(edges);
        
        for (i=0;i<2;i++) 
            free(edgesLCP[i]); free(edgesLCP);
        
        free(Rnk);
        free(local_seeds);
        for (i=0; i<nb_labels-1; i++)  free(local_labels[i]);
        free(local_labels);
        
        free(LCVP);
        free(Es);
        free(NEs);
        free(indic_E);
        free(indic_VP);
        free(indic_P);
        free(Fth);
        for (i=0; i<nb_labels-1; i++)  free(proba[i]);
        free(proba);
        free(sorted_weights);
        
        //return temp;
    }
    
    
    
    
    
    
    
    
    ////////////////////////////////////////////////////////////////////////
    
    
    double PWSegmentation::GetPixelValue(int vertIdx)
    {
        CHECK_GE(vertIdx, 2);
        CHECK_LT(vertIdx, this->totalVoxelNum + 2);
        itk_3d_double::IndexType idx;
        int nx, ny, nz;
        Coordinate(vertIdx, nx, ny, nz);
        idx[0] = nx; idx[1] = ny; idx[2] = nz;
        double pix = this->inputImage->GetPixel(idx);
        return pix;
    }
    
    double PWSegmentation::UnarySource(int vertIdx)
    {
        CHECK_GE(vertIdx, 2);
        CHECK_LT(vertIdx, this->totalVoxelNum + 2);
        
        itk_3d_double::IndexType idx;
        int nx, ny, nz;
        Coordinate(vertIdx, nx, ny, nz);
        idx[0] = nx; idx[1] = ny; idx[2] = nz;
        double pix = this->inputImage->GetPixel(idx);
        
        //CV model, modification here.
        double magFactor = 100000 / (255.0*255.0);
        double C1 = 50;
        return std::pow(pix-C1, 2)*magFactor;
    }
    
    double PWSegmentation::UnarySink(int vertIdx)
    {
        CHECK_GE(vertIdx, 2);
        CHECK_LT(vertIdx, this->totalVoxelNum + 2);
        
        itk_3d_double::IndexType idx;
        int nx, ny, nz;
        Coordinate(vertIdx, nx, ny, nz);
        idx[0] = nx; idx[1] = ny; idx[2] = nz;
        double pix = this->inputImage->GetPixel(idx);
        
        //CV model, modification here.
        double magFactor = 100000 / (255.0*255.0);
        double C2 = 200;
        return std::pow(pix-C2, 2)*magFactor;
    }
    
    void PWSegmentation::GenerateWeight()
    {
        double beta = 1.0 / (255.0*255);
        double magFactor = 100000;
        
        //for test, we use CV model, thus we have C1 and C2 avaliable.
        //double C1 = 200;
        //double C2 = 50;
        
        //first, unary weight.
        this->weights = new WEIGHT_TYPE[this->edgeNum];
        for (int i=0; i<this->totalVoxelNum; i++) {
            int vertIdx = edges[1][i]; //(0,i)=0是source.
            this->weights[i] = UnarySource(vertIdx);
        }
        for (int i=totalVoxelNum; i<2*totalVoxelNum; i++) {
            int vertIdx = edges[1][i]; //(0,i)=1是sink.
            this->weights[i] = UnarySource(vertIdx);
        }
        
        //next, binary weight.
        for (int i=2*totalVoxelNum; i<this->edgeNum; i++) {
            int e1 = edges[0][i];
            int e2 = edges[1][i];
            
            double p1Val = GetPixelValue(e1);
            double p2Val = GetPixelValue(e2);
            
            this->weights[i] = exp(-1.0 * beta * std::pow(p1Val-p2Val, 2)) * magFactor;
        }
        
        this->max_weight = 1;
        for (int i=0; i<this->edgeNum; i++) {
            if(this->weights[i] > max_weight)
            {
                this->max_weight = this->weights[i];
            }
        }
    }
    
    void PWSegmentation::GetNeighborhood(int vertIdx, int forbidEdgeIdx, vector<int> &neighborhoods)
    {
        
        CHECK_LT(vertIdx, this->totalVoxelNum+2);
        CHECK_GE(vertIdx, 0);
        int e1 = this->edges[0][forbidEdgeIdx];
        int e2 = this->edges[1][forbidEdgeIdx];
        CHECK( (e1 == vertIdx || e2 == vertIdx) );
        
        neighborhoods.clear();
        
        if(vertIdx == 0)
        {
            //if source, add everything except forbidEdgeIdxes;
            for (int i=0; i<totalVoxelNum; i++) {
                if( i != forbidEdgeIdx)
                    neighborhoods.push_back(i);
            }
            return;
        }
        
        if(vertIdx == 1)
        {
            //if sink, add everything except forbidEdgeIdxes;
            for (int i=totalVoxelNum; i<2*totalVoxelNum; i++) {
                if (i != forbidEdgeIdx)
                {
                    neighborhoods.push_back(i);
                }
            }
            return;
        }
        
        //ordinary node.
        int nx, ny, nz;
        Coordinate(vertIdx, nx, ny, nz);
        
        CHECK_GE(nx, 0);
        CHECK_GE(ny, 0);
        CHECK_GE(nz, 0);
        CHECK_LT(nx, imgDim[0]);
        CHECK_LT(ny, imgDim[1]);
        CHECK_LT(nz, imgDim[2]);
        
        //ordinary links.
        if(nx > 0 && Idx(nx-1, ny, nz) != forbidEdgeIdx) neighborhoods.push_back(Idx(nx-1, ny, nz));
        if(nx < imgDim[0] - 1 && Idx(nx+1, ny, nz) != forbidEdgeIdx) neighborhoods.push_back(Idx(nx+1, ny, nz));
        
        if(ny > 0 && Idx(nx, ny-1, nz) != forbidEdgeIdx) neighborhoods.push_back(Idx(nx, ny-1, nz));
        if(ny < imgDim[1] - 1 && Idx(nx, ny+1, nz) != forbidEdgeIdx) neighborhoods.push_back(Idx(nx, ny+1, nz));
        
        if(nz > 0 && Idx(nx, ny, nz-1) != forbidEdgeIdx) neighborhoods.push_back(Idx(nx, ny, nz-1));
        if(nz < imgDim[2] - 1 && Idx(nx, ny, nz+1) != forbidEdgeIdx) neighborhoods.push_back(Idx(nx, ny, nz+1));

        //source/sink links.
        if(vertIdx - 2 != forbidEdgeIdx)
        {
            neighborhoods.push_back(vertIdx-2);
        }
        
        if(vertIdx - 2 + totalVoxelNum != forbidEdgeIdx)
        {
            neighborhoods.push_back(vertIdx - 2 + totalVoxelNum);
        }
    }
    
    void PWSegmentation::GenerateNeighborhood(int edgeIdx, vector<int> &neighborhood)
    {
        CHECK_LT(edgeIdx, this->edgeNum);
        CHECK_GE(edgeIdx, 0);
        
        int e1 = this->edges[0][edgeIdx];
        int e2 = this->edges[1][edgeIdx];
        
        neighborhood.clear();
        
        vector<int> neiE1;
        GetNeighborhood(e1, edgeIdx, neiE1);
        vector<int> neiE2;
        GetNeighborhood(e2, edgeIdx, neiE2);
        
        neighborhood.insert(neighborhood.end(), neiE1.begin(), neiE1.end());
        neighborhood.insert(neighborhood.end(), neiE2.begin(), neiE2.end());
    }
    
    PWSegmentation::~PWSegmentation()
    {
//        delete [] seed;
//        for (int i=0; i<2; i++) {
//            delete [] edges[i];
//        }
        //delete [] edges;
    }
    
    int PWSegmentation::Idx(int x, int y, int z)
    {
        CHECK_GE(x, 0);
        CHECK_GE(y, 0);
        CHECK_GE(z, 0);
        CHECK_LT(x, imgDim[0]);
        CHECK_LT(y, imgDim[1]);
        CHECK_LT(z, imgDim[2]);
        
        return x + y * imgDim[0] + z * imgDim[1] * imgDim[0] + 2;
    }
    
    void PWSegmentation::Coordinate(int Idx, int &x, int &y, int &z)
    {
        CHECK_GE(Idx, 2);
        CHECK_LT(Idx, 2+this->totalVoxelNum);
        int actualIdx = Idx - 2;// = x + y * imgDim[2] + z * imgDim[1] * imgDim[2]
        int dimXY = imgDim[0] * imgDim[1];
        z = actualIdx / dimXY;
        y = (actualIdx - z * dimXY) / imgDim[0];
        x = actualIdx - z * dimXY - y * imgDim[0];
        
        CHECK_GE(x, 0);
        CHECK_GE(y, 0);
        CHECK_GE(z, 0);
        CHECK_LT(x, imgDim[0]);
        CHECK_LT(y, imgDim[1]);
        CHECK_LT(z, imgDim[2]);
    }
    
    void PWSegmentation::BuildSeed()
    {
        this->seed_size = 2+this->fgSeeds.size() + this->bgSeeds.size();
        this->seed = new int[seed_size];
        this->labels = new unsigned char[seed_size];
        
        //first hepatic.
        int nCounter = 0;
        seed[nCounter] = 0;//Source.
        labels[nCounter++] = 1;
        for (int i=0; i<fgSeeds.size(); i++) {
            int fgIdx = Idx(fgSeeds[i].x, fgSeeds[i].y, fgSeeds[i].z);
            seed[nCounter] = fgIdx;
            labels[nCounter++] = 1;
        }
        
        seed[nCounter] = 1;//Sink.
        labels[nCounter++] = 2;
        for (int i=0; i<bgSeeds.size(); i++) {
            int bgIdx = Idx(bgSeeds[i].x, bgSeeds[i].y, bgSeeds[i].z);
            seed[nCounter] = bgIdx;
            labels[nCounter++] = 2;
        }
        
        CHECK_EQ(nCounter, this->seed_size);
    }
    
    PWSegmentation::PWSegmentation(itk_3d_double::Pointer& _inputImage,
                                   vector<cv::Point3i> _fgSeeds,
                                   vector<cv::Point3i> _bgSeeds)
    {
        this->inputImage = itk_3d_double::New();
        DeepCopy(_inputImage, inputImage);
        
        this->fgSeeds = fgSeeds;
        this->bgSeeds = bgSeeds;
        itk_3d_double::SizeType sz = _inputImage->GetLargestPossibleRegion().GetSize();
        imgDim[0] = sz[0];
        imgDim[1] = sz[1];
        imgDim[2] = sz[2];
        this->nbLabels = 2;
        this->totalVoxelNum = sz[0]*sz[1]*sz[2];
        
        //debug check.
        CHECK_GE(imgDim[0], 1);
        CHECK_GE(imgDim[1], 1);
        CHECK_GE(imgDim[2], 1);
        
        //build graph.
        BuildVertexesAndEdges();
        
        //build seed
        BuildSeed();
        
        //build weight
        GenerateWeight();
        
        std::cout<<"max weight = "<<this->max_weight<<std::endl;
    }
    
    void PWSegmentation::BuildVertexesAndEdges()
    {
        this->edgeNum = (imgDim[0]-1)*imgDim[1]*imgDim[2]+ imgDim[0]*(imgDim[1]-1)*imgDim[2] + imgDim[0]*imgDim[1]*(imgDim[2] -1)+ 2*totalVoxelNum; // edges inside image plus two virtual points connecting all voxels inside image.
        this->edgeXOffset = (imgDim[0]-1)*imgDim[1]*imgDim[2];
        this->edgeYOffset = imgDim[0]*(imgDim[1]-1)*imgDim[2];
        
        edges = new int*[2];
        for (int i=0; i<2; i++) {
            edges[i] = new int[this->edgeNum];
        }
        
        int nCounter = 0;
        //allocate edges.
        //first links with the source
        for (int z=0; z<imgDim[2]; z++) {
            for (int y=0; y<imgDim[1]; y++) {
                for (int x=0; x<imgDim[0]; x++) {
                    edges[0][nCounter] = 0;//Source
                    edges[1][nCounter] = Idx(x, y, z);
                    nCounter++;
                }
            }
        }
        
        //next links with the sink.
        for (int z=0; z<imgDim[2]; z++) {
            for (int y=0; y<imgDim[1]; y++) {
                for (int x=0; x<imgDim[0]; x++) {
                    edges[0][nCounter] = 1;//Source
                    edges[1][nCounter] = Idx(x, y, z);
                    nCounter++;
                }
            }
        }
        
        //finally directional links.
        //first x
        for (int z=0; z<imgDim[2]; z++) {
            for (int y=0; y<imgDim[1]; y++) {
                for (int x=0; x<imgDim[0]-1; x++) {
                    edges[0][nCounter] = Idx(x, y, z);
                    edges[1][nCounter] = Idx(x+1, y, z);
                    nCounter++;
                }
            }
        }
        
        //then y.
        for (int z=0; z<imgDim[2]; z++) {
            for (int y=0; y<imgDim[1]-1; y++) {
                for (int x=0; x<imgDim[0]; x++) {
                    edges[0][nCounter] = Idx(x, y, z);
                    edges[1][nCounter] = Idx(x, y+1, z);
                    nCounter++;
                }
            }
        }
        //finally z.
        for (int z=0; z<imgDim[2]-1; z++) {
            for (int y=0; y<imgDim[1]; y++) {
                for (int x=0; x<imgDim[0]; x++) {
                    edges[0][nCounter] = Idx(x, y, z);
                    edges[1][nCounter] = Idx(x, y, z+1);
                    nCounter++;
                }
            }
        }

        CHECK_EQ(nCounter, edgeNum);
    }
    
    void PWSegmentation::DeepCopy(itk_3d_double::Pointer &input, itk_3d_double::Pointer &output)
    {
        output->SetRegions(input->GetLargestPossibleRegion());
        output->Allocate();
        output->SetSpacing(input->GetSpacing());
        output->SetOrigin(input->GetOrigin());
        
        itk::ImageRegionConstIterator<itk_3d_double> inputIterator(input, input->GetLargestPossibleRegion());
        itk::ImageRegionIterator<itk_3d_double> outputIterator(output, output->GetLargestPossibleRegion());
        while(!inputIterator.IsAtEnd())
        {
            outputIterator.Set(inputIterator.Get());
            ++inputIterator;
            ++outputIterator;
        }
    }
}