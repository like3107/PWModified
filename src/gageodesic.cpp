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
#include <stdint.h>
#include <stdlib.h>

#include "../include/YZXHeader.h"

#include "../PINK/mccodimage.h"
#include "../PINK/mclifo.h"
#include "../PINK/mcutil.h"
#include "../include/gageodesic.h"
#include "../include/ccsort.h"
#include "../include/union_find.h"
#include "../include/cccodimage.h"
#include "../include/powerwatsegm.h"

uint32_t MAX;

/*================================================*/
void element_link_geod_dilate( int n,
			       int p,
			       int *Fth, 
			       uint32_t *G,
			       uint32_t  *O)
/*================================================*/  
{ 
  int r = element_find(n, Fth);
  
    //std::cout<<"r = "<<r<<" n = "<<n<<std::endl;
    
  if (r != p)
    {
      if((G[r] == G[p])||(G[p]>=O[r]))
	{
	  Fth[r] = p;
	  O[p] = mcmax(O[r],O[p]);
	}
      else O[p] = MAX;
    } 
}


/* ==================================================================================================== */
/**
 gageodilate_union_find(seeds_function, normal_weights, weights, edges, rs, cs, ds, 255, quicksort);
 
 normal_weights[e_ij] = 255 - |I_i - I_j|
 
 uint32_t * seeds_function =(uint32_t *)calloc(M,sizeof(uint32_t));
 注意这边这个seeds_function在非seed点附近都是0，在seed点的6邻域是normal_weights，就是255- |\nabla I|
 
 */

void gageodilate_union_find( uint32_t  *F,  /* f : image seeds */ 
			     uint32_t  *G,  /* g : image weights */
			     uint32_t * O,  /* O : result of the reconstruction by dilation of g under f */
			     int ** edges,  /* list of couple of vertices forming edges*/
			     int32_t rs,    /* row size */
			     int32_t cs,    /* col size */
			     int32_t ds,    /* depth size */
			     int max_weight, 
			     bool quicksort)  
/* ===================================================================================================== */
/* reconstruction by dilation of g under f.  Union-find method described by Luc Vicent. */
{
 
    int k, p,i,n;
    int M = ds*rs*(cs-1)+ds*(rs-1)*cs+(ds-1)*cs*rs;      /*number of edges*/
    bool * Mrk = (bool*)calloc(M,sizeof(bool));
    int * Fth = (int*)malloc(M*sizeof(int));
    
    // Es : E sorted by decreasing weights
    int * Es = (int*)malloc(M*sizeof(int));
    
    for(k=0;k<M;k++)
    {
        Fth[k]=k;
        O[k]=F[k];
        F[k]=G[k];
        Es[k] = k;
    }
    
    //Fth, Es似乎都是index的矩阵
    //O就是output一开始初始化成seed周围的一些梯度值。而F就是Seed image被改成了初始的image weight?
    //就是255-梯度绝对值那个。
    
    cv::Mat tmpdata;
    Data2CVMat(O, cs, rs, tmpdata);
    cv::imshow("O img", tmpdata);
    
    
    //将F从小到大排序（就是把初始的weight image从小到大排序）
    MAX = max_weight;
    if(quicksort) BucketSortCroiss (F, Es, M, MAX+1);
    else TriRapideStochastique_inc(F, Es, 0, M-1);
    
//    for (int i=0;i<M-1;i++)
//    {
//        std::cout<<Es[i]<<",";
//    }
//    std::cout<<std::endl;
    
//    for (int i=0; i<M; i++)
//    {
//        std::cout<<Mrk[i]<<",";
//    }
//    std::cout<<std::endl;
    
    /* first pass */
    if(ds==1)//2D
    {
        //从大到小遍历，这个值越大说明当前点的梯度就越小。
        for(k=M-1;k>=0;k--)
        {
            p = Es[k];
            for (i = 1; i <= 6; i += 1) // The 6 neighbors
            {
                n = neighbor_edge(p, i, rs, cs, ds);
                
                //p是当前点，n是当前的edge, Fth是index, 而且Fth(x) = x;
                if (n != -1 && Mrk[n])
                {
                    element_link_geod_dilate(n,p, Fth, G, O);
                }
                
                Mrk[p]=true;
            }
        }
    }
    else // 3D
    {
        for(k=M-1;k>=0;k--)
        {
            p = Es[k];
            for (i = 1; i <= 12; i += 1) // The 12 neighbors
            {
                n = neighbor_edge_3D(edges[0][p], edges[1][p], p, i, rs, cs, ds);
                if (n != -1)
                    if(Mrk[n])
                        element_link_geod_dilate(n,p, Fth, G, O);
                Mrk[p]=true;
            }
        }
    }
    
    /* second pass */
    
    for(k=0;k<M;k++)
    {
        p = Es[k];
        if (Fth[p]==p) // p is root
        {
            if (O[p]==MAX) O[p]=G[p];
        }
        else O[p]= O[Fth[p]];
    }
    
    free(Es);
    free(Mrk);
    free(Fth);
}

