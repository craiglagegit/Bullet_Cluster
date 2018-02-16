

#include "allvars.h"


 int     NumPart;



 double TreeAllocFactor;

 int *Head, *Next, *Tail, *Len;

 int *GroupLen, *GroupOffset;

 FILE *Logfile;


 double BoxSize, BoxHalf;

 int    DesDensNgb;

 int Axis1, Axis2, Axis3;

 float *Hsml;
 float *Quantity1, *Quantity2, *Quantity3, *Quantity4;
 float *Mass;
 float *Rho;

 float Xmin, Ymin, Xmax, Ymax, Zmin, Zmax, Hmax;

 float Xc, Yc, Zc;

 int Xpixels, Ypixels, Zpixels;

 float *Value, *ValueQuantity1, *ValueQuantity2, *ValueQuantity3, *ValueQuantity4;

 float Softening;


 struct particle_data
*P;                            /*!< points to particles on this processor */





 struct r2data
*R2list;






 int    AnzNodes;
 int    MaxNodes;



 int *Nextnode;
 int *Father;

 struct NODE
*Nodes_base,                    /*!< points to the actual memory allocted for the nodes */
*Nodes;                         /*!< this is a pointer used to access the nodes which is shifted such that Nodes[All.MaxPart] 
				  gives the first allocated node */











