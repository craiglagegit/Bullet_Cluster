/*
  ------------------------------------------------------------------------------
  Author: Craig Lage, NYU CCPP
  Date: Feb 6, 2014

*/

//****************** triaxial.cpp **************************

#include "triaxial.h"

//***************MAIN PROGRAM***********************

int main(int argc, char *argv[])
{
  if(argc != 2)
    {
      printf("\nwrong number of arguments\n");
      printf("Only argument should be name of .cfg file");
      exit(0);
    }

  Collision* Coll = new Collision(argv[1]);

  delete Coll;
  return 0;
}
//***************END MAIN PROGRAM***********************

