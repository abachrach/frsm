#include <stdio.h>
#include <stdlib.h>
#include <frsm/frsm.hpp>

using namespace std;
using namespace frsm;

void runTest()
{
  //hardcoded scan matcher params
  double metersPerPixel = .02; //translational resolution for the brute force search
  double thetaResolution = .01; //angular step size for the brute force search
  frsm_incremental_matching_modes_t matchingMode = FRSM_GRID_COORD; //use gradient descent to improve estimate after brute force search
  int useMultires = 3; // low resolution will have resolution metersPerPixel * 2^useMultiRes

  int useThreads = 1;

  //create the actual scan matcher object
  ScanMatcher * sm = new ScanMatcher(metersPerPixel, thetaResolution, useMultires, useThreads, true);

  if (sm->isUsingIPP())
    fprintf(stderr, "Using IPP\n");
  else
    fprintf(stderr, "NOT using IPP\n");

  //TODO: load two fake scans, and match them... this won't run
  ScanTransform T;
  sm->addScan(NULL, 0, &T, FRSM_HOKUYO_UTM, 0, true);
  sm->gridMatch(NULL, 0, &T, .5, .5, .1);



}

int main(int argc, char *argv[])
{

  runTest();

  return 0;
}

