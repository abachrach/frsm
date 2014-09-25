#include "ScanMatcher.hpp"
#include <stdio.h>
#include <unistd.h>
#include "RasterLookupTable.hpp"
#include "Contour.hpp"
#include "ScanMatchingUtils.hpp"
#include <signal.h>
#include <assert.h>

using namespace std;
namespace frsm {

void * ScanMatcher_thread_wrapper_func(void * user)
{
  ScanMatcher*parent = (ScanMatcher*) user;
  fprintf(stderr, "starting rebuild thread\n");
  parent->rebuildThreadFunc();
  return NULL;
}

ScanMatcher::ScanMatcher(double metersPerPixel_, double thetaResolution_, int useMultires_, bool useThreads_,
    bool verbose_) :

    metersPerPixel(metersPerPixel_), thetaResolution(thetaResolution_), useMultiRes(useMultires_), verbose(verbose_), contour_extractor(
        NULL), rlt(NULL), rltTmp(NULL), rlt_low_res(NULL), rltTmp_low_res(NULL), useThreads(useThreads_), killThread(1)
{

  hitThresh = 100;

  downsampleFactor = (1 << useMultiRes);

  //line drawing kernel stuff
  double sigma = .0675 / metersPerPixel; //sigma is in pixels
  //  sigma = 21;
  draw_kernels = new LutKernel(sigma);

  lutSq_first_zero = -1;
  lutSq = NULL;

  //set reasonable defaults for successive matching
  memset(&currentPose, 0, sizeof(currentPose));
  memset(&prevPose, 0, sizeof(prevPose));
  maxNumScans = 30;
  initialSearchRangeXY = .13;
  maxSearchRangeXY = .25;
  initialSearchRangeTheta = .1;
  maxSearchRangeTheta = .2;
  matchingMode = FRSM_GRID_COORD;
  addScanHitThresh = .8;

  //threading stuff
  if (useThreads) {
    //remember to make sure that frsm_tictoc gets initialized
    killThread = 0;

    /* Initialize mutex and condition variable objects */
    scans_mutex = g_mutex_new();
    rlt_mutex = g_mutex_new();
    toBeProcessed_mutex = g_mutex_new();
    toBeProcessed_cv = g_cond_new();

    //create rebuilder thread
    rebuilder = g_thread_create(ScanMatcher_thread_wrapper_func, (void *) this, 1, NULL);
  }
}

void ScanMatcher::initSuccessiveMatchingParams(unsigned int maxNumScans_, double initialSearchRangeXY_,
    double maxSearchRangeXY_, double initialSearchRangeTheta_, double maxSearchRangeTheta_,
    frsm_incremental_matching_modes_t matchingMode_, double addScanHitThresh_, bool stationaryMotionModel_,
    double motionModelPriorWeight_, ScanTransform * startPose)
{

  maxNumScans = maxNumScans_;
  initialSearchRangeXY = initialSearchRangeXY_;
  maxSearchRangeXY = maxSearchRangeXY_;
  initialSearchRangeTheta = initialSearchRangeTheta_;
  maxSearchRangeTheta = maxSearchRangeTheta_;
  matchingMode = matchingMode_;
  addScanHitThresh = addScanHitThresh_;
  stationaryMotionModel = stationaryMotionModel_;
  motionModelPriorWeight = motionModelPriorWeight_;

  if (maxSearchRangeXY < initialSearchRangeXY) {
    fprintf(stderr, "WARNING: maxSearchRangeXY of %f is less than initialSearchRangeXY of %f, setting max=initial\n",
        maxSearchRangeXY, initialSearchRangeXY);
    maxSearchRangeXY = initialSearchRangeXY;
  }
  if (maxSearchRangeTheta < initialSearchRangeTheta) {
    fprintf(stderr,
        "WARNING: maxSearchRangeTheta of %f is less than initialSearchRangeTheta of %f, setting max=initial\n",
        maxSearchRangeTheta, initialSearchRangeTheta);
    maxSearchRangeTheta = initialSearchRangeTheta;
  }
  if (startPose != NULL) {
    memcpy(&currentPose, startPose, sizeof(currentPose));
    memcpy(&prevPose, startPose, sizeof(prevPose));
  }
  else {
    memset(&currentPose, 0, sizeof(currentPose));
    memset(&prevPose, 0, sizeof(prevPose));
  }

}

ScanMatcher::~ScanMatcher()
{
  clearScans(true);

  //TODO: delete kernels/LUT

  if (useThreads) {
    //aquire all locks so we can destroy them
    g_mutex_lock(scans_mutex);
    g_mutex_lock(rlt_mutex);
    //kill the rebuilder thread
    while (killThread != -1) {
      killThread = 1;
      g_cond_broadcast(toBeProcessed_cv);
      usleep(10000);

    }

    // destroy mutex and condition variable objects
    g_mutex_free(scans_mutex);
    g_mutex_free(rlt_mutex);
    g_mutex_free(toBeProcessed_mutex);
    g_cond_free(toBeProcessed_cv);

  }

  if (lutSq != NULL
  )
    free(lutSq);
  if (draw_kernels != NULL
  )
    delete draw_kernels;

}

void ScanMatcher::clearScans(bool deleteScans)
{
  if (useThreads) {
    g_mutex_lock(scans_mutex);
    g_mutex_lock(rlt_mutex);
    g_mutex_lock(toBeProcessed_mutex);
  }
  if (deleteScans) {
    while (scans.size() > 1) {
      delete scans.front();
      scans.pop_front();
    }
    while (scansToBeProcessed.size() > 1) {
      delete scansToBeProcessed.front();
      scansToBeProcessed.pop_front();
    }
  }
  else {
    scans.clear();
    scansToBeProcessed.clear();
  }
  //rlt is no longer valid
  if (rlt != NULL
  )
    delete rlt;
  if (rltTmp != NULL
  )
    delete rltTmp;

  if (rlt_low_res != NULL
  )
    delete rlt_low_res;
  if (rltTmp_low_res != NULL
  )
    delete rltTmp_low_res;

  rlt = NULL;
  rltTmp = NULL;
  rlt_low_res = NULL;
  rltTmp_low_res = NULL;

  if (useThreads) {
    g_mutex_unlock(scans_mutex);
    g_mutex_unlock(rlt_mutex);
    g_mutex_unlock(toBeProcessed_mutex);
  }
}

void ScanMatcher::rebuildThreadFunc()
{
  if (!useThreads) {
    fprintf(stderr, "ERROR, threading is disabled! this shouldn't be being called\n");
    exit(1);
  }

  while (true) {
    g_mutex_lock(toBeProcessed_mutex);
    while (scansToBeProcessed.empty()) {
      g_cond_wait(toBeProcessed_cv, toBeProcessed_mutex);
      if (killThread) {
        g_mutex_unlock(toBeProcessed_mutex);
        killThread = -1;
        fprintf(stderr, "rebuild thread exiting\n");
        return;
      }
    }
    std::list<Scan *> scansBeingProcessed;
    //there are scans to be processed... lets take care of business
    //first swap out the lists so we can let go of the lock...
    scansBeingProcessed = scansToBeProcessed;
    scansToBeProcessed.clear();
    g_mutex_unlock(toBeProcessed_mutex);

    if (scansBeingProcessed.size() > 1) {
      if (verbose) {
        fprintf(stderr, "there are %zu scans to be processed, discarding all but most recent\n",
            scansBeingProcessed.size());
      }
      while (scansBeingProcessed.size() > 1) {
        delete scansBeingProcessed.front();
        scansBeingProcessed.pop_front();
      }
    }

    if (scansBeingProcessed.size() > maxNumScans) {
      fprintf(stderr, "Scan Matcher is way behind!, there are %zu scans to be added\n", scansBeingProcessed.size());
      while (scansBeingProcessed.size() > maxNumScans) {
        delete scansBeingProcessed.front();
        scansBeingProcessed.pop_front();
      }
    }

    g_mutex_lock(scans_mutex);
    //make space for new scans
    while (numScans() + scansBeingProcessed.size() > maxNumScans) {
      delete scans.front();
      scans.pop_front();
    }

    list<Scan *>::iterator it;
    it = scans.end(); //save the current end for drawing the new scans
    while (!scansBeingProcessed.empty()) {
      Scan * s = scansBeingProcessed.front();
      frsm_tictoc("findContours");

      if (contour_extractor != NULL && contour_extractor->laserType != s->laser_type) {
        delete contour_extractor;
        contour_extractor = NULL;
      }
      if (contour_extractor == NULL
      )
        contour_extractor = new ContourExtractor(s->laser_type);
      contour_extractor->findContours(s->ppoints, s->numPoints, s->contours);

      frsm_tictoc("findContours");
      //            s->drawContours(1000,100);

      scans.push_back(s);
      scansBeingProcessed.pop_front();
    }

    frsm_tictoc("rebuildRaster");
    rebuildRaster(&rltTmp); //rebuild the raster in the tmp
    frsm_tictoc("rebuildRaster");

    if (useMultiRes > 0) {
      frsm_tictoc("rebuild_lowRes");
      rltTmp_low_res = new RasterLookupTable(rltTmp, downsampleFactor);
      frsm_tictoc("rebuild_lowRes");
    }

    //swap rltTmp with rlt to put it in use
    g_mutex_lock(rlt_mutex);
    if (cancelAdd) {
      if (verbose)
        fprintf(stderr, "Scan add was canceled!\n");
      delete rltTmp;
      rltTmp = NULL;
      if (useMultiRes > 0) {
        delete rltTmp_low_res;
        rltTmp_low_res = NULL;
      }
    }
    else {
      delete rlt;
      rlt = rltTmp;
      rltTmp = NULL;

      if (useMultiRes > 0) {
        delete rlt_low_res;
        rlt_low_res = rltTmp_low_res;
        rltTmp_low_res = NULL;
      }
      if (verbose)
        fprintf(stderr, "rlt swapped!\n");
    }
    g_mutex_unlock(rlt_mutex);
    g_mutex_unlock(scans_mutex);

    //    //clear out scans that got added in the interim
    //    g_mutex_lock(toBeProcessed_mutex);
    //    scansToBeProcessed.clear();
    //    g_mutex_unlock(toBeProcessed_mutex);

  }

}

void ScanMatcher::computeBounds(double *minx, double *miny, double *maxx, double *maxy)
{
  *minx = DBL_MAX;
  *maxx = -DBL_MAX;
  *miny = DBL_MAX;
  *maxy = -DBL_MAX;

  frsm_tictoc("rebuild_bounds");
  // Compute bounds of the scans.
  list<Scan *>::iterator it;
  for (it = scans.begin(); it != scans.end(); ++it) {
    Scan * s = *it;
    for (unsigned cidx = 0; cidx < s->contours.size(); cidx++) {
      for (unsigned i = 0; i < s->contours[cidx]->points.size(); i++) {
        frsmPoint p = s->contours[cidx]->points[i];

        *minx = fmin(*minx, p.x);
        *maxx = fmax(*maxx, p.x);
        *miny = fmin(*miny, p.y);
        *maxy = fmax(*maxy, p.y);
      }
    }
  }
  frsm_tictoc("rebuild_bounds");

}

void ScanMatcher::rebuildRaster(RasterLookupTable ** rasterTable)
{
  //dispatcher for which rebuildRaster version to use
  rebuildRaster_blurLine(rasterTable);
}

void ScanMatcher::rebuildRaster_olson(RasterLookupTable ** rasterTable)
{
  double minx, maxx;
  double miny, maxy;
  computeBounds(&minx, &miny, &maxx, &maxy);

  // create the LUT, with a bit more space than we currently
  // think we'll need (the thought being that we could add new
  // scans without having to rebuild the whole damn thing
  double margin = fmax(1, 0.1 * fmax(maxx - minx, maxy - miny));

  RasterLookupTable * rt = new RasterLookupTable(minx - margin, miny - margin, maxx + margin, maxy + margin,
      metersPerPixel, downsampleFactor);

  //initialize the lookup table if it hasn't ben already done
  if (lutSq_first_zero < 0 || lutSq == NULL) {
    /**
     * Make a lookup table with 'sz' entries. The lookup table will have an
     * exponential fall-off with maximum value 255.
     *
     * lut[0] = weight*255 ... lut[i] = weight*255*e^{-i*alpha} ... lut[sz] =
     * weight*255*e^{-maxChiSq}
     *
     * Weight should be [0,1].
     */
    int sz = 512;
    double maxChiSq = 10;
    double weight = 1.0;

    lutSq = (uint8_t *) malloc(sz * sizeof(char));
    bool set_lutSq_first_zero = false;
    for (int i = 0; i < sz; i++) {
      int v = (int) round(255.0 * weight * exp(-maxChiSq * i / (sz - 1)));
      if (v == 0 && !set_lutSq_first_zero) {
        set_lutSq_first_zero = true;
        lutSq_first_zero = i;
      }
      assert(v >= 0);
      assert(v <= 255);
      lutSq[i] = (uint8_t) v;
    }
    if (!set_lutSq_first_zero) {
      set_lutSq_first_zero = true;
      lutSq_first_zero = sz;
    }

  }

  double lutSqRange = pow(0.25, 2);

  // draw each scan.
  list<Scan *>::iterator it;
  for (it = scans.begin(); it != scans.end(); ++it) {
    Scan * s = *it;
    for (unsigned cidx = 0; cidx < s->contours.size(); cidx++) {
      for (unsigned i = 0; i + 1 < s->contours[cidx]->points.size(); i++) {
        frsmPoint p0 = s->contours[cidx]->points[i];
        frsmPoint p1 = s->contours[cidx]->points[i + 1];

        double length = frsm_dist(&p0, &p1);

        rt->drawRectangle((p0.x + p1.x) / 2, (p0.y + p1.y) / 2, length, 0, atan2(p1.y - p0.y, p1.x - p0.x), lutSq, 512,
            lutSq_first_zero, lutSqRange);

      }
    }
  }

  //  free(lutSq);

  if (*rasterTable != NULL
  )
    delete *rasterTable;
  *rasterTable = rt;
  return;
}

void ScanMatcher::drawBlurredScan(RasterLookupTable * rt, Scan * s)
{
  for (unsigned cidx = 0; cidx < s->contours.size(); cidx++) {
    if (s->contours[cidx]->points.size() >= 2) {
      for (unsigned i = 0; i < s->contours[cidx]->points.size() - 1; i++) {
        const frsmPoint &p0 = s->contours[cidx]->points[i];
        const frsmPoint &p1 = s->contours[cidx]->points[i + 1];
        rt->drawBlurredPoint(&p0, draw_kernels);
        rt->drawBlurredLine(&p0, &p1, draw_kernels);
      }
    }
    //draw the last point in the contour
    rt->drawBlurredPoint(&s->contours[cidx]->points.back(), draw_kernels);
  }
}

void ScanMatcher::rebuildRaster_blurLine(RasterLookupTable ** rasterTable)
{
  double minx, maxx;
  double miny, maxy;
  computeBounds(&minx, &miny, &maxx, &maxy);

  double margin = fmax(.5, 0.1 * fmax(maxx - minx, maxy - miny));

  frsm_tictoc("rebuild_alloc");
  RasterLookupTable * rt = new RasterLookupTable(minx - margin, miny - margin, maxx + margin, maxy + margin,
      metersPerPixel, downsampleFactor);
  frsm_tictoc("rebuild_alloc");

  frsm_tictoc("drawBlurLines");
  // draw each scan.
  list<Scan *>::iterator it;
  for (it = scans.begin(); it != scans.end(); ++it) {
    Scan * s = *it;
    drawBlurredScan(rt, s);
  }
  frsm_tictoc("drawBlurLines");

  if (*rasterTable != NULL
  )
    delete *rasterTable;
  *rasterTable = rt;
  return;
}

ScanTransform ScanMatcher::gridMatch(frsmPoint * points, unsigned numPoints, ScanTransform * prior, double xRange,
    double yRange, double thetaRange, int * xSat, int *ySat, int * thetaSat)
{
  if (useThreads)
    g_mutex_lock(rlt_mutex);

  //  fprintf(stderr,"matching has rlt lock\n");

  if (rlt == NULL) {
    fprintf(stderr, "ERROR: raster lookup table is NULL\n");
    exit(1);
  }

  if (useMultiRes > 0 && rlt_low_res == NULL) {
    fprintf(stderr, "ERROR: low res raster lookup table is NULL, and multiRes is %d\n", useMultiRes);
    exit(1);
  }

  //////////////////////////////////////////////
  // Here's where we actually do a scan match!
  ScanTransform r;
  if (useMultiRes <= 0) {
    frsm_tictoc("evaluate3D");
    r = rlt->evaluate3D(points, numPoints, prior, xRange, yRange, thetaRange, thetaResolution, hitThresh, xSat, ySat,
        thetaSat);
    frsm_tictoc("evaluate3D");
  }
  else {
    frsm_tictoc("evaluate3D_multiRes");
    r = rlt_low_res->evaluate3D_multiRes(rlt, points, numPoints, prior, xRange, yRange, thetaRange, thetaResolution,
        hitThresh, xSat, ySat, thetaSat);
    frsm_tictoc("evaluate3D_multiRes");
  }
  r.theta = frsm_normalize_theta(r.theta);

  //    fprintf(stderr,"r.x=%f \t r.y=%f \t r.t=%f\t r.score=%f\n", r.x, r.y, r.theta, r.score);
  //  fprintf(stderr,"r1.x=%f \t r1.y=%f \t r1.t=%f\t r1.score=%f\n",r1.x,r1.y,r1.theta,r1.score);
  //  fprintf(stderr,"r2.x=%f \t r2.y=%f \t r2.t=%f\t r2.score=%f\n",r2.x,r2.y,r2.theta,r2.score);

  if (useThreads)
    g_mutex_unlock(rlt_mutex);
  //  fprintf(stderr,"matching released rlt lock\n");
  return r;

}

//TODO: gauss newton or ESM optimizaiton?
ScanTransform ScanMatcher::coordAscentMatch(frsmPoint * points, unsigned numPoints, ScanTransform * startTrans)
{
  typedef enum {
    Front = 0, Back = 1, Right = 2, Left = 3, TurnLeft = 4, TurnRight = 5, Done = 6
  } move_type_t;
  if (useThreads)
    g_mutex_lock(rlt_mutex);

  //halve the step size this many times.
  int numRounds = 5;
  //at each round, try this many steps
  int max_steps = 100;
  if (matchingMode > FRSM_RUN_GRID) {
    max_steps = 1000;
  }
  //start with a step size that is half the brute force search's resolution
  float ldelta = metersPerPixel / 2.0;
  //  float ldelta_diag = sqrt(2) / 2* metersPerPixel / 2.0;
  float adelta = thetaResolution / 2.0;
  //  float adelta_diag = sqrt(2) / 2* thetaResolution / 2.0;

  //if matching mode is Y only, start with the "Right" step
  move_type_t startMove = Front;
  if (matchingMode == FRSM_Y_COORD_ONLY || matchingMode == FRSM_Y_GRID_COORD)
    startMove = Right;

  ScanTransform currentTrans = *startTrans;
  int totalNumSteps = 0; //keep track of total number of tests
  int numSteps = 0;
  for (int i = 0; i < numRounds; i++) {
    bool hadImprovement = true;
    numSteps = 0;
    while (hadImprovement) {
      if (numSteps > max_steps) {
        fprintf(stderr, "gradient ascent didn't converge on round %d\n", i);
        break;
      }
      numSteps++;
      totalNumSteps++;

      hadImprovement = false; //at least one direction should improve
      //try each direction
      for (int move = startMove; move < Done; move++) {
        ScanTransform testTrans = currentTrans;
        switch (move) {
        case Front:
          testTrans.x += ldelta;
          break;
        case Back:
          testTrans.x -= ldelta;
          break;
        case Right:
          testTrans.y += ldelta;
          break;
        case Left:
          testTrans.y -= ldelta;
          break;
        case TurnLeft:
          testTrans.theta += adelta;
          break;
        case TurnRight:
          testTrans.theta -= adelta;
          break;
        default:
          assert(false);
          break;
        }
        testTrans.score = rlt->getScore(points, numPoints, &testTrans);
        if (testTrans.score > currentTrans.score) {
          //we got an improvement, keep it
          currentTrans = testTrans;
          hadImprovement = true;
        }
      }
    }
    //    fprintf(stderr,"didn't have improvment after %d steps, total steps is %d\n",numSteps,totalNumSteps);
    //half the step size again
    ldelta /= 2.0;
    adelta /= 2.0;
    //    ldelta_diag /= 2.0;
    //    adelta_diag /= 2.0;

  }
  //  fprintf(stderr, "totalNumSteps= %d, numSteps at lowest level=%d\n", totalNumSteps, numSteps);
  //  if (fabs(currentTrans.x - startTrans->x) > metersPerPixel || fabs(currentTrans.y - startTrans->y) > metersPerPixel
  //      || fabs(currentTrans.theta - startTrans->theta) > thetaResolution) {
  //    fprintf(stderr,
  //        "WARNING:gradient descent moved more than 1 cell away!, totalNumSteps= %d, numSteps at lowest level=%d\n",
  //        totalNumSteps, numSteps);
  //    fprintf(
  //        stderr,
  //        "polished: score=%f  (%.6f,%.6f,%.6f),  \t normal: score=%f  (%.6f,%.6f,%.6f) \t delta: score=%f  (%.6f,%.6f,%.6f)\n",
  //        currentTrans.score, currentTrans.x, currentTrans.y, currentTrans.theta, startTrans->score, startTrans->x,
  //        startTrans->y, startTrans->theta, currentTrans.score - startTrans->score, currentTrans.x - startTrans->x,
  //        currentTrans.y - startTrans->y, currentTrans.theta - startTrans->theta);
  //  }

  //get the number of hits
  currentTrans.hits = rlt->getNumHits(points, numPoints, &currentTrans, hitThresh);

  if (useThreads)
    g_mutex_unlock(rlt_mutex);

  return currentTrans;
}

void ScanMatcher::addScanToBeProcessed(frsmPoint * points, unsigned numPoints, ScanTransform * T,
    frsm_laser_type_t laser_type, int64_t utime)
{
  if (!useThreads) {
    addScan(points, numPoints, T, laser_type, utime);
    if (verbose)
      fprintf(stderr, "done\n");
    return;
  }
  else {
    Scan * s = new Scan(numPoints, points, *T, laser_type, utime);
    g_mutex_lock(toBeProcessed_mutex);
    scansToBeProcessed.push_back(s);
    cancelAdd = false;
    g_mutex_unlock(toBeProcessed_mutex);
    g_cond_broadcast(toBeProcessed_cv);
  }
}

void ScanMatcher::addScan(frsmPoint * points, unsigned numPoints, ScanTransform * T, frsm_laser_type_t laser_type,
    int64_t utime, bool rebuildNow)
{
  Scan * s = new Scan(numPoints, points, *T, laser_type, utime, true);
  addScan(s, rebuildNow);
}

void ScanMatcher::addScan(Scan *s, bool rebuildNow)
{

  if (useThreads)
    g_mutex_lock(scans_mutex);
  if (s != NULL) {
    scans.push_back(s);
    if (scans.size() > maxNumScans) {
      delete scans.front();
      scans.pop_front();
    }
  }

  if (rebuildNow) {
    list<Scan *>::iterator it;
    for (it = scans.begin(); it != scans.end(); ++it) {
      Scan * scan = *it;
      if (scan->contours.size() == 0) {
        frsm_tictoc("findContours");
        if (contour_extractor == NULL && contour_extractor->laserType != s->laser_type) {
          delete contour_extractor;
          contour_extractor = NULL;
        }
        if (contour_extractor == NULL
        )
          contour_extractor = new ContourExtractor(s->laser_type);
        contour_extractor->findContours(scan->ppoints, scan->numPoints, scan->contours);
        frsm_tictoc("findContours"); //      s->drawContours();
      }
    }

    if (useThreads)
      g_mutex_lock(rlt_mutex);

    //rebuild now
    //    frsm_tictoc("rebuildRaster_olson");
    //    rebuildRaster_olson(&rlt);
    //    frsm_tictoc("rebuildRaster_olson");

    //  frsm_tictoc("rebuildRaster_blur");
    //  rebuildRaster_blur(&rlt);
    //  frsm_tictoc("rebuildRaster_blur");

    //        frsm_tictoc("rebuildRaster_blurLine");
    //        rebuildRaster_blurLine(&rlt);
    //        frsm_tictoc("rebuildRaster_blurLine");

    frsm_tictoc("rebuildRaster");
    rebuildRaster(&rlt);
    frsm_tictoc("rebuildRaster");

    if (useMultiRes > 0) {
      frsm_tictoc("rebuild_lowRes");
      if (rlt_low_res != NULL)
        delete rlt_low_res;
      rlt_low_res = new RasterLookupTable(rlt, downsampleFactor);
      frsm_tictoc("rebuild_lowRes");
    }

    if (useThreads) {
      g_mutex_unlock(rlt_mutex);
    }
  }
  if (useThreads) {
    g_mutex_unlock(scans_mutex);
  }

}

ScanTransform ScanMatcher::matchSuccessive(frsmPoint * points, unsigned numPoints, frsm_laser_type_t laser_type,
    int64_t utime, bool preventAddScan, ScanTransform * prior)
{
  int xSat = 0, ySat = 0, thetaSat = 0;
  if (numScans() > 0) {
    frsm_tictoc("scanMatch");

    double xRange1 = initialSearchRangeXY;
    double yRange1 = initialSearchRangeXY;
    double thetaRange1 = initialSearchRangeTheta;

    ScanTransform poseGuess;
    if (prior != NULL) {
      poseGuess = *prior;
    }
    else {
      if (stationaryMotionModel) {
        //guess Velocity is 0
        poseGuess.x = currentPose.x;
        poseGuess.y = currentPose.y;
        poseGuess.theta = currentPose.theta;
      }
      else {
        //guess that velocity is constant, but capped at the maxSearchRange...
        double dx = (currentPose.x - prevPose.x);
        if (fabs(dx) > maxSearchRangeXY)
          dx = frsm_fsign(dx) * maxSearchRangeXY;
        poseGuess.x = currentPose.x + dx;

        double dy = (currentPose.y - prevPose.y);
        if (fabs(dy) > maxSearchRangeXY)
          dy = frsm_fsign(dy) * maxSearchRangeXY;
        poseGuess.y = currentPose.y + dy;

        double dt = frsm_normalize_theta(currentPose.theta - prevPose.theta);
        if (fabs(dt) > maxSearchRangeTheta)
          dt = frsm_fsign(dt) * maxSearchRangeTheta;
        poseGuess.theta = frsm_normalize_theta(currentPose.theta + dt);
      }
      poseGuess.score = motionModelPriorWeight;
    }
    //make sure that the search range is large enough for the current velocity!
    xRange1 = fmax(xRange1, fabs(1.5 * (currentPose.x - prevPose.x)));
    yRange1 = fmax(yRange1, fabs(1.5 * (currentPose.y - prevPose.y)));
    thetaRange1 = fmax(thetaRange1, fabs(1.5 * frsm_normalize_theta(currentPose.theta - prevPose.theta)));

    prevPose = currentPose;

    int numTries = 3;
    int tr = 0;
    if (matchingMode == FRSM_Y_GRID_COORD) {
      xRange1 = 0; //gets rounded up to 1 pixel
      numTries = 1;
    }
    if (matchingMode < FRSM_RUN_GRID) {
      for (tr = 0; tr < numTries; tr++) {
        xRange1 = fmin(xRange1, maxSearchRangeXY);
        yRange1 = fmin(yRange1, maxSearchRangeXY);
        thetaRange1 = fmin(thetaRange1, maxSearchRangeTheta);
        if (tr > 0)
          fprintf(stderr, "xRange=%f, yRange=%f, thetaRange=%f\n", xRange1, yRange1, thetaRange1);
        currentPose = gridMatch(points, numPoints, &poseGuess, xRange1, yRange1, thetaRange1, &xSat, &ySat, &thetaSat);
        if (matchingMode == FRSM_Y_GRID_COORD)
          xSat = 0; //x doesn't matter
        if (!(xSat || ySat || thetaSat))
          break;
        else {
          if (thetaSat) {
            fprintf(stderr, "Hit Theta!");
            thetaRange1 *= 2;
          }
          else {
            poseGuess = currentPose; //hopefully all are close
            //          poseGuess.theta = currentPose.theta; //other ones are probably crap
            thetaRange1 = thetaResolution * 6; //theta should be almost correct!
            if (xSat && ySat) {
              fprintf(stderr, "hit BOTH x and y!");
              xRange1 *= 2;
              yRange1 *= 2;
            }
            else if (xSat) {
              fprintf(stderr, "just hit X");
              xRange1 *= 3;
              yRange1 /= 1.5;
              poseGuess.y = currentPose.y;
            }
            else if (ySat) {
              fprintf(stderr, "just hit Y");
              xRange1 /= 3;
              yRange1 *= 1.5;
              poseGuess.x = currentPose.x;
            }
            fprintf(stderr, " ... trying again... ");

          }
        }

      }

      if (matchingMode == FRSM_GRID_COORD) {
        frsm_tictoc("GradientAscent Polish");
        ScanTransform polished = coordAscentMatch(points, numPoints, &currentPose);
        currentPose = polished;
        //        }
        frsm_tictoc("GradientAscent Polish");
      }
    }
    else {
      frsm_tictoc("GradientAscent Match");
      ScanTransform polished = coordAscentMatch(points, numPoints, &poseGuess);
      memset(polished.sigma, 0, 9 * sizeof(double));
      polished.sigma[0] = polished.sigma[4] = polished.sigma[8] = .0001;
      currentPose = polished;
      frsm_tictoc("GradientAscent Match");

    }
    //make sure that the variances are nonzero
    currentPose.sigma[0] = fmax(1e-12, currentPose.sigma[0]);
    currentPose.sigma[4] = fmax(1e-12, currentPose.sigma[4]);
    currentPose.sigma[8] = fmax(1e-12, currentPose.sigma[8]);

    //increase variances if we hit a rail and had to try more than once...
    if (tr > 1) {
      currentPose.sigma[0] *= tr;
      currentPose.sigma[1] *= tr;
      currentPose.sigma[4] *= tr;
      if (xSat || ySat || thetaSat)
        fprintf(stderr, "Warning!hit a rail in the end!\n");
      else
        fprintf(stderr, "Found large enough search region\n");
    }
    if (xSat) {
      fprintf(stderr, "hit rail in X direction %d\n", xSat);
      currentPose.sigma[0] *= 3;
      currentPose.sigma[1] *= 2;
      currentPose.sigma[4] *= 1.5;
    }
    if (ySat) {
      fprintf(stderr, "hit rail in Y direction %d\n", ySat);
      currentPose.sigma[0] *= 1.5;
      currentPose.sigma[1] *= 2;
      currentPose.sigma[4] *= 3;
    }
    if (thetaSat) {
      fprintf(stderr, "hit rail in theta direction %d\n", thetaSat);
      currentPose.sigma[0] *= 3;
      currentPose.sigma[1] *= 3;
      currentPose.sigma[4] *= 3;
    }

    frsm_tictoc("scanMatch");

  }
  if (matchingMode == FRSM_Y_COORD_ONLY || matchingMode == FRSM_Y_GRID_COORD) {
    currentPose.x = prevPose.x; //don't let things move in x...
  }

  if (numScans() > 0) {
    bool addScan = false;
    if ((double) currentPose.hits / (double) numPoints < addScanHitThresh) {
      if (verbose)
        fprintf(stderr, "hits = %.1f%% should add scan...  ", (double) currentPose.hits / (double) numPoints * 100);
      addScan = true;
    }
    else {
      if (useThreads) { //dip in hits was temporary... don't add the map if we haven't yet
        g_mutex_lock(rlt_mutex);
        cancelAdd = true;
        g_mutex_unlock(rlt_mutex);
      }
    }

    if (addScan && preventAddScan) {
      if (verbose)
        fprintf(stderr, "NOT adding scan due to preventAddScan! \n");
      addScan = false;
    }
    if (addScan && (xSat || ySat || thetaSat)) {
      if (verbose)
        fprintf(stderr, "NOT adding scan due to hitting a rail! \n");
      addScan = false;
    }

    if (addScan) {
      frsm_tictoc("addScan");
      addScanToBeProcessed(points, numPoints, &currentPose, laser_type, utime);
      frsm_tictoc("addScan");
    }
  }
  else {
    if (verbose)
      fprintf(stderr, "adding first scan\n");
    if (prior != NULL
    )
      memcpy(&currentPose, prior, sizeof(currentPose));
    addScan(points, numPoints, &currentPose, laser_type, utime); //do a blocking add...

    if (verbose)
      fprintf(stderr, "done adding first scan\n");

    //set the covariance to be nonzero so it doesn't crash isam
    currentPose.sigma[0] = .1;
    currentPose.sigma[4] = .1;
    currentPose.sigma[8] = .1;
  }

  frsm_tictoc("laser_handler");

  return currentPose;

}

int ScanMatcher::isUsingIPP()
{
#ifdef USE_IPP
  return 1;
#else
  return 0;
#endif
}

#ifndef NO_BOT_LCMGL

void ScanMatcher::draw_state_lcmgl(bot_lcmgl_t * lcmgl)
{
  if (useThreads) {
    g_mutex_lock(scans_mutex);
    g_mutex_lock(rlt_mutex);
  }

  //draw the gridmap
  rlt->draw_lcmgl(lcmgl);

  // draw each scan.
  bot_lcmgl_line_width(lcmgl, 2);
  bot_lcmgl_point_size(lcmgl, 4);
  list<Scan *>::iterator it;
  for (it = scans.begin(); it != scans.end(); ++it) {
    Scan * s = *it;
    for (unsigned cidx = 0; cidx < s->contours.size(); cidx++) {
      bot_lcmgl_begin(lcmgl, LCMGL_LINE_STRIP);
      bot_lcmgl_color3f(lcmgl, 0, 1, 0);
      for (unsigned i = 0; i < s->contours[cidx]->points.size(); i++) {
        frsmPoint p0 = s->contours[cidx]->points[i];
        bot_lcmgl_vertex2f(lcmgl, p0.x, p0.y);
      }
      bot_lcmgl_end(lcmgl);
      bot_lcmgl_color3f(lcmgl, 0, 1, 1);
      bot_lcmgl_begin(lcmgl, LCMGL_POINTS);
      for (unsigned i = 0; i < s->contours[cidx]->points.size(); i++) {
        frsmPoint p0 = s->contours[cidx]->points[i];
        bot_lcmgl_vertex2f(lcmgl, p0.x, p0.y);
      }
      bot_lcmgl_end(lcmgl);
    }
  }

  if (useThreads) {
    g_mutex_unlock(scans_mutex);
    g_mutex_unlock(rlt_mutex);
  }
}

void ScanMatcher::draw_scan_lcmgl(bot_lcmgl_t * lcmgl, frsmPoint * points, unsigned numPoints, const ScanTransform * T)
{
  //draw the current scan
  bot_lcmgl_point_size(lcmgl, 4);
  bot_lcmgl_color3f(lcmgl, 1, 0, 0);
  Scan * s = new Scan(numPoints, points, *T, FRSM_DUMMY_LASER, 0, false);
  bot_lcmgl_begin(lcmgl, LCMGL_POINTS);
  for (int i = 0; i < numPoints; i++) {
    bot_lcmgl_vertex3f(lcmgl, s->ppoints[i].x, s->ppoints[i].y, 0);
  }
  bot_lcmgl_end(lcmgl);
  delete s;

  //draw the robot location

  bot_lcmgl_line_width(lcmgl, 4);
  double xyz[3] = { T->x, T->y, 0 };
  bot_lcmgl_color3f(lcmgl, 0, 0, 1);
  bot_lcmgl_circle(lcmgl, xyz, .5);
  double heading[3] = { T->x + cos(T->theta), T->y + sin(T->theta), 0 };
  bot_lcmgl_begin(lcmgl, LCMGL_LINES);
  bot_lcmgl_color3f(lcmgl, 1, 1, 0);
  bot_lcmgl_vertex3f(lcmgl, xyz[0], xyz[1], xyz[2]);
  bot_lcmgl_vertex3f(lcmgl, heading[0], heading[1], heading[2]);
  bot_lcmgl_end(lcmgl);
}

#endif

} //namespace frsm
