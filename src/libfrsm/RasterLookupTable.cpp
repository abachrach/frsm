#ifdef USE_IPP
#include "ipp.h"
#endif

#include "RasterLookupTable.hpp"
#include "ScanMatchingUtils.hpp"
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include "ScanMatchingUtils.hpp"
#include <vector>
#include <algorithm>
#include <iostream>
#include <iomanip>

#define DRAW_COST_SURFACE 0
#if DRAW_COST_SURFACE
//for debugging/visualization
#include <bot_lcmgl_client/lcmgl.h>
#include <GL/gl.h>
#endif

#define VAR_COMP_THRESH .80

using namespace std;
namespace frsm {

typedef struct {
  double x[3];
  double p;
} cov_sample_t;

RasterLookupTable::~RasterLookupTable()
{
  delete table;
}

RasterLookupTable::RasterLookupTable(double x0i, double y0i, double x1i, double y1i, double mPP, int pixelDivisor,
    uint8_t initialValue)
{
  pixelsPerMeter = 1.0 / mPP;

  int widthI = ceil((pixelsPerMeter) * (x1i - x0i));
  int heightI = ceil((pixelsPerMeter) * (y1i - y0i));

  //make it divisible by pixelDivisor
  widthI = ceil(widthI / (double) pixelDivisor) * pixelDivisor;
  heightI = ceil(heightI / (double) pixelDivisor) * pixelDivisor;

  double xy0[2] = { x0i, y0i };
  double xy1[2] = { x0i + widthI * mPP, y0i + heightI * mPP };

  table = new occ_map::Uint8PixelMap(xy0, xy1, mPP, initialValue, true, false);

  //fill out conveniance variables
  x0 = table->xy0[0];
  y0 = table->xy0[1];
  x1 = table->xy1[0];
  y1 = table->xy1[1];
  width = table->dimensions[0];
  height = table->dimensions[1];
  metersPerPixel = table->metersPerPixel;
  pixelsPerMeter = 1.0 / table->metersPerPixel;
}

RasterLookupTable::RasterLookupTable(RasterLookupTable * hi_res, int downsampleFactor)
{
  double x0i = hi_res->table->xy0[0];
  double y0i = hi_res->table->xy0[1];
  double x1i = hi_res->table->xy1[0];
  double y1i = hi_res->table->xy1[1];
  metersPerPixel = hi_res->table->metersPerPixel * downsampleFactor;
  pixelsPerMeter = 1.0 / metersPerPixel;

  width = round((pixelsPerMeter) * (x1i - x0i));
  height = round((pixelsPerMeter) * (y1i - y0i));
  assert(width == hi_res->width / downsampleFactor);
  assert(height == hi_res->height / downsampleFactor);

  double xy0[2] = { x0i, y0i };
  double xy1[2] = { x0i + width * metersPerPixel, y0i + height * metersPerPixel };
  table = new occ_map::Uint8PixelMap(xy0, xy1, metersPerPixel, 0, true, false);

  //TODO: Is this too slow?
  frsm_tictoc("downsample_exp");
  //downsample the high res table
  //each pixel is set to the max of all its high-res counterparts
  for (int i = 0; i < hi_res->table->dimensions[1]; i++) {
    for (int j = 0; j < hi_res->table->dimensions[0]; j++) {
      int lind = i / downsampleFactor * width + j / downsampleFactor;
      int hind = i * hi_res->table->dimensions[0] + j;
      table->data[lind] = frsm_ucmax(table->data[lind], hi_res->table->data[hind]);
    }
  }
  frsm_tictoc("downsample_exp");

  //fill out conveniance variables
  x0 = table->xy0[0];
  y0 = table->xy0[1];
  x1 = table->xy1[0];
  y1 = table->xy1[1];
  width = table->dimensions[0];
  height = table->dimensions[1];
  metersPerPixel = table->metersPerPixel;
  pixelsPerMeter = 1.0 / table->metersPerPixel;
}

void RasterLookupTable::drawRectangle(double cx, double cy, double x_size, double y_size, double theta,
    uint8_t * lutSq, int lutSq_size, int lutSq_first_zero, double lutSqRange)
{
  double ux = cos(theta), uy = sin(theta);

  double lutRange = sqrt(lutSqRange);
  double invLutSqRange = 1.0 / lutSqRange;

  double x_bound = (x_size / 2.0 * fabs(ux) + y_size / 2.0 * fabs(uy)) + lutRange;
  double y_bound = (x_size / 2.0 * fabs(uy) + y_size / 2.0 * fabs(ux)) + lutRange;

  int ix0 = frsm_clamp((int) ((cx - x_bound - x0) * pixelsPerMeter), 0, width - 1);
  int ix1 = frsm_clamp((int) ((cx + x_bound - x0) * pixelsPerMeter), 0, width - 1);

  int iy0 = frsm_clamp((int) ((cy - y_bound - y0) * pixelsPerMeter), 0, height - 1);
  int iy1 = frsm_clamp((int) ((cy + y_bound - y0) * pixelsPerMeter), 0, height - 1);

  double y = y0 + (iy0 + .5) * metersPerPixel;

  for (int iy = iy0; iy <= iy1; iy++) {

    double x = x0 + (ix0 + .5) * metersPerPixel;

    for (int ix = ix0; ix <= ix1; ix++) {

      // distances from query point to center of rectangle
      double dx = x - cx, dy = y - cy;

      // how long are the projections of the vector (dx,dy) onto the
      // two principle
      // components of the rectangle? How much longer are they than
      // the dimensions
      // of the rectangle?
      double c1 = fabs(dx * ux + dy * uy) - (x_size / 2);
      double c2 = fabs(-dx * uy + dy * ux) - (y_size / 2);

      // if the projection length is < 0, we're *inside* the
      // rectangle.
      c1 = fmax(0, c1);
      c2 = fmax(0, c2);

      double distSq = frsm_sq(c1) + frsm_sq(c2);

      int lutSqIdx = (int) (lutSq_size * distSq * invLutSqRange + .5);
      //      printf("dist = %f, lutSqIdx=%d, lutSq_first_zero=%d\n",distSq,lutSqIdx,lutSq_first_zero);

      if (lutSqIdx < lutSq_first_zero) {
        int idx = iy * width + ix;
        table->data[idx] = frsm_ucmax(table->data[idx], lutSq[lutSqIdx]);
      }

      x += metersPerPixel;
    }

    y += metersPerPixel;
  }

}

void RasterLookupTable::drawKernel(int ix, int iy, const uint8_t*kernel, int kernel_width, int kernel_height)
{
  //TODO: compare against opencv/IPP?
  ix -= kernel_width / 2;
  iy -= kernel_height / 2;
  for (int j = 0; j < kernel_height; j++) {
    int table_iy = iy + j;
    int krow = j * kernel_width;
    for (int i = 0; i < kernel_width; i++) {
      int table_ix = ix + i;
      if (readTable(table_ix, table_iy) < kernel[krow + i])
        writeTable(table_ix, table_iy, kernel[krow + i]);
    }
  }

}

void RasterLookupTable::drawBlurredPoint(const frsmPoint * p, const LutKernel * kern)
{
  int ix, iy;
  worldToTable(p->x, p->y, &ix, &iy);
  drawKernel(ix, iy, kern->square_kernel, kern->kernel_size, kern->kernel_size);

}
void RasterLookupTable::drawBlurredLine(const frsmPoint * p1, const frsmPoint * p2, const LutKernel * kern)
{
  frsmPoint d = { p2->x - p1->x, p2->y - p1->y };
  if (d.x == 0 && d.y == 0) {
    return;
  }

  double slope = atan2(d.y, d.x);
  int kernel_width = -1;
  int kernel_height = -1;
  uint8_t * kernel_pointer;
  if (fabs(d.x) > fabs(d.y)) { //vertical kernel
    int slope_ind = kern->slope_step * fabs(slope);
    kernel_width = 1;
    kernel_height = kern->line_kernel_sizes[slope_ind];
    kernel_pointer = kern->line_kernels + (slope_ind * kern->line_kernel_stride);
  }
  else {
    int slope_ind = kern->slope_step * fabs(frsm_angle_subtract(M_PI / 2.0, slope));
    kernel_width = kern->line_kernel_sizes[slope_ind];
    kernel_height = 1;
    kernel_pointer = kern->line_kernels + (slope_ind * kern->line_kernel_stride);
  }

  int start[2], end[2];
  worldToTable(p1->x, p1->y, &start[0], &start[1]);
  worldToTable(p2->x, p2->y, &end[0], &end[1]);

  int curr[2] = { start[0], start[1] };

  // normalize
  int xstep = 1;
  int ystep = 1;
  int dx = end[0] - start[0];
  if (dx < 0) {
    dx = -dx;
    xstep = -1;
  }
  int dy = end[1] - start[1];
  if (dy < 0) {
    dy = -dy;
    ystep = -1;
  }

  if (dx == 0) {
    // vertical
    for (int i = 0; i <= dy; i++) {
      int wasHit = curr[0] == end[0] && curr[1] == end[1];
      drawKernel(curr[0], curr[1], kernel_pointer, kernel_width, kernel_height);
      curr[1] = curr[1] + ystep;
    }
  }
  else if (dy == 0) {
    // horizontal
    for (int i = 0; i <= dx; i++) {
      int wasHit = curr[0] == end[0] && curr[1] == end[1];
      drawKernel(curr[0], curr[1], kernel_pointer, kernel_width, kernel_height);
      curr[0] += xstep;
    }
  }
  else if (dx > dy) {
    // bresenham, horizontal slope
    int n = dx;
    dy += dy;
    int e = dy - dx;
    dx += dx;

    for (int i = 0; i <= n; i++) {
      int wasHit = curr[0] == end[0] && curr[1] == end[1];
      drawKernel(curr[0], curr[1], kernel_pointer, kernel_width, kernel_height);
      if (e >= 0) {
        curr[1] += ystep;
        e -= dx;
      }
      e += dy;
      curr[0] += xstep;
    }
  }
  else {
    // bresenham, vertical slope
    int n = dy;
    dx += dx;
    int e = dx - dy;
    dy += dy;

    for (int i = 0; i <= n; i++) {
      int wasHit = curr[0] == end[0] && curr[1] == end[1];
      drawKernel(curr[0], curr[1], kernel_pointer, kernel_width, kernel_height);
      if (e >= 0) {
        curr[0] += xstep;
        e -= dy;
      }
      e += dx;
      curr[1] += ystep;
    }
  }

}

static int computeKernelSize(double sigma, double slope, float cutoff)
{
  unsigned int sz = 0;
  float val = 255;
  while (val > cutoff) {
    sz++;
    double d = sz * cos(slope);
    val = round(255.0 * exp(-d / sigma));
  }
  return 2 * sz + 1;
}

int RasterLookupTable::getNumHits(const frsmPoint * points, const unsigned numPoints, const ScanTransform * XYT0,
    int hitThresh)
{
  int hits = 0;
  double ct = cos(XYT0->theta), st = sin(XYT0->theta);
  // Evaluate each point for a fixed transform
  for (unsigned pidx = 0; pidx < numPoints; pidx++) {

    // project the point
    frsmPoint p = points[pidx];
    double x = p.x * ct - p.y * st + XYT0->x;
    double y = p.x * st + p.y * ct + XYT0->y;
    int ix, iy;
    worldToTable(x, y, &ix, &iy);
    if (table->data[iy * width + ix] >= hitThresh)
      hits++;
  }
  return hits;
}

float RasterLookupTable::getScore(const frsmPoint * points, const unsigned numPoints, const ScanTransform * XYT0)
{
  float score = 0.0;
  double ct = cos(XYT0->theta), st = sin(XYT0->theta);
  // Evaluate each point for a fixed transform
  for (unsigned pidx = 0; pidx < numPoints; pidx++) {
    // project the point
    frsmPoint p = points[pidx];
    double x = p.x * ct - p.y * st + XYT0->x;
    double y = p.x * st + p.y * ct + XYT0->y;
    int ix, iy;
    worldToTable(x, y, &ix, &iy);
    score += table->data[iy * width + ix];
  }
  return score;
}

float RasterLookupTable::getScoreDump(const frsmPoint * points, const unsigned numPoints, const ScanTransform * XYT0,
    const char * name)
{
  FILE * f = fopen(name, "w");
  float score = 0.0;
  double ct = cos(XYT0->theta), st = sin(XYT0->theta);
  // Evaluate each point for a fixed transform
  for (unsigned pidx = 0; pidx < numPoints; pidx++) {
    // project the point
    frsmPoint p = points[pidx];
    double x = p.x * ct - p.y * st + XYT0->x;
    double y = p.x * st + p.y * ct + XYT0->y;
    int ix, iy;
    worldToTable(x, y, &ix, &iy);
    score += table->data[iy * width + ix];
    fprintf(f, "blah %d %d %d %f %f %d %f \n", pidx, ix, iy, x, y, table->data[iy * width + ix], score);
  }
  fclose(f);
  return score;
}

ScanTransform RasterLookupTable::evaluate2D(const frsmPoint * points, const unsigned numPoints,
    const ScanTransform * XYT0, const ScanTransform * prior, int ixrange, int iyrange, int ixdim, int iydim,
    float * scores, int * bestScoreIndX, int *bestScoreIndY)
{

#ifdef USE_IPP
  int scoresStep = ixdim * sizeof(float);
  int dataStep = width * sizeof(uint8_t);
#endif
  //zero out the scores array!
  memset(scores, 0, ixdim * iydim * sizeof(float));

  double ct = cos(XYT0->theta), st = sin(XYT0->theta);

  // Evaluate each point for a fixed rotation but variable
  // translation
  frsm_tictoc("forloop_acum");
  for (unsigned pidx = 0; pidx < numPoints; pidx++) {

    // project the point
    frsmPoint p = points[pidx];
    double x = p.x * ct - p.y * st + XYT0->x;
    double y = p.x * st + p.y * ct + XYT0->y;

    // (ix0, iy0) are the coordinates in table->data that
    // correspond to scores.x. It's the (nominal) upper-left
    // corner of our search window.
    int ix, iy;
    worldToTable(x, y, &ix, &iy);
    int ix0 = ix - ixrange;
    int iy0 = iy - iyrange;

    // compute the intersection of the box
    // (ix0,iy0)->(ix0+ixdim-1,iy0+iydim-1) and the box
    // (0,0)->(width-1, height-1). This will be our actual
    // search window.
    int bx0 = frsm_imax(ix0, 0);
    int by0 = frsm_imax(iy0, 0);

    int bx1 = frsm_imin(ix0 + ixdim - 1, width - 1);
    int by1 = frsm_imin(iy0 + iydim - 1, height - 1);

    if (by1 < by0 || bx1 < bx0)
      continue; //point is way off map!

#ifdef USE_IPP
    //ipp
    IppiSize roiSize = {bx1 - bx0 + 1, by1 - by0 + 1}; //+1 due to <= below
    uint8_t * pSrc = table->data+ width*by0+bx0;
    float * pScoreAcum = scores + (by0 - iy0)*ixdim + (bx0 - ix0);
    ippiAdd_8u32f_C1IR(pSrc, dataStep, pScoreAcum, scoresStep, roiSize);
#else

    for (int iy = by0; iy <= by1; iy++) {
      int sy = iy - iy0; // y coordinate in scores[]
      for (int ix = bx0; ix <= bx1; ix++) {

        int lutval = table->data[iy * width + ix];
        int sx = ix - ix0;

        int sidx = sy * ixdim + sx;

        if (lutval > 0) {
          scores[sidx] += lutval;
        }
      }
    }
#endif

  }
  frsm_tictoc("forloop_acum");

  //factor in wide gaussian prior
  if (prior->score > .1) {
    double x, y;
    for (int sy = 0; sy < iydim; sy++) {
      for (int sx = 0; sx < ixdim; sx++) {
        int sidx = sy * ixdim + sx;
        scoresToWorld(XYT0, ixrange, iyrange, sx, sy, &x, &y);
        //assuming diag covariance
        scores[sidx] *= exp(-1.0 / prior->score * (frsm_sq(x - prior->x) + frsm_sq(y - prior->y)));
      }
    }
  }

  // Look for the best score
  float bestscore = -1;
  int bestsx = 0, bestsy = 0;
  for (int sy = 0; sy < iydim; sy++) {
    for (int sx = 0; sx < ixdim; sx++) {
      int sidx = sy * ixdim + sx;
      if (scores[sidx] > bestscore) {
        bestscore = scores[sidx];
        bestsx = sx;
        bestsy = sy;
      }
    }
  }
  *bestScoreIndX = bestsx;
  *bestScoreIndY = bestsy;

  int bsidx = bestsy * ixdim + bestsx;

  ScanTransform result;
  scoresToWorld(XYT0, ixrange, iyrange, bestsx, bestsy, &result.x, &result.y);
  result.theta = XYT0->theta;

  result.score = scores[bsidx];

  return result;
}

class score_entry {
public:
  score_entry(float score_, int it_, int ix_, int iy_) :
      score(score_), it(it_), ix(ix_), iy(iy_)
  {
  }
  float score;
  int it;
  int ix;
  int iy;
};

bool score_entry_comp(const score_entry & i, const score_entry & j)
{
  return i.score > j.score;
}

/** Perform a brute-force search in 3DOFs. * */
ScanTransform RasterLookupTable::evaluate3D_multiRes(RasterLookupTable * rlt_high_res, const frsmPoint * points,
    const unsigned numPoints, const ScanTransform * prior, double xrange, double yrange, double thetarange,
    double thetastep, int hitThresh, int * xSat, int *ySat, int * thetaSat)
{
  ScanTransform XYT0;

  //round the prior to a cell
  int ix, iy;
  worldToTable(prior->x, prior->y, &ix, &iy);
  tableToWorld(ix, iy, &XYT0.x, &XYT0.y);

  XYT0.theta = prior->theta;

  ScanTransform bestResult;
  bestResult.score = -1;

  // search range in pixels
  int ixrange = (int) (xrange * pixelsPerMeter) + 1; // always round up,
  // make sure > 0
  int iyrange = (int) (yrange * pixelsPerMeter) + 1; // ditto...

  // allocate our scoring arrays.
  int ixdim = (2 * ixrange + 1);
  int iydim = (2 * iyrange + 1);

  int itdim = (int) (2 * thetarange / thetastep) + 1;

  float * allScores;
  int numScores = itdim * ixdim * iydim;
  allScores = (float *) malloc(numScores * sizeof(float));
  int allScores_step = ixdim * iydim;

  //compute the scores for the full 3D voxel grid at the low resolution
  ScanTransform xyt;
  xyt.x = XYT0.x;
  xyt.y = XYT0.y;
  double dtheta = -thetarange;
  for (int it = 0; it < itdim; it++) {
    xyt.theta = XYT0.theta + dtheta;
    int xInd = 0, yInd = 0;
    ScanTransform thisResult = evaluate2D(points, numPoints, &xyt, prior, ixrange, iyrange, ixdim, iydim, allScores
        + (it * allScores_step), &xInd, &yInd);

    //    //sanity check
    //    float sc = getScore(points, numPoints, &xyt);
    //    float sc2;
    //    ScanTransform tc = evaluate2D(points, numPoints, &xyt, 0, 0, 1, 1, &sc2, &xInd, &yInd);
    //    fprintf(stderr, "e2dscore=%f, scorecheck=%f\n", sc2, sc);
    //
    //    for (int iy = 0; iy < iydim; iy++) {
    //      for (int ix = 0; ix < ixdim; ix++) {
    //
    //        ScanTransform t;
    //        scoresToWorld(&xyt, ixrange, iyrange, ix, iy, &t.x, &t.y);
    //        t.theta = xyt.theta;
    //        float sc = getScore(points, numPoints, &t);
    //        fprintf(stderr, "(%f,%f,%f) (%d,%d,%d) score=%f, scorecheck=%f\n", t.x, t.y,t.theta,ix,iy,it,
    //            allScores[it * allScores_step + iy * ixdim + ix], sc);
    //        assert(fabs(sc-allScores[it * allScores_step + iy * ixdim + ix])<1e-6);
    //      }
    //    }

    dtheta += thetastep;
  }

  //have upper bounds on scores from low res, now find best
  std::vector<score_entry> sorted_scores;
  sorted_scores.reserve(numScores);
  for (int it = 0; it < itdim; it++) {
    for (int iy = 0; iy < iydim; iy++) {
      for (int ix = 0; ix < ixdim; ix++) {
        float score = allScores[it * allScores_step + iy * ixdim + ix];
        sorted_scores.push_back(score_entry(score, it, ix, iy));
      }
    }
  }
  std::sort(sorted_scores.begin(), sorted_scores.end(), score_entry_comp);

  // high_res search range in pixels
  int ixrange_high_res = (int) (metersPerPixel / 2.0 * rlt_high_res->pixelsPerMeter);
  // make sure > 0
  int iyrange_high_res = (int) (metersPerPixel / 2.0 * rlt_high_res->pixelsPerMeter);

  // allocate our scoring arrays.
  int ixdim_high_res = (2 * ixrange_high_res + 1);
  int iydim_high_res = (2 * iyrange_high_res + 1);

  float * scores_high_res;
  scores_high_res = (float *) calloc(ixdim_high_res * iydim_high_res, sizeof(float));

  //initialize the variables for the covariance computation
  vector<cov_sample_t> cov_samples;

  int i;
  int canBeDone = 0;
  int bestIx = -1, bestIy = -1, bestIt = -1;
  int cov_samples_ind = 0;
  double maxProb = 0;
  for (i = 0; i < numScores; i++) {
    //make sure we do at least the best 6 voxels for covariance and robustness...
    if (sorted_scores[i].score < bestResult.score && i > 6)
      canBeDone = 1; //scores from low_res are an upper bound, so we're done...
    //make sure we have checked the voxels around the best
    if (canBeDone && (abs(sorted_scores[i].ix - bestIx) + abs(sorted_scores[i].iy - bestIy) + abs(sorted_scores[i].it
        - bestIt)) != 1)
      continue;

    ScanTransform t;
    //convert from the scores index back to world...
    scoresToWorld(&XYT0, ixrange, iyrange, sorted_scores[i].ix, sorted_scores[i].iy, &t.x, &t.y);
    t.theta = XYT0.theta - thetarange + sorted_scores[i].it * thetastep;

    int xInd, yInd;
    ScanTransform result = rlt_high_res->evaluate2D(points, numPoints, &t, prior, ixrange_high_res, iyrange_high_res,
        ixdim_high_res, iydim_high_res, scores_high_res, &xInd, &yInd);

    //add these "samples" to the covariance computation
    assert(cov_samples_ind == cov_samples.size());
    cov_samples.resize(cov_samples.size() + ixdim_high_res * iydim_high_res);
    double wx, wy, p, score;
    for (int sy = 0; sy < iydim_high_res; sy++) {
      for (int sx = 0; sx < ixdim_high_res; sx++) {
        score = scores_high_res[sy * ixdim_high_res + sx];
        p = score / ((double) numPoints * 255.0);
        p = exp(p) / 2.71828183;
        cov_samples[cov_samples_ind].p = p;
        maxProb = fmax(maxProb, p);

        rlt_high_res->scoresToWorld(&t, ixrange_high_res, iyrange_high_res, sx, sy, &wx, &wy);
        cov_samples[cov_samples_ind].x[0] = wx;
        cov_samples[cov_samples_ind].x[1] = wy;
        cov_samples[cov_samples_ind].x[2] = t.theta;
        cov_samples_ind++;
      }
    }
    if (result.score > bestResult.score) {
      bestResult = result;
      bestIx = sorted_scores[i].ix;
      bestIy = sorted_scores[i].iy;
      bestIt = sorted_scores[i].it;
    }

    //    fprintf(stderr, "low_res_score=%f, hi_res_score=%f, it=%d iy=%d ix=%d\n", sorted_scores[i].score, bestResult.score,
    //        sorted_scores[i].it, sorted_scores[i].iy, sorted_scores[i].ix);
    if (result.score > sorted_scores[i].score) { //it should be an upper bound!
      //    TODO:Not sure why this is happening sometimes!
      //            fprintf(
      //                    stderr,
      //                    "WARNING! score num %d (%d,%d,%d) hi_res score of %f higher than low resolution of %f\n",
      //                    i, sorted_scores[i].ix, sorted_scores[i].iy,
      //                    sorted_scores[i].it, result.score, sorted_scores[i].score);
      //      dumpTable("low_res");
      //      rlt_high_res->dumpTable("hi_res");
      //      ScanTransform tc;
      //      scoresToWorld(&XYT0, ixrange, iyrange, sorted_scores[i].ix, sorted_scores[i].iy, &tc.x, &tc.y);
      //      tc.theta = XYT0.theta - thetarange + sorted_scores[i].it * thetastep;
      //      float lrs = getScoreDump(points, numPoints, &tc, "low_res_dump");
      //      float lrs2 = getScoreDump(points, numPoints, &result, "low_res_dump2");
      //      float lrs3 = getScoreDump(points, numPoints, &t, "low_res_dump3");
      //
      //      int ix, iy;
      //      int ix2, iy2;
      //      int ix3, iy3;
      //      worldToTable(tc.x, tc.y, &ix, &iy);
      //      worldToTable(result.x, result.y, &ix2, &iy2);
      //      worldToTable(t.x, t.y, &ix3, &iy3);
      //      fprintf(stderr, "t=(%f,%f),tc=(%f,%f),r=(%f,%f) t_i=(%d,%d),tc_i=(%d,%d), r_i=(%d,%d)\n", t.x, t.y, tc.x, tc.y,
      //          result.x, result.y, ix3, iy3, ix, iy, ix2, iy2);
      //
      //      float hrs = rlt_high_res->getScoreDump(points, numPoints, &tc, "high_res_dump");
      //      float hrs2 = rlt_high_res->getScoreDump(points, numPoints, &result, "high_res_dump2");
      //      float hrs3 = rlt_high_res->getScoreDump(points, numPoints, &t, "high_res_dump3");
      //      fprintf(stderr, "dump done, lrs=%f,lrs2=%f,lrs3=%f  hrs=%f,hrs2=%f\n", lrs, lrs2, lrs3, hrs, hrs2, hrs3);
      //      exit(1);
    }
  }

  //  fprintf(stderr, "only checked %d of %d voxels at high_res\n",i,ixdim*iydim*itdim);
  //get the number of hits for this transform in the high res table
  bestResult.hits = rlt_high_res->getNumHits(points, numPoints, &bestResult, hitThresh);

  //actually compute the covariance
  double K[9] = { 0 };
  double u[3] = { 0 };
  double s = 0;

  double x[3];
  double Ktmp[9];
  double prob;
  //    FILE * f = fopen("vardump.log","w");
  for (unsigned int i = 0; i < cov_samples.size(); i++) {
    memcpy(x, cov_samples[i].x, 3 * sizeof(double));
    prob = cov_samples[i].p;
    //      fprintf(f,"p %f %f %f %f\n",x[0],x[1],x[2],prob);
    if (prob / maxProb < VAR_COMP_THRESH)
      continue;

    //        K = K+ x*x'*p;
    frsm_vector_vector_outer_product_3d(x, x, Ktmp);
    frsm_vector_scale_nd(Ktmp, 9, prob);
    frsm_vector_add_nd(K, Ktmp, 9, K);
    //        u = u+x*p;
    frsm_vector_scale_3d(x, prob);
    frsm_vector_add_3d(u, x, u);
    //        s = s+p;
    s += prob;
  }
  //    fclose(f);

  //sig_hat = 1/s*K-1/s^2*u*u'
  frsm_vector_scale_nd(K, 9, 1.0 / s);

  double u_outer[9];
  frsm_vector_vector_outer_product_3d(u, u, u_outer);
  frsm_vector_scale_nd(u_outer, 9, 1.0 / (s * s));

  double cov[9];
  frsm_vector_sub_nd(K, u_outer, 9, cov);
  cov[8] = fmax(cov[8], 1e-3); //often times the top 85% will all be the same theta
  memcpy(bestResult.sigma, cov, 9 * sizeof(double));

#if DRAW_COST_SURFACE
  //draw score cloud with lcmgl
  static lcm_t * lcm = lcm_create(NULL);
  static bot_lcmgl_t *lcmgl = bot_lcmgl_init(lcm, "frsm_scores");
  bot_lcmgl_push_attrib(lcmgl,GL_DEPTH_BUFFER_BIT | GL_POINT_BIT | GL_CURRENT_BIT);
  bot_lcmgl_enable(lcmgl,GL_DEPTH_TEST);
  bot_lcmgl_depth_func(lcmgl,GL_LESS);

  bot_lcmgl_point_size(lcmgl, 6);
  bot_lcmgl_begin(lcmgl, GL_POINTS);
  for (unsigned int i = 0; i < cov_samples.size(); i++) {
    if (fabs(frsm_angle_subtract(cov_samples[i].x[2],bestResult.theta))>.01 ||cov_samples[i].p/maxProb <VAR_COMP_THRESH)
    continue;
    float * color3fv = frsm_color_util_jet(cov_samples[i].p/maxProb);
    bot_lcmgl_color3f(lcmgl, color3fv[0], color3fv[1], color3fv[2]);
    //      double x = cov_samples[i].x[0]-bestResult.x;
    //      double y = cov_samples[i].x[1]-bestResult.y;
    double x = cov_samples[i].x[0];
    double y = cov_samples[i].x[1];
    bot_lcmgl_vertex3f(lcmgl, x,y ,cov_samples[i].p/maxProb);
  }

  for (unsigned int i = 0; i < cov_samples.size(); i++) {
    if (fabs(frsm_angle_subtract(cov_samples[i].x[2],bestResult.theta))>.01)
    continue;
    float * color3fv = frsm_color_util_jet(cov_samples[i].p/maxProb);
    bot_lcmgl_color3f(lcmgl, color3fv[0], color3fv[1], color3fv[2]);
    //      double x = cov_samples[i].x[0]-bestResult.x;
    //      double y = cov_samples[i].x[1]-bestResult.y;
    double x = cov_samples[i].x[0];
    double y = cov_samples[i].x[1];
    bot_lcmgl_vertex3f(lcmgl, x+1,y ,cov_samples[i].p/maxProb);
  }

  bot_lcmgl_end(lcmgl);
  bot_lcmgl_pop_attrib(lcmgl);
  bot_lcmgl_switch_buffer(lcmgl);

#endif

  //set the sat variables to zero, cuz we're ignoring it for now...
  if (xSat != NULL)
    *xSat = 0;
  if (ySat != NULL)
    *ySat = 0;
  if (thetaSat != NULL)
    *thetaSat = 0;

  //cleanup
  free(scores_high_res);
  free(allScores);

  return bestResult;
}

/** Perform a brute-force search in 3DOFs. * */
ScanTransform RasterLookupTable::evaluate3D(const frsmPoint * points, const unsigned numPoints,
    const ScanTransform * prior, double xrange, double yrange, double thetarange, double thetastep, int hitThresh,
    int * xSat, int *ySat, int * thetaSat)
{
  ScanTransform xyt = *prior;
  ScanTransform bestResult;
  bestResult.score = -1;

  // search range in pixels
  int ixrange = (int) (xrange * pixelsPerMeter) + 1; // always round up,
  // make sure > 0
  int iyrange = (int) (yrange * pixelsPerMeter) + 1; // ditto...

  // allocate our scoring arrays.
  int ixdim = (2 * ixrange + 1);
  int iydim = (2 * iyrange + 1);

  int itdim = (int) (2 * thetarange / thetastep) + 1;

  float ** allScores;
  allScores = (float **) calloc(itdim, sizeof(float *));
  float * bestThetaScores;
  bestThetaScores = (float *) calloc(itdim, sizeof(float));

  for (int i = 0; i < itdim; i++) {
    allScores[i] = (float *) calloc(ixdim * iydim, sizeof(float));
  }

  int bestThetaInd = 0, bestXInd = 0, bestYInd = 0;
  double dtheta = -thetarange;
  //TODO: this could/should be threaded...
  for (int it = 0; it < itdim; it++) {
    xyt.theta = prior->theta + dtheta;

    int xInd = 0, yInd = 0;
    ScanTransform thisResult = evaluate2D(points, numPoints, &xyt, prior, ixrange, iyrange, ixdim, iydim,
        allScores[it], &xInd, &yInd);
    bestThetaScores[it] = thisResult.score;
    if (bestResult.score == -1 || thisResult.score > bestResult.score) {
      bestResult = thisResult;
      bestThetaInd = it;
      bestXInd = xInd;
      bestYInd = yInd;
    }

    dtheta += thetastep;
  }

  frsm_tictoc("compute_Variance");
  //compute covariance by fitting multivariate gaussian directly
  //initialize the variables for the covariance computation
  double K[9] = { 0 };
  double u[3] = { 0 };
  double s = 0;

  double x[3];
  double Ktmp[9];
  double prob;
  double score;
  double maxProb = exp(bestResult.score / ((double) numPoints * 255.0)) / 2.71828183;
  //    FILE * f = fopen("vardump.log","w");
  for (int st = 0; st < itdim; st++) {
    for (int sy = 0; sy < iydim; sy++) {
      for (int sx = 0; sx < ixdim; sx++) {
        score = allScores[st][sy * ixdim + sx];
        prob = score / ((double) numPoints * 255.0);
        prob = exp(prob) / 2.71828183;
        if (prob / maxProb < VAR_COMP_THRESH)
          continue;

        scoresToWorld(prior, ixrange, iyrange, sx, sy, &x[0], &x[1]);
        x[2] = prior->theta - thetarange + st * thetastep;
        //                fprintf(f,"p %f %f %f %f\n",x[0],x[1],x[2],prob);

        //        K = K+ x*x'*p;
        frsm_vector_vector_outer_product_3d(x, x, Ktmp);
        frsm_vector_scale_nd(Ktmp, 9, prob);
        frsm_vector_add_nd(K, Ktmp, 9, K);
        //        u = u+x*p;
        frsm_vector_scale_3d(x, prob);
        frsm_vector_add_3d(u, x, u);
        //        s = s+p;
        s += prob;
      }
    }
  }
  //    fclose(f);

  //sig_hat = 1/s*K-1/s^2*u*u'
  frsm_vector_scale_nd(K, 9, 1.0 / s);

  double u_outer[9];
  frsm_vector_vector_outer_product_3d(u, u, u_outer);
  frsm_vector_scale_nd(u_outer, 9, 1.0 / (s * s));

  double cov[9];
  frsm_vector_sub_nd(K, u_outer, 9, cov);
  cov[8] = fmax(cov[8], 1e-3); //often times the top 85% will all be the same theta
  memcpy(bestResult.sigma, cov, 9 * sizeof(double));

  frsm_tictoc("compute_Variance");

#if DRAW_COST_SURFACE
  //draw score cloud with lcmgl
  static lcm_t * lcm = lcm_create(NULL);
  static bot_lcmgl_t *lcmgl = bot_lcmgl_init(lcm, "frsm_scores2");
  bot_lcmgl_push_attrib(lcmgl,GL_DEPTH_BUFFER_BIT | GL_POINT_BIT | GL_CURRENT_BIT);
  bot_lcmgl_enable(lcmgl,GL_DEPTH_TEST);
  bot_lcmgl_depth_func(lcmgl,GL_LESS);

  bot_lcmgl_point_size(lcmgl, 6);
  bot_lcmgl_begin(lcmgl, GL_POINTS);

  for (int st = 0; st < itdim; st++) {
    for (int sy = 0; sy < iydim; sy++) {
      for (int sx = 0; sx < ixdim; sx++) {
        score = allScores[st][sy * ixdim + sx];
        prob = score / ((double) numPoints * 255.0);
        prob = exp(prob)/2.71828183;
        scoresToWorld(prior, ixrange, iyrange, sx, sy, &x[0], &x[1]);
        x[2] = prior->theta - thetarange + st * thetastep;
        if (fabs(frsm_angle_subtract(x[2],bestResult.theta))>.01)
        continue;
        float * color3fv = frsm_color_util_jet(prob/maxProb);
        bot_lcmgl_color3f(lcmgl, color3fv[0], color3fv[1], color3fv[2]);
        bot_lcmgl_vertex3f(lcmgl, x[0]+1,x[1] ,prob/maxProb);

      }
    }
  }

  bot_lcmgl_end(lcmgl);
  bot_lcmgl_pop_attrib(lcmgl);
  bot_lcmgl_switch_buffer(lcmgl);

#endif

  //check whether we hit the edge of our ranges...

  if (xSat != NULL) {
    *xSat = 0;
    if (bestXInd == 0)
      *xSat = -1;
    else if (bestXInd == ixdim - 1)
      *xSat = 1;
  }

  if (ySat != NULL) {
    *ySat = 0;
    if (bestYInd == 0)
      *ySat = -1;
    else if (bestYInd == iydim - 1)
      *ySat = 1;
  }

  if (thetaSat != NULL) {
    *thetaSat = 0;
    if (bestThetaInd == 0)
      *thetaSat = -1;
    else if (bestThetaInd == itdim - 1)
      *thetaSat = 1;
  }

  //get the number of hits for this transform
  bestResult.hits = getNumHits(points, numPoints, &bestResult, hitThresh);

  //cleanup
  for (int i = 0; i < itdim; i++) {
    free(allScores[i]);
  }
  free(allScores);
  free(bestThetaScores);
  //  free(costMap);
  //  free(tmap);

  //  printf("br.x=%f \t br.y=%f \t br.t=%f\t br.score=%f\n", bestResult.x, bestResult.y, bestResult.theta, bestResult.score);
  return bestResult;
}

LutKernel::LutKernel(double sigma)
{
  //make the square kernel for the endpoints
  kernel_size = computeKernelSize(sigma, 0, 32);
  square_kernel = (uint8_t *) calloc(frsm_sq(kernel_size), sizeof(uint8_t));
  uint8_t * square_kernel_p = square_kernel;
  for (int i = 0; i < kernel_size; i++) {
    for (int j = 0; j < kernel_size; j++) {
      int center = kernel_size / 2;
      double x = i - center;
      double y = j - center;
      double d = sqrt(frsm_sq(x) + frsm_sq(y));
      int val = round(255.0 * exp(-d / sigma));
      square_kernel_p[j * kernel_size + i] = (uint8_t) val;
    }
  }

  //make the line kernels
  line_kernel_stride = computeKernelSize(sigma, M_PI / 4.0, 32);
  slope_step = M_PI / 180.0;
  int num_line_kernels = (45.0 * M_PI / 180.0) / slope_step;
  line_kernels = (uint8_t*) calloc(line_kernel_stride * num_line_kernels, sizeof(uint8_t));
  line_kernel_sizes = (int *) calloc(num_line_kernels, sizeof(int));
  for (int sl = 0; sl < num_line_kernels; sl++) {
    double slope = slope_step * sl;
    line_kernel_sizes[sl] = computeKernelSize(sigma, slope, 32);
    uint8_t * line_kernel_p = line_kernels + (sl * line_kernel_stride);
    for (int i = 0; i < line_kernel_sizes[sl]; i++) {
      int center = line_kernel_sizes[sl] / 2;
      double d = fabs(i - center) * cos(slope);
      int val = round(255.0 * exp(-d / sigma));
      line_kernel_p[i] = (uint8_t) val;
    }
  }
}

LutKernel::~LutKernel()
{
  if (square_kernel != NULL)
    free(square_kernel);
  if (line_kernel_sizes != NULL)
    free(line_kernel_sizes);
  if (line_kernels != NULL)
    free(line_kernels);
}

#ifndef NO_BOT_LCMGL
/**
 * Draw the table via LCMGL in libbot2
 */
void RasterLookupTable::draw_lcmgl(bot_lcmgl_t * lcmgl)
{
  int texid = bot_lcmgl_texture2d(lcmgl, table->data, width, height, width * sizeof(uint8_t),
      BOT_LCMGL_LUMINANCE, BOT_LCMGL_UNSIGNED_BYTE, BOT_LCMGL_COMPRESS_ZLIB);

  bot_lcmgl_enable(lcmgl, 0x0BE2); // #define 0x0BE2 GL_BLEND
  bot_lcmgl_enable(lcmgl, 0x0B71); //#define GL_DEPTH_TEST 0x0B71

  bot_lcmgl_texture_draw_quad(lcmgl, texid,
      x0, y0, 0,
      x0, y1, 0,
      x1, y1, 0,
      x1, y0, 0);

  bot_lcmgl_disable(lcmgl, 0x0BE2); //GL_BLEND
  bot_lcmgl_disable(lcmgl, 0x0B71); //GL_DEPTH_TEST

}
#endif

#ifndef NO_LCM
/**
 * Save the table to a file
 */
void RasterLookupTable::save_to_file(const std::string & channel)
{
  //TODO:
}
/**
 * Load the table from a file
 */
void RasterLookupTable::load_from_file(const std::string & channel)
{
  //TODO:
}

/**
 * Publish the table over LCM
 */
occ_map_pixel_map_t RasterLookupTable::to_lcm_msg()
{

}

/**
 * Fill the table with data from an LCM message
 */
void RasterLookupTable::from_lcm_msg(const occ_map_pixel_map_t & msg)
{
  //TODO:
}

/**
 * Publish the table over LCM
 */
void RasterLookupTable::lcm_publish(lcm_t * lcm, const std::string & channel)
{

}
#endif

}  //namespace frsm
