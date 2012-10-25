/**
 * SECTION:frsm
 * @title: RasterLookupTable
 * @short_description: Local occupancy gridmap for scan matching
 * @include: bot/frsm/RasterLookupTable.hpp
 *
 * Maintain an occupancy gridmap for laser scan-matching.
 *
 * The internal functions for matching a scan are contained in this class.
 *
 *
 * Linking: -lfrsm
 * namespace: frsm
 */

#ifndef RASTERLOOKUPTABLE_H_
#define RASTERLOOKUPTABLE_H_

#include "ScanMatchingUtils.hpp"
#include <stdint.h>
#include <vector>

#include <occ_map/PixelMap.hpp>

#ifndef NO_BOT_LCMGL
#include <bot_core/bot_core.h>
#include <bot_lcmgl_client/lcmgl.h>
#endif

#ifndef NO_LCM
#include <lcm/lcm.h>
#include <lcmtypes/frsm/pixel_map_t.hpp>
#endif

namespace frsm {

class LutKernel {
public:
  LutKernel(double sigma);
  ~LutKernel();
  int kernel_size;
  uint8_t * square_kernel;

  float slope_step;
  int * line_kernel_sizes;
  uint8_t * line_kernels;
  int line_kernel_stride;

};

class RasterLookupTable {
public:
  /**
   * RasterLookupTable:
   * Constructor for a normal table.
   * @x0i: min x position in table in meters
   * @y0i: min y position in table in meters
   * @x1i: max x position in table in meters
   * @y1i: max y position in table in meters
   * @mPP: resolution of the table (meters per pixel)
   * @pixelDivisor: round size of table to ensure that the height/width are a factor of pixelDivisor
   * @initialValue: initialize entries in the table to this value
   *
   */
  RasterLookupTable(double x0i, double y0i, double x1i, double y1i, double mPP, int pixelDivisor = 1,
      uint8_t initialValue = 0);

  /**
   * RasterLookupTable:
   * Constructor for the low resolution table, corresponding to the hi_res table.
   * @hi_res: pointer to the corresponding high resolution table
   * @downsampleFactor: downsample the hi_res table by this factor
   */
  RasterLookupTable(RasterLookupTable * hi_res, int downsampleFactor);
  virtual
  ~RasterLookupTable();

  /**
   * evaluate2D:
   * Perform an exhaustive search over a grid of possible translational transforms of
   *  the scan to find the best one.(searches X,Y)
   * @points: laser points
   * @numPoints: number of laser points
   * @XYT0: center of the search window
   * @prior: initial guess for the aligning position of the scan (normally XTY0==prior, except when
   * performing a multi-resolution search)
   * @ixrange: x-range to be searched in pixels
   * @iyrange: y-range to be searched in pixels
   * @ixdim: x-size of the scores array
   * @iydim: y-size of the scores array
   * @bestScoreIndX: return value for x-index of the best score in the scores array
   * @bestScoreIndY: return value for y-index of the best score in the scores array
   * Returns: the ScanTransform corresponding to the best alignment
   *
   */
  ScanTransform
  evaluate2D(const frsmPoint * points, const unsigned numPoints, const ScanTransform * XYT0,
      const ScanTransform * prior,
      int ixrange, int iyrange, int ixdim, int iydim, float * scores, int * bestScoreIndX, int *bestScoreIndY);

  /**
   * evaluate3D:
   * Perform an exhaustive search over a grid of possible translational and rotational transforms of the
   *  scan to find the best one. (searches X,Y,theta)
   * @points: laser points
   * @numPoints: number of laser points
   * @prior: initial guess for the aligning position of the scan. Used as the center of the window
   * @xrange: x-range to be searched in meters
   * @yrange: y-range to be searched in meters
   * @thetarange: theta-range to be searched in radians
   * @dthetastep: step size to use in theta
   * @hitThresh: threshold for considering a point a "hit"
   * @xSat: return value for whether the best transform was on the x-edge of the search window
   * @ySat: return value for whether the best transform was on the y-edge of the search window
   * @thetaSat: return value for whether the best transform was on the theta-edge of the search window
   * Returns: the ScanTransform corresponding to the best alignment
   *
   */
  ScanTransform
  evaluate3D(const frsmPoint * points, const unsigned numPoints, const ScanTransform * prior, double xrange,
      double yrange, double thetarange, double dthetastep, int hitThresh, int * xSat, int *ySat, int * thetaSat);

  /**
   * evaluate3D_multiRes:
   * Perform a multi-resolution search over the grid of possible translational and rotational transforms of the
   *  scan to find the best one. (searches X,Y,theta)
   * @points: laser points
   * @numPoints: number of laser points
   * @prior: initial guess for the aligning position of the scan. Used as the center of the window
   * @xrange: x-range to be searched in meters
   * @yrange: y-range to be searched in meters
   * @thetarange: theta-range to be searched in radians
   * @dthetastep: step size to use in theta
   * @hitThresh: threshold for considering a point a "hit"
   * @xSat: return value for whether the best transform was on the x-edge of the search window
   * @ySat: return value for whether the best transform was on the y-edge of the search window
   * @thetaSat: return value for whether the best transform was on the theta-edge of the search window
   * Returns: the ScanTransform corresponding to the best alignment
   *
   */
  ScanTransform
  evaluate3D_multiRes(RasterLookupTable * hi_res, const frsmPoint * points, const unsigned numPoints,
      const ScanTransform * prior, double xrange, double yrange, double thetarange, double thetastep, int hitThresh,
      int * xSat, int *ySat, int * thetaSat);

  /**
   * getScore:
   * evaluate the score of the transform XYT0
   */
  float
  getScore(const frsmPoint * points, const unsigned numPoints, const ScanTransform * XYT0);
  /**
   * getScoreDump:
   * evaluate the score of the transform XYT0, and dump individual point scores to a file
   */
  float
  getScoreDump(const frsmPoint * points, const unsigned numPoints, const ScanTransform * XYT0, const char * name);

  /**
   * getNumHits:
   * get number of points that are given a likelihood above the @hitThresh threshold
   */
  int
  getNumHits(const frsmPoint * points, const unsigned numPoints, const ScanTransform * XYT, int hitThresh);

  /**
   * drawRectangle:
   * function used by Ed's method for rendering the lookup table.
   * renders the cost function for a line segment in the rectangular region around the segment
   */
  void
  drawRectangle(double cx, double cy, double x_size, double y_size, double theta, uint8_t *lutSq, int lutSq_size,
      int lutSq_first_zero, double lutSqRange);

  void drawKernel(int ix, int iy, const uint8_t*kernel, int kernel_width, int kernel_height);
  void drawBlurredPoint(const frsmPoint *p, const LutKernel * kern);
  void drawBlurredLine(const frsmPoint *p1, const frsmPoint *p2, const LutKernel * kern);

  /**
   * worldToTable:
   * accessor function that converts from world coordinates in meters to the table coordinates in pixels
   */
  inline void worldToTable(double x, double y, int * ix, int * iy)
  {
    double xy[2] = { x, y };
    int ixy[2];
    table->worldToTable(xy, ixy);
    *ix = ixy[0];
    *iy = ixy[1];
  }

  /**
   * worldToTable:
   * accessor function that converts from table coordinates in pixels to world coordinates in meters
   */
  inline void tableToWorld(int ix, int iy, double * x, double * y)
  {
    int ixy[2] = { ix, iy };
    double xy[2];
    table->tableToWorld(ixy, xy);
    *x = xy[0];
    *y = xy[1];
  }

  /**
   * scoresToWorld:
   * convert from the scoring matrix indices back to world coordinates in meters
   */
  inline void scoresToWorld(const ScanTransform * XYT0, int ixrange, int iyrange, int sx, int sy, double * x,
      double * y)
  {
    int iXYT0[2];
    double dXYTO[2] = { XYT0->x, XYT0->y };
    table->worldToTable(dXYTO, iXYT0);
    int ixy[2];
    ixy[0] = iXYT0[0] - ixrange + sx;
    ixy[1] = iXYT0[1] - iyrange + sy;
    double xy[2];
    table->tableToWorld(ixy, xy);
    *x = xy[0];
    *y = xy[1];
  }

  /**
   * readTable:
   * returns the value of the table at pixel ix,iy
   */
  inline uint8_t readTable(int ix, int iy)
  {
    int ixy[2] = { ix, iy };
    return table->readValue(ixy);
  }

  /**
   * readTable:
   * returns the value of the table at position x,y (in meters)
   */
  inline uint8_t readTable(double x, double y)
  {
    double xy[2] = { x, y };
    return table->readValue(xy);
  }

  /**
   * writeTable:
   * explicitly set the value of pixel ix,iy to v
   */
  inline void writeTable(int ix, int iy, uint8_t v)
  {
    int ixy[2] = { ix, iy };
    table->writeValue(ixy, v);
  }

  /**
   * writeTable:
   * explicitly set the value at position x,y (in meters)to v
   */
  inline void writeTable(double x, double y, uint8_t v)
  {
    double xy[2] = { x, y };
    table->writeValue(xy, v);
  }

#ifndef NO_BOT_LCMGL
  /**
   * Draw the table via LCMGL in libbot2
   */
  void draw_lcmgl(bot_lcmgl_t * lcmgl);

#endif

#ifndef NO_LCM
  /**
   * Save the table to a file
   */
  void save_to_file(const std::string & channel); //TODO: IMPLIMENT ME!
  /**
   * Load the table from a file
   */
  void load_from_file(const std::string & channel); //TODO: IMPLIMENT ME!

  /**
   * Publish the table over LCM
   */
  occ_map_pixel_map_t to_lcm_msg();

  /**
   * Fill the table with data from an LCM message
   */
  void from_lcm_msg(const occ_map_pixel_map_t & msg); //TODO: IMPLIMENT ME!

  /**
   * Publish the table over LCM
   */
  void lcm_publish(lcm_t * lcm, const std::string & channel);

#endif

private:
  occ_map::Uint8PixelMap * table;

  //conveniance variables
  //extremum of the table in meters
  double x0, y0, x1, y1;
  //resolution of the table
  double pixelsPerMeter, metersPerPixel;
  //pixel dimensions
  int width, height;
};
}
#endif /*RASTERLOOKUPTABLE_H_*/
