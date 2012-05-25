#include <stdio.h>
#include <stdlib.h>
#include <signal.h>
#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif
#include <getopt.h>

#include <lcm/lcm.h>
#include <frsm/FRSM.hpp>
#include <lcmtypes/frsm_rigid_transform_2d_t.h>
#include <lcmtypes/frsm_planar_lidar_t.h>
#include <lcmtypes/frsm_pose_t.h>
#include <lcmtypes/frsm_pixel_map_t.h>
#include <zlib.h>
#include <deque>
#include <pthread.h>
#include <bot_core/bot_core.h>
#include <bot_lcmgl_client/lcmgl.h>
#include <GL/gl.h>

using namespace std;
using namespace frsm;

typedef struct {
  lcm_t * lcm;
  ScanMatcher * sm;
  bot_lcmgl_t * lcmgl;
  int do_drawing;
  int publish_relative;
  int publish_pose;
  char * lidar_chan;
  char * odom_chan;
  char * pose_chan;
  frsm_laser_type_t laser_type;
  int beam_skip; //downsample ranges by only taking 1 out of every beam_skip points
  double spatialDecimationThresh; //don't discard a point if its range is more than this many std devs from the mean range (end of hallway)
  double maxRange; //discard beams with reading further than this value
  float validBeamAngles[2]; //valid part of the field of view of the laser in radians, 0 is the center beam
  frsm_rigid_transform_2d_t prev_odom;
  int verbose;

  //lcm reading thread stuff
  pthread_t processor_thread;
  pthread_mutex_t lcm_data_mutex;
  pthread_cond_t newLcmData_cv;
  deque<frsm_planar_lidar_t *> * laser_queue;
  int noDrop;

} app_t;

////////////////////////////////////////////////////////////////////
//where all the work is done
////////////////////////////////////////////////////////////////////

namespace frsm {

void RasterLookupTable_draw_func(void * user)
{
  app_t * app = (app_t *) user;

  bot_lcmgl_t * lcmgl = app->lcmgl;
  RasterLookupTable * rlt = app->sm->rlt;
  int texid = bot_lcmgl_texture2d(lcmgl, rlt->distdata, rlt->width, rlt->height, rlt->width * sizeof(uint8_t),
      BOT_LCMGL_LUMINANCE, BOT_LCMGL_COMPRESS_ZLIB);

  bot_lcmgl_enable(lcmgl, GL_BLEND);
  bot_lcmgl_enable(lcmgl, GL_DEPTH_TEST);

  bot_lcmgl_texture_draw_quad(lcmgl, texid, rlt->x0, rlt->y0, 0, rlt->x0, rlt->y1, 0, rlt->x1, rlt->y1, 0, rlt->x1,
      rlt->y0, 0);

  bot_lcmgl_disable(lcmgl, GL_BLEND);
  bot_lcmgl_disable(lcmgl, GL_DEPTH_TEST);

}



void ScanMatcher_draw_func(void * user)
{
  app_t * app = (app_t *) user;
  ScanMatcher * sm = app->sm;
  bot_lcmgl_t * lcmgl = app->lcmgl;
  if (sm->useThreads) {
    pthread_mutex_lock(&sm->scans_mutex);
    pthread_mutex_lock(&sm->rlt_mutex);
  }

  //draw the gridmap
  RasterLookupTable_draw_func(app);

  // draw each scan.
  bot_lcmgl_line_width(lcmgl, 2);
  bot_lcmgl_point_size(lcmgl, 4);
  list<Scan *>::iterator it;
  for (it = sm->scans.begin(); it != sm->scans.end(); ++it) {
    Scan * s = *it;
    for (unsigned cidx = 0; cidx < s->contours.size(); cidx++) {
      bot_lcmgl_begin(lcmgl, GL_LINE_STRIP);
      bot_lcmgl_color3f(lcmgl, 0, 1, 0);
      for (unsigned i = 0; i < s->contours[cidx]->points.size(); i++) {
        frsmPoint p0 = s->contours[cidx]->points[i];
        bot_lcmgl_vertex3f(lcmgl, p0.x, p0.y, 0);
      }
      bot_lcmgl_end(lcmgl);
      bot_lcmgl_color3f(lcmgl, 0, 1, 1);
      bot_lcmgl_begin(lcmgl, GL_POINTS); //TODO: Is there a way to do this all at once?
      for (unsigned i = 0; i < s->contours[cidx]->points.size(); i++) {
        frsmPoint p0 = s->contours[cidx]->points[i];
        bot_lcmgl_vertex3f(lcmgl, p0.x, p0.y, 0);
      }
      bot_lcmgl_end(lcmgl);
    }
  }

  if (sm->useThreads) {
    pthread_mutex_unlock(&sm->scans_mutex);
    pthread_mutex_unlock(&sm->rlt_mutex);
  }
}

void * RasterLookupTable_to_msg(void * user)
{
  RasterLookupTable * rlt = (RasterLookupTable *) user;
  frsm_pixel_map_t * msg = (frsm_pixel_map_t*) calloc(1, sizeof(frsm_pixel_map_t));
  msg->xy0[0] = rlt->x0;
  msg->xy0[1] = rlt->y0;
  msg->xy1[0] = rlt->x1;
  msg->xy1[1] = rlt->y1;
  msg->mpp = rlt->metersPerPixel;
  msg->dimensions[0] = rlt->width;
  msg->dimensions[1] = rlt->height;

  int num_cells = rlt->width * rlt->height;
  uLong uncompressed_size = num_cells * sizeof(uint8_t);

  uLong compress_buf_size = uncompressed_size * 1.01 + 12; //with extra space for zlib
  msg->mapData = (uint8_t *) malloc(compress_buf_size);
  int compress_return = compress2((Bytef *) msg->mapData, &compress_buf_size, (Bytef *) rlt->distdata,
      uncompressed_size, Z_BEST_SPEED);
  if (compress_return != Z_OK) {
    fprintf(stderr, "ERROR: Could not compress voxel map!\n");
    exit(1);
  }
  //    fprintf(stderr, "uncompressed_size=%ld compressed_size=%ld\n", uncompressed_size, compress_buf_size);
  msg->datasize = compress_buf_size;
  msg->compressed = 1;

  //set the data_type
  msg->data_type = FRSM_PIXEL_MAP_T_TYPE_UINT8;

  //  msg->utime = utime;

  return msg;
}

}

void draw(app_t * app, frsmPoint * points, unsigned numPoints, const ScanTransform * T)
{
  ScanMatcher_draw_func(app);

  //draw the current scan
  bot_lcmgl_t * lcmgl = app->lcmgl;
  bot_lcmgl_point_size(lcmgl, 4);
  bot_lcmgl_color3f(lcmgl, 1, 0, 0);
  Scan * s = new Scan(numPoints, points, *T, SM_DUMMY_LASER, 0, false);
  bot_lcmgl_begin(lcmgl, GL_POINTS);
  for (int i = 0; i < numPoints; i++) {
    bot_lcmgl_vertex3f(lcmgl, s->ppoints[i].x, s->ppoints[i].y, 0);
  }
  bot_lcmgl_end(lcmgl);
  delete s;

  //draw the robot location

  bot_lcmgl_line_width(lcmgl,4);
  double xyz[3] = { T->x, T->y, 0 };
  bot_lcmgl_color3f(lcmgl, 0, 0, 1);
  bot_lcmgl_circle(lcmgl, xyz, .5);
  double heading[3] = { T->x + cos(T->theta), T->y + sin(T->theta), 0 };
  bot_lcmgl_begin(lcmgl, GL_LINES);
  bot_lcmgl_color3f(lcmgl, 1, 1, 0);
  bot_lcmgl_vertex3f(lcmgl, xyz[0], xyz[1], xyz[2]);
  bot_lcmgl_vertex3f(lcmgl, heading[0], heading[1], heading[2]);
  bot_lcmgl_end(lcmgl);

  bot_lcmgl_switch_buffer(app->lcmgl);
}

static void laser_handler(const lcm_recv_buf_t *rbuf __attribute__((unused)), const char * channel __attribute__((unused)), const frsm_planar_lidar_t * msg,
    void * user  __attribute__((unused)))
{
  app_t * app = (app_t *) user;
  pthread_mutex_lock(&app->lcm_data_mutex);
  app->laser_queue->push_back(frsm_planar_lidar_t_copy(msg));
  pthread_mutex_unlock(&app->lcm_data_mutex);
  pthread_cond_broadcast(&app->newLcmData_cv);
}

static void process_laser(const frsm_planar_lidar_t * msg, void * user  __attribute__((unused)))
{
  app_t * app = (app_t *) user;
  frsm_tictoc("process_laser");
  frsm_tictoc("recToSend");

  ////////////////////////////////////////////////////////////////////
  //Project ranges into points, and decimate points so we don't have too many
  ////////////////////////////////////////////////////////////////////
  frsmPoint * points = (frsmPoint *) calloc(msg->nranges, sizeof(frsmPoint));
  int numValidPoints = frsm_projectRangesAndDecimate(app->beam_skip, app->spatialDecimationThresh, msg->ranges,
      msg->nranges, msg->rad0, msg->radstep, points, app->maxRange, app->validBeamAngles[0], app->validBeamAngles[1]);
  if (numValidPoints < 30) {
    fprintf(stderr, "WARNING! NOT ENOUGH VALID POINTS! numValid=%d\n", numValidPoints);
    return;
  }

  ////////////////////////////////////////////////////////////////////
  //Actually do the matching
  ////////////////////////////////////////////////////////////////////
  ScanTransform r = app->sm->matchSuccessive(points, numValidPoints, app->laser_type, msg->utime, NULL); //don't have a better estimate than prev, so just set prior to NULL
  //utime is ONLY used to tag the scans that get added to the map, doesn't actually matter

  ////////////////////////////////////////////////////////////////////
  //Publish
  ////////////////////////////////////////////////////////////////////
  //publish integrated absolute position instead of delta to last scan
  frsm_rigid_transform_2d_t odom;
  odom.utime = msg->utime;
  memset(&odom, 0, sizeof(odom));
  odom.pos[0] = r.x;
  odom.pos[1] = r.y;
  odom.theta = r.theta;

  //publish integrated absolute position instead of delta to last scan
  frsm_rigid_transform_2d_t cur_odom;
  memset(&cur_odom, 0, sizeof(cur_odom));
  cur_odom.utime = msg->utime;
  cur_odom.pos[0] = r.x;
  cur_odom.pos[1] = r.y;
  cur_odom.theta = r.theta;
  memcpy(cur_odom.cov, r.sigma, 9 * sizeof(double));

  if (!app->publish_relative)
    frsm_rigid_transform_2d_t_publish(app->lcm, app->odom_chan, &cur_odom);
  else {
    //compute the relative odometry
    if (app->prev_odom.utime > 0) {
      frsm_rigid_transform_2d_t rel_odom;
      rel_odom.utime = cur_odom.utime;
      rel_odom.utime_prev = app->prev_odom.utime;
      double delta[2];
      frsm_vector_sub_2d(cur_odom.pos, app->prev_odom.pos, delta);

      frsm_rotate2D(delta, -cur_odom.theta, rel_odom.pos);
      rel_odom.theta = frsm_angle_subtract(cur_odom.theta, app->prev_odom.theta);
      //rotate the covariance estimate to body frame
      frsm_rotateCov2D(cur_odom.cov, -cur_odom.theta, rel_odom.cov);
      frsm_rigid_transform_2d_t_publish(app->lcm, app->odom_chan, &rel_odom);
    }

  }

  if (app->publish_pose) {
    frsm_pose_t pose;
    memset(&pose, 0, sizeof(pose));
    pose.utime = cur_odom.utime;

    memcpy(pose.pos, cur_odom.pos, 2 * sizeof(double));

    double rpy[3] = { 0, 0, cur_odom.theta };
    bot_roll_pitch_yaw_to_quat(rpy, pose.orientation);

    frsm_pose_t_publish(app->lcm, app->pose_chan, &pose);
  }
  frsm_tictoc("recToSend");

  app->prev_odom = cur_odom;

  ////////////////////////////////////////////////////////////////////
  //Print current position periodically!
  ////////////////////////////////////////////////////////////////////
  static double lastPrintTime = 0;
  if (frsm_get_time() - lastPrintTime > 2.0) {
    lastPrintTime = frsm_get_time();
    //print out current state
    fprintf(stderr, "x=%+7.3f y=%+7.3f t=%+7.3f\t score=%f hits=%.2f sx=%.2f sxy=%.2f sy=%.2f st=%.2f, numValid=%d\n",
        r.x, r.y, r.theta, r.score, (double) r.hits / (double) numValidPoints, r.sigma[0], r.sigma[1], r.sigma[4],
        r.sigma[8], numValidPoints);
  }

  ////////////////////////////////////////////////////////////////////
  //Do drawing periodically!
  ////////////////////////////////////////////////////////////////////
  static double lastDrawTime = 0;
  if (app->do_drawing && frsm_get_time() - lastDrawTime > .2) {
    lastDrawTime = frsm_get_time();
    frsm_tictoc("drawing");
    draw(app, points, numValidPoints, &r);
    frsm_tictoc("drawing");
  }

  ////////////////////////////////////////////////////////////////////
  //cleanup!
  ////////////////////////////////////////////////////////////////////

  free(points);
  frsm_tictoc("process_laser");

}

void app_destroy(app_t *app)
{
  // dump timing stats
  frsm_tictoc(NULL);

  if (app) {
    if (app->lidar_chan)
      free(app->lidar_chan);
    if (app->odom_chan)
      free(app->odom_chan);
    if (app->pose_chan)
      free(app->pose_chan);

    lcm_destroy(app->lcm);
    //TODO: kill thread

    free(app);
  }

}

static void usage(const char *name)
{
  fprintf(stderr, "usage: %s [options]\n"
    "\n"
    "  -h, --help                      Shows this help text and exits\n"
    "  -l, --lidar <LCM CHANNEL>       Input lcm channel default:\"LASER\"\n"
    "  -o, --odometry <LCM CHANNEL>    Output odometry channel default:\"<lidar chan>_ODOMETRY\"\n"
    "  -r, --relative                  Publish relative\n"
    "  -p, --pose   <LCM CHANNEL>      Publish pose message for easy viewing in viewer default:\"POSE\"\n"
    "  -d, --draw                      Show window with scan matches \n"
    "  -n, --nodrop                    don't drop laser messages if we're getting behind \n"
    "  -v, --verbose                   Be verbose\n"
    "  -m, --mode  \"HOKUYO_UTM\"|\"SICK\" configures low-level options.\n"
    "\n"
    "Low-level options:\n"
    "  -M, --mask <min,max>            Mask min max angles in (radians)\n"
    "  -B, --beamskip <n>              Skip every n beams \n"
    "  -D, --decimation <value>        Spatial decimation threshold (meters?)\n"
    "  -R, --range <range>             Maximum range (meters)\n", name);
}

sig_atomic_t still_groovy = 1;

static void sig_action(int signal, siginfo_t *s, void *user)
{
  still_groovy = 0;
}

//dispatcher for new laser data
static void *
processingFunc(void * user)
{
  app_t * app = (app_t *) user;
  frsm_planar_lidar_t * local_laser = NULL;
  pthread_mutex_lock(&app->lcm_data_mutex);
  while (1) {
    if (app->laser_queue->empty()) {
      pthread_cond_wait(&app->newLcmData_cv, &app->lcm_data_mutex);
      continue;
    }

    //there is new data to be processed
    //copy shared data to local storage
    frsm_planar_lidar_t * laser_msg;
    if (app->noDrop)
      laser_msg = app->laser_queue->front();
    else
      laser_msg = app->laser_queue->back();

    if (local_laser != NULL)
      frsm_planar_lidar_t_destroy(local_laser);
    local_laser = frsm_planar_lidar_t_copy(laser_msg);

    //process the data
    pthread_mutex_unlock(&app->lcm_data_mutex);
    process_laser(local_laser, (void *) app);
    pthread_mutex_lock(&app->lcm_data_mutex);

    //remove data from the queue
    int numRemoved = 0;
    while (!app->laser_queue->empty() && app->laser_queue->front()->utime <= local_laser->utime) {
      frsm_planar_lidar_t_destroy(app->laser_queue->front());
      app->laser_queue->pop_front();
      numRemoved++;
    }
    if (numRemoved > 1) {
      fprintf(stderr, "dropped %d laser messages\n", numRemoved - 1);
    }
  }
  return NULL;
}

int main(int argc, char *argv[])
{
  setlinebuf(stdout);

  app_t *app = (app_t *) calloc(1, sizeof(app_t));

  app->lidar_chan = strdup("LASER");
  app->pose_chan = strdup("POSE");
  app->verbose = 0;
  app->do_drawing = 0;
  app->publish_relative = 0;
  app->publish_pose = 0;

  // set to default values
  app->laser_type = SM_HOKUYO_UTM;
  //TODO:these should be read from command line or something...
  //parameters for a hokuyo with the helicopters mirror's attached
  app->validBeamAngles[0] = -2.1;
  app->validBeamAngles[1] = 2.1;
  app->beam_skip = 3;
  app->spatialDecimationThresh = .2;
  app->maxRange = 29.7;
  app->laser_queue = new deque<frsm_planar_lidar_t *> ();

  const char *optstring = "hl:o:drvm:R:B:D:M:p::n:";
  char c;
  struct option long_opts[] = { { "help", no_argument, 0, 'h' },
      { "lidar", required_argument, 0, 'l' },
      { "odometry", required_argument, 0, 'o' },
      { "draw", no_argument, 0, 'd' },
      { "relative", no_argument, 0, 'r' },
      { "nodrop", no_argument, 0, 'n' },
      { "mode", required_argument, 0, 'm' },
      { "range", required_argument, 0, 'R' },
      { "beamskip", required_argument, 0, 'B' },
      { "decimation", required_argument, 0, 'D' },
      { "mask", required_argument, 0, 'M' },
      { "pose", optional_argument, 0, 'p' },
      { "verbose", no_argument, 0, 'v' },
      { 0, 0, 0, 0 } };

  while ((c = getopt_long(argc, argv, optstring, long_opts, 0)) >= 0) {
    switch (c) {
    case 'l':
      free(app->lidar_chan);
      app->lidar_chan = strdup(optarg);
      break;
    case 'o':
      app->odom_chan = strdup(optarg);
      break;
    case 'r':
      app->publish_relative = 1;
      printf("INFO: Publish relative enabled\n");
      break;
    case 'p':
      app->publish_pose = 1;
      printf("INFO: Publish pose message enabled\n");
      if (optarg != NULL) {
        app->pose_chan = strdup(optarg);
      }
      break;
    case 'n':
      app->noDrop = 1;
      break;
    case 'v':
      app->verbose = 1;
      break;
    case 'd':
      app->do_drawing = 1;
      printf("INFO: Drawing enabled\n");

      break;
    case 'm':
      if (optarg && optarg[0] == 'S') {
        // assume we want sick mode:
        app->laser_type = SM_SICK_LMS;
        app->maxRange = 79.0;
        app->beam_skip = 0;
        //                app->spatialDecimationThresh =0;
        // anything else required?
      }
      break;
    case 'R':
      app->maxRange = strtod(optarg, 0);
      break;
    case 'B':
      app->beam_skip = atoi(optarg);
      break;
    case 'D':
      app->spatialDecimationThresh = strtod(optarg, 0);
      break;
    case 'M':
      sscanf(optarg, "%f,%f", &app->validBeamAngles[0], &app->validBeamAngles[1]);
      break;
    case 'h':
    default:
      usage(argv[0]);
      return 1;
    }
  }

  if (!app->odom_chan) {
    int max_len = strlen(app->lidar_chan) + 16;
    //printf("max_len:%d\n",max_len);
    app->odom_chan = (char *) calloc(1, max_len * sizeof(char));
    snprintf(app->odom_chan, max_len, "%s_ODOMETRY", app->lidar_chan);
  }

  if (app->verbose) {
    printf("INFO: Listening to:%s\n", app->lidar_chan);
    printf("INFO: Publishing on:%s\n", app->odom_chan);
    printf("INFO: Do Draw:%d\n", app->do_drawing);
    printf("INFO: Publish relative:%d\n", app->publish_relative);
    printf("INFO: Max Range:%lf\n", app->maxRange);
    printf("INFO: SpatialDecimationThresh:%lf\n", app->spatialDecimationThresh);
    printf("INFO: Beam Skip:%d\n", app->beam_skip);
    printf("INFO: validRange:%f,%f\n", app->validBeamAngles[0], app->validBeamAngles[1]);
  }

  //initialize tictoc for threading
  frsm_tictoc_init();

  //hardcoded scan matcher params
  double metersPerPixel = .02; //translational resolution for the brute force search
  double thetaResolution = .01; //angular step size for the brute force search
  frsm_incremental_matching_modes_t matchingMode = SM_COORD_ONLY; //use gradient descent to improve estimate after brute force search
  int useMultires = 3; // low resolution will have resolution metersPerPixel * 2^useMultiRes

  double initialSearchRangeXY = .15; //nominal range that will be searched over
  double initialSearchRangeTheta = .1;

  //SHOULD be set greater than the initialSearchRange
  double maxSearchRangeXY = .3; //if a good match isn't found I'll expand and try again up to this size...
  double maxSearchRangeTheta = .2; //if a good match isn't found I'll expand and try again up to this size...


  int maxNumScans = 30; //keep around this many scans in the history
  double addScanHitThresh = .80; //add a new scan to the map when the number of "hits" drops below this

  bool stationaryMotionModel = false; //use constant velocity model
  double motionModelPriorWeight = 0; //don't use the prior for anything other than centering the window

  int useThreads = 1;

  //create the actual scan matcher object
  app->sm = new ScanMatcher(metersPerPixel, thetaResolution, useMultires, useThreads, true);

  if (app->sm->isUsingIPP())
    fprintf(stderr, "Using IPP\n");
  else
    fprintf(stderr, "NOT using IPP\n");

  ScanTransform startPose;
  memset(&startPose, 0, sizeof(startPose));
  app->sm->initSuccessiveMatchingParams(maxNumScans, initialSearchRangeXY, maxSearchRangeXY, initialSearchRangeTheta,
      maxSearchRangeTheta, matchingMode, addScanHitThresh, stationaryMotionModel, motionModelPriorWeight, &startPose);

  //setup lcm reading thread
  pthread_mutex_init(&app->lcm_data_mutex, NULL);
  pthread_cond_init(&app->newLcmData_cv, NULL);
  //create processing thread
  pthread_create(&app->processor_thread, 0, (void *
  (*)(void *)) processingFunc, (void *) app);

  /* LCM */
  app->lcm = lcm_create(NULL);
  app->lcmgl = bot_lcmgl_init(app->lcm, "FRSM");
  if (!app->lcm) {
    fprintf(stderr, "ERROR: lcm_create() failed\n");
    return 1;
  }

  frsm_planar_lidar_t_subscribe(app->lcm, app->lidar_chan, laser_handler, app);

  // setup sigaction();
  struct sigaction new_action;
  new_action.sa_sigaction = sig_action;
  sigemptyset(&new_action.sa_mask);
  new_action.sa_flags = 0;

  sigaction(SIGINT, &new_action, NULL);
  sigaction(SIGTERM, &new_action, NULL);
  sigaction(SIGKILL, &new_action, NULL);
  sigaction(SIGHUP, &new_action, NULL);

  /* sit and wait for messages */
  while (still_groovy)
    lcm_handle(app->lcm);

  app_destroy(app);

  return 0;
}

