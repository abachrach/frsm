#include <stdio.h>
#include <stdlib.h>
#include <signal.h>
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

#include "ConciseArgs.hpp"

using namespace std;
using namespace frsm;

class App {
public:
  lcm_t * lcm;
  ScanMatcher * sm;
  bot_lcmgl_t * lcmgl;
  bool do_drawing;
  bool publish_relative;
  bool publish_pose;
  string lidar_chan;
  string odom_chan;
  string pose_chan;
  frsm_laser_type_t laser_type;
  int scan_skip; //Skip this many scans between processed ones
  int beam_skip; //downsample ranges by only taking 1 out of every beam_skip points
  double spatialDecimationThresh; //don't discard a point if its range is more than this many std devs from the mean range (end of hallway)
  double maxRange; //discard beams with reading further than this value
  float validBeamAngles[2]; //valid part of the field of view of the laser in radians, 0 is the center beam
  frsm_rigid_transform_2d_t prev_odom;
  bool verbose;

  //lcm reading thread stuff
  pthread_t processor_thread;
  pthread_mutex_t lcm_data_mutex;
  pthread_cond_t newLcmData_cv;
  deque<frsm_planar_lidar_t *> * laser_queue;
  int noDrop;

};

////////////////////////////////////////////////////////////////////
//where all the work is done
////////////////////////////////////////////////////////////////////

namespace frsm {

}

void draw(App * app, frsmPoint * points, unsigned numPoints, const ScanTransform * T)
{
  app->sm->draw_state_lcmgl(app->lcmgl);
  app->sm->draw_scan_lcmgl(app->lcmgl, points, numPoints, T);
  bot_lcmgl_switch_buffer(app->lcmgl);
}

static void laser_handler(const lcm_recv_buf_t *rbuf __attribute__((unused)),
    const char * channel __attribute__((unused)), const frsm_planar_lidar_t * msg,
    void * user)
{
  App * app = (App *) user;

  pthread_mutex_lock(&app->lcm_data_mutex);
  app->laser_queue->push_back(frsm_planar_lidar_t_copy(msg));
  pthread_mutex_unlock(&app->lcm_data_mutex);
  pthread_cond_broadcast(&app->newLcmData_cv);
}

static void process_laser(const frsm_planar_lidar_t * msg, void * user __attribute__((unused)))
{
  App * app = (App *) user;
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
    frsm_rigid_transform_2d_t_publish(app->lcm, app->odom_chan.c_str(), &cur_odom);
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
      frsm_rigid_transform_2d_t_publish(app->lcm, app->odom_chan.c_str(), &rel_odom);
    }

  }

  if (app->publish_pose) {
    frsm_pose_t pose;
    memset(&pose, 0, sizeof(pose));
    pose.utime = cur_odom.utime;

    memcpy(pose.pos, cur_odom.pos, 2 * sizeof(double));

    double rpy[3] = { 0, 0, cur_odom.theta };
    bot_roll_pitch_yaw_to_quat(rpy, pose.orientation);

    frsm_pose_t_publish(app->lcm, app->pose_chan.c_str(), &pose);
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

sig_atomic_t still_groovy = 1;

static void sig_action(int signal, siginfo_t *s, void *user)
{
  still_groovy = 0;
}

//dispatcher for new laser data
static void *
processingFunc(void * user)
{
  App * app = (App *) user;
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

  App *app = new App;

  app->lidar_chan = "LASER";
  app->pose_chan = "POSE";
  app->verbose = 0;
  app->do_drawing = 0;
  app->publish_relative = 0;
  app->publish_pose = 0;

  // set to default values
  app->laser_type = FRSM_HOKUYO_UTM;
  app->validBeamAngles[0] = -2.1;
  app->validBeamAngles[1] = 2.1;
  std::stringstream ss;
  ss << app->validBeamAngles[0] << "," << app->validBeamAngles[1];
  string fov_string = ss.str();
  app->beam_skip = 3;
  app->spatialDecimationThresh = .2;
  app->maxRange = 29.7;
  app->laser_queue = new deque<frsm_planar_lidar_t *>();

  bool isUtm = true;
  bool isUrg = false, isSick = false;

  ConciseArgs opt(argc, argv);
  opt.add(app->lidar_chan, "l", "lidar_chan", "Input lidar lcm channel");
  opt.add(app->odom_chan, "o", "odom_chan", "Output odometry lcm channel");
  opt.add(app->publish_relative, "r", "relative", "Publish relative odometry");
  opt.add(app->publish_pose, "p", "publish_pose", "Enable publishing of POSE lcm message");
  opt.add(app->scan_skip, "s", "scan_skip", "Skip this many scans between ones that get processed");
  opt.add(app->noDrop, "n", "no_drop", "Don't drop laser messages if we're getting behind");
  opt.add(app->do_drawing, "d", "do_drawing", "Draw the map and pose using LCMGL");
  opt.add(app->verbose, "v", "verbose");
  opt.addUsageSeperator("\nLow-level options");
  opt.add(isUtm, "T", "UTM", "lidar type is a Hokuyo UTM (best tested/default)");
  opt.add(isUrg, "G", "URG", "lidar type is a Hokuyo URG");
  opt.add(isSick, "S", "SICK", "lidar type is a SICK LMS");

  opt.addUsageSeperator("\nLow-level options (override settings from lidar type)");
  opt.add(app->maxRange, "R", "max_range", "Max range of the lidar\n");
  opt.add(app->beam_skip, "B", "beam_skip", "Skipe every n beams");
  opt.add(app->spatialDecimationThresh, "D", "spatial_decimation", "Spatial decimation threshold in meters");
  opt.add(fov_string, "F", "fov", "Valid portion of the field of view <min,max> in radians");
  opt.parse();

  int numModes = (int) isUtm + (int) isUrg + (int) isSick;
  if (numModes > 1) {
    fprintf(stderr, "Only one lidar type may be specified at once\n");
    opt.usage(true);
  }

  if (opt.wasParsed("fov")) {
    if (sscanf(fov_string.c_str(), "%f,%f", &app->validBeamAngles[0], &app->validBeamAngles[1]) != 2) {
      fprintf(stderr, "Invalid FOV string %s\n", fov_string.c_str());
      opt.usage(true);
    }
  }

  if (isUtm) {
    if (!opt.wasParsed("max_range"))
      app->maxRange = 29.7;
    if (!opt.wasParsed("beam_skip"))
      app->beam_skip = 3;
    if (!opt.wasParsed("fov")) {
      app->validBeamAngles[0] = -2.1;
      app->validBeamAngles[1] = 2.1;
    }
  }
  if (isUrg) {
    app->laser_type = FRSM_HOKUYO_URG;
    if (!opt.wasParsed("max_range"))
      app->maxRange = 4.0;
    if (!opt.wasParsed("beam_skip"))
      app->beam_skip = 3;
    if (!opt.wasParsed("fov")) {
      app->validBeamAngles[0] = -2.1;
      app->validBeamAngles[1] = 2.1;
    }
  }
  if (isSick) {
    app->laser_type = FRSM_SICK_LMS;
    if (!opt.wasParsed("max_range"))
      app->maxRange = 79.0;
    if (!opt.wasParsed("beam_skip"))
      app->beam_skip = 0;
    if (!opt.wasParsed("fov")) {
      app->validBeamAngles[0] = -1.57;
      app->validBeamAngles[1] = 1.57;
    }
  }

  if (app->verbose) {
    printf("INFO: Listening to:%s\n", app->lidar_chan.c_str());
    printf("INFO: Publishing on:%s\n", app->odom_chan.c_str());
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
  frsm_incremental_matching_modes_t matchingMode = FRSM_GRID_COORD; //use gradient descent to improve estimate after brute force search
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

  app->lcm = lcm_create(NULL);
  if (!app->lcm) {
    fprintf(stderr, "ERROR: lcm_create() failed\n");
    return 1;
  }
  if (app->do_drawing)
    app->lcmgl = bot_lcmgl_init(app->lcm, "FRSM");

  //subscribe to lidar messages
  frsm_planar_lidar_t_subscribe(app->lcm, app->lidar_chan.c_str(), laser_handler, app);

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

  return 0;
}

