/*
 * Scan.cpp
 *
 *  Created on: Aug 23, 2009
 *      Author: abachrac
 */
#include "Scan.hpp"

using namespace std;
namespace frsm{

Scan::Scan() :
    points(NULL), ppoints(NULL), numPoints(0), utime(-1), laser_type(
            SM_DUMMY_LASER)
{
    memset(&T, 0, sizeof(T));
}

Scan::Scan(int numP, frsmPoint * points_, ScanTransform T_,
        frsm_laser_type_t laser_type_, int64_t utime_, bool buildContours) :
    numPoints(numP), utime(utime_), laser_type(laser_type_)
{
    points = (frsmPoint *) malloc(numPoints * sizeof(frsmPoint));
    memcpy(points, points_, numPoints * sizeof(frsmPoint));
    ppoints = (frsmPoint *) malloc(numPoints * sizeof(frsmPoint));
    memset(&T, 0, sizeof(T));
    applyTransform(T_);
    if (buildContours) {
        ContourExtractor * cextractor = new ContourExtractor(laser_type);
        frsm_tictoc("findContours");
        cextractor->findContours(ppoints, numPoints, contours);
        frsm_tictoc("findContours");
        //      s->drawContours();
        delete cextractor;
    }

}

Scan::~Scan()
{
    for (unsigned i = 0; i < contours.size(); i++)
        delete (contours[i]);
    contours.clear();
    if (ppoints != NULL)
        free(ppoints);
    if (points != NULL)
        free(points);
}

}//namespace frsm
