/*
 * Scan.cpp
 *
 *  Created on: Aug 23, 2009
 *      Author: abachrac
 */
#include "Scan.hpp"

using namespace std;
namespace scanmatch{

Scan::Scan() :
    points(NULL), ppoints(NULL), numPoints(0), utime(-1), laser_type(
            SM_DUMMY_LASER)
{
    memset(&T, 0, sizeof(T));
}

Scan::Scan(int numP, smPoint * points_, ScanTransform T_,
        sm_laser_type_t laser_type_, int64_t utime_, bool buildContours) :
    numPoints(numP), utime(utime_), laser_type(laser_type_)
{
    points = (smPoint *) malloc(numPoints * sizeof(smPoint));
    memcpy(points, points_, numPoints * sizeof(smPoint));
    ppoints = (smPoint *) malloc(numPoints * sizeof(smPoint));
    memset(&T, 0, sizeof(T));
    applyTransform(T_);
    if (buildContours) {
        ContourExtractor * cextractor = new ContourExtractor(laser_type);
        sm_tictoc("findContours");
        cextractor->findContours(ppoints, numPoints, contours);
        sm_tictoc("findContours");
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

}//namespace scanmatch
