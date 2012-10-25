/**
 * SECTION:frsm
 * @title: Scan
 * @short_description: Storage class for keeping track of relevant info associated with a scan.
 * @include: bot/frsm/Scan.hpp
 *
 * These parameters work pretty well with a hokuyo UTM, but others might be better.
 *
 * Linking: -lfrsm
 * namespace: frsm
 */

#ifndef SCAN_H_
#define SCAN_H_

#include "Contour.hpp"
#include <vector>
#include <list>
#include "ScanMatchingUtils.hpp"
#include <float.h>

namespace frsm
{

class Scan
{
public:

    Scan(int numP, frsmPoint * points_, ScanTransform T_,
            frsm_laser_type_t laser_type_, int64_t utime_, bool buildContours =
                    false);
    Scan();
    ~Scan();
    /**
     * number of points that make up this scan.
     */
    unsigned numPoints;
    /**
     * ScanTransform (x,y,theta) of scan origin, presumably to the local frame
     */
    ScanTransform T;
    /**
     * the raw points in the laser's coordinate frame
     */
    frsmPoint * points;
    /**
     * the projected points aftering being transformed by T.
     */
    frsmPoint * ppoints;
    /**
     * timestamp for laser scan.
     */
    int64_t utime;

    /**
     * The type of laser that was used to create this scan. (needed for finding contours)
     */
    frsm_laser_type_t laser_type;
    /**
     *  The set of piecewise linear contours, (projected by T)
     */
    std::vector<Contour*> contours;

    /**
     * applyTransform:
     * @T_: ScanTransform to be applied to the set of points
     *
     * apply the ScanTransform T_ to the raw points in this scan and update
     *  its internal notion of T.
     *
     */
    inline void
    applyTransform(ScanTransform &T_)
    {
        frsm_transformPoints(&T_, points, numPoints, ppoints);
        if (!contours.empty()) {
            double ct_diff = cos(T_.theta - T.theta), st_diff = sin(T_.theta
                    - T.theta);
            for (unsigned int i = 0; i < contours.size(); i++) {
                for (unsigned int j = 0; j < contours[i]->points.size(); j++) {
                    frsmPoint p = contours[i]->points[j];
                    p.x -= T.x;
                    p.y -= T.y;
                    contours[i]->points[j].x = T_.x + ct_diff * p.x - st_diff
                            * p.y;
                    contours[i]->points[j].y = T_.y + st_diff * p.x + ct_diff
                            * p.y;
                }
            }
        }
        T = T_;

    }

};

}

#endif /* SCAN_H_ */
