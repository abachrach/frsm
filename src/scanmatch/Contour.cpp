#include "Contour.hpp"
#include <stdlib.h>
#include <vector>
#include <iostream>
#include <queue>
#include <string>
#include <math.h>
#include <values.h>
#include <assert.h>
#include <opencv/cxcore.h>
#include <opencv/cv.h>

using namespace std;
using namespace scanmatch;

ContourExtractor::ContourExtractor(sm_laser_type_t laser_type)
{
    switch (laser_type)
        {
    case SM_HOKUYO_UTM:
        maxAdjacentDistance = 6;
        maxAdjacentAngularDistance = 3.1 * PI / 180; // radians
        minPointsPerContour = 2;
        alwaysOkayDistance = 0.33; // 0.20;
        maxDistanceRatio = 2;
        maxFirstDistance = 1;
        alwaysAcceptDistance = 0.1;
        searchRange = 10;
        minDistForDistRatio = .08;
        maxDistForSkippingDistRatio = .4;
        simplifyContourThresh = 0;
        break;

    case SM_HOKUYO_URG:
        //These used to work well for a URG, haven't tested in a while.
        maxAdjacentDistance = 0.3;
        maxAdjacentAngularDistance = 10.0 * PI / 180; // radians
        minPointsPerContour = 5;
        alwaysOkayDistance = 0.15; // 0.20;
        maxDistanceRatio = 100;
        maxFirstDistance = .2;
        alwaysAcceptDistance = 0.05;
        searchRange = 8;
        minDistForDistRatio = .08;
        maxDistForSkippingDistRatio = .4;
        simplifyContourThresh = 0;
        break;

    case SM_SICK_LMS:
        //THESE ARE WHAT ED USED... PRESUMABLY WITH A SICK
        maxAdjacentDistance = 5;
        maxAdjacentAngularDistance = 3.1 * PI / 180; // radians
        minPointsPerContour = 1;
        alwaysOkayDistance = 0.33; // 0.20;
        maxDistanceRatio = 1.8;
        maxFirstDistance = 1;
        alwaysAcceptDistance = 0.1;
        searchRange = 4;
        minDistForDistRatio = .08;
        maxDistForSkippingDistRatio = .4;
        simplifyContourThresh = 0;
        break;

    default:
        fprintf(stderr,
                "ERROR: unknown laser type (%d) for ContourExtraction!\n",
                laser_type);
        exit(1);

        }
}

ContourExtractor::~ContourExtractor()
{

}

void
ContourExtractor::findContours(smPoint * points, unsigned numValidPoints,
        std::vector<Contour*> &contours)
{

    priority_queue<Join*, std::vector<Join*>, JoinCompare> joins;

    std::vector<PointRecord> pointrecs;

    int range = searchRange;

    pointrecs.resize(numValidPoints);
    //  int i = 0;
    // create the point record table, one for every point.
    for (unsigned int parent = 0; parent < numValidPoints; parent++) {
        pointrecs[parent].point = points[parent];
    }

    // build the joins...
    // everybody gets their first pick to begin with.
    int numJoins = 0;
    for (unsigned int parent = 0; parent < pointrecs.size(); parent++) {
        Join *j;
        j = findClosestPoint(pointrecs, parent, parent - range, parent - 1);
        if (j != NULL) {
            joins.push(j);
            numJoins++;
        }

        j = findClosestPoint(pointrecs, parent, parent + 1, parent + range);
        if (j != NULL) {
            joins.push(j);
            numJoins++;
        }
    }

    // now we start plucking the best joins. If someone's best choice is
    // taken, we let them
    // pick a new partner and reinsert them into the queue.
    Join j;
    while (joins.size() > 0) {
        Join *jtmp = joins.top();
        joins.pop();
        j = *jtmp;
        delete jtmp;

        //    printf("JOIN cost=%f, parent=%d, victim=%d\n",j.cost,j.parent,j.victim);
        //    continue;

        // victim <--> parent
        if (j.parent > j.victim) {
            // if this parent has already been joined, we're done.
            if (pointrecs[j.parent].left >= 0)
                continue;

            // is the victim still available?
            if (pointrecs[j.victim].right < 0) {
                if (j.cost > alwaysOkayDistance) {
                    if (pointrecs[j.victim].leftDistance > minDistForDistRatio
                            && j.cost / pointrecs[j.victim].leftDistance
                                    > maxDistanceRatio)
                        continue;
                    else if (pointrecs[j.victim].leftDistance
                            < minDistForDistRatio && j.cost
                            > maxDistForSkippingDistRatio)
                        continue;
                    if (pointrecs[j.parent].rightDistance > minDistForDistRatio
                            && j.cost / pointrecs[j.parent].rightDistance
                                    > maxDistanceRatio)
                        continue;
                    else if (pointrecs[j.parent].rightDistance
                            < minDistForDistRatio && j.cost
                            > maxDistForSkippingDistRatio)
                        continue;

                    if (pointrecs[j.victim].leftDistance < 0
                            && pointrecs[j.parent].rightDistance < 0 && j.cost
                            > maxFirstDistance)
                        continue;
                }

                // yes, join.
                pointrecs[j.victim].right = j.parent;
                pointrecs[j.parent].left = j.victim;
                pointrecs[j.parent].leftDistance = j.cost;
                pointrecs[j.victim].rightDistance = j.cost;
            } else {
                // search for a new point.
                jtmp = findClosestPoint(pointrecs, j.parent, j.parent - range,
                        j.parent - 1);
                if (jtmp != NULL)
                    joins.push(jtmp);
            }

            continue;
        }

        // parent <--> victim. Same as above, roles reversed.
        if (j.parent < j.victim) {
            if (pointrecs[j.parent].right >= 0)
                continue;

            if (pointrecs[j.victim].left < 0) {
                if (j.cost > alwaysOkayDistance) {

                    if (pointrecs[j.parent].leftDistance > minDistForDistRatio
                            && j.cost / pointrecs[j.parent].leftDistance
                                    > maxDistanceRatio)
                        continue;
                    else if (pointrecs[j.parent].leftDistance
                            < minDistForDistRatio && j.cost
                            > maxDistForSkippingDistRatio)
                        continue;
                    if (pointrecs[j.victim].rightDistance > minDistForDistRatio
                            && j.cost / pointrecs[j.victim].rightDistance
                                    > maxDistanceRatio)
                        continue;
                    else if (pointrecs[j.victim].rightDistance
                            < minDistForDistRatio && j.cost
                            > maxDistForSkippingDistRatio)
                        continue;

                    if (pointrecs[j.parent].leftDistance < 0
                            && pointrecs[j.victim].rightDistance < 0 && j.cost
                            > maxFirstDistance)
                        continue;
                }

                pointrecs[j.victim].left = j.parent;
                pointrecs[j.parent].right = j.victim;
                pointrecs[j.parent].rightDistance = j.cost;
                pointrecs[j.victim].leftDistance = j.cost;

            } else {
                jtmp = findClosestPoint(pointrecs, j.parent, j.parent + 1,
                        j.parent + range);
                if (jtmp != NULL)
                    joins.push(jtmp);
            }
            continue;
        }

    }

    int contour = 0;

    // pull out the contours.
    for (unsigned int i = 0; i < pointrecs.size(); i++) {
        // we have a new contour.
        if (pointrecs[i].contour < 0) {
            //      if (debug)
            //        System.out.print("contour " + contour + " " + pointrecs[i].left + ": ");

            assert(pointrecs[i].left == -1);

            Contour * c = new Contour();
            int p = i;
            while (p >= 0) {
                //        if (debug)
                //          System.out.print(p + " ");
                c->points.push_back(pointrecs[p].point);
                pointrecs[p].contour = contour;
                p = pointrecs[p].right;
            }

            //      if (debug)
            //        System.out.println("");
            contour++;

            if (c->points.size() == 0) {
                printf("*Adding empty contour!?\n");
            }

            if (c->points.size() >= minPointsPerContour) {
                //TODO: this should be a param
                c->allPoints = c->points;
                if (simplifyContourThresh > 0 && c->points.size() > 2) {
                    //use cvApproxPoly to simplify the contour
                    CvMemStorage* storage = cvCreateMemStorage(0);
                    CvSeq* orig_seq = cvCreateSeq(CV_SEQ_KIND_CURVE | CV_32FC2,
                            sizeof(CvContour), sizeof(CvPoint2D32f), storage);
                    CvSeqWriter writer;
                    cvStartAppendToSeq(orig_seq, &writer);
                    for (unsigned i = 0; i < c->points.size(); i++) {
                        CvPoint2D32f pt =
                            { c->points[i].x, c->points[i].y };
                        CV_WRITE_SEQ_ELEM(pt, writer);
                    }
                    cvEndWriteSeq(&writer);
                    //tolerance of within 2.5cm seems reasonable
                    CvSeq* simple_seq = cvApproxPoly(orig_seq, 0, NULL,
                            CV_POLY_APPROX_DP, simplifyContourThresh);
                    CvSeqReader reader;
                    cvStartReadSeq(simple_seq, &reader, 0);
                    c->points.resize(simple_seq->total);
                    CvPoint2D32f pt;
                    for (int i = 0; i < simple_seq->total; i++) {
                        CV_READ_SEQ_ELEM(pt, reader);
                        c->points[i].x = pt.x;
                        c->points[i].y = pt.y;
                    }
                    //for some reason the last entry is sometimes a duplicate?!?
                    for (unsigned i =0;i<c->points.size()-1;i++){
                      if (sm_dist(&c->points.back(),&c->points[i])<1e-4)
                      {
                        //last one is a duplicate
                        c->points.resize(c->points.size()-1);
                      }
                    }
                    cvReleaseMemStorage(&storage);
		}
		if (c->points.size()>1)
		  contours.push_back(c);
		else
		  delete c;
            } else
                delete c;
        }
    }

    pointrecs.clear();

}

Join *
ContourExtractor::findClosestPoint(std::vector<PointRecord> &pointrecs,
        int parent, int a, int b)
{
    if (a < 0)
        a = 0;
    if (a >= (int) pointrecs.size())
        a = pointrecs.size() - 1;

    if (b < 0)
        b = 0;
    if (b >= (int) pointrecs.size())
        b = pointrecs.size() - 1;

    if (a == parent || b == parent)
        return NULL;

    int bestvictim = -1;
    double bestcost = MAXDOUBLE;

    // how far a span was there in this contour between the last two points?

    // parent better not already have children on both left & right.
    assert(!(pointrecs[parent].left >= 0 && pointrecs[parent].right >= 0));

    smPoint parentp = pointrecs[parent].point;
    for (int i = a; i <= b; i++) {
        PointRecord victim = pointrecs[i];
        if (i == parent)
            continue; //shouldn't happen
        if (victim.left >= 0 && parent < i)
            continue;
        if (victim.right >= 0 && parent > i)
            continue;

        double cost = sm_dist(&parentp, &victim.point);

        if (cost <= alwaysAcceptDistance) {// && i == parent + 1) {
            // stop looking.
            //      bestcost = cost / 4;
            bestcost = cost; //ed had this division by 4... not sure  why.
            bestvictim = i;
            break;
        }

        if (cost < bestcost && cost < maxAdjacentDistance) {
            //      double angularDistance = fabs(sm_angle_subtract(atan2(parentp.y, parentp.x), atan2(victim.point.y, victim.point.x)));
            //
            //      if (angularDistance > maxAdjacentAngularDistance) {
            //        printf("rejected due to Angular distance\n");
            //        continue;
            //      }

            bestcost = cost;
            bestvictim = i;
        }
    }

    // ///////////////////////////////////////////////////
    if (bestvictim < 0)
        return NULL;

    return new Join(parent, bestvictim, bestcost);
}
