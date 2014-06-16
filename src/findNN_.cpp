#include <algorithm>
#include <cmath>


/*

$HeadURL: http://fgcz-svn.unizh.ch/repos/fgcz/testing/proteomics/R/protViz/src/findNN_.cpp $
$Id: findNN_.cpp 6179 2014-02-27 09:34:04Z cpanse $

*/

extern "C" {
    void findNN_ (int *m_, int *n_, double *q_, double *vec_, int *NN_) {

        size_t dist;
        double d1, d2;

        for (int i = 0; i < *m_; i++){

            
            dist =  std::distance (vec_, std::lower_bound(vec_, vec_ + *n_, q_[i]) );
            NN_[i] = dist ;

            if (dist > 0){
                d1 = std::fabs(q_[i] - vec_[dist - 1]);
                d2 = std::fabs(q_[i] - vec_[dist]);

                if (d1 < d2)
                    NN_[i] = dist - 1 ;
            }
        }
    }
}