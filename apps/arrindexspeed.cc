/// @file arrindexspeed.cc
///
/// @copyright (c) 2023 CSIRO
/// Australia Telescope National Facility (ATNF)
/// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
/// PO Box 76, Epping NSW 1710, Australia
/// atnf-enquiries@csiro.au
///
/// This file is part of the ASKAP software distribution.
///
/// The ASKAP software distribution is free software: you can redistribute it
/// and/or modify it under the terms of the GNU General Public License as
/// published by the Free Software Foundation; either version 2 of the License,
/// or (at your option) any later version.
///
/// This program is distributed in the hope that it will be useful,
/// but WITHOUT ANY WARRANTY; without even the implied warranty of
/// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
/// GNU General Public License for more details.
///
/// You should have received a copy of the GNU General Public License
/// along with this program; if not, write to the Free Software
/// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA

#include <casacore/casa/Arrays/Cube.h>
#include <casacore/casa/OS/Timer.h>


#include <stdexcept>
#include <iostream>
#include <iomanip>

using namespace std;
using namespace casacore;

void tout( const string& tag, const Timer &timer ) {
    cout << setw(24) << tag
         << "  user: "   << setw(4) << timer.user()
         << "  system: " << setw(4) << timer.system()
         << "  real: "   << setw(4) << timer.real()
         << endl;
}


// Main function
int main(int argc, const char** argv) {
    try {
        const int n = 1300;
        Timer timer;
        Cube<Float> c(n,n,n,0.0f);
        tout("init c",timer);

        timer.mark();
        c = 0.0f;
        tout("c=0",timer);
        timer.mark();
        c = 1.0f;
        tout("c=1",timer);
        timer.mark();
        c.set(1.1f);
        tout("c.set(1.1)",timer);

        timer.mark();
        for (int k=0; k<n; k++) {
            c(Slice(k),Slice(),Slice()) = 2.1f;
        }
        tout("c(k,,)=2.1 (Slice)",timer);

        timer.mark();
        for (int k=0; k<n; k++) {
            c(Slice(),Slice(k),Slice()) = 2.2f;
        }
        tout("c(,k,)=2.2 (Slice)",timer);
        timer.mark();
        for (int k=0; k<n; k++) {
            c(Slice(),Slice(),Slice(k)) = 2.3f;
        }
        tout("c(,,k)=2.3 (Slice)",timer);

        timer.mark();
        for (int k=0; k<n; k++) {
            c.xyPlane(k) = 3.1f;
        }
        tout("C.xyPlane(k)=3.1",timer);

        timer.mark();
        for (int k=0; k<n; k++) {
            c.yzPlane(k) = 3.2f;
        }
        tout("C.yzPlane(k)=3.2",timer);

        timer.mark();
        for (int k=0; k<n; k++) {
            c.xzPlane(k) = 3.3f;
        }
        tout("C.xzPlane(k)=3.3",timer);

        timer.mark();
        for (int k=0; k<n; k++) {
            for (int j = 0; j<n; j++) {
                for (int i=0; i<n; i++) {
                    c(i,j,k) = 4.1f;
                }
            }
        }
        tout("c(i,j,k)=4.1 (loop)",timer);

        timer.mark();
        for (int k=0; k<n; k++) {
            for (int j = 0; j<n; j++) {
                for (int i=0; i<n; i++) {
                    c(j,i,k) = 4.2f;
                }
            }
        }
        tout("c(j,i,k)=4.2 (loop)",timer);

        timer.mark();
        for (int k=0; k<n; k++) {
            for (int j = 0; j<n; j++) {
                for (int i=0; i<n; i++) {
                    c(k,j,i) = 4.3f;
                }
            }
        }
        tout("c(k,j,i)=4.3 (loop)",timer);
    }

    catch (const exception& x) {
        cerr << "Unexpected exception in " << argv[0] << ": " << x.what()
        << endl;
        exit(1);
    }
    exit(0);
}
