
#include "casacore/casa/OS/Timer.h"
#include "casacore/casa/BasicSL/Complex.h"

#include <iostream>
#include <iomanip>
#include <string>
#include <vector>

using namespace std;

void tout( const std::string tag, const casa::Timer &timer ) {
    cout << setw(24) << tag
         << "  user: "   << setw(4) << timer.user()
         << "  system: " << setw(4) << timer.system()
         << "  real: "   << setw(4) << timer.real()
         << endl;
}

int main(int argc, const char** argv)
{

    const int N = int(2e8);
    casa::Timer timer;
    cout << left;

    timer.mark();
    std::vector<std::complex<float> > a(N), b(N);
    for ( int i=0; i<N; i++ ) {
        a[i] = std::complex<float>(i,i);
        b[i] = std::complex<float>(1.0/(i+1),-0.5/(i+1));
    }
    tout( "init:", timer );

    cout << "std::complex<float>:" << endl;
    std::complex<float> y1 = 0;
    {
        timer.mark();
        std::complex<float> x = 0;
        for ( int i=0; i<N; i++ ) {
            x = x + a[i]*b[i];
        }
        y1 += x;
        tout( "x = x + a*b", timer );
    }
    {
        timer.mark();
        std::complex<float> x = 0;
        for ( int i=0; i<N; i++ ) {
            x += a[i]*b[i];
        }
        y1 += x;
        tout( "x += a*b", timer );
    }

    std::complex<float> y2 = 0;
    {
        timer.mark();
        std::complex<float> x = 0;
        for ( int i=0; i<N; i++ ) {
            x = x + std::complex<float>( a[i].real()*b[i].real() - a[i].imag()*b[i].imag(),
                                         a[i].real()*b[i].imag() + a[i].imag()*b[i].real() );
        }
        y2 += x;
        tout( "x = x + rr-ii,ri+ir", timer );
    }
    {
        timer.mark();
        std::complex<float> x = 0;
        for ( int i=0; i<N; i++ ) {
            x += std::complex<float>( a[i].real()*b[i].real() - a[i].imag()*b[i].imag(),
                                      a[i].real()*b[i].imag() + a[i].imag()*b[i].real() );
        }
        y2 += x;
        tout( "x += rr-ii,ri+ir", timer );
    }

/*
    float k1, k2, k3;
    std::complex<float> y3 = 0;
    {
        timer.mark();
        std::complex<float> x = 0;
        for ( int i=0; i<N; i++ ) {
            k1 = a[i].real()*b[i].real();
            k2 = a[i].imag()*b[i].imag();
            k3 = (a[i].real() + a[i].imag()) * (b[i].real() + b[i].imag());
            x = x + std::complex<float>( k1 - k2, k3 - k1 - k2 );
        }
        y3 += x;
        tout( "x = x + k1-k2,k1+k3", timer );
    }
    {
        timer.mark();
        std::complex<float> x = 0;
        for ( int i=0; i<N; i++ ) {
            k1 = a[i].real()*b[i].real();
            k2 = a[i].imag()*b[i].imag();
            k3 = (a[i].real() + a[i].imag()) * (b[i].real() + b[i].imag());
            x += std::complex<float>( k1 - k2, k3 - k1 - k2 );
        }
        y3 += x;
        tout( "x += k1-k2,k1+k3", timer );
    }
*/

    cout << "casa::Complex:" << endl;
    casa::Complex z1 = 0;
    {
        timer.mark();
        casa::Complex x = 0;
        for ( int i=0; i<N; i++ ) {
            x = x + a[i]*b[i];
        }
        z1 += x;
        tout( "x = x + a*b", timer );
    }
    {
        timer.mark();
        casa::Complex x = 0;
        for ( int i=0; i<N; i++ ) {
            x += a[i]*b[i];
        }
        z1 += x;
        tout( "x += a*b", timer );
    }

    casa::Complex z2 = 0;
    {
        timer.mark();
        casa::Complex x = 0;
        for ( int i=0; i<N; i++ ) {
            x = x + casa::Complex( a[i].real()*b[i].real() - a[i].imag()*b[i].imag(),
                                   a[i].real()*b[i].imag() + a[i].imag()*b[i].real() );
        }
        z2 += x;
        tout( "x = x + rr-ii,ri+ir", timer );
    }
    {
        timer.mark();
        casa::Complex x = 0;
        for ( int i=0; i<N; i++ ) {
            x += casa::Complex( a[i].real()*b[i].real() - a[i].imag()*b[i].imag(),
                                a[i].real()*b[i].imag() + a[i].imag()*b[i].real() );
        }
        z2 += x;
        tout( "x += rr-ii,ri+ir", timer );
    }

    cout << endl;
    cout << "results:" << endl;
    cout << y1 << endl;
    cout << y2 << endl;
    cout << z1 << endl;
    cout << z2 << endl;
    cout << endl;

    return(0);

}

