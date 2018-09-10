// You can use one or the other of the following, but not both!!
/////////////////////////////////////////////////
#ifdef USE_C_TYPE_COMPLEX
#include <complex.h>  // fails w/ std=C++11
#endif
#ifdef USE_CPP_TYPE_COMPLEX
#include <complex>
#endif
/////////////////////////////////////////////////

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <vector>   // C++
#include "utils.h"

#define BIGINT long long

//#define double float   // use to test single-prec speed (only 15% faster).


int main(int argc, char* argv[])
/* C complex type version speed test for mult (plus add) two big 1D arrays.

 See do_complextimings.sh for compilation.

 If you change the array size M, be aware that 8*M doubles of RAM are needed,
 ie 64*M bytes.
*/
{
  BIGINT M=1e8;
  CNTime timer;

  ///////////////////////////////////////////////////////////////////////////////////////
  // Create the random data that will be used throughout
  double *data1_re = (double *)malloc(sizeof(double)*M);
  double *data1_im = (double *)malloc(sizeof(double)*M);
  double *data2_re = (double *)malloc(sizeof(double)*M);
  double *data2_im = (double *)malloc(sizeof(double)*M);
  for(BIGINT i=0;i<M;++i) { data1_re[i] = rand01(); data1_im[i] = rand01(); data2_re[i] = rand01(); data2_im[i] = rand01(); }


  #ifdef USE_C_TYPE_COMPLEX
  ///////////////////////////////////////////////////////////////////////////////////////
  // C-type double
  {
    double *x = (double *)malloc(sizeof(double)*M);
    double *x2 = (double *)malloc(sizeof(double)*M);
    for(BIGINT i=0;i<M;++i) { x[i] = data1_re[i]; x2[i] = data2_re[i]; }
    timer.start();
    for(BIGINT i=1;i<M;++i)     x[i] = x[i-1] + x[i] * x2[i];
    double t=timer.elapsedsec();
    printf("%lld C-type double mults in %.3g s\n\tresult = %.12f\n",M,t,x[M-1]);
    free(x); free(x2);
  }
  

  ///////////////////////////////////////////////////////////////////////////////////////
  // C++ vector double
  {
    std::vector<double> y(M),y2(M);
    for(BIGINT i=0;i<M;++i) { y[i] = data1_re[i]; y2[i] = data2_re[i]; }
    timer.restart();
    for(BIGINT i=1;i<M;++i)     y[i] = y[i-1] + y[i] * y2[i];
    double t=timer.elapsedsec();
    printf("%lld C++ std vector double mults in %.3g s\n\tresult = %.12f\n",M,t,y[M-1]);
    y.clear(); y2.clear();
  }

  ///////////////////////////////////////////////////////////////////////////////////////
  // C-type double complex by hand, separated storage
  {
    double *x = (double *)malloc(sizeof(double)*M); 
    double *x2 = (double *)malloc(sizeof(double)*M);
    double *xi = (double *)malloc(sizeof(double)*M);
    double *xi2 = (double *)malloc(sizeof(double)*M);
    for(BIGINT i=0;i<M;++i) {
      x[i] = data1_re[i]; xi[i] = data1_im[i]; x2[i] = data2_re[i]; xi2[i] = data2_im[i];
    }
    timer.restart();
    for(BIGINT i=1;i<M;++i) {
      double r = x[i]*x2[i] - xi[i]*xi2[i];
      double j = x[i]*xi2[i] + xi[i]*x2[i];
      x[i] = x[i-1]+r;
      xi[i] = xi[i-1]+j;
    }
    double t=timer.elapsedsec();
    printf("%lld complex as doubles by hand (sep storage) mults in %.3g s\n\tresult = %.12f + %.12fI\n",M,t,x[M-1],xi[M-1]);
    free(x); free(xi); free(x2); free(xi2);
  }

  ///////////////////////////////////////////////////////////////////////////////////////
  // C-type double complex by hand
  {
    double *x = (double *)malloc(sizeof(double)*2*M); 
    double *x2 = (double *)malloc(sizeof(double)*2*M);
    for(BIGINT i=0;i<M;++i) {
      x[2*i] = data1_re[i]; x[2*i+1] = data1_im[i]; x2[2*i] = data2_re[i]; x2[2*i+1] = data2_im[i];
    }
    timer.restart();
    for(BIGINT i=2;i<2*M;i+=2) {
      double r = x[i]*x2[i] - x[i+1]*x2[i+1];
      double j = x[i]*x2[i+1] + x[i+1]*x2[i];
      x[i] = x[i-2]+r;
      x[i+1] = x[i-1]+j;
    }
    double t=timer.elapsedsec();
    printf("%lld complex as doubles by hand mults in %.3g s\n\tresult = %.12f + %.12fI\n",M,t,x[2*M-2],x[2*M-1]);
    free(x); free(x2);
  }

  ///////////////////////////////////////////////////////////////////////////////////////
  // C-type complex
  {
    double complex *z = (double complex *)malloc(sizeof(double complex)*M);
    double complex *z2 = (double complex *)malloc(sizeof(double complex)*M);
    for(BIGINT i=0;i<M;++i) {
      z[i] = data1_re[i] + I*data1_im[i];
      z2[i] = data2_re[i] + I*data2_im[i];
    }
    timer.restart();
    for(BIGINT i=1;i<M;++i) {
      z[i] = z[i-1] + z[i] * z2[i];
    }
    double t=timer.elapsedsec();
    printf("%lld C-type complex mults in %.3g s\n\tresult = %.12f + %.12fI\n",M,t,creal(z[M-1]),cimag(z[M-1]));
    free(z); free(z2);
  }
  #endif

  ///////////////////////////////////////////////////////////////////////////////////////
  // C++-type complex
  #ifdef USE_CPP_TYPE_COMPLEX
  {
    const std::complex<double> I(0.0,1.0);
    std::complex<double> *z = (std::complex<double> *)malloc(sizeof(std::complex<double>)*M);
    std::complex<double> *z2 = (std::complex<double> *)malloc(sizeof(std::complex<double>)*M);
    for(BIGINT i=0;i<M;++i) {
      z[i] = data1_re[i] + I*data1_im[i];
      z2[i] = data2_re[i] + I*data2_im[i];
    }
    timer.restart();
    for(BIGINT i=1;i<M;++i) {
      z[i] = z[i-1] + z[i] * z2[i];
    }
    double t=timer.elapsedsec();
    printf("%lld C++-type complex mults in %.3g s\n\tresult = %.12f + %.12fI\n",M,t,z[M-1].real(),z[M-1].imag());
    free(z); free(z2);
  }
  #endif
  
  free(data1_re);
  free(data1_im);
  free(data2_re);
  free(data2_im);

  return 0;
}
