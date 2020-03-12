
// demo of manual SIMD via vectorclass library, edited slightly from Libin Lu's
// ccm_simd_demo github repo of 12/10/19.
// Coords are in 3D, and kernel is 1/r (not the exp of abs(r) in 1D that Libin used).

// Note VCL version 1 includes rsqrt in single prec (eg 8f, 16f), but not double that we test here. Hence we added custom double rsqrt to the VCL in this repo.

// TO DO: single prec tests.

// Alex Barnett, Libin Lu   3/11/20

// Compile & run: (better: use GCC 7 or 9)

// g++ -fPIC -g -O3 -march=native -funroll-loops -fopenmp -std=c++17 -DVCL -I./VCL_AgnerFog_version1+rsqrt -fopenmp lap3dkernel.cpp; ./a.out

// Libin's old compilation commands:
// icpc -fPIC -g -O3 -march=native -funroll-loops -fopenmp -std=c++17 -DVCL -I./version1 main.cc; ./a.out
// g++ -fPIC -g -O3 -march=native -funroll-loops -fopenmp -std=c++17 -DVCL -I./version1 main.cc; ./a.out
// /cm/shared/sw/pkg/vendor/intel-pstudio/2019-3/vtune_amplifier/bin64/amplxe-cl -collect hotspots ./a.out

#include <omp.h>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <vector>

// choose double prec throughout
typedef double Real;
typedef long Integer;

void test_kernel_ts(const Real* Xs, const Real* Ys, const Real* Zs,
                    int ns, const Real* charge,
                    const Real* Xt, const Real* Yt, const Real* Zt, int nt, Real* pot){
#pragma omp parallel for
    for(int t = 0; t < nt; t++){
        for(int s = 0; s < ns; s++){

            Real dx = Xs[s] - Xt[t];
            Real r2 = dx*dx;
            dx = Ys[s] - Yt[t];
            r2 += dx*dx;
            dx = Zs[s] - Zt[t];
            r2 += dx*dx;
            pot[t] += charge[s] / sqrt(r2);

        }
    }
}

void test_kernel_st(const Real* Xs, const Real* Ys, const Real* Zs,
                    int ns, const Real* charge,
                 const Real* Xt, const Real* Yt, const Real* Zt, int nt, Real* pot){
  for(int s = 0; s < ns; s++){
    // note we can't omp over sources without atomic update to pot[t] which is slow
#pragma omp parallel for
    for(int t = 0; t < nt; t++){

            Real dx = Xs[s] - Xt[t];
            Real r2 = dx*dx;
            dx = Ys[s] - Yt[t];
            r2 += dx*dx;
            dx = Zs[s] - Zt[t];
            r2 += dx*dx;
            //             #pragma omp atomic update   // if want out omp loop, but slows
            pot[t] += charge[s] / sqrt(r2);

        }
    }
}

#ifdef VCL
// Use Agner Fog's vector class library
#include "vectorclass.h"
#ifdef __INTEL_COMPILER
#include "vectormath_lib.h"
#else
#include "vectormath_exp.h"
#endif

void test_kernel_vec512(const Real* Xs, const Real* Ys, const Real* Zs,
                    int ns, const Real* charge,
                 const Real* Xt, const Real* Yt, const Real* Zt, int nt, Real* pot){
#pragma omp parallel for
    for(int t = 0; t < nt; t += 8){

      Vec8d Xt_vec, Yt_vec, Zt_vec, pot_vec;
        // load
        Xt_vec.load(Xt + t);
        Yt_vec.load(Yt + t);
        Zt_vec.load(Zt + t);
        pot_vec.load(pot + t);

        for(int s = 0; s < ns; s++){
            Vec8d dx = Vec8d(Xs[s]) - Xt_vec;
            Vec8d dy = Vec8d(Ys[s]) - Yt_vec;
            Vec8d dz = Vec8d(Zs[s]) - Zt_vec;
            Vec8d r2 = dx*dx+dy*dy+dz*dz;
            pot_vec += Vec8d(charge[s]) / sqrt(r2);
        }

        // store
        pot_vec.store(pot + t);
    }
}

void test_kernel_vec512_fast_rsqrt(const Real* Xs, const Real* Ys, const Real* Zs,
                    int ns, const Real* charge,
                 const Real* Xt, const Real* Yt, const Real* Zt, int nt, Real* pot){
#pragma omp parallel for
    for(int t = 0; t < nt; t += 8){

      Vec8d Xt_vec, Yt_vec, Zt_vec, pot_vec;
        // load
        Xt_vec.load(Xt + t);
        Yt_vec.load(Yt + t);
        Zt_vec.load(Zt + t);
        pot_vec.load(pot + t);

        for(int s = 0; s < ns; s++){
            Vec8d dx = Vec8d(Xs[s]) - Xt_vec;
            Vec8d dy = Vec8d(Ys[s]) - Yt_vec;
            Vec8d dz = Vec8d(Zs[s]) - Zt_vec;
            Vec8d r2 = dx*dx+dy*dy+dz*dz;

            // Newton iterations for rsqrt double precision
            Vec8d rinv = approx_rsqrt(r2);
            // Note the trick here changes the Newton iteration constants
            // the standard two Newton iterations are,
            // rinv *= ((3.0) - r2 * rinv * rinv) * 0.5;
            // rinv *= ((3.0) - r2 * rinv * rinv) * 0.5;
            // we can save one simd multiplication " * 0.5 " by using the following trick
            // first Newton iteration
            rinv *= ((3.0) - r2 * rinv * rinv);
            // second Newton iteration
            rinv *= ((12.0) - r2 * rinv * rinv) * (0.0625);

            pot_vec += Vec8d(charge[s]) * rinv;
        }

        // store
        pot_vec.store(pot + t);
    }
}

void test_kernel_vec256(const Real* Xs, const Real* Ys, const Real* Zs,
                    int ns, const Real* charge,
                 const Real* Xt, const Real* Yt, const Real* Zt, int nt, Real* pot){
#pragma omp parallel for
    for(int t = 0; t < nt; t += 4){

      Vec4d Xt_vec, Yt_vec, Zt_vec, pot_vec;
        // load
        Xt_vec.load(Xt + t);
        Yt_vec.load(Yt + t);
        Zt_vec.load(Zt + t);
        pot_vec.load(pot + t);

        for(int s = 0; s < ns; s++){
            Vec4d dx = Vec4d(Xs[s]) - Xt_vec;
            Vec4d dy = Vec4d(Ys[s]) - Yt_vec;
            Vec4d dz = Vec4d(Zs[s]) - Zt_vec;
            Vec4d r2 = dx*dx+dy*dy+dz*dz;
            pot_vec += Vec4d(charge[s]) / sqrt(r2);
        }

        // store
        pot_vec.store(pot + t);
    }
}

void test_kernel_vec256_fast_rsqrt(const Real* Xs, const Real* Ys, const Real* Zs,
                    int ns, const Real* charge,
                 const Real* Xt, const Real* Yt, const Real* Zt, int nt, Real* pot){
#pragma omp parallel for
    for(int t = 0; t < nt; t += 4){

      Vec4d Xt_vec, Yt_vec, Zt_vec, pot_vec;
        // load
        Xt_vec.load(Xt + t);
        Yt_vec.load(Yt + t);
        Zt_vec.load(Zt + t);
        pot_vec.load(pot + t);

        for(int s = 0; s < ns; s++){
            Vec4d dx = Vec4d(Xs[s]) - Xt_vec;
            Vec4d dy = Vec4d(Ys[s]) - Yt_vec;
            Vec4d dz = Vec4d(Zs[s]) - Zt_vec;
            Vec4d r2 = dx*dx+dy*dy+dz*dz;

            // Newton iterations for rsqrt double precision
            Vec4d rinv = approx_rsqrt(r2);
            // Note the trick here changes the Newton iteration constants
            // the standard two Newton iterations are,
            // rinv *= ((3.0) - r2 * rinv * rinv) * 0.5;
            // rinv *= ((3.0) - r2 * rinv * rinv) * 0.5;
            // we can save one simd multiplication " * 0.5 " by using the following trick
            // first Newton iteration
            rinv *= ((3.0) - r2 * rinv * rinv);
            // second Newton iteration
            rinv *= ((12.0) - r2 * rinv * rinv) * (0.0625);

            pot_vec += Vec4d(charge[s]) * rinv;
        }

        // store
        pot_vec.store(pot + t);
    }
}
#endif   // VCL

// Return random real-valued standard vector of specified size
std::vector<Real> random_vec(size_t N) {
  std::vector<Real> a(N);
  std::generate(a.begin(), a.end(), drand48);
  return a;
}

// Test function for specified number of iterations at chargeand positions
// and positions
template <typename F>
void test(const F& f, const std::string& desc,
          size_t ntest, std::vector<Real>& charge,
          std::vector<Real>& Xs, std::vector<Real>& Ys, std::vector<Real>& Zs,
          std::vector<Real>& Xt, std::vector<Real>& Yt, std::vector<Real>& Zt) {
  std::vector<Real> pot(Zt.size(), 0);
  size_t ns = Xs.size();
  size_t nt = Xt.size();
  double t_start = omp_get_wtime();
  for(int i = 0; i < ntest; ++i)
    f(Xs.data(), Ys.data(), Zs.data(), ns, charge.data(),
      Xt.data(), Yt.data(), Zt.data(), nt, pot.data());
  double t_elapsed = omp_get_wtime() - t_start;
  Real ans = std::accumulate(pot.begin(), pot.end(), 0.0);
  double Gpair_per_sec = ntest * ns * nt / t_elapsed / 1e9;
  std::cout << std::setprecision(16)
            << "N=" << ns << ", M=" << nt << ". "
            << desc << ", ans: " << ans
            << std::endl;
  std::cout << std::setprecision(3)
            << "\t\t\t\ttime: " << t_elapsed << " s      \t"
            << Gpair_per_sec << " Gpair/sec"
            << std::endl;
}

// ------------------------------------------------------------------------------------------
int main (void)
{
  int ns = 10000, nt = 10000, ntest = 20;
  std::vector<Real> Xs(random_vec(ns)), Ys(random_vec(ns)), Zs(random_vec(ns)),
      Xt(random_vec(nt)), Yt(random_vec(nt)), Zt(random_vec(nt)),
      charge(random_vec(ns));

  std::cout << ntest << " repetitions each run..."
            << std::endl;

  // avx512 tests, with and without fast sqrt
  test(&test_kernel_vec512, "manual VCL SIMD avx512",
       ntest, charge, Xs, Ys, Zs, Xt, Yt, Zt);
  test(&test_kernel_vec512_fast_rsqrt, "manual VCL SIMD avx512 fast rsqrt",
       ntest, charge, Xs, Ys, Zs, Xt, Yt, Zt);

  // avx2 tests, with and without fast sqrt
  test(&test_kernel_vec256, "manual VCL SIMD avx2 (256)",
       ntest, charge, Xs, Ys, Zs, Xt, Yt, Zt);
  test(&test_kernel_vec256_fast_rsqrt, "manual VCL SIMD avx2 (256) fast rsqrt",
       ntest, charge, Xs, Ys, Zs, Xt, Yt, Zt);

  // non-avx tests, with differing loop orders
  test(&test_kernel_ts, "target outer loop",
       ntest, charge, Xs, Ys, Zs, Xt, Yt, Zt);
  test(&test_kernel_st, "source outer loop",
       ntest, charge, Xs, Ys, Zs, Xt, Yt, Zt);

  return 0;
}
