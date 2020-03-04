
// demo of manual SIMD via vectorclass library, edited slightly from Libin Lu's
// ccm_simd_demo github repo of 12/10/19.
// Coords are in 3D, and kernel is 1/r (not the exp of abs(r) in 1D that Libin used).

// Note VCL includes rsqrt in single prec (eg 8f, 16f), but not double that we test here.

// TO DO: single prec tests, including rsqrt

// Alex Barnett 3/4/20

// Compile & run:
// g++-9 -fPIC -g -O3 -march=native -funroll-loops -fopenmp -std=c++17 -DVCL -I./version1 -fopenmp main.cc; ./a.out

// Libin's compilation commands:
// icpc -fPIC -g -O3 -march=native -funroll-loops -fopenmp -std=c++17 -DVCL -I./version1 main.cc; ./a.out
// g++ -fPIC -g -O3 -march=native -funroll-loops -fopenmp -std=c++17 -DVCL -I./version1 main.cc; ./a.out
// /cm/shared/sw/pkg/vendor/intel-pstudio/2019-3/vtune_amplifier/bin64/amplxe-cl -collect hotspots ./a.out

#include <math.h> 
#include <vector>
#include <iostream>
#include <iomanip>
#include <omp.h>

// choose double prec throughout
typedef double Real;

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
#endif   // VCL


// ------------------------------------------------------------------------------------------
int main (void)
{
  int ns = 10000, nt = 10000, ntest = 5;
  std::vector<Real> Xs(ns), Ys(ns), Zs(ns), Xt(nt), Yt(nt), Zt(nt), charge(ns), pot(nt);
  std::cout << ntest << " repetitions each run...\n";
  
    // init
  for (auto& x : Xs) x = drand48();   // sources
  for (auto& x : Ys) x = drand48();
  for (auto& x : Zs) x = drand48();
  for (auto& x : Xt) x = drand48();  // targets
  for (auto& x : Yt) x = drand48();
  for (auto& x : Zt) x = drand48();
  for (auto& x : charge) x = drand48();
  
  // kernel run vec 8d  ... only useful if have AVX512
  for (auto& x : pot) x = 0.0;            // needed since funcs add to pot
  Real ts = omp_get_wtime();
  for(int i=0; i<ntest; i++) {
    test_kernel_vec512(&Xs[0], &Ys[0], &Zs[0], ns, &charge[0], 
                   &Xt[0], &Yt[0], &Zt[0], nt, &pot[0]);
  }
  Real t = omp_get_wtime()-ts;
  Real ans = 0;
  for(auto& x : pot) ans+=x;    // this answer checks it agrees btw methods, to 15 digits
  std::cout << std::setprecision(15);
  std::cout<<"N="<<ns<<", M="<<nt<<". manual VCL SIMD avx512, ans: "<<ans<<"\n";
  std::cout << std::setprecision(3);
  std::cout<<"\t\t\t\ttime: "<<t<<" s      \t"<<ntest*ns*nt/t/1e9<<" Gpair/sec\n";

  // kernel run vec 4d   ... tests avx2
  for (auto& x : pot) x = 0.0;
  ts = omp_get_wtime();
  for(int i=0; i<ntest; i++) {
    test_kernel_vec256(&Xs[0], &Ys[0], &Zs[0], ns, &charge[0], 
                   &Xt[0], &Yt[0], &Zt[0], nt, &pot[0]);
  }
  t = omp_get_wtime()-ts;
  ans = 0;
  for(auto& x : pot) ans+=x;
  std::cout << std::setprecision(15);
  std::cout<<"N="<<ns<<", M="<<nt<<". manual VCL SIMD avx2 (256), ans: "<<ans<<"\n";
  std::cout << std::setprecision(3);
  std::cout<<"\t\t\t\ttime: "<<t<<" s      \t"<<ntest*ns*nt/t/1e9<<" Gpair/sec\n";
  
  // kernel run ts
  for (auto& x : pot) x = 0.0;
  ts = omp_get_wtime();
  for(int i=0; i<ntest; i++) {
    test_kernel_ts(&Xs[0], &Ys[0], &Zs[0], ns, &charge[0], 
                   &Xt[0], &Yt[0], &Zt[0], nt, &pot[0]);
  }
  t = omp_get_wtime()-ts;
  ans = 0;
  for(auto& x : pot) ans+=x;
    std::cout << std::setprecision(15);
  std::cout<<"N="<<ns<<", M="<<nt<<". target outer loop, ans: "<<ans<<"\n";
    std::cout << std::setprecision(3);
  std::cout<<"\t\t\t\ttime: "<<t<<" s      \t"<<ntest*ns*nt/t/1e9<<" Gpair/sec\n";
    
  // kernel run st
  for (auto& x : pot) x = 0.0;
  ts = omp_get_wtime();
  for(int i=0; i<ntest; i++) {
    test_kernel_st(&Xs[0], &Ys[0], &Zs[0], ns, &charge[0], 
                   &Xt[0], &Yt[0], &Zt[0], nt, &pot[0]);
  }
  t = omp_get_wtime()-ts;
  ans = 0;
  for(auto& x : pot) ans+=x;
    std::cout << std::setprecision(15);
  std::cout<<"N="<<ns<<", M="<<nt<<". source outer loop, ans: "<<ans<<"\n";
  std::cout << std::setprecision(3);
  std::cout<<"\t\t\t\ttime: "<<t<<" s      \t"<<ntest*ns*nt/t/1e9<<" Gpair/sec\n";
  
  return 0;
}

