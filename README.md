# floatingspeed

Various tests of numerical performance, comparing languages, compilers, flags, etc.

Alex Barnett 2016--2026

Contributions by:  
Jeremy Magland  
Libin Lu  
many others  

## Contents

* ``arraysinglethread``: compare C/C++/Fortran speeds for complex and real double-precision arithmetic, on a 1D array.
* ``lap3dkernel``: compare various languages speed in a multithreaded implementation of direct summation of the 1/r kernel.


## Open in a containerized development environment using vscode

You can use [Visual Studio Code](https://code.visualstudio.com/) to open a containerized development environment. This creates a Docker container with all of the compilers and tools pre-installed, including Fortran, Python, Octave, g++, and Julia. You can then run the tests and seamlessly edit the code all within this environment, and the source files are automatically sync'd with the host environment.

Note: MATLAB is not installed within the development container.

### Prerequisites

* [Docker](https://docs.docker.com/install/)
* [Visual Studio Code](https://code.visualstudio.com/)

### Instructions

```bash
# Clone this repo
git clone [repo-url] floatingspeed
cd floatingspeed

### Open in Visual Studio Code
code .
```

Install the [Remote-Containers vscode extension](https://marketplace.visualstudio.com/items?itemName=ms-vscode-remote.remote-containers) (if not already installed)

In the lower-left corner of the vscode window you should see a little green icon. Click that and select "Reopen in Container"

Now the container will be built on your system (may take several minutes). And then your project will appear -- inside the development environment container.

If you are curious, the Docker recipe for installing all the compilers can be found at [.devcontainer/Dockerfile](.devcontainer/Dockerfile)

You can now either run the tests from the built-in vscode terminal, or you can run them using the built-in tasks for this project.

To run a task, launch the vscode command menu (Ctrl+Shift+P) and search for "Run Task" and press enter. Then select one of the two tasks that appear (e.g., "Run arraysinglethread test") and press enter. Then enter again to continue. A new vscode terminal will open and run your test.

If you are curious, the tasks are configured at [.vscode/tasks.json](.vscode/tasks.json).


## Summary of this repo produced by codex CLI using ChatGPT 5.2

 I'm preparing a concise high-level explanation covering the two main
  benchmarks—arraysinglethread testing complex arithmetic in a single thread,
  and lap3dkernel performing 3D 1/r sums across languages with SIMD and OpenMP.
  I'll describe key files, script entry points, output locations, and note VCL
  inclusion, while briefly addressing the non-critical display errors. I'll also
─ Worked for 1m 11s ────────────────────────────────────────────────────────────

• Here’s the shape of the repo and what each part does, based on the files in /
  home/alex/numerics/floatingspeed:

  1) arraysinglethread/ — single‑core array arithmetic

  - Purpose: compare real vs complex arithmetic performance in C, C++, and
    Fortran on large 1D arrays, mostly memory‑bound. See arraysinglethread/
    README.
  - Runner script: arraysinglethread/do_complextimings.sh compiles and runs
    three binaries (Fortran, C++ with C complex, C++ with std::complex). You
    tweak FLAGS, GF, and GXX here.
  - C++ test: arraysinglethread/complexmulttiming.cpp allocates large arrays
    (default M=1e8), fills with fixed random data, and times multiply‑add loops
    across several representations:
      - real doubles
      - “by hand” complex (separate real/imag arrays and interleaved)
      - C double complex (if USE_C_TYPE_COMPLEX)
      - C++ std::complex<double> (if USE_CPP_TYPE_COMPLEX)
  - Fortran test: arraysinglethread/complexmulttiming.f does the same for real*8
    and complex*16.
  - Timing utility: arraysinglethread/utils.h + arraysinglethread/utils.cpp
    provide a simple CNTime timer.
  - Results: arraysinglethread/results/*.txt are saved benchmark outputs for
    various machines/flags.

  2) lap3dkernel/ — multithreaded 3D 1/r kernel

  - Purpose: benchmark dense direct N‑body summation (1/r kernel) from ns=1e4
    sources to nt=1e4 targets. See lap3dkernel/README.
  - Runner script: lap3dkernel/run_all_tests.sh compiles/runs C++ SIMD + OpenMP,
    Fortran + OpenMP, Julia, Python+Numba, MATLAB, Octave. It expects those
    toolchains to be installed.
  What the benchmarks actually measure

  - arraysinglethread: raw array traversal + multiply‑add with accumulation;
    highlights compiler flags and complex number representation overhead.
  - lap3dkernel: “source‑target pairs per second” for 3D 1/r kernel; highlights
    threading, SIMD, and compiler differences across languages.

  If you want, tell me which part you want to go deeper on (e.g., “walk through
  the C++ SIMD kernel,” “what flags matter most,” or “how to run just one
  benchmark”), and I’ll zoom in.


