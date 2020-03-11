# floatingspeed

Various tests of numerical performance, comparing languages, compilers, flags, etc.

Alex Barnett 2016--2020

### contents

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

# Open in Visual Studio Code
code .
```

Install the [Remote-Containers vscode extension](https://marketplace.visualstudio.com/items?itemName=ms-vscode-remote.remote-containers) (if not already installed)

In the lower-left corner of the vscode window you should see a little green icon. Click that and select "Reopen in Container"

Now the container will be built on your system (may take several minutes). And then your project will appear -- inside the development environment container.

If you are curious, the Docker recipe for installing all the compilers can be found at [.devcontainer/Dockerfile](.devcontainer/Dockerfile)

You can now either run the tests from the built-in vscode terminal, or you can run them using the built-in tasks for this project.

To run a task, launch the vscode command menu (Ctrl+Shift+P) and search for "Run Task" and press enter. Then select one of the two tasks that appear (e.g., "Run arraysinglethread test") and press enter. Then enter again to continue. A new vscode terminal will open and run your test.

If you are curious, the tasks are configured at [.vscode/tasks.json](.vscode/tasks.json).



