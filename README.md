      _   _                                   ______
     | | | |                                 |___   |
     | |_| | ___ _ __ _ __ ___   ___  ___        / /
     |  _  |/ _ \ '__| '_ ` _ \ / _ \/ __|      / /
     | | | |  __/ |  | | | | | |  __/\__ \     / /___
     \_| |_/\___|_|  |_| |_| |_|\___||___/    |______| 


Hermes plasma edge simulation model. Uses BOUT++ framework, adds finite volume
operators and neutral gas models.

This is Hermes-2, a hot ion drift-reduced model.

Author: Ben Dudson, University of York <benjamin.dudson@york.ac.uk>

Released under the GPL license

## License

Full text of the license is in the file LICENSE. If you are using Hermes-2,
please cite the relevant papers.

    Copyright B.Dudson, J.Leddy, University of York, September 2017-2019
              email: benjamin.dudson@york.ac.uk

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

## Installing

This version works with the BOUT++ v4.4, currently the master branch.
Either CMake or Autotools can be used to build it.

### CMake

This is probably the most straightforward method to use now.  First
configure BOUT++ and Hermes-2. To use the default options and minimal
dependencies just run:

    $ cmake . -B build

Alternatively the CMake build can be customised: See the [BOUT++
documentation](https://bout-dev.readthedocs.io/en/latest/user_docs/installing.html#cmake)
for examples of using `cmake` arguments, or edit the compile options
interactively before building:

    $ ccmake . -B build

During configuration
[BOUT++](https://github.com/boutproject/BOUT-dev/) will be
automatically downloaded as a submodule, together with some
dependencies (NetCDF and FFTW are assumed to be installed already,
along with optional dependencies like SUNDIALS and PETSc if they are
requested).  Once configured, run build to compile BOUT++ and then
Hermes-3:

    $ cmake --build build

### Autotools

To build and run tests with GNU autoconf and make, first install BOUT++:

    git clone https://github.com/boutproject/BOUT-dev.git
    cd BOUT-dev

To run this model, preconditioning is strongly recommended, and requires the CVODE solver, part of [SUNDIALS](http://computation.llnl.gov/projects/sundials).
Tested with version 2.6.0. To enable CVODE, BOUT++ should be configured using

    ./configure --with-cvode

or

    ./configure --with-sundials

(which then also enables the IDA solver). Compile BOUT++ with

    make

Then clone the Hermes-2 repository

    git clone https://github.com/bendudson/hermes-2

    cd hermes-2


To compile, run "make" and specify the location of BOUT++

    make BOUT_TOP=/path/to/BOUT-next

This path should be the full path, not relative path, to avoid
problems with compilation in subdirectories.
