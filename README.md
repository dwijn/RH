# RH radiative transfer code

See the installation and other documentation in the
[PDF manual](docs/rhmanual.pdf). A [changelog](docs/changelog) is also
available.

The canonical reference for RH is Uitenbroek, H. 2001, ApJ,
557, 389 (ADS: [2001ApJ...557..389U](https://ui.adsabs.harvard.edu/abs/2001ApJ...557..389U/abstract),
DOI: [10.1086/321659](http://dx.doi.org/10.1086/321659)).

## Building

RH needs a C11 compiler, CMake 3.20 or newer, POSIX threads and an XDR
library. On Linux, XDR comes from libtirpc: install `libtirpc-devel`
(RHEL, Fedora, SUSE) or `libtirpc-dev` (Debian, Ubuntu). On macOS it is part
of the system.

```
cmake -B build
cmake --build build -j
ctest --test-dir build                    # quick check that the programs start
cmake --install build --prefix ~/rh       # optional: ~/rh/bin/rhf1d, ...
```

The programs are built in `build/rhf1d/rhf1d`, `build/rhsc2d/rhsc2d`,
`build/rhsc3d/rhsc3d` and `build/rhsphere/rhsphere`, next to the auxiliary
programs of each geometry (`solveray`, `backgrcontr`, ...); the atomic-data
tools are in `build/tools/`. Run a program from a run directory, e.g.
`cd rhf1d/run && ../../build/rhf1d/rhf1d`. Installed, the auxiliary
programs that exist for several geometries get a suffix (`solveray_1d`,
`solveray_2d`, ...).

The default build type is `RelWithDebInfo` (`-O2 -g`). Options, given as
`-D<option>=<value>` when configuring:

| Option | Effect |
|---|---|
| `CMAKE_BUILD_TYPE=Release` | `-O3`, no debug information |
| `CMAKE_C_COMPILER=icx` | another compiler (on a Cray system: `cc`) |
| `RH_NATIVE=ON` | tune for the build machine's CPU (`-march=native`); the programs may not run on other machines |
| `RH_SANITIZE=address,undefined` | AddressSanitizer and UndefinedBehaviorSanitizer, for debugging |
| `RH_WARNINGS=ON` | compile with `-Wall -Wextra` |
| `RH_BUILD_TOOLS=OFF` | skip the tools in `tools/` |

Intel compilers are run with `-fp-model=precise`: RH needs IEEE
floating-point semantics.
