# WASM interface code to diagonalize matrix using LAPACK(E).
We use WebAssembly to speed up the diagonalization routine central to the
normal modes calculation. There are several libraries available to this end
in pure JavaScript but we were unsatisfied with their performance.

So, and because we had time to spare, we decided to try and link the good
old LAPACK routines that many in the scientific community use. To this end,
we wrote a small C code that uses the LAPACK C interface (LAPACKE) to run
our diagonalization and made it work, through WebAssembly, with the rest of
our Javascript code.

What follows is a very simple and probably incomplete guide to help others
interested in getting their FORTRAN/C code running in the browser. We did
manage to call the actual FORTRAN routines from C and transpile them to a
WASM module, but in the end we were not masochistic enough to run with that
as our final choice. Happy to share the code though, just email me or write
me on Twitter.

# How to compile the C/FORTRAN code to WASM/JS.
The end goal here is to get a WASM module (.wasm) and the JavaScript glue
code produced by Emscripten, starting from a C program. Our C program does
not include a `main()` function on purpose, as we ran into some issues with
WASM running it automatically. Probably there is a nice way to prevent this
from happening but we couldn't figure it out (or if we did, we don't remember
it anymore). Also, less code, simpler to read/maintain.

## Setup the compiler/transpiler toolchain
The very first step is to setup the toolchain to compile C/FORTRAN to WASM.
We followed the instructions in [this](https://github.com/StarGate01/Full-Stack-Fortran)
amazing repository and [this](https://chrz.de/2020/04/21/fortran-in-the-browser/)
blog post, with a few changes. In short, and just in case those links rot,
here are a short version of those instructions for an Ubuntu machine 
(correct/adapt as necessary):

```bash
# Make a 'base' folder to download install all these tools
mkdir app
cd app

# Setup compiler toolchain with old compatible compilers
echo "deb http://archive.ubuntu.com/ubuntu/ trusty main restricted universe" >> /etc/apt/sources.list
echo "deb http://security.ubuntu.com/ubuntu/ trusty-security main" >> /etc/apt/sources.list
apt-get update
apt-get -y install --no-install-recommends \
    build-essential python wget git gnutls-bin bash make ca-certificates xz-utils \
    gfortran-4.6 g++-4.6 gcc-4.6 gcc-4.6-plugin-dev llvm-3.3-dev
rm -rf /var/lib/apt/lists/*

# dragonegg
git clone --single-branch --branch="release_33" --depth=1 https://github.com/llvm-mirror/dragonegg.git && \
cd dragonegg
LLVM_CONFIG=llvm-config-3.3 GCC=gcc-4.6 CC=gcc-4.6 CXX=g++-4.6 make -j$(nproc)
mkdir ../bin
mv dragonegg.so ../bin/

# emscripten
git clone --depth=1 https://github.com/emscripten-core/emsdk.git
cd emsdk
./emsdk install "1.39.12"  # We actually upgraded to 1.39.20
./emsdk activate "1.39.12"
mkdir ../scripts/
cp -r scripts ../scripts/
# Adjust these paths as necessary
export PATH="app/bin:app/scripts:app/emsdk:app/emsdk/node/12.9.1_64bit/bin:app/emsdk/upstream/emscripten:${PATH}"

# libgfortran
git clone --single-branch --branch="releases/gcc-4.6" --depth=1 https://github.com/gcc-mirror/gcc.git
cp -r vendor/gfortran gcc-build  # vendor is the folder we provide here, copied from those original sources
cd gcc-build
./patch.sh && make -j$(nproc) build

# lapack
git clone --depth=1 https://github.com/Reference-LAPACK/lapack.git
cp vendor/LAPACK/make.inc lapack/make.inc
cd lapack/SRC
emmake make -j$(nproc) single  # Make sure to add double if necessary for your code.
```

In addition to these instructions, we had to add an include to stdlib.h to both GCC
getlog.c and ctime.c to actually get the code to compile.

After all this was done, we copied `libgfortran.a`, `liblapack.a`, `liblapacke.a`,
and `librefblas.a` to a `f2wasm_lib` folder, which you will see referenced below.

## Compile the C/FORTRAN code to WASM
Building a WASM module is not as trivial as the documentation makes it to be. We ran
into a _lot_ of trouble to get the code to execute properly. We found that the easiest
way was to actually follow the 'hard way'. We also ran into an issue with using 64-bit
code on the browser, because apparently browsers don't like that. That's why you will
see all those `m32` in the following commands.

To compile our code on your machine, run these three lines sequentially. In short, they
first translate the C code to bytecode, then link all the libraries, and finally make
the WASM code (`diag.wasm`) and the JavaScript glue code (`diag.js`). We ran all these
commands on a Windows 10 machine with WSL2 and Ubuntu 18.04, so I will assume they work
equally well on a regular Linux machine.

On the last step, we export the diagonalize function only because that's what we will call.
We don't use the `ccall` and `cwrap` functions that are often mentioned in the tutorials
because we got into a lot of trouble getting them to behave properly. Probably ignorance
from us but hey, at least this works! Ah, also, running with `-O3` was discouraged, but
we cannot remember if it was LAPACK or EMCC, so we left it as `-O2`.

```bash
emcc -O2 --target=wasm32-unknown-emscripten -c -flto -emit-llvm -m32 -I/home/joaor/lapack/LAPACKE/include -Wall -o xyzdiag.bc xyzdiag.c

emcc -O2 -m32 -Wall -flto -r -o assembly.bc xyzdiag.bc ../f2wasm_lib/liblapacke.a ../f2wasm_lib/liblapack.a ../f2wasm_lib/librefblas.a ../f2wasm_lib/libgfortran.a

emcc -O2 -m32 -Wall -flto -s ERROR_ON_UNDEFINED_SYMBOLS=1 -s EXPORTED_FUNCTIONS='["_diagonalize"]' -o diag.js assembly.bc -s ALLOW_MEMORY_GROWTH=1
```

## Additional Notes
Apparently all this will be made much simpler when LLVM supports gfortran out of the box,
which apparently is in the works. Then, it might be truly much simpler than getting
all this weird toolchain setup. Keep an eye open for that! It's called FLANG.

Good luck and let us know if you have questions, suggestions, etc.
