#### Mac
First build MCL
```
mkdir build
cd build
cmake ..
make
```

Then build VFHE (optimised)
```
build_opt_mac.sh
```

Then build VFHE (debug)
```
build_debug_mac.sh
```

```
g++ -Wall -Wextra -Wpedantic -std=gnu++23 \
-Ofast -march=armv8.6-a -mtune=apple-m3 \
-I/Users/jw/Projects/mcl/include \
-L/Users/jw/Projects/mcl/build/lib \
vfhe.cpp \
/Users/jw/Projects/mcl/build/lib/libmcl.a \
-DMCL_USE_OMP=1 -fopenmp \
-pthread -lntl \
-lgmpxx -lgmp \
-o vfhe

g++ -Wall -Wextra -Wpedantic -std=gnu++23 \
-Og -g3 -march=armv8.6-a -mtune=apple-m3 \
-I/Users/jw/Projects/mcl/include \
-L/Users/jw/Projects/mcl/build/lib \
vfhe.cpp \
/Users/jw/Projects/mcl/build/lib/libmcl.a \
-DMCL_USE_OMP=1 -fopenmp \
-pthread -lntl \
-lgmpxx -lgmp \
-o vfhe
```