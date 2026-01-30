# -Ofast is deprecated, use: -O3 -ffast-math
# -fsanitize=address -g3 \
# -lprofiler \
clang++ -Wall -Wextra -Wpedantic -std=gnu++23 \
  -O3 -ffast-math -march=armv8.6-a -mtune=apple-m3 \
  -I/Users/jw/Projects/mcl/include \
  -I/opt/homebrew/opt/libomp/include \
  vfhe.cpp /Users/jw/Projects/mcl/build/lib/libmcl.a \
  -DTIMING_ON -DN_THREADS=4 -DASSERT_ON \
  -DMCL_USE_OMP=1 \
  -Xpreprocessor -fopenmp \
  -L/opt/homebrew/opt/libomp/lib -lomp \
  -pthread -lntl -lgmpxx -lgmp \
  -Wl,-rpath,/opt/homebrew/opt/libomp/lib \
  -o vfhe
