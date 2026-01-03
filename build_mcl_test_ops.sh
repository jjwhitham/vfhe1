clang++ -Wall -Wextra -Wpedantic -std=gnu++23 \
  -fno-omit-frame-pointer \
  -Og -g3 -march=armv8.6-a -mtune=apple-m3 \
  -I/Users/jw/Projects/mcl/include \
  -I/opt/homebrew/opt/libomp/include \
  mcl_test_ops.cpp /Users/jw/Projects/mcl/build/lib/libmcl.a \
  -DTIMING_ON -DN_THREADS=4 \
  -DMCL_USE_OMP=1 \
  -Xpreprocessor -fopenmp \
  -L/opt/homebrew/opt/libomp/lib -lomp \
  -lprofiler \
  -pthread -lntl -lgmpxx -lgmp \
  -Wl,-rpath,/opt/homebrew/opt/libomp/lib \
  -o mcl_test_ops
