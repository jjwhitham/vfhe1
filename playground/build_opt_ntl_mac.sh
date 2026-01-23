# -Ofast is deprecated, use: -O3 -ffast-math
clang++ -Wall -Wextra -Wpedantic -std=gnu++23 \
  -O3 -ffast-math -march=armv8.6-a -mtune=apple-m3 \
  test_ntl.cpp \
  -pthread -lntl \
  -Wl,-rpath,/opt/homebrew/opt/libomp/lib \
  -o test_ntl