g++ -Wall -Wextra -Wpedantic -std=gnu++23 \
    -march=znver3 -mtune=znver3 -Ofast \
    -flto -fomit-frame-pointer -fstrict-aliasing \
    -I/home/jw/Projects/mcl/include \
    vfhe.cpp /home/jw/Projects/mcl/build/lib/libmcl.a \
    -DTIMING_ON -DN_THREADS=16 -DASSERT_ON \
    -DMCL_USE_OMP=1 \
    -fopenmp \
    -pthread -lntl -lgmpxx -lgmp \
    -o vfhe