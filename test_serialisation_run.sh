    -L/home/jw/Projects/mcl/build/lib \

g++ -Wall -Wextra -Wpedantic -std=gnu++23 \
    -march=znver3 -mtune=znver3 -Ofast \
    -flto -fomit-frame-pointer -fstrict-aliasing \
    -I/home/jw/Projects/mcl/include \
    -DMCL_USE_OMP=1 -fopenmp -lmcl \
    -pthread -lntl \
    -lgmpxx -lgmp \
    -o test_ntl_mcl test_ntl_mcl.cpp \
&& LD_LIBRARY_PATH=/home/jw/Projects/mcl/build/lib:$LD_LIBRARY_PATH ./test_ntl_mcl