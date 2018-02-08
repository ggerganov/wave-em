#!/bin/bash

echo "static const char * BUILD_TIMESTAMP=\"`date`\";" > build_timestamp.h

em++ -O3 -std=c++11 -s USE_SDL=2 ./main.cpp -o wave.js -I ./fftw-3.3.3/api ./lib/libfftw3f.a \
    -s EXPORTED_FUNCTIONS='["_getText", "_getSampleRate", "_setText", "_getAverageRxTime_ms", "_setParameters", "_main"]' \
    -s EXTRA_EXPORTED_RUNTIME_METHODS='["ccall", "cwrap"]'
