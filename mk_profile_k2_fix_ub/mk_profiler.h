#ifndef MK_PROFILER_H
#define MK_PROFILER_H
#include <chrono>
#include <iostream>
struct MKProfiler {
    double time_input_packing; double time_blind_rotate_control;
    double time_external_product; double time_sample_extract;
    double time_keygen; double time_encrypt;
    MKProfiler() { time_input_packing=0; time_blind_rotate_control=0; time_external_product=0; time_sample_extract=0; time_keygen=0; time_encrypt=0; }
};
extern MKProfiler global_profiler;
#endif
