#include "time.h"
double
get_time ()
{
    struct timeval t;
    gettimeofday (&t, 0);
    return (double)t.tv_sec + (double)t.tv_usec * 1e-6;
}
double
get_cpu_time ()
{
    struct rusage t;
    getrusage (RUSAGE_THREAD, &t);
    return (double)t.ru_utime.tv_sec + (double)t.ru_utime.tv_usec * 1e-6;
}
