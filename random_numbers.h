#ifndef RANDOM_NUMBERS_H
#define RANDOM_NUMBERS_H

#include <time.h>
#include <stdlib.h>

int IntRandom(int target)
{
    return int(rand()/double(RAND_MAX+1)*double(target));
}
double Random(double min,double max)
{
    return rand()/double(RAND_MAX+1)*(max-min)+min;
}

#endif // RANDOM_NUMBERS_H
