#include "sample.h"

sample::sample()
{
}
sample::~sample()
{
}
void sample::ReadFileClassification(char* filename, int NewSize, int NewNCols,
                                    int NewNVars, int NewNOuts)
{
   std::ifstream fin(filename);
   Size = NewSize;
   NCols = NewNCols;
   ProblemType = 0;
   NVars = NewNVars;
   NOuts = NewNOuts;


}
