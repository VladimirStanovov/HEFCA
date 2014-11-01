#include "sample.h"

sample::sample(int NewSize, int NewNCols, int NewNVars, int NewNClasses)
{
  Size = NewSize;
  NCols = NewNCols;
  NClasses = NewNClasses;
  NVars = NewNVars;
  NOuts = 1;

}
sample::sample(int NewSize, int NewNCols, int NewNVars, int NewNOuts)
{
  Size = NewSize;
  NCols = NewNCols;
  NVars = NewNVars;
  NOuts = NewNOuts;

}
sample::~sample()
{
}
void sample::ReadFileClassification(char* filename)
{
   std::ifstream fin(filename);
}
