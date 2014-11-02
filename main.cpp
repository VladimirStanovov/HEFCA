#include <iostream>
#include <random_numbers.h>
#include <sample.cpp>

using namespace std;

int main()
{
  cout << "Hello World!" << endl;
  cout << "I'm trying to use GitHub"<<endl;

  sample S(150,4,3,10,0.3);
  S.ReadFileClassification((char*)"iris_.txt");
  S.ShowSampleClassification();
  S.ClassPatternsCalc();
  S.SplitCVStratified();
  S.SplitCVRandom();
  S.SplitStratified();
  S.SplitRandom();
  return 0;
}

