#include <iostream>
#include <random_numbers.h>
#include <sample.cpp>

using namespace std;

int main()
{
  cout << "Hello World!" << endl;
  cout << "I'm trying to use GitHub"<<endl;

  sample S(150,4,3,10,0.3,1);
  S.ReadFileClassification((char*)"iris_.txt");
  S.ShowSampleClassification();
  return 0;
}

