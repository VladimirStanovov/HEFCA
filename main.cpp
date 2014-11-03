#include <iostream>
#include <random_numbers.h>
#include <sample.cpp>

using namespace std;

int main()
{
  cout << "Hello World!" << endl;
  cout << "I'm trying to use GitHub"<<endl;

  sample S(150,4,3,10,0.5);
  S.ReadFileClassification((char*)"iris_.txt");
  //S.ShowSampleClassification();
  S.ClassPatternsCalc();
  S.SplitCVStratified();
  int FoldOnTest = 0;
  sample S_CVLearn(S.GetCVLearnSize(FoldOnTest),S.NVars,S.NClasses,10,0.9);
  sample S_CVTest(S.GetCVTestSize(FoldOnTest),S.NVars,S.NClasses,10,0.9);

  S.SetCVLearn(S_CVLearn,FoldOnTest);
  S.SetCVTest(S_CVTest,FoldOnTest);

  S_CVLearn.ShowSampleClassification();
  cout<<endl;
  S_CVTest.ShowSampleClassification();

  cout<<endl<<endl;
  S.SplitStratified();
  sample S_Learn(S.GetLearnSize(),S.NVars,S.NClasses,10,0.9);
  sample S_Test(S.GetTestSize(),S.NVars,S.NClasses,10,0.9);

  S.SetLearn(S_Learn);
  S.SetTest(S_Test);

  S_Learn.ShowSampleClassification();
  cout<<endl;
  S_Test.ShowSampleClassification();

  return 0;
}

