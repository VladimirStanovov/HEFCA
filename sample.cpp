#include "sample.h"

sample::sample(int NewSize, int NewNVars, int NewNClasses, int NewNFolds,
               double NewSplitRate)
{
  Size = NewSize;
  NClasses = NewNClasses;
  NVars = NewNVars;
  ProblemType = 0;  //классификация
  NFolds = NewNFolds;
  SplitRate = NewSplitRate;

  ClassPerFold = new int*[NClasses];
  ClassPositions = new int*[NClasses];
  CVFoldNum = new int[Size];
  FoldSize = new int[NFolds];
  CVSplit = new int[Size];
  NClassInst = new int[NClasses];
  Classes = new int[Size];
  Inputs = new double*[Size];
  MissingInputs = new bool*[Size];
  for(int i=0;i!=Size;i++)
  {
    Inputs[i] = new double[NVars];
    MissingInputs[i] = new bool[NVars];
    for(int j=0;j!=NVars;j++)
    {
        MissingInputs[i][j] = false;
    }
  }
  for(int i=0;i!=NClasses;i++)
  {
      ClassPositions[i] = new int[Size];
      ClassPerFold[i] = new int[NFolds];
  }
}
sample::sample(int NewSize, int NewNCols, int NewNVars, int NewNOuts,
               int NewNFolds, double NewSplitRate)
{
  Size = NewSize;
  NCols = NewNCols;
  NVars = NewNVars;
  NOuts = NewNOuts;
  ProblemType = 1;  //регрессия
  NFolds = NewNFolds;
  SplitRate = NewSplitRate;

  CVFoldNum = new int[Size];
  FoldSize = new int[NFolds];
  CVSplit = new int[Size];
  Inputs = new double*[Size];
  Outputs = new double*[Size];
  MissingInputs = new bool*[Size];
  MissingOutputs = new bool*[Size];
  for(int i=0;i!=Size;i++)
  {
    Inputs[i] = new double[NVars];
    Outputs[i] = new double[NOuts];
    MissingInputs[i] = new bool[NVars];
    MissingOutputs[i] = new bool[NOuts];
    for(int j=0;j!=NVars;j++)
    {
        MissingInputs[i][j] = false;
    }
    for(int j=0;j!=NOuts;j++)
    {
        MissingOutputs[i][j] = false;
    }
  }
}
sample::~sample()
{
  delete FoldSize;
  delete CVSplit;
  for(int i=0;i!=Size;i++)
  {
    delete Inputs[i];
    delete MissingInputs[i];
  }
  if(ProblemType == 1)
    for(int i=0;i!=Size;i++)
    {
      delete Outputs[i];
      delete MissingOutputs[i];
    }
  if(ProblemType == 0)
  {
    delete Classes;
    delete NClassInst;
    for(int i=0;i!=NClasses;i++)
    {
      delete ClassPositions[i];
      delete ClassPerFold[i];
    }
    delete ClassPositions;
    delete ClassPerFold;
  }
}
void sample::ReadFileClassification(char* filename)
{
   std::ifstream fin(filename);
   char tempstring [80];
   for(int i=0;i!=Size;i++)
   {
       for(int j=0;j!=NVars;j++)
       {
          fin>>tempstring;
          if(strcmp(tempstring,"?") == 0)
          {
              MissingInputs[i][j] = true;
              Inputs[i][j] = 0;
          }
          else
          {
              Inputs[i][j] = atof(tempstring);
          }
       }
       fin>>Classes[i];
   }
}
void sample::ReadFileRegression(char* filename)
{
   std::ifstream fin(filename);
   char tempstring [80];
   for(int i=0;i!=Size;i++)
   {
       for(int j=0;j!=NVars;j++)
       {
          fin>>tempstring;
          if(strcmp(tempstring,"?") == 0)
          {
              MissingInputs[i][j] = true;
              Inputs[i][j] = 0;
          }
          else
          {
              Inputs[i][j] = atof(tempstring);
          }
       }
       for(int j=0;j!=NOuts;j++)
       {
          fin>>tempstring;
          if(strcmp(tempstring,"?") == 0)
          {
              MissingOutputs[i][j] = true;
              Outputs[i][j] = 0;
          }
          else
          {
              Outputs[i][j] = atof(tempstring);
          }
       }
   }
}
void sample::ShowSampleRegression()
{
  for(int i=0;i!=Size;i++)
  {
      for(int j=0;j!=NVars;j++)
      {
          if(MissingInputs[i][j])
            cout<<"?\t";
          else
            cout<<Inputs[i][j]<<"\t";
      }
      cout<<"->\t";
      for(int j=0;j!=NOuts;j++)
      {
          if(MissingOutputs[i][j])
            cout<<"?\t";
          else
            cout<<Outputs[i][j]<<"\t";
      }
      cout<<endl;
  }
}
void sample::ShowSampleClassification()
{
  for(int i=0;i!=Size;i++)
  {
      for(int j=0;j!=NVars;j++)
      {
          if(MissingInputs[i][j])
            cout<<"?\t";
          else
            cout<<Inputs[i][j]<<"\t";
      }
      cout<<"->\t";
      cout<<Classes[i]<<endl;
  }
}
double sample::GetValue(int Num,int Var)
{
  return Inputs[Num][Var];
}
double sample::GetOutput(int Num,int Var)
{
  return Outputs[Num][Var];
}
int sample::GetClass(int Num)
{
  return Classes[Num];
}
void sample::SplitCVStratified()
{
  int counter=0;
  int RandomPattern;
  for(int i=0;i!=Size;i++)
  {
    CVFoldNum[i] = -1;
  }
  for(int i=0;i!=NFolds;i++)
  {
      for(int j=0;j!=NClasses;j++)
      {
          for(int k=0;k!=ClassPerFold[j][i];k++)
          {
              while(counter < (int) ((double)i+1.0)*(double)Size/(double)NFolds )
              {
                RandomPattern = IntRandom(NClassInst[j]);
                while(CVFoldNum[ClassPositions[j][RandomPattern]] != -1)
                {
                  RandomPattern++;
                  if(RandomPattern == NClassInst[j])
                    RandomPattern=0;
                }

                CVFoldNum[ClassPositions[j][RandomPattern]] = i;
                counter++;
              }
          }
      }
  }
}
void sample::SplitCVRandom()
{
  int counter=0;
  int RandomPattern;
  for(int i=0;i!=Size;i++)
  {
    CVFoldNum[i] = -1;
  }
  for(int i=0;i!=NFolds;i++)
  {
    while(counter < (int) ((double)i+1.0)*(double)Size/(double)NFolds )
    {
      RandomPattern = IntRandom(Size);
      while(CVFoldNum[RandomPattern] != -1) //Если измерение уже взято в
      {                                     //одну из выборок, то увеличиваем
        RandomPattern++;                    //номер, пока не найдем не взятое
        if(RandomPattern == Size)           //таким образом избегаем смещения
          RandomPattern=0;                  //если достигли конца, начинаем
      }                                     //сначала.

      CVFoldNum[RandomPattern] = i;
      counter++;
    }
  }
}
void sample::ClassPatternsCalc()
{
  for(int i=0;i!=NClasses;i++)
  {
      NClassInst[i] = 0;
  }
  for(int i=0;i!=Size;i++)
  {
      NClassInst[GetClass(i)]++;
  }
  for(int i=0;i!=NClasses;i++)
  {
      int counter = 0;
      for(int j=0;j!=Size;j++)
      {
          if(GetClass(j) == i)
          {
              ClassPositions[i][counter] = j;
              counter++;
          }
      }
      counter = 0;
      for(int j=0;j!=NFolds;j++)
      {
          ClassPerFold[i][j] = 0;
          while(counter < int( ((double)j+1.0)*(double)NClassInst[i]/(double)NFolds ) )
          {
              ClassPerFold[i][j] ++;
              counter++;
          }
      }
  }
}
void sample::SplitRandom()
{
  int RandomPattern;
  for(int i=0;i!=Size;i++)
  {
    CVFoldNum[i] = -1;  //используется здесь как индикатор
                        //использования измерения
  }
  LearnSize = int(SplitRate*(double)Size);
  TestSize = Size-LearnSize;
  for(int i=0;i!=LearnSize;i++)
  {
    RandomPattern = IntRandom(LearnSize);
    while(CVFoldNum[RandomPattern] != -1)
    {
      RandomPattern++;
      if(RandomPattern == LearnSize)
        RandomPattern=0;
    }
    CVFoldNum[RandomPattern] = 0;
  }
  for(int i=0;i!=TestSize;i++)
  {
    RandomPattern = IntRandom(TestSize);
    while(CVFoldNum[RandomPattern] != -1)
    {
      RandomPattern++;
      if(RandomPattern == TestSize)
        RandomPattern=0;
    }
    CVFoldNum[RandomPattern] = 1;
  }
}
void sample::SplitStratified()
{
  int counter=0;
  int RandomPattern;
  for(int i=0;i!=Size;i++)
  {
    CVFoldNum[i] = -1;  //используется здесь как индикатор
                        //использования измерения
  }
  for(int i=0;i!=NClasses;i++)
  {
    for(int j=0;j!=round(SplitRate*NClassInst[i]);j++)
    {
      RandomPattern = IntRandom(NClassInst[j]);
      while(CVFoldNum[ClassPositions[j][RandomPattern]] != -1)
      {
        RandomPattern++;
        if(RandomPattern == NClassInst[j])
          RandomPattern=0;
      }
      CVFoldNum[ClassPositions[j][RandomPattern]] = 0;
      counter++;
    }
  }
  LearnSize = counter;
  TestSize = Size-LearnSize;
  for(int i=0;i!=Size;i++)
  {
    if(CVFoldNum[i] == -1)
      CVFoldNum[i] = 1;
  }
}
