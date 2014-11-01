#include "sample.h"

sample::sample(int NewSize, int NewNVars, int NewNClasses)
{
  Size = NewSize;
  NClasses = NewNClasses;
  NVars = NewNVars;
  ProblemType = 0;  //классификация

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

}
sample::sample(int NewSize, int NewNCols, int NewNVars, int NewNOuts)
{
  Size = NewSize;
  NCols = NewNCols;
  NVars = NewNVars;
  NOuts = NewNOuts;
  ProblemType = 1;  //регрессия

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
    delete Classes;
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
