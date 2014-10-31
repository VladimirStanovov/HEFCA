#ifndef SAMPLE_H
#define SAMPLE_H

#include <fstream>

using namespace std;

class sample
{
public:
  sample();
  ~sample();
  // первичное считывание с файла
  void ReadFileClassification(char* filename, int Size, int NCols,
                              int NVars, int NOuts);
  // общие параметры
  int Size;
  int NCols;
  int NVars;
  int NOuts;
  int ProblemType;

  // параметры для задач классификации
  int NClasses;
};

#endif // SAMPLE_H
