#ifndef SAMPLE_H
#define SAMPLE_H

#include <fstream>
#include <string.h>

using namespace std;

class sample
{
public:
  //конструктор для задач регрессии
  sample(int NewSize, int NewNCols, int NewNVars, int NewNOuts);
  //конструктор для задач классификации
  sample(int NewSize, int NewNVars, int NewNClasses);
  ~sample();
  // первичное считывание с файла
  void ReadFileClassification(char* filename);
  void ReadFileRegression(char* filename);
  //вывод на экран всей выборки
  void ShowSampleClassification();
  void ShowSampleRegression();
  // общие параметры
  int Size;         //объем выборки
  int NCols;        //общее число столбцов в выборке
  int NVars;        //число столбцов входных параметров
  int NOuts;        //число столбцов выходных параметров
  int ProblemType;  //тип задачи

  double** Inputs;  //входы задачи
  double** Outputs; //выходы задачи
  bool** MissingInputs;   //массив пропущенных входных значений
  bool** MissingOutputs;   //массив пропущенных входных значений

  // параметры для задач классификации
  int NClasses;     //число классов в задаче
  int* Classes;     //массив номеров классов

};

#endif // SAMPLE_H
