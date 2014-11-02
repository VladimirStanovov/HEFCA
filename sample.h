#ifndef SAMPLE_H
#define SAMPLE_H

#include <fstream>
#include <string.h>

using namespace std;

class sample
{
public:
  //конструктор для задач регрессии
  sample(int NewSize, int NewNCols, int NewNVars, int NewNOuts, int NewNFolds, double NewSplitRate);
  //конструктор для задач классификации
  sample(int NewSize, int NewNVars, int NewNClasses, int NewNFolds, double NewSplitRate);
  ~sample();
  // первичное считывание с файла
  void ReadFileClassification(char* filename);
  void ReadFileRegression(char* filename);
  //вывод на экран всей выборки
  void ShowSampleClassification();
  void ShowSampleRegression();
  // взять значение переменной из выборки
  double GetValue(int Num,int Var);
  // взять значение выхода из выборки
  double GetOutput(int Num,int Var);
  // получить номер класса для измерения
  int GetClass(int Num);
  // общие параметры
  int Size;         //объем выборки
  int NCols;        //общее число столбцов в выборке
  int NVars;        //число столбцов входных параметров
  int NOuts;        //число столбцов выходных параметров
  int ProblemType;  //тип задачи
  int NFolds;       //число частей при кросс-валидации
  double SplitRate; //доля обучающей выборки, например
                    //0.7 => разбиение 70/30

  double** Inputs;  //входы задачи
  double** Outputs; //выходы задачи
  bool** MissingInputs;   //массив пропущенных входных значений
  bool** MissingOutputs;   //массив пропущенных входных значений
  int* CVSplit;     //определяет к какой части относится
                    //измерение при кросс-валидации

  // параметры для задач классификации
  int NClasses;     //число классов в задаче
  int* Classes;     //массив номеров классов
  int* NClassInst;  //количество объектов в классах

};

#endif // SAMPLE_H
