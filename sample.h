#ifndef SAMPLE_H
#define SAMPLE_H

#include <fstream>
#include <string.h>
#include <math.h>

using namespace std;

class sample
{
public:
  //конструктор для задач регрессии
  sample(int NewSize, int NewNCols, int NewNVars, int NewNOuts,
         int NewNFolds, double NewSplitRate);
  //конструктор для задач классификации
  sample(int NewSize, int NewNVars, int NewNClasses, int NewNFolds,
         double NewSplitRate);
  ~sample();
  //первичное считывание с файла
  void ReadFileClassification(char* filename);
  void ReadFileRegression(char* filename);
  //вывод на экран всей выборки
  void ShowSampleClassification();
  void ShowSampleRegression();
  //взять значение переменной из выборки
  double GetValue(int Num,int Var);
  //взять значение выхода из выборки
  double GetOutput(int Num,int Var);
  //получить номер класса для измерения
  int GetClass(int Num);
  //разбиение выборки, кросс-валидация
  void SplitCVRandom();
  void SplitCVStratified();
  //считает число объектов по классам
  void ClassPatternsCalc();
  //простое разбиение выборки
  void SplitRandom();
  void SplitStratified();
  //общие параметры
  int Size;         //объем выборки
  int NCols;        //общее число столбцов в выборке
  int NVars;        //число столбцов входных параметров
  int NOuts;        //число столбцов выходных параметров
  int ProblemType;  //тип задачи
  int NFolds;       //число частей при кросс-валидации
  double SplitRate; //доля обучающей выборки, например
                    //0.7 => разбиение 70/30
  int LearnSize;    //размер обучающей выборки
  int TestSize;     //размер тестовой выборки

  double** Inputs;  //входы задачи
  double** Outputs; //выходы задачи
  bool** MissingInputs;    //массив пропущенных входных значений
  bool** MissingOutputs;   //массив пропущенных выходных значений
  int* CVSplit;     //определяет к какой части относится
                    //измерение при кросс-валидации
  int* FoldSize;    //размеры частей, на которые разбивается
                    //выборка при кросс-валидации
  int* CVFoldNum;   //номер части, к которой относится измерение

  //параметры для задач классификации
  int NClasses;     //число классов в задаче
  int* Classes;     //массив номеров классов
  int* NClassInst;  //количество объектов в классах
  int** ClassPositions;  //Номера объектов, принадлежащих разным классам
  int** ClassPerFold;    //число объектов классов для каждой части

};

#endif // SAMPLE_H
