#include "main.h"


int main(int argc, char* argv[])
{
    ofstream foutResult("Finished.txt");
    foutResult<<0;
    foutResult.close();
    //cout<<__FILE__<<" "<<__LINE__;
    srand((unsigned)time(NULL));
    clock_t StartTime = clock();
    int NumChanges = 100;
    float LearnPercent = 1.0;
    int SamplingType = 0;
    Params PRS;

    PRS.NumInds = 100;
    PRS.NumGens = 500;
    PRS.MaxNRules = 40;
    PRS.NFSets = 15;    //Change NPartitions if change this!
    PRS.NPartitions = 4;
    PRS.CritType = 0;
    PRS.TourSize = 2;
    PRS.NumSel = 2;
    PRS.NumCrs = 3;
    PRS.NumMut = 3;
    PRS.NumHeu = 2;
    PRS.HasGoodRules = false;
    PRS.DelProb = 0.5;
    PRS.AddProb = 0.5;
    PRS.SelType=0;
    PRS.CrossType=2;
    PRS.MutType=2;
    PRS.MichType=0;
    PRS.SADJ=1;
    PRS.K1=50;
    PRS.K2=8;
    PRS.CritType = 0;           //0 = Acc, 1 = Ave;
    PRS.Unbiased_init = 0;
    PRS.PDC = 0.5;
    PRS.formType = 1;
    int TotalNumRuns = 1;
    ifstream par("Params.txt");
    if(argc != 11)
    {
        par>>PRS.NumInds;
        par>>PRS.NumGens;
        par>>PRS.FNum;
        par>>PRS.MaxNRules;
        par>>PRS.CritType;
        par>>PRS.Unbiased_init;
        par>>NumChanges;
        par>>LearnPercent;
        par>>SamplingType;
        par>>TotalNumRuns;
        par>>ArtMis;
        cout<<"Read parameters from file, because argc = "<<argc<<endl;
    }
    else
    {
        PRS.NumInds = atoi(argv[1]);
        PRS.NumGens = atoi(argv[2]);
        PRS.FNum = atoi(argv[3]);
        PRS.MaxNRules = atoi(argv[4]);
        PRS.CritType = atoi(argv[5]);
        PRS.Unbiased_init = atoi(argv[6]);
        NumChanges = atoi(argv[7]);
        LearnPercent = atof(argv[8]);
        SamplingType = atoi(argv[9]);
        TotalNumRuns = atoi(argv[10]);
        ArtMis = atoi(argv[11]);
        cout<<"Handled parameters"<<endl;
        for(int i=0;i!=1;i++)
        {
            cout<<argv[i]<<endl;
        }
    }
    cout<<PRS.NumInds<<endl;
    cout<<PRS.NumGens<<endl;
    cout<<PRS.FNum<<endl;
    cout<<PRS.MaxNRules<<endl;
    cout<<PRS.CritType<<endl;
    cout<<PRS.Unbiased_init<<endl;
    cout<<NumChanges<<endl;
    cout<<LearnPercent<<endl;
    cout<<SamplingType<<endl;
    cout<<TotalNumRuns<<endl;
    cout<<"Parameters set!"<<endl;

    if(SamplingType == 0)
        LearnPercent = 1.0;

    ofstream fouterr("Errors.txt");
    ofstream fouterrU("ErrorsMass.txt");
    ofstream foutconf("ConfMatrix.txt");
    ofstream foutconf2("ConfMatrix2.txt");
    ofstream foutRS("RS.txt");

    ofstream fouterrnums("ErrNums.txt",ios::trunc);

    //ofstream foutbalance("BalanceRates.txt");

    //int TestN=0;
    //char buff[40];
    //sprintf(buff,"%d",TestN+1);
    //strcat(&buff[0],"_Results.txt");
	//ofstream fout(buff);

    for(int Run=0;Run!=TotalNumRuns;Run++)
    {
        sample S;
        switch(PRS.FNum)
        {
            case -19:
            {
                S.Init(5500,40,11,10,0.9);
                S.ReadFileClassification((char*)"texture_type.txt");
                break;
            }
            case -18:
            {
                S.Init(2310,19,7,10,0.9);
                S.ReadFileClassification((char*)"segment_n_type.txt");
                break;
            }
            case -17:
            {
                S.Init(5404,5,2,10,0.9);
                S.ReadFileClassification((char*)"phoneme_type.txt");
                break;
            }
            case -16:
            {
                S.Init(5472,10,5,10,0.9);
                S.ReadFileClassification((char*)"page-blocks_n_type.txt");
                break;
            }
            case -15:
            {
                S.Init(690,15,2,10,0.9);
                S.ReadFileClassification((char*)"crx.txt");
                break;
            }
            case -14:
            {
                S.Init(366,34,6,10,0.9);
                S.ReadFileClassification((char*)"dermatology.txt");
                break;
            }
            case -13:
            {
                S.Init(303,13,5,10,0.9);
                S.ReadFileClassification((char*)"cleveland.txt");
                break;
            }
            case -12:
            {
                S.Init(435,16,2,10,0.9);
                S.ReadFileClassification((char*)"housevotes.txt");
                break;
            }
            case -11:
            {
                S.Init(48788,1,7,10,0.9);
                S.ReadFileClassification((char*)"beeline_train_x8.txt");
                break;
            }
            case -10:
            {
                S.Init(50000,43,7,10,0.9);
                S.ReadFileClassification((char*)"beeline_train_cut.txt");
                break;
            }
            case -9:
            {
                S.Init(50000,62,7,10,0.9);
                S.ReadFileClassification((char*)"beeline_train.txt");
                break;
            }
            case -8:
            {
                S.Init(8993,13,9,10,0.9);
                S.ReadFileClassification((char*)"marketing.txt");
                break;
            }
            case -7:
            {
                S.Init(690,14,2,10,0.9);
                S.ReadFileClassification((char*)"aust_real.txt");
                break;
            }
            case -6:
            {
                S.Init(699,9,2,10,0.9);
                S.ReadFileClassification((char*)"breast_mis.txt");
                break;
            }
            case -5:
            {
                S.Init(10228,378,2,10,0.9);
                S.ReadFileClassification((char*)"Data.dat");
                break;
            }
            case -4:
            {
                S.Init(500000,28,2,10,0.9);
                S.ReadFileClassification((char*)"Spreadsheet2.txt");
                break;
            }
            case -3:
            {
                S.Init(5000,28,2,10,0.9);
                S.ReadFileClassification((char*)"Spreadsheet3.txt");
                break;
            }
            case -2:
            {
                S.Init(2308,19,2,10,0.9);
                S.ReadFileClassification((char*)"segment0.txt");
                break;
            }
            case -1:
            {
                S.Init(539,19,2,10,0.9);
                S.ReadFileClassification((char*)"bands.txt");
                break;
            }
            case 1:
            {
                S.Init(690,14,2,10,0.9);
                S.ReadFileClassification((char*)"aust_n.txt");
                break;
            }
            case 2:
            {
                S.Init(1000,24,2,10,0.9);
                S.ReadFileClassification((char*)"germ_n.txt");
                break;
            }
            case 3:
            {
                S.Init(2310,19,7,10,0.9);
                S.ReadFileClassification((char*)"segment_n.txt");
                break;
            }
            case 4:
            {
                S.Init(5404,5,2,10,0.9);
                S.ReadFileClassification((char*)"phoneme.txt");
                break;
            }
            case 5:
            {
                S.Init(5472,10,5,10,0.9);
                S.ReadFileClassification((char*)"page-blocks_n.txt");
                break;
            }
            case 6:
            {
                S.Init(6435,36,6,10,0.9);
                S.ReadFileClassification((char*)"satimage_n.txt");
                break;
            }
            case 7:
            {
                S.Init(7400,20,2,10,0.9);
                S.ReadFileClassification((char*)"twonorm.txt");
                break;
            }
            case 8:
            {
                S.Init(7400,20,2,10,0.9);
                S.ReadFileClassification((char*)"ring.txt");
                break;
            }
            case 9:
            {
                S.Init(10992,16,10,10,0.9);
                S.ReadFileClassification((char*)"penbased.txt");
                break;
            }
            case 10:
            {
                S.Init(19020,10,2,10,0.9);
                S.ReadFileClassification((char*)"magic.txt");
                break;
            }
            case 11:
            {
                S.Init(5500,40,11,10,0.9);
                S.ReadFileClassification((char*)"texture.txt");
                break;
            }
            case 12:
            {
                S.Init(2308,19,2,10,0.9);
                S.ReadFileClassification((char*)"segment0.txt");
                break;
            }
            case 13:
            {
                S.Init(2536,72,2,10,0.9);
                S.ReadFileClassification((char*)"ozone.txt");
                break;
            }
            case 14:
            {
                S.Init(5500000,28,2,10,0.9);
                S.ReadFileClassification((char*)"train_kaggle.txt");
                break;
            }
        }
        //S.ShowSampleClassification();
        S.ClassPatternsCalc();
        S.SplitCVStratified();

        //ofstream foutSha("ring_CV.txt");
        //for(int i=0;i!=S.Size;i++)
        //{
          //  foutSha<<S.GetCVFoldNum(i)<<"\t";
            //for(int j=0;j!=S.NVars;j++)
            //{
              //  foutSha<<S.GetValue(i,j)<<"\t";
            //}
          //  foutSha<<S.GetClass(i)<<endl;
        //}
        //exit(0);

        for(int FoldOnTest = 0;FoldOnTest!=10;FoldOnTest++)
        {
            StartTime = clock();
            S.NormalizeCV_01(FoldOnTest);
            //S.ShowNormSampleClassification();
            //S.ShowNormSampleClassification();
            //cout<<endl;
            PRS.NVars = S.NVars;

            float*** MRS; InitMRS(MRS,S,PRS);

            // задаем номера используемых из выборки измерений
            // и их число
            float OverBestError = S.Size;
            float OverBestAve = 1.0;

            float SavedError=-1;
            float SavedAve = -1;

            int NumInst;
            int* InstUsed;
            int* Inst;
            int* InstClass;
            int** CanBeUsedNums;
            //float LearnPercent = 1.0;
            int CVLearnSize = S.GetCVLearnSize(FoldOnTest);
            int CVTestSize = S.GetCVTestSize(FoldOnTest);
            float OverTestError = CVLearnSize;
            float OverTestAve = 1;
            NumInst = (int)S.GetCVLearnSize(FoldOnTest)*LearnPercent;
            //int NumChanges = (float)PRS.NumGens/20.;
            //int PrevLearnSize;
            //NumChanges = 100;
            //int SampIncreaseStep = S.GetCVLearnSize(FoldOnTest)*(1.-LearnPercent) / (PRS.NumGens/NumChanges);
            CanBeUsedNums = new int*[S.GetCVLearnSize(FoldOnTest)];
            for(int i=0;i!=S.GetCVLearnSize(FoldOnTest);i++)
            {
                CanBeUsedNums[i] = new int[S.NClasses];
            }
            //cout<<endl<<NumInst<<endl<<endl;
            Inst = new int[CVLearnSize];
            InstClass = new int[S.NClasses];
            for(int i=0;i!=S.NClasses;i++)
                InstClass[i] = 0;
            InstUsed = new int[S.Size];
            int* counterTestClass;
            counterTestClass = new int[S.NClasses];
            int* counterClass;
            counterClass = new int[S.NClasses];
            for(int i=0;i!=S.NClasses;i++)
            {
                counterClass[i] = 0;
                counterTestClass[i] = 0;
            }
            for(int i=0;i!=S.Size;i++)
            {
                InstUsed[i] = 1;
                if(S.CVFoldNum[i] != FoldOnTest)
                {
                    CanBeUsedNums[counterClass[S.Classes[i]]] [S.Classes[i]] = i;
                    //cout<<CanBeUsedNums[counterClass[S.Classes[i]]] [S.Classes[i]]<<" ";
                    counterClass[S.Classes[i]] ++;
                }
                if(S.CVFoldNum[i] == FoldOnTest)
                {
                    counterTestClass[S.Classes[i]] ++;
                }
            }
            //cout<<endl;
            if(SamplingType == 0)
                SetAllInst(S,PRS,NumInst,InstUsed,Inst,InstClass,FoldOnTest,CanBeUsedNums,counterClass,CVLearnSize);
            else if(SamplingType == 1)
                SetAvailableInst(S,PRS,NumInst,InstUsed,Inst,InstClass,FoldOnTest,CanBeUsedNums,counterClass,CVLearnSize);
            else if(SamplingType == 2)
                SetAvailableInstBalanced(S,PRS,NumInst,InstUsed,Inst,InstClass,FoldOnTest,CanBeUsedNums,counterClass,CVLearnSize);

            cout<<"Sample has been set"<<endl;
            for(int i=0;i!=NumInst;i++)
            {
                //cout<<InstUsed[i]<<"\t";
            }
            //cout<<endl;

            //////////////////////////

            float* RuleClassified;
            float* RuleP;
            float* tempRuleC;
            float* BadRNums;
            int* MisClassifiedNums;
            int* RandomNumbers;
            int* CurNumRules;
            float* CurRLength;
            int** GoodRNums;
            float** ChosenRules;
            float* FitMass;
            float* FitMassCopy;
            float* PropSetSel;
            float* ClassConf;
            float** PopulClass;
            float** PopulCF;
            float*** PopulMu;
            int** BestInd;
            float* BestClasses;
            float* BestCF;
            float** BestMu;
            int** OverBestInd;
            int** OverBestInd2;
            float* OverBestClasses;
            float* OverBestClasses2;
            float* OverBestCF;
            float* OverBestCF2;
            float** OverBestMu;
            float* OverErrorClass;
            float** RUPD;
            float* ErrMass;
            int* SelectedNums1;
            int* SelectedNums2;
            float* NumMass;
            float* PMass;
            float* NumMass2;
            float* RankMass;
            float* RMass;
            float* SavedFit;
            float* NumClassInst;
            float* NumClassInstWrong;
            float* NumClassInstTest;
            float* NumClassInstTestWrong;
            float** ConfMatrix;
            float** ConfMatrixTest;
            float** ErrorClass;
            float* AveMass;
            float* TestErrorClass;
            int*** Popul;
            float* psel;
            float* pcrs;
            float* pmut;
            float* pheu;

            float* usedsel;
            float* usedcrs;
            float* usedmut;
            float* usedheu;

            float* sumfitsel;
            float* sumfitcrs;
            float* sumfitmut;
            float* sumfitheu;

            float* numselused;
            float* numcrsused;
            float* nummutused;
            float* numheuused;

            float* sucsel;
            float* succrs;
            float* sucmut;
            float* sucheu;

            OverBestInd = new int*[PRS.MaxNRules];
            for(int i=0;i!=PRS.MaxNRules;i++)
            {
                OverBestInd[i] = new int[S.NVars];
            }
            OverBestClasses = new float[PRS.MaxNRules];
            OverBestCF = new float[PRS.MaxNRules];
            OverBestInd2 = new int*[PRS.MaxNRules];
            for(int i=0;i!=PRS.MaxNRules;i++)
            {
                OverBestInd2[i] = new int[S.NVars];
            }
            OverBestClasses2 = new float[PRS.MaxNRules];
            OverBestCF2 = new float[PRS.MaxNRules];
            OverBestMu = new float*[PRS.MaxNRules];
            for(int i=0;i!=PRS.MaxNRules;i++)
            {
                OverBestMu[i] = new float[S.Size];
            }
            OverErrorClass = new float[S.NClasses];

            Init(S,PRS,
                 RuleClassified, RuleP, tempRuleC, BadRNums, MisClassifiedNums, RandomNumbers, CurNumRules,CurRLength,
                 GoodRNums, ChosenRules, FitMass, FitMassCopy, PropSetSel, ClassConf, PopulClass, PopulCF, PopulMu,BestInd,
                 BestClasses, BestCF, BestMu, RUPD, ErrMass, SelectedNums1, SelectedNums2, NumMass, PMass, NumMass2, RankMass,
                 RMass, SavedFit, NumClassInst, NumClassInstWrong, NumClassInstTest, NumClassInstTestWrong, ConfMatrix,
                 ConfMatrixTest, ErrorClass, AveMass, TestErrorClass, Popul, psel, pcrs, pmut, pheu,usedsel,usedcrs,usedmut,
                 usedheu,sumfitsel,sumfitcrs,sumfitmut,sumfitheu,numselused,numcrsused,nummutused,numheuused,sucsel,succrs,
                 sucmut,sucheu,NumInst,InstUsed,Inst,InstClass,FoldOnTest,MRS,CVLearnSize);

            FindNSaveBest(0,PRS,NumMass,FitMassCopy,Popul,PopulClass,PopulMu,PopulCF,RUPD,AveMass,ErrMass,FitMass,CurNumRules,
                          BestInd,BestClasses,BestCF,BestMu,NumInst,CurRLength,true,OverBestInd,OverBestClasses,OverBestCF,
                          OverBestError,OverBestInd2,OverBestClasses2,OverBestCF2,S,OverBestMu,ClassConf,
                          CVLearnSize,MRS,OverBestAve,OverErrorClass,FoldOnTest,NumClassInst,ErrorClass,ConfMatrix,
                          Inst,InstUsed,InstClass,counterClass);

            for(int i=0;i!=PRS.NumGens+1;i++)
            {
                if(i%NumChanges != 0 || i == 0)
                //if(false)
                {
                    CreateNewPopul(i,PRS,S,RuleClassified, RuleP, tempRuleC, BadRNums, MisClassifiedNums, RandomNumbers, CurNumRules,
                               CurRLength,GoodRNums, ChosenRules, FitMass, FitMassCopy, PropSetSel, ClassConf, PopulClass,
                               PopulCF, PopulMu,BestInd, BestClasses, BestCF, BestMu, RUPD, ErrMass, SelectedNums1,
                               SelectedNums2, NumMass, PMass, NumMass2, RankMass, RMass, SavedFit, NumClassInst,
                               NumClassInstWrong, NumClassInstTest, NumClassInstTestWrong, ConfMatrix, ConfMatrixTest,
                               ErrorClass, AveMass, TestErrorClass, Popul, psel, pcrs, pmut, pheu,usedsel,usedcrs,usedmut,
                               usedheu,sumfitsel,sumfitcrs,sumfitmut,sumfitheu,numselused,numcrsused,nummutused,numheuused,
                               sucsel,succrs, sucmut,sucheu,NumInst,InstUsed,Inst,InstClass,FoldOnTest,MRS,false,
                               OverBestInd,OverBestClasses,OverBestCF,OverBestError,OverBestInd2,OverBestClasses2,OverBestCF2,
                               OverBestMu,CVLearnSize,OverBestAve,OverErrorClass,counterClass);

                }
                else
                {
                    CreateNewPopul(i,PRS,S,RuleClassified, RuleP, tempRuleC, BadRNums, MisClassifiedNums, RandomNumbers, CurNumRules,
                               CurRLength,GoodRNums, ChosenRules, FitMass, FitMassCopy, PropSetSel, ClassConf, PopulClass,
                               PopulCF, PopulMu,BestInd, BestClasses, BestCF, BestMu, RUPD, ErrMass, SelectedNums1,
                               SelectedNums2, NumMass, PMass, NumMass2, RankMass, RMass, SavedFit, NumClassInst,
                               NumClassInstWrong, NumClassInstTest, NumClassInstTestWrong, ConfMatrix, ConfMatrixTest,
                               ErrorClass, AveMass, TestErrorClass, Popul, psel, pcrs, pmut, pheu,usedsel,usedcrs,usedmut,
                               usedheu,sumfitsel,sumfitcrs,sumfitmut,sumfitheu,numselused,numcrsused,nummutused,numheuused,
                               sucsel,succrs, sucmut,sucheu,NumInst,InstUsed,Inst,InstClass,FoldOnTest,MRS,true,
                               OverBestInd,OverBestClasses,OverBestCF,OverBestError,OverBestInd2,OverBestClasses2,OverBestCF2,
                               OverBestMu,CVLearnSize,OverBestAve,OverErrorClass,counterClass);
                    for(int i=0;i!=CVLearnSize;i++)
                    {
                        //if(InstUsed[i] < 1000)
                            //cout<<InstUsed[i]<<" ";
                    }
                    //cout<<endl;
                    FitCalc(0,S,ErrorClass,ConfMatrix,NumInst,PRS,PopulCF,PopulMu,PopulClass,CurNumRules,CurRLength,
                        ErrMass,AveMass,FitMass,NumClassInst,Popul,Inst,1,InstUsed,InstClass,RUPD,MRS,ClassConf);
                    for(int i=0;i!=CVLearnSize;i++)
                    {
                        //if(InstUsed[i] < 1000)
                            //cout<<InstUsed[i]<<" ";
                    }
                    //cout<<endl;
                    /*for(int L=0;L!=PRS.NumInds;L++)
                    {
                        float SavedOverBestError = OverBestError;
                        for(int j=0;j!=PRS.MaxNRules;j++)
                        {
                            for(int k=0;k!=PRS.NVars;k++)
                            {
                                OverBestInd2[j][k] = OverBestInd[j][k];
                            }
                            OverBestClasses2[j] = OverBestClasses[j];
                            OverBestCF2[j] = OverBestCF[j];
                        }
                        for(int j=0;j!=PRS.MaxNRules;j++)
                        {
                            for(int k=0;k!=PRS.NVars;k++)
                            {
                                OverBestInd[j][k] = Popul[L][j][k];
                            }
                            OverBestClasses[j] = PopulClass[L][j];
                            OverBestCF[j] = PopulCF[L][j];
                        }
                        LearnErrorCalc(S,PRS,NumClassInst,OverBestInd,OverBestClasses,OverBestCF,OverBestMu,ClassConf,
                                       CVLearnSize,MRS,OverBestAve,OverErrorClass,OverBestError,FoldOnTest);
                        //cout<<OverBestError<<"\t";
                        //cout<<endl;
                        if(OverBestError > SavedOverBestError)
                        {
                            OverBestError = SavedOverBestError;
                            for(int j=0;j!=PRS.MaxNRules;j++)
                            {
                                for(int k=0;k!=PRS.NVars;k++)
                                {
                                    OverBestInd[j][k] = OverBestInd2[j][k];
                                }
                                OverBestClasses[j] = OverBestClasses2[j];
                                OverBestCF[j] = OverBestCF2[j];
                            }
                        }
                    }*/
                    if(PRS.CritType == 1)
                    {
                        if(SavedAve == -1)
                        {
                            SavedAve = OverBestAve;
                        }
                        else
                        {
                            if(OverBestAve >= SavedAve)
                            {
                                //PrevLearnSize = NumInst;
                                //NumInst+=CVLearnSize*0.02;
                                if(NumInst > CVLearnSize)
                                    NumInst = CVLearnSize;
                            }
                            else
                            {
                                SavedAve = OverBestAve;
                            }
                        }
                    }
                    else
                    {
                        if(SavedError == -1)
                        {
                            SavedError = OverBestError;
                        }
                        else
                        {
                            if(OverBestError >= SavedError)
                            {
                                //PrevLearnSize = NumInst;
                                //NumInst+=CVLearnSize*0.02;
                                if(NumInst > CVLearnSize)
                                    NumInst = CVLearnSize;
                            }
                            else
                            {
                                SavedError = OverBestError;
                            }
                        }
                    }

                    //PrevLearnSize = NumInst;
                    //NumInst+=SampIncreaseStep;
                    //if(NumInst > CVLearnSize)
                        //NumInst = CVLearnSize;
                    //NumInst = (int)S.GetCVLearnSize(FoldOnTest)*LearnPercent;
                    //cout<<endl<<NumInst<<endl<<endl;
                    cout<<endl<<"NewNumInst = "<<NumInst<<endl;
                    //if(NumInst != PrevLearnSize)
                    if(SamplingType == 0)
                        SetAllInst(S,PRS,NumInst,InstUsed,Inst,InstClass,FoldOnTest,CanBeUsedNums,counterClass,CVLearnSize);
                    else if(SamplingType == 1)
                        SetAvailableInst(S,PRS,NumInst,InstUsed,Inst,InstClass,FoldOnTest,CanBeUsedNums,counterClass,CVLearnSize);
                    else if(SamplingType == 2)
                        SetAvailableInstBalanced(S,PRS,NumInst,InstUsed,Inst,InstClass,FoldOnTest,CanBeUsedNums,counterClass,CVLearnSize);
                    for(int j=0;j!=PRS.NumInds;j++)
                    {
                        for(int k=0;k!=PRS.MaxNRules;k++)
                        {
                            RUPD[j][k] = 0;
                            CFClassRUPD(j,k,PRS,NumInst,Inst,Popul,RUPD,MRS,S,ClassConf,PopulMu,
                                        PopulClass,PopulCF,InstClass);
                        }
                        CFClassRB(j,S,PRS,RUPD,Popul,ClassConf,PopulMu,MRS,PopulClass,PopulCF,NumInst,Inst,InstClass);
                        FitCalc(j,S,ErrorClass,ConfMatrix,NumInst,PRS,PopulCF,PopulMu,PopulClass,CurNumRules,CurRLength,
                                ErrMass,AveMass,FitMass,NumClassInst,Popul,Inst,0,InstUsed,InstClass,RUPD,MRS,ClassConf);
                    }
                    FindNSaveBest(0,PRS,NumMass,FitMassCopy,Popul,PopulClass,PopulMu,PopulCF,RUPD,AveMass,ErrMass,FitMass,CurNumRules,
                                BestInd,BestClasses,BestCF,BestMu,NumInst,CurRLength,true,OverBestInd,OverBestClasses,OverBestCF,
                                OverBestError,OverBestInd2,OverBestClasses2,OverBestCF2,S,OverBestMu,ClassConf,
                                CVLearnSize,MRS,OverBestAve,OverErrorClass,FoldOnTest,NumClassInst,ErrorClass,ConfMatrix,
                                Inst,InstUsed,InstClass,counterClass);
                    for(int j=0;j!=NumInst;j++)
                    {
                        //cout<<InstUsed[Inst[j]]<<" ";
                    }
                    //cout<<endl;
                    //qSort(Inst,0,CVLearnSize-1);
                    //cout<<endl;
                    for(int j=0;j!=CVLearnSize;j++)
                    {
                        //cout<<i<<" "<<Inst[i]<<"\n";
                    }
                    //cout<<endl;
                }
                if(PRS.NumGens-1 == i)
                {
                    havetosaveerrors = true;
                    for(int KL=0;KL!=S.Size;KL++)
                        S.ErrOnMiss[KL] = 0;
                    S.NMisErr = 0;
                    S.NNormErr = 0;
                }
                LearnErrorCalc(S,PRS,NumClassInst,OverBestInd,OverBestClasses,OverBestCF,OverBestMu,ClassConf,
                               CVLearnSize,MRS,OverBestAve,OverErrorClass,OverBestError,FoldOnTest,counterClass,ConfMatrix);
                //OverTestError = 0;
                //OverTestAve = 0;
                TestErrorCalc(S,PRS,NumClassInst,OverBestInd,OverBestClasses,OverBestCF,OverBestMu,ClassConf,
                               CVTestSize,MRS,OverTestAve,OverErrorClass,OverTestError,FoldOnTest,counterTestClass,ConfMatrixTest);
                if(PRS.NumGens-1 == i)
                {
                    havetosaveerrors = false;
                    for(int KL=0;KL!=S.Size;KL++)
                    {
                        if(S.ErrOnMiss[KL] == 1 && S.HasMiss[KL] == 1)
                            S.NMisErr++;
                        if(S.ErrOnMiss[KL] == 1 && S.HasMiss[KL] == 0)
                            S.NNormErr++;
                    }
                    fouterrnums<<S.NMisErr<<"\t"<<S.NNormErr<<endl;
                }

                if(PRS.CritType == 0)
                    cout<<OverBestError<<" "<<OverTestError<<endl;
                else
                    cout<<OverBestAve<<" "<<OverTestAve<<endl;
                //cout<<endl;
                //foutbalance<<ErrorClass[0][0]<<"\t"<<ErrorClass[0][1]<<endl;
                //_getch();
                fouterrU<<i<<"\t";
                fouterrU<<1.-(float)PRS.BestErr/(float)NumInst<<"\t"<<PRS.BestErr<<"\t"<<NumInst<<"\t"<<PRS.BestNumRules<<"\t"<<PRS.BestRLength<<"\t";
                fouterrU<<1.-(float)OverBestError/(float)CVLearnSize<<"\t"<<OverBestError<<"\t"<<CVLearnSize<<"\t"<<1.-(float)OverTestError/(float)CVTestSize<<"\t";
                fouterrU<<OverTestError<<"\t"<<CVTestSize<<"\t"<<1.-OverBestAve<<"\t"<<1.-OverTestAve<<"\t";
                fouterrU<<float(clock()-StartTime)/(float)CLOCKS_PER_SEC<<endl;

                foutconf<<i<<endl;
                for(int L=0;L!=S.NClasses;L++)
                {
                    for(int LL=0;LL!=S.NClasses+1;LL++)
                    {
                        foutconf<<ConfMatrix[L][LL]<<"\t";
                    }
                    foutconf<<endl;
                }
                foutconf2<<i<<endl;
                for(int L=0;L!=S.NClasses;L++)
                {
                    for(int LL=0;LL!=S.NClasses+1;LL++)
                    {
                        foutconf2<<ConfMatrixTest[L][LL]<<"\t";
                    }
                    foutconf2<<endl;
                }
                foutRS<<i<<endl;
                for(int i=0;i!=PRS.MaxNRules;i++)
                {
                    if(OverBestClasses[i] != -1)
                    {
                        for(int j=0;j!=PRS.NVars;j++)
                        {
                            foutRS<<OverBestInd[i][j]<<"\t";
                        }
                        foutRS<<"->\t"<<OverBestClasses[i]<<"\t"<<OverBestCF[i];
                        foutRS<<endl;
                    }
                }
            }
            FitCalc(0,S,ErrorClass,ConfMatrix,NumInst,PRS,PopulCF,PopulMu,PopulClass,CurNumRules,CurRLength,
                        ErrMass,AveMass,FitMass,NumClassInst,Popul,Inst,1,InstUsed,InstClass,RUPD,MRS,ClassConf);
            for(int i=0;i!=NumInst;i++)
            {
                cout<<InstUsed[Inst[i]]<<"\t";
            }
            cout<<endl;

            fouterr<<1.-(float)PRS.BestErr/(float)NumInst<<"\t"<<PRS.BestErr<<"\t"<<NumInst<<"\t"<<PRS.BestNumRules<<"\t"<<PRS.BestRLength<<"\t";
            fouterr<<1.-(float)OverBestError/(float)CVLearnSize<<"\t"<<OverBestError<<"\t"<<CVLearnSize<<"\t"<<1.-(float)OverTestError/(float)CVTestSize<<"\t";
            fouterr<<OverTestError<<"\t"<<CVTestSize<<"\t"<<1.-OverBestAve<<"\t"<<1.-OverTestAve<<"\t";
            //fouterr<<1.-(float)TestErr/(float)TestSize<<"\t"<<TestErr<<"\t"<<TestSize<<"\t"<<TestNumRules<<"\t"<<TestRLength<<"\t";
            //fouterr<<1.-BestAve<<"\t"<<1.-TestAve<<"\t";
            /*for(int L=0;L!=CVS.NClasses;L++)
            {
                fouterr<<NumClassInstTest[L]<<"\t"<<NumClassInstTestWrong[L]<<"\t";
                fouterr<<NumClassInst[L]<<"\t"<<ErrorClass[0][L]*NumClassInst[L]<<"\t";
            }*/
            fouterr<<float(clock()-StartTime)/(float)CLOCKS_PER_SEC<<endl;
            //avgerr+=1.-(float)BestErr/(float)LearnSize;
            //avgerrtest+=1.-(float)TestErr/(float)TestSize;
            //cout<<endl;
            for(int i=0;i!=PRS.MaxNRules;i++)
            {
                if(OverBestClasses[i] != -1)
                {
                    for(int j=0;j!=PRS.NVars;j++)
                    {
                        cout<<OverBestInd[i][j]<<" ";
                    }
                    cout<<"-> "<<OverBestClasses[i]<<" "<<OverBestCF[i];
                    cout<<endl;
                }
            }
            cout<<1.-(float)PRS.BestErr/(float)NumInst<<"\t"<<endl;

            Clean(S,PRS,
                 RuleClassified, RuleP, tempRuleC, BadRNums, MisClassifiedNums, RandomNumbers, CurNumRules,CurRLength,
                 GoodRNums, ChosenRules, FitMass, FitMassCopy, PropSetSel, ClassConf, PopulClass, PopulCF, PopulMu,BestInd,
                 BestClasses, BestCF, BestMu, RUPD, ErrMass, SelectedNums1, SelectedNums2, NumMass, PMass, NumMass2, RankMass,
                 RMass, SavedFit, NumClassInst, NumClassInstWrong, NumClassInstTest, NumClassInstTestWrong, ConfMatrix,
                 ConfMatrixTest, ErrorClass, AveMass, TestErrorClass, Popul, psel, pcrs, pmut, pheu,usedsel,usedcrs,usedmut,
                 usedheu,sumfitsel,sumfitcrs,sumfitmut,sumfitheu,numselused,numcrsused,nummutused,numheuused,sucsel,succrs,
                 sucmut,sucheu,NumInst,InstUsed,Inst,InstClass,FoldOnTest,MRS);

            for(int i=0;i!=PRS.MaxNRules;i++)
            {
                delete OverBestMu;
            }
            delete OverBestMu;
            for(int i=0;i!=PRS.MaxNRules;i++)
            {
                delete OverBestInd[i];
            }
            delete OverBestInd;
            delete OverBestClasses;
            delete OverBestCF;
            for(int i=0;i!=PRS.MaxNRules;i++)
            {
                delete OverBestInd2[i];
            }
            delete OverBestInd2;
            delete OverBestClasses2;
            delete OverBestCF2;
            delete OverErrorClass;

            //////////////////////////

            for(int i=0;i!=S.GetSize();i++)
            {
                for(int j=0;j!=S.GetNVars();j++)
                {
                    delete MRS[i][j];
                }
                delete MRS[i];
            }
            delete MRS;
            delete InstClass;
            delete Inst;
            delete InstUsed;
            for(int i=0;i!=S.GetCVLearnSize(FoldOnTest);i++)
                delete CanBeUsedNums[i];
            delete CanBeUsedNums;
            delete counterClass;
            delete counterTestClass;
        }
        S.CleanSamp();
        //fouterrU<<endl;
    }

    cout<<endl<<"Total time: "<<float(clock()-StartTime)/(float)CLOCKS_PER_SEC<<" seconds"<<endl;
    //_getch();
    foutResult.open("Finished.txt");
    foutResult<<1;
    foutResult.close();
    return 0;
}
