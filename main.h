#include "sample.cpp"
#include <conio.h>

struct Params
{
    int NumInds;
    int NumGens;
    int MaxNRules;
    int NumVars;
    float PDC;
    float BestFit;
    float BestErr;
    float BestNumRules;
    float BestRLength;
    float TestFit;
    float TestErr;
    float TestNumRules;
    float TestRLength;
    bool formType;
    bool AddDelete;
    float DelProb;
    float AddProb;
    float BestAve;
    float TestAve;
    int CritType;
    int Unbiased_init;
    int TourSize;
    int NCols;
    int NVars;
    int NOuts;
    int LearnSize;
    int TestSize;
    int SampleSize;
    int NFSets;
    int NPartitions;
    int NMisclassified;
    bool HasGoodRules;
    int SelType;
    int CrossType;
    int MutType;
    int MichType;
    int NumSel;
    int NumCrs;
    int NumMut;
    int NumHeu;
    float K1;
    float K2;
    int savednumheu;
    int savednumgen;
    //CRS = MICHIGAN
    int SADJ;
    bool HasBestFit;
    int FNum;
    float pselall;
    float pcrsall;
    float pmutall;
    float pheuall;
};
float SumEl(float* Mass,int N)
{
	float Sum=0;
	for(int i=0;i!=N;i++)
	{
		Sum=Sum+Mass[i];
	}
	return Sum;
}
float MinElVal(float* Mass,int N)
{
	float min = Mass[0];
	for(int i=1;i!=N;i++)
	{
		if(Mass[i]<min)
			min=Mass[i];
	}
	return min;
}
float MaxElVal(float* Mass,int N)
{
	float max = Mass[0];
	for(int i=1;i!=N;i++)
	{
		if(Mass[i]>max)
			max=Mass[i];
	}
	return max;
}
int SumEl(int* Mass,int N)
{
	int Sum=0;
	for(int i=0;i!=N;i++)
	{
		Sum=Sum+Mass[i];
	}
	return Sum;
}
int MinElVal(int* Mass,int N)
{
    int min = Mass[0];
	for(int i=1;i!=N;i++)
	{
		if(Mass[i]<min)
			min=Mass[i];
	}
	return min;
}
int MaxElVal(int* Mass,int N)
{
	int max = Mass[0];
	for(int i=1;i!=N;i++)
	{
		if(Mass[i]>max)
			max=Mass[i];
	}
	return max;
}
void qSort(float* Mass,float* Mass2, int low, int high)
{
    int i=low;
    int j=high;
    float x=Mass[(low+high)>>1];
    do
    {
        while(Mass[i]<x)    ++i;
        while(Mass[j]>x)    --j;
        if(i<=j)
        {
            float temp=Mass[i];
            Mass[i]=Mass[j];
            Mass[j]=temp;
            float temp2=Mass2[i];
            Mass2[i]=Mass2[j];
            Mass2[j]=temp2;
            i++;    j--;
        }
    } while(i<=j);
    if(low<j)   qSort(Mass,Mass2,low,j);
    if(i<high)  qSort(Mass,Mass2,i,high);
}
void qSort(int* Mass, int low, int high)
{
    int i=low;
    int j=high;
    float x=Mass[(low+high)>>1];
    do
    {
        while(Mass[i]<x)    ++i;
        while(Mass[j]>x)    --j;
        if(i<=j)
        {
            int temp=Mass[i];
            Mass[i]=Mass[j];
            Mass[j]=temp;
            i++;    j--;
        }
    } while(i<=j);
    if(low<j)   qSort(Mass,low,j);
    if(i<high)  qSort(Mass,i,high);
}
void QS(int type, float* &RankMass, float* &NumMass2, Params &PRS)
{
    if(type == 0)
        qSort(RankMass,NumMass2,0,PRS.NumInds-1);
    else
        qSort(NumMass2,RankMass,0,PRS.NumInds-1);
}
float GetMR(float X, int FSNum)
{
    switch(FSNum)
    {
        case 0:
            {
                return -1.;
            }
        case 1:
            {
                if(X<0.)
                    return 1;
                if(X>1.)
                    return 0;
                return 1.-X;
            }
        case 2:
            {
                if(X<0.)
                    return 0;
                if(X>1.)
                    return 1;
                return X;
            }
        case 3:
            {
                if(X<0.)
                    return 1;
                if(X>0.5)
                    return 0.;
                else
                    return 1.-X*2.;
            }
        case 4:
            {
                if(X<0.)
                    return 0;
                if(X>1.)
                    return 0;
                if(X>0.5)
                    return 2.-X*2.;
                else
                    return X*2.;
            }
        case 5:
            {
                if(X>1.)
                    return 1;
                if(X>0.5)
                    return 2.*X-1.;
                else
                    return 0.;
            }
        case 6:
            {
                if(X<0.)
                    return 1;
                if(X>1./3.)
                    return 0.;
                else
                    return 1.-3.*X;
            }
        case 7:
            {
                if(X<0.)
                    return 0;
                if(X>1.)
                    return 0;
                if(X<1./3.)
                    return 3.*X;
                if(X>2./3.)
                    return 0.;
                else
                    return 2.-3.*X;
            }
        case 8:
            {
                if(X<0.)
                    return 0;
                if(X>1.)
                    return 0;
                if(X<1./3.)
                    return 0.;
                if(X>2./3.)
                    return 3.-3.*X;
                else
                    return 3*X-1.;
            }
        case 9:
            {
                if(X>1.)
                    return 1;
                if(X<2./3.)
                    return 0.;
                else
                    return 3.*X-2.;
            }
        case 10:
            {
                if(X<0.)
                    return 1;
                if(X>0.25)
                    return 0.;
                else
                    return 1.-4.*X;
            }
        case 11:
            {
                if(X<0.)
                    return 0;
                if(X>1.)
                    return 0;
                if(X<0.25)
                    return 4.*X;
                if(X>0.5)
                    return 0.;
                else
                    return 2.-4.*X;
            }
        case 12:
            {
                if(X<0.)
                    return 0;
                if(X>1.)
                    return 0;
                if(X<0.25)
                    return 0.;
                if(X<0.5)
                    return 4.*X-1.;
                if(X>0.75)
                    return 0.;
                else
                    return 3.-4.*X;
            }
        case 13:
            {
                if(X<0.)
                    return 0;
                if(X>1.)
                    return 0;
                if(X<0.5)
                    return 0.;
                if(X>0.75)
                    return 4.-4.*X;
                else
                    return 4.*X-2.;
            }
        case 14:
            {
                if(X>1.)
                    return 1;
                if(X<0.75)
                    return 0.;
                else
                    return 4.*X-3.;
            }
        default:
            {
                cout<<"Something went wrong...";
            }
    }
    return -1;
}
void InitMRS(float*** &MRS, sample &Samp, Params &PRS)
{
    MRS = new float**[Samp.GetSize()];
    for(int i=0;i!=Samp.GetSize();i++)
    {
        MRS[i] = new float*[Samp.GetNVars()];
        for(int j=0;j!=Samp.GetNVars();j++)
        {
            MRS[i][j] = new float[PRS.NFSets];
            for(int k=0;k!=PRS.NFSets;k++)
            {
                MRS[i][j][k] = GetMR(Samp.GetNormValue(i,j),k);
                //cout<<MRS[i][j][k]<<"\t";
            }
            //cout<<endl;
        }
        //cout<<endl;
    }
}
/*void InitMRS2(float*** &MRS1,float*** &MRS2,float*** &MRS3,float*** &MRS4,float*** &MRS5, sample &Samp, Params &PRS)
{
    int SampPartSize = Samp.GetSize()/5;
    MRS1 = new float**[SampPartSize];
    MRS2 = new float**[SampPartSize];
    MRS3 = new float**[SampPartSize];
    MRS4 = new float**[SampPartSize];
    MRS5 = new float**[SampPartSizes];
    for(int i=0;i!=SampPartSize;i++)
    {
        MRS1[i] = new float*[Samp.GetNVars()];
        for(int j=0;j!=Samp.GetNVars();j++)
        {
            MRS1[i][j] = new float[PRS.NFSets];
            for(int k=0;k!=PRS.NFSets;k++)
            {
                MRS1[i][j][k] = GetMR(Samp.GetNormValue(i,j),k);
                //cout<<MRS[i][j][k]<<"\t";
            }
            //cout<<endl;
        }
        //cout<<endl;
    }
    for(int i=SampPartSize;i!=SampPartSize*2;i++)
    {
        MRS2[i] = new float*[Samp.GetNVars()];
        for(int j=0;j!=Samp.GetNVars();j++)
        {
            MRS2[i][j] = new float[PRS.NFSets];
            for(int k=0;k!=PRS.NFSets;k++)
            {
                MRS2[i][j][k] = GetMR(Samp.GetNormValue(i,j),k);
                //cout<<MRS[i][j][k]<<"\t";
            }
            //cout<<endl;
        }
        //cout<<endl;
    }
    for(int i=SampPartSize*2;i!=SampPartSize*3;i++)
    {
        MRS3[i] = new float*[Samp.GetNVars()];
        for(int j=0;j!=Samp.GetNVars();j++)
        {
            MRS3[i][j] = new float[PRS.NFSets];
            for(int k=0;k!=PRS.NFSets;k++)
            {
                MRS3[i][j][k] = GetMR(Samp.GetNormValue(i,j),k);
                //cout<<MRS[i][j][k]<<"\t";
            }
            //cout<<endl;
        }
        //cout<<endl;
    }
    for(int i=SampPartSize*3;i!=SampPartSize*4;i++)
    {
        MRS4[i] = new float*[Samp.GetNVars()];
        for(int j=0;j!=Samp.GetNVars();j++)
        {
            MRS4[i][j] = new float[PRS.NFSets];
            for(int k=0;k!=PRS.NFSets;k++)
            {
                MRS4[i][j][k] = GetMR(Samp.GetNormValue(i,j),k);
                //cout<<MRS[i][j][k]<<"\t";
            }
            //cout<<endl;
        }
        //cout<<endl;
    }
    for(int i=SampPartSize*4;i!=SampPartSize*5;i++)
    {
        MRS5[i] = new float*[Samp.GetNVars()];
        for(int j=0;j!=Samp.GetNVars();j++)
        {
            MRS5[i][j] = new float[PRS.NFSets];
            for(int k=0;k!=PRS.NFSets;k++)
            {
                MRS5[i][j][k] = GetMR(Samp.GetNormValue(i,j),k);
                //cout<<MRS[i][j][k]<<"\t";
            }
            //cout<<endl;
        }
        //cout<<endl;
    }
}*/
void SetAllInst(sample &Samp,Params &PRS, int &NumInst, int* &InstUsed, int* &Inst, int* &InstClass,
                      int &FoldOnTest, int** &CanBeUsedNums, int* &counterClass, int& CVLearnSize)
{
    int counter = 0;
    for(int i=0;i!=Samp.NClasses;i++)
    {
        InstClass[i] = counterClass[i];
        for(int j=0;j!=counterClass[i];j++)
        {
            Inst[counter] = CanBeUsedNums[j][i];
            counter++;
        }
    }
}
void SetAvailableInstBalanced(sample &Samp,Params &PRS, int &NumInst, int* &InstUsed, int* &Inst, int* &InstClass,
                      int &FoldOnTest, int** &CanBeUsedNums, int* &counterClass, int& CVLearnSize)
{
    //cout<<endl;
    int TotalNumInst = 0;
    //int BigestClassNum = 0;
    for(int i=0;i!=Samp.NClasses;i++)
    {
        InstClass[i] = 0;
    }
    do
    {
        for(int i=0;i!=Samp.NClasses;i++)
        {
            if(InstClass[i] != counterClass[i])
            {
                InstClass[i] ++;
                TotalNumInst ++;
            }
        }
    } while(TotalNumInst < NumInst);

    /*
    for(int i=0;i!=Samp.NClasses;i++)
    {
        InstClass[i] = (int)round((float)NumInst * (float)counterClass[i] / (float)CVLearnSize);
        TotalNumInst += InstClass[i];
    }
    while(TotalNumInst > NumInst)
    {
        BigestClassNum = 0;
        for(int i=0;i!=Samp.NClasses;i++)
            if(InstClass[i] >= BigestClassNum)
                BigestClassNum = i;
        InstClass[BigestClassNum] --;
        TotalNumInst --;
    }
    while(TotalNumInst < NumInst)
    {
        int NumMin = InstClass[0];
        int MinPos = 0;
        for(int i=0;i!=Samp.NClasses;i++)
            if(NumMin > InstClass[i])
            {
                NumMin = InstClass[i];
                MinPos = i;
            }
        InstClass[MinPos]++;
        TotalNumInst++;
    }
    */
    for(int i=0;i!=Samp.NClasses;i++)
    {
        //cout<<InstClass[i]<<"\t";
    }
    float RandomPoint;
    float PSumm;
    float* PMass;
    PMass = new float[CVLearnSize];
    int TotalInstCounter = 0;
    int PreviousClassEnd = 0;
    for(int i=0;i!=Samp.NClasses;i++)
    {
        PMass[0] = 1./(float)InstUsed[CanBeUsedNums[0][i]];
        PSumm = PMass[0];
        for(int j=1;j!=counterClass[i];j++)
        {
            PMass[j] = PMass[j-1] + 1./(float)InstUsed[CanBeUsedNums[j][i]];
            PSumm += 1./(float)InstUsed[CanBeUsedNums[j][i]];
        }
        for(int j=0;j!=counterClass[i];j++)
        {
            PMass[j] = PMass[j] / PSumm;
        }
        for(int j=0;j!=InstClass[i];j++)
        {
            bool repeatFlag = true;
            while(repeatFlag)
            {
                repeatFlag = false;
                RandomPoint = Random(0,1);
                if(RandomPoint <= PMass[0])
                {
                    Inst[TotalInstCounter] = CanBeUsedNums[0][i];
                    for(int k=PreviousClassEnd;k!=TotalInstCounter;k++)
                    //for(int k=0;k!=TotalInstCounter;k++)
                    {
                        if(Inst[TotalInstCounter] == Inst[k])
                        {
                            repeatFlag = true;
                            break;
                        }
                    }
                }
                else
                {
                    for(int k=1;k!=counterClass[i];k++)
                    {
                        if(RandomPoint > PMass[k-1] && RandomPoint <=PMass[k])
                        {
                            Inst[TotalInstCounter] = CanBeUsedNums[k][i];
                            for(int k=PreviousClassEnd;k!=TotalInstCounter;k++)
                            //for(int k=0;k!=TotalInstCounter;k++)
                            {
                                if(Inst[TotalInstCounter] == Inst[k])
                                {
                                    repeatFlag = true;
                                    break;
                                }
                            }
                            break;
                        }
                    }
                }
            }
            //cout<<Inst[TotalInstCounter]<<" "<<Samp.Classes[Inst[TotalInstCounter]]<<endl;
            TotalInstCounter ++;
        }
        PreviousClassEnd = TotalInstCounter;
    }
}
void SetAvailableInst(sample &Samp,Params &PRS, int &NumInst, int* &InstUsed, int* &Inst, int* &InstClass,
                      int &FoldOnTest, int** &CanBeUsedNums, int* &counterClass, int& CVLearnSize)
{
    //cout<<endl;
    int TotalNumInst = 0;
    int BigestClassNum = 0;
    for(int i=0;i!=Samp.NClasses;i++)
    {
        InstClass[i] = (int)round((float)NumInst * (float)counterClass[i] / (float)CVLearnSize);
        TotalNumInst += InstClass[i];
    }
    while(TotalNumInst > NumInst)
    {
        BigestClassNum = 0;
        for(int i=0;i!=Samp.NClasses;i++)
            if(InstClass[i] >= BigestClassNum)
                BigestClassNum = i;
        InstClass[BigestClassNum] --;
        TotalNumInst --;
    }
    while(TotalNumInst < NumInst)
    {
        int NumMin = InstClass[0];
        int MinPos = 0;
        for(int i=0;i!=Samp.NClasses;i++)
            if(NumMin > InstClass[i])
            {
                NumMin = InstClass[i];
                MinPos = i;
            }
        InstClass[MinPos]++;
        TotalNumInst++;
    }
    float RandomPoint;
    float PSumm;
    float* PMass;
    PMass = new float[CVLearnSize];
    int TotalInstCounter = 0;
    int PreviousClassEnd = 0;
    for(int i=0;i!=Samp.NClasses;i++)
    {
        PMass[0] = 1./(float)InstUsed[CanBeUsedNums[0][i]];
        PSumm = PMass[0];
        for(int j=1;j!=counterClass[i];j++)
        {
            PMass[j] = PMass[j-1] + 1./(float)InstUsed[CanBeUsedNums[j][i]];
            PSumm += 1./(float)InstUsed[CanBeUsedNums[j][i]];
        }
        for(int j=0;j!=counterClass[i];j++)
        {
            PMass[j] = PMass[j] / PSumm;
        }
        for(int j=0;j!=InstClass[i];j++)
        {
            bool repeatFlag = true;
            while(repeatFlag)
            {
                repeatFlag = false;
                RandomPoint = Random(0,1);
                if(RandomPoint <= PMass[0])
                {
                    Inst[TotalInstCounter] = CanBeUsedNums[0][i];
                    for(int k=PreviousClassEnd;k!=TotalInstCounter;k++)
                    //for(int k=0;k!=TotalInstCounter;k++)
                    {
                        if(Inst[TotalInstCounter] == Inst[k])
                        {
                            repeatFlag = true;
                            break;
                        }
                    }
                }
                else
                {
                    for(int k=1;k!=counterClass[i];k++)
                    {
                        if(RandomPoint > PMass[k-1] && RandomPoint <=PMass[k])
                        {
                            Inst[TotalInstCounter] = CanBeUsedNums[k][i];
                            for(int k=PreviousClassEnd;k!=TotalInstCounter;k++)
                            //for(int k=0;k!=TotalInstCounter;k++)
                            {
                                if(Inst[TotalInstCounter] == Inst[k])
                                {
                                    repeatFlag = true;
                                    break;
                                }
                            }
                            break;
                        }
                    }
                }
            }
            //cout<<Inst[TotalInstCounter]<<" "<<Samp.Classes[Inst[TotalInstCounter]]<<endl;
            TotalInstCounter ++;
        }
        PreviousClassEnd = TotalInstCounter;
    }
    //qSort(Inst,0,CVLearnSize-1);
    //cout<<endl;
    //for(int i=0;i!=CVLearnSize;i++)
    //{
        //cout<<Inst[i]<<" ";
    //}
    //cout<<endl;
    /*
    for(int i=0;i!=NumInst;i++)
    {
        int MinEl = InstClass[0];
        int MinNum = 0;
        for(int k=1;k!=Samp.NClasses;k++)
        {
            if((float)InstClass[k]/(float)Samp.NClassInst[k] < (float)MinEl/(float)Samp.NClassInst[MinNum])
            {
                MinEl = InstClass[k];
                MinNum = k;
            }

        }
        int RandomPoint;



        do
        {
            RandomPoint = IntRandom(Samp.NClassInst[MinNum]);
            //if(FoldOnTest == 3 && i > 618)
                //cout<<RandomPoint<<" ";
        } while(Samp.CVFoldNum[Samp.ClassPositions[MinNum][RandomPoint]] == FoldOnTest
           || InstUsed[Samp.ClassPositions[MinNum][RandomPoint]] != 0);



        Inst[i] = Samp.ClassPositions[MinNum][RandomPoint];
        InstUsed[Inst[i]] = 0;
        InstClass[MinNum]++;
        //cout<<i<<"\t"<<Inst[i]<<"\t"<<Samp.Classes[Inst[i]]<<endl;
        //if(FoldOnTest == 3 && i > 618)
            //for(int j=0;j!=Samp.Size;j++)
            //{
                //if(InstUsed[j] == 0)
                    //cout<<j<<" "<<InstUsed[j]<<" "<<Samp.CVFoldNum[j]<<endl;
            //}
        //int j=0;
    }
    */
}
void GenerateRB(int Num, Params &PRS, int NumInst, int* Inst, int*** &Popul, float** &RUPD, float*** &MRS, sample &Samp)
{
    int RandomPattern;
    float RandomPoint;
    float SummOR;
    for(int k=0;k!=PRS.MaxNRules>>1;k++)
    {
        RandomPattern = IntRandom(NumInst);
        for(int i=0;i!=Samp.NVars;i++)
        {
            Popul[Num][k][i] = 0;
            if(Random(0,1) > PRS.PDC )
            {
                SummOR=PRS.NPartitions;                               //Depends on NFSets!
                RandomPoint = Random(0,SummOR);
                do
                {
                    Popul[Num][k][i]++;
                    //SummOR-=LORS[RandomPattern][i][Popul[Num][k][i]];
                    SummOR-=MRS[Inst[RandomPattern]][i][Popul[Num][k][i]];
                } while(SummOR > RandomPoint && Popul[Num][k][i] != 14);

            }
            //cout<<Popul[Num][k][i]<<" ";
        }
        //cout<<"\t";
        //cout<<endl;
    }
    for(int k=PRS.MaxNRules>>1;k!=PRS.MaxNRules;k++)
    {
        for(int i=0;i!=Samp.NVars;i++)
        {
            Popul[Num][k][i] = 0;
        }
    }
    for(int i=0;i!=PRS.MaxNRules;i++)
    {
        RUPD[Num][i] = 0;
    }
    //cout<<endl;
}
void CFClassRUPD(int Num,int k, Params &PRS, int NumInst, int* Inst, int*** &Popul, float** &RUPD, float*** &MRS,
                 sample &Samp, float* &ClassConf, float*** &PopulMu, float** &PopulClass,    float** &PopulCF,
                 int* &InstClass)
{
    float CurConf=0;
    float MaxConf=0;
    float tempConf=0;
    int maxnum=0;
    bool DCRule;
    int SummAc=0;
    DCRule=true;
    for(int i=0;i!=Samp.NVars;i++)
        SummAc+=Popul[Num][k][i];
    for(int i=0;i!=Samp.NClasses;i++)
        ClassConf[i]=0;
    if(SummAc != 0)
    {
        for(int i=0;i!=NumInst;i++)
        {
            tempConf=1;
            for(int j=0;j!=Samp.NVars;j++)
            {
                if(Popul[Num][k][j] > 0)
                {
                    DCRule = false;
                    //tempConf *= LORS[i][j][Popul[Num][k][j]];
                    tempConf *= MRS[Inst[i]][j][Popul[Num][k][j]];
                }
            }
            if(DCRule == true)
            {
                tempConf=0;
            }
            //PopulMu[Num][k][i] = tempConf;
            //cout<<PopulMu[Num][k][i]<<" ";
            ClassConf[ Samp.Classes[Inst[i]] ] += tempConf;
            CurConf += tempConf;
        }
    }
    //cout<<endl;
    if(DCRule == true)
        CurConf=1;
    if(CurConf != 0)
    {
        for(int i=0;i!=Samp.NClasses;i++)
        {
            if(PRS.Unbiased_init == 1)
                ClassConf[i] = ClassConf[i] / InstClass[i] * Samp.Size;
        }
        ///////////////////////////////////////probabilistic rule numbering
        //float SumClassConf = 0;
        for(int i=0;i!=Samp.NClasses;i++)
        {
            ClassConf[i] /= CurConf;
            //SumClassConf += ClassConf[i];
        }
        maxnum = 0;
        MaxConf = ClassConf[0];
        for(int i=1;i!=Samp.NClasses;i++)
        {
            //ClassConf[i] /= CurConf;
            if(ClassConf[i] > MaxConf)
            {
                MaxConf = ClassConf[i];
                maxnum = i;
            }
        }
        //float totalConf = ClassConf[0];
        //float RandomPoint = Random(0,SumClassConf);
        //for(int i=1;i!=Samp.NClasses;i++)
        //{
            //if(RandomPoint > totalConf && RandomPoint <= totalConf + ClassConf[i])
            //{
                //MaxConf = ClassConf[i];
                //maxnum = i;
                //break;
            //}
            //totalConf += ClassConf[i];
        //}
        ////////////////////////////////////////////////////////////////
        if(PRS.Unbiased_init == 1)
            MaxConf = MaxConf * InstClass[maxnum] / Samp.Size;
        //cout<<MaxConf<<endl;
        if(MaxConf < 0.5)
        {
            PopulClass[Num][k] = -1;
            PopulCF[Num][k] = 0;
        }
        else
        {
            PopulClass[Num][k] = maxnum;
            PopulCF[Num][k] = 2.*MaxConf - 1.;
            PRS.HasGoodRules = true;
        }
    }
    else
    {
        PopulClass[Num][k] = -1;
        PopulCF[Num][k] = 0;
    }
    RUPD[Num][k] = 1;
    //cout<<PopulClass[Num][k]<<" ";
    //cout<<PopulCF[Num][k]<<endl;
}
void CFClassRB(int Num, sample &Samp, Params &PRS, float** &RUPD, int*** &Popul, float* &ClassConf, float*** &PopulMu,
               float*** &MRS, float** &PopulClass, float** &PopulCF, int NumInst, int* Inst, int* &InstClass)
{
    PRS.HasGoodRules=false;
    //cout<<endl;
    float CurConf=0;
    float MaxConf=0;
    float tempConf=0;
    int maxnum=0;
    bool DCRule=false;
    int SummAc;
    for(int k=0;k!=PRS.MaxNRules;k++)
    {
        if(RUPD[Num][k] == 0)
        {
            CurConf=0;
            MaxConf=0;
            maxnum=0;
            tempConf=0;
            SummAc=0;
            DCRule=true;
            for(int i=0;i!=Samp.NVars;i++)
                SummAc+=Popul[Num][k][i];
            for(int i=0;i!=Samp.NClasses;i++)
                ClassConf[i]=0;
            if(SummAc != 0)
            {
                for(int i=0;i!=NumInst;i++)
                {
                    tempConf=1;
                    for(int j=0;j!=Samp.NVars;j++)
                    {
                        if(Popul[Num][k][j] > 0)
                        {
                            DCRule = false;
                            //tempConf *= LORS[i][j][Popul[Num][k][j]];
                            tempConf *= MRS[Inst[i]][j][Popul[Num][k][j]];
                        }
                    }
                    if(DCRule == true)
                    {
                        tempConf=0;
                    }
                    //PopulMu[Num][k][i] = tempConf;
                    //cout<<PopulMu[Num][k][i]<<" ";
                    ClassConf[ Samp.Classes[Inst[i]] ] += tempConf;
                    CurConf += tempConf;
                }
            }
            //cout<<endl;
            if(DCRule == true)
                CurConf=1;
            if(CurConf != 0)
            {
                for(int i=0;i!=Samp.NClasses;i++)
                {
                    if(PRS.Unbiased_init == 1)
                        ClassConf[i] = ClassConf[i] / InstClass[i] * Samp.Size;
                }
                ///////////////////////////////////////probabilistic rule numbering
                //float SumClassConf = 0;
                for(int i=0;i!=Samp.NClasses;i++)
                {
                    ClassConf[i] /= CurConf;
                    //SumClassConf += ClassConf[i];
                }
                maxnum = 0;
                MaxConf = ClassConf[0];
                for(int i=1;i!=Samp.NClasses;i++)
                {
                    //ClassConf[i] /= CurConf;
                    if(ClassConf[i] > MaxConf)
                    {
                        MaxConf = ClassConf[i];
                        maxnum = i;
                    }
                }
                //float totalConf = ClassConf[0];
                //float RandomPoint = Random(0,SumClassConf);
                //for(int i=1;i!=Samp.NClasses;i++)
                //{
                    //if(RandomPoint > totalConf && RandomPoint <= totalConf + ClassConf[i])
                    //{
                        //MaxConf = ClassConf[i];
                        //maxnum = i;
                        //break;
                    //}
                    //totalConf += ClassConf[i];
                //}
                ////////////////////////////////////////////////////////////////
                if(PRS.Unbiased_init == 1)
                    MaxConf = MaxConf * InstClass[maxnum] / Samp.Size;
                //cout<<MaxConf<<endl;
                if(MaxConf < 0.5)
                {
                    PopulClass[Num][k] = -1;
                    PopulCF[Num][k] = 0;
                }
                else
                {
                    PopulClass[Num][k] = maxnum;
                    PopulCF[Num][k] = 2.*MaxConf - 1.;
                    PRS.HasGoodRules = true;
                }
            }
            else
            {
                PopulClass[Num][k] = -1;
                PopulCF[Num][k] = 0;
            }
            RUPD[Num][k] = 1;
            //cout<<PopulClass[Num][k]<<" ";
            //cout<<PopulCF[Num][k]<<endl;
        }
    }
    //cout<<endl;
}

void FitCalc(int Num, sample &Samp, float** &ErrorClass, float** &ConfMatrix, int NumInst, Params &PRS,
             float** &PopulCF, float*** &PopulMu, float** &PopulClass,int* &CurNumRules,
             float* &CurRLength, float* &ErrMass, float* &AveMass, float* &FitMass, float* &NumClassInst,
             int*** &Popul, int* &Inst, int ChangeInst, int* &InstUsed, int* &InstClass,float*** &MRS)
{
    if(ChangeInst == 1)
        for(int i=0;i!=NumInst;i++)
        {
            InstUsed[Inst[i]] += 1;
            //cout<<Inst[i]<<"\t"<<InstUsed[Inst[i]]<<endl;
        }
    float maxCFMu;
    int maxnum;
    float Error=0;
    for(int i=0;i!=Samp.NClasses;i++)
    {
        ErrorClass[Num][i] = 0;
    }
    for(int i=0;i!=Samp.NClasses;i++)
        for(int j=0;j!=Samp.NClasses+1;j++)
        {
            ConfMatrix[i][j] = 0;
        }
    /////////////////////////////

    float CurConf=0;
    //float MaxConf=0;
    float tempConf=0;
    bool DCRule=false;
    int SummAc;

    for(int k=0;k!=PRS.MaxNRules;k++)
    {
        CurConf=0;
        //MaxConf=0;
        tempConf=0;
        SummAc=0;
        DCRule=true;
        for(int i=0;i!=Samp.NVars;i++)
            SummAc+=Popul[Num][k][i];
        if(SummAc != 0)
        {
            for(int i=0;i!=NumInst;i++)
            {
                tempConf=1;
                for(int j=0;j!=Samp.NVars;j++)
                {
                    if(Popul[Num][k][j] > 0)
                    {
                        DCRule = false;
                        //tempConf *= LORS[i][j][Popul[Num][k][j]];
                        tempConf *= MRS[Inst[i]][j][Popul[Num][k][j]];
                    }
                }
                if(DCRule == true)
                {
                    tempConf=0;
                }
                PopulMu[0][k][i] = tempConf;
                //cout<<PopulMu[Num][k][i]<<" ";
                CurConf += tempConf;
            }
        }
    }

    /////////////////////////////
    for(int i=0;i!=NumInst;i++)
    {
        maxCFMu=0;
        maxnum=-1;
        for(int j=0;j!=PRS.MaxNRules;j++)
        {
            if(PopulCF[Num][j]*PopulMu[0][j][i] > maxCFMu && PopulCF[Num][j]*PopulMu[0][j][i] > 0)
            {
                maxCFMu = PopulCF[Num][j]*PopulMu[0][j][i];
                maxnum = j;
                continue;
            }
            if(PopulCF[Num][j]*PopulMu[0][j][i] == maxCFMu && PopulCF[Num][j]*PopulMu[0][j][i] > 0)
            {
                maxnum = -1;
                break;
            }
        }
        if(maxnum != -1)
        {
            if(PopulClass[Num][maxnum] != Samp.Classes[Inst[i]])
            {
                //cout<<Inst[i]<<" ";
                Error++;
                ErrorClass[Num][ Samp.Classes[Inst[i]] ] ++;
                if(ChangeInst == 1)
                {
                    InstUsed[Inst[i]] = 1;
                    //InstUsed[Inst[i]] = int((float)InstUsed[Inst[i]]/2.);
                    //cout<<Inst[i]<<" ";
                }
            }
            ConfMatrix[ Samp.Classes[Inst[i]] ][ (int)PopulClass[Num][maxnum] ] ++;
        }
        else
        {
            //cout<<Inst[i]<<" ";
            Error++;
            ErrorClass[Num][ Samp.Classes[Inst[i]] ] ++;
            ConfMatrix[ Samp.Classes[Inst[i]] ][ Samp.NClasses ] ++;
            if(ChangeInst == 1)
            {
                InstUsed[Inst[i]] = 1;
                //InstUsed[Inst[i]] = int((float)InstUsed[Inst[i]]/2.);
                //cout<<Inst[i]<<" ";
            }
        }
    }
    //cout<<endl;
    CurNumRules[Num]=0;
    CurRLength[Num] =0;
    for(int i=0;i!=PRS.MaxNRules;i++)
    {
        if(PopulClass[Num][i] >= 0)
        {
            CurNumRules[Num]++;
            for(int j=0;j!=Samp.NVars;j++)
            {
                if(Popul[Num][i][j] > 0)
                {
                    CurRLength[Num]++;
                }
                //cout<<Popul[Num][i][j]<<" ";
            }
            //cout<<"-> "<<PopulClass[Num][i]<<" "<<(PopulCF[Num][i])<<endl;
        }
    }
    if(CurNumRules[Num] != 0)
        CurRLength[Num] /= CurNumRules[Num];
    else
        CurRLength[Num] = -1;
    ErrMass[Num] = Error;
    AveMass[Num]=0;
    if(ChangeInst == 1)
        for(int i=0;i!=NumInst;i++)
        {
            //cout<<Inst[i]<<"\t"<<InstUsed[Inst[i]]<<endl;
        }
    for(int i=0;i!=Samp.NClasses;i++)
    {
        ErrorClass[Num][i] /= InstClass[i];
        AveMass[Num] += ErrorClass[Num][i]/(float)Samp.NClasses;
    }
    if(PRS.CritType == 0)
        FitMass[Num] = Error/NumInst*10000 + 1*CurNumRules[Num] + 1*CurNumRules[Num]*CurRLength[Num];
    else
        FitMass[Num] = AveMass[Num]*10000 + 1*CurNumRules[Num] + 1*CurNumRules[Num]*CurRLength[Num];
    //FitMass[Num] = 1./(1.+FitMass[Num]);

    if(CurNumRules[Num] == 1)
        FitMass[Num] *= 2;
    if(CurRLength[Num] == 1)
        FitMass[Num] *= 2;
}

void Init(sample &Samp,Params &PRS, float* &RuleClassified,
    float* &RuleP,    float* &tempRuleC,    float* &BadRNums,    int* &MisClassifiedNums,    int* &RandomNumbers,
    int* &CurNumRules,    float* &CurRLength,    int** &GoodRNums,    float** &ChosenRules,    float* &FitMass,
    float* &FitMassCopy,    float* &PropSetSel,    float* &ClassConf,    float** &PopulClass,    float** &PopulCF,
    float*** &PopulMu,    int** &BestInd,    float* &BestClasses,    float* &BestCF,    float** &BestMu,
    float** &RUPD, float* &ErrMass,    int* &SelectedNums1,    int* &SelectedNums2,    float* &NumMass,
    float* &PMass,    float* &NumMass2, float* &RankMass,    float* &RMass,    float* &SavedFit,
    float* &NumClassInst,    float* &NumClassInstWrong, float* &NumClassInstTest,    float* &NumClassInstTestWrong,
    float** &ConfMatrix,    float** &ConfMatrixTest, float** &ErrorClass,    float* &AveMass,
    float* &TestErrorClass,    int*** &Popul,    float* &psel,    float* &pcrs,  float* &pmut,
    float* &pheu,    float* &usedsel,    float* &usedcrs,    float* &usedmut,    float* &usedheu,
    float* &sumfitsel,    float* &sumfitcrs,    float* &sumfitmut,    float* &sumfitheu,    float* &numselused,
    float* &numcrsused,    float* &nummutused,    float* &numheuused,    float* &sucsel,    float* &succrs,
    float* &sucmut,    float* &sucheu, int &NumInst, int* &InstUsed, int* &Inst, int* &InstClass, int &FoldOnTest,
    float*** &MRS, int CVLearnSize)
{
    TestErrorClass = new float[Samp.NClasses];
    AveMass=new float[PRS.NumInds*2];
    ErrorClass = new float*[PRS.NumInds*2];
    for(int i=0;i!=PRS.NumInds*2;i++)
    {
        ErrorClass[i] = new float[Samp.NClasses];
    }
    NumClassInst=new float[Samp.NClasses];
    for(int i=0;i!=Samp.NClasses;i++)
    {
        NumClassInst[i]=InstClass[i];
    }
    NumClassInstWrong=new float[Samp.NClasses];
    for(int i=0;i!=Samp.NClasses;i++)
    {
        NumClassInstWrong[i]=0;
    }
    NumClassInstTest=new float[Samp.NClasses];
    for(int i=0;i!=Samp.NClasses;i++)
    {
        NumClassInstTest[i]=0;
    }
    //for(int i=0;i!=S.TestSize;i++)
    {
        //NumClassInstTest[ (int)S.GetClass(i) ] ++;
    }
    NumClassInstTestWrong=new float[Samp.NClasses];
    for(int i=0;i!=Samp.NClasses;i++)
    {
        NumClassInstTestWrong[i]=0;
    }
    ConfMatrix = new float*[Samp.NClasses];
    for(int i=0;i!=Samp.NClasses;i++)
    {
        ConfMatrix[i] = new float[Samp.NClasses+1];
        for(int j=0;j!=Samp.NClasses+1;j++)
        {
            ConfMatrix[i][j] = 0;
        }
    }
    ConfMatrixTest = new float*[Samp.NClasses];
    for(int i=0;i!=Samp.NClasses;i++)
    {
        ConfMatrixTest[i] = new float[Samp.NClasses+1];
        for(int j=0;j!=Samp.NClasses+1;j++)
        {
            ConfMatrixTest[i][j] = 0;
        }
    }
    psel=new float[PRS.NumSel];
	pcrs=new float[PRS.NumCrs];
	pmut=new float[PRS.NumMut];
	pheu=new float[PRS.NumHeu];
	usedsel=new float[PRS.NumSel];
	usedcrs=new float[PRS.NumCrs];
	usedmut=new float[PRS.NumMut];
	usedheu=new float[PRS.NumHeu];
	sumfitsel=new float[PRS.NumSel];
	sumfitcrs=new float[PRS.NumCrs];
	sumfitmut=new float[PRS.NumMut];
	sumfitheu=new float[PRS.NumHeu];
	numselused=new float[PRS.NumSel];
	numcrsused=new float[PRS.NumCrs];
	nummutused=new float[PRS.NumMut];
	numheuused=new float[PRS.NumHeu];
	sucsel=new float[PRS.NumSel];
	succrs=new float[PRS.NumMut];
	sucmut=new float[PRS.NumMut];
	sucheu=new float[PRS.NumHeu];
    for(int i=0;i!=PRS.NumSel;i++)
    {
        psel[i]=1./PRS.NumSel;
        usedsel[i] = 1;
        sucsel[i] = 1;
    }
    for(int i=0;i!=PRS.NumCrs;i++)
    {
        pcrs[i]=1./PRS.NumCrs;
        usedcrs[i] = 1;
        succrs[i] = 1;
    }
    for(int i=0;i!=PRS.NumMut;i++)
    {
        pmut[i]=1./PRS.NumMut;
        usedmut[i] = 1;
        sucmut[i] = 1;
    }
    for(int i=0;i!=PRS.NumHeu;i++)
    {
        pheu[i]=1./PRS.NumHeu;
        usedheu[i] = 1;
        sucheu[i] = 1;
    }
    PRS.pselall=0.2/PRS.NumSel;
    PRS.pcrsall=0.2/PRS.NumCrs;
    PRS.pmutall=0.2/PRS.NumMut;
    PRS.pheuall=0.2/PRS.NumHeu;
    PRS.HasBestFit=false;
    SavedFit = new float[PRS.NumInds*2];
    RMass = new float[PRS.NumInds*2];
    RankMass = new float[PRS.NumInds*2];
    NumMass2 = new float[PRS.NumInds*2];
    PMass = new float[PRS.NumInds*2];
    RUPD = new float*[PRS.NumInds*2];
    for(int i=0;i!=PRS.NumInds*2;i++)
    {
        RUPD[i] = new float[PRS.MaxNRules];
    }
    FitMassCopy = new float[PRS.NumInds*2];
    NumMass = new float[PRS.NumInds*2];
    SelectedNums1 = new int[PRS.NumInds];
    SelectedNums2 = new int[PRS.NumInds];
    MisClassifiedNums = new int[CVLearnSize];
    ErrMass = new float[PRS.NumInds*2];
    RuleClassified = new float[PRS.MaxNRules];
    RuleP = new float[PRS.MaxNRules*2];
    tempRuleC = new float[PRS.MaxNRules];
    BadRNums = new float[PRS.MaxNRules];
    CurRLength = new float[PRS.NumInds*2];
    BestInd = new int*[PRS.MaxNRules];
    for(int i=0;i!=PRS.MaxNRules;i++)
    {
        BestInd[i] = new int[Samp.NVars];
    }
    BestClasses = new float[PRS.MaxNRules];
    BestCF = new float[PRS.MaxNRules];
    BestMu = new float*[PRS.MaxNRules];
    for(int i=0;i!=PRS.MaxNRules;i++)
    {
        BestMu[i] = new float[CVLearnSize];
    }
    GoodRNums = new int*[PRS.NumInds*2];
    for(int i=0;i!=PRS.NumInds*2;i++)
    {
        GoodRNums[i] = new int[PRS.MaxNRules];
    }
    ChosenRules = new float*[PRS.MaxNRules*2];
    for(int i=0;i!=PRS.MaxNRules*2;i++)
    {
        ChosenRules[i] = new float[PRS.NVars+1];
    }
    CurNumRules = new int[PRS.NumInds*2];
    RandomNumbers = new int[PRS.NumInds];
    FitMass = new float[PRS.NumInds*2];
    PopulCF = new float*[PRS.NumInds*2];
    PopulClass = new float*[PRS.NumInds*2];
    PopulMu = new float**[1];

    PopulMu[0] = new float*[PRS.MaxNRules];
    for(int j=0;j!=PRS.MaxNRules;j++)
    {
        PopulMu[0][j] = new float[NumInst];
    }
    for(int i=0;i!=PRS.NumInds*2;i++)
    {
        PopulCF[i] = new float[PRS.MaxNRules];
        PopulClass[i] = new float[PRS.MaxNRules];
    }
    ClassConf = new float[Samp.NClasses];
    PropSetSel = new float[PRS.NFSets];
    Popul = new int**[PRS.NumInds*2];
    for(int i=0;i!=PRS.NumInds*2;i++)
    {
        Popul[i] = new int*[PRS.MaxNRules];
        for(int j=0;j!=PRS.MaxNRules;j++)
        {
            Popul[i][j] = new int[Samp.NVars];
            for(int k=0;k!=Samp.NVars;k++)
            {
                Popul[i][j][k] = -1;
            }
        }
        if(i<PRS.NumInds)
        {
            PRS.HasGoodRules=false;
            do
            {
                GenerateRB(i,PRS,NumInst,Inst,Popul,RUPD,MRS,Samp);
                CFClassRB(i,Samp,PRS,RUPD,Popul,ClassConf,PopulMu,MRS,PopulClass,PopulCF,NumInst,Inst,InstClass);

            } while(PRS.HasGoodRules == false); //HasGoodRules == false &&
            FitCalc(i,Samp,ErrorClass,ConfMatrix,NumInst,PRS,PopulCF,PopulMu,PopulClass,CurNumRules,CurRLength,ErrMass,
                    AveMass,FitMass,NumClassInst,Popul,Inst,0,InstUsed,InstClass,MRS);
            //cout<<ErrMass[i]<<endl;
        }
    }
}
void Clean(sample &Samp,Params &PRS, float* &RuleClassified,
    float* &RuleP,    float* &tempRuleC,    float* &BadRNums,    int* &MisClassifiedNums,    int* &RandomNumbers,
    int* &CurNumRules,    float* &CurRLength,    int** &GoodRNums,    float** &ChosenRules,    float* &FitMass,
    float* &FitMassCopy,    float* &PropSetSel,    float* &ClassConf,    float** &PopulClass,    float** &PopulCF,
    float*** &PopulMu,    int** &BestInd,    float* &BestClasses,    float* &BestCF,    float** &BestMu,
    float** &RUPD, float* &ErrMass,    int* &SelectedNums1,    int* &SelectedNums2,    float* &NumMass,
    float* &PMass,    float* &NumMass2, float* &RankMass,    float* &RMass,    float* &SavedFit,
    float* &NumClassInst,    float* &NumClassInstWrong, float* &NumClassInstTest,    float* &NumClassInstTestWrong,
    float** &ConfMatrix,    float** &ConfMatrixTest, float** &ErrorClass,    float* &AveMass,
    float* &TestErrorClass,    int*** &Popul,    float* &psel,    float* &pcrs,  float* &pmut,
    float* &pheu,    float* &usedsel,    float* &usedcrs,    float* &usedmut,    float* &usedheu,
    float* &sumfitsel,    float* &sumfitcrs,    float* &sumfitmut,    float* &sumfitheu,    float* &numselused,
    float* &numcrsused,    float* &nummutused,    float* &numheuused,    float* &sucsel,    float* &succrs,
    float* &sucmut,    float* &sucheu, int &NumInst, int* &InstUsed, int* &Inst, int* &InstClass, int &FoldOnTest,
    float*** &MRS)
{
    //cout<<"Cleaning..."<<endl;
    for(int i=0;i!=Samp.NClasses;i++)
    {
        delete ConfMatrix[i];
        delete ConfMatrixTest[i];
    }
    delete ConfMatrix;
    delete ConfMatrixTest;
    delete TestErrorClass;
    delete NumClassInstTest;
    delete NumClassInstTestWrong;
    delete AveMass;
    for(int i=0;i!=PRS.NumInds*2;i++)
    {
        delete ErrorClass[i];
    }
    delete ErrorClass;
    delete NumClassInst;
    delete SavedFit;
    delete RMass;
    delete RankMass;
    delete NumMass2;
    delete PMass;
    for(int i=0;i!=PRS.NumInds*2;i++)
    {
        delete RUPD[i];
    }
    delete RUPD;
    delete FitMassCopy;
    delete NumMass;
    delete SelectedNums1;
    delete SelectedNums2;
    delete MisClassifiedNums;
    delete ErrMass;
    delete RuleClassified;
    delete RuleP;
    delete tempRuleC;
    delete BadRNums;
    delete CurRLength;
    for(int i=0;i!=PRS.MaxNRules;i++)
    {
        delete BestInd[i];
        delete BestMu[i];
    }
    delete BestInd;
    delete BestMu;
    delete BestCF;
    delete BestClasses;
    for(int i=0;i!=PRS.NumInds*2;i++)
    {
        delete GoodRNums[i];
    }
    delete GoodRNums;
    for(int i=0;i!=PRS.MaxNRules*2;i++)
    {
        delete ChosenRules[i];
    }
    delete ChosenRules;
    delete CurNumRules;
    delete RandomNumbers;
    delete FitMass;
    for(int j=0;j!=PRS.MaxNRules;j++)
    {
        delete PopulMu[0][j];
    }
    delete PopulMu[0];
    for(int i=0;i!=PRS.NumInds*2;i++)
    {
        delete PopulCF[i];
        delete PopulClass[i];
    }
    delete PopulCF;
    delete PopulClass;
    delete PopulMu;
    delete ClassConf;
    delete PropSetSel;
    for(int i=0;i!=PRS.NumInds*2;i++)
    {
        for(int j=0;j!=PRS.MaxNRules;j++)
        {
            delete Popul[i][j];
        }
        delete Popul[i];
    }
    delete Popul;
    //cout<<"Done cleaning!"<<endl;
}
void PropSelection(int Num, Params &PRS, float* &PMass, int* &SelectedNums1, int* &SelectedNums2)
{
    int Taken1=0;
    int Taken2=0;
    float RN1=Random(0,1);
    float RN2=Random(0,1);
    if(RN1 >= PMass[0])
    {
        for(int j=1;j!=PRS.NumInds;j++)
        {
            if( (RN1>=PMass[j-1]) && (RN1<=PMass[j]) )
                Taken1=j;
        }
    }
    if(RN2 >= PMass[0])
    {
        for(int j=1;j!=PRS.NumInds;j++)
        {
            if( (RN2>=PMass[j-1]) && (RN2<=PMass[j]) )
                Taken2=j;
        }
    }
    Taken1=1;
    Taken2=1;
    SelectedNums1[Num] = Taken1;
    SelectedNums2[Num] = Taken2;
}
void RankSelection(int Num, Params &PRS, float* &PMass, float* RMass, int* &SelectedNums1, int* &SelectedNums2)
{
    int Taken1=0;
    int Taken2=0;
    float RN1=Random(0,1);
    float RN2=Random(0,1);

    if(RN1 >= RMass[0])
    {
        for(int j=1;j!=PRS.NumInds;j++)
        {
            if( (RN1>=RMass[j-1]) && (RN1<=RMass[j]) )
                Taken1=j;
        }
    }
    if(RN2 >= RMass[0])
    {
        for(int j=1;j!=PRS.NumInds;j++)
        {
            if( (RN2>=RMass[j-1]) && (RN2<=RMass[j]) )
                Taken2=j;
        }
    }
    SelectedNums1[Num] = Taken1;
    SelectedNums2[Num] = Taken2;
}
void TourSelection(int Num, Params &PRS, float* &FitMass, int* &SelectedNums1, int* &SelectedNums2, int* &RandomNumbers)
{
    float maxfit=0;
    int Taken1;
    int Taken2;
    Taken1=0;
    for(int j=0;j!=PRS.TourSize;j++)
    {
        RandomNumbers[j]=IntRandom(PRS.NumInds);
        if(j==0)
            maxfit=FitMass[RandomNumbers[0]];
        else
        {
            if(FitMass[RandomNumbers[j]] < maxfit)
            {
                Taken1 = RandomNumbers[j];
                maxfit = FitMass[RandomNumbers[j]];
            }
        }
    }
    Taken2=0;
    for(int j=0;j!=PRS.TourSize;j++)
    {
        RandomNumbers[j]=IntRandom(PRS.NumInds);
        if(j==0)
            maxfit=FitMass[RandomNumbers[0]];
        else
        {
            if(FitMass[RandomNumbers[j]] < maxfit)
            {
                Taken2 = RandomNumbers[j];
                maxfit = FitMass[RandomNumbers[j]];
            }
        }
    }
    SelectedNums1[Num] = Taken1;
    SelectedNums2[Num] = Taken2;
    /*for(int i=0;i!=MaxNRules;i++)
    {
        for(int j=0;j!=NVars;j++)
        {
            Popul[NumInds+Num][i][j]=Popul[Taken1][i][j];
            Popul[NumInds*2+Num][i][j]=Popul[Taken2][i][j];
        }
        PopulClass[NumInds+Num][i] = PopulClass[Taken1][i];
        PopulClass[NumInds*2+Num][i] = PopulClass[Taken2][i];
        PopulCF[NumInds+Num][i] = PopulCF[Taken1][i];
        PopulCF[NumInds*2+Num][i] = PopulCF[Taken2][i];
        for(int j=0;j!=LearnSize;j++)
        {
            PopulMu[NumInds+Num][i][j] = PopulMu[Taken1][i][j];
            PopulMu[NumInds*2+Num][i][j] = PopulMu[Taken2][i][j];
        }
    }
    //cout<<Taken1<<" "<<Taken2<<endl;
    //ofstream foutrules("Rulestest.txt",ios::app);
    for(int i=0;i!= MaxNRules;i++)
    {
        if(PopulClass[NumInds+Num][i] != -1)
        {
            //foutrules<<i<<" ";
            for(int j=0;j!=NVars;j++)
            {
                //foutrules<<Popul[NumInds+Num][i][j]<<" ";
            }
            //foutrules<<"-> "<<PopulClass[NumInds+Num][i]<<endl;
        }
    }
    //foutrules<<endl;
    for(int i=0;i!= MaxNRules;i++)
    {
        if(PopulClass[NumInds*2+Num][i] != -1)
        {
            //foutrules<<i<<" ";
            for(int j=0;j!=NVars;j++)
            {
                //foutrules<<Popul[NumInds*2+Num][i][j]<<" ";
            }
            //foutrules<<"-> "<<PopulClass[NumInds*2+Num][i]<<endl;
        }
    }
    //foutrules<<endl;*/
}

void PittsCrossover(int Num, Params &PRS, int* &SelectedNums1, int* &SelectedNums2, int** &GoodRNums,
                    float* &FitMass, float** &ChosenRules,int* &CurNumRules,float** &RUPD,
                    int*** &Popul,float*** &PopulMu,float** &PopulClass,    float** &PopulCF, int NumInst)
{
    int CurNRules = 0;
    for(int i=0;i!=PRS.MaxNRules;i++)
    {
        GoodRNums[SelectedNums1[Num]][i] = 0;
        if(PopulClass[SelectedNums1[Num]][i] >= 0)
        {
            GoodRNums[SelectedNums1[Num]][CurNRules] = i;
            CurNRules++;
        }
    }
    CurNumRules[SelectedNums1[Num]] = CurNRules;
    CurNRules = 0;
    for(int i=0;i!=PRS.MaxNRules;i++)
    {
        GoodRNums[SelectedNums2[Num]][i] = 0;
        if(PopulClass[SelectedNums2[Num]][i] >= 0)
        {
            GoodRNums[SelectedNums2[Num]][CurNRules] = i;
            CurNRules++;
        }
    }
    CurNumRules[SelectedNums2[Num]] = CurNRules;
    if(Random(0,1) < 0.9)
    {
        for(int i=0;i!=PRS.MaxNRules;i++)
        {
            //cout<<RuleP[i]<<"\t";
        }
        /////////
        int TotalNRules = CurNumRules[SelectedNums1[Num]] + CurNumRules[SelectedNums2[Num]];
        int NewTotalNRules = 0;
        int CurRandom = 0;
        for(int i=0;i!=CurNumRules[SelectedNums1[Num]];i++)
        {
            ChosenRules[CurRandom][PRS.NVars] = 0;
            CurRandom++;
        }
        for(int i=0;i!=CurNumRules[SelectedNums2[Num]];i++)
        {
            ChosenRules[CurRandom][PRS.NVars] = 0;
            CurRandom++;
        }
        CurRandom=0;
        NewTotalNRules = IntRandom(TotalNRules-1)+1;
        if(NewTotalNRules > PRS.MaxNRules)
            NewTotalNRules = PRS.MaxNRules;
        for(int i=0;i!=NewTotalNRules;i++)
        {
            do
            {
                if(PRS.MichType == 1)
                {
                    if(Random(0,1) < FitMass[SelectedNums1[Num]]/(FitMass[SelectedNums1[Num]]+FitMass[SelectedNums2[Num]]))
                        CurRandom = IntRandom(CurNumRules[SelectedNums1[Num]]);
                    else
                        CurRandom = IntRandom(CurNumRules[SelectedNums2[Num]])+CurNumRules[SelectedNums1[Num]];
                }
                if(PRS.MichType == 2)
                {
                    if(Random(0,1) < 0.5)
                        CurRandom = IntRandom(CurNumRules[SelectedNums1[Num]]);
                    else
                        CurRandom = IntRandom(CurNumRules[SelectedNums2[Num]])+CurNumRules[SelectedNums1[Num]];
                }
                else
                {
                    CurRandom = IntRandom(TotalNRules);
                }

            } while(ChosenRules[CurRandom][PRS.NVars] != 0);
            ChosenRules[CurRandom][PRS.NVars] = 1;

            //cout<<CurRandom<<" ";
            if(CurRandom >= CurNumRules[SelectedNums1[Num]])
            {
                CurRandom -= CurNumRules[SelectedNums1[Num]];
                for(int j=0;j!=PRS.NVars;j++)
                {
                    Popul[Num+PRS.NumInds][i][j] = Popul[SelectedNums2[Num]][ GoodRNums[SelectedNums2[Num]][CurRandom] ][j];
                }
                PopulClass[Num+PRS.NumInds][i] = PopulClass[SelectedNums2[Num]][GoodRNums[SelectedNums2[Num]][CurRandom]];
                PopulCF[Num+PRS.NumInds][i] = PopulCF[SelectedNums2[Num]][GoodRNums[SelectedNums2[Num]][CurRandom]];
                for(int j=0;j!=1;j++)
                {
                    //PopulMu[Num+PRS.NumInds][i][j] = PopulMu[SelectedNums2[Num]][GoodRNums[SelectedNums2[Num]][CurRandom]][j];
                }
                RUPD[Num+PRS.NumInds][i] = 1;
                //if(PopulClass[Num][i] == -1)
                  //  continue;
                //CurRandom += CurNumRules[NumInds+Num];
            }
            else
            {
                for(int j=0;j!=PRS.NVars;j++)
                {
                    Popul[Num+PRS.NumInds][i][j] = Popul[SelectedNums1[Num]][ GoodRNums[SelectedNums1[Num]][CurRandom] ][j];
                }
                PopulClass[Num+PRS.NumInds][i] = PopulClass[SelectedNums1[Num]][GoodRNums[SelectedNums1[Num]][CurRandom]];
                PopulCF[Num+PRS.NumInds][i] = PopulCF[SelectedNums1[Num]][GoodRNums[SelectedNums1[Num]][CurRandom]];
                for(int j=0;j!=1;j++)
                {
                    //PopulMu[Num+PRS.NumInds][i][j] = PopulMu[SelectedNums1[Num]][GoodRNums[SelectedNums1[Num]][CurRandom]][j];
                }
                RUPD[Num+PRS.NumInds][i] = 1;
                //if(PopulClass[Num][i] == -1)
                  //  continue;
            }
            //cout<<endl;
        }
        CurNumRules[Num+PRS.NumInds] = NewTotalNRules;
        for(int i=NewTotalNRules;i!=PRS.MaxNRules;i++)
        {
            for(int j=0;j!=PRS.NVars;j++)
            {
                Popul[Num+PRS.NumInds][i][j] = 0;
            }
            PopulClass[Num+PRS.NumInds][i] = -1;
            PopulCF[Num+PRS.NumInds][i] = 0;
            for(int j=0;j!=NumInst;j++)
            {
                //PopulMu[Num+PRS.NumInds][i][j] = 0;
            }
            RUPD[Num+PRS.NumInds][i] = 1;
        }
    }
    else
    {
        for(int i=0;i!=CurNumRules[SelectedNums1[Num]];i++)
        {
            for(int j=0;j!=PRS.NVars;j++)
            {
                Popul[Num+PRS.NumInds][i][j] = Popul[SelectedNums1[Num]][ GoodRNums[SelectedNums1[Num]][i] ][j];
            }
            PopulClass[Num+PRS.NumInds][i] = PopulClass[SelectedNums1[Num]][GoodRNums[SelectedNums1[Num]][i]];
            PopulCF[Num+PRS.NumInds][i] = PopulCF[SelectedNums1[Num]][GoodRNums[SelectedNums1[Num]][i]];
            for(int j=0;j!=NumInst;j++)
            {
                //PopulMu[Num+PRS.NumInds][i][j] = PopulMu[SelectedNums1[Num]][GoodRNums[SelectedNums1[Num]][i]][j];
            }
            RUPD[Num+PRS.NumInds][i] = 1;
        }
        for(int i=CurNumRules[SelectedNums1[Num]];i!=PRS.MaxNRules;i++)
        {
            for(int j=0;j!=PRS.NVars;j++)
            {
                Popul[Num+PRS.NumInds][i][j] = 0;
            }
            PopulClass[Num+PRS.NumInds][i] = -1;
            PopulCF[Num+PRS.NumInds][i] = 0;
            for(int j=0;j!=NumInst;j++)
            {
                //PopulMu[Num+PRS.NumInds][i][j] = 0;
            }
            RUPD[Num+PRS.NumInds][i] = 1;
        }
    }
    for(int i=0;i!=PRS.MaxNRules;i++)
    {
        for(int j=0;j!=PRS.NVars;j++)
        {
            //cout<<Popul[Num][i][j]<<" ";
        }
        //cout<<PopulClass[Num][i]<<" "<<PopulCF[Num][i];
        //cout<<endl;
    }
    //cout<<endl;
    //_getch();
}
void PittsMutation(int Num, Params &PRS, float** &PopulClass,int*** &Popul,float** &RUPD,int* &CurNumRules,
                   float*** &PopulMu, float** &PopulCF, int NumInst, int* Inst, float*** MRS, sample &Samp,
                   float* &ClassConf, int* &InstClass)
{
    if(PRS.MutType == 0)
    {
        for(int k=0;k!=PRS.MaxNRules;k++)
        {
            if(PopulClass[Num][k] != -1)
            {
                for(int L=0;L!=PRS.NVars;L++)
                {
                    if(Random(0,1) < 1./(float)CurNumRules[Num]/(int)PRS.NVars/3.)
                    {
                        Popul[Num][k][L] = IntRandom(PRS.NFSets);
                        RUPD[Num][k] = 0;
                        CFClassRUPD(Num,k,PRS,NumInst,Inst,Popul,RUPD,MRS,Samp,ClassConf,PopulMu,PopulClass,PopulCF,InstClass);
                    }
                }
            }
        }
    }
    if(PRS.MutType == 1)
    {
        for(int k=0;k!=PRS.MaxNRules;k++)
        {
            if(PopulClass[Num][k] != -1)
            {
                for(int L=0;L!=PRS.NVars;L++)
                {
                    if(Random(0,1) < 1./(float)CurNumRules[Num]/(int)PRS.NVars)
                    {
                        Popul[Num][k][L] = IntRandom(PRS.NFSets);
                        RUPD[Num][k] = 0;
                        CFClassRUPD(Num,k,PRS,NumInst,Inst,Popul,RUPD,MRS,Samp,ClassConf,PopulMu,PopulClass,PopulCF,InstClass);
                    }
                }
            }
        }
    }
    if(PRS.MutType == 2)
    {
        for(int k=0;k!=PRS.MaxNRules;k++)
        {
            if(PopulClass[Num][k] != -1)
            {
                for(int L=0;L!=PRS.NVars;L++)
                {
                    if(Random(0,1) < 3./(float)CurNumRules[Num]/(int)PRS.NVars)
                    {
                        Popul[Num][k][L] = IntRandom(PRS.NFSets);
                        RUPD[Num][k] = 0;
                        CFClassRUPD(Num,k,PRS,NumInst,Inst,Popul,RUPD,MRS,Samp,ClassConf,PopulMu,PopulClass,PopulCF,InstClass);
                    }
                }
            }
        }
    }
    return;
}
void MichOperators(int Num, Params &PRS,float** &PopulClass,int* &CurNumRules,  float* &RuleClassified,
                   int** &GoodRNums, int* &RandomNumbers,int*** &Popul,float** &RUPD,int NumInst, int* Inst,
                   float*** MRS, sample &Samp,float* &ClassConf, float*** &PopulMu, float** &PopulCF,
                   int* &InstClass)
{
    int NewRulePlace=-1;
    for(int i=0;i!=PRS.MaxNRules;i++)
    {
        if(PopulClass[Num][i] == -1)
        {
            NewRulePlace = i;
        }
    }
    if(NewRulePlace == -1)
        return;
    float minfit=0;
    int Taken1;
    int Taken2;
    Taken1=0;
    for(int j=0;j!=PRS.TourSize;j++)
    {
        RandomNumbers[j]=IntRandom(CurNumRules[Num]);
        if(j==0)
            minfit=RuleClassified[ GoodRNums[Num][RandomNumbers[0]] ];
        else
        {
            if(RuleClassified[ GoodRNums[Num][RandomNumbers[j]] ] > minfit)
            {
                Taken1 = RandomNumbers[j];
                minfit = RuleClassified[ GoodRNums[Num][RandomNumbers[j]] ];
            }
        }
    }
    Taken2=0;
    for(int j=0;j!=PRS.TourSize;j++)
    {
        RandomNumbers[j]=IntRandom(CurNumRules[Num]);
        if(j==0)
            minfit=RuleClassified[ GoodRNums[Num][RandomNumbers[0]] ];
        else
        {
            if(RuleClassified[ GoodRNums[Num][RandomNumbers[j]] ] > minfit)
            {
                Taken2 = RandomNumbers[j];
                minfit = RuleClassified[ GoodRNums[Num][RandomNumbers[j]] ];
            }
        }
    }
    if(Random(0,1) < 0.9)
    {
        for(int j=0;j!=PRS.NVars;j++)
        {
            if(Random(0,1) < 0.5)
            {
                Popul[Num][NewRulePlace][j] = Popul[Num][GoodRNums[Num][Taken2]][j];
            }
            else
            {
                Popul[Num][NewRulePlace][j] = Popul[Num][GoodRNums[Num][Taken1]][j];
            }
        }
    }
    for(int L=0;L!=PRS.NVars;L++)
    {
        if(Random(0,1) < 1./(float)PRS.NVars)
        {
            Popul[Num][NewRulePlace][L] = IntRandom(PRS.NFSets);
            RUPD[Num][NewRulePlace] = 0;
            CFClassRUPD(Num,NewRulePlace,PRS,NumInst,Inst,Popul,RUPD,MRS,Samp,ClassConf,PopulMu,PopulClass,PopulCF,InstClass);
        }
    }
}
void HeuristicRG(int Num, Params &PRS, int*** &Popul,float** &RUPD,int NumInst, int* Inst, float*** &MRS,
                 sample &Samp, float* &ClassConf, float*** &PopulMu, float** &PopulCF,float** &PopulClass,
                 int* &InstClass)
{
    int NewRulePlace=-1;
    for(int i=0;i!=PRS.MaxNRules;i++)
    {
        if(PopulClass[Num][i] == -1)
        {
            NewRulePlace = i;
        }
    }
    if(NewRulePlace == -1)
        return;
    int RandomPattern = IntRandom(NumInst);
    float SummOR;
    float RandomPoint;
    for(int i=0;i!=PRS.NVars;i++)
    {
        Popul[Num][NewRulePlace][i] = 0;
        if(Random(0,1) > PRS.PDC )
        {
            SummOR=PRS.NPartitions;                               //Depends on NFSets!
            RandomPoint = Random(0,SummOR);
            do
            {
                Popul[Num][NewRulePlace][i]++;
                //SummOR-=LORS[RandomPattern][i][Popul[Num][NewRulePlace][i]];
                SummOR-=MRS[Inst[RandomPattern]][i][Popul[Num][NewRulePlace][i]];
            } while(SummOR > RandomPoint && Popul[Num][NewRulePlace][i] != 14);
        }
        //cout<<Popul[Num][k][i]<<" ";
    }
    RUPD[Num][NewRulePlace] = 0;
    CFClassRUPD(Num,NewRulePlace,PRS,NumInst,Inst,Popul,RUPD,MRS,Samp,ClassConf,PopulMu,PopulClass,PopulCF,InstClass);
}
void MichiganPart(int Num, Params &PRS,float** &PopulClass, int*** &Popul,int* &CurNumRules, sample &Samp,
                  float* &RuleClassified, int NumInst, int* Inst, float*** &PopulMu, float** &PopulCF,
                  int* &MisClassifiedNums, float* &BadRNums, int** &GoodRNums, float* &tempRuleC,
                  float** &RUPD, float*** &MRS, float* &ClassConf, float* &CurRLength, int* &RandomNumbers,
                  int* &InstClass)
{
    if(PRS.SADJ == 0)
    {
        if(Random(0,1) < PRS.AddProb)
            PRS.CrossType=0;
        else
            PRS.CrossType=3;
    }
    if(PRS.CrossType == 0 || PRS.CrossType == 2)
    {

        float CurConf=0;
        //float MaxConf=0;
        float tempConf=0;
        bool DCRule=false;
        int SummAc;

        for(int k=0;k!=PRS.MaxNRules;k++)
        {
            CurConf=0;
            //MaxConf=0;
            tempConf=0;
            SummAc=0;
            DCRule=true;
            for(int i=0;i!=Samp.NVars;i++)
                SummAc+=Popul[Num][k][i];
            if(SummAc != 0)
            {
                for(int i=0;i!=NumInst;i++)
                {
                    tempConf=1;
                    for(int j=0;j!=Samp.NVars;j++)
                    {
                        if(Popul[Num][k][j] > 0)
                        {
                            DCRule = false;
                            //tempConf *= LORS[i][j][Popul[Num][k][j]];
                            tempConf *= MRS[Inst[i]][j][Popul[Num][k][j]];
                        }
                    }
                    if(DCRule == true)
                    {
                        tempConf=0;
                    }
                    PopulMu[0][k][i] = tempConf;
                    //cout<<PopulMu[Num][k][i]<<" ";
                    CurConf += tempConf;
                }
            }
        }
            int CurNRules = 0;
            float CurRLength1=0;
            for(int i=0;i!=PRS.MaxNRules;i++)
            {
                GoodRNums[Num][i] = 0;
                if(PopulClass[Num][i] >= 0)
                {
                    GoodRNums[Num][(int)CurNRules] = i;
                    CurNRules++;
                    for(int j=0;j!=PRS.NVars;j++)
                    {
                        if(Popul[Num][i][j] > 0)
                        {
                            CurRLength1++;
                        }
                        //cout<<Popul[Num][i][j]<<" ";
                    }
                }
            }
            CurNumRules[Num] = CurNRules;
            for(int i=0;i!=PRS.MaxNRules;i++)
            {
                RuleClassified[i] = 0;
            }
            float maxCFMu;
            int maxnum;
            float Error=0;
            CurNRules=0;
            float SavedError=0;
            for(int i=0;i!=NumInst;i++)
            {
                maxCFMu=0;
                maxnum=-1;
                for(int j=0;j!=PRS.MaxNRules;j++)
                {
                    if(PopulCF[Num][j]*PopulMu[0][j][i] > maxCFMu && PopulCF[Num][j]*PopulMu[0][j][i] > 0)
                    {
                        maxCFMu = PopulCF[Num][j]*PopulMu[0][j][i];
                        maxnum = j;
                        continue;
                    }
                }
                SavedError=Error;
                if(maxnum != -1)
                {
                    if(PopulClass[Num][maxnum] != Samp.Classes[Inst[i]])
                    {
                        MisClassifiedNums[(int)Error] = i;
                        Error++;
                    }
                }
                else
                {
                    MisClassifiedNums[(int)Error] = i;
                    Error++;
                }
                if(Error==SavedError)
                {
                    RuleClassified[maxnum]++;
                }
            }
            for(int i=0;i!=CurNumRules[Num];i++)
            {
                tempRuleC[i] = RuleClassified[GoodRNums[Num][i]];
                BadRNums[i] = GoodRNums[Num][i];
            }
            if(PRS.AddDelete)
                qSort(tempRuleC,BadRNums,0,CurNumRules[Num]-1);

            int NumNewRules=0;
            int NumHeuristic=0;
            int NumGenetic = 0;
            //NumNewRules = IntRandom( ((int)CurNumRules[Num]/5+1)*2 );
            NumNewRules = (int)CurNumRules[Num]/5+1;
            if(PRS.AddDelete)
            {
                for(int i=0;i!=NumNewRules;i++)
                {
                    PopulClass[Num][(int)BadRNums[i]] = -1;
                    for(int j=0;j!=PRS.NVars;j++)
                    {
                        Popul[Num][(int)BadRNums[i]][j] = 0;
                    }
                    RUPD[Num][(int)BadRNums[i]] = 0;
                    CFClassRUPD(Num,(int)BadRNums[i],PRS,NumInst,Inst,Popul,RUPD,MRS,Samp,ClassConf,PopulMu,PopulClass,PopulCF,InstClass);
                }
            }
            CurNRules=0;
            for(int i=0;i!=PRS.MaxNRules;i++)
            {
                GoodRNums[Num][i] = 0;
                if(PopulClass[Num][i] >= 0)
                {
                    GoodRNums[Num][(int)CurNRules] = i;
                    CurNRules++;
                    for(int j=0;j!=PRS.NVars;j++)
                    {
                        if(Popul[Num][i][j] > 0)
                        {
                            CurRLength1++;
                        }
                        //cout<<Popul[Num][i][j]<<" ";
                    }
                }
            }
            CurNumRules[Num] = CurNRules;
            if(CurNRules != 0)
                CurRLength1 /= CurNRules;
            else
                CurRLength1 = -1;
            CurRLength[Num] = CurRLength1;

            if(PRS.MaxNRules - CurNumRules[Num] < NumNewRules)
                NumNewRules = PRS.MaxNRules - CurNumRules[Num];
            NumHeuristic = int((float)NumNewRules/2.)+IntRandom(1);
            PRS.NMisclassified = Error;

            /*NumGenetic=0;
            NumHeuristic=0;
            for(int i=0;i!=NumNewRules;i++)
            {
                if(Random(0,1) < pheu[0])
                    NumHeuristic++;
            }*/
            if(PRS.NMisclassified < NumHeuristic)
            {
                NumHeuristic = PRS.NMisclassified;
            }
            NumGenetic = NumNewRules - NumHeuristic;
            PRS.savednumgen=NumGenetic;
            PRS.savednumheu=NumHeuristic;

            for(int i=0;i!=NumGenetic;i++)
                MichOperators(Num,PRS,PopulClass,CurNumRules,RuleClassified,GoodRNums,RandomNumbers,Popul,
                              RUPD,NumInst,Inst,MRS,Samp,ClassConf,PopulMu,PopulCF,InstClass);
            for(int i=0;i<NumHeuristic;i++)
                HeuristicRG(Num,PRS,Popul,RUPD,NumInst,Inst,MRS,Samp,ClassConf,PopulMu,PopulCF,PopulClass,InstClass);

    }

    if(PRS.SADJ == 0)
    {
        if(Random(0,1) < PRS.DelProb)
            PRS.CrossType=1;
        else
            PRS.CrossType=3;
    }
    if(PRS.CrossType == 1 || PRS.CrossType == 2)
    {

        float CurConf=0;
        //float MaxConf=0;
        float tempConf=0;
        bool DCRule=false;
        int SummAc;

        for(int k=0;k!=PRS.MaxNRules;k++)
        {
            CurConf=0;
            //MaxConf=0;
            tempConf=0;
            SummAc=0;
            DCRule=true;
            for(int i=0;i!=Samp.NVars;i++)
                SummAc+=Popul[Num][k][i];
            if(SummAc != 0)
            {
                for(int i=0;i!=NumInst;i++)
                {
                    tempConf=1;
                    for(int j=0;j!=Samp.NVars;j++)
                    {
                        if(Popul[Num][k][j] > 0)
                        {
                            DCRule = false;
                            //tempConf *= LORS[i][j][Popul[Num][k][j]];
                            tempConf *= MRS[Inst[i]][j][Popul[Num][k][j]];
                        }
                    }
                    if(DCRule == true)
                    {
                        tempConf=0;
                    }
                    PopulMu[0][k][i] = tempConf;
                    //cout<<PopulMu[Num][k][i]<<" ";
                    CurConf += tempConf;
                }
            }
        }
            int CurNRules = 0;
            float CurRLength1=0;
            for(int i=0;i!=PRS.MaxNRules;i++)
            {
                GoodRNums[Num][i] = 0;
                if(PopulClass[Num][i] >= 0)
                {
                    GoodRNums[Num][(int)CurNRules] = i;
                    CurNRules++;
                    for(int j=0;j!=PRS.NVars;j++)
                    {
                        if(Popul[Num][i][j] > 0)
                        {
                            CurRLength1++;
                        }
                        //cout<<Popul[Num][i][j]<<" ";
                    }
                }
            }
            CurNumRules[Num] = CurNRules;
            for(int i=0;i!=PRS.MaxNRules;i++)
            {
                RuleClassified[i] = 0;
            }
            float maxCFMu;
            int maxnum;
            float Error=0;
            CurNRules=0;
            float SavedError=0;
            for(int i=0;i!=NumInst;i++)
            {
                maxCFMu=0;
                maxnum=-1;
                for(int j=0;j!=PRS.MaxNRules;j++)
                {
                    if(PopulCF[Num][j]*PopulMu[0][j][i] > maxCFMu && PopulCF[Num][j]*PopulMu[0][j][i] > 0)
                    {
                        maxCFMu = PopulCF[Num][j]*PopulMu[0][j][i];
                        maxnum = j;
                        continue;
                    }
                }
                SavedError=Error;
                if(maxnum != -1)
                {
                    if(PopulClass[Num][maxnum] != Samp.Classes[Inst[i]])
                    {
                        MisClassifiedNums[(int)Error] = i;
                        Error++;
                    }
                }
                else
                {
                    MisClassifiedNums[(int)Error] = i;
                    Error++;
                }
                if(Error==SavedError)
                {
                    RuleClassified[maxnum]++;
                }
            }
            for(int i=0;i!=CurNumRules[Num];i++)
            {
                tempRuleC[i] = RuleClassified[GoodRNums[Num][i]];
                BadRNums[i] = GoodRNums[Num][i];
            }
            qSort(tempRuleC,BadRNums,0,CurNumRules[Num]-1);

            int NumNewRules = (int)CurNumRules[Num]/5+1;

            for(int i=0;i!=NumNewRules;i++)
            {
                PopulClass[Num][(int)BadRNums[i]] = -1;
                for(int j=0;j!=PRS.NVars;j++)
                {
                    Popul[Num][(int)BadRNums[i]][j] = 0;
                }
                RUPD[Num][(int)BadRNums[i]] = 0;
                CFClassRUPD(Num,(int)BadRNums[i],PRS,NumInst,Inst,Popul,RUPD,MRS,Samp,ClassConf,PopulMu,
                            PopulClass,PopulCF,InstClass);
            }
            CurNRules=0;
            for(int i=0;i!=PRS.MaxNRules;i++)
            {
                GoodRNums[Num][i] = 0;
                if(PopulClass[Num][i] >= 0)
                {
                    GoodRNums[Num][(int)CurNRules] = i;
                    CurNRules++;
                    for(int j=0;j!=PRS.NVars;j++)
                    {
                        if(Popul[Num][i][j] > 0)
                        {
                            CurRLength1++;
                        }
                        //cout<<Popul[Num][i][j]<<" ";
                    }
                }
            }
            CurNumRules[Num] = CurNRules;
            if(CurNRules != 0)
                CurRLength1 /= CurNRules;
            else
                CurRLength1 = -1;
            CurRLength[Num] = CurRLength1;
    }
    if(PRS.SADJ == 0)
        PRS.CrossType=2;
}
void LearnErrorCalc(sample &Samp, Params &PRS, float* &NumClassInst, int** &OverBestInd, float* &OverBestClasses,
                    float* &OverBestCF, float** &OverBestMu, float* &ClassConf, int CVLearnSize, float*** &MRS,
                    float &OverBestAve,float* &OverErrorClass, float &OverBestError, int FoldOnTest,
                    int* &InstClass, float** &ConfMatrix)
{
    //float CurConf=0;
    //float MaxConf=0;
    //float savedOverBestMu;
    float tempConf=0;
    int maxnum=0;
    bool DCRule=false;
    int* InstNums;
    InstNums = new int[CVLearnSize];
    int counter = 0;
    for(int i=0;i!=Samp.Size;i++)
    {
        if(Samp.CVFoldNum[i] != FoldOnTest)
        {
            InstNums[counter] = i;
            counter++;
        }
    }
    for(int k=0;k!=PRS.MaxNRules;k++)
    {
        //CurConf=0;
        //MaxConf=0;
        maxnum=0;
        tempConf=0;
        DCRule=true;
        for(int i=0;i!=Samp.NClasses;i++)
            ClassConf[i]=0;
        //cout<<endl;
        for(int i=0;i!=CVLearnSize;i++)
        {
            tempConf=1;
            for(int j=0;j!=PRS.NVars;j++)
            {
                if(OverBestInd[k][j] > 0)
                {
                    DCRule = false;
                    //tempConf *= TORS[i][j][Popul[Num][k][j]];
                    tempConf *= MRS[InstNums[i]][j][OverBestInd[k][j]];
                }
            }
            if(DCRule == true)
            {
                tempConf=0;
            }
            //OverBestMu[k][i] = tempConf;
            OverBestMu[k][InstNums[i]] = tempConf;
            //cout<<OverBestMu[k][i]<<" ";
        }
    }

    float maxCFMu;
    float Error=0;
    for(int i=0;i!=Samp.NClasses;i++)
    {
        OverErrorClass[i] = 0;
    }
    for(int i=0;i!=Samp.NClasses;i++)
        for(int j=0;j!=Samp.NClasses+1;j++)
        {
            ConfMatrix[i][j] = 0;
        }
    for(int i=0;i!=CVLearnSize;i++)
    {
        maxCFMu=0;
        maxnum=-1;
        for(int j=0;j!=PRS.MaxNRules;j++)
        {
            if(OverBestCF[j]*OverBestMu[j][InstNums[i]] > maxCFMu && OverBestCF[j]*OverBestMu[j][InstNums[i]] > 0)
            {
                maxCFMu = OverBestCF[j]*OverBestMu[j][InstNums[i]];
                maxnum = j;
                continue;
            }
            if(OverBestCF[j]*OverBestMu[j][InstNums[i]] == maxCFMu && OverBestCF[j]*OverBestMu[j][InstNums[i]] > 0)
            {
                maxnum = -1;
                break;
            }
        }
        if(maxnum != -1)
        {
            if(OverBestClasses[maxnum] != Samp.Classes[InstNums[i]])
            {
                //cout<<InstNums[i]<<" ";
                Error++;
                OverErrorClass[ Samp.Classes[InstNums[i]] ] ++;
            }
            ConfMatrix[ Samp.Classes[InstNums[i]] ][ (int)OverBestClasses[maxnum] ] ++;
        }
        else
        {
            //cout<<InstNums[i]<<" ";
            Error++;
            OverErrorClass[ Samp.Classes[InstNums[i]] ] ++;
            ConfMatrix[ Samp.Classes[InstNums[i]] ][ Samp.NClasses ] ++;
        }
    }
    //cout<<endl;
    //CurNumRules[0]=0;
    //CurRLength[0] =0;
    float CurNumRules = 0;
    float CurRLength = 0;
    for(int i=0;i!=PRS.MaxNRules;i++)
    {
        if(OverBestClasses[i] >= 0)
        {
            CurNumRules++;
            for(int j=0;j!=Samp.NVars;j++)
            {
                if(OverBestInd[i][j] > 0)
                {
                    CurRLength++;
                }
                //cout<<Popul[Num][i][j]<<" ";
            }
            //cout<<"-> "<<PopulClass[Num][i]<<" "<<(PopulCF[Num][i])<<endl;
        }
    }
    if(CurNumRules != 0)
        CurRLength /= CurNumRules;
    else
        CurRLength = -1;
    //ErrMass[0] = Error;
    //AveMass[0]=0;
    OverBestAve = 0;
    OverBestError = Error;
    for(int i=0;i!=Samp.NClasses;i++)
    {
        OverErrorClass[i] /= InstClass[i];
        OverBestAve += OverErrorClass[i]/(float)Samp.NClasses;
    }
    delete InstNums;
    //if(PRS.CritType == 0)
        //FitMass[0] = Error/NumInst*10000 + 1*CurNumRules[0] + 1*CurNumRules[0]*CurRLength[0];
    //else
        //FitMass[0] = AveMass[0]*10000 + 1*CurNumRules[0] + 1*CurNumRules[0]*CurRLength[0];
    //FitMass[Num] = 1./(1.+FitMass[Num]);
}
void TestErrorCalc(sample &Samp, Params &PRS, float* &NumClassInst, int** &OverBestInd, float* &OverBestClasses,
                    float* &OverBestCF, float** &OverBestMu, float* &ClassConf, int CVTestSize, float*** &MRS,
                    float &OverBestAve,float* &OverErrorClass, float &OverBestError, int FoldOnTest,
                    int* &InstClass, float** &ConfMatrixTest)
{
    //float CurConf=0;
    //float MaxConf=0;
    float tempConf=0;
    int maxnum=0;
    bool DCRule=false;
    int* InstNums;
    InstNums = new int[CVTestSize];
    int counter = 0;
    for(int i=0;i!=Samp.Size;i++)
    {
        if(Samp.CVFoldNum[i] == FoldOnTest)
        {
            InstNums[counter] = i;
            counter++;
        }
    }
    if(counter != CVTestSize)
        cout<<"sdfhsdfgdsgfndfgnDFgnD";
    for(int k=0;k!=PRS.MaxNRules;k++)
    {
        //CurConf=0;
        //MaxConf=0;
        maxnum=0;
        tempConf=0;
        DCRule=true;
        for(int i=0;i!=Samp.NClasses;i++)
            ClassConf[i]=0;
        //cout<<endl;
        for(int i=0;i!=CVTestSize;i++)
        {
            tempConf=1;
            for(int j=0;j!=PRS.NVars;j++)
            {
                if(OverBestInd[k][j] > 0)
                {
                    DCRule = false;
                    //tempConf *= TORS[i][j][Popul[Num][k][j]];
                    tempConf *= MRS[InstNums[i]][j][OverBestInd[k][j]];
                }
            }
            if(DCRule == true)
            {
                tempConf=0;
            }
            OverBestMu[k][InstNums[i]] = tempConf;
        }
    }

    float maxCFMu;
    float Error=0;
    for(int i=0;i!=Samp.NClasses;i++)
    {
        OverErrorClass[i] = 0;
    }
    for(int i=0;i!=Samp.NClasses;i++)
        for(int j=0;j!=Samp.NClasses+1;j++)
        {
            ConfMatrixTest[i][j] = 0;
        }
    for(int i=0;i!=CVTestSize;i++)
    {
        maxCFMu=0;
        maxnum=-1;
        for(int j=0;j!=PRS.MaxNRules;j++)
        {
            if(OverBestCF[j]*OverBestMu[j][InstNums[i]] > maxCFMu && OverBestCF[j]*OverBestMu[j][InstNums[i]] > 0)
            {
                maxCFMu = OverBestCF[j]*OverBestMu[j][InstNums[i]];
                maxnum = j;
                continue;
            }
            if(OverBestCF[j]*OverBestMu[j][InstNums[i]] == maxCFMu && OverBestCF[j]*OverBestMu[j][InstNums[i]] > 0)
            {
                maxnum = -1;
                break;
            }
        }
        if(maxnum != -1)
        {
            if(OverBestClasses[maxnum] != Samp.Classes[InstNums[i]])
            {
                Error++;
                OverErrorClass[ Samp.Classes[InstNums[i]] ] ++;
            }
            ConfMatrixTest[ Samp.Classes[InstNums[i]] ][ (int)OverBestClasses[maxnum] ] ++;
        }
        else
        {
            Error++;
            OverErrorClass[ Samp.Classes[InstNums[i]] ] ++;
            ConfMatrixTest[ Samp.Classes[InstNums[i]] ][ Samp.NClasses ] ++;
        }
    }
    //CurNumRules[0]=0;
    //CurRLength[0] =0;
    float CurNumRules = 0;
    float CurRLength = 0;
    for(int i=0;i!=PRS.MaxNRules;i++)
    {
        if(OverBestClasses[i] >= 0)
        {
            CurNumRules++;
            for(int j=0;j!=Samp.NVars;j++)
            {
                if(OverBestInd[i][j] > 0)
                {
                    CurRLength++;
                }
                //cout<<Popul[Num][i][j]<<" ";
            }
            //cout<<"-> "<<PopulClass[Num][i]<<" "<<(PopulCF[Num][i])<<endl;
        }
    }
    if(CurNumRules != 0)
        CurRLength /= CurNumRules;
    else
        CurRLength = -1;
    //ErrMass[0] = Error;
    //AveMass[0]=0;
    OverBestAve = 0;
    OverBestError = Error;
    for(int i=0;i!=Samp.NClasses;i++)
    {
        OverErrorClass[i] /= InstClass[i];
        OverBestAve += OverErrorClass[i]/(float)Samp.NClasses;
    }
    delete InstNums;
    //if(PRS.CritType == 0)
        //FitMass[0] = Error/NumInst*10000 + 1*CurNumRules[0] + 1*CurNumRules[0]*CurRLength[0];
    //else
        //FitMass[0] = AveMass[0]*10000 + 1*CurNumRules[0] + 1*CurNumRules[0]*CurRLength[0];
    //FitMass[Num] = 1./(1.+FitMass[Num]);
}
void FindNSaveBest(bool HasBestFit,Params &PRS, float* &NumMass,float* &FitMassCopy, int*** &Popul,
                   float** &PopulClass,float*** &PopulMu, float** &PopulCF,float** &RUPD, float* &AveMass,
                   float* &ErrMass,float* &FitMass,int* &CurNumRules, int** &BestInd, float* &BestClasses,
                   float* &BestCF, float** &BestMu, int NumInst, float* &CurRLength, bool UpdateBest,
                   int** &OverBestInd, float* &OverBestClasses, float* &OverBestCF, float &OverBestError,
                   int** &OverBestInd2, float* &OverBestClasses2, float* &OverBestCF2, sample &Samp,
                   float** &OverBestMu, float* &ClassConf, int &CVLearnSize, float*** &MRS,
                   float &OverBestAve,float* &OverErrorClass,int FoldOnTest,float* &NumClassInst,
                   float** &ErrorClass, float** &ConfMatrix, int* &Inst, int* &InstUsed, int* &InstClass,
                   int* &counterClass)
{
    if(HasBestFit == 1)
    {
        if(UpdateBest)
        {
            //PRS.BestFit = FitMass[0];
            //PRS.BestAve = AveMass[0];
            //HasBestFit=true;
        }
        if(PRS.formType == 1)
        {
            for(int i=0;i!=PRS.NumInds*2;i++)
            {
                NumMass[i] = i;
                FitMassCopy[i] = FitMass[i];
            }
            qSort(FitMassCopy,NumMass,0,PRS.NumInds*2-1);
        }
        else
        {
            for(int i=0;i!=PRS.NumInds;i++)
            {
                NumMass[i] = PRS.NumInds+i;
            }
        }
        qSort(NumMass,FitMassCopy,0,PRS.NumInds-1);

        for(int i=0;i!=PRS.NumInds;i++)
        {
            for(int j=0;j!=PRS.MaxNRules;j++)
            {
                for(int k=0;k!=PRS.NVars;k++)
                {
                    Popul[i][j][k] = Popul[(int)NumMass[i]][j][k];
                }
                PopulCF[i][j] = PopulCF[(int)NumMass[i]][j];
                PopulClass[i][j] = PopulClass[(int)NumMass[i]][j];
                for(int k=0;k!=NumInst;k++)
                {
                    //PopulMu[i][j][0] = PopulMu[(int)NumMass[i]][j][0];
                }
                RUPD[i][j] = RUPD[(int)NumMass[i]][j];
            }
            FitMass[i] = FitMass[(int)NumMass[i]];
            AveMass[i] = AveMass[(int)NumMass[i]];
            ErrMass[i] = ErrMass[(int)NumMass[i]];
        }

        for(int Num=0;Num!=PRS.NumInds;Num++)
        {
            float CurNRules=0;
            float CurRLength1=0;
            for(int i=0;i!=PRS.MaxNRules;i++)
            {
                if(PopulClass[Num][i] >= 0)
                {
                    CurNRules++;
                    for(int j=0;j!=PRS.NVars;j++)
                    {
                        if(Popul[Num][i][j] > 0)
                        {
                            CurRLength1++;
                        }
                    }
                }
            }
            CurNumRules[Num] = CurNRules;
            if(CurNRules != 0)
                CurRLength1 /= CurNRules;
            else
                CurRLength1 = -1;
            CurRLength[Num] = CurRLength1;
        }
    }
    //float MinFit;
    if(HasBestFit == false)
    {
         PRS.BestFit = FitMass[0];
         PRS.BestAve = AveMass[0];
         //HasBestFit=true;
    }
    for(int i=0;i!=PRS.NumInds;i++)
    {
        if(FitMass[i] <= PRS.BestFit)
        {
            PRS.BestFit = FitMass[i];
            PRS.BestAve = AveMass[i];
            PRS.BestNumRules = CurNumRules[i];
            PRS.BestRLength = CurRLength[i];
            PRS.BestErr = ErrMass[i];
            for(int j=0;j!=PRS.MaxNRules;j++)
            {
                for(int k=0;k!=PRS.NVars;k++)
                {
                    BestInd[j][k] = Popul[i][j][k];
                }
                for(int k=0;k!=NumInst;k++)
                {
                    //BestMu[j][k] = PopulMu[i][j][0];
                }
                BestCF[j] = PopulCF[i][j];
                BestClasses[j] = PopulClass[i][j];
            }
        }
    }
    FitMass[0] = PRS.BestFit;
    AveMass[0] = PRS.BestAve;
    for(int j=0;j!=PRS.MaxNRules;j++)
    {
        for(int k=0;k!=PRS.NVars;k++)
        {
            Popul[0][j][k] = BestInd[j][k];
            //cout<<BestInd[j][k]<<" ";
        }
        for(int k=0;k!=NumInst;k++)
        {
            //PopulMu[0][j][0] = BestMu[j][k];
        }
        PopulCF[0][j] = BestCF[j];
        //cout<<BestCF[j]<<" "<<BestClasses[j]<<endl;
        PopulClass[0][j] = BestClasses[j];
    }

    //cout<<OverBestError<<"\t";
    float SavedOverBestAve = OverBestAve;
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
            OverBestInd[j][k] = BestInd[j][k];
        }
        OverBestClasses[j] = BestClasses[j];
        OverBestCF[j] = BestCF[j];
    }
    //cout<<"asdfa";
    LearnErrorCalc(Samp,PRS,NumClassInst,OverBestInd,OverBestClasses,OverBestCF,OverBestMu,ClassConf,
                   CVLearnSize,MRS,OverBestAve,OverErrorClass,OverBestError,FoldOnTest,counterClass,ConfMatrix);
    //cout<<OverBestError<<"\t";
    //cout<<endl;
    if(PRS.CritType == 0)
    {
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
    }
    else
    {
        if(OverBestAve > SavedOverBestAve)
        {
            OverBestAve = SavedOverBestAve;
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
    }

    for(int j=0;j!=PRS.MaxNRules;j++)
    {
        for(int k=0;k!=PRS.NVars;k++)
        {
            Popul[1][j][k] = OverBestInd[j][k];
            //cout<<BestInd[j][k]<<" ";
        }
        PopulCF[1][j] = OverBestCF[j];
        //cout<<BestCF[j]<<" "<<BestClasses[j]<<endl;
        PopulClass[1][j] = OverBestClasses[j];
        RUPD[1][j] = 0;
        CFClassRUPD(1,j,PRS,NumInst,Inst,Popul,RUPD,MRS,Samp,ClassConf,PopulMu,PopulClass,PopulCF,InstClass);
    }
    CFClassRB(1,Samp,PRS,RUPD,Popul,ClassConf,PopulMu,MRS,PopulClass,PopulCF,NumInst,Inst,InstClass);
    FitCalc(1,Samp,ErrorClass,ConfMatrix,NumInst,PRS,PopulCF,PopulMu,PopulClass,CurNumRules,CurRLength,ErrMass,AveMass,
            FitMass,NumClassInst,Popul,Inst,0,InstUsed,InstClass,MRS);
}
void CreateNewPopul(int Gen,Params &PRS,sample &Samp, float* &RuleClassified,
    float* &RuleP,    float* &tempRuleC,    float* &BadRNums,    int* &MisClassifiedNums,    int* &RandomNumbers,
    int* &CurNumRules,    float* &CurRLength,    int** &GoodRNums,    float** &ChosenRules,    float* &FitMass,
    float* &FitMassCopy,    float* &PropSetSel,    float* &ClassConf,    float** &PopulClass,    float** &PopulCF,
    float*** &PopulMu,    int** &BestInd,    float* &BestClasses,    float* &BestCF,    float** &BestMu,
    float** &RUPD, float* &ErrMass,    int* &SelectedNums1,    int* &SelectedNums2,    float* &NumMass,
    float* &PMass,    float* &NumMass2, float* &RankMass,    float* &RMass,    float* &SavedFit,
    float* &NumClassInst,    float* &NumClassInstWrong, float* &NumClassInstTest,    float* &NumClassInstTestWrong,
    float** &ConfMatrix,    float** &ConfMatrixTest, float** &ErrorClass,    float* &AveMass,
    float* &TestErrorClass,    int*** &Popul,    float* &psel,    float* &pcrs,  float* &pmut,
    float* &pheu,    float* &usedsel,    float* &usedcrs,    float* &usedmut,    float* &usedheu,
    float* &sumfitsel,    float* &sumfitcrs,    float* &sumfitmut,    float* &sumfitheu,    float* &numselused,
    float* &numcrsused,    float* &nummutused,    float* &numheuused,    float* &sucsel,    float* &succrs,
    float* &sucmut,    float* &sucheu, int &NumInst, int* &InstUsed, int* &Inst, int* &InstClass, int &FoldOnTest,
    float*** &MRS, bool UpdateBest, int** &OverBestInd, float* &OverBestClasses, float* &OverBestCF,
    float &OverBestError,int** &OverBestInd2, float* &OverBestClasses2, float* &OverBestCF2,
    float** OverBestMu, int CVLearnSize,float &OverBestAve,float* &OverErrorClass, int* &counterClass)
{
    for(int i=0;i!=PRS.NumInds;i++)
    {
        SavedFit[i]=FitMass[i];
    }
    float MinFit=MinElVal(FitMass,PRS.NumInds);
    float MaxFit=MaxElVal(FitMass,PRS.NumInds);
    if(MinFit != MaxFit)
    {
        for(int i=0;i!=PRS.NumInds;i++)
            PMass[i]=(FitMass[i]-MinFit)/(MaxFit-MinFit);
            //PMass[i]=-FitMass[i]+MaxFit+MinFit;
        float Sum=SumEl(PMass,PRS.NumInds);
        for(int i=0;i!=PRS.NumInds;i++)
        {
            PMass[i]=PMass[i]/Sum;
            if(i>0)
                PMass[i]=PMass[i]+PMass[i-1];
        }
    }
    else
        for(int i=0;i!=PRS.NumInds;i++)
        {
            PMass[i]=1.0/PRS.NumInds;
            if(i>0)
                PMass[i]=PMass[i]+PMass[i-1];
        }
    //ofstream foutRanks("Ranks.txt",ios::trunc);
    for(int i=0;i!=PRS.NumInds;i++)
    {
        SavedFit[i]=FitMass[i];
        RankMass[i]=FitMass[i];
        NumMass2[i]=(float)i;
    }
    QS(0,RankMass,NumMass2,PRS);
    float NumEq=1;
    float SumEq=0;
    float Prev=RankMass[0];
    //cout<<endl<<endl;
    for(int i=1;i!=PRS.NumInds+1;i++)
    {
        if(i<PRS.NumInds)
        {
            if(RankMass[i] == Prev)
            {
                NumEq++;
            }
            else
            {
                for(int j=i-NumEq;j!=i;j++)
                {
                    SumEq+=(j);
                }
                SumEq=(SumEq)/NumEq;
                for(int j=i-NumEq;j!=i;j++)
                {
                    RankMass[j]=SumEq;
                }
                SumEq=0;
                NumEq=1;
                Prev=RankMass[i];
            }
        }
        else
        {
            for(int j=i-NumEq;j!=i;j++)
            {
                SumEq+=(j);
            }
            SumEq=(SumEq)/NumEq;
            for(int j=i-NumEq;j!=i;j++)
            {
                RankMass[j]=SumEq;
            }
        }
    }
    //for(int i=0;i!=NumInds;i++)
        //cout<<RankMass[i]<<endl;
    //for(int i=0;i!=NumInds;i++)
        //RankMass[i]=(float)i+1;
    QS(1,RankMass,NumMass2,PRS);

    float MinRank=MinElVal(RankMass,PRS.NumInds);
    float MaxRank=MaxElVal(RankMass,PRS.NumInds);

    for(int i=0;i!=PRS.NumInds;i++)
      //  RMass[i]=(RankMass[i]-MinRank)/(MaxRank-MinRank);
      RMass[i]=-RankMass[i]+MaxRank+MinRank;
    if(MinRank != MaxRank)
    {
        float Sum=SumEl(RMass,PRS.NumInds);
        for(int i=0;i!=PRS.NumInds;i++)
        {
            RMass[i]=RMass[i]/Sum;
            if(i>0)
                RMass[i]=RMass[i]+RMass[i-1];
        }
    }
    else
    {
        for(int i=0;i!=PRS.NumInds;i++)
        {
            RMass[i]=1.0/PRS.NumInds;
            if(i>0)
                RMass[i]=RMass[i]+RMass[i-1];
        }
    }

    for(int i=0;i!=PRS.NumInds;i++)
    {
        //foutRanks<<FitMass[i]<<"\t"<<PMass[i]<<"\t"<<RMass[i]<<"\t"<<RankMass[i]<<endl;
    }

    //float MinR=MinElVal(RMass,NumInds);
    //float MaxR=MaxElVal(RMass,NumInds);
    //for(int i=0;i!=NumInds;i++)
      //  RMass[i]=(RMass[i]-MinR)/(MaxR-MinR);

    //_getch();

    //ofstream foutSA("SA.txt",ios::app);

    for(int i=0;i!=PRS.NumSel;i++)
    {
        sumfitsel[i] = 0;
        numselused[i] = 0;
    }
    for(int i=0;i!=PRS.NumCrs;i++)
    {
        sumfitcrs[i] = 0;
        numcrsused[i] = 0;
    }
    for(int i=0;i!=PRS.NumMut;i++)
    {
        sumfitmut[i] = 0;
        nummutused[i] = 0;
    }
    for(int i=0;i!=PRS.NumHeu;i++)
    {
        sumfitheu[i] = 0;
        numheuused[i] = 0;
    }

    for(int i=0;i!=PRS.NumInds;i++)
    {
        PRS.savednumgen=0;
        PRS.savednumheu=0;
        if(PRS.SADJ == 1 || PRS.SADJ == 2)
        {
            float tempRand=Random(0,1);
            if(tempRand <= psel[0])
                PRS.SelType=0;
            else PRS.SelType=1;
            /*
            if(tempRand <= psel[0])
                SelType=0;
            else if(tempRand <= psel[0] + psel[1])
                SelType=1;
            else
                SelType=2;*/

            tempRand=Random(0,1);
            if(tempRand <= pcrs[0])
                PRS.CrossType=0;
            else if(tempRand <= pcrs[0] + pcrs[1])
                PRS.CrossType=1;
            else
                PRS.CrossType=2;
            tempRand=Random(0,1);
            if(tempRand <= pmut[0])
                PRS.MutType=0;
            else if(tempRand <= pmut[0] + pmut[1])
                PRS.MutType=1;
            else
                PRS.MutType=2;
            tempRand=Random(0,1);
        }
        switch(PRS.SelType)
        {
            case 0:
            {
                //PropSelection(i);
                RankSelection(i,PRS,PMass,RMass,SelectedNums1,SelectedNums2);
                break;
            }
            case 1:
            {
                //RankSelection(i);
                PRS.TourSize=5;
                TourSelection(i,PRS,FitMass,SelectedNums1,SelectedNums2,RandomNumbers);
                break;
            }
            case 2:
            {
                PRS.TourSize=5;
                TourSelection(i,PRS,FitMass,SelectedNums1,SelectedNums2,RandomNumbers);
                break;
            }
        }
        PittsCrossover(i,PRS,SelectedNums1,SelectedNums2,GoodRNums,FitMass,ChosenRules,CurNumRules,RUPD,
                       Popul,PopulMu,PopulClass,PopulCF,NumInst);
        PittsMutation(i+PRS.NumInds,PRS,PopulClass,Popul,RUPD,CurNumRules,PopulMu,PopulCF,NumInst,Inst,
                      MRS,Samp,ClassConf,InstClass);
        MichiganPart(i+PRS.NumInds,PRS,PopulClass,Popul,CurNumRules,Samp,RuleClassified,NumInst,Inst,PopulMu,
                     PopulCF,MisClassifiedNums,BadRNums,GoodRNums,tempRuleC,RUPD,MRS,ClassConf,CurRLength,RandomNumbers,InstClass);
        FitCalc(PRS.NumInds+i,Samp,ErrorClass,ConfMatrix,NumInst,PRS,PopulCF,PopulMu,PopulClass,CurNumRules,CurRLength,
                ErrMass,AveMass,FitMass,NumClassInst,Popul,Inst,0,InstUsed,InstClass,MRS);

        if(PRS.SADJ == 1)
        {
            //sumfitsel[SelType]+=1./(1.-FitMass[i]);
            //sumfitcrs[CrossType]+=1./(1.-FitMass[i]);
            //sumfitmut[MutType]+=1./(1.-FitMass[i]);
            numselused[PRS.SelType]++;
            numcrsused[PRS.CrossType]++;
            nummutused[PRS.MutType]++;
            sumfitsel[PRS.SelType]+=(FitMass[i]);
            sumfitcrs[PRS.CrossType]+=(FitMass[i]);
            sumfitmut[PRS.MutType]+=(FitMass[i]);
            if(PRS.savednumheu > PRS.savednumgen)
            {
                numheuused[0]++;
                sumfitheu[0]+=(FitMass[i]);
            }
            else if(PRS.savednumgen > PRS.savednumheu)
            {
                numheuused[1]++;
                sumfitheu[1]+=(FitMass[i]);
            }
            else if(PRS.savednumgen != 0 && PRS.savednumheu != 0)
            {
                numheuused[0]++;
                numheuused[1]++;
                sumfitheu[0]+=(FitMass[i]);
                sumfitheu[1]+=(FitMass[i]);
            }
        }
        if(PRS.SADJ == 2)
        {
            usedsel[PRS.SelType]++;
            usedcrs[PRS.CrossType]++;
            usedmut[PRS.MutType]++;
            if(FitMass[i] > SavedFit[i])
            {
                sucsel[PRS.SelType]++;
                succrs[PRS.CrossType]++;
                sucmut[PRS.MutType]++;
            }
        }
    }
    if(PRS.SADJ == 1)
    {
        for(int i=0;i!=PRS.NumSel;i++)
        {
            sumfitsel[i] /= numselused[i];
            //sumfitsel[i] /= psel[i];
        }
        for(int i=0;i!=PRS.NumCrs;i++)
        {
            sumfitcrs[i] /= numcrsused[i];
            //sumfitcrs[i] /= pcrs[i];
        }
        for(int i=0;i!=PRS.NumMut;i++)
        {
            sumfitmut[i] /= nummutused[i];
            //sumfitmut[i] /= pmut[i];
        }
        for(int i=0;i!=PRS.NumHeu;i++)
        {
            sumfitheu[i] /= numheuused[i];
        }
        int winsel;
        int wincrs;
        int winmut;
        int winheu;
        if( sumfitsel[0] > sumfitsel[1] )
            winsel = 0;
        else if(sumfitsel[1] > sumfitsel[0])
            winsel = 1;
        else
            winsel = 1;
        //if( (sumfitsel[0] > sumfitsel[1] && sumfitsel[0] >= sumfitsel[2]) || (sumfitsel[0] >= sumfitsel[1] && sumfitsel[0] > sumfitsel[2]))
        //    winsel = 0;
        //else if((sumfitsel[1] > sumfitsel[0] && sumfitsel[1] >= sumfitsel[2]) || (sumfitsel[1] >= sumfitsel[0] && sumfitsel[1] > sumfitsel[2]))
        //    winsel = 1;
        //else if((sumfitsel[2] > sumfitsel[1] && sumfitsel[2] >= sumfitsel[0]) || (sumfitsel[2] >= sumfitsel[1] && sumfitsel[2] > sumfitsel[0]))
        //    winsel = 2;
        //else winsel=-1;

        if( (sumfitcrs[0] > sumfitcrs[1] && sumfitcrs[0] >= sumfitcrs[2]) || (sumfitcrs[0] >= sumfitcrs[1] && sumfitcrs[0] > sumfitcrs[2]))
            wincrs = 0;
        else if((sumfitcrs[1] > sumfitcrs[0] && sumfitcrs[1] >= sumfitcrs[2]) || (sumfitcrs[1] >= sumfitcrs[0] && sumfitcrs[1] > sumfitcrs[2]))
            wincrs = 1;
        else if((sumfitcrs[2] > sumfitcrs[1] && sumfitcrs[2] >= sumfitcrs[0]) || (sumfitcrs[2] >= sumfitcrs[1] && sumfitcrs[2] > sumfitcrs[0]))
            wincrs = 2;
        else wincrs=-1;

        if( (sumfitmut[0] > sumfitmut[1] && sumfitmut[0] >= sumfitmut[2]) || (sumfitmut[0] >= sumfitmut[1] && sumfitmut[0] > sumfitmut[2]))
            winmut = 0;
        else if((sumfitmut[1] > sumfitmut[0] && sumfitmut[1] >= sumfitmut[2]) || (sumfitmut[1] >= sumfitmut[0] && sumfitmut[1] > sumfitmut[2]))
            winmut = 1;
        else if((sumfitmut[2] > sumfitmut[1] && sumfitmut[2] >= sumfitmut[0]) || (sumfitmut[2] >= sumfitmut[1] && sumfitmut[2] > sumfitmut[0]))
            winmut = 2;
        else winmut=-1;

        if(sumfitheu[0] > sumfitheu[1])
            winheu = 0;
        else if(sumfitheu[1] > sumfitheu[0])
            winheu = 1;
        else winheu=-1;

        bool breakcur=0;
        if(winsel != -1)
        {
            for(int i=0;i!=PRS.NumSel;i++)
            {
                if(i!=winsel)
                    if(psel[i] <= 1./PRS.K2)
                        breakcur=1;
            }
            if(breakcur == 0)
            {
                for(int i=0;i!=PRS.NumSel;i++)
                {
                    if(i==winsel)
                        psel[i] = psel[i] + (PRS.NumSel-1.)/(PRS.NumSel*PRS.K1);
                    else
                        psel[i] = psel[i] - (1.)/(PRS.NumSel*PRS.K1);
                }
            }
        }
        breakcur=0;
        if(wincrs != -1)
        {
            for(int i=0;i!=PRS.NumCrs;i++)
            {
                if(i!=wincrs)
                    if(pcrs[i] <= 1./PRS.K2)
                        breakcur=1;
            }
            if(breakcur == 0)
            {
                for(int i=0;i!=PRS.NumCrs;i++)
                {
                    if(i==wincrs)
                        pcrs[i] = pcrs[i] + (PRS.NumCrs-1.)/(PRS.NumCrs*PRS.K1);
                    else
                        pcrs[i] = pcrs[i] - (1.)/(PRS.NumCrs*PRS.K1);
                }
            }
        }
        breakcur=0;
        if(winmut != -1)
        {
            for(int i=0;i!=PRS.NumMut;i++)
            {
                if(i!=winmut)
                    if(pmut[i] <= 1./PRS.K2)
                        breakcur=1;
            }
            if(breakcur == 0)
            {
                for(int i=0;i!=PRS.NumMut;i++)
                {
                    if(i==winmut)
                        pmut[i] = pmut[i] + (PRS.NumMut-1.)/(PRS.NumMut*PRS.K1);
                    else
                        pmut[i] = pmut[i] - (1.)/(PRS.NumMut*PRS.K1);
                }
            }
        }
        breakcur=0;
        if(winheu != -1)
        {
            for(int i=0;i!=PRS.NumHeu;i++)
            {
                if(i!=winheu)
                    if(pheu[i] <= 1./PRS.K2)
                        breakcur=1;
            }
            if(breakcur == 0)
            {
                for(int i=0;i!=PRS.NumHeu;i++)
                {
                    if(i==winheu)
                        pheu[i] = pheu[i] + (PRS.NumHeu-1.)/(PRS.NumHeu*PRS.K1);
                    else
                        pheu[i] = pheu[i] - (1.)/(PRS.NumHeu*PRS.K1);
                }
            }
        }
    }
    float summrsel=0;
    float summrcrs=0;
    float summrmut=0;
    if(PRS.SADJ == 2)
    {
        for(int i=0;i!=PRS.NumSel;i++)
        {
            summrsel+=sucsel[i]*sucsel[i]/usedsel[i];
        }
        for(int i=0;i!=PRS.NumCrs;i++)
        {
            summrcrs+=succrs[i]*succrs[i]/usedcrs[i];
        }
        for(int i=0;i!=PRS.NumMut;i++)
        {
            summrmut+=sucmut[i]*sucmut[i]/usedmut[i];
        }
        for(int i=0;i!=PRS.NumSel;i++)
        {
            psel[i] = PRS.pselall + sucsel[i]*sucsel[i]/usedsel[i]*(1.-PRS.NumSel*PRS.pselall)/summrsel;
        }
        for(int i=0;i!=PRS.NumCrs;i++)
        {
            pcrs[i] = PRS.pcrsall + succrs[i]*succrs[i]/usedcrs[i]*(1.-PRS.NumCrs*PRS.pcrsall)/summrcrs;
        }
        for(int i=0;i!=PRS.NumMut;i++)
        {
            pmut[i] = PRS.pmutall + sucmut[i]*sucmut[i]/usedmut[i]*(1.-PRS.NumMut*PRS.pmutall)/summrmut;
        }
    }
    //foutSA<<psel[0]<<"\t"<<psel[1]<<"\t"<<psel[2]<<"\t"<<pcrs[0]<<"\t"<<pcrs[1]<<"\t"<<pcrs[2]<<"\t"<<pmut[0]<<"\t"<<pmut[1]<<"\t"<<pmut[2]<<"\t"<<pheu[0]<<"\t"<<pheu[1]<<endl;
    //FindNSaveBest(1,PRS,NumMass,FitMassCopy,Popul,PopulClass,PopulMu,PopulCF,RUPD,AveMass,ErrMass,FitMass,CurNumRules,
      //            BestInd,BestClasses,BestCF,BestMu,NumInst,CurRLength,UpdateBest,OverBestInd,OverBestClasses,OverBestCF,
        //          OverBestError,OverBestInd2,OverBestClasses2,OverBestCF2,Samp,OverBestMu,ClassConf,
          //        CVLearnSize,MRS,OverBestAve,ErrorClass,OverErrorClass,FoldOnTest,NumClassInst);
    FindNSaveBest(1,PRS,NumMass,FitMassCopy,Popul,PopulClass,PopulMu,PopulCF,RUPD,AveMass,ErrMass,FitMass,CurNumRules,
                  BestInd,BestClasses,BestCF,BestMu,NumInst,CurRLength,UpdateBest,OverBestInd,OverBestClasses,OverBestCF,
                  OverBestError,OverBestInd2,OverBestClasses2,OverBestCF2,Samp,OverBestMu,ClassConf,CVLearnSize,MRS,
                  OverBestAve,OverErrorClass,FoldOnTest,NumClassInst,ErrorClass,ConfMatrix,Inst,InstUsed,InstClass,
                  counterClass);

    //if(Gen%(int)(PRS.NumGens/1000.) == 0)
    {
        //TestErrCalc();
        cout<<Gen<<"\t"<<PRS.BestFit<<"\t"<<PRS.BestErr<<"\t"<<PRS.BestErr/NumInst<<"\t"<<PRS.BestNumRules<<"\t"<<PRS.BestRLength<<"\t"<<PRS.BestAve<<"\t";
        //cout<<endl;
        //FitCalc(0);
    }
}
