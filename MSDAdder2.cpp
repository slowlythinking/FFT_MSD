//
//  MSD_Adder.cpp
//

#include <fstream>
#include <sstream>
#include <iostream>
#include <map>
#include <cmath>
#include <vector>
#include <iomanip>
using namespace std;
const int DFT_N=32;//N点：输入向量为N维
const int MSD_m=32;//m位：输入向量每一维为m位
const int MSD_p=24;//小数部分位数（MSD precision）
static unsigned long G_YEJING = 0;//
static unsigned long G_TIME = 0;//
static int FLAGforYEJING= 0;//
static int FLAGforTIME= 0;//
static int FLAGforMultiYEJING= 0;//
static int FLAGforMultiTIME= 0;//


typedef string TERNARY_DATA; //用0、1和u符号串表示的三值数据
typedef map<char, map<char, char> > TERNARY_LOGICAL_TRUTH_TABLE; //存放三值逻辑变换真值表的二维映射表


//模拟三值逻辑运算器功能的函数
TERNARY_DATA TernaryLogicalOperator(TERNARY_DATA operandA, TERNARY_DATA operandB, TERNARY_LOGICAL_TRUTH_TABLE &truthTable)
{
    //操作数位数对齐(位数小的操作数补零)
    if(operandA.length() > operandB.length()){
	TERNARY_DATA zeroPadding(operandA.length() - operandB.length(), '0');
	operandB = zeroPadding + operandB;
    }else if(operandB.length() > operandA.length()){
	TERNARY_DATA zeroPadding(operandB.length() - operandA.length(), '0');
	operandA = zeroPadding + operandA;
    }

    //按位进行三值逻辑变换
    TERNARY_DATA resultC(operandA.length(), '0');

    for(unsigned long uIndex = 0; uIndex < operandA.length(); ++uIndex){
	resultC[uIndex] = truthTable[operandA[uIndex]][operandB[uIndex]];
    }

    return resultC;
}


//模拟MSD加法器功能的函数
TERNARY_DATA MSDAdder(TERNARY_DATA operandA, TERNARY_DATA operandB)
{
    //三值逻辑T变换真值表
    TERNARY_LOGICAL_TRUTH_TABLE truthTableT;
    truthTableT['u']['u']='u';
    truthTableT['u']['0']='u';
    truthTableT['u']['1']='0';
    truthTableT['0']['u']='u';
    truthTableT['0']['0']='0';
    truthTableT['0']['1']='1';
    truthTableT['1']['u']='0';
    truthTableT['1']['0']='1';
    truthTableT['1']['1']='1';

    //三值逻辑W变换真值表
    TERNARY_LOGICAL_TRUTH_TABLE truthTableW;
    truthTableW['u']['u']='0';
    truthTableW['u']['0']='1';
    truthTableW['u']['1']='0';
    truthTableW['0']['u']='1';
    truthTableW['0']['0']='0';
    truthTableW['0']['1']='u';
    truthTableW['1']['u']='0';
    truthTableW['1']['0']='u';
    truthTableW['1']['1']='0';

    //三值逻辑T'变换真值表
    TERNARY_LOGICAL_TRUTH_TABLE truthTableT1;
    truthTableT1['u']['u']='u';
    truthTableT1['u']['0']='0';
    truthTableT1['u']['1']='0';
    truthTableT1['0']['u']='0';
    truthTableT1['0']['0']='0';
    truthTableT1['0']['1']='0';
    truthTableT1['1']['u']='0';
    truthTableT1['1']['0']='0';
    truthTableT1['1']['1']='1';

    //三值逻辑W'变换真值表
    TERNARY_LOGICAL_TRUTH_TABLE truthTableW1;
    truthTableW1['u']['u']='0';
    truthTableW1['u']['0']='u';
    truthTableW1['u']['1']='0';
    truthTableW1['0']['u']='u';
    truthTableW1['0']['0']='0';
    truthTableW1['0']['1']='1';
    truthTableW1['1']['u']='0';
    truthTableW1['1']['0']='1';
    truthTableW1['1']['1']='0';

    //操作数位数对齐(位数小的操作数补零)
    if(operandA.length() > operandB.length()){
	TERNARY_DATA zeroPadding(operandA.length() - operandB.length(), '0');
	operandB = zeroPadding + operandB;
    }else if(operandB.length() > operandA.length()){
	TERNARY_DATA zeroPadding(operandB.length() - operandA.length(), '0');
	operandA = zeroPadding + operandA;
    }
    //执行三值逻辑T变换，结果补零
    TERNARY_DATA resultT = TernaryLogicalOperator(operandA, operandB, truthTableT) + '0';
    //执行三制逻辑W变换，结果补零
    TERNARY_DATA resultW = '0' + TernaryLogicalOperator(operandA, operandB, truthTableW);
    //执行三值逻辑T'变换，结果补零
    TERNARY_DATA resultT1 = TernaryLogicalOperator(resultT, resultW, truthTableT1) + '0';
    //执行三值逻辑W'变换，结果补零
    TERNARY_DATA resultW1 = '0' + TernaryLogicalOperator(resultT, resultW, truthTableW1);
    //三值逻辑T2变换真值表
    TERNARY_LOGICAL_TRUTH_TABLE& truthTableT2 = truthTableT;
    //执行三值逻辑T2变换，结果返回
    TERNARY_DATA resultC = TernaryLogicalOperator(resultT1, resultW1, truthTableT2);

    if (FLAGforYEJING)
//	G_YEJING += (unsigned long)(resultC.length());
	G_YEJING += 5*(unsigned long)(resultC.length());

    return resultC;
}

//MSD乘法实现
TERNARY_DATA MultiplicationRoutine(TERNARY_DATA operandA, TERNARY_DATA operandB){
    
    unsigned long PartialNum = operandB.length();
    unsigned long AddTreeLayers = ceil(log2(PartialNum));
    
    //复制乘数B的每一位到QB数据队列
    TERNARY_DATA QueueB[100];
    TERNARY_DATA temp;
    for(unsigned long uIndex = 0; uIndex < operandB.length(); ++uIndex)
	QueueB[uIndex].insert(0, operandA.length(), operandB[operandB.length() - 1 - uIndex]);

    //三值逻辑M变换真值表
    TERNARY_LOGICAL_TRUTH_TABLE truthTableM;
    truthTableM['u']['u']='1';
    truthTableM['u']['0']='0';
    truthTableM['u']['1']='u';
    truthTableM['0']['u']='0';
    truthTableM['0']['0']='0';
    truthTableM['0']['1']='0';
    truthTableM['1']['u']='u';
    truthTableM['1']['0']='0';
    truthTableM['1']['1']='1';

    //将乘数A与QB数据队列中每一项进行三值逻辑M变换，结果补零存入QS数据队列
    TERNARY_DATA QueueS[100];
    for(unsigned long uIndex = 0; uIndex < operandB.length(); ++uIndex){ //将A和QB[i]送入M变换器进行运算，结果补零为和数项，存入QS数据队列
	QueueS[uIndex] = TernaryLogicalOperator(operandA, QueueB[uIndex], truthTableM);
	QueueS[uIndex].insert(QueueS[uIndex].length(), uIndex, '0');
    }
    if (FLAGforTIME == 1)
	G_TIME += 1;


    //QS数据队列中的和数项或部分和两两相加迭代，结果存入QH数据队列，迭代结束返回最终乘积
    TERNARY_DATA *QueueH = QueueS;
    for(unsigned long turnK = 1; turnK <= ceil(log2(operandB.length())); ++turnK){
	if (FLAGforTIME == 1)
	    G_TIME += 3;
	if (FLAGforMultiYEJING)
	    FLAGforYEJING = 1;
	//求和项个数为奇数时，增加一个值为0的偶数项
	if(1 == (unsigned long)ceil(operandB.length() / pow(2, turnK - 1)) % 2)
	    QueueH[(unsigned long)ceil(operandB.length() / pow(2, turnK - 1))] = "0";

	//第k轮两两相加迭代
	for(unsigned long uIndex = 0; uIndex < ceil(operandB.length() / pow(2, turnK)); ++uIndex){
	    QueueH[uIndex] = MSDAdder(QueueH[2 * uIndex], QueueH[2 * uIndex + 1]);
	}
	FLAGforYEJING = 0;
    }
    //
    temp = QueueH[0].substr(0,QueueH[0].length()-MSD_p);
   
//    return QueueH[0];
    return temp;
}

//将三角函数sin和cos的double结果转为MSD数
TERNARY_DATA ConvertDoubleToMSD(double x)
{
    TERNARY_DATA msd;
    double number = (double)x;
    double decimal;
    int integer;
    TERNARY_DATA temp(MSD_p,'0');

    integer = (int)number;
    decimal = number - (double)integer;

    if (integer == 1) 
	msd = "1" + temp;
    else if (integer == -1) 
	msd = "u" + temp;
    else
    {
	msd = "0";
	int i=0;
	int count=0;
	while (count < MSD_p) //小数部分取MSD_p位
	{
	    decimal*=2;
	    if ((int)decimal == 1)
		msd += '1';
	    else if ((int)decimal == -1)
		msd += 'u';
	    else 
		msd += '0';
	    if (decimal>=1.0)
		decimal-=1.0;
	    else if (decimal <= -1.0)
		decimal+=1.0;
	    ++count;
	}
    }
    return msd;
}

//将整数MSD数转化为十进制数
long ConTernaryNumberInt(TERNARY_DATA a)
{
    long temp = 0;
    int index = 0;
    for(int iThBit =(a.length() - 1);iThBit != -1;--iThBit)
    {
	int i;
	if (a[iThBit] == '0')
	    i = 0;
	else if (a[iThBit] == '1')
	    i = 1;
	else
	    i = -1;
	temp += i*pow(2,index);
	++index;
    }
}
//将含有MSD_p位小数的MSD数转化为十进制数
double ConTernaryNumberdouble(string a)
{
    double temp = 0.0;
    int index1 = -1;
    int index2 = 0;
    for(int iThBit =(a.length() - MSD_p);iThBit != a.length();++iThBit)
    {
	int i;
	if (a[iThBit] == '0')
	    i = 0;
	else if (a[iThBit] == '1')
	    i = 1;
	else
	    i = -1;
	temp += i*pow(2,index1);
	--index1;
    }
    for(int iThBit =(a.length() - MSD_p - 1);iThBit != -1;--iThBit)
    {
	int i;
	if (a[iThBit] == '0')
	    i = 0;
	else if (a[iThBit] == '1')
	    i = 1;
	else
	    i = -1;
	temp += i*pow(2,index2);
	++index2;
    }
    return temp;
}

//求一个MSD数的相反数
TERNARY_DATA MSDinverse(TERNARY_DATA temp)
{
    TERNARY_DATA tempResult;
    for(int iThBit = 0;iThBit != temp.length();++iThBit)
    {
	if (temp[iThBit] == '0')
	    tempResult += '0';
	else if (temp[iThBit] == '1')
	    tempResult += 'u';
	else
	    tempResult += '1';
    }
    return tempResult;
}

//DFT十进制实现
void ConDFTDec(double (&dataInreal)[DFT_N],double (&dataInimag)[DFT_N],double (&dataOutreal)[DFT_N],double (&dataOutimag)[DFT_N])
{
    for (unsigned int k = 0; k != DFT_N; k++)
    {  /* For each output element */
	cout << "DFT-decimble computing the " << k << "th data"<< endl;
	double sumreal = 0;
	double sumimag = 0;
	for (unsigned int t = 0; t != DFT_N; t++) {  /* For each input element */
	    double angle = 2 * M_PI * t * k / DFT_N;
	    sumreal +=  dataInreal[t] * cos(angle) + dataInimag[t] * sin(angle);
	    sumimag += -dataInreal[t] * sin(angle) + dataInimag[t] * cos(angle);
	}
	dataOutreal[k] = sumreal;
	dataOutimag[k] = sumimag;
    }
}
//DFT MSD实现
void ComDFTMSD(TERNARY_DATA (&dataInreal)[DFT_N],TERNARY_DATA (&dataInimag)[DFT_N],TERNARY_DATA (&dataOutreal)[DFT_N],TERNARY_DATA (&dataOutimag)[DFT_N])
{
    FLAGforTIME = 1;
    for (unsigned int k = 0; k != DFT_N; k++) {  /* For each output element */
	cout << "DFT-MSD computing the " << k << "th data"<< endl;
	TERNARY_DATA sumreal;
	TERNARY_DATA sumimag;
	int tempFlag = 1;

	//every dimension
	for (unsigned int t = 0; t != DFT_N; t++) {  /* For each input element */
	    double angle = 2 * M_PI * t * k / DFT_N; 
	    TERNARY_DATA tempreal,tempreal1,tempreal2;
	    TERNARY_DATA tempimag,tempimag1,tempimag2;
	    TERNARY_DATA cosMSD = ConvertDoubleToMSD(cos(angle));
	    TERNARY_DATA sinMSD = ConvertDoubleToMSD(sin(angle));

	    //parallel computing 
	    //1.tempreal and tempimag can be parallel computed 
	    //2.MultiplicationRoutine(dataInreal[t],cosMSD) and MultiplicationRoutine(dataInimag[t],sinMSD) can be parallel computed 

	    //计算矩阵一项与输入向量一项做复数乘法后结果实部的一部分
	    FLAGforMultiYEJING = 1;
	    tempreal1 = MultiplicationRoutine(dataInreal[t],cosMSD);
	    FLAGforTIME = 0;
	    tempreal2 = MultiplicationRoutine(dataInimag[t],sinMSD);
	    tempreal = MSDAdder(tempreal1,tempreal2);


	    //计算矩阵一项与输入向量一项做复数乘法后结果虚部的一部分
	    tempimag1 = MultiplicationRoutine(MSDinverse(dataInreal[t]),sinMSD);
	    tempimag2 = MultiplicationRoutine(dataInimag[t],cosMSD);
	    tempimag = MSDAdder(tempimag1,tempimag2);

	    sumreal = MSDAdder(sumreal,tempreal);
	    sumimag = MSDAdder(sumimag,tempimag);
	}
	dataOutreal[k] = sumreal;
	dataOutimag[k] = sumimag;
    }
    G_TIME += (unsigned long)(ceil(log2(DFT_N))*3);
}

//基4FFT-MSD实现
void ComRadix_4_FFT_MSD(TERNARY_DATA (&dataInreal)[DFT_N],TERNARY_DATA (&dataInimag)[DFT_N],TERNARY_DATA (&dataOutreal)[DFT_N],TERNARY_DATA (&dataOutimag)[DFT_N])
{
    FLAGforTIME = 1;

    //输入N维数据（实部或虚部）转化为四个N/4维数组（实部或虚部）
    int inputIndex[] = { 0, 2, 1, 3 };
    TERNARY_DATA inputReal[4][DFT_N/4]; 
    TERNARY_DATA inputImag[4][DFT_N/4]; 

    //基4FFT计算的中间结果
    TERNARY_DATA F_Real[4][DFT_N/4]; 
    TERNARY_DATA F_Imag[4][DFT_N/4]; 

    //需要用到的由旋转矩阵构成的参数矩阵
    TERNARY_DATA TwiddleFactorMatrixReal[4][DFT_N/4][DFT_N/4] ; 
    TERNARY_DATA TwiddleFactorMatrixImag[4][DFT_N/4][DFT_N/4] ; 

    //初始化，准备数据
    for (unsigned int k = 0; k != DFT_N/4; k++) {  /* For each output element */
	//把输入数据分为四组,0到3依次存放X1,X3,X2,X4
	for(int i = 0;i != 4;++i)
	{
	    inputReal[i][k] = dataInreal[inputIndex[i]];
	    inputImag[i][k] = dataInimag[inputIndex[i]];
	    inputIndex[i] += 4; 
	}

	//计算相关旋转因子组成的参数矩阵的值(顺序为X1,X3*W2k,X2*Wk,X4*W3k)，八个矩阵，四个实部矩阵，四个虚部矩阵。
	for (unsigned int t = 0; t != DFT_N/4; t++) {  /* For each input element */
	    double angle[4];
	    angle[0] = 2 * M_PI * 4 * t * k / DFT_N; 
	    angle[1] = 2 * M_PI * (4 * t * k + 2*k) / DFT_N; 
	    angle[2] = 2 * M_PI * (4 * t * k + k) / DFT_N; 
	    angle[3] = 2 * M_PI * (4 * t * k + 3*k) / DFT_N; 

	    for(int j = 0;j != 4;++j)
	    {
		TERNARY_DATA cosMSD = ConvertDoubleToMSD(cos(angle[j]));
		TERNARY_DATA sinMSD = ConvertDoubleToMSD(sin(angle[j]));
		TwiddleFactorMatrixReal[j][k][t] = cosMSD;
		TwiddleFactorMatrixImag[j][k][t] = MSDinverse(sinMSD);
	    }
	}
    }

    //计算中间结果，即文中F1,F2，F3,F4四个复数矩阵
    for (int j = 0;j != 4;++j)
    {
	cout << "Radix4FFT-MSD computing the " << j << " F matix"<< endl;

	for (unsigned int k = 0; k != DFT_N/4; k++) {  /* For each output element */
	    TERNARY_DATA sumreal;
	    TERNARY_DATA sumimag;

	    //every dimension
	    for (unsigned int t = 0; t != DFT_N/4; t++) {  /* For each input element */

		TERNARY_DATA tempreal,tempreal1,tempreal2;
		TERNARY_DATA tempimag,tempimag1,tempimag2;
		//计算矩阵一项与输入向量一项做复数乘法后结果实部的一部分
		tempreal1 = MultiplicationRoutine(inputReal[j][t],TwiddleFactorMatrixReal[j][k][t]);
		tempreal2 = MultiplicationRoutine(inputImag[j][t],TwiddleFactorMatrixImag[j][k][t]);
		tempreal = MSDAdder(tempreal1,MSDinverse(tempreal2));


		//计算矩阵一项与输入向量一项做复数乘法后结果虚部的一部分
		tempimag1 = MultiplicationRoutine(inputReal[j][t],TwiddleFactorMatrixImag[j][k][t]);
		tempimag2 = MultiplicationRoutine(inputImag[j][t],TwiddleFactorMatrixReal[j][k][t]);
		tempimag = MSDAdder(tempimag1,tempimag2);

		sumreal = MSDAdder(sumreal,tempreal);
		sumimag = MSDAdder(sumimag,tempimag);
	    }
	F_Real[j][k] = sumreal;
	F_Imag[j][k] = sumimag;
	}
    }

    //计算最终结果
   for (int i = 0;i != DFT_N/4;++i)
   {
        dataOutreal[i] = MSDAdder( MSDAdder(F_Real[0][i],F_Real[1][i]) , MSDAdder(F_Real[2][i],F_Real[3][i]) );
        dataOutimag[i] =MSDAdder( MSDAdder(F_Imag[0][i],F_Imag[1][i]) , MSDAdder(F_Imag[2][i],F_Imag[3][i]) );//F_Imag[0][i] + F_Imag[1][i] + F_Imag[2][i] + F_Imag[3][i];
   }
   for (int i = DFT_N/4, j = 0;i != DFT_N/2;++i,++j)
   {
        dataOutreal[i] = MSDAdder( MSDAdder(F_Real[0][j],MSDinverse(F_Real[1][j])) , MSDAdder(F_Imag[2][j],MSDinverse(F_Imag[3][j])) );//F_Real[0][i] - F_Real[1][i] - F_Imag[2][i] + F_Imag[3][i];
        dataOutimag[i] = MSDAdder(MSDAdder(F_Imag[0][j],MSDinverse(F_Imag[1][j])),MSDAdder(MSDinverse(F_Real[2][j]),F_Real[3][j]));//F_Imag[0][i] - F_Imag[1][i] - F_Real[2][i] + F_Real[3][i];
   }
   for (int i = DFT_N/2, j = 0;i != 3*DFT_N/4;++i,++j)
   {
        dataOutreal[i] = MSDAdder( MSDAdder(F_Real[0][j],F_Real[1][j]) , MSDAdder(MSDinverse(F_Real[2][j]),MSDinverse(F_Real[3][j])));//F_Real[0][i] + F_Real[1][i] - F_Real[2][i] - F_Real[3][i];
        dataOutimag[i] = MSDAdder(MSDAdder(F_Imag[0][j],F_Imag[1][j]),MSDAdder(MSDinverse(F_Imag[2][j]),MSDinverse(F_Imag[3][j])));//F_Imag[0][i] + F_Imag[1][i] - F_Imag[2][i] - F_Imag[3][i];
   }
   for (int i = 3*DFT_N/4, j = 0;i != DFT_N;++i,++j)
   {
        dataOutreal[i] = MSDAdder(MSDAdder(F_Real[0][j],MSDinverse(F_Real[1][j])),MSDAdder(MSDinverse(F_Imag[2][j]),F_Imag[3][j]));//F_Real[0][i] - F_Real[1][i] + F_Imag[2][i] - F_Imag[3][i];
        dataOutimag[i] = MSDAdder(MSDAdder(F_Imag[0][j],MSDinverse(F_Imag[1][j])),MSDAdder(F_Real[2][j],MSDinverse(F_Real[3][j])));//F_Imag[0][i] - F_Imag[1][i] + F_Real[2][i] - F_Real[3][i];
   }
}

//基8FFT-MSD实现
TERNARY_DATA MSDADD8(TERNARY_DATA t1,TERNARY_DATA t2,TERNARY_DATA t3,TERNARY_DATA t4,TERNARY_DATA t5,TERNARY_DATA t6,TERNARY_DATA t7,TERNARY_DATA t8)
{
    return MSDAdder(MSDAdder(MSDAdder(t1,t2),MSDAdder(t3,t4)),MSDAdder(MSDAdder(t5,t6),MSDAdder(t7,t8)));
}


void ComRadix_8_FFT_MSD(TERNARY_DATA (&dataInreal)[DFT_N],TERNARY_DATA (&dataInimag)[DFT_N],TERNARY_DATA (&dataOutreal)[DFT_N],TERNARY_DATA (&dataOutimag)[DFT_N])
{
    FLAGforTIME = 1;

    //输入N维数据（实部或虚部）转化为八个N/8维数组（实部或虚部）
    int inputIndex[] = { 0, 1, 2, 3, 4, 5, 6, 7 };
    TERNARY_DATA inputReal[8][DFT_N/8]; 
    TERNARY_DATA inputImag[8][DFT_N/8]; 

    //基8FFT计算的中间结果
    TERNARY_DATA First_Real[16][DFT_N/8]; 
    TERNARY_DATA First_Imag[16][DFT_N/8]; 
//
//    TERNARY_DATA Second_Real[8][DFT_N/8]; 
//    TERNARY_DATA Second_Imag[8][DFT_N/8]; 
//
//    TERNARY_DATA Second_Real[8][DFT_N/8]; 
//    TERNARY_DATA Second_Imag[8][DFT_N/8]; 
    //需要用到的由旋转矩阵构成的参数矩阵
    TERNARY_DATA TwiddleFactorMatrixReal[16][DFT_N/8][DFT_N/8] ; 
    TERNARY_DATA TwiddleFactorMatrixImag[16][DFT_N/8][DFT_N/8] ; 

    //初始化，准备数据
    for (unsigned int k = 0; k != DFT_N/8; k++) {  /* For each output element */
	//把输入数据分为八组
	for(int i = 0;i != 8;++i)
	{
	    inputReal[i][k] = dataInreal[inputIndex[i]];
	    inputImag[i][k] = dataInimag[inputIndex[i]];
	    inputIndex[i] += 8; 
	}

	//计算相关旋转因子组成的参数矩阵的值(顺序为X1,X3*W2k,X2*Wk,X4*W3k)，二十四个矩阵，十六个实部矩阵，十六个虚部矩阵。
	for (unsigned int t = 0; t != DFT_N/8; t++) {  /* For each input element */
	    double angle[8];
	    for (int j = 0;j != 8;++j)
	    {
	    angle[j] = 2 * M_PI * (8 * t * k + j*k) / DFT_N; 
	    }
//	    angle[8] = 2 * M_PI * ((8 * t * k + k) / DFT_N + 1/8);
//	    angle[9] = 2 * M_PI * ((8 * t * k + k) / DFT_N + 3/8);
//
//	    angle[10] = 2 * M_PI * ((8 * t * k +3*k) / DFT_N + 3/8);
//	    angle[11] = 2 * M_PI * ((8 * t * k +3*k) / DFT_N + 9/8);
//
//	    angle[12] = 2 * M_PI * ((8 * t * k +5*k) / DFT_N + 5/8);
//	    angle[13] = 2 * M_PI * ((8 * t * k +5*k) / DFT_N + 15/8);
//
//	    angle[14] = 2 * M_PI * ((8 * t * k +7*k) / DFT_N + 7/8);
//	    angle[15] = 2 * M_PI * ((8 * t * k +7*k) / DFT_N + 21/8);

//	    for(int j = 0;j != 16;++j)
	    for(int j = 0;j != 8;++j)
	    {
		TERNARY_DATA cosMSD = ConvertDoubleToMSD(cos(angle[j]));
		TERNARY_DATA sinMSD = ConvertDoubleToMSD(sin(angle[j]));
		TwiddleFactorMatrixReal[j][k][t] = cosMSD;
		TwiddleFactorMatrixImag[j][k][t] = MSDinverse(sinMSD);
	    }
	}
    }
    
    for (unsigned int k = 0; k != DFT_N/8; k++) {  /* For each output element */
	for (unsigned int t = 0; t != DFT_N/8; t++) {
	    double angle8[16];
	    angle8[8] = 2*M_PI*((k+DFT_N/8)*(8*t+1))/DFT_N;
	    angle8[9] = 2*M_PI*((k+3*DFT_N/8)*(8*t+1))/DFT_N;
	    angle8[10] = 2*M_PI*((k+DFT_N/8)*(8*t+3))/DFT_N;
	    angle8[11] = 2*M_PI*((k+3*DFT_N/8)*(8*t+3))/DFT_N;
	    angle8[12] = 2*M_PI*((k+DFT_N/8)*(8*t+5))/DFT_N;
	    angle8[13] = 2*M_PI*((k+3*DFT_N/8)*(8*t+5))/DFT_N;
	    angle8[14] = 2*M_PI*((k+DFT_N/8)*(8*t+7))/DFT_N;
	    angle8[15] = 2*M_PI*((k+3*DFT_N/8)*(8*t+7))/DFT_N;
	    for(int j = 8;j != 16;++j)
	    {
		TERNARY_DATA cosMSD = ConvertDoubleToMSD(cos(angle8[j]));
		TERNARY_DATA sinMSD = ConvertDoubleToMSD(sin(angle8[j]));
		TwiddleFactorMatrixReal[j][k][t] = cosMSD;
		TwiddleFactorMatrixImag[j][k][t] = MSDinverse(sinMSD);
	    }
	}
    }

    //计算中间结果，即十六个复数参数矩阵与输入的八个八维向量数据的乘积
    for (int j = 0;j != 16;++j)
    {
	cout << "Radix8FFT-MSD computing the " << j << " F matix"<< endl;
	int tempj;
	if(j < 8)
	    tempj = j;
	else if (j == 8 || j == 9)
	    tempj = 1;
	else if (j == 10 || j == 11)
	    tempj = 3;
	else if (j == 12 || j == 13)
	    tempj = 5;
	else
	    tempj = 7;

	for (unsigned int k = 0; k != DFT_N/8; k++) {  /* For each output element */
	    TERNARY_DATA sumreal;
	    TERNARY_DATA sumimag;

	    //every dimension
	    for (unsigned int t = 0; t != DFT_N/8; t++) {  /* For each input element */

		TERNARY_DATA tempreal,tempreal1,tempreal2;
		TERNARY_DATA tempimag,tempimag1,tempimag2;
		//计算矩阵一项与输入向量一项做复数乘法后结果实部的一部分
		tempreal1 = MultiplicationRoutine(inputReal[tempj][t],TwiddleFactorMatrixReal[j][k][t]);
		tempreal2 = MultiplicationRoutine(inputImag[tempj][t],TwiddleFactorMatrixImag[j][k][t]);
		tempreal = MSDAdder(tempreal1,MSDinverse(tempreal2));


		//计算矩阵一项与输入向量一项做复数乘法后结果虚部的一部分
		tempimag1 = MultiplicationRoutine(inputReal[tempj][t],TwiddleFactorMatrixImag[j][k][t]);
		tempimag2 = MultiplicationRoutine(inputImag[tempj][t],TwiddleFactorMatrixReal[j][k][t]);
		tempimag = MSDAdder(tempimag1,tempimag2);

		sumreal = MSDAdder(sumreal,tempreal);
		sumimag = MSDAdder(sumimag,tempimag);
	    }
	First_Real[j][k] = sumreal;
	First_Imag[j][k] = sumimag;
	}
    }


    //计算最终结果
    //1-0到N/8维 0+4+2+6+1+3+5+7
    cout << "第1部分" << endl;
    for (int i = 0;i != DFT_N/8;++i)
    {
	dataOutreal[i] = MSDADD8(First_Real[0][i],First_Real[4][i],First_Real[2][i],First_Real[6][i],First_Real[1][i],First_Real[3][i],First_Real[5][i],First_Real[7][i]);
	dataOutimag[i] = MSDADD8(First_Imag[0][i],First_Imag[4][i],First_Imag[2][i],First_Imag[6][i],First_Imag[1][i],First_Imag[3][i],First_Imag[5][i],First_Imag[7][i]);
	cout << "第" << i << "维，dataOutreal[i]:" << ConTernaryNumberdouble(dataOutreal[i]) << "------------dataOutImag[i]:" << ConTernaryNumberdouble(dataOutimag[i]) << endl;
    }
	cout << endl;
    //2-N/8维到2N/8维
    cout << "第2部分" << endl;
    for (int i = DFT_N/8, j = 0;i != DFT_N/4;++i,++j)
   {
	dataOutreal[i] = MSDADD8(First_Real[0][j],MSDinverse(First_Real[4][j]),First_Imag[2][j],MSDinverse(First_Imag[6][j]),First_Real[8][j],First_Real[10][j],First_Real[12][j],First_Real[14][j]);
	dataOutimag[i] = MSDADD8(First_Imag[0][j],MSDinverse(First_Imag[4][j]),MSDinverse(First_Real[2][j]),First_Real[6][j],First_Imag[8][j],First_Imag[10][j],First_Imag[12][j],First_Imag[14][j]);
	cout << "第" << i << "维，dataOutreal[i]:" << ConTernaryNumberdouble(dataOutreal[i]) << "------------dataOutImag[i]:" << ConTernaryNumberdouble(dataOutimag[i]) << endl;
   }
	cout << endl;
    //3-2N/8维到3N/8维 0+4-2-6-j(1-3+5-7)
    cout << "第3部分" << endl;
   for (int i = DFT_N/4, j = 0;i != 3*DFT_N/8;++i,++j)
   {
	dataOutreal[i] = MSDADD8(First_Real[0][j],First_Real[4][j],MSDinverse(First_Real[2][j]),MSDinverse(First_Real[6][j]),First_Imag[1][j],MSDinverse(First_Imag[3][j]),First_Imag[5][j],MSDinverse(First_Imag[7][j]));
	dataOutimag[i] = MSDADD8(First_Imag[0][j],First_Imag[4][j],MSDinverse(First_Imag[2][j]),MSDinverse(First_Imag[6][j]),MSDinverse(First_Real[1][j]),First_Real[3][j],MSDinverse(First_Real[5][j]),First_Real[7][j]);
	cout << "第" << i << "维，dataOutreal[i]:" << ConTernaryNumberdouble(dataOutreal[i]) << "------------dataOutImag[i]:" << ConTernaryNumberdouble(dataOutimag[i]) << endl;
   }
	cout << endl;
    //4-3N/8维到4N/8维
    cout << "第4部分" << endl;
   for (int i = 3*DFT_N/8, j = 0;i != DFT_N/2;++i,++j)
   {
dataOutreal[i] = MSDADD8(First_Real[0][j],MSDinverse(First_Real[4][j]),MSDinverse(First_Imag[2][j]),First_Imag[6][j],First_Real[9][j],First_Real[11][j],First_Real[13][j],First_Real[15][j]);
dataOutimag[i] = MSDADD8(First_Imag[0][j],MSDinverse(First_Imag[4][j]),First_Real[2][j],MSDinverse(First_Real[6][j]),First_Imag[9][j],First_Imag[11][j],First_Imag[13][j],First_Imag[15][j]);  
	cout << "第" << i << "维，dataOutreal[i]:" << ConTernaryNumberdouble(dataOutreal[i]) << "------------dataOutImag[i]:" << ConTernaryNumberdouble(dataOutimag[i]) << endl;
    }
	cout << endl;
    //5-4N/8维到5N/8维 0+4+2+6-(1+3+5+7)
    cout << "第5部分" << endl;
   for (int i = DFT_N/2, j = 0;i !=5*DFT_N/8;++i,++j)
   {
	dataOutreal[i] = MSDADD8(First_Real[0][j],First_Real[4][j],First_Real[2][j],First_Real[6][j],MSDinverse(First_Real[1][j]),MSDinverse(First_Real[3][j]),MSDinverse(First_Real[5][j]),MSDinverse(First_Real[7][j]));
	dataOutimag[i] = MSDADD8(First_Imag[0][j],First_Imag[4][j],First_Imag[2][j],First_Imag[6][j],MSDinverse(First_Imag[1][j]),MSDinverse(First_Imag[3][j]),MSDinverse(First_Imag[5][j]),MSDinverse(First_Imag[7][j]));
	cout << "第" << i << "维，dataOutreal[i]:" << ConTernaryNumberdouble(dataOutreal[i]) << "------------dataOutImag[i]:" << ConTernaryNumberdouble(dataOutimag[i]) << endl;
   }
	cout << endl;
    //6-5N/8维到6N/8维
    cout << "第6部分" << endl;
   for (int i = 5*DFT_N/8, j = 0;i != 3*DFT_N/4;++i,++j)
   {
	dataOutreal[i] = MSDADD8(First_Real[0][j],MSDinverse(First_Real[4][j]),First_Imag[2][j],MSDinverse(First_Imag[6][j]),MSDinverse(First_Real[8][j]),MSDinverse(First_Real[10][j]),MSDinverse(First_Real[12][j]),MSDinverse(First_Real[14][j]));
	dataOutimag[i] = MSDADD8(First_Imag[0][j],MSDinverse(First_Imag[4][j]),MSDinverse(First_Real[2][j]),First_Real[6][j],MSDinverse(First_Imag[8][j]),MSDinverse(First_Imag[10][j]),MSDinverse(First_Imag[12][j]),MSDinverse(First_Imag[14][j]));
	cout << "第" << i << "维，dataOutreal[i]:" << ConTernaryNumberdouble(dataOutreal[i]) << "------------dataOutImag[i]:" << ConTernaryNumberdouble(dataOutimag[i]) << endl;
   }
	cout << endl;
    //7-6N/8维到7N/8维 0+4-(2+6)+j(1-3+5-7)
    cout << "第7部分" << endl;
   for (int i = 3*DFT_N/4, j = 0;i != 7*DFT_N/8;++i,++j)
   {
	dataOutreal[i] = MSDADD8(First_Real[0][j],First_Real[4][j],MSDinverse(First_Real[2][j]),MSDinverse(First_Real[6][j]),MSDinverse(First_Imag[1][j]),First_Imag[3][j],MSDinverse(First_Imag[5][j]),First_Imag[7][j]);
	dataOutimag[i] = MSDADD8(First_Imag[0][j],First_Imag[4][j],MSDinverse(First_Imag[2][j]),MSDinverse(First_Imag[6][j]),First_Real[1][j],MSDinverse(First_Real[3][j]),First_Real[5][j],MSDinverse(First_Real[7][j]));
	cout << "第" << i << "维，dataOutreal[i]:" << ConTernaryNumberdouble(dataOutreal[i]) << "------------dataOutImag[i]:" << ConTernaryNumberdouble(dataOutimag[i]) << endl;
   }
	cout << endl;
   //8-7N/8维到第N维
    cout << "第8部分" << endl;
   for (int i = 7*DFT_N/8, j = 0;i != DFT_N;++i,++j)
   {
       dataOutreal[i] = MSDADD8(First_Real[0][j],MSDinverse(First_Real[4][j]),MSDinverse(First_Imag[2][j]),First_Imag[6][j],MSDinverse(First_Real[9][j]),MSDinverse(First_Real[11][j]),MSDinverse(First_Real[13][j]),MSDinverse(First_Real[15][j]));
       dataOutimag[i] = MSDADD8(First_Imag[0][j],MSDinverse(First_Imag[4][j]),First_Real[2][j],MSDinverse(First_Real[6][j]),MSDinverse(First_Imag[9][j]),MSDinverse(First_Imag[11][j]),MSDinverse(First_Imag[13][j]),MSDinverse(First_Imag[15][j]));  
       cout << "第" << i << "维，dataOutreal[i]:" << ConTernaryNumberdouble(dataOutreal[i]) << "------------dataOutImag[i]:" << ConTernaryNumberdouble(dataOutimag[i]) << endl;
   }
}


// 电子计算机实现快速傅立叶变换 Fast Fourier Transform
// 改进了《算法导论》的算法，旋转因子取 ωn-kj  (ωnkj 的共轭复数)
// 且只计算 n / 2 次，而未改进前需要计算 (n * lg n) / 2 次。
inline void swap (double &a, double &b)
{
    double t;
    t = a;
    a = b;
    b = t;
}

void bitrp (double xreal [], double ximag [], int n)
{
    // 位反转置换 Bit-reversal Permutation
    int i, j, a, b, p;

    for (i = 1, p = 0; i < n; i *= 2)
        {
        p ++;
        }
    for (i = 0; i < n; i ++)
        {
        a = i;
        b = 0;
        for (j = 0; j < p; j ++)
            {
            b = (b << 1) + (a & 1);    // b = b * 2 + a % 2;
            a >>= 1;        // a = a / 2;
            }
        if ( b > i)
            {
            swap (xreal [i], xreal [b]);
            swap (ximag [i], ximag [b]);
            }
        }
}

//void FFT(double (&xreal)[DFT_N], double (&ximag)[DFT_N])
void FFT(double xreal[], double ximag[])
{
    // 快速傅立叶变换，将复数 x 变换后仍保存在 x 中，xreal, ximag 分别是 x 的实部和虚部
    double wreal [DFT_N / 2], wimag [DFT_N / 2], treal, timag, ureal, uimag, arg;
    int m, k, j, t, index1, index2;
    int n = DFT_N;
    
    bitrp (xreal, ximag, n);

    // 计算 1 的前 n / 2 个 n 次方根的共轭复数 W'j = wreal [j] + i * wimag [j] , j = 0, 1, ... , n / 2 - 1
    arg = - 2 * M_PI / n;
    treal = cos (arg);
    timag = sin (arg);
    wreal [0] = 1.0;
    wimag [0] = 0.0;
    for (j = 1; j < n / 2; j ++)
        {
        wreal [j] = wreal [j - 1] * treal - wimag [j - 1] * timag;
        wimag [j] = wreal [j - 1] * timag + wimag [j - 1] * treal;
        }

    for (m = 2; m <= n; m *= 2)
        {
        for (k = 0; k < n; k += m)
            {
            for (j = 0; j < m / 2; j ++)
                {
                index1 = k + j;
                index2 = index1 + m / 2;
                t = n * j / m;    // 旋转因子 w 的实部在 wreal [] 中的下标为 t
                treal = wreal [t] * xreal [index2] - wimag [t] * ximag [index2];
                timag = wreal [t] * ximag [index2] + wimag [t] * xreal [index2];
                ureal = xreal [index1];
                uimag = ximag [index1];
                xreal [index1] = ureal + treal;
                ximag [index1] = uimag + timag;
                xreal [index2] = ureal - treal;
                ximag [index2] = uimag - timag;
                }
            }
        }
}

/*
void  IFFT (double (&xreal)[DFT_N], double (&ximag)[DFT_N])
{
    // 快速傅立叶逆变换
    double wreal [DFT_N / 2], wimag [DFT_N / 2], treal, timag, ureal, uimag, arg;
    int m, k, j, t, index1, index2;
    int n = DFT_N;
    
    bitrp (xreal, ximag, n);

    // 计算 1 的前 n / 2 个 n 次方根 Wj = wreal [j] + i * wimag [j] , j = 0, 1, ... , n / 2 - 1
    arg = 2 * M_PI / n;
    treal = cos (arg);
    timag = sin (arg);
    wreal [0] = 1.0;
    wimag [0] = 0.0;
    for (j = 1; j < n / 2; j ++)
        {
        wreal [j] = wreal [j - 1] * treal - wimag [j - 1] * timag;
        wimag [j] = wreal [j - 1] * timag + wimag [j - 1] * treal;
        }

    for (m = 2; m <= n; m *= 2)
        {
        for (k = 0; k < n; k += m)
            {
            for (j = 0; j < m / 2; j ++)
                {
                index1 = k + j;
                index2 = index1 + m / 2;
                t = n * j / m;    // 旋转因子 w 的实部在 wreal [] 中的下标为 t
                treal = wreal [t] * xreal [index2] - wimag [t] * ximag [index2];
                timag = wreal [t] * ximag [index2] + wimag [t] * xreal [index2];
                ureal = xreal [index1];
                uimag = ximag [index1];
                xreal [index1] = ureal + treal;
                ximag [index1] = uimag + timag;
                xreal [index2] = ureal - treal;
                ximag [index2] = uimag - timag;
                }
            }
        }

    for (j=0; j < n; j ++)
        {
        xreal [j] /= n;
        ximag [j] /= n;
        }
}
*/

//主函数从实验用例文件读取加法操作数，运算结果输出到实验结果文件中
int main(int argc, const char * argv[])
{
//    cout << M_PI;
//    return 0;
    stringstream SDFT_N,SMSD_m,SMSD_p;
    SDFT_N << DFT_N;
    SMSD_m << MSD_m;
    SMSD_p << MSD_p;
    string INPUTFILE_s = "UseCase" + SDFT_N.str() + "_" + SMSD_m.str() + "_" + SMSD_p.str() + ".test";
    string OUTPUTFILE_s = "Result" + SDFT_N.str() + "_" + SMSD_m.str() + "_" + SMSD_p.str() + ".test";
    const char* INPUTFILE = INPUTFILE_s.c_str();
    const char* OUTPUTFILE = OUTPUTFILE_s.c_str();


    ifstream inputFile(INPUTFILE);
    ofstream outputFile(OUTPUTFILE);
    outputFile.precision(20);
    TERNARY_DATA InputReal[DFT_N];
    TERNARY_DATA InputImag[DFT_N];
    TERNARY_DATA OutputReal[DFT_N];
    TERNARY_DATA OutputImag[DFT_N];
    TERNARY_DATA Radix4OutputReal[DFT_N];
    TERNARY_DATA Radix4OutputImag[DFT_N];
    TERNARY_DATA Radix8OutputReal[DFT_N];
    TERNARY_DATA Radix8OutputImag[DFT_N];
    double inreal[DFT_N];
    double inimag[DFT_N];
    double outreal[DFT_N];
    double outimag[DFT_N];
    double CheckReal[DFT_N];
    double CheckImag[DFT_N];

    int times = 0;

    //读取输入数据
    for (int iThNum = 0;iThNum != DFT_N;++iThNum)
    { 
	double A,B;
	inputFile >> InputReal[iThNum];
	inputFile >> InputImag[iThNum];
	A = ConTernaryNumberdouble(InputReal[iThNum]);  
	B = ConTernaryNumberdouble(InputImag[iThNum]);  
	inreal[iThNum] = A;
	inimag[iThNum] = B;  
	CheckReal[iThNum] = A;
	CheckImag[iThNum] = B;  

	outputFile << "第" << iThNum << "个输入数MSD表示： " <<  InputReal[iThNum]<< endl;
	outputFile << "第" << iThNum << "个输入数MSD表示： " <<  InputImag[iThNum]<< endl;
	outputFile << "第" << iThNum << "个输入数十进制表示： "<< A << endl;
	outputFile << "第" << iThNum << "个输入数十进制表示： "<< B << endl;
    }
    cout << "Input done ...." << endl;

    //十进制数计算DFT
    cout << "Decimble-DFT computing ..." << endl;
    ConDFTDec(inreal,inimag,outreal,outimag);
    cout << "Decimble-DFT done ..." << endl;

    //MSD数计算DFT
    cout << "MSD-DFT computing ..." << endl;
    ComDFTMSD(InputReal,InputImag,OutputReal,OutputImag);
    cout << "MSD-DFT done ..." << endl;

    //用一般方法计算FFT作为检验值
    cout << "check FFT computing ..." << endl;
    FFT(CheckReal,CheckImag);
    cout << "check FFT done ..." << endl;

    //基4-FFT计算
    cout << "Radix-4 FFT computing ..." << endl;
    ComRadix_4_FFT_MSD(InputReal,InputImag,Radix4OutputReal,Radix4OutputImag);
    cout << "Radix-4 FFT done ..." << endl;

    //基8-FFT计算
    cout << "Radix-8 FFT computing ..." << endl;
    ComRadix_8_FFT_MSD(InputReal,InputImag,Radix8OutputReal,Radix8OutputImag);
    cout << "Radix-8 FFT done ..." << endl;

    //输出结果
    for(int i = 0;i != DFT_N;++i)
    {

	outputFile << "MSD数DFT计算结果(MSD数)：\t"; 
	outputFile << "第 " << i << " 维实部：\t" << OutputReal[i] << endl; 
	outputFile << "\t\t\t\t第 " << i << " 维虚部：\t" << OutputImag[i] << endl; 

	outputFile << "基4FFT算法计算结果(MSD数)：" ; 
	outputFile << "\t第 " << i << " 维实部：\t" << Radix4OutputReal[i] << endl; 
	outputFile << "\t\t\t\t第 " << i << " 维虚部：\t" << Radix4OutputImag[i] << endl; 

	outputFile << "基8FFT算法计算结果(MSD数)：" ; 
	outputFile << "\t第 " << i << " 维实部：\t" << Radix8OutputReal[i] << endl; 
	outputFile << "\t\t\t\t第 " << i << " 维虚部：\t" << Radix8OutputImag[i] << endl; 

	outputFile << "MSD数DFT计算结果(十进制数)：\t"; 
	outputFile << "第 " << i << " 维实部：\t" << ConTernaryNumberdouble(OutputReal[i]); 
	outputFile << "\t\t第 " << i << " 维虚部：\t" << ConTernaryNumberdouble(OutputImag[i]) << endl; 

	outputFile << "基4FFT算法计算结果(十进制数)：\t" ; 
	outputFile << "第 " << i << " 维实部：\t" << ConTernaryNumberdouble(Radix4OutputReal[i]); 
	outputFile << "\t\t第 " << i << " 维虚部：\t" << ConTernaryNumberdouble(Radix4OutputImag[i]) << endl; 
	
	outputFile << "基8FFT算法计算结果(十进制数)：\t" ; 
	outputFile << "第 " << i << " 维实部：\t" << ConTernaryNumberdouble(Radix8OutputReal[i]); 
	outputFile << "\t\t第 " << i << " 维虚部：\t" << ConTernaryNumberdouble(Radix8OutputImag[i]) << endl; 

	outputFile << "十进制数DFT计算结果：\t\t" ; 
	outputFile << "第 " << i << " 维实部：\t" << outreal[i]; 
	outputFile << "\t\t第 " << i << " 维虚部：\t" << outimag[i] << endl; 

	outputFile << "验证FFT算法计算结果：\t\t" ; 
	outputFile << "第 " << i << " 维实部：\t" << CheckReal[i]; 
	outputFile << "\t\t第 " << i << " 维虚部：\t" << CheckImag[i] << endl << endl; 



//	outputFile << "MSD数表示：" << endl; 
//	outputFile << "第 " << i << " 维实部：" << endl; 
//	outputFile << OutputReal[i] << endl; 
//	outputFile << ConTernaryNumberdouble(OutputReal[i]) << endl;
//	outputFile << "第 " << i << " 维虚部：" << endl; 
//	outputFile << OutputImag[i] << endl; 
//	outputFile << ConTernaryNumberdouble(OutputImag[i]) << endl;
//
//
//	outputFile << "十进制数表示：" << endl; 
//	outputFile << "第 " << i << " 维实部：" << endl; 
//	outputFile << outreal[i] << endl; 
//	outputFile << "第 " << i << " 维虚部：" << endl; 
//	outputFile << outimag[i] << endl << endl; 

    }
    outputFile << "所需液晶数目： " << G_YEJING << endl;
    outputFile << "所需时钟周期数： " << G_TIME << endl;

    return 0;
}
