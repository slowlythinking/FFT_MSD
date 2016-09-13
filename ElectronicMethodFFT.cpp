#include <fstream>
#include <sstream>
#include <iomanip>
#include <iostream>
#include <time.h>
#include <cmath>
using namespace std;
typedef string TERNARY_DATA;

const int DFT_N=1024;//N点：输入向量为N维
const int MSD_m=16;//m位：输入向量每一维为m位
const int MSD_p=8;//小数部分位数（MSD precision）

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

int main()
{
    stringstream SDFT_N,SMSD_m,SMSD_p;
    SDFT_N << DFT_N;
    SMSD_m << MSD_m;
    SMSD_p << MSD_p;
    string INPUTFILE_s = "UseCase" + SDFT_N.str() + "_" + SMSD_m.str() + "_" + SMSD_p.str() + ".test";
    string OUTPUTFILE_s = "ElectronicResult" + SDFT_N.str() + "_" + SMSD_m.str() + "_" + SMSD_p.str() + ".test";
    const char* INPUTFILE = INPUTFILE_s.c_str();
    const char* OUTPUTFILE = OUTPUTFILE_s.c_str();


    ifstream inputFile(INPUTFILE);
    ofstream outputFile(OUTPUTFILE);
    outputFile.precision(20);

    TERNARY_DATA InputReal[DFT_N];
    TERNARY_DATA InputImag[DFT_N];
    double inreal[DFT_N];
    double inimag[DFT_N];
    double outreal[DFT_N];
    double outimag[DFT_N];

    int times = 0;

    //读取输入数据
    cout << "Read data from input file" << endl;
    for (int iThNum = 0;iThNum != DFT_N;++iThNum)
    { 
	double A,B;
	inputFile >> InputReal[iThNum];
	inputFile >> InputImag[iThNum];
	A = ConTernaryNumberdouble(InputReal[iThNum]);  
	B = ConTernaryNumberdouble(InputImag[iThNum]);  
	inreal[iThNum] = A;
	inimag[iThNum] = B;  

	outputFile << "第" << iThNum << "个输入数MSD表示： " <<  InputReal[iThNum]<< endl;
	outputFile << "第" << iThNum << "个输入数MSD表示： " <<  InputImag[iThNum]<< endl;
	outputFile << "第" << iThNum << "个输入数十进制表示： "<< A << endl;
	outputFile << "第" << iThNum << "个输入数十进制表示： "<< B << endl;
    }
    cout << "Input done ...." << endl;

    cout << "Computing FFT..." << endl;
    
    //电子方法计算FFT
    clock_t StartTime = clock();
    FFT(inreal,inimag);
    clock_t EndTime = clock();

    cout << "Computing FFT done ..." << endl;

    //输出结果
    for(int i = 0;i != DFT_N;++i)
    {
	outputFile << "验证FFT算法计算结果：\t\t" ; 
	outputFile << "第 " << i << " 维实部：\t" << inreal[i]; 
	outputFile << "\t\t第 " << i << " 维虚部：\t" << inimag[i] << endl << endl; 
    }
    outputFile << "Number of running cycles is:" << EndTime - StartTime << endl;
    cout << "Number of running cycles is:" << EndTime - StartTime << endl;

    return 0;
}
