//
//  CreatUseCase.cpp
//  

#include <stdlib.h>
#include <fstream>
#include <sstream>
using namespace std;

//用0、1和u符号串表示的三值数据
typedef string  TERNARY_DATA;

//主函数用于生成加法的数据文件
int main()
{
    const int DFT_N = 32;
    const int MSD_m = 32;
    const int MSD_p = 24;
    int temp = MSD_m - MSD_p;
    
    stringstream SDFT_N,SMSD_m,SMSD_p;
    SDFT_N << DFT_N;
    SMSD_m << MSD_m;
    SMSD_p << MSD_p;
    string OUTPUTFILE_s = "UseCase" + SDFT_N.str() + "_" + SMSD_m.str() + "_" + SMSD_p.str() + ".test";
    const char* OUTPUTFILE = OUTPUTFILE_s.c_str();

    //加法数据文件
    ofstream outputFile(OUTPUTFILE);


    //随机产生10组50位验证数据和10组400位测试数据
    char ternarySymbol[3] = {'0', '1', 'u'};
    TERNARY_DATA operandA50(DFT_N, '0'), operandB50(DFT_N, '0');

    for(unsigned uIndexGroup = 0; uIndexGroup < DFT_N; ++uIndexGroup)
    {
	for(unsigned uIndexBit = 0; uIndexBit < MSD_m; ++uIndexBit)
//	for(unsigned uIndexBit = MSD_m - MSD_p - 8; uIndexBit < MSD_m; ++uIndexBit)
	{
	    operandA50[uIndexBit] = ternarySymbol[rand() % 3];
	    operandB50[uIndexBit] = ternarySymbol[rand() % 3];
	}
	outputFile << operandA50 << endl;
	outputFile << operandB50 << endl << endl;
    }


    return 0;
}
