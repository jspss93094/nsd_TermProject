#include <iostream>
#include <fstream>
#include <memory.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <unistd.h>
#include <cmath>

#define WIDTH 352
#define HEIGHT 288

using namespace std;

typedef unsigned char BYTE;

FILE *f_in, *f_out;
int fsize = WIDTH*HEIGHT, fnum = 100;
int uvsize = fsize / 4;  //yuv 4:2:0
int onesize = fsize + uvsize*2;
int dataLen = onesize * fnum;

float Qstep[52] = {0.625, 0.6875, 0.8125, 0.875,   1, 1.125,
				    1.25,  1.375,  1.625,  1.75,   2,  2.25,
				     2.5,   2.75,   3.25,   3.5,   4,   4.5,
				       5,    5.5,    6.5,     7,   8,     9,
				      10,     11,     13,    14,  16,    18,
				      20,     22,     26,    28,  32,    36,
				      40,     44,     52,    56,  64,    72,
				      80,     88,    104,   112, 128,   144,
				     160,    176,    208,   224};

void Input(char filename[], BYTE data[]){
	f_in = fopen(filename,"rb");
	fread(data,dataLen,1,f_in);
	fclose(f_in);
	
//	index = 0;
//	for(int i = 0;i < outsize;++i){
//		if(i >= 176*144 + index * outone && i < outone * (index + 1)) dataout[i] = 128;  //u,v = 128
//		else dataout[i] = int((adata[i] - min) / (max - min) * 255);  //0 ~ 255
//		if(i == outone * (index + 1) - 1) ++index;
//	}
	
	return;
}

void Output(char filename[],BYTE data[]){
	char fileout[50] = "";
	for(int i = 0;i < strlen(filename)-15;++i){
		fileout[i] = filename[i];
	}
	strncat (fileout, "_after.yuv", 49);
	
	f_out = fopen(fileout,"wb");
	fwrite(data, dataLen, 1, f_out);
	fclose(f_out);
	
	return;
}

void Quant(double data[], double out[], int qp, int size){
	double qs = Qstep[qp];
	for(int i = 0;i < size;++i){
		out[i] = round(data[i] / qs);
	}
	return;
}

void IQuant(double data[], double out[], int qp, int size){
double qs = Qstep[qp];
	for(int i = 0;i < size;++i){
		out[i] = data[i] * qs;
	}
	return;
}

int main(){
	//char filename1[30] = "AKIYO_352x288_10.yuv";
	//char filename2[30] = "MOBILE_352x288_10.yuv";
	//BYTE *data = new BYTE[dataLen];
	//Input(filename1,data);
	//Output(filename1,data);
	double test[16] = {0,1,-2,3,-4,5,-6,7,-8,9,-10,11,-12,13,-14,255};
	double test2[4] = {0,1,2,3};
	double p[16],p2[16],p3[4],p4[4];
	Quant(test,p, 4,16);
	IQuant(p,p2, 4,16);
	for(int i = 0;i < 16;++i){
		cout<<i<<" = "<<p2[i]<<endl;
	}

	return 0;
}
