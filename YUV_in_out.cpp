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

int main(){
	char filename1[30] = "AKIYO_352x288_10.yuv";
	char filename2[30] = "MOBILE_352x288_10.yuv";
	BYTE *data = new BYTE[dataLen];
	Input(filename1,data);
	Output(filename1,data);

	return 0;
}
