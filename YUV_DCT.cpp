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

void DCT4(double data[], double out[]){
	double temp1[16],temp2[16];
	for(int i = 0;i < 16;++i){
		if(data[i] <= 0) temp1[i] = data[i] - 1.0; ////solve the round-off error
		else temp1[i] = data[i];
	} 
	double a = 0.5, b = sqrt(2.0/5.0);
	double Cf[4][4] = {{ 1,  1,  1,  1},
					   { 2,  1, -1, -2},
					   { 1, -1, -1,  1},
					   { 1, -2,  2, -1}};
	double Ef[4][4] = {{  a*a, a*b/2,   a*a, a*b/2},
					   {a*b/2, b*b/4, a*b/2, b*b/4},
					   {  a*a, a*b/2,   a*a, a*b/2},
					   {a*b/2, b*b/4, a*b/2, b*b/4}};  
	for(int i = 0;i < 4;++i){
		for(int j = 0;j < 4;++j){
			double sum = 0.0;
			for(int k = 0;k < 4;++k){
				sum += Cf[i][k]*temp1[j+k*4];
			}
			temp2[i*4+j] = sum;
			//cout<<sum<<endl;
		}
	}
	for(int i = 0;i < 4;++i){
		for(int j = 0;j < 4;++j){
			double sum = 0.0;
			for(int k = 0;k < 4;++k){
				sum += temp2[i*4+k]*Cf[j][k];
			}
			out[i*4+j] = sum;
			//cout<<sum<<endl;
		}
	}
	for(int i = 0;i < 16;++i){
		out[i] = out[i]*Ef[i/4][i%4];
	}
	return;
}

void IDCT4(double data[], double out[]){
	double temp[16];
	double a = 0.5, b = sqrt(2.0/5.0);
	double Ci[4][4] = {{  1,   1,    1,    1},
					  {  1, 0.5, -0.5,   -1},
					  {  1,  -1,   -1,    1},
					  {0.5,  -1,    1, -0.5}};
	double E[4][4] = {{a*a, a*b, a*a, a*b},
					  {a*b, b*b, a*b, b*b},
					  {a*a, a*b, a*a, a*b},
					  {a*b, b*b, a*b, b*b}};
	for(int i = 0;i < 16;++i){
		out[i] = data[i]*E[i/4][i%4];
		//cout<<out[i]<<endl;
	}
	for(int i = 0;i < 4;++i){
		for(int j = 0;j < 4;++j){
			double sum = 0.0;
			for(int k = 0;k < 4;++k){
				sum += Ci[k][i]*out[j+k*4];
			}
			temp[i*4+j] = sum;
			//cout<<sum<<endl;
		}
	}
	for(int i = 0;i < 4;++i){
		for(int j = 0;j < 4;++j){
			double sum = 0.0;
			for(int k = 0;k < 4;++k){
				//if(i == 0 && j == 0) cout<<sum<<endl;
				sum += temp[i*4+k]*Ci[k][j];
			}
			out[i*4+j] = sum;
			//cout<<sum<<endl;
		}
	}
	for(int i = 0;i < 16;++i){
		//cout<<out[i]<<endl;
		if(out[i] < 0.0){
			cout<<out[i]<<endl;
			if(out[i] <= -1.00000000001){ //solve the round-off error
				out[i] += 1.0;
			}
			else{
				out[i] = 0.0;
			}
			//cout<<out[i]<<endl;
		}
		//cout<<out[i]<<endl;
	}
	return;
}

void DCT2(double data[], double out[]){
	double temp[4];
	double Cf[2][2] = {{ 1,  1},
					   { 1, -1}};		  
	for(int i = 0;i < 2;++i){
		for(int j = 0;j < 2;++j){
			double sum = 0.0;
			for(int k = 0;k < 2;++k){
				sum += Cf[i][k]*data[j+k*2];
			}
			temp[i*2+j] = sum;
			//cout<<sum<<endl;
		}
	}
	for(int i = 0;i < 2;++i){
		for(int j = 0;j < 2;++j){
			double sum = 0.0;
			for(int k = 0;k < 2;++k){
				sum += temp[i*2+k]*Cf[j][k];
			}
			out[i*2+j] = sum;
			//cout<<sum<<endl;
		}
	}
	return;
}

void IDCT2(double data[], double out[]){
	double temp[4];
	double Ci[2][2] = {{ 1,  1},
					   { 1, -1}};
	for(int i = 0;i < 2;++i){
		for(int j = 0;j < 2;++j){
			double sum = 0.0;
			for(int k = 0;k < 2;++k){
				sum += Ci[k][i]*data[j+k*2];
			}
			temp[i*2+j] = sum;
			//cout<<sum<<endl;
		}
	}
	for(int i = 0;i < 2;++i){
		for(int j = 0;j < 2;++j){
			double sum = 0.0;
			for(int k = 0;k < 2;++k){
				sum += temp[i*2+k]*Ci[k][j];
			}
			out[i*2+j] = sum/4;
			//cout<<sum<<endl;
		}
	}
	return;
}

int main(){
	//char filename1[30] = "AKIYO_352x288_10.yuv";
	//char filename2[30] = "MOBILE_352x288_10.yuv";
	//BYTE *data = new BYTE[dataLen];
	//Input(filename1,data);
	//Output(filename1,data);
	double test[16] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15};
	double test2[4] = {0,1,2,3};
	double p[16],p2[16],p3[4],p4[4];
	//DCT4(test, p);
	//IDCT4(p, p2);
	DCT2(test2, p3);
	IDCT2(p3, p4);
	for(int i = 0;i < 16;++i){
		//cout<<i<<" o= "<<test[i]<<endl;
		//cout<<i<<" = "<<p2[i]<<endl;
	}
	for(int i = 0;i < 4;++i){
		cout<<i<<" = "<<p4[i]<<endl;
	}

	return 0;
}
