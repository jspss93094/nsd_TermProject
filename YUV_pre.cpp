#include <iostream>
#include <fstream>
#include <memory.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <unistd.h>
#include <cmath>

#define WIDTH 8
#define HEIGHT 8
//352*288
using namespace std;

typedef unsigned char BYTE;

FILE *f_in, *f_out;
int fsize = WIDTH*HEIGHT, fnum = 10;
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

double Pre_ver(double data[]){
	double diff = 0.0;
	float pred[fsize];
	int nrow = HEIGHT / 4;
	int ncol = WIDTH / 4;
	for(int i = 0;i < fsize;++i){
		pred[i] = 128; //edge
	}
	for(int y = 1;y < nrow;++y){ //y = 0 is edge
		for(int x = 0;x < ncol;++x){
			for(int i = 0;i < 4;++i){
				double temp = data[(y*4-1)*WIDTH + x*4+i];
				for(int j = 0;j < 4;++j){
					pred[(y*4+j)*WIDTH + x*4+i] = temp;
				}
			}
		}
	}
	for(int i = 0;i < fsize;++i){
		diff += abs(data[i] - pred[i]);
		//cout<<data[i] - pred[i]<<" "<<abs(data[i] - pred[i])<<endl;
	}
	return diff;
}

double Pre_hor(double data[]){
	double diff = 0.0;
	float pred[fsize];
	int nrow = HEIGHT / 4;
	int ncol = WIDTH / 4;
	for(int i = 0;i < fsize;++i){
		pred[i] = 128; //edge
	}
	for(int y = 0;y < nrow;++y){
		for(int x = 1;x < ncol;++x){ //x = 0 is edge
			for(int i = 0;i < 4;++i){
				double temp = data[(y*4+i)*WIDTH + x*4-1];
				for(int j = 0;j < 4;++j){
					pred[(y*4+i)*WIDTH + x*4+j] = temp;
				}
			}
		}
	}
//	for(int i = 0;i < HEIGHT;++i){
//		for(int j = 0;j < WIDTH;++j){
//			cout<<data[i*WIDTH+j]<<" ";
//		}
//		cout<<endl;
//	}
//	cout<<endl;
//	for(int i = 0;i < HEIGHT;++i){
//		for(int j = 0;j < WIDTH;++j){
//			cout<<pred[i*WIDTH+j]<<" ";
//		}
//		cout<<endl;
//	}
	for(int i = 0;i < fsize;++i){
		diff += abs(data[i] - pred[i]);
		//cout<<data[i] - pred[i]<<" "<<abs(data[i] - pred[i])<<endl;
	}
	return diff;
}

double Pre_DC(double data[]){
	double diff = 0.0;
	float pred[fsize];
	int nrow = HEIGHT / 4;
	int ncol = WIDTH / 4;
	for(int y = 0;y < nrow;++y){
		for(int x = 0;x < ncol;++x){
			for(int i = 0;i < 4;++i){
				double temp1;
				if(x != 0) temp1 = data[(y*4+i)*WIDTH + x*4-1];
				else temp1 = 128;
				for(int j = 0;j < 4;++j){
					double temp2;
					if(y != 0) temp2 = data[(y*4-1)*WIDTH + x*4+j];
					else temp2 = 128;
					pred[(y*4+i)*WIDTH + x*4+j] = (temp1 + temp2)/2;
				}
			}
		}
	}
//	for(int i = 0;i < HEIGHT;++i){
//		for(int j = 0;j < WIDTH;++j){
//			cout<<data[i*WIDTH+j]<<" ";
//		}
//		cout<<endl;
//	}
//	cout<<endl;
//	for(int i = 0;i < HEIGHT;++i){
//		for(int j = 0;j < WIDTH;++j){
//			cout<<pred[i*WIDTH+j]<<" ";
//		}
//		cout<<endl;
//	}
	for(int i = 0;i < fsize;++i){
		diff += abs(data[i] - pred[i]);
		//cout<<data[i] - pred[i]<<" "<<abs(data[i] - pred[i])<<endl;
	}
	return diff;
}

double Pre_ref(double data[],double refdata[],int max_ran){
	int min_ran = 0 - max_ran;
	double diff = 0.0;
	float pred[fsize];
	int nrow = HEIGHT / 4;
	int ncol = WIDTH / 4;
	for(int y = 0;y < nrow;++y){
		for(int x = 0;x < ncol;++x){
			float mindif = -1.0;
			int mv_x, mv_y;
			for(int sy = min_ran;sy <= max_ran;++sy){
				for(int sx = min_ran;sx <= max_ran;++sx){
					if(y*4+sy < 0 || x*4+sx < 0 || y*4+sy+3 >= HEIGHT || x*4+sx+3 >= WIDTH) continue;   //detect edge
					float dif = 0.0;
					for(int ry = 0;ry < 4;++ry){
						for(int rx = 0;rx < 4;++rx){
							dif += abs(data[(y*4+ry)*WIDTH + (x*4+rx)] - refdata[(y*4+ry+sy)*WIDTH+(x*4+rx+sx)]);  //compute difference
							//cout<<data[(y*4+ry)*WIDTH + (x*4+rx)] - refdata[(y*4+ry+sy)*WIDTH+(x*4+rx+sx)]<<" ";
						}
					}
					cout<<endl;
					if(dif < mindif || mindif < 0.0){
						mv_x = sx;
						mv_y = sy;
						mindif = dif;
					}
				}
			}
			for(int ry = 0;ry < 4;++ry){
				for(int rx = 0;rx < 4;++rx){
					pred[(y*4+ry)*WIDTH + x*4+rx] = refdata[(y*4+ry+mv_y)*WIDTH+(x*4+rx+mv_x)];
					//cout<<mv_x<<" "<<mv_y<<" "<<refdata[(y*4+ry+mv_y)*WIDTH+(x*4+rx+mv_x)]<<endl;
				}
			}
		}
	}
	for(int i = 0;i < HEIGHT;++i){
		for(int j = 0;j < WIDTH;++j){
			cout<<refdata[i*WIDTH+j]<<" ";
		}
		cout<<endl;
	}
	cout<<endl;
	for(int i = 0;i < HEIGHT;++i){
		for(int j = 0;j < WIDTH;++j){
			cout<<data[i*WIDTH+j]<<" ";
		}
		cout<<endl;
	}
	cout<<endl;
	for(int i = 0;i < HEIGHT;++i){
		for(int j = 0;j < WIDTH;++j){
			cout<<pred[i*WIDTH+j]<<" ";
		}
		cout<<endl;
	}
	for(int i = 0;i < fsize;++i){
		diff += abs(data[i] - pred[i]);
		//cout<<data[i] - pred[i]<<" "<<abs(data[i] - pred[i])<<endl;
	}
	return diff;
}

int Pre_compare(double data[]){
	double dif_ver, dif_hor, dif_DC, dif_ref;
	return 0;
}

int main(){
	//char filename1[30] = "AKIYO_352x288_10.yuv";
	//char filename2[30] = "MOBILE_352x288_10.yuv";
	//BYTE *data = new BYTE[dataLen];
	//Input(filename1,data);
	//Output(filename1,data);
	double test[16*4] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,255,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,255,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,255,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,255};
	//double test1[16] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,255};
	double reftest[16*4] = {0,1,2,3,4,5,6,7,0,1,2,3,4,5,6,7,0,1,2,3,4,5,6,7,0,1,2,3,4,5,6,7,0,1,2,3,4,5,6,7,0,1,2,3,4,5,6,7,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,255};
	//double reftest1[16] = {8,9,10,11,12,13,14,255,0,1,2,3,4,5,6,7};
	double test2[4] = {0,1,2,3};
	double p[16],p2[16] = {0},p3[4],p4[4];
	double dif;
	//dif = Pre_ver(test);
	//dif = Pre_hor(test);
	//dif = Pre_DC(test);
	dif = Pre_ref(test,reftest,8);
	//cout<<dif<<endl;
	//for(int i = 0;i < 16;++i){
		//cout<<i<<" o= "<<test[i]<<endl;
		//cout<<i<<" = "<<p2[i]<<endl;
	//}
//	for(int i = 0;i < 4;++i){
//		cout<<i<<" = "<<p4[i]<<endl;
//	}

	return 0;
}
