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
//352*288
using namespace std;

typedef unsigned char BYTE;

FILE *f_in, *f_out;
int fsize = WIDTH*HEIGHT, fnum = 10;
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

//input and output
void Input(char filename[], BYTE data[]){
	f_in = fopen(filename,"rb");
	fread(data,dataLen,1,f_in);
	fclose(f_in);	
	return;
}

void Output(char filename[],BYTE data[], int dsize, bool iscp){
	char fileout[50] = "";
	for(int i = 0;i < strlen(filename)-4;++i){
		fileout[i] = filename[i];
	}
	if(iscp) strncat (fileout, "_cp.yuv", 49);
	else strncat (fileout, "_de.yuv", 49);
	
	f_out = fopen(fileout,"wb");
	fwrite(data, dsize, 1, f_out);
	fclose(f_out);
	
	return;
}

//frame to block
void ftob(BYTE data[],BYTE bdata[]){
	int nrow = HEIGHT / 4;
	int ncol = WIDTH / 4;
	for(int y = 0;y < nrow;++y){
		for(int x = 0;x < ncol;++x){
			for(int j = 0;j < 4;++j){
				for(int i = 0;i < 4;++i){
					bdata[y*ncol*16+x*16+j*4+i] = data[(y*4+j)*WIDTH + x*4+i];
				}
			}
		}
	}
	return;
}

void uvtob(BYTE data[], BYTE udata[], BYTE vdata[]){
	int nrow = HEIGHT / 4;
	int ncol = WIDTH / 4;
	for(int y = 0;y < nrow;++y){
		for(int x = 0;x < ncol;++x){
			for(int j = 0;j < 2;++j){
				for(int i = 0;i < 2;++i){
					udata[y*ncol*4+x*4+j*2+i] = data[(y*2+j)*WIDTH/2 + x*2+i];
					vdata[y*ncol*4+x*4+j*2+i] = data[uvsize + (y*2+j)*WIDTH/2 + x*2+i];
				}
			}
		}
	}
	return;
}

//block to frame
void btof(BYTE data[],BYTE fdata[]){
	int nrow = HEIGHT / 4;
	int ncol = WIDTH / 4;
	for(int y = 0;y < nrow;++y){
		for(int x = 0;x < ncol;++x){
			for(int j = 0;j < 4;++j){
				for(int i = 0;i < 4;++i){
					fdata[(y*4+j)*WIDTH + x*4+i] = data[y*ncol*16+x*16+j*4+i];
				}
			}
		}
	}
	return;
}

void btouv(BYTE udata[], BYTE vdata[], BYTE fdata[]){
	int nrow = HEIGHT / 4;
	int ncol = WIDTH / 4;
	for(int y = 0;y < nrow;++y){
		for(int x = 0;x < ncol;++x){
			for(int j = 0;j < 2;++j){
				for(int i = 0;i < 2;++i){
					fdata[(y*2+j)*WIDTH/2 + x*2+i] = udata[y*ncol*4+x*4+j*2+i];
					fdata[uvsize + (y*2+j)*WIDTH/2 + x*2+i] = vdata[y*ncol*4+x*4+j*2+i];
				}
			}
		}
	}
	return;
}

//prediction
double Pre_ver(BYTE data[], BYTE ref[], int x, int y){
	double diff = 0.0;
	BYTE pred[16];
	for(int i = 0;i < 16;++i){
		pred[i] = 128; //edge
	}
	if(y != 0){
		for(int i = 0;i < 4;++i){
			for(int j = 0;j < 4;++j){
				pred[i*4+j] = ref[(y*4-1)*WIDTH + x*4+j];
			}
		}
	}
	for(int i = 0;i < 16;++i){
		diff += abs(data[i] - pred[i]);
		//cout<<data[i] - pred[i]<<" "<<abs(data[i] - pred[i])<<endl;
	}
	return diff;
}

void Pre_ver(BYTE data[], BYTE ref[], int x, int y, BYTE resid[]){
	BYTE pred[16];
	for(int i = 0;i < 16;++i){
		pred[i] = 128; //edge
	}
	if(y != 0){
		for(int i = 0;i < 4;++i){
			for(int j = 0;j < 4;++j){
				pred[i*4+j] = ref[(y*4-1)*WIDTH + x*4+j];
			}
		}
	}
	for(int i = 0;i < 16;++i){
		resid[i] = data[i] - pred[i];
	}
	return;
}

void IPre_ver(BYTE data[], BYTE ref[], int x, int y, BYTE sum[]){
	BYTE pred[16];
	for(int i = 0;i < 16;++i){
		pred[i] = 128; //edge
	}
	if(y != 0){
		for(int i = 0;i < 4;++i){
			for(int j = 0;j < 4;++j){
				pred[i*4+j] = ref[(y*4-1)*WIDTH + x*4+j];
			}
		}
	}
	for(int i = 0;i < 16;++i){
		sum[i] = data[i] + pred[i];
	}
	return;
}

void Pre_ver_uv(BYTE data[], BYTE ref[], int x, int y, BYTE resid[]){
	BYTE pred[4];
	for(int i = 0;i < 4;++i){
		pred[i] = 128; //edge
	}
	if(y != 0){
		for(int i = 0;i < 2;++i){
			for(int j = 0;j < 2;++j){
				pred[i*2+j] = ref[(y*2-1)*WIDTH/2 + x*2+j];
			}
		}
	}
	for(int i = 0;i < 4;++i){
		resid[i] = data[i] - pred[i];
	}
	return;
}

void IPre_ver_uv(BYTE data[], BYTE ref[], int x, int y, BYTE sum[]){
	BYTE pred[4];
	for(int i = 0;i < 4;++i){
		pred[i] = 128; //edge
	}
	if(y != 0){
		for(int i = 0;i < 2;++i){
			for(int j = 0;j < 2;++j){
				pred[i*2+j] = ref[(y*2-1)*WIDTH/2 + x*2+j];
			}
		}
	}
	for(int i = 0;i < 4;++i){
		sum[i] = data[i] + pred[i];
	}
	return;
}

double Pre_hor(BYTE data[], BYTE ref[], int x, int y){
	double diff = 0.0;
	BYTE pred[16];
	for(int i = 0;i < 16;++i){
		pred[i] = 128; //edge
	}
	if(x != 0){
		for(int i = 0;i < 4;++i){
			for(int j = 0;j < 4;++j){
				pred[i*4+j] = ref[(y*4+i)*WIDTH + x*4-1];
			}
		}
	}
	for(int i = 0;i < 16;++i){
		diff += abs(data[i] - pred[i]);
		//cout<<data[i] - pred[i]<<" "<<abs(data[i] - pred[i])<<endl;
	}
	return diff;
}

void Pre_hor(BYTE data[], BYTE ref[], int x, int y,BYTE resid[]){
	BYTE pred[16];
	for(int i = 0;i < 16;++i){
		pred[i] = 128; //edge
	}
	if(x != 0){
		for(int i = 0;i < 4;++i){
			for(int j = 0;j < 4;++j){
				pred[i*4+j] = ref[(y*4+i)*WIDTH + x*4-1];
			}
		}
	}
	for(int i = 0;i < 16;++i){
		resid[i] = data[i] - pred[i];
	}
	return;
}

void IPre_hor(BYTE data[], BYTE ref[], int x, int y,BYTE sum[]){
	BYTE pred[16];
	for(int i = 0;i < 16;++i){
		pred[i] = 128; //edge
	}
	if(x != 0){
		for(int i = 0;i < 4;++i){
			for(int j = 0;j < 4;++j){
				pred[i*4+j] = ref[(y*4+i)*WIDTH + x*4-1];
			}
		}
	}
	for(int i = 0;i < 16;++i){
		sum[i] = data[i] + pred[i];
	}
	return;
}

void Pre_hor_uv(BYTE data[], BYTE ref[], int x, int y,BYTE resid[]){
	BYTE pred[4];
	for(int i = 0;i < 4;++i){
		pred[i] = 128; //edge
	}
	if(x != 0){
		for(int i = 0;i < 2;++i){
			for(int j = 0;j < 2;++j){
				pred[i*2+j] = ref[(y*2+i)*WIDTH/2 + x*2-1];
			}
		}
	}
	for(int i = 0;i < 4;++i){
		resid[i] = data[i] - pred[i];
	}
	return;
}

void IPre_hor_uv(BYTE data[], BYTE ref[], int x, int y,BYTE sum[]){
	BYTE pred[4];
	for(int i = 0;i < 4;++i){
		pred[i] = 128; //edge
	}
	if(x != 0){
		for(int i = 0;i < 2;++i){
			for(int j = 0;j < 2;++j){
				pred[i*2+j] = ref[(y*2+i)*WIDTH/2 + x*2-1];
			}
		}
	}
	for(int i = 0;i < 4;++i){
		sum[i] = data[i] + pred[i];
	}
	return;
}

double Pre_DC(BYTE data[], BYTE ref[], int x, int y){
	double diff = 0.0;
	float pred[16];
	for(int i = 0;i < 4;++i){
		BYTE temp1;
		if(x != 0) temp1 = ref[(y*4+i)*WIDTH + x*4-1];
		else temp1 = 128;
		for(int j = 0;j < 4;++j){
			BYTE temp2;
			if(y != 0) temp2 = ref[(y*4-1)*WIDTH + x*4+j];
			else temp2 = 128;
			pred[i*4+j] = (temp1 + temp2)/2.0;
		}
	}
	for(int i = 0;i < 16;++i){
		diff += abs(data[i] - pred[i]);
		//cout<<data[i] - pred[i]<<" "<<abs(data[i] - pred[i])<<endl;
	}
	return diff;
}

void Pre_DC(BYTE data[], BYTE ref[], int x, int y, BYTE resid[]){
	BYTE pred[16];
	for(int i = 0;i < 4;++i){
		BYTE temp1;
		if(x != 0) temp1 = ref[(y*4+i)*WIDTH + x*4-1];
		else temp1 = 128;
		for(int j = 0;j < 4;++j){
			BYTE temp2;
			if(y != 0) temp2 = ref[(y*4-1)*WIDTH + x*4+j];
			else temp2 = 128;
			pred[i*4+j] = (temp1 + temp2)/2;
		}
	}
	for(int i = 0;i < 16;++i){
		resid[i] = data[i] - pred[i];
	}
	return;
}

void IPre_DC(BYTE data[], BYTE ref[], int x, int y, BYTE sum[]){
	BYTE pred[16];
	for(int i = 0;i < 4;++i){
		BYTE temp1;
		if(x != 0) temp1 = ref[(y*4+i)*WIDTH + x*4-1];
		else temp1 = 128;
		for(int j = 0;j < 4;++j){
			BYTE temp2;
			if(y != 0) temp2 = ref[(y*4-1)*WIDTH + x*4+j];
			else temp2 = 128;
			pred[i*4+j] = (temp1 + temp2)/2;
		}
	}
	for(int i = 0;i < 16;++i){
		sum[i] = data[i] + pred[i];
	}
	return;
}

void Pre_DC_uv(BYTE data[], BYTE ref[], int x, int y, BYTE resid[]){
	BYTE pred[4];
	for(int i = 0;i < 2;++i){
		BYTE temp1;
		if(x != 0) temp1 = ref[(y*2+i)*WIDTH/2 + x*2-1];
		else temp1 = 128;
		for(int j = 0;j < 2;++j){
			BYTE temp2;
			if(y != 0) temp2 = ref[(y*2-1)*WIDTH/2 + x*2+j];
			else temp2 = 128;
			pred[i*2+j] = (temp1 + temp2)/2;
		}
	}
	for(int i = 0;i < 4;++i){
		resid[i] = data[i] - pred[i];
	}
	return;
}

void IPre_DC_uv(BYTE data[], BYTE ref[], int x, int y, BYTE sum[]){
	BYTE pred[4];
	for(int i = 0;i < 2;++i){
		BYTE temp1;
		if(x != 0) temp1 = ref[(y*2+i)*WIDTH/2 + x*2-1];
		else temp1 = 128;
		for(int j = 0;j < 2;++j){
			BYTE temp2;
			if(y != 0) temp2 = ref[(y*2-1)*WIDTH/2 + x*2+j];
			else temp2 = 128;
			pred[i*2+j] = (temp1 + temp2)/2;
		}
	}
	for(int i = 0;i < 4;++i){
		sum[i] = data[i] + pred[i];
	}
	return;
}

double Pre_ref(BYTE data[], BYTE refdata[], int x, int y, int max_ran){
	max_ran = abs(max_ran);
	int min_ran = 0 - max_ran;
	double diff = 0.0;
	BYTE pred[16];
	float mindif = -1.0;
	int mv_x, mv_y;
	for(int sy = min_ran;sy <= max_ran;++sy){
		for(int sx = min_ran;sx <= max_ran;++sx){
			if(y*4+sy < 0 || x*4+sx < 0 || y*4+sy+3 >= HEIGHT || x*4+sx+3 >= WIDTH) continue;   //detect edge
			float dif = 0.0;
			for(int ry = 0;ry < 4;++ry){
				for(int rx = 0;rx < 4;++rx){
					dif += abs(data[ry*4+rx] - refdata[(y*4+ry+sy)*WIDTH+(x*4+rx+sx)]);  //compute difference
					//cout<<data[(y*4+ry)*WIDTH + (x*4+rx)] - refdata[(y*4+ry+sy)*WIDTH+(x*4+rx+sx)]<<" ";
				}
			}
			//cout<<endl;
			if(dif < mindif || mindif < 0.0){
				mv_x = sx;
				mv_y = sy;
				mindif = dif;
			}
		}
	}
	for(int ry = 0;ry < 4;++ry){
		for(int rx = 0;rx < 4;++rx){
			pred[ry*4+rx] = refdata[(y*4+ry+mv_y)*WIDTH+(x*4+rx+mv_x)];
		}
	}
	for(int i = 0;i < 16;++i){
		diff += abs(data[i] - pred[i]);
		//cout<<data[i] - pred[i]<<" "<<abs(data[i] - pred[i])<<endl;
	}
	return diff;
}

void Pre_ref(BYTE data[], BYTE refdata[], int x, int y, int max_ran, BYTE resid[], int &mv_x, int &mv_y){
	max_ran = abs(max_ran);
	int min_ran = 0 - max_ran;
	BYTE pred[16];
	int mindif = -1;
	for(int sy = min_ran;sy <= max_ran;++sy){
		for(int sx = min_ran;sx <= max_ran;++sx){
			if(y*4+sy < 0 || x*4+sx < 0 || y*4+sy+3 >= HEIGHT || x*4+sx+3 >= WIDTH) continue;   //detect edge
			int dif = 0;
			for(int ry = 0;ry < 4;++ry){
				for(int rx = 0;rx < 4;++rx){
					dif += abs(data[ry*4+rx] - refdata[(y*4+ry+sy)*WIDTH+(x*4+rx+sx)]);  //compute difference
					//cout<<data[(y*4+ry)*WIDTH + (x*4+rx)] - refdata[(y*4+ry+sy)*WIDTH+(x*4+rx+sx)]<<" ";
				}
			}
			//cout<<endl;
			if(dif < mindif || mindif < 0){
				mv_x = sx;
				mv_y = sy;
				mindif = dif;
			}
		}
	}
	for(int ry = 0;ry < 4;++ry){
		for(int rx = 0;rx < 4;++rx){
			pred[ry*4+rx] = refdata[(y*4+ry+mv_y)*WIDTH+(x*4+rx+mv_x)];
		}
	}
	for(int i = 0;i < 16;++i){
		resid[i] = data[i] - pred[i];
	}
	return;
}

void IPre_ref(BYTE data[], BYTE refdata[], int x, int y, BYTE sum[], int mv_x, int mv_y){
	BYTE pred[16];
	for(int ry = 0;ry < 4;++ry){
		for(int rx = 0;rx < 4;++rx){
			pred[ry*4+rx] = refdata[(y*4+ry+mv_y)*WIDTH+(x*4+rx+mv_x)];
		}
	}
	for(int i = 0;i < 16;++i){
		sum[i] = data[i] + pred[i];
	}
	return;
}

void Pre_ref_uv(BYTE data[], BYTE refdata[], int x, int y, BYTE resid[], int mv_x, int mv_y){
	BYTE pred[4];
	BYTE refexp[fsize];
	int nrow = HEIGHT / 2;
	int ncol = WIDTH / 2;
	for(int ny = 0;ny < nrow;++ny){
		for(int nx = 0;nx < ncol;++nx){
			for(int j = 0;j < 2;++j){
				for(int i = 0;i < 2;++i){
					refexp[(ny*2+j)*WIDTH + nx*2+i] = refdata[ny*ncol + nx];
				}
			}
		}
	}
	for(int ry = 0;ry < 2;++ry){
		for(int rx = 0;rx < 2;++rx){
			pred[ry*2+rx] = refexp[(y*4+ry+mv_y)*WIDTH+(x*4+rx+mv_x)];
		}
	}
	for(int i = 0;i < 4;++i){
		resid[i] = data[i] - pred[i];
	}
	return;
}

void IPre_ref_uv(BYTE data[], BYTE refdata[], int x, int y, BYTE sum[], int mv_x, int mv_y){
	BYTE pred[4];
	BYTE refexp[fsize];
	int nrow = HEIGHT / 2;
	int ncol = WIDTH / 2;
	for(int ny = 0;ny < nrow;++ny){
		for(int nx = 0;nx < ncol;++nx){
			for(int j = 0;j < 2;++j){
				for(int i = 0;i < 2;++i){
					refexp[(ny*2+j)*WIDTH + nx*2+i] = refdata[ny*ncol + nx];
				}
			}
		}
	}
	for(int ry = 0;ry < 2;++ry){
		for(int rx = 0;rx < 2;++rx){
			pred[ry*2+rx] = refexp[(y*4+ry+mv_y)*WIDTH+(x*4+rx+mv_x)];
		}
	}
	for(int i = 0;i < 4;++i){
		sum[i] = data[i] + pred[i];
	}
	return;
}

//find the best prediction mode
int Pre_compare(BYTE data[], BYTE ref[], int x, int y){
	int ind = 0;
	double dif_ver, dif_hor, dif_DC;
	dif_ver = Pre_ver(data, ref, x, y);
	double min_val = dif_ver;
	dif_hor = Pre_hor(data, ref, x, y);
	if(dif_hor < min_val) ind = 1;
	
	dif_DC = Pre_DC(data, ref, x, y);
	if(dif_DC < min_val) ind = 2;
	
	//cout<<dif_ver<<" "<<dif_hor<<" "<<dif_DC<<endl;
	return ind;
}

int Pre_compare(BYTE data[], BYTE ref[], BYTE refb[], int x, int y, int range){
	int ind = 0;
	double dif_ver, dif_hor, dif_DC, dif_ref;
	dif_ver = Pre_ver(data, ref, x, y);
	double min_val = dif_ver;
	dif_hor = Pre_hor(data, ref, x, y);
	if(dif_hor < min_val) ind = 1;
	
	dif_DC = Pre_DC(data, ref, x, y);
	if(dif_DC < min_val) ind = 2;
	
	dif_ref = Pre_ref(data, refb, x, y, range);
	if(dif_ref < min_val) ind = 3;
	
	//cout<<dif_ver<<" "<<dif_hor<<" "<<dif_DC<<" "<<dif_ref<<endl;
	return ind;
}

//DCT
/*void DCT4(BYTE data[], BYTE out[]){
	BYTE temp1[16],temp2[16];
	for(int i = 0;i < 16;++i){
		if(data[i] < 0) temp1[i] = data[i] - 1; ////solve the round-off error
		else temp1[i] = data[i];
	} 
	float a = 0.5, b = sqrt(2.0/5.0);
	BYTE Cf[4][4] = {{ 1,  1,  1,  1},
					   { 2,  1, -1, -2},
					   { 1, -1, -1,  1},
					   { 1, -2,  2, -1}};
	BYTE Ef[4][4] = {{  a*a, a*b/2,   a*a, a*b/2},
					   {a*b/2, b*b/4, a*b/2, b*b/4},
					   {  a*a, a*b/2,   a*a, a*b/2},
					   {a*b/2, b*b/4, a*b/2, b*b/4}};  
	for(int i = 0;i < 4;++i){
		for(int j = 0;j < 4;++j){
			int sum = 0;
			for(int k = 0;k < 4;++k){
				sum += Cf[i][k]*temp1[j+k*4];
			}
			temp2[i*4+j] = sum;
			//cout<<sum<<endl;
		}
	}
	for(int i = 0;i < 4;++i){
		for(int j = 0;j < 4;++j){
			int sum = 0;
			for(int k = 0;k < 4;++k){
				sum += temp2[i*4+k]*Cf[j][k];
			}
			out[i*4+j] = (char)sum;
			//cout<<sum<<endl;
		}
	}
	for(int i = 0;i < 16;++i){
		out[i] = out[i]*Ef[i/4][i%4];
	}
	return;
}

void IDCT4(BYTE data[], BYTE out[]){
	BYTE temp[16];
	BYTE outtemp[16];
	float a = 0.5, b = sqrt(2.0/5.0);
	BYTE Ci[4][4] = {{  1,   1,    1,    1},
					  {  1, 0.5, -0.5,   -1},
					  {  1,  -1,   -1,    1},
					  {0.5,  -1,    1, -0.5}};
	BYTE E[4][4] = {{a*a, a*b, a*a, a*b},
					  {a*b, b*b, a*b, b*b},
					  {a*a, a*b, a*a, a*b},
					  {a*b, b*b, a*b, b*b}};
	for(int i = 0;i < 16;++i){
		outtemp[i] = data[i]*E[i/4][i%4];
		//cout<<out[i]<<endl;
	}
	for(int i = 0;i < 4;++i){
		for(int j = 0;j < 4;++j){
			int sum = 0;
			for(int k = 0;k < 4;++k){
				sum += Ci[k][i]*outtemp[j+k*4];
			}
			temp[i*4+j] = sum;
			//cout<<sum<<endl;
		}
	}
	for(int i = 0;i < 4;++i){
		for(int j = 0;j < 4;++j){
			int sum = 0;
			for(int k = 0;k < 4;++k){
				//if(i == 0 && j == 0) cout<<sum<<endl;
				sum += temp[i*4+k]*Ci[k][j];
			}
			outtemp[i*4+j] = (char)sum;
			//cout<<sum<<endl;
		}
	}
	for(int i = 0;i < 16;++i){
		//cout<<out[i]<<endl;
		if(outtemp[i] < 0.0){
			if(outtemp[i] <= -1.00000000001){ //solve the round-off error
				outtemp[i] += 1.0;
				out[i] = (char)outtemp[i];
			}
			else{
				//outtempt[i] = 0.0;
				out[i] = 0;
			}
			//cout<<out[i]<<endl;
		}
		//cout<<out[i]<<endl;
	}
	return;
}

void DCT2(BYTE data[], BYTE out[]){
	float temp[4];
	float Cf[2][2] = {{ 1,  1},
					   { 1, -1}};		  
	for(int i = 0;i < 2;++i){
		for(int j = 0;j < 2;++j){
			float sum = 0.0;
			for(int k = 0;k < 2;++k){
				sum += Cf[i][k]*data[j+k*2];
			}
			temp[i*2+j] = sum;
			//cout<<sum<<endl;
		}
	}
	for(int i = 0;i < 2;++i){
		for(int j = 0;j < 2;++j){
			float sum = 0.0;
			for(int k = 0;k < 2;++k){
				sum += temp[i*2+k]*Cf[j][k];
			}
			out[i*2+j] = sum;
			//cout<<sum<<endl;
		}
	}
	return;
}

void IDCT2(BYTE data[], BYTE out[]){
	float temp[4];
	float Ci[2][2] = {{ 1,  1},
					   { 1, -1}};
	for(int i = 0;i < 2;++i){
		for(int j = 0;j < 2;++j){
			float sum = 0.0;
			for(int k = 0;k < 2;++k){
				sum += Ci[k][i]*data[j+k*2];
			}
			temp[i*2+j] = sum;
			//cout<<sum<<endl;
		}
	}
	for(int i = 0;i < 2;++i){
		for(int j = 0;j < 2;++j){
			float sum = 0.0;
			for(int k = 0;k < 2;++k){
				sum += temp[i*2+k]*Ci[k][j];
			}
			out[i*2+j] = (char)(sum/4);
			//cout<<sum<<endl;
		}
	}
	return;
}*/

//Quantization
void Quant(BYTE data[], BYTE out[], int qp, int size){
	float qs = Qstep[qp];
	for(int i = 0;i < size;++i){
		out[i] = (char)(round(data[i] / qs));
	}
	return;
}

void IQuant(BYTE data[], BYTE out[], int qp, int size){
	float qs = Qstep[qp];
	for(int i = 0;i < size;++i){
		out[i] = (float)data[i] * qs;
	}
	return;
}

//Zigzag
int Zig(BYTE data[], BYTE after[]){
	int order[16] = {0, 1, 4, 8, 5, 2, 3, 6, 9, 12, 13, 10, 7, 11, 14, 15};
	int count = 0, zsize = 0;
	for(int i = 0;i < 16;++i){
		if(data[order[i]] == 0) ++count;
		else{
			after[zsize++] = data[order[i]];
			after[zsize++] = count;
			count = 0;
		}
	}
	after[zsize++] = 0;
	return zsize;
}

int Izig(BYTE data[], BYTE after[]){
	int order[16] = {0, 1, 4, 8, 5, 2, 3, 6, 9, 12, 13, 10, 7, 11, 14, 15};
	int index = 0;
	int zsize = 0;
	for(int i = 0;i < 33;i += 2){
		if(data[i] == 0){
			for(int j = index;j < 16;++j){
				after[order[index++]] = 0;
			}
			zsize++;
			break;
		}
		for(int j = 0;j < (int)data[i+1];++j){
			after[order[index++]] = 0;
		}
		zsize++;
		after[order[index++]] = data[i];
		zsize++;
	}
	return zsize;
}

//main function to call
void YUV_compress(char file[], int mode, int QP){
	BYTE *data = new BYTE[dataLen];
	BYTE *output = new BYTE[100000000];//block size | width/4 | height/4 | frame number | QP ...
	int osize = 0;
	Input(file,data);
	output[osize++] = 4;//block size
	output[osize++] = WIDTH/4;//width/4
	output[osize++] = HEIGHT/4;//height/4
	output[osize++] = (char)fnum;//frame number
	output[osize++] = (char)QP;//QP
	int nrow = HEIGHT / 4;
	int ncol = WIDTH / 4;
	int bmode;
	BYTE *bdata = new BYTE[fsize];
	BYTE *udata = new BYTE[uvsize];
	BYTE *vdata = new BYTE[uvsize];
	for(int f = 0;f < fnum;++f){
		cout<<f<<endl;
		ftob(&data[onesize*f],bdata);
		uvtob(&data[onesize*f+fsize], udata, vdata);
		for(int y = 0;y < nrow;++y){
			for(int x = 0;x < ncol;++x){
				//select mode
				if(mode == 4){
					if(f != 0){
						bmode = Pre_compare(&bdata[y*ncol*16+x*16], 
											&data[onesize*f], 
											&data[onesize*(f-1)], 
											x, 
											y, 
											0);
					}
					else{
						bmode = Pre_compare(&bdata[y*ncol*16+x*16], 
											&data[onesize*f], 
											x, 
											y);
					}
				}
				else bmode = mode;
				output[osize++] = (char)bmode;
				//use mode to compute the block residual
				BYTE *resid = new BYTE[24];
				if(bmode == 0){
					Pre_ver(&bdata[y*ncol*16+x*16], &data[onesize*f], x, y, &resid[0]);
					Pre_ver_uv(&udata[y*ncol*4+x*4], &data[onesize*f+fsize], x, y, &resid[16]);
					Pre_ver_uv(&vdata[y*ncol*4+x*4], &data[onesize*f+fsize+uvsize], x, y, &resid[20]);
				}
				else if(bmode == 1){
					Pre_hor(&bdata[y*ncol*16+x*16], &data[onesize*f], x, y, &resid[0]);
					Pre_hor_uv(&udata[y*ncol*4+x*4], &data[onesize*f+fsize], x, y, &resid[16]);
					Pre_hor_uv(&vdata[y*ncol*4+x*4], &data[onesize*f+fsize+uvsize], x, y, &resid[20]);
				}
				else if(bmode == 2){
					Pre_DC(&bdata[y*ncol*16+x*16], &data[onesize*f], x, y, &resid[0]);
					Pre_DC_uv(&udata[y*ncol*4+x*4], &data[onesize*f+fsize], x, y, &resid[16]);
					Pre_DC_uv(&vdata[y*ncol*4+x*4], &data[onesize*f+fsize+uvsize], x, y, &resid[20]);
				}
				else{
					if(f != 0){
						int mv_x, mv_y;
						Pre_ref(&bdata[y*ncol*16+x*16], &data[onesize*(f-1)], x, y, 0, &resid[0], mv_x, mv_y);
						Pre_ref_uv(&udata[y*ncol*4+x*4], &data[onesize*(f-1)+fsize], x, y, &resid[16], mv_x, mv_y);
						Pre_ref_uv(&vdata[y*ncol*4+x*4], &data[onesize*(f-1)+fsize+uvsize], x, y, &resid[20], mv_x, mv_y);
						output[osize++] = (char)mv_x;
						output[osize++] = (char)mv_y;
					}
					else{
						for(int i = 0;i < 16;++i){
							resid[i] = bdata[y*ncol*16+x*16+i];
						}
						for(int i = 0;i < 4;++i){
							resid[16+i] = udata[y*ncol*4+x*4+i];
						}
						for(int i = 0;i < 4;++i){
							resid[20+i] = vdata[y*ncol*4+x*4+i];
						}
					}
				}
				
				//DCT
				BYTE *residf = new BYTE[24];
//				DCT4(&resid[0], &residf[0]);
//				DCT2(&resid[16], &residf[16]);
//				DCT2(&resid[20], &residf[20]);
				
				//Quantization
				Quant(&resid[0], &resid[0], QP, 24);
				delete [] residf;
				
				//Zigzag
				BYTE *zdata = new BYTE[33];
				int zsize;
				zsize = Zig(&resid[0], zdata);
				
				for(int i = 0;i < zsize;++i){
					output[osize++] = zdata[i];
				}
				delete [] zdata;
				for(int i = 0;i < 8;++i){
					output[osize++] = resid[16+i];
				}
				delete [] resid;
			}
		}
	}
	delete [] bdata;
	delete [] udata;
	delete [] vdata;
	delete [] data;
	Output(file, output, osize, true);
	cout<<dataLen<<endl;
	cout<<osize<<endl;
	delete [] output;
	return;
}

void YUV_decompress(char file[]){
	BYTE *inputd = new BYTE[100000000];
	Input(file,inputd);
	int index = 0;
	int inbsize = (int)inputd[index++];
	int inwidth = (int)inputd[index++]*4;
	int inheight = (int)inputd[index++]*4;
	int infnum = (int)inputd[index++];
	int QP = (int)inputd[index++];
	int inLen = (inwidth*inheight + inwidth*inheight/2)*infnum;
	BYTE *data = new BYTE[inLen];
	int nrow = HEIGHT / 4;
	int ncol = WIDTH / 4;
	for(int f = 0;f < infnum;++f){
		cout<<f<<endl;
		BYTE *bdata = new BYTE[fsize];
		BYTE *udata = new BYTE[uvsize];
		BYTE *vdata = new BYTE[uvsize];
		for(int y = 0;y < nrow;++y){
			for(int x = 0;x < ncol;++x){
				int mode = (int)inputd[index++];
				int mv_x = 0, mv_y = 0;
				if(mode == 3 && f != 0){
					mv_x = (int)inputd[index++];
					mv_y = (int)inputd[index++];
				}
				//Inverse Zigzag
				BYTE *tempd = new BYTE[24];
				int zsize;
				zsize = Izig(&inputd[index], tempd);
				index += zsize;
				for(int i = 0;i < 8;++i){
					tempd[16+i] = inputd[index++];
				}
				
				//Inverse Quantization
				BYTE *tempf = new BYTE[24];
				IQuant(tempd, tempd, QP, 24);
				
				//IDCT
//				IDCT4(&tempf[0], &tempd[0]);
//				IDCT2(&tempf[16], &tempd[16]);
//				IDCT2(&tempf[20], &tempd[20]);
				delete [] tempf;
				
				//mode
				if(mode == 0){
					IPre_ver(&tempd[0], &data[onesize*f], x, y, &bdata[y*ncol*16+x*16]);
					IPre_ver_uv(&tempd[16], &data[onesize*f+fsize], x, y, &udata[y*ncol*4+x*4]);
					IPre_ver_uv(&tempd[20], &data[onesize*f+fsize+uvsize], x, y, &vdata[y*ncol*4+x*4]);
				}
				else if(mode == 1){
					IPre_hor(&tempd[0], &data[onesize*f], x, y, &bdata[y*ncol*16+x*16]);
					IPre_hor_uv(&tempd[16], &data[onesize*f+fsize], x, y, &udata[y*ncol*4+x*4]);
					IPre_hor_uv(&tempd[20], &data[onesize*f+fsize+uvsize], x, y, &vdata[y*ncol*4+x*4]);
				}
				else if(mode == 2){
					IPre_DC(&tempd[0], &data[onesize*f], x, y, &bdata[y*ncol*16+x*16]);
					IPre_DC_uv(&tempd[16], &data[onesize*f+fsize], x, y, &udata[y*ncol*4+x*4]);
					IPre_DC_uv(&tempd[20], &data[onesize*f+fsize+uvsize], x, y, &vdata[y*ncol*4+x*4]);
				}
				else{
					if(f != 0){
						IPre_ref(&tempd[0], &data[onesize*(f-1)], x, y, &bdata[y*ncol*16+x*16], mv_x, mv_y);
						IPre_ref_uv(&tempd[16], &data[onesize*(f-1)+fsize], x, y, &udata[y*ncol*4+x*4], mv_x, mv_y);
						IPre_ref_uv(&tempd[20], &data[onesize*(f-1)+fsize+uvsize], x, y, &vdata[y*ncol*4+x*4], mv_x, mv_y);
					}
					else{
						for(int i = 0;i < 16;++i){
							bdata[y*ncol*16+x*16+i] = tempd[i];
						}
						for(int i = 0;i < 4;++i){
							udata[y*ncol*4+x*4+i] = tempd[16+i];
						}
						for(int i = 0;i < 4;++i){
							vdata[y*ncol*4+x*4+i] = tempd[20+i];
						}
					}
				}
				delete [] tempd;
				
				//block to frame
				btof(bdata,&data[onesize*f]);
				btouv(udata, vdata, &data[onesize*f+fsize]);
			}
		}
		delete [] bdata;
		delete [] udata;
		delete [] vdata;
	}
	delete [] inputd;
	Output(file, data, inLen, false);
	delete [] data;
	return;
}

int main(){
	char filename1[30] = "AKIYO_352x288_10.yuv";
	char filename2[30] = "AKIYO_352x288_10_cp.yuv";
	//char filename2[30] = "MOBILE_352x288_10.yuv";
	YUV_compress(filename1, 4, 4);
	YUV_decompress(filename2);

	return 0;
}
