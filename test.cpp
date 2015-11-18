#include <iostream>
#include "PCA.h"
#include <stdio.h>
#include <stdlib.h>

using namespace std;
int main()
{
	int M=2, N=10;
	double** data = new double*[N];
	for(int i=0;i<N;i++)
		data[i] = new double[M];
	
	data[0][0]=1.5;  data[0][1]=2.3;
	data[1][0]=3.0;  data[1][1]=1.7; 
	data[2][0]=1.2;  data[2][1]=2.9; 
	data[3][0]=2.1;  data[3][1]=2.2; 
	data[4][0]=3.1;  data[4][1]=3.1; 
	data[5][0]=1.3;  data[5][1]=2.7; 
	data[6][0]=2.0;  data[6][1]=1.7; 
	data[7][0]=1.0;  data[7][1]=2.0; 
	data[8][0]=0.5;  data[8][1]=0.6; 
	data[9][0]=1.0;  data[9][1]=0.9;
		
//	for(int i=0;i<N+1;i++) result[i][0]=result[i][1]=0;
		
	double** result=PCA(data,M,N,true); 
	
	cout<<result[0][0]<<" "<<result[0][1]<<" "<<result[0][2]<<endl;
	cout<<result[1][0]<<" "<<result[1][1]<<" "<<result[1][2]<<endl;
	
	for(int i=0;i<M;++i)
		delete [] result[i];
	delete [] result;
	return 0;
}
