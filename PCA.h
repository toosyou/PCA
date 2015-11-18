#ifndef PCA_H
#define PCA_H

#include <iostream>
#include <cmath>

using namespace std;

template <class T>
T** toCreateJacobian(T** Input,const int N,const int i,const int j,const T rPhi)
{
    T rSp = sin(rPhi);
    T rCp = cos(rPhi);
    T** Temp = new T*[N];
    for(int ii=0;ii<N;++ii)
        Temp[ii] = new T[N];
    for(int ii=0;ii<N;++ii)
        for(int jj=0;jj<N;++jj)
        {
            Temp[ii][jj]=Input[ii][jj];
            Input[ii][jj]=(T)0;
        }

    for(long ii=0; ii<N; ++ii){
        for(long jj=0; jj<N; ++jj){
            if( ii==i ){     // row i
                if( jj==i ) Input[ii][jj] = Temp[i][i]*rCp*rCp + Temp[j][j]*rSp*rSp + 2*Temp[i][j]*rCp*rSp;
                else if( jj==j ) Input[ii][jj] = (Temp[j][j]-Temp[i][i])*rSp*rCp + Temp[i][j]*(rCp*rCp-rSp*rSp);
                else Input[ii][jj] = Temp[i][jj]*rCp + Temp[j][jj]*rSp;
            } else if ( ii==j ) {// row j
                if( jj==i ) Input[ii][jj] = (Temp[j][j]-Temp[i][i])*rSp*rCp + Temp[i][j]*(rCp*rCp-rSp*rSp);
                else if( jj==j ) Input[ii][jj] = Temp[i][i]*rSp*rSp + Temp[j][j]*rCp*rCp - 2*Temp[i][j]*rCp*rSp;
                else Input[ii][jj] = Temp[j][jj]*rCp - Temp[i][jj]*rSp;
            } else {            // row l ( l!=i,j )
                if( jj==i ) Input[ii][jj] = Temp[i][ii]*rCp + Temp[j][ii]*rSp;
                else if( jj==j ) Input[ii][jj] = Temp[j][ii]*rCp - Temp[i][ii]*rSp;
                else Input[ii][jj] = Temp[ii][jj];
            }
        }
    }

    for(int ii=0;ii<N;++ii)
        delete [] Temp[ii] ;
    delete [] Temp;

    return Input;

}


template <class T>
T getMax(T** Input,const int N,int &nRow,int &nCol)
{
    T rMax=Input[0][1];
    nRow = 0;
    nCol = 1;

    for(int i=0;i<N;++i)
        for(int j=0;j<N;++j)
            if(i!=j)
                if( abs(Input[i][j]) > rMax )
                {
                    rMax = abs(Input[i][j]);
                    nRow = i ;
                    nCol = j ;
                }
    return rMax;
}


//Data[N][M] ; Return EigenMatrix[M][M+1] (Last Column of Matrix is for EigenValue
//If the data are formed by points Points_Vectors = true , as to by vectors = false
template <class T>
T** PCA(T** Data,const int M,const int N, const bool Points_Vectors = true ,const T Err=0.00001)
{
    T*  Average     = new T[M] ;
    T** TempData    = new T*[N]; //[N][M]
    T** S           = new T*[M]; //[M][M]
    T** EigenVector = new T*[M]; //[M][M+1]
    T*  EigenValue  = new T[M];
    T** eTemp       = new T*[M]; //[M][M]
    T** eVec        = new T*[M]; //[M][M]
    T** eC          = new T*[M]; //[M][M]

    //Creat TempData[N][M] , S[M][M] , EigenVector[M][N] ,eTemp[M][M] , eVec[M][M],eC[M][M]
    for(int i=0;i<N;++i)
        TempData[i] =new T[M];
    for(int i=0;i<M;++i)
    {
        S[i]            = new T[M];
        EigenVector[i]  = new T[M+1];
        eTemp[i]        = new T[M];
        eVec[i]         = new T[M];
        eC[i]           = new T[M];
    }

    //Reset
    for(int i=0;i<M;++i)
        Average[i]=(T)0;
    for(int i=0;i<N;++i)
        for(int j=0;j<M;++j)
            TempData[i][j]=Data[i][j];
    for(int i=0;i<M;++i)
        for(int j=0;j<M;++j)
            S[i][j]=EigenVector[i][j+1]=eVec[i][j]=eC[i][j]=(T)0;
    for(int i=0;i<M;++i)
        for(int j=0;j<M;++j)
        {
            eTemp[i][j]=(T)0;
            if(i==j)
                eTemp[i][j]=(T)1;
        }
    for(int i=0;i<M;++i)
        EigenValue[i]=(T)0;


    if(Points_Vectors)
    {
        //Calculate Average Point
        for(int i=0;i<N;++i)
            for(int j=0;j<M;++j)
                Average[j]+=Data[i][j]/(T)N;

        //Calculate TempData
        for(int i=0;i<N;++i)
            for(int j=0;j<M;++j)
                TempData[i][j] -= Average[j];
    }

    //Calculate S
    for(int i=0;i<M;++i)
        for(int j=0;j<M;++j)
            for(int k=0;k<N;++k)
                S[i][j]+=TempData[k][i] * TempData[k][j] / (T)N;

    //Calculate EigenVector and EigenValue of S

    while(1)
    {
        int i=0,j=0;
        T rMax = getMax(S,M,i,j);
        if(rMax <= Err) break ;

        T rPhi = atan2((T)2*S[i][j], S[i][i] - S[j][j]) / (T)2;
        S = toCreateJacobian(S,M,i,j,rPhi);
        for(int x=0;x<M;++x)
            eC[x][x] = (T)1;
        eC[j][j] = eC[i][i] = cos(rPhi);
        eC[j][i] = sin(rPhi);
        eC[i][j] = -eC[j][i];

        for(int x=0; x<M; ++x) // for eigenvectors
            for(int y=0; y<M; ++y)
                for(int z=0; z<M; ++z)
                    eVec[x][y] = eVec[x][y] + eTemp[x][z] * eC[z][y];

        for(int x=0; x<M; ++x)
            for(int y=0; y<M; ++y)
            {
                eTemp[x][y] = eVec[x][y];
                eVec[x][y]  = (T)0;
                eC[x][y]    = (T)0;
            }
    }
    for(int i=0;i<M;++i)
        EigenValue[i] = EigenVector[i][M] = S[i][i];
    for(int i=0;i<M;++i)
        for(int j=0;j<M;++j)
            EigenVector[i][j] = eTemp[j][i];

    //Sort the EigenVector order of the EigenValue
    T* TempEigenVector = new T[M];

    for(int i=0;i<M;++i)
        for(int j=i;j<M-1;++j)
            if(EigenValue[j]<EigenValue[j+1])
            {
                T TempEigenValue = EigenValue[j+1];
                EigenValue[j+1]=EigenVector[j+1][M]=EigenValue[j];
                EigenValue[j]  =EigenVector[j][M]  =TempEigenValue;
                for(int k=0;k<M;++k)
                    TempEigenVector[k]  =EigenVector[j+1][k];
                for(int k=0;k<M;++k)
                    EigenVector[j+1][k] =EigenVector[j][k];
                for(int k=0;k<M;++k)
                    EigenVector[j][k]   =TempEigenVector[k];
            }

    delete [] TempEigenVector;

    //Delete All the stuffs
    delete Average ;
    delete EigenValue ;

    for(int i=0;i<N;++i)
        delete [] TempData[i];
    delete [] TempData;

    for(int i=0;i<M;++i)
    {
        delete [] S[i];
        delete [] eTemp[i];
        delete [] eVec[i];
        delete [] eC[i];
    }
    delete [] S;
    delete [] eTemp;
    delete [] eVec;
    delete [] eC;

    return EigenVector;


}

#endif // PCA_H
