# PCA
Principal component analysis in C++  
wiki : https://en.wikipedia.org/wiki/Principal_component_analysis

##Usage
###input
* `template <class T>`
    * T could be float or double
* `const int N, M`
    * size of data
    * N : number of vector
    * M : dimension of vector
* `T**Data`
    * size : [N][M]
* `const bool Points_Vectors`
    * `true` : the data are formed by coordinates of points
        * this will cause a shift to the center of data
    * `false` : the data are formed by vectors
* `const T Err`
    * precision of answer
    * default with 0.00001
        * OK with T = `float`

###output
* return : pca result
    * `T** EigenMatrix`
        * size : [M][M+1]
        * [Eigenvector[0],Eigenvector[1],...,Eigenvector[m-1],Eigenvalue] * M
        * sorted by eigenvalues

###test
* original data : 
* expected result :

##Contact Me
###Me
A silly student of NCTU. Coding for life :)
###Email
ycc.cs03@nctu.edu.tw
