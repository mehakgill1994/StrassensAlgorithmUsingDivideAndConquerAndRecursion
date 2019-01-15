#include <iostream>
#include <fstream>
#include "simple_matrix.hpp"

using namespace std;

//function to ADD two matrices
template <class T>
SMatrix<T> add(SMatrix<T> P, SMatrix<T> Q)
{
	int s = P.size();
	SMatrix<T> result(s);
	for(int i=0; i<s; i++)
	{
		for(int j=0; j<s; j++)
		{
			result(i,j) = P(i,j) + Q(i,j);
		}
	}
	return result;
}

//function to SUBTRACT two matrices
template <class T>
SMatrix<T> subtract(SMatrix<T> P, SMatrix<T> Q)
{
	int s = P.size();
	SMatrix<T> result(s);
	for(int i=0; i<s; i++)
	{
		for(int j=0; j<s; j++)
		{
			result(i,j) = P(i,j) - Q(i,j);
		}
	}
	return result;
}

//function to divide any matrix M into a square submatrix staring from indices (startx, starty) to (endx, endy)
template <class T>
SMatrix<T> partition(SMatrix<T> M, int startx, int starty, int endx, int endy)
{
	int s = M.size();
	SMatrix<T> P (s/2);
	for(int i=startx; i<endx; i++)
	{
		for(int j=starty; j<endy; j++)
		{
			P(i-startx,j-starty) = M(i,j);
		}	
	}	
	return P;		
}

//function to join 4 submartices of size 'n/2' to make a single matrix of size 'n'
template <class T>
SMatrix<T> combine(SMatrix<T> P, SMatrix<T> Q, SMatrix<T> R, SMatrix<T> S)
{
	int s = P.size();
	SMatrix<T> M(2*s);
	for(int i=0; i<s; i++)
	{
		for(int j=0; j<s; j++)
		{
			M(i,j) = P(i,j);
			M(i,j+s) = Q(i,j);
			M(i+s,j) = R(i,j);
			M(i+s,j+s) = S(i,j);
		}	
	}	
	return M;		
}

//function that implements strassen's matrix multiplication using divide and conquer 
template <class T>
SMatrix<T> function(SMatrix<T> A,SMatrix<T> B)
{
	int n = A.size();
	int r;
	SMatrix<T> c(n), S1(n/2), S2(n/2), S3(n/2), S4(n/2), S5(n/2), S6(n/2), S7(n/2), S8(n/2), S9(n/2), S10(n/2), P1(n/2), P2(n/2), P3(n/2), P4(n/2), P5(n/2), P6(n/2), P7(n/2), c00(n/2), c01(n/2), c10(n/2), c11(n/2);	
	if(A.size() == 1)
	{
		c(0,0) = A(0,0) * B(0,0);
	}	
	else
	{
		/*
		The 4 partitions of matrix M are as follows:
		M11 -> (M, 0, 0, n/2, n/2)
		M12 -> (M, 0, n/2, n/2, n)
		M21 -> (M, n/2, 0, n, n/2)
		M22 -> (M, n/2, n/2, n, n)
		*/
		S1 = subtract(partition(B, 0, n/2, n/2, n), partition(B, n/2, n/2, n, n));
		S2 = add(partition(A, 0, 0, n/2, n/2), partition(A, 0, n/2, n/2, n));
		S3 = add(partition(A, n/2, 0, n, n/2), partition(A, n/2, n/2, n, n));
		S4 = subtract(partition(B, n/2, 0, n, n/2), partition(B, 0, 0, n/2, n/2));
		S5 = add(partition(A, 0, 0, n/2, n/2), partition(A, n/2, n/2, n, n));
		S6 = add(partition(B, 0, 0, n/2, n/2), partition(B, n/2, n/2, n, n));
		S7 = subtract(partition(A, 0, n/2, n/2, n), partition(A, n/2, n/2, n, n));
		S8 = add(partition(B, n/2, 0, n, n/2), partition(B, n/2, n/2, n, n));
		S9 = subtract(partition(A, 0, 0, n/2, n/2), partition(A, n/2, 0, n, n/2));
		S10 = add(partition(B, 0, 0, n/2, n/2), partition(B, 0, n/2, n/2, n));
		P1 = function(partition(A, 0, 0, n/2, n/2), S1);
		P2 = function(S2, partition(B, n/2, n/2, n, n));
		P3 = function(S3, partition(B, 0, 0, n/2, n/2));
		P4 = function(partition(A, n/2, n/2, n, n), S4);
		P5 = function(S5, S6);
		P6 = function(S7, S8);
		P7 = function(S9, S10);
		
		c00 = add(subtract(add(P5, P4), P2), P6);
		c01 = add(P1, P2);
		c10 = add(P3, P4);
		c11 = subtract(subtract(add(P5, P1), P3), P7);
		
		c = combine(c00, c01, c10, c11);
	}
	return c;
}


int main(int argc, char* argv[])
{
	ifstream fi(argv[1]);
	if (!fi)
	{
		cout << "cant open fi" << endl;
		return -1;
	}

	ofstream fo(argv[2]);
	if (!fo)
	{
		cout << "cant open fo" << endl;
		return -1;
	}

	int m, n, i, count;
	
	//reading the first value (i.e. the number of datasets) in the input set
	fi >> m;
	
	//check if the input file is empty
	if (m == 0)
	{
		cout << "the file has no dataset";
	}
	else
	{
		for (int j = 1; j <= m; j++)
		{
			//reading the size of the matrices to multiplied
			fi >> n;
			if (n == 0)
			{
				fo << "EMPTY" << endl;
			}
			else
			{
				SMatrix<float> A(n), B(n), C(n);
				
				//reading in the first matrix
				for(int k=0; k<n; k++)
				{
					for(int l=0; l<n; l++)
					{
						fi>>A(k,l);
					}
				}
				
				//reading in the second matrix
				for(int k=0; k<n; k++)
				{
					for(int l=0; l<n; l++)
					{
						fi>>B(k,l);
					}
				}
					
				//Strassen's Algorithm			
				C = function(A,B);

				//writing the solution to the output file
				for(int k=0; k<n; k++){
					for(int l=0; l<n; l++){
						fo<<C(k,l);
						if(k+l != 2*n-2)
							fo<<" ";							
					}
				}
				fo<<endl;
			}
		}
	}

return 1;
}