#include <iostream>
#include <libstructural.h>

using namespace std;
using namespace LIB_STRUCTURAL;


int main(int argc, char*argv[])
{
	if (argc != 2)
	{
		cerr << "Need one argument, full path to SBML file." << endl;
		return -1;
	}
	
	// get an instance of the library
	LibStructural* instance = LibStructural::getInstance();

	// print out model overview
	cout << instance->loadSBMLFromFile(argv[1]);
	
	// print test details
	cout << instance->getTestDetails();
	
	return 0;
}

// Running the program for BorisEJB.xml (available from the SBW 
// distribution) yields:
//
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//STRUCTURAL ANALYSIS MODULE : Results 
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//Size of Stochiometric Matrix: 8 x 10 (Rank is  5)
//Nonzero entries in Stochiometric Matrix: 20  (25% full)
//
//Independent Species (5) :
//MKK_P, MAPK_P, MKKK, MKK, MAPK
//
//Dependent Species (3) :
//MKK_PP, MKKK_P, MAPK_PP
//
//L0 : There are 3 dependencies. L0 is a 3x5 matrix.
//
//Conserved Entities
//1:  + MKK_P + MKK + MKK_PP
//2:  + MKKK + MKKK_P
//3:  + MAPK_P + MAPK + MAPK_PP
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//Developed by the Computational Systems Biology Group at Keck Graduate Institute 
//and the Saurolab at the Bioengineering Departmant at  University of Washington.
//Contact : Frank T. Bergmann (fbergman@u.washington.edu) or Herbert M. Sauro.   
//
//(previous authors) Ravishankar Rao Vallabhajosyula                   
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//
//Testing Validity of Conservation Laws.
//
//Passed Test 1 : Gamma*N = 0 (Zero matrix)
//Passed Test 2 : Rank(N) using SVD (5) is same as m0 (5)
//Passed Test 3 : Rank(NR) using SVD (5) is same as m0 (5)
//Passed Test 4 : Rank(NR) using QR (5) is same as m0 (5)
//Passed Test 5 : L0 obtained with QR matches Q21*inv(Q11)
//Passed Test 6 : N*K = 0 (Zero matrix)

