#include "FCrand.h"
#include "Tools1.h"

using namespace std;


int main(){
	int N=30;
	//double x[N];
	int x[N];

	FCrand R;

	for(int i=0; i<N; i++){
		x[i] = R.GetRandom(0,N);
		cout << x[i] << endl;
	}
	cout << "/////////////////////////////////////" << endl;

	clock_t tStart = clock();
	quicksort(x,0,N-1);
	printf("Time taken to run this code: %.7fs\n", (double)(clock() - tStart)/CLOCKS_PER_SEC);

	for(int i=0; i<N; i++)
		cout << x[i] << endl;


	//binary search
	tStart = clock();
	int vintesete = binarySearch(x,0,N-1,27);
	printf("Time taken to run this code: %.7fs\n", (double)(clock() - tStart)/CLOCKS_PER_SEC);

	if(vintesete!=-1)
		cout << "The number 27 was found in index " << vintesete << 
		" [proof: " << x[vintesete] << "]" << endl;

	return 0;
}