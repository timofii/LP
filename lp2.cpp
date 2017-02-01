#include "qif"

#include <iostream>
#include <math.h>
#include <list>

using namespace qif;
using namespace std;

const static int n = 4;

typedef qif::MatrixEntry<double> ME;

struct point2 {
		int x;
		int y;
		};

static double Pi[n];
static point2 listOfPoints[n];

static list<ME> entries; 
static LinearProgram<double> lp;


double distance(point2 i, point2 j){
	
	return sqrt(((i.x - j.x)*(i.x - j.x)) + ((i.y - j.y)*(i.y - j.y)));

}

void initPriors(){
    cout<< "Priors:" << endl;
    for(int i = 0; i< n; i++)
    {
        Pi[i] = 1.0/n;
        cout<< Pi[i] << endl;
    }

}

void initSetOfLocations(){

//The size of a grid
int size = sqrt(n);
cout << size << endl;
	for(int i = 0; i < size; i++){
		for(int j = 0; j < size; j++){
			listOfPoints[(size*i) + j].x = i;
			listOfPoints[(size*i) + j].y = j;
			cout <<listOfPoints[(size*i) + j].x<<" "<<listOfPoints[(size*i) + j].y<< endl; 
		}
	}

}


void assignVariables(){

  int n2 = n*n;
  int n3 = n2*n;

  lp.c.set_size(n2);
  lp.b.set_size(n3 + 1);  
  lp.sense.set_size(n3 + 1);
/*
  point2 p1, p2; //, p3;
  p1.x=0; p2.x=3; //p3.x=4;
  p1.y=0; p2.y=4; //p3.y=3;
  listOfPoints[0]=p1; listOfPoints[1]=p2;// listOfPoints[2]=p3;
*/
	for(int i = 0; i < n; i++){
		for(int j = 0; j < n; j++){

			lp.c(n*i + j) = Pi[i] * distance(listOfPoints[i],listOfPoints[j]);
			entries.push_back(ME( n3 ,(n*i) + j, 1));
			for(int k = 0; k < n; k++){
				if(i!=k){
					entries.push_back(ME((n2*k) + (n*i) + j, (n*i) + j, 1));
					entries.push_back(ME((n2*k) + (n*i) + j, (n*k) + j, -exp(distance(listOfPoints[i],listOfPoints[k]))));
				}	
				lp.sense((n2*k)+ (n*i) + j) = '<';
				lp.b((n2*k)+ (n*i) + j) = 0;
			
			}
			
		}
	}
 
  lp.b(n3) = 1;
  lp.sense(n3) = '=';

  lp.fill_A(entries);

}


void solve(){

lp.method = LinearProgram<double>::method_t::simplex_primal;
bool solved = lp.solve();

cout << lp.A << endl;
cout << lp.status << endl;
cout << "Solve:" << solved << endl;
    cout << endl;
    cout << lp.x;

}


int main() {

cout<<"Launching function main"<<endl;

initSetOfLocations();
initPriors();
assignVariables();
solve();

cout<<"END!"<<endl;

}

