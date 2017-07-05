#include "qif"

#include <iostream>
#include <math.h>
#include <list>
#include <time.h>

using namespace qif;
using namespace std;

const static int n = 49;
const static double epsilon = 2.0;

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
//    cout<< "Priors:" << endl;
    for(int i = 0; i< n; i++)
    {
        Pi[i] = 1.0/n;
//        cout<< Pi[i] << endl;
    }

}

void initSetOfLocations(){

//The size of a grid
int size = sqrt(n);
//cout << "SIZE: " << size << ",
cout << "NUMBER OF POINTS:" << n << endl;
    
	for(int i = 0; i < size; i++){
		for(int j = 0; j < size; j++){
			listOfPoints[(size*i) + j].x = i;
			listOfPoints[(size*i) + j].y = j;
			//cout << listOfPoints[(size*i) + j].x << " " << listOfPoints[(size*i) + j].y << endl;
		}
	}

}


void assignVariables(){

  int number = 0;
  int counter = 0;
  int n2 = n*n;
  int n3 = n2*n;

    
    lp.c.set_size(n2);
    lp.b.set_size(n3 + n);
    lp.sense.set_size(n3 + n);
    
	for(int i = 0; i < n; i++){
		for(int j = 0; j < n; j++){

            lp.c(n*i + j) = Pi[i] * distance(listOfPoints[i],listOfPoints[j]);
            //cout << "vvvv:" << Pi[i] * distance(listOfPoints[i],listOfPoints[j]) << " :: "<<n*i + j<< endl;
			//entries.push_back(ME( n3 ,(n*i) + j, 1));
			for(int k = 0; k < n; k++){
                
                if(i!=k){
                    counter++;
                    if(distance(listOfPoints[i],listOfPoints[j]) < log(100.0))//Threshold
                    {
                        entries.push_back(ME((n2*k) + (n*i) + j, (n*i) + j, 1));
                    //cout << "EXP:" << distance(listOfPoints[i],listOfPoints[k]) << endl;
                        entries.push_back(ME((n2*k) + (n*i) + j, (n*k) + j, -exp(epsilon * distance(listOfPoints[i],listOfPoints[k]))));//entries.push_back(ME((n2*k) + (n*i) + j, (n*k) + j, -exp(distance(listOfPoints[i],listOfPoints[k]))));
                    }
                    else
                    {
                        number++;
                        //entries.push_back(ME((n2*k) + (n*i) + j, (n*i) + j, 0));
                        //entries.push_back(ME((n2*k) + (n*i) + j, (n*k) + j, 0));
                    }
                    lp.sense((n2*k)+ (n*i) + j) = '<';
                    lp.b((n2*k)+ (n*i) + j) = 0;
                }
                else
                {
                    counter++;
                    entries.push_back(ME((n2*k) + (n*i) + j, (n*i) + j, 0));
                    lp.sense((n2*k)+ (n*i) + j) = '=';
                    lp.b((n2*k)+ (n*i) + j) = 0;
                }
			}
			
		}
	}
    cout << "____number___ " << number << endl;
    for ( int i = n3; i < n3 + n ; i++) {
        counter++;
        for (int j = 0; j < n; ++j) {
            entries.push_back(ME( i , (n * (i - n3)) + j , 1));
        }
        lp.b(i) = 1;
        lp.sense(i) = '=';
    }
 
    cout << "NUMBER OF CONSTRAINTS: " << counter<< endl;
    cout << "VARIABLES: " << n2 << endl;
    lp.fill_A(entries);

}


void solve(){

    lp.maximize = false;
    // methods: simplex_primal, simplex_dual, interior
    lp.method = LinearProgram<double>::method_t::simplex_primal;
    cout <<" Method: "<< "interior" << endl;
    
    cout <<" A.n_rows: "<< lp.A.n_rows << endl;
    cout <<" b.n_elem: "<< lp.b.n_elem << endl;
    cout <<" sense.n_elem: "<< lp.sense.n_elem << endl;
    cout <<" A.n_cols: "<< lp.A.n_cols << endl;
    cout <<" c.n_elem: "<< lp.c.n_elem << endl;
    

    clock_t tStart = clock();
    bool solved = lp.solve();
    printf(" Time: %.2fs\n", (double)(clock() - tStart)/CLOCKS_PER_SEC);
    cout <<" Solved: " << solved << endl;
    cout <<" Status: "<< lp.status << endl;
    cout <<" lp.x.n_elem: "<< lp.x.n_elem <<endl;
    //cout << lp.A << endl;
    //cout << "Solve:" << solved << endl;
    //cout << endl;
    
    //cout <<" lp.x.at(5):"<< lp.x.at(5) << endl;
    //cout <<" lp.x.at(5):"<< lp.x.at(345) << endl;
    
    cout << " Optimum: "<< lp.optimum() << endl;

}


int main() {
   
    cout<<"-----------++++++++++++++++----------------++++++++++++++++---------------------"<<endl;

    initSetOfLocations();
    initPriors();
    assignVariables();
    solve();

    cout<<"END!"<<endl;
    cout<<"================================================================================"<<endl;

}

