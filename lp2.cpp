#include "qif"

#include <iostream>
#include <math.h>
#include <list>
#include <time.h>

using namespace qif::lp;
using namespace std;

const static int n = 49;
const static int n2 = n*n;
const static double epsilon = 1.0;

typedef MatrixEntry<double> ME;

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

int counter = 0;
    lp.c.set_size( n2 );
    
	for(int i = 0; i < n; i++){
		for(int j = 0; j < n; j++){
        
            lp.c( n*i  + j) = Pi[i] * distance(listOfPoints[i],listOfPoints[j]);
            for(int k = 0; k < n; k++){
                if(i!=k){
                    //if((epsilon * distance( listOfPoints[i],listOfPoints[j])) < std::log(1e200))
                    //{
                        counter++;
                        entries.push_back(ME(counter, (n*i) + j, 1));
                        entries.push_back(ME(counter, (n*k) + j, -exp(epsilon * distance(listOfPoints[i],listOfPoints[k]))));
                    //}
                    
                }
			}
			
		}
	}
    
    for ( int i = counter; i < counter + n ; ++i) {
        for (int j = 0; j < n; ++j) {
            entries.push_back(ME( i , (n * (i - counter)) + j , 1));
        }
    }

    lp.b.set_size(counter + n);
    lp.sense.set_size(counter + n);
    
    for (int i = 0; i < counter; ++i) {
        lp.sense(i) = '<';
        lp.b(i) = 0;
    }
    
    for (int i = counter; i < counter + n; ++i) {
        lp.b(i) = 1;
        lp.sense(i) = '=';
    }
    
    lp.fill_A(entries);

}


void solve(){

    cout << " TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT " << endl;
    lp.maximize = false;
    // methods: simplex_primal, simplex_dual (default), interior
    lp.method = method_t::simplex_primal;
	lp.glp_msg_level = msg_level_t::all;
	lp.glp_presolve = true;
    cout << "===================================" << endl;
    cout <<" Method: "<< lp.method << endl;
    cout << "===================================" << endl;
    
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
    cout << " TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT " << endl;


//	qif::chan C(rows, cols);


}

int kostas() {
	std::cout << "KKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKK" << endl;
	
	Defaults::glp_msg_level = msg_level_t::all;
		//Timo: set the method for solver
	Defaults::method = method_t::simplex_primal;
	std::cout << "===================================" << endl;
    	std::cout <<" Method: "<< Defaults::method << endl;
    	std::cout << "===================================" << endl;
	
	uint width = std::sqrt(n);

	auto d_euclid = qif::metric::grid<double>(width);
	auto d_priv = epsilon * d_euclid;
	auto d_loss = d_euclid;
	double threshold = 1e200;

	qif::prob pi = qif::probab::uniform<double>(n);

	qif::chan opt = qif::mechanism::optimal_utility(pi, n, d_priv, d_loss, threshold);

	std::cout << "size: " << opt.n_elem << "\n";
	std::cout << "Pr(0 | 0) = " << opt.at(0, 0) << "\n";
	std::cout << "utility: " << qif::utility::expected_distance(d_loss, pi, opt) << "\n";
	cout << "actual epsilon: " << qif::mechanism::smallest_epsilon(opt, d_priv) << "\n";
	cout << "is proper: " << qif::channel::is_proper(opt) << "\n";


    	std::cout << "KKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKK" << endl;
	return 0;
}


int main() {
return kostas();
   
    cout<<"-----------++++++++++++++++----------------++++++++++++++++---------------------"<<endl;

    initSetOfLocations();
    initPriors();
    assignVariables();
    solve();

    cout<<"END!"<<endl;
    cout<<"================================================================================"<<endl;

}

