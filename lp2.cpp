#include "qif"

#include <iostream>
#include <math.h>
#include <list>
#include <time.h>

using namespace qif::lp;
using namespace std;

const static int n = 64;
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

double computing_delta(double ro, double R)
{
    if(ro != 0.0){
        cout << "R = " << R << " ro =" << ro << endl;
        return (R/ro)/((R/ro) + 2) ;
        //return ((R/ro) + 2) / (R/ro) ;
    }
    else{
        cout << "Can not compute delta: ro = 0" << endl;
        return 0.0;
    }
}

void initPriors(){
    for(int i = 0; i< n; i++){
        Pi[i] = 1.0/n;
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

//variable KostasMethod make a switch between Direct and Kostas method with parameters ro and R
//true to run Kostas method false to run Direct
//As Catuscia mentioned I have put delta = sqrt(2)/2 in her email to all of us
//before was computation of delta according the formula therefore parameter ro is not in use 
void assignVariables(bool KostasMethod , double ro, double R){

    int counter = 0;
    lp.c.set_size( n2 );

    double delta = sqrt(2)/2;//computing_delta(ro,R);
    cout << "DELTA * EPS: "<< delta * epsilon << endl;
    
    for(int i = 0; i < n; i++){
        for(int j = 0; j < n; j++){

            lp.c( n*i  + j) = Pi[i] * distance(listOfPoints[i],listOfPoints[j]);
            for(int k = 0; k < n; k++){
                if(i!=k){
                    //if((epsilon * distance( listOfPoints[i],listOfPoints[j])) < std::log(1e200))
                    //{
                    if( KostasMethod )
                    {
                        if(distance(listOfPoints[i],listOfPoints[k]) <= R){
                            entries.push_back(ME(counter, (n*i) + j, 1));
                            entries.push_back(ME(counter, (n*k) + j, -exp( delta * epsilon * distance(listOfPoints[i],listOfPoints[k]))));
                            counter++;
                        }
                    }
                    else
                    {
                        entries.push_back(ME(counter, (n*i) + j, 1));
                        entries.push_back(ME(counter, (n*k) + j, -exp(epsilon * distance(listOfPoints[i],listOfPoints[k]))));
                        counter++;
                    }
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
    
    lp.maximize = false;
    // methods: simplex_primal, simplex_dual (default), interior
    lp.method = method_t::simplex_dual;
    lp.glp_msg_level = msg_level_t::all;
    lp.glp_presolve = true;
    cout << "===================================" << endl;
    cout <<"Method: "<< lp.method << endl;
    cout << "===================================" << endl;
    
    cout <<" CONSTRAINTS: "<< lp.A.n_rows << endl;
    cout <<" VARIABLES: "<< lp.A.n_cols << endl;
    //cout <<" c.n_elem: "<< lp.c.n_elem << endl;
    //cout <<" b.n_elem: "<< lp.b.n_elem << endl;
    //cout <<" sense.n_elem: "<< lp.sense.n_elem << endl;
    
    cout << endl;
    clock_t tStart = clock();
    bool solved = lp.solve();
    printf(" Time: %.2fs\n", (double)(clock() - tStart)/CLOCKS_PER_SEC);
    cout <<" Solved: " << solved << endl;
    cout <<" Status: "<< lp.status << endl;
    
    //cout <<" lp.x.n_elem: "<< lp.x.n_elem <<endl;
    //cout << lp.A << endl;
    //cout << "Solve:" << solved << endl;
    //cout << endl;
    
    //Only to check the correctness of the optimal solution
    qif::chan opt(n,n);
    for(int i = 0; i < n; ++i){
        for(int j = 0; j < n; ++j){
            opt(i,j) = lp.x.at( i * n + j);
        }
    }
    
    //add probability distr that kostas use
    //add parameters that kostas use
    uint width = std::sqrt(n);
    auto d_euclid = qif::metric::grid<double>(width);
    auto d_priv = epsilon * d_euclid;
    auto d_loss = d_euclid;
    qif::prob pi = qif::probab::uniform<double>(n);
    
    cout << endl << "CHECK THE CORRECTNESS" << endl;
    cout << "Optimum(c,x): "<< lp.optimum() << endl;
    std::cout << "utility: " << qif::utility::expected_distance(d_loss, pi, opt) << "\n";
    cout << "epsilon in a programm: " << epsilon << endl;
    cout << "actual epsilon: " << qif::mechanism::smallest_epsilon(opt, d_priv) << "\n";
    cout << "is proper: " << qif::channel::is_proper(opt) << "\n";
    
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
    
    //std::cout << "size: " << opt.n_elem << "\n";
    //std::cout << "Pr(0 | 0) = " << opt.at(0, 0) << "\n";
    
            //////////////
    std::cout << "utility: " << qif::utility::expected_distance(d_loss, pi, opt) << "\n";
    cout << "actual epsilon: " << qif::mechanism::smallest_epsilon(opt, d_priv) << "\n";
    cout << "is proper: " << qif::channel::is_proper(opt) << "\n";
    
    std::cout << "KKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKK" << endl;
    return 0;
}


int main() {
    //return kostas();
    
    cout<<"-----------++++++++++++++++----------------++++++++++++++++---------------------"<<endl;
    
    initSetOfLocations();
    initPriors();
    assignVariables(true, 0.0 ,6.5);//,4);//////
    
    solve();
    
    cout<<"END!"<<endl;
    cout<<"================================================================================"<<endl;
    
}

