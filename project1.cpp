#include <iostream>
#include <armadillo>
#include <vector>
#include <cmath>
#include <fstream>
#include <string>
#include <sstream>

using namespace std;
using namespace arma;

double g(double x){

    return 100. * exp(-10 * x);

}


void relative_error(vec &v, vec &v_exact, int n){

    ofstream myfile;
    myfile.open("relative_error.dat", ios_base::app);

    double error = 0;
    double tmp_error = 0;
    for(int i = 2; i < n+1; i++){
        tmp_error = abs(((v(i) - v_exact(i))))/abs(v_exact(i));
        if (tmp_error > error) {
            error = tmp_error;
        } 
    }


    myfile << n << " " << error << "\n";
    myfile.close();
}  


vec own_general_solver(vec &a, vec &b, vec &c, vec &f, int n){ 
    vec v = zeros<vec>(n + 2);  //the unknown vector
    vec b_tilde(n + 2); //don't need this as vector because only interested in V  
    vec f_tilde(n + 2); //updated f after row reduction
    double tmp; //temporary value used for row reduction

    a(1) = 0; //outside matrix
    c(n) = 0; //outside matrix
    b_tilde(1) = b(1);//for first iteration 
    f_tilde(1) = f(1);

    /*Foreward substitution */
    for(int i = 2; i < n+1; i++){
        tmp = a(i)/b_tilde(i-1); 
        b_tilde(i) = b(i) - tmp*c(i-1);//updated mid diagonal 
        f_tilde(i) = f(i) - tmp*f_tilde(i-1); 
    }

    v(n) = f_tilde(n)/b_tilde(n); 

    /*Backward substitution*/
    for(int i = n-1; i >= 1; i--){
        v(i) = (f_tilde(i) - c(i)*v(i+1))/b_tilde(i);
    }
      
    return v;

}  

vec special_solver(vec &f, int n){ 
    
    vec v = zeros<vec>(n+2); //the unknown vector, TODO change name from v when not dead tired. 
    vec f_tilde = zeros<vec>(n+2);   
    vec b_tilde = zeros<vec>(n+2);  

    b_tilde(1) = 2.;
    f_tilde(1) = f(1);

    /*foreward substitution*/
    for(int i = 2; i < n+1; i++){
        b_tilde(i) = (i+1.)/i;
        f_tilde(i) = f(i) + f_tilde(i-1)/b_tilde(i-1); 
    }             
         
    v(n) = f_tilde(n)/b_tilde(n); 

    for(int i = n-1; i >= 1; i--){
        v(i) = (f_tilde(i) + v(i+1))/ b_tilde(i);
    }
    return v;
}      

vec lu_solver(mat A,vec f){
    mat L, U, P;
    lu(L,U, P, A);
    vec v = solve(trimatu(U), solve(trimatl(L), P*f));

    return v;
}


/*writing data to file */
void write_file(string filename, vec x, vec v, vec u, int n){

    ofstream myfile;
    myfile.open(filename);

    for(int i = 0; i < n+2; i++){

        myfile << x(i) <<" " << v(i) <<" " << u(i)<< " \n";
    }
    myfile.close();
}

int main(int argc, char* argv[]){

    int n; // array size, will be taken from command line
    int decide_solver; //decides which solver to use

    /*checking for at least one command-line arguments */
    if(argc > 2) {
        n = atoi(argv[1]); 
        decide_solver = atoi(argv[2]);
    } else {
        cout << "you need to give the array size" <<endl;
        exit(1);
    }     

    int N = 1; //end value
    double h = N/(n+1.0); //step length
    
    //allocate three diagonals
    
    vec a = zeros<vec>(n+2); //bottom diagonal
    vec b = zeros<vec>(n+2); //middle diagonal
    vec c = zeros<vec>(n+2); //top diagonal
   
    //filling the three diagonal vectors
    a.fill(-1.);
    b.fill(2.);
    c.fill(-1.); 
    
    //setting up vectors for x and f
    vec x = zeros<vec>(n + 2); //x-grid
    vec f = zeros <vec>(n + 2); // the known vector (called b_tilde in exercise)
    vec u = zeros <vec>(n + 2); // closed form solution
   

    for(int i = 0; i < n + 2; i++){

        x(i) = i*h;
        f(i) = g(x(i)) * h * h;
        u(i) = 1 - (1 - exp(-10))*x(i)- exp(-10*x(i));
    }  

    switch(decide_solver) {
        case 1:{
               cout<<"using my own general solver"<<endl;

               clock_t start, finish; //declare start and finish time
               start = clock();

               vec v_general = own_general_solver(a, b, c, f, n);  

               finish = clock();

               cout << ((double) (finish - start) /((double)CLOCKS_PER_SEC)) << endl;

               ostringstream os;
               os << "ownsolver_n_" << n << ".dat";
               string filename = os.str();
               write_file(filename, x, v_general, u,  n);
               relative_error(v_general, u, n);
               break;
           }
        case 2:{
                cout<<"own solver for special case problem"<<endl;

                clock_t start, finish; //declare start and finish time
                start = clock();

                vec v_special = special_solver(f, n); 

                finish = clock();


                cout << ((double) (finish - start) /((double)CLOCKS_PER_SEC)) << endl;
                                    
                ostringstream os;
                os << "specialsolver_n_" << n << ".dat";
                string filename = os.str();
                write_file(filename, x, v_special, u, n);

                break;
             }
        case 3:{
                cout<<"LU-decom>p to file here!"<<endl;
                mat A = zeros(n + 2, n + 2);
                A.diag(-1) += -1; 
                A.diag(0) += 2;
                A.diag(1) += -1;

                clock_t start, finish; //declare start and finish time
                start = clock();     

                vec v_lu = lu_solver(A, f);

                finish = clock();  

                cout << ((double) (finish - start) /((double)CLOCKS_PER_SEC)) << endl;
               break;
           }
    }
}

