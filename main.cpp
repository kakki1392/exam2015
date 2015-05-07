#include "koch.h"
#include <iostream>
#include <vector>
#include <unistd.h>
#include <string>
#include <sstream>
using namespace arma;
using namespace std;

int main(){
	cout.precision(7);

	Koch k;
	/*
	k.gplt.cmd("set term pdfcairo");
	k.gplt.cmd("set output 'data/koch3.pdf'");
	k.gplt.cmd("set size square");
	k.gplt.cmd("set xlabel 'x'");
	k.gplt.cmd("set ylabel 'y'");
	k.gplt.cmd("set title 'l = 3, L = 1'");
	*/
	//k.gplt.cmd("set size square");

	k.fill_x();
	//k.x.print("x: ");
	//cout << k.x(k.origin) << endl;
	k.initialize_line();
	k.initialize_interior();
	k.update_l();
	k.draw_lines();
	k.update_l();
	k.draw_lines();
	//k.update_l();
	//k.draw_lines();
	k.fill_corner();
//	k.plot_corners();
//	k.update_l();
//	k.draw_lines();
//	k.update_l();
//	k.draw_lines();
//	k.gplt.cmd("set size square");
	//k.fill_A_eff();
	//k.solve_A_eff();
	k.fill_A_eff();
	k.solve_A_eff();
	/*
	for(size_t it=0; it<10; it++){
		k.extract_eigvec_A_eff_biharmonic(it);
		k.plot_u_eff_biharmonic();
		sleep(3);
	}
	*/
	/*
	for(size_t it=0; it<50; it++){
		k.gplt.cmd("set term pdfcairo");
		k.gplt.cmd("set size square");
		stringstream ss1;
		ss1 << "set output 'data/drum_l2_highres_" << it << ".pdf'";
		string str;
		str = ss1.str();
		k.gplt.cmd(str);
		k.gplt.cmd("set xlabel 'x/L'");
		k.gplt.cmd("set ylabel 'y/L'");
		stringstream ss2;
		ss2 << "set title 'membrane, omega = " << k.omega_A_eff(50-1-it) << "'";
		string str2;
		str2 = ss2.str();
		k.gplt.cmd(str2);
		k.extract_eigvec_A_eff(it);
		k.plot_u_eff();
		k.gplt.cmd("unset output");
	}
	*/
	/*
	k.fill_A_eff_biharmonic();
	k.solve_A();
	k.solve_A_eff();
	k.solve_A_eff_biharmonic();

	for(size_t it=0; it<50; it++){
		k.extract_eigvec_A_eff_biharmonic(it);
		k.plot_u_eff_biharmonic();
		sleep(3);
	}
	*/

	//k.extract_eigvec_A_eff_biharmonic(2);
	//k.plot_u_eff_biharmonic();
	//k.extract_eigvec_A_eff(3);
	//k.extract_eigvec(9);
	//k.plot_u();
	//k.plot_u_eff();
	
/*	
	for(size_t it=0; it<69; it++){
		k.extract_eigvec_A_eff(it);
		k.plot_u_eff();
		sleep(1);
	}
*/

//	k.x.print("x: ");

	cout << "grid max: " << k.grid_max << endl;
	cout << "test grid max: " << ((double) (k.N-1)/2)*k.delta << endl;
//	cout << "grid_max*N/2: " << k.grid_max*((double) k.N/2.0) << endl;
	cout << "delta: " << k.delta << endl;
	cout << "N: " << k.N << endl;
	cout << "l: " << k.l << endl;
	cout << "matrix size: " << k.N*k.N << endl;
	cout << "number of interior points: " << accu(k.interior) << endl;

	return 0;
}
