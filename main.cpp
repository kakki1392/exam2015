#include "koch.h"
#include <iostream>
#include <vector>
#include <unistd.h>
using namespace arma;
using namespace std;

int main(){
	cout.precision(15);

	Koch k;
	/*
	k.gplt.cmd("set term pdfcairo");
	k.gplt.cmd("set output 'data/koch3.pdf'");
	k.gplt.cmd("set size square");
	k.gplt.cmd("set xlabel 'x'");
	k.gplt.cmd("set ylabel 'y'");
	k.gplt.cmd("set title 'l = 3, L = 1'");
	*/


	k.fill_x();
	//k.x.print("x: ");
	//cout << k.x(k.origin) << endl;
	k.initialize_line();
	k.initialize_interior();
	k.update_l();
	k.draw_lines();
	k.update_l();
	k.draw_lines();
//	k.update_l();
//	k.draw_lines();
	//k.update_l();
//	k.draw_lines();
	k.gplt.cmd("set size square");
	k.plot_interior_boundary();
	//k.plot_lines();
	k.fill_A();
	k.solve_A();
	cout << "grid max: " << k.grid_max << endl;
	cout << "test grid max: " << ((double) (k.N-1)/2)*k.delta << endl;
//	cout << "grid_max*N/2: " << k.grid_max*((double) k.N/2.0) << endl;
	cout << "delta: " << k.delta << endl;
	cout << "N: " << k.N << endl;
	cout << "l: " << k.l << endl;
	cout << "matrix size: " << k.N*k.N << endl;
	cout << "number of interior points: " << accu(k.interior) << endl;
	cout << "points in B: " << accu(k.B) << endl;
	cout << "points in C: " << accu(k.C) << endl;
	cout << "points in D: " << accu(k.D) << endl;
	//k.plot_lines();
	//k.draw_lines();
	//k.plot_boundary();
	//k.plot_lines();
	//k.draw_lines();
	//sleep(4);
	//k.plot_boundary();
	//k.plot_lines();
	//k.draw_lines();
	//sleep(4);
	//k.plot_boundary();
	/*
	k.update_l();
	sleep(2);
	k.plot_lines();
	k.update_l();
	sleep(2);
	k.plot_lines();
	*/
	//k.draw_lines();
	//sleep(4);
	//k.plot_lines();
	//k.plot_single_line_corners(3);
	/*
	k.plot_lines();
//	k.draw_lines();
	k.update_l();
	sleep(4);
	k.plot_lines();
	k.update_l();
	sleep(4);
	k.plot_lines();
	sleep(4);
	k.update_l();
	k.plot_lines();
	cout << k.lines.size() << endl;
	for(size_t i=0; i<k.lines.size(); i++){
		k.lines[i].print();
	}
	*/

	return 0;
}
