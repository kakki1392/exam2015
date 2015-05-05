#include <armadillo>
#include <cstdio>
#include <vector>
#include "gnuplotting.h"

using namespace arma;
using namespace std;

class Line{
	public:
		Line();
		~Line();
		Line(size_t i_1, size_t j_1, size_t i_2, size_t j_2, int d);
		void print();
		size_t i1;
		size_t j1;
		size_t i2;
		size_t j2;
		int dir;
};

class Koch{
	public:
		Koch();
		~Koch();
		vector<Line> lines;
		size_t n;         //number of lines
		int l_max;        //maximum allowed l
		int l;            //current l
		double s;         //length of line
		size_t s_index;   //length of line in indeces
		double L;         //length of line at l=0
		int m;            //
		double delta_min; //minimum grid step
		double delta;     //used grid step
		double grid_max;

		size_t lstp;
		size_t stp;
		size_t N;
		size_t origin;
		umat boundary;
		umat interior;
		mat grid;
		vec x;

		void fill_x();

		void initialize_line();
		void initialize_interior();

		void test_lines();
		size_t intpower(size_t base, size_t exponent);
		void update_l();
		void draw_lines();
		void plot_single_line_corners(size_t it);
		void plot_lines();
		void plot_boundary();
		void plot_interior();
		void plot_interior_boundary();
		void fill_interior();
		Gnuplotting gplt;

		sp_mat A;
		void fill_A();
		void solve_A();

		umat B;
		umat C;
		umat D;

};
/*
class Eigensolv{
	public:
		Eigensolv();
		~Eigensolv();

		Koch curve;

		void initialize();
		void set_l();
		void fill_mat();
		spmat A;

};

*/

	



