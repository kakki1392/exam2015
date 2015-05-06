#include <armadillo>
#include <cstdio>
#include <vector>
#include "gnuplotting.h"

using namespace arma;
using namespace std;

class Point{
	public:
		Point();
		Point(size_t i_add, size_t j_add);
		~Point();
		bool equal(Point &p);
		size_t i;
		size_t j;
};


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
		double m;            //
		int m_int;
		double delta_min; //minimum grid step
		double delta;     //used grid step
		double grid_max;

		size_t lstp;
		size_t stp;
		size_t N;
		size_t origin;
		umat boundary;
		umat interior;
		umat corner;
		vec x;

		void fill_x();
		void fill_corner();

		void initialize_line();
		void initialize_interior();

		void test_lines();
		size_t intpower(size_t base, size_t exponent);
		void update_l();

		void draw_lines();
		void plot_single_line_corners(size_t it);
		void plot_corners();
		void plot_lines();
		void plot_boundary();
		void plot_interior();
		void plot_interior_boundary();
		void plot_u();

		void fill_interior();
		Gnuplotting gplt;

		void fill_A();
		void fill_A_eff();
		void fill_A_eff_biharmonic();

		void solve_A();
		void solve_A_eff();
		void solve_A_eff_biharmonic();

		void extract_eigvec(size_t k);
		size_t number_eig_val;
		vec eigval;
		vec omega;
		vec coeff_omega;
		mat eigvec;
		mat u;

		vec omega_A;
		vec omega_A_eff;
		vec omega_A_eff_biharmonic;

		vec coeff_A;
		vec coeff_A_eff;
		vec coeff_A_eff_biharmonic;
		
		mat eigvec_A;
		mat eigvec_A_eff;
		mat eigvec_A_eff_biharmonic;

		vec eigval_A;
		vec eigval_A_eff;
		vec eigval_A_eff_biharmonic;

		sp_mat A;
		sp_mat A_eff;
		sp_mat A_eff_biharmonic;

		umat B;
		umat B_eff;
		umat B_eff_biharmonic;

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

	



