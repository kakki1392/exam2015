#include <armadillo>
#include <cstdio>
#include <vector>
#include "gnuplotting.h"

using namespace arma;
using namespace std;

/* POINT CLASS
 * This class represents a point on the grid.
 * The point is defined by two indices (i,j).
 * One can check if a point is equal to another point.
 */
class Point{
	public:
		Point();
		Point(size_t i_add, size_t j_add);
		~Point();
		bool equal(Point &p);
		size_t i;
		size_t j;
};


/* LINE CLASS
 * This class represents a line on the Koch curve.
 * The line is defined by start and end indices
 * (i1,j1) ----- (i2,j2)
 * It also has a normal vector "dir".
 * 	dir = 1 is \hat y
 * 	dir = 2 is \hat x
 * 	dir = 3 is -\hat y
 * 	dir = 4 is -\hat x
 */
class Line{
	public:
		Line();
		~Line();
		Line(size_t i_1, size_t j_1, size_t i_2, size_t j_2, int d);
		size_t i1;
		size_t j1;
		size_t i2;
		size_t j2;
		int dir;
};

/* KOCH CLASS
 * This class represents the physical system.
 * It contains:
 * 	A set of <class> Lines that defines the Koch curve
 * 	Variables that define the size of the grid
 * 	Matrices that define if point on grid is interior/corner/boundary
 * 	Functions that update grid for each fractal level "l"
 * 	Matrices in eigenvalue problem that can be solved
 *
 * 	
 *
 */
class Koch{
	public:
		Koch();
		~Koch();
		vector<Line> lines;
		size_t n;           //Number of lines in Koch curve
		int l_max;          //Maximum allowed l
		int l;              //Current l
		size_t s_index;     //Length of line in indeces
		double L;           //Length of line at l=0
		double m;           //Defines resolution of grid
		int m_int;          //Defines resolution of grid
		double delta_min;   //Minimum grid step
		double delta;       //Used grid step, delta=delta_min/m
		double grid_max;    //Max value of x (y) on grid

		size_t stp;
		size_t N;           //Size of grid
		size_t origin;      //Index of origin
		
		//Matrices that defines boundary, interior and corners.
		umat boundary;
		umat interior;
		umat corner;

		//vectors that defines (x,y) on grid
		vector<Point> points;
		vec x;

		void fill_x();
		void fill_corner();

		void initialize_line();
		void initialize_interior();

		void update_l();
		void draw_lines();

		size_t intpower(size_t base, size_t exponent);

		//Plotting functions. No external data storage
		Gnuplotting gplt;
		void plot_single_line_corners(size_t it);
		void plot_corners();
		void plot_lines();
		void plot_boundary();
		void plot_interior();
		void plot_interior_boundary();
		void plot_u();
		void plot_u_eff();
		void plot_u_eff_biharmonic();

		//Initialize eigenvalue matrices
		void fill_A();
		void fill_A_eff();
		void fill_A_eff_biharmonic();

		//Solve eigenvalue matrices
		size_t number_eig_val;    //Number of eigenvalues wanted
		void solve_A();
		void solve_A_eff();
		void solve_A_eff_biharmonic();

		//Exctract eigenvectors
		void extract_eigvec(size_t k);
		void extract_eigvec_A_eff(size_t k);
		void extract_eigvec_A_eff_biharmonic(size_t k);

		vec eigval;
		vec omega;
		vec coeff_omega;
		mat eigvec;
		mat u;

		//The matrices that will be solved
		sp_mat A;
		sp_mat A_eff;
		sp_mat A_eff_biharmonic;

		//Clones of A matrices, only used for visualizing
		//big matrices
		umat B;
		umat B_eff;
		umat B_eff_biharmonic;
		
		//Grid defining u(x,y) for k'th eigenvalue
		mat u_A;
		mat u_A_eff;
		mat u_A_eff_biharmonic;
		
		//Storage of all eigenvalues
		vec eigval_A;
		vec eigval_A_eff;
		vec eigval_A_eff_biharmonic;

		//Storage of all frequencies omega
		vec omega_A;
		vec omega_A_eff;
		vec omega_A_eff_biharmonic;

		//Storage of all frequency coefficients
		vec coeff_A;
		vec coeff_A_eff;
		vec coeff_A_eff_biharmonic;
		
		//Storage of all eigenvectors
		mat eigvec_A;
		mat eigvec_A_eff;
		mat eigvec_A_eff_biharmonic;

};
	
