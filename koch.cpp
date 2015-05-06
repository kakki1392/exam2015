#include "koch.h"
#include <cmath>
#include <iostream>

Line::Line(){
	i1 = 0;
	i2 = 0;
	j1 = 0;
	j2 = 0;
	dir = 0;
}

Line::Line(size_t i_1, size_t j_1, size_t i_2, size_t j_2, int d){
	i1 = i_1;
	j1 = j_1;
	i2 = i_2;
	j2 = j_2;
	dir = d;
}

Line::~Line(){}

void Line::print(){
	cout << "i1 :" << i1 << endl
	     << "j1 :" << j1 << endl
	     << "i2 :" << i2 << endl
	     << "j2 :" << j2 << endl
	     << "dir: " << dir << endl << endl;
}


Koch::Koch(){
	number_eig_val = 10;
	n = 4;
	l_max = 2;
	l = 0;
	L = 1.0;
	s = 1.0;
	m = 2.0;
	m_int = 2;
	delta_min = L*pow(0.25,l_max);
	delta = delta_min/m;
	grid_max = L/2.0;
	for(int i = 1; i<(l_max+1);i++){
		grid_max = grid_max + L*pow(0.25,i);
	}
	lines = vector<Line>(n);	
	double N_double = grid_max/delta;
	N = (size_t) N_double;
	N = 2*N + 1; //+3
	origin = N/2; 
	//boundary = zeros<umat>(N,N);
	//interior = zeros<umat>(N,N);
	//grid = zeros<mat>(N,N);
}

Koch::~Koch(){}

void Koch::fill_x(){
	x = zeros<vec>(N);
	for(size_t i=(origin+1); i<N; i++){
		x(i) = x(i-1) + delta;
	}
	for(size_t i=0; i<(origin); i++){
		x(i) = -x(N-1-i);
	}
}

size_t Koch::intpower(size_t base, size_t exponent){
	size_t out = base;
	if(exponent==0){
		return base;
	}
	for(size_t i=0; i<(exponent-1); i++){
		out = out*base;
	}
	return out;
}

void Koch::plot_single_line_corners(size_t it){
	vec u = zeros<vec>(2);
	vec v = zeros<vec>(2);
	size_t length = 2;
	u(0) = x(lines[it].i1);
	u(1) = x(lines[it].i2);
	v(0) = x(lines[it].j1);
	v(1) = x(lines[it].j2);
	gplt.xystream(length,u,v);
}

void Koch::plot_lines(){
	vector<double> u ; 
	vector<double> v;
	for(size_t it=0; it<n; it++){
		u.push_back(x(lines[it].i1));
		u.push_back(x(lines[it].i2));
		v.push_back(x(lines[it].j1));
		v.push_back(x(lines[it].j2));
	}
	gplt.xystream(u,v);
}

void Koch::plot_interior(){
	vector<double> u; 
	vector<double> v;
	for(size_t i=0; i<N; i++){
		for(size_t j=0; j<N; j++){
			if(interior(i,j)==1){
				u.push_back(x(i));
				v.push_back(x(j));
			}
		}
	}
	gplt.xystream(u,v);
}


void Koch::plot_boundary(){
	vector<double> u; 
	vector<double> v;
	for(size_t i=0; i<N; i++){
		for(size_t j=0; j<N; j++){
			if(boundary(i,j)==0){
				u.push_back(x(i));
				v.push_back(x(j));
			}
		}
	}
	gplt.xystream(u,v);
}

void Koch::plot_interior_boundary(){
	vector<double> u; 
	vector<double> v;
	vector<double> u_out; 
	vector<double> v_out;
	for(size_t i=0; i<N; i++){
		for(size_t j=0; j<N; j++){
			if(boundary(i,j)==0){
				u_out.push_back(x(i));
				v_out.push_back(x(j));
			}
			if(interior(i,j)==1){
				u.push_back(x(i));
				v.push_back(x(j));
			}
		}
	}
	gplt.two_xystream(u_out,v_out,"boundary",u,v,"interior");
}

void Koch::plot_u(){
	vector<double> X ; 
	vector<double> Y;
	for(size_t it=0; it<n; it++){
		X.push_back(x(lines[it].i1));
		X.push_back(x(lines[it].i2));
		Y.push_back(x(lines[it].j1));
		Y.push_back(x(lines[it].j2));
	}
	gplt.heatmap_coords(x,x,u,N,X,Y);
}

void Koch::initialize_line(){
	size_t steps = m_int*intpower(4,l_max)/4;
	stp = steps;
	s_index = 4*steps + 1;
	cout << "steps: " << steps << endl;
	cout << "s_index: " << s_index << endl;
	size_t longstep = 2*steps;
	/*
	stp = steps;
	lstp = longstep;
	cout << "longstep: " << longstep << endl;
	interior.submat(origin-longstep, origin-longstep, origin+longstep, origin+longstep) = ones<umat>(2*longstep+1,2*longstep+1);
	boundary = interior;
	boundary.submat(origin-longstep+1, origin-longstep+1, origin+longstep-1, origin+longstep-1) = zeros<umat>(2*(longstep-1)+1,2*(longstep-1)+1);
	interior.save("interior.dat", arma_ascii);
	boundary.save("boundary.dat", arma_ascii);
	*/
	
	//line 1
	Line line1;
	line1.i1 = origin - longstep;
	line1.i2 = origin + longstep;
	line1.j1 = origin + longstep;
	line1.j2 = line1.j1;
	line1.dir = 1;
	lines[0] = line1;
	
	//line 2
	Line line2;
	line2.i1 = origin + longstep;
	line2.i2 = line2.i1;
	line2.j1 = origin + longstep;
	line2.j2 = origin - longstep;
	line2.dir = 2;
	lines[1] = line2;
	
	//line 3
	Line line3;
	line3.i1 = origin + longstep;
	line3.i2 = origin - longstep;
	line3.j1 = origin - longstep; line3.j2 = line3.j1; line3.dir = 3;
	lines[2] = line3;
	
	//line 4
	Line line4;
	line4.i1 = origin - longstep;
	line4.i2 = line4.i1;
	line4.j1 = origin - longstep;
	line4.j2 = origin + longstep;
	line4.dir = 4;
	lines[3] = line4;
}

void Koch::initialize_interior(){
	interior = zeros<umat>(N,N);
	umat box = ones<umat>(s_index-2,s_index-2);
	size_t low = origin-2*stp+1;
	size_t high = origin+2*stp-1;
	//interior.submat(low,low,high,high) = box;
	interior.submat(lines[3].i1+1,lines[3].j1+1,lines[0].i2-1,lines[0].j2-1) = box;
	interior.save("data/interior_init.dat",arma_ascii);
}

void Koch::draw_lines(){
	cout << "entering draw_lines" << endl;
	boundary = ones<umat>(N,N);
	uvec ins_col = zeros<uvec>(s_index);
	urowvec ins_row = zeros<urowvec>(s_index);
	cout << "s_index: " << s_index << endl;
	for(size_t it=0; it<lines.size(); it++){
		if(lines[it].dir == 1){
			boundary.col(lines[it].j1).subvec(lines[it].i1, lines[it].i2) = ins_col;
			interior.col(lines[it].j1).subvec(lines[it].i1, lines[it].i2) = ins_col;
		}
		else if(lines[it].dir == 2){
			boundary.row(lines[it].i1).subvec(lines[it].j2, lines[it].j1) = ins_row;
			interior.row(lines[it].i1).subvec(lines[it].j2, lines[it].j1) = ins_row;
		}
		else if(lines[it].dir == 3){
			boundary.col(lines[it].j1).subvec(lines[it].i2, lines[it].i1) = ins_col;
			interior.col(lines[it].j1).subvec(lines[it].i2, lines[it].i1) = ins_col;
		}
		else{
			boundary.row(lines[it].i1).subvec(lines[it].j1, lines[it].j2) = ins_row;
			interior.row(lines[it].i1).subvec(lines[it].j1, lines[it].j2) = ins_row;
		}
	}
	//boundary.save("data/boundaryBigBig5.dat", raw_ascii);
	interior.save("data/interior_draw.dat",arma_ascii);
	cout << "leaving draw_lines" << endl;
}

/*
void Koch::fill_interior(){
	for(size_t j=0; j<N; j++){
		size_t start = 0;
		size_t stop = 0;
		size_t i = 0;
		bool hitwall = false;
		bool outsideBoundary = true;
		for(size_t i=0; i<N, i++){
			if(outsideBoundary){
				boundary(i,j) = 0;
			}
			size_t peek = i+1;
			if(boundary(peek,j)==0){
				hitwall = true;

		while(!hitwall){
			boundary(i,j) = 0;
			i = i+1;
			if(boundary(i,j)==0){
				hitwall = true;
				start = i;
*/




void Koch::test_lines(){
	umat test = zeros<umat>(N,N);

	uvec del_col = zeros<uvec>(stp-1);
	uvec ins_col = ones<uvec>(2*lstp+1);
	urowvec del_row = zeros<urowvec>(stp-1);
	urowvec ins_row = ones<urowvec>(2*lstp+1);

	ins_col.print();
	cout << endl;
	ins_row.print();
	cout << endl;

	for(size_t it=0; it<lines.size(); it++){
		if(lines[it].dir == 1){
			test.col(lines[it].j1).subvec(lines[it].i1, lines[it].i2) = ins_col;
		}
		else if(lines[it].dir == 2){
			test.row(lines[it].i1).subvec(lines[it].j1, lines[it].j2) = ins_row;
		}
		else if(lines[it].dir == 3){
			test.col(lines[it].j1).subvec(lines[it].i2, lines[it].i1) = ins_col;
		}
		else{
			test.row(lines[it].i1).subvec(lines[it].j2, lines[it].j1) = ins_row;
		}
	}
	test.save("test_lines.dat",arma_ascii);
}

void Koch::update_l(){
	if((l+1) > l_max){
		return;
	}
	cout << "entering update_l()" << endl;
	int new_l = l+1;
	size_t new_n = 8*n;
	size_t s = m_int*intpower(4,l_max)/intpower(4,new_l); //CAREFUL
	size_t longstep = 2*s;
	cout << "new_l: " << new_l << endl;
	cout << "steps: " << s << endl;
	cout << "longstep: " << longstep << endl;
	vector<Line> new_lines;
	s_index = s + 1;
	l = new_l;
	n = new_n;
	cout << "s_index:" << s_index << endl;

	//update interior first
	umat add = ones<umat>(s_index,s_index);
	umat del = zeros<umat>(s_index,s_index);
	for(size_t it=0; it<lines.size(); it++){
		if(lines[it].dir == 1){
			//cout << "first line" << endl;
			interior.submat(lines[it].i1+s,lines[it].j1,lines[it].i1+2*s,lines[it].j1+s) = add;
			interior.submat(lines[it].i1+2*s,lines[it].j1-s,lines[it].i1+3*s,lines[it].j1) = del;
		}else if(lines[it].dir == 2){
			//cout << "second line" << endl;
			interior.submat(lines[it].i1,lines[it].j1-2*s,lines[it].i1+s,lines[it].j1-s) = add;
			interior.submat(lines[it].i1-s,lines[it].j1-3*s,lines[it].i1,lines[it].j1-2*s) = del;
		}else if(lines[it].dir == 3){
			//cout << "third line" << endl;
			interior.submat(lines[it].i1-2*s,lines[it].j1-s,lines[it].i1-s,lines[it].j1) = add;
			interior.submat(lines[it].i1-3*s,lines[it].j1,lines[it].i1-2*s,lines[it].j1+s) = del;
		}else{
			//cout << "fourth line" << endl;
			interior.submat(lines[it].i1-s,lines[it].j1+s,lines[it].i1,lines[it].j1+2*s) = add;
			interior.submat(lines[it].i1,lines[it].j1+2*s,lines[it].i1+s,lines[it].j1+3*s) = del;
		}
	}
	interior.save("data/interior_update.dat",arma_ascii);

	//Update lines
	for(size_t it=0; it<lines.size(); it++){
		if(lines[it].dir == 1){ 
			size_t i1 = lines[it].i1;
			size_t j1 = lines[it].j1;
			Line line1(i1,j1,i1+s,j1,1);
			new_lines.push_back(line1);
			Line line2(i1+s,j1,i1+s,j1+s,4);
			new_lines.push_back(line2);
			Line line3(i1+s,j1+s,i1+2*s,j1+s,1);
			new_lines.push_back(line3);
			Line line4(i1+2*s,j1+s,i1+2*s,j1,2);
			new_lines.push_back(line4);
			Line line5(i1+2*s,j1,i1+2*s,j1-s,2);
			new_lines.push_back(line5);
			Line line6(i1+2*s,j1-s,i1+3*s,j1-s,1);
			new_lines.push_back(line6);
			Line line7(i1+3*s,j1-s,i1+3*s,j1,4);
			new_lines.push_back(line7);
			Line line8(i1+3*s,j1,i1+4*s,j1,1);
			new_lines.push_back(line8);
		}
		else if(lines[it].dir == 2){
			size_t i1 = lines[it].i1;
			size_t j1 = lines[it].j1;
			Line line1(i1,j1,i1,j1-s,2);
			new_lines.push_back(line1);
			Line line2(i1,j1-s,i1+s,j1-s,1);
			new_lines.push_back(line2);
			Line line3(i1+s,j1-s,i1+s,j1-2*s,2);
			new_lines.push_back(line3);
			Line line4(i1+s,j1-2*s,i1,j1-2*s,3);
			new_lines.push_back(line4);
			Line line5(i1,j1-2*s,i1-s,j1-2*s,3);
			new_lines.push_back(line5);
			Line line6(i1-s,j1-2*s,i1-s,j1-3*s,2);
			new_lines.push_back(line6);
			Line line7(i1-s,j1-3*s,i1,j1-3*s,1);
			new_lines.push_back(line7);
			Line line8(i1,j1-3*s,i1,j1-4*s,2);
			new_lines.push_back(line8);
		}
		else if(lines[it].dir == 3){
			size_t i1 = lines[it].i1;
			size_t j1 = lines[it].j1;
			Line line1(i1,j1,i1-s,j1,3);
			new_lines.push_back(line1);
			Line line2(i1-s,j1,i1-s,j1-s,2);
			new_lines.push_back(line2);
			Line line3(i1-s,j1-s,i1-2*s,j1-s,3);
			new_lines.push_back(line3);
			Line line4(i1-2*s,j1-s,i1-2*s,j1,4);
			new_lines.push_back(line4);
			Line line5(i1-2*s,j1,i1-2*s,j1+s,4);
			new_lines.push_back(line5);
			Line line6(i1-2*s,j1+s,i1-3*s,j1+s,3);
			new_lines.push_back(line6);
			Line line7(i1-3*s,j1+s,i1-3*s,j1,2);
			new_lines.push_back(line7);
			Line line8(i1-3*s,j1,i1-4*s,j1,3);
			new_lines.push_back(line8);
		}
		else{
			size_t i1 = lines[it].i1;
			size_t j1 = lines[it].j1;
			Line line1(i1,j1,i1,j1+s,4);
			new_lines.push_back(line1);
			Line line2(i1,j1+s,i1-s,j1+s,3);
			new_lines.push_back(line2);
			Line line3(i1-s,j1+s,i1-s,j1+2*s,4);
			new_lines.push_back(line3);
			Line line4(i1-s,j1+2*s,i1,j1+2*s,1);
			new_lines.push_back(line4);
			Line line5(i1,j1+2*s,i1+s,j1+2*s,1);
			new_lines.push_back(line5);
			Line line6(i1+s,j1+2*s,i1+s,j1+3*s,4);
			new_lines.push_back(line6);
			Line line7(i1+s,j1+3*s,i1,j1+3*s,3);
			new_lines.push_back(line7);
			Line line8(i1,j1+3*s,i1,j1+4*s,4);
			new_lines.push_back(line8);
		}
	}
	lines = new_lines;
	cout << "leaving update_l()" << endl;
}

void Koch::fill_A(){
	cout << "entering fill_A()" << endl;
	A = sp_mat(N*N,N*N);
	B = zeros<umat>(N*N,N*N);
	size_t i = 0;
	for(size_t row=0; row<(N*N); row++){
		if(i==N){
			i=0;
		}
		size_t max = N*N - 1;
		size_t i = row % N;
		size_t j = row/N;
		if(interior(i,j)==0){
			continue;
		}
		//A(row,row) = 0.0;
		B(row,row) = 0;
		//left
		if((i-1) >= 0){
			if(interior(i-1,j)==1){
				A(row,row-1) = 1.0;
				B(row,row-1) = 1;
			}
		}
		//right
		if((i+1) <= max){
			if(interior(i+1,j)==1){
				A(row,row+1) = 1.0;
				B(row,row+1) = 1;
			}
		}
		//up
		if((j+1) <= max){
			if(interior(i,j+1)==1){
				A(row,row+N) = 1.0;
				B(row,row+N) = 1;
			}
		}
		//down
		if((j-1) >= 0){
			if(interior(i,j-1)==1){
				A(row,row-N) = 1.0;
				B(row,row-N) = 1;
			}
		}
		i++;
	}
	C = trans(B);
	D = (B==C);
	cout << "leaving fill_A()" << endl;
	//B.save("data/matrix.dat", arma_ascii);
	//size_t length = N*N;
	//gplt.show_matrix(length,B);
}

void Koch::solve_A(){
	eigs_sym(eigval, eigvec, A, number_eig_val,"la", 1.0e-6);
	//eigval = eigs_sym(A,number_eig_val,"la",1.0e-10);
	omega = eigval;
//	vec coeff_omega = (1.0/(M_PI*sqrt(2)))*omega;
	for(size_t i=0; i<number_eig_val; i++){
		omega(i) = sqrt((4.0 - eigval(i))/(delta*delta));
	}
	coeff_omega = (1.0/(M_PI*sqrt(2)))*omega;
	eigval.raw_print(cout,"eigvals: ");
	omega.raw_print(cout,"omega: ");
	coeff_omega.raw_print(cout,"coeffs: ");
	eigval.save("data/eigval_l2.dat", raw_ascii);
	omega.save("data/omega_l2.dat",raw_ascii);
}

void Koch::extract_eigvec(size_t k){
	vec u_vec = eigvec.col(k);
	u_vec.save("data/eigvec.dat",raw_ascii);
	u = zeros<mat>(N,N);
	size_t start = 0;
	size_t stop = 0;
	for(size_t it=0; it<N; it++){
		stop = start + N - 1;
		u.col(it) = u_vec.subvec(start,stop);
		start = stop + 1;
	}
}


/*
Eigensolv::Eigensolv(){
	curve = Koch();
}

Eigensolv::~Eigensolv(){}

void Eigensolv::initialize(){
	curve.initialize_line();
	curve.initialize_interior();
}

*/




