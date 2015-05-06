#ifndef GNUPLOTTING_H
#define GNUPLOTTING_H

#include <string>
#include <cstdio>
#include <vector>
#include <armadillo>
using namespace std;
using namespace arma;

/*
Gnuplotting & operator << (Gnuplotting & outstream, char * command){
	outstream.cmd(command);
	return outstream;
}
*/

class Gnuplotting {
	public:
		Gnuplotting();
		Gnuplotting(string filename);
		Gnuplotting(char *filename);
		~Gnuplotting();

		void xrange(double xmin, double xmax);
		void yrange(double ymin, double ymax);
		void title(string & titlename);
		void title(const char * titlename);
		void xlabel(string & x);
		void ylabel(string & y);
		void xlabel(const char * x);
		void ylabel(const char * y);
		void cmd(string & command);
		void cmd(const char *  command);

		Gnuplotting & operator<<(const char * command);

		void show_matrix(size_t &N, arma::mat & x);
		void xystream(vector<double> & x, vector<double> & y);
		void xystream(size_t & N, arma::vec & x, arma::vec & y);
		void xystream_replot(size_t & N, arma::vec & x, arma::vec & y, const char* title);
		void two_xystream(size_t &N1, vec &x1, vec &y1, const char* title1, size_t &N2, vec &x2, vec &y2, const char* title2);
		void two_xystream(vector<double> &x1, vector<double> &y1, const char* title1, vector<double> &x2, vector<double> &y2, const char* title2);
		void xyzstream(size_t & N, arma::vec & x, arma::vec & y, arma::mat & z);
		void heatmap(size_t & Nx, size_t & Ny, arma::mat & z);
		void show_matrix(size_t &N, umat & x);
		void heatmap_coords(vec &x, vec &y, mat &u, size_t &N, vector<double> &X, vector<double> &Y);
	private:
		string filename;
		FILE * pipe;
};

#endif

