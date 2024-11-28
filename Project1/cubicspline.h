#pragma once

#include <iostream>
#include <Eigen/Dense>

//using Eigen::MatrixXd;
using namespace Eigen;
using namespace Eigen::internal;
using namespace Eigen::Architecture;

using namespace std;

class CubicSpline {

private:

	vector<vector<double>> point2D;
	vector<vector<double>> CubicSplineParameter;//a, b, c, d.
	vector<double> h;
	vector<double> m;

public:

	void CubicSpline_init(vector<vector<double>> point2D_input) {

		point2D = point2D_input;

		//init h
		h.clear();
		h.resize(point2D.size() - 1);
		for (int i = 0; i < point2D.size() - 1; i++) {
			double x1 = point2D[i][0];
			double x2 = point2D[i + 1][0];
			double h_i = abs(x2 - x1);
			h[i] = h_i;
		}

		//init m. m.size = point2D.size()
		//1, compute yh coefficient
		vector<double> yh(point2D_input.size());
		for (int i = 0; i < yh.size(); i++) {
			if (i == 0 || i == yh.size() - 1) {
				yh[i] = 0;
			}
			else {
				yh[i] = 6 * ((point2D[i + 1][1] - point2D[i][1]) / h[i] - (point2D[i][1] - point2D[i - 1][1]) / h[i - 1]);
			}
		}

		MatrixXf A(point2D.size(), point2D.size());
		MatrixXf B(point2D.size(), 1);
		MatrixXf m;

		//2, init A, B
		B(0, 0) = yh[0];
		B(point2D.size() - 1, 0) = yh[point2D.size() - 1];

		for (int i = 0; i < point2D.size() - 1; i++) {

			A(0, i) = 0;
			A(point2D.size() - 1, i) = 0;

		}
		A(0, 0) = 1;
		A(point2D.size() - 1, point2D.size() - 2) = 1;

		for (int i = 1; i < point2D.size() - 1; i++) {

			B(i, 0) = yh[i];

			for (int j = 0; j < point2D.size(); j++) {

				if (j == i) {
					A(i, j) = 2 * (h[i - 1] + h[i]);
				}
				else if (j == i - 1) {
					A(i, j) = h[i - 1];
				}
				else if (j == i + 1) {
					A(i, j) = h[i];
				}
				else {
					A(i, j) = 0;
				}

			}

		}

		m = A.llt().solve(B);
		vector<double> mV(point2D.size());
		for (int i = 0; i < point2D.size(); i++) {
			mV[i] = m(i, 0);
		}

		for (int i = 0; i < m.size() - 1; i++) {

			vector<double> CubicSplineParameter_i;
			double a = point2D[i][1];
			double b = (point2D[i + 1][1] - point2D[i][1]) / h[i] - h[i] / 2 * mV[i] - h[i] / 6 * (mV[i + 1] - mV[i]);
			double c = mV[i] / 2;
			double d = (mV[i + 1] - mV[i]) / (6 * h[i]);
			CubicSplineParameter_i.push_back(a);
			CubicSplineParameter_i.push_back(b);
			CubicSplineParameter_i.push_back(c);
			CubicSplineParameter_i.push_back(d);
			CubicSplineParameter.push_back(CubicSplineParameter_i);

		}

	}

	vector<vector<double>> CubicSpline_Insert(int step) {

		vector<vector<double>> insertList;

		for (int i = 0; i < CubicSplineParameter.size(); i++) {
			double h_i = h[i] / (double)step;
			insertList.push_back(point2D[i]);
			double a = CubicSplineParameter[i][0];
			double b = CubicSplineParameter[i][1];
			double c = CubicSplineParameter[i][2];
			double d = CubicSplineParameter[i][3];
			for (int j = 1; j < step; j++) {
				double x_new = point2D[i][0] + h_i * j;
				double y_new = a + b * (x_new - point2D[i][0])
					+ c * (x_new - point2D[i][0]) * (x_new - point2D[i][0])
					+ d * (x_new - point2D[i][0]) * (x_new - point2D[i][0]) * (x_new - point2D[i][0]);
				vector<double> p_new_ij;
				p_new_ij.push_back(x_new);
				p_new_ij.push_back(y_new);
				insertList.push_back(p_new_ij);
			}
		}

		insertList.push_back(point2D[point2D.size() - 1]);
		return insertList;

	}

	vector<vector<double>> CubicSpline_Insert(double step) {

		vector<vector<double>> insertList;

		for (int i = 0; i < CubicSplineParameter.size(); i++) {
			int h_i = h[i] / (double)step;
			insertList.push_back(point2D[i]);
			double a = CubicSplineParameter[i][0];
			double b = CubicSplineParameter[i][1];
			double c = CubicSplineParameter[i][2];
			double d = CubicSplineParameter[i][3];
			for (int j = 1; j < h_i; j++) {
				double x_new = point2D[i][0] + step * j;
				double y_new = a + b * (x_new - point2D[i][0])
					+ c * (x_new - point2D[i][0]) * (x_new - point2D[i][0])
					+ d * (x_new - point2D[i][0]) * (x_new - point2D[i][0]) * (x_new - point2D[i][0]);
				vector<double> p_new_ij;
				p_new_ij.push_back(x_new);
				p_new_ij.push_back(y_new);
				insertList.push_back(p_new_ij);
			}
		}

		insertList.push_back(point2D[point2D.size() - 1]);
		return insertList;

	}

};
 