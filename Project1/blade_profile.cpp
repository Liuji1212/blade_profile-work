#include <cmath>
#include "matplotlibcpp.h"
#define M_PI 3.1415926
#include <math.h>
#include <Eigen/Dense> 
#include <utility>
//#include "cubicspline.h"



using namespace Eigen;
using namespace Eigen::internal;
using namespace Eigen::Architecture;
using namespace std;
namespace plt = matplotlibcpp;


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






struct arc_para
{
    double arc_x;
    double arc_y;
    double arc_r;
    //double test;
};

struct poly_para
{
    double poly3;
    double poly2;
    double poly1;
    double poly0;
};

struct point
{
    double x;
    double y; 
};

poly_para poly(double v1, double v2, double vk, double w1, double w2, double wk);
arc_para arc(double x1, double y1, double x2, double y2, double angle);
vector<std::pair<double,double>> centerline(double endpoint_x, double intersection_x, 
                                                     double intersection_y, double lead_angle, double trail_angle);
//vector<std::pair<double, double>> suctionline(int index_maxcircle, double r_suc, vector<std::pair<double, double>> centerline_temp);
vector<double> the_suction_point(double mid_x, double mid_y, double left_x, double left_y, double right_x, double right_y, double r);
vector<double> the_pressure_point(double mid_x_p, double mid_y_p, double left_x_p, double left_y_p, double right_x_p, double right_y_p, double r_p);
//void polyfit(const vector<double>& x, const vector<double>& y, int n, vector<double>& a);
void PlotDraw_int(vector<vector<double>> plistInsert);




/*class profile {
private:
    pair<double, double> thick_point; 
    pair<double, double> thick_point_left;
    pair<double, double> thick_point_right;
    double r_circle;
     
public:
    profile(pair<double, double>in_thick_point, pair<double, double>in_thick_point_left,
        pair<double, double>in_thick_point_right, double in_r_circle) :
        thick_point(in_thick_point), thick_point_left(in_thick_point_left), thick_point_right(in_thick_point_right), r_circle(in_r_circle)
    {
        setparameter(in_thick_point,in_thick_point_left,in_thick_point_right,in_r_circle);
    }

    void setparameter(pair<double, double>set_thick_point, pair<double, double>set_thick_point_left,
        pair<double, double>set_thick_point_right, double set_r_circle)
    {
        thick_point = set_thick_point;
        thick_point_left = set_thick_point_left;
        thick_point_right = set_thick_point_right;
        r_circle = set_r_circle;
    }
    double circle_slope() const{
        double slope = (thick_point_left.second - thick_point_right.second) / (thick_point_left.first - thick_point_right.first);
        double slope_ver = -1.0 / slope;
        return slope_ver;
    } 

    point suction_point() const{
        double point_suc_y = thick_point.second + 
                          r_circle * abs(circle_slope() / sqrt(pow(circle_slope(), 2) + 1));
        double point_suc_x = (point_suc_y - thick_point.second) / circle_slope() + thick_point.first;

        point point_suc = { point_suc_x,point_suc_y };
        return point_suc;
    }

    point pressure_point() const {
        double point_pres_y = thick_point.second -
            r_circle * abs(circle_slope() / sqrt(pow(circle_slope(), 2) + 1));
        double point_pres_x = (point_pres_y - thick_point.second) / circle_slope() + thick_point.first;

        point point_pres = { point_pres_x, point_pres_y };
        return point_pres;
    }
}; 
*/

int main()
{
    
    double endpoint_x, intersection_x, intersection_y, lead_angle, trail_angle, max_circle_location;
    std::cout << "input endpoint_x: ";
    cin >> endpoint_x;
    std::cout << "input intersection_x: ";
    cin >> intersection_x;
    std::cout << "input intersection_y: ";
    cin >> intersection_y;
    std::cout << "input lead_angle: ";
    cin >> lead_angle;
    std::cout << "input trail_angle: ";
    cin >> trail_angle;
    std::cout << "input location of maximum circle: ";
    cin >> max_circle_location;
    
    double max_circle_r = endpoint_x * 1.5 / 78.0;
    double max_circle_x;
    double max_circle_y;
    int circle_quantity=8;
    vector<double> x_label;
    vector<double> y_label;
    vector<double> x_label_suc;
    vector<double> y_label_suc;
    vector<double> x_label_press;
    vector<double> y_label_press;



    auto centerline_imp = centerline( endpoint_x, intersection_x, intersection_y, lead_angle, trail_angle);
    for (int i = 0; i <= 10001; i++)
    {
        x_label.push_back(centerline_imp[i].first); 
        y_label.push_back(centerline_imp[i].second);
    }
    
    max_circle_x = max_circle_location;
    int index_maxcircle = max_circle_x / endpoint_x * 10000;
    //std::cout << index_maxcircle;
    //vector<vector<double>> profile_sucpoint;
    //vector<vector<double>> profile_presspoint;
    //profile_sucpoint.push_back(vector<double>());
    //profile_presspoint.push_back(vector<double>());
    //profile_sucpoint[0].push_back(0);
    //profile_sucpoint[0].push_back(0);
    //profile_presspoint[0].push_back(0);
    //profile_presspoint[0].push_back(0);

    double profile_sucpoint[19][2];
    double profile_presspoint[19][2];
    profile_sucpoint[0][0] = 0;
    profile_sucpoint[0][1] = 0;
    profile_presspoint[0][0] = 0;
    profile_presspoint[0][1] = 0;


    for (double m = 0.1; m <= 0.92; m = m + 0.05)
    {
        static int i = 1;
        int index_circle_x = m * 10000;
        double circle_r;

        if (index_circle_x < index_maxcircle)
            circle_r = max_circle_r * m * 2.2; 
        else
            circle_r = max_circle_r * (9.66 * pow(m, 4) - 20.46 * pow(m, 3) + 9.82 * pow(m, 2) + 0.975 * m);

        double midx = centerline_imp[index_circle_x].first;
        double midy = centerline_imp[index_circle_x].second;
        double leftx = centerline_imp[index_circle_x - 1].first;
        double lefty = centerline_imp[index_circle_x - 1].second;
        double rightx = centerline_imp[index_circle_x + 1].first;
        double righty = centerline_imp[index_circle_x + 1].second;

        //vector<double> loop_suc_point = the_suction_point(midx, midy, leftx, lefty, rightx, righty, circle_r);
        //profile_sucpoint[i].push_back(loop_suc_point[0]);
        //profile_sucpoint[i].push_back(loop_suc_point[1]);

        //vector<double> loop_press_point = the_pressure_point(midx, midy, leftx, lefty, rightx, righty, circle_r);
        //profile_presspoint[i].push_back(loop_press_point[0]);
        //profile_presspoint[i].push_back(loop_press_point[1]);
        

        vector<double> loop_suc_point = the_suction_point(midx, midy, leftx, lefty, rightx, righty, circle_r);
        vector<double> loop_press_point = the_pressure_point(midx, midy, leftx, lefty, rightx, righty, circle_r);
        profile_sucpoint[i][0] = loop_suc_point[0];
        profile_sucpoint[i][1] = loop_suc_point[1];
        profile_presspoint[i][0] = loop_press_point[0];
        profile_presspoint[i][1] = loop_press_point[1];

        std::cout << circle_r;
        i = i + 1;
        //std::cout << i;
        if (i == 18)
        {
            //profile_sucpoint[i].push_back(endpoint_x);
            //profile_sucpoint[i].push_back(0);
            //profile_presspoint[i].push_back(endpoint_x);
            //profile_presspoint[i].push_back(0);

            profile_sucpoint[i][0] = endpoint_x;
            profile_sucpoint[i][1] = 0;
            profile_presspoint[i][0] = endpoint_x;
            profile_presspoint[i][1] = 0;




        }
        //profile profile_p(pair<double, double> p_mid (midx,midy), pair<double, double>in_thick_point_left,
            //pair<double, double>in_thick_point_right, double in_r_circle)
    };

    vector<vector<double>> vec_profile_sucpoint;
    std::vector<double>vec_p_s;
    for (int i = 0; i < 19; i++) {
        for (int j = 0; j < 2; j++) {
         
            vec_p_s.push_back(profile_sucpoint[i][j]);
        }
        vec_profile_sucpoint.push_back(vec_p_s);
        vec_p_s.clear();
    }

    vector<vector<double>> vec_profile_presspoint;
    std::vector<double>vec_p_p;
    for (int i = 0; i < 19; i++) {
        for (int j = 0; j < 2; j++) {   

            vec_p_p.push_back(profile_presspoint[i][j]);
        }
        vec_profile_presspoint.push_back(vec_p_p);
        vec_p_p.clear();
    }

    for (int i = 0; i < 19; i++) 
    {
       x_label_suc.push_back(profile_sucpoint[i][0]);
       y_label_suc.push_back(profile_sucpoint[i][1]);
   
    }

    for (int i = 0; i < 19; i++)
    {
        x_label_press.push_back(profile_presspoint[i][0]);
        y_label_press.push_back(profile_presspoint[i][1]);

    }

    //std::cout << vec_profile_sucpoint[6][1] << "\n";



    CubicSpline cssuc;
    cssuc.CubicSpline_init(vec_profile_sucpoint);
    vector<vector<double>> pointInsert_suc = cssuc.CubicSpline_Insert(3);

    CubicSpline cspress;
    cspress.CubicSpline_init(vec_profile_presspoint);
    vector<vector<double>> pointInsert_press = cspress.CubicSpline_Insert(3);

    //std::cout << pointInsert_suc[20][0];
    //std::cout << pointInsert_suc[20][1];

    plt::xlim(0, 80); 
    plt::ylim(0, 30); 

    plt::plot(x_label, y_label);
    plt::plot(x_label_suc, y_label_suc);
    plt::plot(x_label_press, y_label_press);
    //PlotDraw_int(pointInsert_suc);
    //PlotDraw_int(pointInsert_press);
    plt::show();
    //std::cout << centerline_imp[10001].second;
    //std::cout << profile_sucpoint[2][0] << "\n";
    //std::cout << profile_sucpoint[5][1] << "\n";
    //std::cout << profile_presspoint[5][0] << "\n";
    //std::cout << profile_presspoint[7][1] << "\n";
}; 

arc_para arc(double x1, double y1, double x2, double y2, double angle) //angle need  minus 90
{
    arc_para arc_ins;
    double center_x, center_y;
    double k1_ver, b1_ver, k2, b2;

    k1_ver = -1.0 / ((y2 - y1) / (x2 - x1));
    center_x = (x1 + x2) / 2.0;
    center_y = (y1 + y2) / 2.0;
    b1_ver = center_y - k1_ver * center_x;

    k2 = tan((90 - angle) * M_PI / 180);
    b2 = y2 - k2 * x2;
    arc_ins.arc_x = (b2 - b1_ver) / (k1_ver - k2);
    arc_ins.arc_y = arc_ins.arc_x * k2 + b2;
    arc_ins.arc_r = sqrt(pow((arc_ins.arc_x - x1), 2) + pow((arc_ins.arc_y - y1), 2));
    //arc_ins.test = y1 +1;
    
    return arc_ins;
}

poly_para poly(double vx, double vy, double vk, double wx, double wy, double wk)
{
    poly_para poly_ins;
    Eigen::Matrix<double, 4, 4> matrix_left;
    Eigen::Vector4d vector_right;

    matrix_left << 3 * pow(vx, 2), 2 * vx, 1.0, 0.0,
                   3 * pow(wx, 2), 2 * wx, 1.0, 0.0,
                   pow(vx, 3), pow(vx, 2), vx, 1.0,
                   pow(wx, 3), pow(wx, 2), wx, 1.0;
    vector_right << vk, wk, vy, wy;

    Eigen::Vector4d solve_coeff = matrix_left.colPivHouseholderQr().solve(vector_right);
    
    poly_ins.poly3 = solve_coeff(0);
    poly_ins.poly2 = solve_coeff(1);
    poly_ins.poly1 = solve_coeff(2);
    poly_ins.poly0 = solve_coeff(3);

    return poly_ins;
}

 
vector<std::pair<double, double>> centerline(double endpoint_x, double intersection_x,
                                  double intersection_y, double lead_angle, double trail_angle)
{
    double startpoint_x = 0.0, startpoint_y = 0.0;
    double endpoint_y = 0.0;
                         //input the end point of airfoil

    double lead_k = tan(M_PI * lead_angle / 180);
    arc_para arc_part = arc(intersection_x, intersection_y, endpoint_x, endpoint_y, trail_angle);
    double intersection_k = -1.0 / ((intersection_y - arc_part.arc_y) / (intersection_x - arc_part.arc_x));
    poly_para poly_part = poly(startpoint_x, startpoint_y, lead_k, intersection_x, intersection_y, intersection_k);


    int n_poly = 5000, n_arc = 5000;
    int n_centerline = n_poly + n_arc+1;
    vector<double> centerline_x(n_centerline+1), centerline_y(n_centerline + 1);
    vector<std::pair<double, double>> centerline_ins;

    for (int i = 0; i <= n_centerline; i++)
    {
        centerline_x.at(i) = i * endpoint_x / n_centerline;

        if (i <= n_poly)
        {
            centerline_y.at(i) = poly_part.poly3 * pow(centerline_x.at(i), 3)
                + poly_part.poly2 * pow(centerline_x.at(i), 2)
                + poly_part.poly1 * centerline_x.at(i) + poly_part.poly0;
        }

        if (n_poly < i <= n_centerline)
        {
            centerline_y.at(i) = arc_part.arc_y 
                + sqrt(pow(arc_part.arc_r, 2) - pow((centerline_x.at(i) - arc_part.arc_x), 2));
        
        }

        centerline_ins.push_back(make_pair(centerline_x.at(i), centerline_y.at(i)));
        //centerline_ins.push_back(make_pair(arc_part.test, 999.0));
    }
    return centerline_ins;

}

//vector<std::pair<double, double>> suctionline(int index_maxcircle, double r_suc, vector<std::pair<double, double>> centerline_temp)
//{
    
    //std::pair<double, double> suction_mid<x, centerline_temp[index_maxcircle].second>;
    //profile suction(centerline_temp.at(index_maxcircle), )
//}

vector<double> the_suction_point(double mid_x, double mid_y, double left_x, double left_y, double right_x, double right_y, double r)

{
    vector<double> coor_suc;
    double tang = (right_y - left_y) / (right_x - left_x);
    double tang_ver = -1.0 / tang;
    coor_suc.reserve(200);
    coor_suc[1] = mid_y + r * abs(tang_ver / sqrt((pow(tang_ver, 2) + 1)));
    coor_suc[0] = (coor_suc[1] - mid_y) / tang_ver + mid_x;
    return coor_suc;
}

vector<double> the_pressure_point(double mid_x_p, double mid_y_p, double left_x_p, double left_y_p, double right_x_p, double right_y_p, double r_p)

{
    vector<double> coor_press;
    coor_press.reserve(200); 
    double tang = (right_y_p - left_y_p) / (right_x_p - left_x_p);
    double tang_ver = -1.0 / tang;
    coor_press[1] = mid_y_p - r_p * abs(tang_ver / sqrt((pow(tang_ver, 2) + 1)));
    coor_press[0] = (coor_press[1] - mid_y_p) / tang_ver + mid_x_p;
    return coor_press;
}


void PlotDraw_int(vector<vector<double>> plistInsert) {

  
    vector<double> xI;
    vector<double> yI;

 

    for (int i = 0; i < plistInsert.size(); i++) {
        xI.push_back(plistInsert[i][0]);
        yI.push_back(plistInsert[i][1]);
    }
    std::cout << xI[4];
    std::cout << yI[4];
    
    // Plot a red dashed line from given x and y data.
    plt::plot(xI, yI);

}