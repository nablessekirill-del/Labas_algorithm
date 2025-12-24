#ifndef ANALYSIS_H
#define ANALYSIS_H

#include <vector>
#include <string>
#include <cmath>
#include <algorithm>


struct PlotData {
    std::vector<double> x_emp;
    std::vector<double> y_emp;
    std::vector<double> x_low;
    std::vector<double> x_est;
    std::vector<double> x_up;
    std::vector<double> y_line;
    double param1; // c_hat
    double param2; // b_hat
    std::vector<std::vector<double>> cov;

    double getXMin() const {
        double xmin = 1e18;
        if(!x_emp.empty()) xmin = std::min(xmin, *std::min_element(x_emp.begin(), x_emp.end()));
        if(!x_low.empty()) xmin = std::min(xmin, *std::min_element(x_low.begin(), x_low.end()));
        return (xmin == 1e18) ? 0 : xmin * 0.8;
    }
    double getXMax() const {
        double xmax = -1e18;
        if(!x_emp.empty()) xmax = std::max(xmax, *std::max_element(x_emp.begin(), x_emp.end()));
        if(!x_up.empty()) xmax = std::max(xmax, *std::max_element(x_up.begin(), x_up.end()));
        return (xmax == -1e18) ? 1000 : xmax * 1.1;
    }
};

void calculate_weibull_intervals(PlotData& pd, const std::vector<double>& p_vec, double beta);

struct OptData {
    int n;
    std::vector<double> x;
    std::vector<int> r;
};
extern OptData globalOptData;

struct Sample {
    std::vector<double> x;
    std::vector<int> r;
    int n = 0;
};

struct EmpiricalKM {
    std::vector<double> x_sorted;
    std::vector<double> F_emp;
};

extern const std::vector<double>* G_X;
extern const std::vector<int>* G_R;


double norm_pdf(double z);
double norm_cdf(double z);
double norm_ppf(double p);
EmpiricalKM kaplan_meier_Itype(const std::vector<double>& x, const std::vector<int>& r);


std::pair<double, double> weibull_mle_2par(const std::vector<double>& x, const std::vector<int>& r);
std::pair<double, double> weibull_regression_fallback(const std::vector<double>& x, const std::vector<int>& r);
std::pair<std::vector<std::vector<double>>, int> cov_weibull_asymp_eff(const std::vector<double>& x, const std::vector<int>& r, double c, double b);


Sample read_input_normal(const std::string& tag);
Sample read_input_weibull(const std::string& tag);

#endif
