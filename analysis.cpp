#include "analysis.h"
#include <boost/math/distributions/normal.hpp>
#include <algorithm>
#include <numeric>
#include <fstream>
#include <iostream>

OptData globalOptData;
const std::vector<double>* G_X = nullptr;
const std::vector<int>* G_R = nullptr;

double norm_pdf(double z) {
    return 0.3989422804014327 * std::exp(-0.5 * z * z);
}

double norm_cdf(double z) {
    boost::math::normal_distribution<> N(0.0, 1.0);
    return boost::math::cdf(N, z);
}

double norm_ppf(double p) {
    if (p <= 0.0) p = 1e-15;
    if (p >= 1.0) p = 1.0 - 1e-15;
    boost::math::normal_distribution<> N(0.0, 1.0);
    return boost::math::quantile(N, p);
}

EmpiricalKM kaplan_meier_Itype(const std::vector<double>& x, const std::vector<int>& r) {
    int n = static_cast<int>(x.size());
    if (n == 0) return {};
    std::vector<std::pair<double, int>> data(n);
    for(int i=0; i<n; ++i) data[i] = {x[i], r[i]};
    std::sort(data.begin(), data.end());

    EmpiricalKM res;
    double S = 1.0;
    for (int i = 0; i < n; ++i) {
        double n_at_risk = static_cast<double>(n - i);
        if (data[i].second == 0) { // 0 - ОТКАЗ
            S *= (1.0 - 1.0 / n_at_risk);
            res.x_sorted.push_back(data[i].first);
            res.F_emp.push_back(1.0 - S);
        }
    }
    return res;
}

// РЕГРЕССИОННЫЙ ФОЛБЭК ДЛЯ ВЕЙБУЛЛА
std::pair<double, double> weibull_regression_fallback(const std::vector<double>& x, const std::vector<int>& r) {
    auto emp = kaplan_meier_Itype(x, r);
    std::vector<double> X_log, Y_log;
    for (size_t i = 0; i < emp.x_sorted.size(); ++i) {
        double F = emp.F_emp[i];
        if (F > 1e-6 && F < 0.999) {
            X_log.push_back(std::log(emp.x_sorted[i]));
            Y_log.push_back(std::log(std::log(1.0 / (1.0 - F))));
        }
    }
    if (X_log.size() < 2) return {10000.0, 1.0};
    double Sx=0, Sy=0, Sxx=0, Sxy=0, nn = X_log.size();
    for(size_t i=0; i<nn; ++i) {
        Sx += X_log[i]; Sy += Y_log[i];
        Sxx += X_log[i]*X_log[i]; Sxy += X_log[i]*Y_log[i];
    }
    double b = (nn*Sxy - Sx*Sy)/(nn*Sxx - Sx*Sx);
    double a = (Sy - b*Sx)/nn;
    return { std::exp(-a/b), b };
}

std::pair<double, double> weibull_mle_2par(const std::vector<double>& x, const std::vector<int>& r) {
    std::vector<double> obs;
    for(size_t i=0; i<x.size(); ++i) if(r[i] == 0) obs.push_back(x[i]);
    if(obs.size() < 2) return weibull_regression_fallback(x, r);

    double b = weibull_regression_fallback(x, r).second;
    for(int it=0; it<100; ++it) {
        double S0=0, S1=0, S2=0, L=0;
        for(double v : obs) {
            double lv = std::log(v), vb = std::pow(v, b);
            S0 += vb; S1 += vb*lv; S2 += vb*lv*lv; L += lv;
        }
        double f = (L/obs.size()) - (S1/S0) + 1.0/b;
        double df = -(S2*S0 - S1*S1)/(S0*S0) - 1.0/(b*b);
        double step = f/df;
        b -= step;
        if(std::abs(step) < 1e-7) break;
    }
    double s0 = 0;
    for(double v : obs) s0 += std::pow(v, b);
    return { std::pow(s0/obs.size(), 1.0/b), b };
}

// КОВАРИАЦИЯ ВЕЙБУЛЛА
std::pair<std::vector<std::vector<double>>, int> cov_weibull_asymp_eff(const std::vector<double>& x, const std::vector<int>& r, double c, double b) {
    int n_eff = 0;
    for (int ri : r) if (ri == 0) n_eff++;
    if (n_eff < 2) return {{{1,0},{0,1}}, n_eff};

    double p = std::log(c), q = 1.0/b, Jpp=0, Jpq=0, Jqq=0;
    for (size_t i=0; i<x.size(); ++i) {
        if (r[i]!=0 || x[i]<=0) continue;
        double z = (std::log(x[i])-p)/q;
        Jpp += 1.0; Jpq += z; Jqq += 1.0 + z*z;
    }
    double det = (Jpp*Jqq - Jpq*Jpq);
    if(std::abs(det) < 1e-15) det = 1e-15;
    std::vector<std::vector<double>> V = {
        {(Jqq/det)*(q*q/n_eff), (-Jpq/det)*(q*q/n_eff)},
        {(-Jpq/det)*(q*q/n_eff), (Jpp/det)*(q*q/n_eff)}
    };
    return {V, n_eff};
}

Sample read_input_normal(const std::string& tag) {
    std::ifstream in("Inp/" + tag + ".inp");
    Sample S; std::string tmp;
    if (!in) return S;
    in >> tmp >> S.n;
    while(in >> tmp && tmp != "X" && tmp != "Data");
    for(int i=0; i<S.n; ++i) {
        std::string val; in >> val;
        if(val == ",") { i--; continue; }
        try { S.x.push_back(std::stod(val)); } catch(...) {}
    }
    while(in >> tmp && (tmp != "R" && tmp != "Censorizes"));
    for(int i=0; i<S.n; ++i) {
        std::string val; in >> val;
        if(val == ",") { i--; continue; }
        try { S.r.push_back(std::stoi(val)); } catch(...) {}
    }
    return S;
}

Sample read_input_weibull(const std::string& tag) {
    return read_input_normal(tag);
}

void calculate_weibull_intervals(PlotData& pd, const std::vector<double>& p_vec, double beta) {
    double c = pd.param1;
    double b = pd.param2;


    double alpha = 1.0 - beta;
    double u_gamma = norm_ppf(1.0 - alpha / 2.0);


    auto cov = pd.cov;
    if (cov.empty() || cov.size() < 2) return;

    pd.x_low.clear();
    pd.x_est.clear();
    pd.x_up.clear();

    for (double p : p_vec) {
        double w = std::log(std::log(1.0 / (1.0 - p)));
        double xp = c * std::pow(-std::log(1.0 - p), 1.0 / b);
        pd.x_est.push_back(xp);

        // Дисперсия логарифма квантиля (метод дельта)
        double var_log_xp = cov[0][0] + w * w * cov[1][1] + 2.0 * w * cov[0][1];
        double std_log_xp = std::sqrt(std::max(0.0, var_log_xp));

        // Границы (логарифмически нормальный интервал)
        pd.x_low.push_back(xp * std::exp(-u_gamma * std_log_xp));
        pd.x_up.push_back(xp * std::exp(u_gamma * std_log_xp));
    }
}
