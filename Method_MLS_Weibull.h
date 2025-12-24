#ifndef METHOD_MLS_WEIBULL_H
#define METHOD_MLS_WEIBULL_H

#include "AbstractMethod.h"
#include <cmath>
#include <numeric>
#include <algorithm>
#include <vector>
#include <QString>

class Method_MLS_Weibull : public AbstractMethod {
private:
    double a_w = 0, b_w = 1;
    double s_res = 0, sum_w = 0, mean_z = 0, SS_z = 0;
    std::vector<double> lastData;
    std::vector<int> lastCens;

    double weibull_z(double p) { return std::log(-std::log(std::max(1e-12, 1.0 - p))); }

    void local_km(const std::vector<double>& x, const std::vector<int>& r,
                  std::vector<double>& ycum, std::vector<double>& fcum) {
        int n = x.size();
        std::vector<int> idx(n); std::iota(idx.begin(), idx.end(), 0);
        std::sort(idx.begin(), idx.end(), [&](int i, int j) { return x[i] < x[j]; });
        double S = 1.0; int at_risk = n;
        ycum.clear(); fcum.clear();
        for (int i = 0; i < n;) {
            double xi = x[idx[i]];
            int d = 0, c = 0, j = i;
            while (j < n && x[idx[j]] == xi) { if (r[idx[j]] == 0) d++; else c++; ++j; }
            if (d > 0) {
                S *= (double)(at_risk - d) / (double)at_risk;
                ycum.push_back(xi); fcum.push_back(1.0 - S);
            }
            at_risk -= (d + c); i = j;
        }
    }

public:
    bool hasGraph() override { return true; }

    QString calculate(const std::vector<double>& data, const std::vector<int>& cens) override {
        lastData = data; lastCens = cens;
        std::vector<double> ycum, fcum;
        local_km(data, cens, ycum, fcum);
        int m = ycum.size();
        if (m < 3) return "Ошибка: мало данных";

        std::vector<double> zi(m), lnx(m);
        sum_w = 0; double sum_wz = 0, sum_wzx = 0, sum_wx = 0;
        for (int i = 0; i < m; ++i) {
            zi[i] = weibull_z((i + 1.0 - 0.3) / (m + 0.4));
            lnx[i] = std::log(ycum[i]);
            sum_w += 1.0; sum_wz += zi[i];
            sum_wx += lnx[i]; sum_wzx += zi[i] * lnx[i];
        }
        mean_z = sum_wz / sum_w;
        SS_z = 0;
        for (int i = 0; i < m; ++i) SS_z += std::pow(zi[i] - mean_z, 2);
        double sigma_hat = (sum_wzx - (sum_wz * sum_wx) / sum_w) / SS_z;
        double mu_hat = (sum_wx / sum_w) - sigma_hat * mean_z;
        a_w = 1.0 / sigma_hat; b_w = std::exp(mu_hat);
        double rss = 0;
        for (int i = 0; i < m; ++i) rss += std::pow(lnx[i] - (mu_hat + sigma_hat * zi[i]), 2);
        s_res = std::sqrt(rss / (m - 2));

        return QString("МНК (Вейбулл)\na: %1, b: %2").arg(a_w).arg(b_w);
    }

    std::vector<GraphSeriesData> getGraphData() override {
        std::vector<GraphSeriesData> series;
        std::vector<double> ycum, fcum;
        local_km(lastData, lastCens, ycum, fcum);

        GraphSeriesData points, censored, line, low, up;
        points.name = "Events"; points.isScatter = true;
        censored.name = "Censored"; censored.isScatter = true;
        line.name = "MLE Линия"; low.name = "Lower CI"; up.name = "Upper CI";

        // Рисуем в логарифмах, чтобы на линейной оси MainWindow была прямая
        for(size_t i=0; i<ycum.size(); ++i) {
            points.x.push_back(std::log(ycum[i]));
            points.y.push_back(weibull_z(fcum[i]));
        }
        double mu = std::log(b_w), sigma = 1.0/a_w;
        for(size_t i=0; i<lastData.size(); ++i) {
            if (lastCens[i] == 1) {
                double lnx_val = std::log(lastData[i]);
                censored.x.push_back(lnx_val);
                censored.y.push_back((lnx_val - mu) / sigma);
            }
        }
        for(int i = 0; i <= 100; ++i) {
            double p = 0.005 + i * 0.99 / 100.0;
            double z = weibull_z(p);
            double lnx_hat = mu + sigma * z;
            double se = s_res * std::sqrt(1.0 + 1.0/sum_w + std::pow(z - mean_z, 2) / SS_z);
            line.x.push_back(lnx_hat); line.y.push_back(z);
            low.x.push_back(lnx_hat - 2.5 * se); low.y.push_back(z);
            up.x.push_back(lnx_hat + 2.5 * se);  up.y.push_back(z);
        }
        series.push_back(points); series.push_back(censored);
        series.push_back(line); series.push_back(low); series.push_back(up);
        return series;
    }
};
#endif
