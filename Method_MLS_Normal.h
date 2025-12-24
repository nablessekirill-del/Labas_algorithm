#ifndef METHOD_MLS_NORMAL_H
#define METHOD_MLS_NORMAL_H

#include "AbstractMethod.h"
#include "analysis.h"
#include <cmath>
#include <numeric>
#include <algorithm>
#include <vector>
#include <QString>

class Method_MLS_Normal : public AbstractMethod {
private:
    double mu_mls = 0, sigma_mls = 1;
    double s_res = 0, sum_w = 0, mean_z = 0, SS_z = 0;
    std::vector<double> lastData;
    std::vector<int> lastCens;

    void local_km(const std::vector<double>& x, const std::vector<int>& r,
                  std::vector<double>& ycum, std::vector<double>& fcum) {
        int n = (int)x.size();
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
        int m = (int)ycum.size();
        if (m < 3) return "Ошибка: мало данных";

        std::vector<double> zi(m), wi(m);
        sum_w = 0; double sum_wz = 0;
        for (int i = 0; i < m; ++i) {
            zi[i] = norm_ppf((i + 1.0 - 0.375) / (m + 0.25));
            wi[i] = 1.0 / std::max(1e-12, fcum[i] * (1.0 - fcum[i]));
            sum_w += wi[i]; sum_wz += wi[i] * zi[i];
        }
        mean_z = sum_wz / sum_w;
        SS_z = 0; double sum_wzx = 0, sum_wx = 0;
        for (int i = 0; i < m; ++i) {
            SS_z += wi[i] * std::pow(zi[i] - mean_z, 2);
            sum_wzx += wi[i] * zi[i] * ycum[i];
            sum_wx += wi[i] * ycum[i];
        }
        sigma_mls = (sum_wzx - (sum_wz * sum_wx) / sum_w) / SS_z;
        mu_mls = (sum_wx / sum_w) - sigma_mls * mean_z;

        double rss = 0;
        for (int i = 0; i < m; ++i) rss += wi[i] * std::pow(ycum[i] - (mu_mls + sigma_mls * zi[i]), 2);
        s_res = std::sqrt(rss / (m - 2));

        return QString("МНК (Нормальное)\nmu: %1\nsigma: %2").arg(mu_mls).arg(sigma_mls);
    }

    std::vector<GraphSeriesData> getGraphData() override {
        std::vector<GraphSeriesData> series;
        std::vector<double> ycum, fcum;
        local_km(lastData, lastCens, ycum, fcum);

        GraphSeriesData points, censored, line, low, up;
        points.name = "Events"; points.isScatter = true;
        censored.name = "Censored"; censored.isScatter = true;
        line.name = "MLE Линия"; low.name = "Lower CI"; up.name = "Upper CI";

        for(size_t i=0; i<ycum.size(); ++i) {
            points.x.push_back(ycum[i]);
            points.y.push_back(5.0 + norm_ppf(fcum[i]));
        }
        for(size_t i=0; i<lastData.size(); ++i) {
            if (lastCens[i] == 1) {
                censored.x.push_back(lastData[i]);
                censored.y.push_back(5.0 + (lastData[i] - mu_mls) / sigma_mls);
            }
        }
        for(int i = 0; i <= 100; ++i) {
            double p = 0.001 + i * 0.998 / 100.0;
            double z = norm_ppf(p);
            double x_hat = mu_mls + sigma_mls * z;
            double se = s_res * std::sqrt(1.0 + 1.0/sum_w + std::pow(z - mean_z, 2) / SS_z);
            line.x.push_back(x_hat); line.y.push_back(5.0 + z);
            low.x.push_back(x_hat - 2.5 * se); low.y.push_back(5.0 + z);
            up.x.push_back(x_hat + 2.5 * se);  up.y.push_back(5.0 + z);
        }
        series.push_back(points); series.push_back(censored);
        series.push_back(line); series.push_back(low); series.push_back(up);
        return series;
    }
};
#endif
