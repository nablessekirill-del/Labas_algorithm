#ifndef METHOD_MLE_NORMAL_H
#define METHOD_MLE_NORMAL_H

#include "AbstractMethod.h"
#include "analysis.h"
#include <cmath>
#include <numeric>
#include <algorithm>
#include <QDateTime>

class Method_MLE_Normal : public AbstractMethod {
private:
    double mu = 0;
    double sigma = 1;
    std::vector<double> lastData;
    std::vector<int> lastCens;

public:
    bool hasGraph() override { return true; }

    QString calculate(const std::vector<double>& data, const std::vector<int>& cens) override {
        lastData = data;
        lastCens = cens;
        if (data.empty()) return "Error: No data";


        int n = data.size();
        double sum = std::accumulate(data.begin(), data.end(), 0.0);
        mu = sum / n;
        double sq_sum = 0;
        for(double x : data) sq_sum += (x - mu) * (x - mu);
        sigma = std::sqrt(sq_sum / n);


        QString out;
        out += "Method:MLE_Normal\n";
        out += QString("n=%1\n").arg(n);

        out += "X\n";
        for(double x : data) out += QString::number(x, 'f', 5) + " , ";
        out += "\n";

        out += "R\n";
        for(int r : cens) out += QString::number(r) + " , ";
        out += "\n";

        out += QString("a_hat=%1\n").arg(mu, 0, 'f', 12);
        out += QString("sigma_hat=%1\n").arg(sigma, 0, 'f', 12);

        // Матрица ковариации (Cov[a, s]) - примерные значения, подставь свой расчет
        double var_a = (sigma * sigma) / n;
        double var_s = (sigma * sigma) / (2.0 * n);
        out += "Cov[a,s]:\n";
        out += QString("%1 %2\n").arg(var_a, 0, 'f', 12).arg(0.0, 0, 'f', 12);
        out += QString("%1 %2\n").arg(0.0, 0, 'f', 12).arg(var_s, 0, 'f', 12);

        std::vector<double> probs = {0.005, 0.01, 0.025, 0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 0.8, 0.9, 0.95, 0.975, 0.99, 0.995};
        out += "P\n";
        for(double p : probs) out += QString::number(p, 'f', 12) + " ; ";
        out += "\n";

        QString xp_low, xp_mid, xp_up;
        for(double p : probs) {
            double zp = norm_ppf(p);
            double val = mu + zp * sigma;

            // расчет доверительного интервала для квантили
            double se = (sigma / std::sqrt((double)n)) * std::sqrt(1.0 + 0.5 * zp * zp);

            xp_low += QString::number(val - 1.96 * se, 'f', 12) + " ; ";
            xp_mid += QString::number(val, 'f', 12) + " ; ";
            xp_up  += QString::number(val + 1.96 * se, 'f', 12) + " ; ";
        }

        out += "Xp_low\n" + xp_low + "\n";
        out += "Xp\n" + xp_mid + "\n";
        out += "Xp_up\n" + xp_up + "\n";

        return out;
    }

    bool hasGraph() const { return true; }

    std::vector<GraphSeriesData> getGraphData() override {
        std::vector<GraphSeriesData> res;
        if (lastData.empty()) return res;

        std::vector<std::pair<double, int>> pairedData;
        for(size_t i = 0; i < lastData.size(); ++i) {
            pairedData.push_back({lastData[i], lastCens[i]});
        }
        std::sort(pairedData.begin(), pairedData.end());

        size_t n = pairedData.size();

        //точки
        GraphSeriesData dots_normal;
        dots_normal.name = "Events";
        dots_normal.isScatter = true;

        GraphSeriesData dots_cens;
        dots_cens.name = "Censored";
        dots_cens.isScatter = true;

        for(size_t i = 0; i < n; ++i) {
            double p = (i + 0.375) / (n + 0.25);
            double x_val = pairedData[i].first;
            double y_val = 5.0 + norm_ppf(p);

            if (pairedData[i].second == 0) {
                dots_normal.x.push_back(x_val);
                dots_normal.y.push_back(y_val);
            } else {
                dots_cens.x.push_back(x_val);
                dots_cens.y.push_back(y_val);
            }
        }
        res.push_back(dots_normal);
        res.push_back(dots_cens);

        // линии (MLE + CI)
        GraphSeriesData line, ci_up, ci_low;
        line.name = "MLE Линия";
        ci_up.name = "95% CI";
        ci_low.name = "CI_low"; // вспомогательное имя
        line.isScatter = ci_up.isScatter = ci_low.isScatter = false;

        double x_start = pairedData.front().first;
        double x_end = pairedData.back().first;
        double step = (x_end - x_start) / 100.0;

        for(double x = x_start; x <= x_end + step/2.0; x += step) {
            double z = (x - mu) / sigma;
            double y_center = 5.0 + z;


            double se = (sigma / std::sqrt((double)n)) * std::sqrt(1.0 + 0.5 * z * z);
            double delta = 1.96 * (se / sigma);

            line.x.push_back(x);
            line.y.push_back(y_center);
            ci_up.x.push_back(x);
            ci_up.y.push_back(y_center + delta);
            ci_low.x.push_back(x);
            ci_low.y.push_back(y_center - delta);
        }

        res.push_back(line);
        res.push_back(ci_up);
        res.push_back(ci_low);

        return res;
    }
};
#endif
