#ifndef METHOD_MLE_WEIBULL_H
#define METHOD_MLE_WEIBULL_H

#include "AbstractMethod.h"
#include "analysis.h"
#include <cmath>
#include <numeric>
#include <algorithm>
#include <vector>

class Method_MLE_Weibull : public AbstractMethod {
private:
    double c_hat = 0;
    double b_hat = 1;
    std::vector<double> lastData;
    std::vector<int> lastCens;

public:
    bool hasGraph() override { return true; }

    QString calculate(const std::vector<double>& data, const std::vector<int>& cens) override {
        lastData = data;
        lastCens = cens;
        if (data.empty()) return "Error: No data";


        auto est = weibull_mle_2par(data, cens);
        c_hat = est.first;
        b_hat = est.second;

        if (!(std::isfinite(c_hat) && std::isfinite(b_hat)) || c_hat <= 0.0 || b_hat <= 0.0) {
            auto est_fb = weibull_regression_fallback(data, cens);
            c_hat = est_fb.first; b_hat = est_fb.second;
        }

        int n = data.size();
        QString out;
        out += "Method:MLE_Weibull\n";
        out += QString("n=%1\n").arg(n);

        out += "X\n";
        for(double x : data) out += QString::number(x, 'f', 5) + " , ";
        out += "\n";

        out += "R\n";
        for(int r : cens) out += QString::number(r) + " , ";
        out += "\n";

        out += QString("c_hat=%1\n").arg(c_hat, 0, 'f', 12);
        out += QString("b_hat=%1\n").arg(b_hat, 0, 'f', 12);

        // Матрица ковариации
        auto cv = cov_weibull_asymp_eff(data, cens, c_hat, b_hat);
        out += "Cov[c,b]:\n";
        out += QString("%1 %2\n").arg(cv.first[0][0], 0, 'f', 12).arg(cv.first[0][1], 0, 'f', 12);
        out += QString("%1 %2\n").arg(cv.first[1][0], 0, 'f', 12).arg(cv.first[1][1], 0, 'f', 12);

        // Блок P и расчет квантилей
        std::vector<double> probs = {0.005, 0.01, 0.025, 0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 0.8, 0.9, 0.95, 0.975, 0.99, 0.995};
        out += "P\n";
        for(double p : probs) out += QString::number(p, 'f', 12) + " ; ";
        out += "\n";

        QString xp_low, xp_mid, xp_up;
        for(double p : probs) {
            double tp = std::log(-std::log(1.0 - p));
            double val = c_hat * std::exp(tp / b_hat); // Квантиль Xp

            // SE для Вейбулла
            double z = tp;
            double se_y = (1.0 / b_hat) * std::sqrt((1.0 + 0.5 * z * z) / (double)n);

            xp_low += QString::number(val * std::exp(-1.96 * se_y), 'f', 12) + " ; ";
            xp_mid += QString::number(val, 'f', 12) + " ; ";
            xp_up  += QString::number(val * std::exp(1.96 * se_y), 'f', 12) + " ; ";
        }

        out += "Xp_low\n" + xp_low + "\n";
        out += "Xp\n" + xp_mid + "\n";
        out += "Xp_up\n" + xp_up + "\n";

        return out;
    }

    std::vector<GraphSeriesData> getGraphData() override {
        std::vector<GraphSeriesData> res;
        if (lastData.empty()) return res;

        std::vector<std::pair<double, int>> pairedData;
        for(size_t i = 0; i < lastData.size(); ++i) pairedData.push_back({lastData[i], lastCens[i]});
        std::sort(pairedData.begin(), pairedData.end());

        size_t n = pairedData.size();

        //точки
        GraphSeriesData dots_ev, dots_cens;
        dots_ev.name = "Events"; dots_ev.isScatter = true;
        dots_cens.name = "Censored"; dots_cens.isScatter = true;

        for(size_t i = 0; i < n; ++i) {
            double p = (i + 0.3) / (n + 0.4); // Агамировское смещение
            double x_val = pairedData[i].first;
            // Спрямляющая ось Y для Вейбулла: 5 + ln(-ln(1-p))
            double y_val = 5.0 + std::log(-std::log(1.0 - p));

            if (pairedData[i].second == 0) {
                dots_ev.x.push_back(x_val); dots_ev.y.push_back(y_val);
            } else {
                dots_cens.x.push_back(x_val); dots_cens.y.push_back(y_val);
            }
        }
        res.push_back(dots_ev); res.push_back(dots_cens);

        // линии
        GraphSeriesData line, ci_up, ci_low;
        line.name = "MLE Линия"; ci_up.name = "95% CI"; ci_low.name = "CI_low";
        line.isScatter = ci_up.isScatter = ci_low.isScatter = false;

        double x_start = pairedData.front().first;
        double x_end = pairedData.back().first;
        double step = (x_end - x_start) / 100.0;

        for(double x = x_start; x <= x_end + step/2.0; x += step) {
            if (x <= 0) continue;
            // z для Вейбулла в спрямляющих координатах
            double z = b_hat * std::log(x / c_hat);
            double y_center = 5.0 + z;

            // Расчет Delta строго по твоей логике из Normal
            // Коэффициент (1 + 0.5*z*z) создает тот самый раструб (расширение к краям)
            double se = (1.0 / std::sqrt((double)n)) * std::sqrt(1.0 + 0.5 * z * z);
            double delta = 1.96 * se;

            line.x.push_back(x);     line.y.push_back(y_center);
            ci_up.x.push_back(x);    ci_up.y.push_back(y_center + delta);
            ci_low.x.push_back(x);   ci_low.y.push_back(y_center - delta);
        }

        res.push_back(line); res.push_back(ci_up); res.push_back(ci_low);
        return res;
    }
};

#endif
