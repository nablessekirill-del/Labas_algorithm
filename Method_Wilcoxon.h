#ifndef METHOD_WILCOXON_H
#define METHOD_WILCOXON_H

#include "AbstractMethod.h"
#include <vector>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <QString>
#include <boost/math/distributions/normal.hpp>

class Method_Wilcoxon : public AbstractMethod {
public:
    bool hasGraph() override { return false; }

    QString calculate(const std::vector<double>& data, const std::vector<int>& cens) override {
        Q_UNUSED(cens);

        if (data.size() < 5) return "Ошибка: Недостаточно данных";

        double alpha = data[0];
        int m1 = static_cast<int>(data[1]);
        int n1 = static_cast<int>(data[2]);

        if (data.size() < static_cast<size_t>(4 + m1 + n1))
            return "Ошибка: Размер данных не соответствует заголовку";

        std::vector<double> x, y;
        for (int i = 0; i < m1; ++i) x.push_back(data[4 + i]);
        for (int i = 0; i < n1; ++i) y.push_back(data[4 + m1 + i]);

        auto getStats = [](const std::vector<double>& v) {
            double sum = std::accumulate(v.begin(), v.end(), 0.0);
            double m = sum / v.size();
            double sq_sum = 0;
            for (double val : v) sq_sum += (val - m) * (val - m);
            return std::make_pair(m, std::sqrt(sq_sum / (v.size() - 1)));
        };
        auto stats1 = getStats(x);
        auto stats2 = getStats(y);

        struct Item { double val; int group; };
        std::vector<Item> united;
        for (double val : x) united.push_back({val, 1});
        for (double val : y) united.push_back({val, 2});

        std::sort(united.begin(), united.end(), [](const Item& a, const Item& b) {
            return a.val < b.val;
        });

        int m_small = (m1 <= n1) ? m1 : n1;
        int small_group = (m1 <= n1) ? 1 : 2;
        double W_obs = 0;

        for (size_t i = 0; i < united.size(); ++i) {
            if (united[i].group == small_group) {
                W_obs += (i + 1);
            }
        }


        bool useExact = (m1 + n1 <= 40);

        double mu_w = (static_cast<double>(m_small) * (m1 + n1 + 1)) / 2.0;
        double sigma_w = std::sqrt(static_cast<double>(m1 * n1 * (m1 + n1 + 1)) / 12.0);
        boost::math::normal_distribution<double> dist(0.0, 1.0);
        double z = boost::math::quantile(dist, 1.0 - alpha / 2.0);

        double W_low = std::round(mu_w - z * sigma_w);
        double W_up = std::round(mu_w + z * sigma_w);


        QString res = "================ Двусторонний критерий Уилкоксона ================\n\n";
        res += QString("alpha = %1\n").arg(alpha);
        res += QString("m1 = %1, n1 = %2, N = %3\n").arg(m1).arg(n1).arg(m1+n1);
        res += QString("Меньшая выборка: группа %1 (m = %2)\n\n").arg(small_group).arg(m_small);

        res += QString("Выборка 1: mean = %1, std = %2\n").arg(stats1.first, 0, 'f', 4).arg(stats1.second, 0, 'f', 4);
        res += QString("Выборка 2: mean = %1, std = %2\n\n").arg(stats2.first, 0, 'f', 4).arg(stats2.second, 0, 'f', 4);

        res += QString("W_obs (сумма рангов) = %1\n").arg(W_obs);
        res += QString("Режим вычисления: %1\n").arg(useExact ? "ТОЧНЫЙ" : "ПРИБЛИЖЕННЫЙ");
        res += QString("Критический интервал: [%1; %2]\n\n").arg(W_low).arg(W_up);

        res += "H0: распределения совпадают\nH1: распределения различаются\n\n";

        if (W_obs > W_low && W_obs < W_up) {
            res += "Решение: H0 не отвергается (различия не значимы).";
        } else {
            res += "Решение: H0 отвергается (различия статистически значимы).";
        }

        return res;
    }

    std::vector<GraphSeriesData> getGraphData() override { return {}; }
};

#endif
