#ifndef METHOD_ANOVA_H
#define METHOD_ANOVA_H

#include "AbstractMethod.h"
#include <vector>
#include <cmath>
#include <numeric>
#include <QString>
#include <iomanip>
#include <boost/math/distributions/fisher_f.hpp>

class Method_Anova : public AbstractMethod {
public:
    bool hasGraph() override { return false; }

    QString calculate(const std::vector<double>& data, const std::vector<int>& cens) override {
        Q_UNUSED(cens);

        if (data.empty()) return "Ошибка: Входные данные пусты";

        int k = static_cast<int>(data[0]);
        if (k < 2) return "Ошибка: Для ANOVA требуется минимум 2 группы";

        struct Group {
            std::vector<double> values;
            double mean;
            double var; // выборочная дисперсия (s^2)
            int n;
        };

        std::vector<Group> groups;
        size_t currentIdx = 1;
        int totalN = 0;
        double sumAll = 0;

        for (int i = 0; i < k; ++i) {
            if (currentIdx >= data.size()) break;
            int n_i = static_cast<int>(data[currentIdx++]);
            Group g;
            g.n = n_i;
            totalN += n_i;

            double groupSum = 0;
            for (int j = 0; j < n_i; ++j) {
                double val = data[currentIdx++];
                g.values.push_back(val);
                groupSum += val;
            }
            g.mean = groupSum / n_i;
            sumAll += groupSum;

            double sqDiffSum = 0;
            for (double val : g.values) {
                sqDiffSum += (val - g.mean) * (val - g.mean);
            }
            g.var = (n_i > 1) ? sqDiffSum / (n_i - 1) : 0;
            groups.push_back(g);
        }

        // ANOVA
        double generalMean = sumAll / totalN;

        // Межгрупповая дисперсия (S_out)
        double ssOut = 0;
        for (const auto& g : groups) {
            ssOut += g.n * (g.mean - generalMean) * (g.mean - generalMean);
        }
        int df1 = k - 1;
        double sOut = ssOut / df1;

        // Внутригрупповая дисперсия (S_in)
        double ssIn = 0;
        for (const auto& g : groups) {
            ssIn += (g.n - 1) * g.var;
        }
        int df2 = totalN - k;
        double sIn = ssIn / df2;

        // F-статистика
        double f_obs = sOut / sIn;
        double alpha = 0.05; // По умолчанию из твоего ТЗ

        boost::math::fisher_f_distribution<double> dist(df1, df2);
        double f_crit = boost::math::quantile(boost::math::complement(dist, alpha));

        // отчет
        QString res = "Результаты однофакторного дисперсионного анализа (ANOVA)\n\n";
        res += QString("Уровень значимости alpha = %1\n").arg(alpha);
        res += QString("Число групп (выборок) k = %1\n").arg(k);
        res += QString("Общая средняя X.. = %1\n\n").arg(generalMean, 0, 'f', 10);

        for (int i = 0; i < groups.size(); ++i) {
            res += QString("Группа %1: n = %2, mean = %3, stdDev = %4\n")
                       .arg(i + 1)
                       .arg(groups[i].n)
                       .arg(groups[i].mean, 0, 'f', 10)
                       .arg(std::sqrt(groups[i].var), 0, 'f', 10);
        }

        res += QString("\nМежгрупповая дисперсия S_out = %1\n").arg(sOut, 0, 'f', 10);
        res += QString("Внутригрупповая дисперсия S_in  = %1\n").arg(sIn, 0, 'f', 10);
        res += QString("Наблюдаемое значение F = %1\n").arg(f_obs, 0, 'f', 10);
        res += QString("Степени свободы: f1 = %1, f2 = %2\n").arg(df1).arg(df2);
        res += QString("Критическое значение F_(1-alpha) = %1\n\n").arg(f_crit, 0, 'f', 10);

        res += "Решение: ";
        if (f_obs <= f_crit) {
            res += "H0 принимается (различия между средними статистически несущественны).";
        } else {
            res += "H0 отвергается (различия между средними статистически значимы).";
        }

        return res;
    }

    std::vector<GraphSeriesData> getGraphData() override { return {}; }
};

#endif
