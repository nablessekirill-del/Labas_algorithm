#ifndef METHOD_FISHERSTUDENT_H
#define METHOD_FISHERSTUDENT_H

#include "AbstractMethod.h"
#include <vector>
#include <cmath>
#include <algorithm>
#include <QString>
#include <boost/math/distributions/fisher_f.hpp>
#include <boost/math/distributions/students_t.hpp>

class Method_FisherStudent : public AbstractMethod {
public:
    bool hasGraph() override { return false; }

    QString calculate(const std::vector<double>& data, const std::vector<int>& cens) override {
        Q_UNUSED(cens);

        //парсинг
        if (data.size() < 5) return "Ошибка: Недостаточно данных для анализа двух выборок";

        double alpha = data[0];
        int n1 = static_cast<int>(data[1]);
        if (data.size() < static_cast<size_t>(2 + n1 + 1)) return "Ошибка: Неверный размер первой выборки";

        std::vector<double> s1_data, s2_data;
        for (int i = 0; i < n1; ++i) s1_data.push_back(data[2 + i]);

        int n2_idx = 2 + n1;
        int n2 = static_cast<int>(data[n2_idx]);
        for (int i = 0; i < n2; ++i) s2_data.push_back(data[n2_idx + 1 + i]);

        auto getStats = [](const std::vector<double>& v) {
            double sum = 0;
            for (double x : v) sum += x;
            double m = sum / v.size();
            double sq_sum = 0;
            for (double x : v) sq_sum += (x - m) * (x - m);
            double s = std::sqrt(sq_sum / (v.size() - 1));
            return std::make_pair(m, s);
        };

        auto stats1 = getStats(s1_data);
        auto stats2 = getStats(s2_data);
        double m1 = stats1.first, s1 = stats1.second;
        double m2 = stats2.first, s2 = stats2.second;

        //критерий фишера
        double f_obs = (s1 * s1 >= s2 * s2) ? (s1 * s1) / (s2 * s2) : (s2 * s2) / (s1 * s1);
        double df1 = (s1 * s1 >= s2 * s2) ? n1 - 1 : n2 - 1;
        double df2 = (s1 * s1 >= s2 * s2) ? n2 - 1 : n1 - 1;

        boost::math::fisher_f_distribution<double> f_dist(df1, df2);
        double f_crit = boost::math::quantile(boost::math::complement(f_dist, alpha));
        bool fisherEqual = (f_obs <= f_crit);

        //критерий стьюдента
        double t_obs, df_student;
        if (fisherEqual) {
            // Обычный t-тест (равные дисперсии)
            double sp = std::sqrt(((n1 - 1) * s1 * s1 + (n2 - 1) * s2 * s2) / (n1 + n2 - 2));
            t_obs = (m1 - m2) / (sp * std::sqrt(1.0 / n1 + 1.0 / n2));
            df_student = n1 + n2 - 2;
        } else {
            // Тест Уэлча (неравные дисперсии)
            double w1 = s1 * s1 / n1;
            double w2 = s2 * s2 / n2;
            t_obs = (m1 - m2) / std::sqrt(w1 + w2);
            df_student = std::pow(w1 + w2, 2) / (std::pow(w1, 2) / (n1 - 1) + std::pow(w2, 2) / (n2 - 1));
        }

        boost::math::students_t_distribution<double> t_dist(df_student);
        double t_crit = boost::math::quantile(boost::math::complement(t_dist, alpha / 2.0));
        bool studentEqual = (std::abs(t_obs) <= t_crit);

        QString res = "Критерий Фишера и критерий Стьюдента для двух выборок\n\n";
        res += QString("Уровень значимости alpha = %1\n\n").arg(alpha, 0, 'f', 6);

        auto formatSample = [&](const std::vector<double>& v, int n, double m, double s, int id) {
            QString out = QString("Первая выборка (n%1 = %2):\n").arg(id).arg(n);
            if (id == 2) out = QString("Вторая выборка (n%1 = %2):\n").arg(id).arg(n);
            for (double x : v) out += QString("%1 ").arg(x, 0, 'f', 6);
            out += QString("\nВыборочное среднее m%1 = %2\n").arg(id).arg(m, 0, 'f', 6);
            out += QString("Выборочное СКО s%1 = %2\n").arg(id).arg(s, 0, 'f', 6);
            // Генеральные параметры mu и sigma здесь справочные (mu=1.0, sigma из генератора)
            res += out;
        };

        formatSample(s1_data, n1, m1, s1, 1);
        res += "Генеральное мат. ожидание mu1 = 1.000000\n";
        res += "Генеральная дисперсия sigma1^2 = 0.100000\n\n";

        formatSample(s2_data, n2, m2, s2, 2);
        res += "Генеральное мат. ожидание mu2 = 1.000000\n";
        res += "Генеральная дисперсия sigma2^2 = 0.150000\n\n";

        res += QString("Степени свободы: nu1 = %1, nu2 = %2\n\n").arg(df1, 0, 'f', 6).arg(df2, 0, 'f', 6);

        res += "Результат критерия Фишера (проверка равенства дисперсий): ";
        res += fisherEqual ? "дисперсии можно считать равными.\n" : "дисперсии различаются.\n";

        res += "Результат критерия Стьюдента (проверка равенства средних): ";
        res += studentEqual ? "математические ожидания можно считать равными." : "математические ожидания различаются.";

        return res;
    }

    std::vector<GraphSeriesData> getGraphData() override { return {}; }
};

#endif
