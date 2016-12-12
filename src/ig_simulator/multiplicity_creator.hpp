#pragma once

#include "repertoire.hpp"

template<class IgVariableRegionPtr>
class BaseMultiplicityCreator {
public:
    virtual size_t AssignMultiplicity(IgVariableRegionPtr ig_variable_region_ptr) = 0;
};

//-------------------------------------------------------------------------------------------

template<class IgVariableRegionPtr>
class ExponentialMultiplicityCreator : public BaseMultiplicityCreator<IgVariableRegionPtr> {
    double lambda_;

public:
    ExponentialMultiplicityCreator(size_t n1, size_t n2) :
            lambda_(double(n1) / n2) { }

    size_t AssignMultiplicity(IgVariableRegionPtr ig_variable_region_ptr) {
        return size_t((-1) * log(1 - double(rand()) / RAND_MAX) / lambda_) + 1;
    }
};

typedef ExponentialMultiplicityCreator<HC_VariableRegionPtr> HC_ExponentialMultiplicityCreator;
typedef ExponentialMultiplicityCreator<LC_VariableRegionPtr> LC_ExponentialMultiplicityCreator;

//-------------------------------------------------------------------------------------------

template<class IgVariableRegionPtr>
class PowerLawMultiplicityCreator : public BaseMultiplicityCreator<IgVariableRegionPtr> {
    double T_;
    double lambda_;
    double C_;
    double mult_;

    double FindLambda(double M, double N) {
        double l = 1;
        double r = 1000;
        for (size_t i = 0; i < 100; i++) {
            double s = (l + r) / 2;
            double tmp = (s - 1) / (s - 2) * (1 - pow(T_, 2 - s)) / (1 - pow(T_, 1 - s));
            //cout << tmp << " " << s << endl;
            if(tmp < M / N)
                r = s;
            else
                l = s;
        }
        //cout << l - 2 << endl;
        return l; //(N + 2) / (N + 1);
    }

    double s(double n) {
        return (1 - pow(n + 1, 1 - lambda_)) * C_ / (lambda_ - 1);
    }

    double p(double n) {
        return pow(n, - lambda_) * C_;
    }

public:
    PowerLawMultiplicityCreator(size_t n1, size_t n2) :
            T_(double(n2) * 0.3),
            lambda_(FindLambda(n2, n1)),
            C_((lambda_ - 1) / (1 - pow(double(T_), 1 - lambda_))),
            mult_(1 / (1 - lambda_)) {
//        srand (time(NULL));
//        cout << lambda_ << " " << C_ << endl;
//        cout << (lambda_ + 1) / lambda_ * (1 - pow(T_, -lambda_)) << endl;
        double sum0 = 0;
        double sum1 = 0;
        double div = 0;
        for(double i = T_; i > .5; i -= 1) {
            sum0 += p(i);
            sum1 += (s(i) - s(i - 1)) * i;
            div += (s(i) - s(i - 1)) * i * i;
        }
//        cout << sum0 << " " << sum1 << " " << sqrt(div - sum1 * sum1) << endl;
    }

    size_t AssignMultiplicity(IgVariableRegionPtr ig_variable_region_ptr) {
        double base = 1 - double(rand()) / RAND_MAX * (-1 + lambda_) / C_; //double(i + .5) / n1_;
        return size_t(pow(base, mult_));
    }
};

typedef PowerLawMultiplicityCreator<HC_VariableRegionPtr> HC_PowerLawMultiplicityCreator;
typedef PowerLawMultiplicityCreator<LC_VariableRegionPtr> LC_PowerLawMultiplicityCreator;