#pragma once
#include "date.hpp"
#include <vector>

struct MarketTermStructure
{
    Date date;
    double yield, div, vol;
};

class TermStructure
{
    public:
    TermStructure();
    TermStructure(Date dates, double rates)
    {
        dates_.push_back(dates);
        rates_.push_back(rates);
    }
    void setPoint(Date dates, double rates);
    std::vector<Date> getDates();
    double getValue(Date dt);
    
    void operator+=(const double d);
    void operator*=(const double d);

    private:
    std::vector<Date> dates_;
    std::vector<double> rates_;

};

class YieldTermStructure : public TermStructure
{
    public:
    YieldTermStructure();
    ~YieldTermStructure();
    YieldTermStructure(Date dates, double rates)
    : TermStructure(dates, rates) {}

    double discountFactor(Date dt);
    double forwardRate(Date dt1, Date dt2);
    double zeroRate(Date dt);

    private:
};

class VolatilityTermStructure : public TermStructure
{
    public:
    VolatilityTermStructure() {};
    VolatilityTermStructure(Date dates, double rates):
    TermStructure(dates, rates) {}
    double getVariance(Date dt);
    
    private:

};

class DividendTermStructure : public TermStructure
{
    public:
    DividendTermStructure() {};
    DividendTermStructure(Date dates, double rates)
    : TermStructure(dates, rates){}

    private:
};
