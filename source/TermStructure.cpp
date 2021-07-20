#include "termStructure.hpp"
#include <iostream>
#include <cmath>
#include <stdexcept>

void TermStructure::setPoint(Date dates, double rates)
{
    dates_.push_back(dates);
    rates_.push_back(rates);
}

std::vector<Date> TermStructure::getDates()
{
    return dates_;
}

double TermStructure::getLinearValue(Date dt)
{
    double inter_rate = 0;
    if(dt < dates_.at(0))
    {
        std::cout << " Out of Range ! " << std::endl;
        throw dt;
    }

    else if(dt == dates_.at(0))
    {
        inter_rate = rates_.at(0);
    }

    else
    {
        int i = 0;
        while (dt>dates_.at(i))
        {
            i++;
            if(i==dates_.size())
            {
                std::cout << " Out of Range ! " << std::endl;
                throw dt;
            }
        }
        Date dt2 = dates_.at(i);
        Date dt1 = dates_.at(i-1);
        double value2 = rates_.at(i);
        double value1 = rates_.at(i-1);

        double x1 = Date_Between(dt1, dt);
        double x2 = Date_Between(dt, dt2);

        inter_rate = x2/(x1+x2)*value1+
        x1 / (x1+x2)*value2;
    }
    return inter_rate;
}

void TermStructure::operator+=(double d)
{
    for(unsigned int i =0;i<rates_.size();i++)
    {
        rates_[i] += d;
    }
}


void TermStructure::operator+=(double d)
{
    for(unsigned int i =0;i<rates_.size();i++)
    {
        rates_[i] *= d;
    }
}

double YieldTermStructure::discountFactor(compounding comp = Continuous,Date dt)
{
    Date start_date = getDates().front();
    double t = Date_Between(start_date, dt);

    if (comp == Continuous)
    {
        
    }
}

double YeldTermStructure::forwardRate(Date dt1, Date dt2)
{
    double t = Date_Between(dt1, dt2);
    double forward_rate=
    365/t*std::log(discountFactor(dt1)/discountFactor(dt2));
    return forward_rate;
}