#include "IRengine.hpp"

double SimpleIRengine::singleSimplePeriod(double value)
{
    double f= value*(1+this->rate_);
    return f;
}

double CompoundIREngine::multiCompoundingPeriod(double value, int numPeriods)
{
    double f = value * std::pow(1+rate_, numPeriods);
    return f;
}

double CompoundIREngine::Continuous(double value, int numPeriods)
{
    double f = value * std::exp(rate_ * numPeriods);
    return f;
}