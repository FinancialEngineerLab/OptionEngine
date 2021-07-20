#ifndef __AdditionalJob__IREngine__
#define __AdiitionalJob__IREngine__

#include "Date.hpp"
#include <iostream>
//#include <ql/quantlib.hpp>

class IREngine
{
    public:
    IREngine() {};
    IREngine(double rate):rate_(rate) {};
    IREngine(const IREngine &v):rate_(v.rate_) {};
    IREngine &operator = (const IREngine &v)
    {
        if(this!=&v)
        {
            this->rate_ = v.rate_;
        }
        return *this;
    }
    ~IREngine(){};

    //double singleSimplePeriod(double value);

    private:
    double rate_;
}

class SimpleIREngine::IREngine
{
    public:
    SimpleIREngine(){};
    SimpleIREngine(double rate):
    IREngine(rate) {};
    SimpleIREngine(const SimpleIREngine &v):rate(v.rate);
    SingleSimpleIREngine &operator = (const SingleSimpleIREngine &v)
    {
        if(&v!=this)
        {
            this->rate = v.rate;
        }
        return *this;
    };
    ~SimpleIREngine(){};

    double singleSimplePeriod(double value);

    private:
}

class CompoundIREngine::IREngine
{
    public:
    AnnualCompoundIREngine();
    AnnualCompoundIREngine(double rate):
    IREngine(rate){};
    AnnualCompoundIREngine(const AnnualCompoundIREngine& v):rate(v.rate);
    AnnualCompoundIREngine &operator = (const AnnualCompoundIREngine& v)
    {
        if(&v!=this)
        {
            this->rate = v.rate;
        }
        return *this;
    };
    ~AnnualCompoundIREngine(){};

    double multiCompoundingPeriod(double value, int numPeriods);
    double Coninuous(double value, int numPeriods);
    
    private:

}



/*
inline double IREngine::singleSimplePeriod(double value)
{
    double f = value * (1+this->rate_);
    return f;
}
*/



#endif