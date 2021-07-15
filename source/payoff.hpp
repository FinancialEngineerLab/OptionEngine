#pragma once

/*
class Payoff
{
public:
virtual ~Payoff() {}
virtual double operator() (double s) =0;
};
*/
#define MAX(x, y) ((x > y) ? x : y)

enum OptionType { Call = 1, Put = -1 };

class payOff
{
public:
	payOff() {}
    virtual double operator() (double spot)=0;
    virtual ~payOff() {}
    

private:

};

class VanillaPayOff : public payOff
{
public:
	VanillaPayOff(double strike, OptionType type) :
		strike_(strike), type_(type) {}

    virtual double operator() (double spot) override
	{
		return MAX(type_*(spot - strike_), 0.0);
	}

	virtual ~VanillaPayOff() {};

private:
	double strike_;
	OptionType type_;
};

class DigitalPayOff : public payOff
{
public:
	DigitalPayOff(double strike, OptionType type) :
		strike_(strike), type_(type) {}

	virtual double operator() (double spot) override
	{
		return (type_*(spot - strike_) >= 0.0) ? 1.0 : 0.0;
	}

	virtual ~DigitalPayOff() {};
private:
	double strike_;
	OptionType type_;
};
