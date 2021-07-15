#pragma once
#include <iostream>
#include <cmath>
#include <random>
#include <vector>

#include "date.hpp"
#include "process.hpp"
#include "payoff.hpp"
#include "normal.hpp"


using namespace std;

//enum OptionType { Call = 1, Put = -1 };
enum OptionProduct { Vanilla = 1, Digital = 0 };
enum CalcMethod { MC=0, FDM=1, BSM=2};
enum Solver { Thomas = 1, SOR = 2};

class optionOutput
{
public:
	optionOutput(double strike, OptionType type) :
		strike_(strike), type_(type) {}

	void setProcess(ProcessGBM p);
	void setValuationDate(Date valuationDate, Date mat);
    void setValuationDate(double t);
	void setOp(OptionProduct op);

	Date getMat() { return mat_; }
    
    virtual double bsmPrice()= 0;
    double mcPrice(unsigned long nSim);
    double fdmPrice(Solver solver, double di, double dj);
 
	double delta(double nSim);
	double delta(Solver  solver, double di, double dj);
	double gamma(double nSim);
	double gamma(Solver  solver, double di, double dj);
	double vega(double nSim);
	double vega(Solver solver, double di, double dj);
	double rho(double nSim);
	double rho(Solver solver, double di, double dj);
	double theta(double nSim);
	double theta(Solver solver, double di, double dj);
    
    virtual double BSdelta()=0;
	virtual double BSgamma()=0;
	virtual double BSvega()=0;
	virtual double BSrho()=0;
	virtual double BStheta()=0;
  
	virtual double impVol(double mktPrice, double init, double tol)=0;

	void linspace(std::vector<double> &x, double start, double end, int step);
	vector<double> thomas(vector<double> a, vector<double> b, vector<double> g, vector<double> f);
    vector<double> sor(int n, vector<double>a, vector<double>b,   double relax, int max_iter);
	virtual void printInfo() = 0;

    virtual ~optionOutput()
	{
		delete payOff_;
	}

protected:
	double getd1();
	double getd2();
    
    Date valuationDate_;
    Date mat_;

	payOff* payOff_;

	double strike_;
	OptionType type_;
    OptionProduct op_;
    
	double s_, r_, q_, t_, sig_;
    
    ProcessGBM p_;

	void setTau(double tau) { t_ = tau; }
	void setSpot(double Spot) { s_ = Spot; }
	void setSig(double sig) { sig_ = sig; }
	void setRf(double rf) { r_ = rf; }

};

class VanillaOption : public optionOutput
{
public:
    
	VanillaOption(double strike, OptionType type) :
		optionOutput(strike, type)
	{
		payOff_ = new VanillaPayOff(strike, type);
	}

	virtual double bsmPrice() override;
	virtual double BSdelta() override ;
    virtual double BSgamma() override;
    virtual double BSrho() override ;
    virtual double BStheta() override;
    virtual double BSvega() override;
	virtual double impVol(double mktPrice, double init, double tol) override;
    virtual void printInfo() override;

	~VanillaOption() {}
    
private:
	void setVol(double sigma)
	{
		sigma_ = sigma;
	}
	double sigma_;
	//OptionType type_;
    //payOff* payOff_;
};

class DigitalOption : public optionOutput
{
public:
	DigitalOption(double  strike, OptionType  type) :
		optionOutput(strike, type)
	{
		payOff_ = new DigitalPayOff(strike, type);
	}
	
	virtual double bsmPrice() override;

	virtual double BSdelta() override { return 0; };
	virtual double BSgamma() override { return 0; };
	virtual double BSvega()  override { return 0; };
	virtual double BStheta() override { return 0; };
	virtual double BSrho()   override { return 0; };
	virtual double impVol(double mktPrice, double init, double tol) override
	{ return 0; } ;

	virtual void printInfo() override;
    
	~DigitalOption() {}

private:
	//payOff* payOff_;
};



struct optionInput
{
    optionOutput* option;
    ProcessGBM gbmp;
    Date valuationDt;
    int quantity;
};
