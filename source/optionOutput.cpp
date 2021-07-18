//#include "stdafx.h"
#include "optionOutput.hpp"
#include "normal.hpp"

using namespace std;

void optionOutput::setValuationDate(Date dt, Date mat)
{
	Date *dateMethod = nullptr;
	valuationDate_ = dt;
	t_ = dateMethod->Date_Between(valuationDate_, mat) / 365.0;
}


void optionOutput::setValuationDate(double t)
{
	t_ = t;
}

void optionOutput::setOp(OptionProduct op)
{
	op_ = op;
}

void optionOutput::setProcess(ProcessGBM p)
{
	p_ = p;
	s_ = p_.getSpot();
	r_ = p_.getRf();
	q_ = p_.getDiv();
	sig_ = p_.getRealVol();
}

double optionOutput::getd1()
{
	return (std::log(s_ / strike_) + (r_ - q_ + 0.5 * sig_ * sig_) * t_) / (sig_*std::sqrt(t_));
}

double optionOutput::getd2()
{
	return getd1() - sig_ * std::sqrt(t_);
}

double optionOutput::mcPrice(unsigned long nSim)
{
	double sumPayoff = 0;
	double df = std::exp(-r_ * t_);

	std::mt19937_64 gen;
	std::normal_distribution<double> engine(0.0, 1.0);
	gen.seed(1);

	double eS = s_ * std::exp((r_ - q_ - 0.5*sig_*sig_)*t_);
	double diffusion = sig_ * std::sqrt(t_);

	for (unsigned int i = 0; i < nSim; ++i)
	{
		double e = engine(gen);
		for (int j = 0; j < 2; ++j)
		{
			// antithetic variates for variance reduction
			double sT = eS * std::exp(diffusion*(1 - j * 2) * e);
			double p = (*payOff_)(sT);
			sumPayoff += df * p;
		}
	}
	return sumPayoff / nSim / 2.0;
}


double optionOutput::delta(double nSim = 0)
{
	double s = s_;
	double s1 = s * 1.01;
	double s2 = s * 0.99;
	double p1 = 0.0; double p2 = 0.0;
	double delta;


	//s_ = s1;
	setSpot(s1);
	p1 = mcPrice(nSim);
	setSpot(s2);
	p2 = mcPrice(nSim);

	delta = (p1 - p2) / (0.02*s);

	setSpot(s);
	return delta;
}

double optionOutput::delta(Solver solver, double di, double dj)
{
	double s = s_;
	double s1 = s * 1.01;
	double s2 = s * 0.99;
	double p1 = 0.0; double p2 = 0.0;
	double delta;

	setSpot(s1);
	p1 = fdmPrice(solver, di, dj);
	setSpot(s2);
	p2 = fdmPrice(solver, di, dj);
	delta = (p1 - p2) / (0.02*s);

	setSpot(s);
	return delta;
}

double optionOutput::gamma(double nSim)
{
	double s = s_;
	double s1 = s * 1.01;
	double s2 = s * 0.99;
	double p0 = 0.0; double p1 = 0.0; double p2 = 0.0;

	p0 = mcPrice(nSim);
	setSpot(s1);
	p1 = mcPrice(nSim);
	setSpot(s2);
	p2 = mcPrice(nSim);
	setSpot(s);

	double gamma = (p1 + p2 - 2 * p0) / std::pow(0.01*s, 2);
	return gamma;
}


double optionOutput::gamma(Solver solver, double di, double dj)
{
	double s = s_;
	double s1 = s * 1.01;
	double s2 = s * 0.99;
	double p0 = 0.0; double p1 = 0.0; double p2 = 0.0;

	p0 = fdmPrice(solver, di, dj);
	setSpot(s1);
	p1 = fdmPrice(solver, di, dj);
	setSpot(s2);
	p2 = fdmPrice(solver, di, dj);
	setSpot(s);

	double gamma = (p1 + p2 - 2 * p0) / std::pow(0.01*s, 2);
	return gamma;
}

double optionOutput::vega(double nSim)
{
	double v = sig_;
	double v1 = v + 0.01;
	double v2 = v - 0.01;
	double p1 = 0.0; double p2 = 0.0;

	setSig(v1);
	p1 = mcPrice(nSim);
	setSig(v2);
	p2 = mcPrice(nSim);
	setSig(v);
	double vega = (p1 - p2) / 0.02;
	return vega;
}


double optionOutput::vega(Solver solver, double di, double dj)
{
	double v = sig_;
	double v1 = v + 0.01;
	double v2 = v - 0.01;
	double p1 = 0.0; double p2 = 0.0;

	setSig(v1);
	p1 = fdmPrice(solver, di, dj);
	setSig(v2);
	p2 = fdmPrice(solver, di, dj);
	setSig(v);
	double vega = (p1 - p2) / 0.02;
	return vega;

}

double optionOutput::rho(double nSim)
{
	double rf = r_;
	double r1 = rf + 0.0001;
	double r2 = rf - 0.0001;
	double p1 = 0.0; double p2 = 0.0;

	setRf(r1);
	p1 = mcPrice(nSim);
	setRf(r2);
	p2 = mcPrice(nSim);
	setRf(rf);

	double rho = (p1 - p2) / 0.0002;
	return rho;

}

double optionOutput::rho(Solver solver, double di, double dj)
{
	double rf = r_;
	double r1 = rf + 0.0001;
	double r2 = rf - 0.0001;
	double p1 = 0.0; double p2 = 0.0;

	setRf(r1);
	p1 = fdmPrice(solver, di, dj);
	setRf(r2);
	p2 = fdmPrice(solver, di, dj);
	setRf(rf);

	double rho = (p1 - p2) / 0.0002;
	return rho;
}

double optionOutput::theta(double nSim)
{
	double t = t_;
	double t1 = t - 0.0001;
	double p0 = 0.0; double p1 = 0.0;

	p0 = mcPrice(nSim);
	setTau(t1);
	p1 = mcPrice(nSim);

	setTau(t);

	double theta = (p1 - p0) / 0.0001;
	return theta;
}

double optionOutput::theta(Solver solver, double di, double dj)
{
	double t = t_;
	double t1 = t - 0.0001;
	double p0 = 0.0; double p1 = 0.0;

	p0 = fdmPrice(solver, di, dj);
	setTau(t1);
	p1 = fdmPrice(solver, di, dj);

	setTau(t);
	double theta = (p1 - p0) / (0.0001);
	return theta;
}


double VanillaOption::impVol(double mktPrice, double init = 0.8, double tol = 1.0e-6)
{
	double init_v = sig_;
	double x = init;
	double e = 9999;
	double dynamicVega;
	//VanillaOption *method = nullptr;
	//method->setProcess(p_);
	//method->setValuationDate(t_);

	while (e > tol)
	{
		setSig(x);
		double diff = bsmPrice() - mktPrice;
		e = abs(diff);
		dynamicVega = BSvega();
		x = x - diff / dynamicVega;
	}
	setVol(init_v);
	return x;
}

double VanillaOption::bsmPrice()
{
	double d1 = getd1();
	double d2 = getd2();
	double nd1 = stdNomCumDist(type_ * d1);
	double nd2 = stdNomCumDist(type_ * d2);
	double price = type_ * (s_*std::exp(-q_ * t_)*nd1 -
		strike_ * std::exp(-r_ * t_)*nd2);
	return price;
}

double VanillaOption::BSdelta()
{
	double d1 = getd1();
	double nd1 = stdNomCumDist(d1);
	double delta = std::exp(-q_ * t_)*(nd1 + (type_ - 1) / 2);
	return delta;
}

double VanillaOption::BSgamma()
{
	double d1 = getd1();
	double npd1 = stdNomCumDist(d1);
	double gamma = (npd1 * exp(-q_ * t_)) / (s_*sig_*std::sqrt(t_));
	return gamma;
}

double VanillaOption::BSrho()
{
	double d2 = getd2();
	double nd2 = stdNomCumDist(type_*d2);
	double rho = type_ * strike_*t_*std::exp(-r_ * t_)*nd2;
	return rho;
}

double VanillaOption::BStheta()
{
	double d1 = getd1();
	double d2 = getd2();
	double npd1 = stdNomPrbDist(d1);
	double nd1 = stdNomPrbDist(type_*d1);
	double nd2 = stdNomPrbDist(type_*d2);
	double theta = -s_ * npd1 * sig_ * std::exp(-q_ * t_) / (2 * std::sqrt(t_)) + type_ * q_*s_*nd1*std::exp(-q_ * t_) -
		type_ * r_* strike_ * std::exp(-r_ * t_)*nd2;
	return theta;
}

double VanillaOption::BSvega()
{
	double d1 = getd1();
	double npd1 = stdNomCumDist(d1);
	double vega = s_ * std::exp(-q_ * t_)*npd1 * std::sqrt(t_);
	return vega;
}

void VanillaOption::printInfo()
{
	std::cout << " Plain Vanilla Option Inputs" << std::endl;

	if (type_ == Call)
	{
		std::cout << " Type : Call " << std::endl;
	}
	else if (type_ == Put)
	{
		std::cout << " Type : Put " << std::endl;
	}
	std::cout << " Spot :" << s_ << std::endl;
	std::cout << " Strike :" << strike_ << std::endl;
	std::cout << " set Time to Maturty :" << t_ << std::endl;
	std::cout << " risk free flat : " << r_ << std::endl;
	std::cout << " dividend rate flat : " << q_ << std::endl;
	std::cout << " sigma : " << sig_ << std::endl;
}

double DigitalOption::bsmPrice()
{
	double d2 = getd2();
	double nd2 = stdNomCumDist(type_*d2);
	return std::exp(-r_ * t_)*nd2;
}


void DigitalOption::printInfo()
{
	std::cout << " Digital Option ($1) Inputs " << std::endl;
	if (type_ == Call)
	{
		std::cout << " Type : Call " << std::endl;
	}
	else
	{
		std::cout << " Type : Put " << std::endl;
	}
	std::cout << " Spot : " << s_ << std::endl;
	std::cout << " Strike : " << strike_ << std::endl;
	std::cout << " set Time to Maturty : " << t_ << std::endl;
	std::cout << " risk free flat : " << r_ << std::endl;
	std::cout << " dividend rate flat : " << q_ << std::endl;
	std::cout << " sigma : " << sig_ << std::endl;
}


vector<double> optionOutput::thomas(vector<double> a, vector<double> b, vector<double> g, vector<double> f)
{
	int n = f.size() - 1;
	for (int i = 2; i <= n; i++)
	{
		double mult = a[i] / b[i - 1];
		b[i] = b[i] - mult * g[i - 1];
		f[i] = f[i] - mult * f[i - 1];
	}
	vector<double> x(n + 1);
	x[n] = f[n] / b[n];
	for (int i = n - 1; i >= 1; i--)
	{
		x[i] = (f[i] - g[i] * x[i + 1]) / b[i];
	}
	return x;
}


void optionOutput::linspace(vector<double> &x, double start, double end, int step)
{
	double h = (end - start) / (step - 1.0);
	x.push_back(0.0);
	for (int i = 1; i <= step; i++)
	{
		x.push_back(start + h * (i - 1.0));
	}
}

double optionOutput::fdmPrice(Solver solver, double di, double dj)
{
	int nL = 2 * s_;
	int nt = 1 / di;
	int ns = 1 / dj;

	int w = di / (dj*dj);

	vector<double> aa(ns * 3);
	aa[0 + 0 * 3] = 0.0;
	aa[1 + 0 * 3] = 1.0;
	aa[0 + 1 * 3] = 0.0;
	for (int i = 1; i < ns - 1; i++)
	{
		aa[2 + (i - 1) * 3] = -w;
		aa[1 + i * 3] = 1.0 + 2.0*w;
		aa[0 + (i + 1) * 3] = -w;
	}
	aa[2 + (ns - 2) * 3] = 0.0;
	aa[1 + (ns - 1) * 3] = 1.0;
	aa[2 + (ns - 1) * 3] = 0.0;


	vector<double> s;
	vector<vector<double>> u(ns + 1, vector<double>(nt + 2));
	linspace(s, 0, nL, ns);

	if (op_ == Vanilla)
	{
		for (int i = 1; i <= ns; i++)
		{
			u[i][1] = (type_*(s[i] - strike_) >= 0.0) ? (type_*(s[i] - strike_)) : 0.0;
		}
	}
	else if (op_ == Digital)
	{
		for (int i = 1; i <= ns; i++)
		{
			u[i][1] = (type_*(s[i] - strike_) >= 0.0) ? 1.0 : 0.0;
		}
	}

	vector<double> d(ns), c(ns), a(ns);

	for (int i = 1; i <= ns - 1; i++)
	{
		// beta, gamma, alpha 
		d[i] = 1 / (t_ / (1 / di)) + std::pow(sig_*i, 2) + r_ - q_;
		c[i] = -(r_ - q_)* i * 0.5 - std::pow(sig_*i, 2)*0.5;
		a[i] = (r_ - q_) * i * 0.5 - std::pow(sig_*i, 2)*0.5;
	}

	// linear boundary conditions
	d[ns - 1] = d[ns - 1] + 2 * c[ns - 1];
	a[ns - 1] = a[ns - 1] - c[ns - 1];

	vector<double> b(ns);
	vector<double> x(ns + 1);
	vector<double> fvec(ns + 1);
	vector<double> temp;
	
	// Solving Algorhitm
	for (int n = 1; n <= nt; n++)
	{
		for (int i = 1; i <= ns - 1; i++)
		{
			b[i] = u[i + 1][n] / (t_ / (1 / di));
		}
		if (solver == Thomas)
		{
			temp = thomas(a, d, c, b);
		}
		//else if (solver == SOR)
		//{
		//	for (int i = 1; i <= ns - 1; i++)
		//	{
		//		b[i] = b[i] * (t_ / (1 / di));
		//	}
		//	temp = sor(ns-1, aa, b,1.2, 100);
		//}

		for (int h = 2; h <= ns; h++)
		{
			u[h][n + 1] = temp[h - 1]; // for nt, ns
			fvec[h - 1] = temp[h - 1]; // for single (ns)
		}

	}


	// above why starting index = 1, shaping 1 ~ ns
	double value = fvec[ns / 2];
	return value;


}



/*
vector<double> optionOutput::sor(int n, vector<double>a, vector<double>b,
	double relax, int max_iter)
{
	//int n = f.size() - 1;
	int iter, i;
	double square_sum = 0.0, tol = 1e-6;


	vector<double> x(n + 1);
	vector<double> x_new(n + 1);

	for (i = 0; i < n; i++)
	{
		x[i] = 0.0;
	}

	for (iter = 1; iter <= max_iter; iter++)
	{
		for (i = 1; i < n - 1; i++)
		{
			x_new[i] = b[i];
			x_new[i] = x_new[i] - a[3 * i] * x_new[i - 1];
			x_new[i] = x_new[i] - a[3 * i + 2] * x[i + 1];
			x_new[i] = x_new[i] / a[3 * i + 1];
			x_new[i] = (1.0 - relax) * x[i] + relax * x_new[i];
			square_sum += (x_new[i] - x[i]) * (x_new[i] - x[i]);
		}

		if (std::sqrt(square_sum) < tol || iter == max_iter)
		{
			x_new.clear();
			return x;
		}
		for (i = 1; i < n - 1; i++)
		{
			x[i] = x_new[i];
		}
		square_sum = 0.0;
	}
}
*/
