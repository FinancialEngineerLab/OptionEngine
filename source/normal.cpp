#define _USE_MATH_DEFINES
#include "pch.h"
#include <cmath>
#include "normal.hpp"

double stdNomPrbDist(double x, double mu, double sigma)
{
	return 1.0 / std::sqrt(2 * M_PI) / sigma * std::exp(-0.5*(x - mu)*(x - mu) / (sigma*sigma));
}

double stdNomCumDist(double x, double mu, double sigma)
{
	double v = (x - mu) / sigma;
	return 0.5 * std::erfc(-v * M_SQRT1_2);
}