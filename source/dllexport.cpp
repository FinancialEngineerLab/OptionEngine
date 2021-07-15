#include "pch.h"
#include "optionFactory.hpp"
#include "optionOutput.hpp"

extern double __stdcall BSPrice(OptionProduct op, OptionType type, double spot, double strike,
	double rf, double div, double vol, double tau);
extern double __stdcall MCPrice(OptionProduct op, OptionType type, double spot, double strike,
	double rf, double div, double vol, double tau, unsigned long nSim);
extern double __stdcall FDMPrice(OptionProduct op, OptionType type, double spot, double strike,
	double rf, double div, double vol, double tau, Solver solver, double di, double dj);

extern double __stdcall Delta(OptionProduct op, OptionType type, double spot, double strike,
	double rf, double div, double vol, double tau, string method,
	double nSim, Solver solver, double di, double dj);
extern double __stdcall BSDelta(OptionProduct op, OptionType type, double spot, double strike,
	double rf, double div, double vol, double tau);



double __stdcall MCPrice(OptionProduct op, OptionType type, double spot, double strike,
	double rf, double div, double vol, double tau, unsigned long nSim)
{
	double price;
	OptionFactory factory;
	ProcessGBM process(spot, rf, div, vol);
	optionOutput* option = factory.CreateOption(op, type, strike);

	option->setValuationDate(tau);
	option->setProcess(process);
	price = option->mcPrice(nSim);
	return price;
}


double __stdcall BSPrice(OptionProduct op, OptionType type, double spot, double strike,
		double rf, double div, double vol, double tau)
{
		double price;
		OptionFactory factory;
		ProcessGBM process(spot, rf, div, vol);
		optionOutput* option = factory.CreateOption(op, type, strike);

		option->setValuationDate(tau);
		option->setProcess(process);
		price = option->bsmPrice();
		return price;
}


double __stdcall FDMPrice(OptionProduct op, OptionType type, double spot, double strike,
		double rf, double div, double vol, double tau, Solver solver, double di, double dj)
{
		double price;
		OptionFactory factory;
		ProcessGBM process(spot, rf, div, vol);
		optionOutput* option = factory.CreateOption(op, type, strike);

		option->setValuationDate(tau);
		option->setProcess(process);
		price = option->fdmPrice(solver,di, dj);
		return price;
}


double __stdcall Delta(OptionProduct op, OptionType type, double spot, double strike,
	double rf, double div, double vol, double tau, string method,
	double nSim=0, Solver solver=Thomas, double di=1, double dj=1)
{
	double value;
	OptionFactory factory;
	ProcessGBM process(spot, rf, div, vol);
	optionOutput* option = factory.CreateOption(op, type, strike);

	option->setValuationDate(tau);
	option->setProcess(process);
	if (method == "MC")
	{
		value = option->delta(nSim);
	}
	else if (method == "FDM")
	{
		value = option->delta(solver, di, dj);
	}
	return value;
}


double __stdcall BSDelta(OptionProduct op, OptionType type, double spot, double strike,
		double rf, double div, double vol, double tau)
{
	double value;
	OptionFactory factory;
	ProcessGBM process(spot, rf, div, vol);
	optionOutput* option = factory.CreateOption(op, type, strike);

	option->setValuationDate(tau);
	option->setProcess(process);
	value = option->BSdelta();
	return value;
}