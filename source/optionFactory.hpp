#pragma once
#include "optionOutput.hpp"

class OptionFactory
{
public:
	OptionFactory() {}
	
	optionOutput* CreateOption(OptionProduct op, OptionType type,
		double strike)
	{
		optionOutput* option = NULL;
		if (op == Vanilla)
		{
			option = new VanillaOption(strike, type);
		}
		else if (op == Digital)
		{
			option = new DigitalOption(strike, type);
		}
		option->setOp(op);
		return option;
	}
	~OptionFactory() {}

protected:
	OptionProduct op_;
	OptionType type_;

};