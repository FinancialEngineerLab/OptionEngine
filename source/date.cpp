#include "pch.h"
#include "date.hpp"
#include <ctime>
#include <iostream>
#include <sstream>

Date::Date(std::string ymd)
{
	y_ = std::atoi(ymd.substr(0, 4).c_str());
	m_ = std::atoi(ymd.substr(4, 2).c_str());
	d_ = std::atoi(ymd.substr(6, 2).c_str());
}

int Date::Date_Diff(Date dt)
{
	std::tm a = { 0,0,0,d_, m_ - 1, y_ - 1900 };
	std::tm b = { 0,0,0, dt.day(), dt.month() - 1, dt.year() - 1900 };
	std::time_t x = std::mktime(&a);
	std::time_t y = std::mktime(&b);
	double difference = std::difftime(y, x) / (24 * 60 * 60);
	return difference;
}

double Date::Date_Between(Date dt1, Date dt2)
{
	return dt2.Date_Diff(dt1);
}

bool Date::operator>(Date dt)
{
	int diff = (*this).Date_Diff(dt);
	if (diff < 0)
	{
		return true;
	}
	else
	{
		return false;
	}
}

bool Date::operator<(Date dt)
{
	int diff = (*this).Date_Diff(dt);
	if (diff > 0)
	{
		return true;
	}
	else
	{
		return false;
	}
}

bool Date::operator==(Date dt)
{
	bool a = !(*this > dt);
	bool b = !(*this < dt);

	if (a && b)
	{
		return true;
	}
	else
	{
		return false;
	}
}
std::string Date::to_str()
{
	std::stringstream ss;
	ss << y_ << "-" << m_ << "-" << d_;
	return ss.str();
}

void Date::Date_Print()
{
	std::cout << y_ << "-" << m_ << "-" << d_ << std::endl;
}