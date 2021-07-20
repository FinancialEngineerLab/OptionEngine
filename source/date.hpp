#pragma once
#include <string>

class Date
{
public:
	Date() {}
	Date(int y, int m, int d) : y_(y), m_(m), d_(d) {}
	Date(std::string ymd);

	int year() { return y_; }
	int month() { return m_; }
	int day() { return d_; }
	int Date_Diff(Date dt);
	double Date_Between(Date dt1, Date dt2);
	void Date_Print();

	bool operator> (const Date dt);
	bool operator< (const Date dt);
	bool operator== (const Date dt);

	std::string to_str();

	~Date() {}

private:
	int y_; int m_; int d_;
};
