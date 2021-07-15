#pragma once

class ProcessGBM
{
public:
	ProcessGBM() {}
	ProcessGBM(double spot, double rf, double div, double realVol) :
		spot_(spot), rf_(rf), div_(div), realVol_(realVol) {}

	double getSpot() { return spot_; }
	double getRf() { return rf_; }
	double getDiv() { return div_; }
	double getRealVol() { return realVol_; }

	void setSpot(double spot) { spot_ = spot; }
	void setRf(double rf) { rf_ = rf; }
	void setDiv(double div) { div_ = div; }
	void setSig(double vol) { realVol_ = vol; }

	~ProcessGBM() {}

protected:
	double spot_, rf_, div_, realVol_;
};