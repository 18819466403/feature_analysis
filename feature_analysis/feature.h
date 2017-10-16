#pragma once
#include"matrix.h"
#include"iostream"
using namespace std;
class Feature{
private:
	Matrix<int>*block;
	// ha
	double har_mean;
	double har_var;

	double cv;

	double hist;
	double hist_mean;
	double hist_variance;
	double hist_skewness;
	double hist_kurtosis;

	double bdip;
	double nmsid;

	double glcm_entropy;
	double glcm_contrast;
	double glcm_uniformity;
	double glcm_idm;
	double glcm_im;
	double glcm_correlation;
	double glcm_energy;
	double glcm_brightness;

	double tamura_coarseness;
	double tamura_contrast;
	double tamura_directionality;
	double tamura_linelikeness;
	double tamura_regularity;
	double tamura_roughness;

protected:
	void cal_har_mean_var();
	void cal_cv();
	void cal_hist();
	void cal_bdip();
	void cal_nmsid();
	void cal_glcm();
public:
	Feature(Matrix<int> &matrix);
	Matrix<int> & getBlock();
	
	double getHar_mean() {
		return har_mean;
	}
	double getHar_var() {
		return har_var;
	}
	double getCv() {
		return cv;
	}
	double getHist(){
		return hist_skewness;
	}
	double getHist_mean() {
		return hist_mean;
	}
	double getHist_variance(){
		return hist_variance;
	}
	double getHist_skewness() {
		return hist_skewness;
	}
	double getHist_kurtosis() {
		return hist_kurtosis;
	}
	double getBdip() {
		return bdip;
	}
	double getNmsid() {
		return nmsid;
	}
	double operator[](int i) {
		switch (i)
		{
		case 0:  return har_mean; break;
		case 1:  return har_var; break;
		case 2:  return cv; break;
		case 3:  return hist; break;
		case 4:  return hist_mean; break;
		case 5:  return hist_variance; break;
		case 6:  return hist_skewness; break;
		case 7:  return hist_kurtosis; break;
		case 8:  return bdip; break;
		case 9:  return nmsid; break;
		case 10:  return glcm_entropy; break;
		case 11:  return glcm_contrast; break;
		case 12:  return glcm_uniformity; break;
		case 13:  return glcm_idm; break;
		case 14:  return glcm_im; break;
		case 15:  return glcm_correlation; break;
		case 16:  return glcm_energy; break;
		case 17:  return glcm_brightness; break;
		case 18:  return tamura_coarseness; break;
		case 19:  return tamura_contrast; break;
		case 20:  return tamura_directionality; break;
		case 21:  return tamura_linelikeness; break;
		case 22:  return tamura_regularity; break;
		case 23:  return tamura_roughness; break;
		default: cerr << i << " is out of range\n";
		}

	}
	
};


Feature::Feature(Matrix<int> &matrix) {
	if (matrix.getRow() != matrix.getColumn())
		return;
	block = new Matrix<int>(matrix.getRow(), matrix.getColumn());
	*block = matrix;

	cal_har_mean_var();
	cal_cv();
	cal_hist();
	cal_bdip();
	cal_nmsid();
}

Matrix<int> & Feature::getBlock() {
	return *block;
}

void Feature::cal_glcm() {

}

void Feature::cal_har_mean_var() {
	Matrix<double> temp(block->getRow()/2, block->getColumn()/2);
	double total = 0;
	for (int i = 0; i < temp.getRow(); i++)
	{
		for (int j = 0; j < temp.getColumn(); j++)
		{
			temp[i][j] = (double)(block->getElement(2 * i, 2 * j) + block->getElement(2 * i + 1, 2 * j) +
				block->getElement(2 * i, 2 * j + 1) + block->getElement(2 * i + 1, 2 * j + 1)) / 4;
			total += temp[i][j];
		}
	}
	har_mean = total / (block->getRow()*block->getColumn());

	double var = 0;
	for(int i=0;i<temp.getRow();i++)
		for (int j = 0; j < temp.getColumn(); j++) {
			var += (temp[i][j] - har_mean)*(temp[i][j] - har_mean);
		}

	har_var = var / (block->getRow()*block->getColumn());	
}

void Feature::cal_cv() {
	double total = 0;
	for(int i=0;i<block->getRow();i++)
		for (int j = 0; j < block->getColumn(); j++) {
			total += block->getElement(i, j);
		}
	double mean = (double)(total) / (block->getRow()*block->getColumn());
	double var_total = 0;
	for (int i = 0; i<block->getRow(); i++)
		for (int j = 0; j < block->getColumn(); j++) {
			var_total += (block->getElement(i, j) - mean)*(block->getElement(i, j) - mean);
		}
	double standard = sqrt(var_total / (block->getRow()*block->getColumn()));
	if (mean == 0) cv = 0;
	else cv = standard / mean;
}

void Feature::cal_hist() {
	hist_mean = 0.0;
	hist_variance=0.0;
	hist_skewness=0.0;
	hist_kurtosis=0.0;
	double hist_standard = 0;
	int hito[256];
	double p[256];
	hist = 0;
	for (int i = 0; i < 256; i++)
	{
		hito[i] = 0;
	}
	for(int i=0;i<block->getRow();i++)
		for (int j = 0; j < block->getColumn(); j++) {
			hito[block->getElement(i, j)]++;
		}
	int max = 0;
	for (int i = 0; i < 256; i++)
	{
		if (hito[i] > max) max = i;
		p[i] =1.0* hito[i] / block->getRow() / block->getColumn();
		hist_mean += i*p[i];
	}
	for (int i = max - 10; i < max + 10; i++) {
		if (i >= 0 && i <= 255 && i!=max) {
			hist += hito[i];
		}
	}

	for (int i = 0; i < 256; i++)
	{
		hist_variance += (i - hist_mean)*(i - hist_mean)*p[i];
	}
	hist_standard = sqrt(hist_variance);
	if (hist_mean == 0) {
		hist_skewness = 0;
		hist_kurtosis = 0;
	}
	else
	{
		for (int i = 0; i < 256; i++) {
			if (hist_standard != 0)
			{
				hist_skewness += pow(hist_standard, -3)*(i - hist_mean)*(i - hist_mean)*(i - hist_mean)*p[i];
				hist_kurtosis += pow(hist_standard, -4)*(i - hist_mean)*(i - hist_mean)*(i - hist_mean)*(i - hist_mean)*p[i];
			}
		}
	}
}

void Feature::cal_bdip() {
	int max = 0;
	int total_pixel = 0;
	for(int i=0;i<block->getRow();i++)
		for (int j = 0; j < block->getColumn(); j++) {
			if (block->getElement(i, j) > max) max = block->getElement(i, j);
			total_pixel += block->getElement(i, j);
	
		}
	if (max == 0) bdip = 0;
	else bdip =1.0*(block->getRow()*block->getColumn() - total_pixel) / max;
}

void Feature::cal_nmsid()
{
	nmsid = 0;
	int M = block->getRow();
	for (int k = 1; k <= M; k++) {
		double a = 0, b = 0, c = 0, d = 0;
		for (int x = 0; x < M; x++)
		{
			for (int y = 0; y < M-k; y++)
			{
				a += (double)(abs(block->getElement(x, y) - block->getElement(x, y + k)))/M/(M-k);
			}
		}
		for(int x=0;x<M-k;x++)
			for (int y = 0; y < M; y++) {
				b += (double)(abs(block->getElement(x, y) - block->getElement(x + k, y)))/ M / (M - k);
			}

		for(int x=0;x<M-k;x++)
			for (int y = 0; y < M - k; y++) {
				c += (double)(abs(block->getElement(x, y) - block->getElement(x + k, y + k)))/ (M - k) / (M - k);
			}

		for(int x=0;x<M-k;x++)
			for (int y = 0; y < M - k; y++) {
				d += (double)(abs(block->getElement(x, M - y -1) - block->getElement(x, M - y - k -1)))/ (M - k) / (M - k);
			}

		nmsid += (a + b + c + d) / 4;
		
	}

}

