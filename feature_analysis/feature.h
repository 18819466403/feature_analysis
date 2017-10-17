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

	double glcm_energy;
	double glcm_contrast;
	double glcm_uniformity;
	double glcm_idm;
	double glcm_im;
	double glcm_correlation;
	double glcm_entropy;
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
	void cal_glcm(int grade, int gap);
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
		case 10:  return glcm_energy; break;
		case 11:  return glcm_contrast; break;
		case 12:  return glcm_uniformity; break;
		case 13:  return glcm_idm; break;
		case 14:  return glcm_im; break;
		case 15:  return glcm_correlation; break;
		case 16:  return glcm_entropy; break;
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
	cal_glcm(16,1);
}

Matrix<int> & Feature::getBlock() {
	return *block;
}

void Feature::cal_glcm(int grade, int gap) {
	Matrix<size_t> mat_0(grade, grade);
	Matrix<size_t> mat_45(grade, grade);
	Matrix<size_t> mat_135(grade, grade);
	Matrix<size_t> mat_90(grade, grade);//GLCM of four directions of 0, 45, 90, 135 degree. 
	Matrix<double> mat(grade, grade);
	size_t grade_factor = 256 / grade;
	size_t row = block->getRow();
	size_t cloumn = block->getColumn();
	for (int i = 0; i < block->getRow(); i++)
	{
		for (int j = 0; j < block->getColumn()-gap; j++) {
			mat_0[block->getElement(i, j) / grade_factor][block->getElement(i, j + gap) / grade_factor]++;
			if(block->getElement(i, j) / grade_factor!= block->getElement(i, j + gap) / grade_factor)
			mat_0[block->getElement(i, j+gap) / grade_factor][block->getElement(i, j ) / grade_factor]++;
		}
	}

	for (int j = 0; j < block->getColumn(); j++)
		for (int i = 0; i < block->getRow() - gap; i++)
		{
			mat_90[block->getElement(i+gap, j) / grade_factor][block->getElement(i, j) / grade_factor]++;
			if (block->getElement(i+gap, j) / grade_factor != block->getElement(i, j) / grade_factor)
			mat_90[block->getElement(i, j) / grade_factor][block->getElement(i+gap, j) / grade_factor]++;
		}

	for (int i = 0; i < block->getRow() - gap; i++)
	{
		for (int j = 0; j < block->getColumn() - gap; j++) {
			mat_45[block->getElement(i, j) / grade_factor][block->getElement(i+gap, j+gap) / grade_factor]++;
			if (block->getElement(i + gap, j+gap) / grade_factor != block->getElement(i, j) / grade_factor)
				mat_45[block->getElement(i+gap, j+gap) / grade_factor][block->getElement(i, j) / grade_factor]++;
		}
	}

	for (int i = 0; i < block->getRow() - gap; i++){
		for (int j = block->getColumn()-1; j >= gap; j--) {
			mat_135[block->getElement(i, j) / grade_factor][block->getElement(i + 1, j - 1) / grade_factor]++;
			if(block->getElement(i, j) / grade_factor!= block->getElement(i + 1, j - 1) / grade_factor)
				mat_135[block->getElement(i+1, j-1) / grade_factor][block->getElement(i, j) / grade_factor]++;
		}
	}
    
	for (int i = 0; i < grade; i++) {
		for (int j = 0; j < grade; j++) {
			mat[i][j] = (mat_0[i][j] + mat_45[i][j] + mat_90[i][j] + mat_135[i][j])*1.0 / 4;
		}
	}

	// calcute energy
	glcm_energy = 0;
	for (int i = 0; i < grade; i++) {
		for (int j = 0; j < grade; j++) {
			glcm_energy += mat[i][j] * mat[i][j];
		}
	}

	//calcute contrast
	glcm_contrast = 0;
	for (int i = 0; i < grade; i++) {
		for (int j = 0; j < grade; j++) {
			glcm_contrast += (i - j)*(i - j)*mat[i][j];
		}
	}

	//calcute uniformity
	glcm_uniformity = 0;
	for (int i = 0; i < grade; i++) {
		for (int j = 0; j < grade; j++) {
			glcm_uniformity += mat[i][j]*1.0 / (1 + abs(i - j));
		}
	}

	//calcute inverse difference moment
	glcm_idm = 0;
	for (int i = 0; i < grade; i++) {
		for (int j = 0; j < grade; j++) {
			glcm_idm += mat[i][j]*1.0 / (1 + (i - j)*(i - j));
		}
	}

	//calcute inertia moment
	glcm_im = 0;
	for (int i = 0; i < grade; i++) {
		for (int j = 0; j < grade; j++) {
			glcm_im += (i - j)*(i - j)*mat[i][j];
       }
	}

	//calcute correlation
	glcm_correlation = 0;
	int miu_x = 0;
	int miu_y = 0;
	int sei_ta_x = 0;
	int sei_ta_y = 0;
	for (int i = 0; i < grade; i++) {
		int tmp_x = 0;
		for (int j = 0; j < grade; j++) {
			tmp_x += mat[i][j];
		}
		miu_x += i*tmp_x;
	}
	for (int j = 0; j < grade; j++) {
		int tmp_y = 0;
		for (int i = 0; i < grade; i++) {
			tmp_y += mat[i][j];
		}
		miu_x += j*tmp_y;
	}
	for (int i = 0; i < grade; i++) {
		int tmp = 0;
		for (int j = 0; j < grade; j++) {
			tmp += mat[i][j];
		}
		sei_ta_x += (i - miu_x)*(i - miu_x)*tmp;
	}
	for (int j = 0; j < grade; j++) {
		int tmp = 0;
		for (int i = 0; i < grade; i++) {
			tmp += mat[i][j];
		}
		sei_ta_y += (j - miu_y)*(j - miu_y)*tmp;
	}
	for (int i = 0; i < grade; i++) {
		for (int j = 0; j < grade; j++) {
			if(sei_ta_x != 0 && sei_ta_y != 0)
			glcm_correlation += (i - miu_x)*(j - miu_y)*mat[i][j] * 1.0 / sei_ta_x / sei_ta_y;
		}
	}
	//calcute entropy
	glcm_entropy = 0;
	for (int i = 0; i < grade; i++) {
		for (int j = 0; j < grade; j++) {
			if(0 != mat[i][j])
			glcm_entropy -= mat[i][j] * log10(mat[i][j]);
		}
	}

	//calcute brightness
	glcm_brightness = 0;
	for (int i = 0; i < grade; i++) {
		for (int j = 0; j < grade; j++) {
			glcm_brightness += (i + j - miu_x - miu_y)*(i + j - miu_x - miu_y)*(i + j - miu_x - miu_y)*mat[i][j];
		}
	}
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

