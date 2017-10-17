#include"feature.h"
#include"image.h"
#include"matrix.h"
#include"imageProcess.h"
#include<iostream>
#include<fstream>
#define OUT 0
#define OUTEND 1
using namespace std;
typedef double type;


const int sampleNumber = 408;//1593; //样本数量
const int data_length = 18;//6;   //每个样本向量长度（特征数量）
					  //const int k = 2 * data_length + 1; //隐藏层节点数目
const int output = 1; //输出层节点数目

Matrix<type> samples(sampleNumber, data_length); //输入样本
Matrix<type> expectOutput(sampleNumber, output); //期望输出
Matrix<type> actualOutput(sampleNumber, output); //输入样本的世界输出 

double feature_means[data_length];
double feature_standards[data_length];
void extractBlock(cv::String imageName, cv::String edgeImageName, int non_thy, int thy) {
	Image edgeimage(edgeImageName);
	Image sourceImage(imageName);
	edgeimage.loadImage(0);
	sourceImage.loadImage(0);
	Matrix<int> sourceMatrix(sourceImage);
	Matrix<int> matrix(edgeimage);
	Matrix<int> flag(matrix.getRow(), matrix.getColumn());
	Matrix<int> filter(matrix.getRow(), matrix.getColumn());
	int *region = NULL;
	locate(sourceMatrix, region);
	cout << non_thy << " " << thy << endl;
	AWMF(sourceMatrix, filter, 9, 0.25, 10);
	morphology(filter);
	compensation(filter);
	region[0] = 0;
	region[1] = matrix.getRow();
	for (int i = 0; i < matrix.getRow(); i++) {
		for (int j = 0; j < matrix.getColumn(); j++) {
			if (matrix[i][j] >= 200) flag[i][j] = 1;
			else flag[i][j] = 0;
		}

	}
	for (int i = 3 + region[0]; i < region[1] - 3;) {
		for (int j = 3; j < matrix.getColumn() - 3;) {
			if (i + 15 <= region[1] && j + 15 <= matrix.getColumn()) {
				int flags = 0;
				for (int m = i; m <= i + 15; m++) {
					for (int n = j; n <= j + 15; n++) {
						flags += flag[m][n];
					}
				}


				if (flags == 0) {
					Matrix<int> non_thyroid_region_block(16, 16);

					for (int a = 0; a < 16; a++)
						for (int b = 0; b < 16; b++)
							non_thyroid_region_block[a][b] = filter[i + a][j + b];

					Feature fea(non_thyroid_region_block);
					if (fea.getHar_mean() != 0)
					{
						Image non_thyroid(non_thyroid_region_block.getArray(), 16, 16);
						cv::String road = "K:\\trainImage\\block\\non_thyroid\\";

						road += to_string(non_thy++);
						road += ".bmp";

						non_thyroid.saveImage(road);
					}
				}
				if (flags == 256) {
					Matrix<int> thyroid_region_block(16, 16);
					for (int a = 0; a < 16; a++)
						for (int b = 0; b < 16; b++)
							thyroid_region_block[a][b] = filter[i + a][j + b];
					Image thyroid(thyroid_region_block.getArray(), 16, 16);
					cv::String road = "K:\\trainImage\\block\\thyroid\\";

					road += to_string(thy++);
					road += ".bmp";

					thyroid.saveImage(road);
				}


			}

			j += 16;
		}
		i += 16;
	}
	cout << --non_thy << " " << --thy << " " << non_thy + thy << endl;
	delete[]region;
	region = NULL;
}
void generateFeature(int non_thyroid_number, int thyroid_number, int way) {
	ofstream outstuf;
	if (way == 1)
		outstuf.open("K:\\trainImage\\feature\\features.txt", ios::out && ios::end);
	else
		outstuf.open("K:\\trainImage\\feature\\features.txt", ios::out);
	for (int i = 0; i < non_thyroid_number; i++) {
		cv::String fileName = "K:\\trainImage\\block\\non_thyroid\\";
		fileName += to_string(i + 1);
		fileName += ".bmp";
		Image image(fileName);
		image.loadImage(0);
		Matrix<int> matrix(image);
		Feature feature(matrix);
		for (int i = 0; i < data_length; i++)
			outstuf << feature[i] << ' ';
		outstuf <<0<< '\n';
	}
	for (int i = 0; i < thyroid_number; i++) {
		cv::String fileName = "K:\\trainImage\\block\\thyroid\\";
		fileName += to_string(i + 1);
		fileName += ".bmp";
		Image image(fileName);
		image.loadImage(0);
		Matrix<int> matrix(image);
		Feature feature(matrix);
		
		for (int i = 0; i < data_length; i++)
			outstuf << feature[i] << ' ';
		outstuf <<1<< '\n';
	}


	outstuf.close();
}
void generateSamples()
{
	double har_mean, har_var, cv, hist, hist_mean, hist_variance, hist_skewness, hist_kurtosis, bdip, nmsid, kind;
	ifstream in("K:\\trainImage\\feature\\features.txt", ios::in);
	in.seekg(0, ios::beg);
	if (!in) {
		cerr << "File could not  be open." << endl;
		abort();
	}
	double total[data_length];

	for (int i = 0; i < data_length; i++)
	{
		total[i] = 0;
		feature_means[i] = 0;
		feature_standards[i] = 0;
	}

	for (int i = 0; i < sampleNumber; i++) {
		for (int j = 0; j < data_length; j++) {
			in >> samples[i][j];
			total[j] += samples[i][j];
		}
		in >> expectOutput[i][0];
	}
	for (int i = 0; i < data_length; i++) {
		feature_means[i] = total[i] / sampleNumber;
		for (int j = 0; j < sampleNumber; j++) {
			feature_standards[i] += (feature_means[i] - samples[j][i])*(feature_means[i] - samples[j][i]);
		}
		if (feature_means[i] != 0)
			feature_standards[i] = sqrt(feature_standards[i] / sampleNumber);
		else
			feature_standards[i] = 0;
	}
	for (int i = 0; i < sampleNumber; i++) {
		for (int j = 0; j < data_length; j++) {
			if (feature_standards[j] != 0)
				samples[i][j] = (samples[i][j] - feature_means[j]) / feature_standards[j];
		}
	}
	in.close();
}

void save_feature() {
	ofstream output;
	output.open("K:\\trainImage\\roc_features\\features.txt", ios::out);
	for (int j = 0; j < data_length; j++) {
		for (int i = 0; i < sampleNumber; i++)
		{
			output << samples[i][j] << ' ';
		}
		output << '\n';
	}
	
	for (int i = 0; i < sampleNumber; i++)
		output << expectOutput[i][0] << ' ';
	output.close();
}

int main() {
	generateFeature(378, 30, OUT);
	generateSamples();
	save_feature();
	//cv::String fileName = "K:\\trainImage\\block\\non_thyroid\\";
	//fileName += to_string(1);
	//fileName += ".bmp";
	//Image image(fileName);
	//image.loadImage(0);
	//Matrix<int> matrix(image);
	//Feature feature(matrix);
	//for (int i = 0; i < data_length; i++)
	//	cout << feature[i] << ' ';
	//cout << 0 << '\n';

}
