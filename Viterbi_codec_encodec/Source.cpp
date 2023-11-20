#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <cstdlib>
#include <ctime>
#include "gnuplot-iostream.h"

using namespace std;

int one(int x) {
	int count = 0;
	for (; x > 0; x = x >> 1) {
		if (x % 2 == 1) {
			count++;
		}
	}
	return count;
}
string canal_p_0_1(string in, double p) {
	string out;
	std::srand(std::time(nullptr));
	for (char c : in) {
		double r = (double)rand() / RAND_MAX;
		if (p >= r) {
			if (c == '1') {
				out.push_back('0');
			}
			else
				out.push_back('1');
		}
		else
			out.push_back(c);
	}
	return out;
}

class Codec_Viterbi
{
public:
	Codec_Viterbi();
	Codec_Viterbi(vector<int> polynoms);
	~Codec_Viterbi();

	string codec(string in);

private:
	vector<int> polynom;
	int kl;
};
Codec_Viterbi::Codec_Viterbi() {
	polynom = { 13,15 };
	kl = 0;
	for (int i = 0; i < polynom.size(); ++i) {
		if ((int)log(polynom[i]) > kl) {
			kl = (int)log(polynom[i]);
		}
	}
	kl++;
}
Codec_Viterbi::Codec_Viterbi(vector<int> pols) {
	polynom = pols;
	kl = 0;
	for (int i = 0; i < polynom.size(); ++i) {
		if ((int)log(polynom[i]) > kl) {
			kl = (int)log(polynom[i]);
		}
	}
	kl++;
}
Codec_Viterbi::~Codec_Viterbi() {
}
string Codec_Viterbi::codec(string in) {
	string res;
	int sc = 0;
	for (char c : in) {
		if (c == '1') {
			sc = sc >> 1;
			sc += (int)pow(2, kl);
		}
		else
			sc = sc >> 1;
		for (int i = 0; i < polynom.size(); ++i) {
			char b = polynom[i] & sc;
			if (one(b) % 2 == 1)
				res.push_back('1');
			else
				res.push_back('0');
		}
	}
	return res;
}

class Encodec_Viterbi
{
public:
	Encodec_Viterbi();
	~Encodec_Viterbi();
	Encodec_Viterbi(vector<int> pols);

	string encodec(string in);
	void filling_matrix_weight();
private:
	vector<int> polynom;
	vector<vector<int>> matrix_weight;
	int kl;
	int max_l_res;
};
Encodec_Viterbi::Encodec_Viterbi() {
	polynom = { 13,15 };
	kl = 0;
	for (int i = 0; i < polynom.size(); ++i) {
		if ((int)log(polynom[i]) > kl) {
			kl = (int)log(polynom[i]);
		}
	}
	kl++;
	max_l_res = pow(2, kl + 2);
}
Encodec_Viterbi::Encodec_Viterbi(vector<int> pols) {
	polynom = pols;
	kl = 0;
	for (int i = 0; i < polynom.size(); ++i) {
		if ((int)log(polynom[i]) > kl) {
			kl = (int)log(polynom[i]);
		}
	}
	kl++;
	max_l_res = pow(2, kl + 2);  // огриначение на кол-во путей в памяти
}
Encodec_Viterbi::~Encodec_Viterbi() {
}
void Encodec_Viterbi::filling_matrix_weight() {

	for (int i = 0; i < pow(2, kl); ++i) {
		vector<int> tmp;
		for (int j = 0; j < pow(2, kl); ++j) {
			tmp.push_back(-1);
		}
		matrix_weight.push_back(tmp);
	}

	for (int i = 0; i < pow(2, kl); ++i) { // перебор всех вариантов регистра
		for (int j = 0; j < 2; ++j) { // перебор вариантов пришедших битов за шаг 
			int u = i;
			u += j << kl;
			matrix_weight[i][u >> 1] = 0;
			for (int x = 0; x < polynom.size(); ++x) { //получение выходных битов
				char b = polynom[x] & u; // выход регистра 
				int count = 0;
				matrix_weight[i][u >> 1] = matrix_weight[i][u >> 1] << 1;
				for (; b > 0; b = b >> 1) { // подсчет 1чек
					if (b % 2 == 1)
						count++;
				}
				if (count % 2 == 1)  // бит выхода
					matrix_weight[i][u >> 1] += 1;
			}
		}
	}
}
string Encodec_Viterbi::encodec(string in) {
	filling_matrix_weight();
	vector<int> input;
	int buff = 0, count = 0;
	for (char c : in) {
		if (c == '1') {
			buff = buff << 1;
			buff += 1;
			count++;
		}
		else {
			buff = buff << 1;
			count++;
		}
		if (count == polynom.size()) {
			input.push_back(buff);
			count = 0;
			buff = 0;
		}
	}

	vector<vector<int>> paths = { {0} };
	vector<int> weight_paths = { 0 };
	int min = INT_MAX, ind = 0;

	// в случае если не все сообщение декодированно
	for (int n_byts = 0; n_byts < input.size(); ++n_byts) {
		vector<int> weight_paths_sort = weight_paths;
		sort(weight_paths_sort.begin(), weight_paths_sort.end());
		// проверка на большое кол-во путей в памяти -> удаление путей с большим весом
		for (; paths.size() > max_l_res;) {
			for (int i = 0; i < weight_paths.size() and paths.size() > max_l_res; ++i) {
				if (weight_paths[i] >= weight_paths_sort[weight_paths_sort.size() / 2]) {
					paths.erase(paths.begin() + i);
					weight_paths.erase(weight_paths.begin() + i);
					i--;
				}
			}
		}
		// добавление новых путей 
		int len_res = paths.size();
		for (int i = 0; i < len_res; ++i) {
			vector<int> tmp_p = paths[i];
			for (int j = 0; j < matrix_weight.size(); ++j) {
				int tmp_w = weight_paths[i];
				int pos = paths[i][paths[i].size() - 1];
				if (matrix_weight[pos][j] >= 0) {
					tmp_w += one(matrix_weight[pos][j] ^ input[n_byts]);
					tmp_p.push_back(j);
					paths.push_back(tmp_p);
					weight_paths.push_back(tmp_w);
					tmp_p.pop_back();
				}
			}
		}
		paths.erase(paths.begin(), paths.begin() + len_res);
		weight_paths.erase(weight_paths.begin(), weight_paths.begin() + len_res);
	}
	//поиск пути с минимальным весом
	for (int i = 0; i < paths.size(); ++i) {
		for (int j = 1; j < paths[i].size(); ++j) {
			weight_paths[i] += one(matrix_weight[paths[i][j - 1]][paths[i][j]] ^ input[j - 1]);
		}
		if (min > weight_paths[i]) {
			min = weight_paths[i];
			ind = i;
		}
	}

	string res;
	for (double x = 1; x < paths[ind].size(); ++x) {
		int r = paths[ind][x];
		r = r >> (kl - 1);
		(r % 2) ? res.push_back('1') : res.push_back('0');
	}
	return res;
}

int main() {
	vector<int> polynoms = { 13,15,11 }; // 1101 1111 1011
	//vector<int> polynoms = { 13,15 }; //1101 1111
	vector<double> probability_error_s;
	vector<double> probability;
	for (double p = 0; p <= 1; p += 0.01) {
		int count_error = 0, count_right = 0;
		for (int i = 0; i < 100; i++) {
			string in = "000000000000000000000000000000";
			//cout << "in - " << in << '\n';
			Codec_Viterbi test_codec(polynoms);
			string incanal = test_codec.codec(in);
			//cout << "codec - " << incanal << '\n';
			string in_plus_error = canal_p_0_1(incanal, p);
			//cout << in_plus_error << '\n';
			Encodec_Viterbi test_encodec(polynoms);
			string out = test_encodec.encodec(in_plus_error);
			//cout << "out - " << out << '\n' << '\n';
			for (int j = 0; j < out.size(); ++j) {
				if (in[j] != out[j])
					count_error++;
				else
					count_right++;
			}
		}
		probability_error_s.push_back((double)count_error / (count_error + count_right));
		probability.push_back(p);
	}
	vector<vector<double>> res;
	res.push_back(probability);
	res.push_back(probability_error_s);
	Gnuplot gp;
	gp << "plot '-' with lines title 'v0' ,"
		<< "'-' with lines title 'v' \n";
	gp.send(res);
	return 0;
}