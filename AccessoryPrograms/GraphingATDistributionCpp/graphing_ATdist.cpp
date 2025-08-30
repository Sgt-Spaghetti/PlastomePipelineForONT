#include <iostream>
#include <fstream>
#include <vector>
#include <typeinfo>

using namespace std;

int main (int argc, char* argv[]) {

	string input_file_path = argv[1];
	ifstream input_file(input_file_path);
	
	int number_of_bins = 100;
	int minimum = 0;
	int maximum = 100;
	double stepsize = double(maximum - minimum) / double(number_of_bins);
	double bins[number_of_bins]{ };
	long int frequency_density [number_of_bins] { };
	for (int i = 0; i <  number_of_bins; i++){
		bins[i] = i * stepsize;
	}

	string line_contents;
	int line_number = 0;
	int section = 0;
	
	while (getline(input_file, line_contents)){
		string value_percent = "0";
		string value_length = "0";
		if (line_number > 0){
			for (int i=0; i < line_contents.length(); i++){
				char digit = line_contents[i];
				if (digit == ','){
					++section;
				}
				else {
					if (section == 1) {value_percent += digit;}
					else if (section == 2) {value_length += digit;}
				}
			}
		}
		section = 0;
		line_number += 1;
		
		double ATpercentage = stod(string(value_percent));		
		int read_length = stoi(value_length);
		bool binned = false;
		int i = 0;
		while (binned == false){
			if (ATpercentage > bins[i]) {
				++i;
			} else {
				frequency_density[i] = frequency_density[i] + read_length;
				binned = true;
			} 
		}
		
		value_percent = "0";
		value_length = "0";
	}
	for (int i=0; i < number_of_bins; i++){
		cout << bins[i] << ',' << frequency_density[i] << "\n";
	}
	return 0;
}
