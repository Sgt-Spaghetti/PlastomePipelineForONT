#include <fstream>
#include <iostream>
#include <vector>

int calculate_AT_percentage(const std::string &read, const int &read_length){
	int ATcount = 0;
	for (int i = 0; i < read_length; i++){
		if (read[i]=='A' || read[i]=='T' || read[i]=='U'){
			++ATcount;
		}
	}
	const int read_AT_percentage = int(((double(ATcount)/double(read_length))*100)+0.5);
	return read_AT_percentage;
}

int calculate_quality(const std::string &quality, const int &read_length){
	const int ascii_code_offset_phred_score = 33;
	int Qscore = 0;
	for (int i=0; i < read_length; i++){
		Qscore += quality[i] - ascii_code_offset_phred_score;
	}
	const int Qmean = int((double(Qscore)/double(read_length))+0.5);
	return Qmean;
}

int main(int argc, char* argv[]){

	const std::string input_file_path = argv[1];
	const int min_length = atoi(argv[2]);
	const int min_quality = atoi(argv[3]);
	//const std::string output_file_path = argv[4];
	
	const int lines_per_entry = 4; 
	int point_in_entry = 0;
	const int sequence_location = 1;
	const int quality_location = 3;

	int parameters[2];
	std::string entry[4] = {"","","",""};

	bool discard_entry = false;

	//vector<int> data_output;

	std::cout << "AT_percentage, Read_Length,"<< "\n";

	std::ifstream input_file(input_file_path);

	std::string line_in_file;

	while (getline(input_file, line_in_file)){
		if (point_in_entry < lines_per_entry-1){
			entry[point_in_entry] = line_in_file;
			if (point_in_entry == sequence_location){
				int sequence_length = line_in_file.length();
				if (sequence_length < min_length) {
					discard_entry = true;
				}
			}
			++point_in_entry;
		} else {
			entry[point_in_entry] = line_in_file;
			int read_quality = calculate_quality(line_in_file, line_in_file.length());
			if (read_quality < min_quality){
				discard_entry = true;
			}
			point_in_entry = 0;
			
			if (discard_entry == false){
				int ATpercentage = calculate_AT_percentage(entry[1], entry[1].length());
				std::cout << ATpercentage << "," << entry[1].length() << "\n";
			} else {discard_entry = false;}

		}
	}

	input_file.close();

	return 0;
}
