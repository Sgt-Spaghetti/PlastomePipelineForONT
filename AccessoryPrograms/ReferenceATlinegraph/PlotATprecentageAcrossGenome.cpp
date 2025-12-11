/* Calculate the AT percentage of "x" bp bins moving along
   a given sequence (in fasta format). Defaults to 100bp
   bins. Outputs results as a CSV directly to stdout, so
   can easily be piped to other programs or saved into a
   file.
*/

#include <iostream>
#include <string>
#include <fstream>

double calculate_AT_percentage(std::string& sequence){
	int ATcount = 0;
	int sequence_length = sequence.length();
	for (int i = 0; i < sequence_length; i++){
		char base = sequence[i];
		// convert lowercase to uppercase
		if (base >= 'a' && base <= 'z'){
			base = base - 32;
		} 
		if (base == 'A' || base == 'T' || base == 'U'){
			++ATcount;
		}
	}
	double ATpercentage = (double(ATcount) / double(sequence_length))*100;
	return ATpercentage;
}

int main(int argc, char* argv[]){

	std::string input_file_path = argv[1];
	int buffer_size = std::atoi(argv[2]);


	std::string* line_in_file = new std::string("");
	line_in_file -> reserve(1000);
	std::string* line_buffer = new std::string("");
	line_buffer -> reserve(buffer_size);

	int counter = 1;
	int position = 0 - buffer_size / 2;

	std::ifstream my_file;
	my_file.open(input_file_path);

	std::cout << "Position,AT Percentage,\n";

	
	while (std::getline(my_file, *line_in_file)){
		if ((*line_in_file)[0] != '>'){	
			for (int i = 0; i < line_in_file -> length(); i++){	
				if (counter < buffer_size){
					line_buffer -> push_back((*line_in_file)[i]);
					++counter;
				} else {
					line_buffer -> push_back((*line_in_file)[i]);
					position = position + buffer_size;
					double buffer_AT_percentage = calculate_AT_percentage(*line_buffer);
					std::cout << position << "," << buffer_AT_percentage <<",\n";
					counter = 1;
					delete line_buffer;
					line_buffer = new std::string("");
					line_buffer -> reserve(buffer_size);	
				}
			}
	
		}
	}

	
	double buffer_AT_percentage = calculate_AT_percentage(*line_buffer);
	position = position + line_buffer -> size();
	std::cout << position << "," << buffer_AT_percentage <<"\n";

	my_file.close();

	if (line_buffer != NULL) {delete line_buffer;}
	if (line_in_file != NULL) {delete line_in_file;}
	return 0;

}
