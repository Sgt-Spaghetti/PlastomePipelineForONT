#include <stdlib.h>
#include <iostream>
#include <fstream>

using namespace std;

class Entry {

	private:
		string * pread_name;
		string * pread_sequence;
		string * pread_adapter;
		string * pread_quality;
		
	public:
		~Entry(){
			if(pread_name!=NULL)  delete pread_name;
			if(pread_sequence!=NULL)  delete pread_sequence;
			if(pread_adapter!=NULL) delete pread_adapter;
			if(pread_quality!=NULL) delete pread_quality;
		}
		void set_name(string * n){
			pread_name = n;
		}
		void set_sequence(string * n){
			pread_sequence = n;
		}
		void set_adapter(string * n){
			pread_adapter = n;
		}
		void set_quality(string * n){
			pread_quality = n;
		}
		string * get_name(){
			return pread_name;
		}
		string * get_sequence(){
			return pread_sequence;
		}
		string * get_adapter(){
			return pread_adapter;
		}
		string * get_quality(){
			return pread_quality;
		}
};

bool CheckReadParameters(const string &Read, const int ATminimum, const int ATmaximum, const int MinimumLength, const unsigned int ReadLength){ // returns true to discard a read, or false to keep it
	if (ReadLength < MinimumLength){
		return true;
	}
	else {
		unsigned int ReadATcount = 0;
		for (unsigned int i=0; i<ReadLength; i++){ // Calculate AT percentage of the sequence
			char Base = Read[i];
			if (Base == 'A' || Base == 'T' || Base == 'U'){
				++ReadATcount;
			}
		}
                const double ReadATPercentage = (double(ReadATcount) / double(ReadLength)) * 100;
                if (ReadATPercentage < ATminimum || ReadATPercentage >= ATmaximum){
					return true;
				}
				else {
					return false;
		}
	}
}

bool CheckReadQuality(const string &ReadQuality, const int MinimumQualityScore, const unsigned int QualityScoreLength){
	const int AsciiCodeOffsetForPhredScore = 33;
	unsigned int Qscore = 0;
	for (unsigned int i=0; i<QualityScoreLength; i++){
		Qscore = Qscore + int(ReadQuality[i]) - AsciiCodeOffsetForPhredScore;
	}
	const double Qmean = double(Qscore)/double(QualityScoreLength);
	if (Qmean < MinimumQualityScore){
		return true;
	}
	else {
		return false;
	}
}

int main(int argc, char* argv[]) {

	const string InputFilePath = argv[1];
	const int ATmin = atoi(argv[2]);
	const int ATmax = atoi(argv[3]);
	const int LengthMin = atoi(argv[4]);
	const int QualityMin = atoi(argv[5]);
	const string OutputFilePath = argv[6];
	
	const int LinesPerEntry = 4;
	int PointInEntry = 0;

	const int LocationOfSequence = 1;
	bool DiscardEntry = false;

	ifstream InputFile(InputFilePath);

	string *LineInFile = new string();
	LineInFile -> reserve (60000);

	Entry * pnew_entry = new Entry();

	unsigned int SequenceLength = 0;

	while (getline(InputFile, *LineInFile)){
		LineInFile -> shrink_to_fit();
		if (PointInEntry<LinesPerEntry-1){
			if (PointInEntry == 0){
				 pnew_entry -> set_name(LineInFile);
			} else if (PointInEntry == 1){
				pnew_entry -> set_sequence(LineInFile);
				SequenceLength = (*(pnew_entry -> get_sequence())).length(); 
				DiscardEntry = CheckReadParameters(*(pnew_entry -> get_sequence()), ATmin, ATmax, LengthMin, SequenceLength);
			} else if (PointInEntry == 2){
				pnew_entry -> set_adapter(LineInFile);
			}

			++PointInEntry;
		}
		else { // We are in the last section of the entry, prepare reset
			if (DiscardEntry == false){
				pnew_entry -> set_quality(LineInFile); // We are at the quality score!
				if (CheckReadQuality(*(pnew_entry -> get_quality()), QualityMin, SequenceLength) == false){ // Can it pass the quality score check?
				
					ofstream OutputFile(OutputFilePath, std::ios_base::app);
					OutputFile << *(pnew_entry -> get_name()) << "\n";
					OutputFile << *(pnew_entry -> get_sequence()) << "\n";
					OutputFile << *(pnew_entry -> get_adapter()) << "\n";
					OutputFile << *(pnew_entry -> get_quality()) << "\n";;
	
					OutputFile.flush();
					OutputFile.close();
				}
			}
			SequenceLength = 0;
			PointInEntry = 0;
			delete pnew_entry;
			pnew_entry = new Entry();
		}
	LineInFile = new string();
	LineInFile -> reserve(60000);
	}

	InputFile.close();

	return 0;	
}
