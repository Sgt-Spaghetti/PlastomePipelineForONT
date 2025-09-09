#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <vector>

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

bool CheckReadParameters(const string &Read, const int &ATminimum, const int &ATmaximum, const int &MinimumLength){ // returns true to discard a read, or false to keep it
	const size_t ReadLength{Read.length()};
	if (ReadLength < MinimumLength){
		return true;
	}
	else {
		int ReadATcount{0};
		for (int i=0; i<ReadLength; i++){ // Calculate AT percentage of the sequence
			if (Read[i] == 'A' || Read[i] == 'T' || Read[i] == 'U'){
				++ReadATcount;
			}
		}
                const int ReadATPercentage{int(((double(ReadATcount) / double(ReadLength)) * 100) + 0.5)};
		// Adding 0.5 allows rounding behaviour when casting
                if (ReadATPercentage < ATminimum || ReadATPercentage >= ATmaximum){
					return true;
				}
				else {
					return false;
		}
	}
}

bool CheckReadQuality(const string &ReadQuality, const int &MinimumQualityScore){
	const int AsciiCodeOffsetForPhredScore{33};
	const size_t QualityScoreLength{ReadQuality.length()};
	int Qscore{0};
	for (int i=0; i<QualityScoreLength; i++){
		Qscore = Qscore + int(ReadQuality[i]) - AsciiCodeOffsetForPhredScore;
	}
	const int Qmean{int((double(Qscore)/double(QualityScoreLength))+0.5)};
	if (Qmean < MinimumQualityScore){
		return true;
	}
	else {
		return false;
	}
}

int main(int argc, char* argv[]) {

	const string InputFilePath{argv[1]};
	const int ATmin{atoi(argv[2])};
	const int ATmax{atoi(argv[3])};
	const int LengthMin{atoi(argv[4])};
	const int QualityMin{atoi(argv[5])};
	const string OutputFilePath{argv[6]};
	
	const int LinesPerEntry{4};
	int PointInEntry{0};

	const int LocationOfSequence{1};
	bool DiscardEntry{false};

	vector<string> EntriesToKeep;

	ifstream InputFile(InputFilePath);

	string *LineInFile = new string();

	Entry * pnew_entry = new Entry();

	while (getline(InputFile, *LineInFile)){
		if (PointInEntry<LinesPerEntry-1){
			if (PointInEntry == 0){
				 pnew_entry -> set_name(LineInFile);
			} else if (PointInEntry == 1){
				pnew_entry -> set_sequence(LineInFile);
			} else if (PointInEntry == 2){
				pnew_entry -> set_adapter(LineInFile);
			}

			if (PointInEntry == LocationOfSequence){
				DiscardEntry = CheckReadParameters(*(pnew_entry -> get_sequence()), ATmin, ATmax, LengthMin);
			}
			++PointInEntry;
		}
		else { // We are in the last section of the entry, prepare reset
			pnew_entry -> set_quality(LineInFile); // We are at the quality score!
			if (DiscardEntry == false){ // It has passed the AT etc checks
				if (CheckReadQuality(*(pnew_entry -> get_quality()), QualityMin) == false){ // Can it pass the quality score check?
					EntriesToKeep.push_back(*(pnew_entry -> get_name()));
					EntriesToKeep.push_back(*(pnew_entry -> get_sequence()));
					EntriesToKeep.push_back(*(pnew_entry -> get_adapter()));
					EntriesToKeep.push_back(*(pnew_entry -> get_quality()));
				}
			}
			PointInEntry = 0;
			delete pnew_entry;
			Entry * pnew_entry = new Entry();
		}
	LineInFile = new string();
	}

	InputFile.close();

	
	ofstream OutputFile(OutputFilePath);
	for (int i=0; i < EntriesToKeep.size(); i++){
		OutputFile << EntriesToKeep[i] << "\n";
	}
	OutputFile.close();


	return 0;	
}
