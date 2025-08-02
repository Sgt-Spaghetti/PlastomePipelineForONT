#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <vector>

using namespace std;

bool CheckReadParameters(const string &Read, const double &ATminimum, const double &ATmaximum, const double &MinimumLength){ // returns true to discard a read, or false to keep it
	const size_t ReadLength{Read.length()};
	if (ReadLength < MinimumLength){
		return true;
	}
	else {
		int ReadATcount{0};
		for (int i=0; i<ReadLength; i++){ // Calculate AT percentage of the sequence
			if (Read[i] == 'A' || Read[i] == 'T'){
				++ReadATcount;
			}
		}
		const double ReadATPercentage{100.0 * (double(ReadATcount) / double(ReadLength))};
		if (ReadATPercentage < ATminimum || ReadATPercentage >= ATmaximum){
			return true;
		}
		else {
			return false;
		}
	}
}

bool CheckReadQuality(const string &ReadQuality, const double &MinimumQualityScore){
	const int AsciiCodeOffsetForPhredScore{33};
	const size_t QualityScoreLength{ReadQuality.length()};
	int Qscore{0};
	for (int i=0; i<QualityScoreLength; i++){
		Qscore = Qscore + int(ReadQuality[i]) - AsciiCodeOffsetForPhredScore;
	}
	const double Qmean{double(Qscore)/double(QualityScoreLength)};
	if (Qmean < MinimumQualityScore){
		return true;
	}
	else {
		return false;
	}
}

int main(int argc, char* argv[]) {

	const string InputFilePath{argv[1]};
	const double ATmin{atof(argv[2])};
	const double ATmax{atof(argv[3])};
	const double LengthMin{atof(argv[4])};
	const double QualityMin{atof(argv[5])};
	const string OutputFilePath{argv[6]};
	
	const int LinesPerEntry{4};
	int PointInEntry{0};
	string Entry[4]{"","","",""};
	const int LocationOfSequence{1};
	bool DiscardEntry{false};

	vector<string> EntriesToKeep;

	ifstream InputFile(InputFilePath);

	string LineInFile;

	while (getline(InputFile, LineInFile)){
		if (PointInEntry<LinesPerEntry-1){
			Entry[PointInEntry] = LineInFile;
			if (PointInEntry == LocationOfSequence){
				DiscardEntry = CheckReadParameters(Entry[PointInEntry], ATmin, ATmax, LengthMin);
			}
			++PointInEntry;
		}
		else { // We are in the last section of the entry, prepare reset
			Entry[PointInEntry] = LineInFile; // We are at the quality score!
			if (DiscardEntry == false){ // It has passed the AT etc checks
				if (CheckReadQuality(Entry[PointInEntry], QualityMin) == false){ // Can it pass the quality score check?
					for (int i = 0; i < LinesPerEntry; i++){
						EntriesToKeep.push_back(Entry[i]);
					}
				}
			}
			PointInEntry = 0;
		}
	}

	InputFile.close();

	
	ofstream OutputFile(OutputFilePath);
	for (int i=0; i < EntriesToKeep.size(); i++){
		OutputFile << EntriesToKeep[i] << "\n";
	}
	OutputFile.close();


	return 0;	
}
