#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <iterator>
#include <sstream>
#include <stdlib.h>
#include <Rcpp.h>

using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
List read_PLY(std::string input,std::string output)
{
	ifstream file(input.c_str());
    string line;
    vector<string> tokens;
    int i;
	int j;
	FILE *ptr_myfile;
	FILE *ptr_mynewfile;
	List out(3);

	float my_record;
	int intmy_record;
	unsigned char character;
    int contador=0;
    for(i=1;i<=11;i++)
    {
            getline(file, line);
            istringstream iss(line);
            copy(istream_iterator<string>(iss),
            istream_iterator<string>(),
            back_inserter<vector<string> >(tokens));
    }
    file.close();
    int vertex=atoi(tokens[9].c_str());
    int faces=atoi(tokens[24].c_str());
	
    for(int j=0;j<tokens.size();j++)
    {
        contador+=tokens[j].size()*sizeof(char);
        contador+=sizeof(char);
    }
		
	ptr_myfile=fopen(input.c_str(),"rb");
    ptr_mynewfile=fopen(output.c_str(),"w");

	fseek(ptr_myfile,contador,SEEK_SET);
	for ( i=1; i <= (vertex); i++)
	{
		for(j=1;j<=4;j++){
			fread(&my_record,sizeof(float),1,ptr_myfile);
			fprintf(ptr_mynewfile, "%f ",my_record);
		}
		fprintf(ptr_mynewfile, "\n");
	}

	for ( i=1; i <= (faces); i++)
	{
		fread(&character,sizeof(char),1,ptr_myfile);
		fprintf(ptr_mynewfile, "%u ",character);
		for(j=1;j<=3;j++){
		    fread(&intmy_record,sizeof(int),1,ptr_myfile);
		    fprintf(ptr_mynewfile, "%i ",intmy_record);
		}
    fprintf(ptr_mynewfile, "\n");
	}
	
	fclose(ptr_myfile);
	fclose(ptr_mynewfile);
	out[0]=vertex;
	out[1]=faces;
	out[2]=output;
	return out;
}
