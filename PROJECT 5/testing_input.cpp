#include<armadillo>
#include<iostream>
#include<string>
#include<fstream>
#include<vector>
#include"Box.hpp"
#include"utils.hpp"
using namespace arma;
using namespace std;


vector<string> split (const string &s, char delim) {
    vector<string> result;
    stringstream ss (s);
    string item;

    while (getline (ss, item, delim)) {
        result.push_back (item);
    }

    return result;
}


int main(int argc, const char* argv[]) {

	string init_conditions = argv[1];
	fstream initial_conditions; 
	initial_conditions.open(init_conditions, ios::in);
	string input; 
	string values; 
	string delimiter = ","; 
	const char* delim = ","; 
	if (!initial_conditions) {
		cout << "No such file";
	}
	else {
		getline(initial_conditions, input); 
		getline(initial_conditions, values); 
	
	}
	vector<string> vals = split(values, ','); 
	cout << vals[0] << endl; 
	cout << input << endl; 
	cout << values << endl; 
	initial_conditions.close(); 
	string v0 = values.substr(0, values.find(delimiter));
	string d_t = values.substr(1, values.find(delimiter)); 
	string len_time = values.substr(2,values.find(delimiter)); 
	string center_x = values.substr(3,values.find(delimiter)); 
	string center_y = values.substr(4,values.find(delimiter)); 
	string spread_x = values.substr(5,values.find(delimiter)); 
	string spread_y = values.substr(6,values.find(delimiter)); 
	string m_x = values.substr(7, values.find(delimiter)); 
	string m_y = values.substr(8, values.find(delimiter)); 
	
	cout << vals[0] << endl; 
	cout << vals[1] << endl; 
	cout << vals[2] << endl; 
	cout << vals[3] + vals[4] + vals[5] + vals[6] << endl; 
	cout << vals[7] + vals[8] << endl;  
	
	
	return 0; 	
}
