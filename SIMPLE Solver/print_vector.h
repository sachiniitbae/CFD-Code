#include <iostream>
#include <vector>
using namespace std;

void print_vector_1d(vector<float> vec)
{
    cout << "Shape: ( " <<vec.size()<< ",1)" << endl;
    for (float d : vec)
        cout << d << "\t";
    cout << endl;
}

void print_vector_2d(vector<vector<float>> vec)
{
    cout << "Shape:(" <<vec.size()<< "," << vec[0].size() <<")" << endl;
    for(int i=0 ;i<vec.size();i++)
    {
        for(int j=0;j<vec[0].size();j++)
        {
            cout << vec[i][j] << "\t";
        }
        cout << endl;
    }
}

