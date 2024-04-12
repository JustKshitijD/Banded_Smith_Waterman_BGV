#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>

using namespace std;

// Function to split a string into individual digits
vector<int> splitDigits(const string& str) {
    vector<int> digits;
    for (char digit : str) {
        if (isdigit(digit)) {
            digits.push_back(digit - '0');
        }
    }
    return digits;
}

int main() {
    ifstream inputFile("number_nucl_file"); // Replace "your_file.txt" with the actual file path
    if (!inputFile) {
        cerr << "Error opening file." << endl;
        return 1;
    }

    int count = 10; // Number of lines to read
    int cnt=count;
    int size=3;
    int t_arr[count][size]; // Final array to store arrays from the first section
    int q_arr[count][size]; // Final array to store arrays from the second section

    string line;
    while (count > 0 && getline(inputFile, line)) {
        if (line.find("-----") != string::npos) {
            break; // Stop reading when encountering a line containing "-----"
        }
        vector<int> digits = splitDigits(line);

        for(int i=0;i<size;i++)
        {
            if(i<digits.size())
                t_arr[cnt-count][i]=digits[i];
            else
                t_arr[cnt-count][i]=5;
        }
        // t_arr.push_back(digits);
        count--;
    }

    while(getline(inputFile, line))
    {
        if (line.find("-----") != string::npos) {
            break; // Stop reading when encountering a line containing "-----"
        }
    }

    count = cnt; // Reset count for the second section
    while (count > 0 && getline(inputFile, line)) {
        vector<int> digits = splitDigits(line);
        for(int i=0;i<size;i++)
        {
            if(i<digits.size())
                q_arr[cnt-count][i]=digits[i];
            else
                q_arr[cnt-count][i]=6;
        }
        // t_arr.push_back(digits);
        count--;
    }

    for(int i=0;i<cnt;i++)
    {
        for(int j=0;j<size;j++)
        {
            cout<<t_arr[i][j]<<" ";
        }
        cout<<endl;
    }

    cout << endl;

    for(int i=0;i<cnt;i++)
    {
        for(int j=0;j<size;j++)
        {
            cout<<q_arr[i][j]<<" ";
        }
        cout<<endl;
    }

    return 0;
}
