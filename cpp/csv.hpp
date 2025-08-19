#pragma once
#include <bits/stdc++.h>
using namespace std;

inline vector<vector<string>> read_csv(const string &filename) {
    ifstream file(filename);
    vector<vector<string>> rows;
    string line;
    while (getline(file, line)) {
        stringstream ss(line);
        string cell;
        vector<string> row;
        while (getline(ss, cell, ',')) row.push_back(cell);
        rows.push_back(row);
    }
    return rows;
}
