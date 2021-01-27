/* C++ script to write a file "ts.remove" containing the indices of all TSs corresponding to a barrier
   greater than a specified threshold, for either of the pair of minima which the TS connects, and
   all TSs with energy greater than a second specified threshold.

   Compile with:
   g++ -std=c++11 -o tidy_maxbarrier tidy_maxbarrier.cpp

   Example usage (e.g. barrier threshold = 50 kcal/mol, raw TS energy threshold = -2050 kcal/mol):
   tidy_maxbarrier min.data ts.data 50.0 -2050.0

   Daniel J. Sharpe */

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <string>
#include <tuple>
#include <vector>
using namespace std;

bool check_barrier(double *min_energies, pair<int,int> conn_min, double ts_en, \
                   double thresh) {

    if ( (ts_en - min_energies[conn_min.first-1] > thresh) || \
       (ts_en - min_energies[conn_min.second-1] > thresh) ) {
        return true; }
    else {
        return false; }
}

// return an array containing minima energies
double *read_min_info(string mindataf) {

    int n_min = 0;
    string line;
    ifstream file(mindataf);
    while (!file.eof()) {
        getline(file,line);
        n_min++; }
    cout << "File " << mindataf << " contains " << n_min-1 << " minima\n";

    double *min_energies = new double[n_min-1];

    int i = 0;
    file.clear(); // clear state of stream
    file.seekg(0,file.beg); // go back to first line of file
    while (getline(file,line)) {
        string firstcol;
        stringstream(line) >> firstcol;
        // cout << firstcol << "\n";
        min_energies[i] = stod(firstcol);
        // cout << std::setprecision(10) << min_energies[i] << "\n";
        i++;
    }   
    cout << "Read in energies of " << i << " minima from file " << mindataf << "\n";

    return min_energies;
}

// write a file ts.remove containing indices of TSs associated with a barrier
// greater than the specified threshold
void parse_ts_info(string tsdataf, double *min_energies, double thresh, \
                   double raw_thresh) {

    vector<int> bad_ts; // list of TSs associated with a barrier or raw energy above the thresholds
    pair<int,int> conn_min; // pair of minima connected by the TS
    double ts_en; // energy of the TS

    // find "bad" TSs
    int n_ts = 0, n_bad_ts = 0;
    string line;
    ifstream file(tsdataf);
    while (getline(file,line)) {
        string firstcol, secondcol, thirdcol, fourthcol, fifthcol;
        stringstream(line) >> firstcol >> secondcol >> thirdcol >> fourthcol >> fifthcol;
        conn_min = make_pair(stod(fourthcol),stod(fifthcol));
        ts_en = stod(firstcol);
        if ( (check_barrier(min_energies,conn_min,ts_en,thresh)) ||
             (ts_en > raw_thresh) ) {
            n_bad_ts++;
            bad_ts.emplace_back(n_ts+1); }
        n_ts++;
    }
    cout << "Finished processing " << n_ts << " transition states from file " << \
            tsdataf << "\n";
    cout << "There are " << n_bad_ts << " transition states associated with barriers above " << \
            "the threshold of " << thresh << " kcal/mol\n";

    // write out ts.remove file
    ofstream tsremovef;
    tsremovef.open("ts.remove",std::ios::app); // append mode
    tsremovef << n_bad_ts << std::endl;
    for (auto ts_idx: bad_ts) {
        tsremovef << ts_idx << std::endl;
    }
    tsremovef.close();
}

int main(int argc, char** argv) {

    // read arguments to script
    string mindataf = argv[1];
    string tsdataf = argv[2];
    double thresh = atof(argv[3]);
    double raw_thresh = atof(argv[4]);

    // read the min.data file
    double *min_energies = read_min_info(mindataf);

    // parse the ts.data file
    parse_ts_info(tsdataf,min_energies,thresh,raw_thresh);

    // cleanup
    delete min_energies;

    return 0;
}
