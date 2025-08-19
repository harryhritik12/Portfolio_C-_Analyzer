#include <bits/stdc++.h>
#include "csv.hpp"
#include "json.hpp"
using namespace std;
using json = nlohmann::json;

// Convert string to double safely
inline double to_double(const string &s) {
    try { return stod(s); } catch (...) { return 0.0; }
}

// Compute daily returns
vector<vector<double>> compute_returns(const vector<vector<string>> &prices) {
    int T = prices.size() - 1;   // exclude header
    int N = prices[0].size() - 1; // exclude date
    vector<vector<double>> R(T-1, vector<double>(N));

    for (int t = 1; t < T; t++) {
        for (int j = 1; j <= N; j++) {
            double p0 = to_double(prices[t][j]);
            double p1 = to_double(prices[t+1][j]);
            R[t-1][j-1] = log(p1 / p0); // log returns
        }
    }
    return R;
}

// Monte Carlo Simulation
vector<double> simulate_portfolio(
    int n_sims, const vector<vector<double>> &R, const vector<double> &weights, double V0)
{
    int T = R.size();
    int N = weights.size();

    // Compute mean and variance of portfolio returns
    vector<double> portR(T);
    for (int t = 0; t < T; t++) {
        double r = 0;
        for (int j = 0; j < N; j++) r += weights[j] * R[t][j];
        portR[t] = r;
    }

    double mu = accumulate(portR.begin(), portR.end(), 0.0) / T;
    double var = 0;
    for (double r : portR) var += (r - mu) * (r - mu);
    var /= (T - 1);
    double sigma = sqrt(var);

    // Monte Carlo simulation
    mt19937_64 rng(random_device{}());
    normal_distribution<double> dist(mu, sigma);

    vector<double> pnl;
    pnl.reserve(n_sims);
    for (int i = 0; i < n_sims; i++) {
        double ret = dist(rng);
        pnl.push_back(V0 * ret);
    }
    return pnl;
}

// Percentile
double percentile(vector<double> x, double p) {
    sort(x.begin(), x.end());
    double idx = p * (x.size()-1);
    size_t i = floor(idx), j = ceil(idx);
    if (i == j) return x[i];
    double w = idx - i;
    return x[i]*(1-w) + x[j]*w;
}

int main() {
    try {
        // Load prices
        auto prices = read_csv("../data/prices_100.csv");
        if (prices.empty()) throw runtime_error("No prices.csv found");

        // Load weights
        auto w_csv = read_csv("../data/weights.csv");
        vector<double> weights;
        for (string s : w_csv[0]) weights.push_back(to_double(s));

        // Compute returns
        auto R = compute_returns(prices);

        // Simulate
        double V0 = 100000.0;
        auto pnl = simulate_portfolio(50000, R, weights, V0);

        // Compute VaR, CVaR
        double var95 = -percentile(pnl, 0.05);
        double var99 = -percentile(pnl, 0.01);

        double q95 = percentile(pnl, 0.05);
        double sum95 = 0; int c95 = 0;
        for (double v : pnl) if (v <= q95) { sum95 += -v; c95++; }
        double cvar95 = (c95>0 ? sum95/c95 : -q95);

        double q99 = percentile(pnl, 0.01);
        double sum99 = 0; int c99 = 0;
        for (double v : pnl) if (v <= q99) { sum99 += -v; c99++; }
        double cvar99 = (c99>0 ? sum99/c99 : -q99);

        // Histogram
        int bins = 30;
        double minv = *min_element(pnl.begin(), pnl.end());
        double maxv = *max_element(pnl.begin(), pnl.end());
        double width = (maxv - minv) / bins;
        vector<int> counts(bins, 0);
        for (double v : pnl) {
            int b = (int)((v - minv) / width);
            if (b >= bins) b = bins - 1;
            if (b < 0) b = 0;
            counts[b]++;
        }

        json hist = json::array();
        for (int i = 0; i < bins; i++) {
            hist.push_back({{"binStart", minv+i*width}, {"binEnd", minv+(i+1)*width}, {"count", counts[i]}});
        }

        // Output JSON
        json out;
        out["var95"] = var95;
        out["var99"] = var99;
        out["cvar95"] = cvar95;
        out["cvar99"] = cvar99;
        out["histogram"] = hist;

        ofstream fo("../web/public/results.json");
        fo << out.dump(2);
        cout << "Results written to ../web/public/results.json\n";
    }
    catch (exception &e) {
        cerr << "Error: " << e.what() << endl;
    }
}
