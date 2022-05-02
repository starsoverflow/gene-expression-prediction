#include <iostream>

#include <string>
#include <vector>

#include <sstream>
#include <iterator>
#include <fstream>

#include <map>
#include <set>
#include <unistd.h>

#include <future>
#include <filesystem>
namespace fs = std::filesystem;
using namespace std;

std::vector<std::string> getNextLineAndSplitIntoTokens(std::istream& str)
{
    std::vector<std::string>   result;
    std::string                line;
    std::getline(str,line);

    std::stringstream          lineStream(line);
    std::string                cell;

    while(std::getline(lineStream,cell, '\t'))
    {
        result.push_back(cell);
    }
    if (!lineStream && cell.empty())
    {
        // If there was a trailing comma then add an empty element.
        result.push_back("");
    }
    return result;
}


const int BIN_COUNT = 200;
const int BIN_LENGTH_BP = 100;
const int BIN_RANGE = 10000;  // +/- 10000 bp

struct rate {
    double met;
    double acc;
};

struct bin {
    string gene;
    int64_t id;
    bool operator==(const bin &o) const {
        return gene == o.gene && id == o.id;
    }

    bool operator<(const bin &o) const {
        return gene < o.gene || (gene == o.gene && id < o.id);
    }
};


int process_sample(bool isMet, ifstream& file2, const map<string, map<int64_t, string>> & tss_data, map<bin, rate>& bin_data, set<string> & appeared_genes)
{
    int count = 0;
    while (file2) {
        count++;
        if (count % 1000000 == 0) {
            printf("Processing line %d\n", count);

        }
        auto data = getNextLineAndSplitIntoTokens(file2);
        string chrom = data[0];
        if (chrom.empty()) {
            continue;
        }
        if (chrom == "chr") {
            continue;
        }
        int64_t position = std::stoll(data[1]);
        double rate = std::stod(data[2]);

        auto iterch = tss_data.find(chrom);
        if (iterch == tss_data.end()) continue;
        auto & tss_chrom = iterch->second;

        // find tss >= position - BIN_RANGE, tss + BIN_RANGE >= pos
        auto iter = tss_chrom.lower_bound(position - BIN_RANGE);
        if (iter == tss_chrom.end()) continue;
        auto first_tss = iter->first;
        while (true) {
            auto tss = iter->first;
            auto gene = iter->second;
            // tss - BIN_RANGE < position < tss + BIN_RANGE
            if (tss - BIN_RANGE >= position) {
                if (iter == tss_chrom.begin()) break;
                iter--;
                continue;
            }

            if (tss + BIN_RANGE < position) break;

            if (! ((tss-BIN_RANGE <= position) && (position <= tss+BIN_RANGE))) {
                printf("error!");
            }


            auto bin_id = max(0LL, position - (tss - BIN_RANGE)) / BIN_LENGTH_BP;
            if (bin_id >= BIN_COUNT) bin_id = BIN_COUNT - 1;
            if (isMet)
                bin_data[{gene, bin_id}].met += rate;
            else
                bin_data[{gene, bin_id}].acc += rate;
            if (bin_data.count({gene, 0}) == 0) {
                bin_data[{gene, 0}] = {0, 0};
            }
            appeared_genes.insert(gene);

            if (iter == tss_chrom.begin()) break;
            iter--;

        }
    }
    return 0;
}

// serum cells only
const string sample_names_tmp[] = {"A02","A03","A04","A05","A06","A07","A09","B01","B02","B03","B04","B05","B06","B07","B09","C01","C02","C03","C04","C05","C07","C09","D01","D05","D06","D07","D09","E01","E02","E03","E05","E07","E09","F02","F03","F04","F05","F06","F07","F09","G01","G02","G03","G05","G06","G07","G09","H02","H03","H09"};
const int sample_c = 50; // omit last two

const string DIR = "../../data/";
const string OUTPUTDIR = "../../processed_data/";

string get_sample_full_path(string sample_name)
{
    std::string path = DIR + "GSE109262_TSV/";
    for (const auto & entry : fs::directory_iterator(path)) {
        size_t pos;
        auto str = entry.path().string();
        if ((pos = str.find(sample_name + "_CpG-met_processed")) != string::npos) {
            return str.substr(0, pos);
        }
    }
    throw std::exception();
}

int write_one(string sample_name, map<bin, rate> & bin_data,  set<string> & intersect_genes)
{

    string output_file = OUTPUTDIR + "Cell_" + sample_name + ".all.csv";

    printf("Writing output to: %s\n", output_file.c_str());

    ofstream ofile(output_file, std::ofstream::trunc);


    string previous_bin = "";
    int64_t previous_number = -1;
    int64_t rows_written = 0;

    for (auto & bin: bin_data) {
        if (intersect_genes.find(bin.first.gene) == intersect_genes.end()) {
            continue;
        }
        if (previous_bin != bin.first.gene && previous_bin != "") {
            for (int64_t i = previous_number + 1; i < BIN_COUNT; i++) {
                string bin_id_str = std::to_string(i);
                string bin_name = previous_bin + "\t" + std::string(3-bin_id_str.length(), '0').append(bin_id_str);
                ofile << bin_name << "\t0\t0\n";
                rows_written++;
            }
            previous_number = -1;
        }

        if (bin.first.id != previous_number + 1) {
            for (int64_t i = previous_number + 1; i < bin.first.id; i++) {
                string bin_id_str = std::to_string(i);
                string bin_name = bin.first.gene + "\t" + std::string(3-bin_id_str.length(), '0').append(bin_id_str);
                ofile << bin_name << "\t0\t0\n";
                rows_written++;
            }
        }
        previous_number = bin.first.id;

        previous_bin = bin.first.gene;
        string bin_id_str = std::to_string(bin.first.id);
        string bin_name = bin.first.gene + "\t" + std::string(3-bin_id_str.length(), '0').append(bin_id_str);

        ofile << bin_name << "\t" << bin.second.met << "\t" << bin.second.acc << "\n";
        rows_written++;
    }
    for (int64_t i = previous_number + 1; i < BIN_COUNT; i++) {
        string bin_id_str = std::to_string(i);
        string bin_name = previous_bin + "\t" + std::string(3-bin_id_str.length(), '0').append(bin_id_str);
        ofile << bin_name << "\t0\t0\n";
        rows_written++;
    }
    previous_number = -1;
    cout << "Write " << rows_written << " rows\n";

    return 0;
}

map<bin, rate> process_one(string sample_name, string prefix,
                const map<string, map<int64_t, string>> & tss_data, set<string> & appeared_genes)
{

    cout << "Processing " << sample_name << endl;
    map<bin, rate> bin_data;

    std::ifstream file2(prefix + sample_name + "_CpG-met_processed.tsv");
    process_sample(true, file2, tss_data, bin_data, appeared_genes);
    std::ifstream file3(prefix + sample_name + "_GpC-acc_processed.tsv");
    process_sample(false, file3, tss_data, bin_data, appeared_genes);

    cout << "Processing done " << sample_name << endl;
    return bin_data;
}

int main()
{

    map<string, map<int64_t, string>> tss_data;

    std::ifstream file(DIR + "TSS_clean.csv");
    std::string cwd = getcwd(0, 0);
    auto ignored_ = getNextLineAndSplitIntoTokens(file);
    while (file) {
        auto data = getNextLineAndSplitIntoTokens(file);
        if (data.size() < 4) continue;
        string gene = data[1];
        if (gene.empty()) {
            continue;
        }
        int64_t tss = std::stoll(data[2]);
        string chrom = data[3];
        if (chrom.empty()) {
            continue;
        }
        tss_data[chrom][tss] = gene;
    }

    vector<map<bin, rate>> all_bin_data(sample_c);
    vector<future<set<string>>> futures(sample_c);

    vector<string> sample_names;
    for (int i = 0; i < sample_c; i++) {
        sample_names.push_back(string("ESC_" + sample_names_tmp[i]));
    }


    for (int i = 0; i < sample_c; i++) {
        futures.at(i) = async([i, &tss_data, &all_bin_data, &sample_names ]{
            set<string> appeared_genes;
            try {
                auto bin_data = process_one(sample_names[i], get_sample_full_path(sample_names[i]), tss_data, appeared_genes);
                all_bin_data.at(i) = std::move(bin_data);
            } catch (...) {
                cout << "exception\n";
            }
            return appeared_genes;
        });
        // futures.at(i).wait();
    }

    set<string> intersect_genes;

    for (int i = 0; i < sample_c; i++) {
        futures[i].wait();
        auto result = futures[i].get();
        if (i == 0)
            intersect_genes.insert(result.begin(), result.end());
        else {
            set<string> intersect;
            set_intersection(intersect_genes.begin(), intersect_genes.end(), result.begin(), result.end(),
                             std::inserter(intersect, intersect.begin()));
            intersect_genes = std::move(intersect);
        }
    }

    cout << intersect_genes.size() << " genes for each cell\n";


    /* add missing genes
    for (int i = 0; i < sample_c; i++) {
        for (auto & gene : all_genes) {
            auto & bin_data = all_bin_data[i];
            if (bin_data.count({gene, 0}) == 0)
                bin_data[{gene, 0}] = {0, 0};
        }
    }
     */

    for (int i = 0; i < sample_c; i++) {
        write_one(sample_names[i], all_bin_data[i], intersect_genes);
    }

    return 0;
}
