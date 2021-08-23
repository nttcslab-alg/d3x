#include <math.h>
#include <unistd.h>

#include <chrono>
#include <unordered_set>

#include "dancing_on_zdd.h"
#include "dp_manager.h"

/**
 * main function
 */

// extern uint64_t num_search_tree_nodes;
// extern uint64_t search_tree_depth;
// extern uint64_t search_tree_max_depth;
// extern uint64_t num_solutions;
// extern uint64_t ZddWithLinks::num_updates;
// extern uint64_t num_inactive_updates;

/**
 * get number of variables from zdd file
 */
int get_num_vars_from_zdd_file(const string& file_name) {
    ifstream ifs(file_name);

    if (!ifs) {
        cerr << "can't open " << file_name << endl;
        exit(1);
    }

    string line;

    unordered_set<int> vars;

    while (getline(ifs, line)) {
        if (line[0] == '.' || line[0] == '\n' || line.size() == 0) continue;

        istringstream iss(line);
        int nid;
        int var;
        string lo_str;
        int lo_id;
        string hi_str;
        int hi_id;
        iss >> nid;
        iss >> var;

        vars.emplace(var);
    }

    return vars.size();
}

void show_help_and_exit() {
    std::cerr << "usage: ./dancing_on_zdd_main -z zdd_file\n" << std::endl;
    exit(1);
}

int main(int argc, char** argv) {
    int opt;
    string zdd_file_name;
    int num_var = -1;
    while ((opt = getopt(argc, argv, "z:h")) != -1) {
        switch (opt) {
            case 'z':
                zdd_file_name = optarg;
                break;
            case 'h':
                show_help_and_exit();
                break;
        }
    }

    if (zdd_file_name.empty()) {
        show_help_and_exit();
    }

    num_var = get_num_vars_from_zdd_file(zdd_file_name);

    ZddWithLinks zdd_with_links(num_var, false);
    zdd_with_links.load_zdd_from_file(zdd_file_name);
    if (zdd_with_links.sanity()) {
        fprintf(stderr, "initial zdd is invalid\n");
    }
    fprintf(stderr, "load files done\n");
    vector<vector<uint16_t>> solution;
    auto start_time = std::chrono::system_clock::now();
    zdd_with_links.search(solution, 0);
    auto end_time = std::chrono::system_clock::now();
    printf("num nodes %llu, num solutions %llu, num updates %llu, "
           "time: %llu msecs\n", ZddWithLinks::num_search_tree_nodes,
           ZddWithLinks::num_solutions, ZddWithLinks::num_updates,
           std::chrono::duration_cast<std::chrono::milliseconds>(end_time -
                                                                 start_time)
               .count());

    return 0;
}