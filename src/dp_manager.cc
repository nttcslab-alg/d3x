#include "dp_manager.h"

using namespace std;

DpManager::DpManager(const vector<Node> &nodes, const int num_var)
    : table_elems_(nodes.size(), 0),
      var_heads_(num_var + 1, 0),
      num_elems_(num_var + 1, 0),
      diff_counter_(nodes.size(), 0),
      diff_counter_hi_(nodes.size(), 0),
      entries_counter_(0),
      num_var_(num_var),
      upper_varorder_pq_(),
      lower_varorder_pq_() {
    int previous_var = -1;

    for (size_t i = 0; i < nodes.size(); i++) {
        const Node &node = nodes[i];
        if (node.var != previous_var) {
            var_heads_[node.var] = i;
            previous_var = node.var;
        }
    }
}
