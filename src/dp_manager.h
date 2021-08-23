#ifndef DP_MANAGER_H_
#define DP_MANAGER_H_

#include <memory>
#include <queue>
#include <vector>

#include "dancing_on_zdd.h"

/**
 * A class having dp tables for cover / uncover operaitons.
 *  When performing cover/uncover operations, the order must be reversed. 
 *  the class stores the order of processed node cell ids.
 */

struct Node;
class DpManager {
   public:
    DpManager(const std::vector<Node> &nodes, const int num_var);
    DpManager(const DpManager &obj) = delete;

    void add_node_diff_count(uint16_t var, int32_t node_id, count_t count) {
        diff_counter_[node_id] += count;
        if (diff_counter_[node_id] > count) {
            return;
        }

        auto var_num_elems = num_elems_[var];
        num_elems_[var]++;
        table_elems_[var_heads_[var] + var_num_elems] = node_id;
        if (!var_num_elems) {
            lower_varorder_pq_.push(var);
        }
    }

    void add_node_diff_count_high(uint16_t var, int32_t node_id,
                                  count_t count) {
        if (diff_counter_[node_id] > 0 || diff_counter_hi_[node_id] > 0) {
            diff_counter_hi_[node_id] += count;
            return;
        }

        diff_counter_hi_[node_id] += count;
        table_elems_[var_heads_[var] + num_elems_[var]] = node_id;
        num_elems_[var]++;
        if (num_elems_[var] == 1) {
            upper_varorder_pq_.push(var);
        }
    }

    void add_node_diff_count_low(uint16_t var, int32_t node_id, count_t count) {
        if (diff_counter_[node_id] > 0 || diff_counter_hi_[node_id] > 0) {
            diff_counter_[node_id] += count;
            return;
        }

        diff_counter_[node_id] += count;
        table_elems_[var_heads_[var] + num_elems_[var]] = node_id;
        num_elems_[var]++;
        if (num_elems_[var] == 1) {
            upper_varorder_pq_.push(var);
        }
    }

    inline int32_t at(uint16_t var, int32_t i) const noexcept {
        return table_elems_[var_heads_[var] + i];
    }

    int32_t num_elems(uint16_t var) const { return num_elems_[var]; }

    count_t count_at(int32_t node_id) const { return diff_counter_[node_id]; }

    count_t low_count_at(int32_t node_id) const {
        return diff_counter_[node_id];
    }

    count_t high_count_at(int32_t node_id) const {
        return diff_counter_hi_[node_id];
    }

    count_t get_count_and_clear(const int32_t node_id) {
        auto c = diff_counter_[node_id];
        diff_counter_[node_id] = 0;
        return c;
    }

    count_t get_low_count_and_clear(const int32_t node_id) {
        auto c = diff_counter_[node_id];
        diff_counter_[node_id] = 0;
        return c;
    }

    count_t get_high_count_and_clear(const int32_t node_id) {
        auto c = diff_counter_hi_[node_id];
        diff_counter_hi_[node_id] = 0;
        return c;
    }

    void clear_var_counter(uint16_t var) { num_elems_[var] = 0; }

    void clear_var_elems(uint16_t var) {
        for (size_t i = 0; i < num_elems_[var]; i++) {
            auto node_id = at(var, i);
            diff_counter_[node_id] = 0UL;
            diff_counter_hi_[node_id] = 0UL;
        }
        num_elems_[var] = 0;
    }

    uint16_t upper_nonzero_var() {
        if (upper_varorder_pq_.empty()) return 0;

        uint16_t next = upper_varorder_pq_.top();
        upper_varorder_pq_.pop();
        return next;
    }


    uint16_t lower_nonzero_var() {
        if (lower_varorder_pq_.empty()) return 0;

        uint16_t next = lower_varorder_pq_.top();
        lower_varorder_pq_.pop();
        return next;

        return 0;
    }

    void add_upper_var(uint16_t var) { upper_varorder_pq_.push(var); }

    void add_lower_var(uint16_t var) { lower_varorder_pq_.push(var); }


   private:
    std::vector<int32_t> table_elems_;
    std::vector<int32_t> var_heads_;
    std::vector<int32_t> num_elems_;
    std::vector<count_t> diff_counter_;
    std::vector<count_t> diff_counter_hi_;

    uint32_t entries_counter_;
    const uint16_t num_var_;
    int var_cache_;
    std::priority_queue<uint16_t, std::vector<uint16_t>, std::greater<uint16_t>>
        lower_varorder_pq_;
    std::priority_queue<uint16_t> upper_varorder_pq_;
};

#endif  // DP_MANAGER_H_