#ifndef HIDDEN_NODE_STACK_H_
#define HIDDEN_NODE_STACK_H_

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <stack>
#include <vector>
/**
 * A stack storing the node cells to hide.
 * this class helps to restore hidden nodes in an appropriate order.
 */

class HiddenNodeStack {
   public:
    enum class HideType { UpperZero, LowerZero, CoverUp, CoverDown };
   
    using stack_value_t = std::pair<int32_t, HideType>;


    HiddenNodeStack();

    stack_value_t top() const { return *node_stack_.crbegin(); }

    void pop() { node_stack_.pop_back(); }
    void push_cover_down(const int32_t val) {
        node_stack_.emplace_back(val, HideType::CoverDown);
    }

    void push_cover_up(const int32_t val) {
        node_stack_.emplace_back(val, HideType::CoverUp);
    }

    void push_upperzero(const int32_t val) {
        node_stack_.emplace_back(val, HideType::UpperZero);
    }

    void push_lowerzero(const int32_t val) {
        node_stack_.emplace_back(val, HideType::LowerZero);
    }

    bool is_empty() const {
        return node_stack_.size() == stack_start_positions_.top();
    }
    void push_checkpoint() { stack_start_positions_.push(node_stack_.size()); }

    void reverse_current_stack() {
        std::reverse(node_stack_.begin() + stack_start_positions_.top(),
                     node_stack_.end());
    }

    std::vector<stack_value_t>::const_iterator stack_cbegin() const {
        return node_stack_.cbegin() + stack_start_positions_.top();
    }

    std::vector<stack_value_t>::const_iterator stack_cend() const {
        return node_stack_.cend();
    }

    void pop_checkpoint() {
        assert(is_empty());
        stack_start_positions_.pop();
    }

   private:
    std::vector<stack_value_t> node_stack_;
    std::stack<int32_t, std::vector<int32_t>> stack_start_positions_;
};

#endif  // HIDDEN_NODE_STACK_H_