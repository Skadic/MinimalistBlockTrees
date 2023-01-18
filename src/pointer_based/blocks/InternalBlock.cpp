#include <algorithm>
#include <iostream>
#include <numeric>
#include <pointer_based/blocks/BackBlock.h>
#include <pointer_based/blocks/InternalBlock.h>
#include <pointer_based/blocks/LeafBlock.h>
#include <ranges>
#include <string_view>

InternalBlock::InternalBlock(Block *parent, int64_t start_index, int64_t end_index, const std::string &input) :
    Block(parent, start_index, end_index, input) {}

InternalBlock::~InternalBlock() {
    for (int i = children_.size() - 1; i >= 0; i--) {
        delete children_[i];
    }
}

int InternalBlock::add_rank_select_support(int character) {
    int cumulative_count = 0;
    // Add rank select support to each child and count the occurrences
    for (Block *child : children_) {
        cumulative_count += child->add_rank_select_support(character);
    }
    pop_counts_[character] = cumulative_count;
    return pop_counts_[character];
}

int InternalBlock::rank(const int character, const int i) const {
    // The cumulative length of the string up to the current child
    int cumulative_length = 0;
    // The number of times the character appeared up to the current child
    int cumulative_count = 0;
    for (Block *child : children_) {
        cumulative_length += child->length();
        if (i < cumulative_length) {
            // If the character we are looking for is inside this child, use rank directly on the child
            return cumulative_count + child->rank(character, i - (cumulative_length - child->length()));
        }
        // If it is not inside the child then read the pre-calculated rank from the child
        cumulative_count += child->pop_counts_[character];
    }
    return 0;
}

int InternalBlock::select(const int character, const int rank) const {
    // The cumulative length of the string up to the current child
    int cumulative_length = 0;
    // The number of times the character appeared up to the current child
    int cumulative_count = 0;
    for (auto it = children_.begin(); it != children_.end(); ++it) {
        if ((it + 1) == children_.end()) {
            // If this is the last child, we try to find it in that child
            return cumulative_length + (*it)->select(character, rank - cumulative_count);
        }
        if (cumulative_count + (*it)->pop_counts_[character] >= rank) {
            // if including the number of times 'character' appears is more than we want, then the character we are
            // looking for is in this child
            return cumulative_length + (*it)->select(character, rank - cumulative_count);
        }
        // Otherwise, update the variables and go to the next child
        cumulative_count += (*it)->pop_counts_[character];
        cumulative_length += (*it)->length();
    }
    return -1;
}

std::vector<Block *> &InternalBlock::children(const int leaf_length, const int arity) {
    // If the children are already calculated, return them
    if (!children_.empty()) {
        return children_;
    }

    // If the children are not calculated yet, then calculate them
    const int child_block_length = length() / arity;
    for (int i = 0; i < arity; ++i) {
        const int start = start_index_ + i * child_block_length;
        const int end   = start_index_ + (i + 1) * child_block_length - 1;
        if (start < input_.size()) {
            // If the child would be small enough to be a leaf or is smaller than its arity (so it wouldn't be able to
            // split anymore) create the children as leaf blocks
            Block *child = (child_block_length <= leaf_length || child_block_length <= arity)
                               ? static_cast<Block *>(new LeafBlock(this, start, end, input_))
                               : static_cast<Block *>(new InternalBlock(this, start, end, input_));
            children_.push_back(child);
        }
    }

    return children_;
}

void InternalBlock::clean_unnecessary_expansions() {
    for (auto it = children_.rbegin(); it != children_.rend(); ++it) {
        (*it)->clean_unnecessary_expansions();
    }

    // check whether all children are children
    bool all_children_leaves = std::ranges::all_of(children_, [&](Block *child) { return child->is_leaf(); });

    // If all children are leaves, no other block has this as its source, and this block points back to an earlier block
    // as its source, we can replace this entire subtree with a back block pointing to its source
    if (all_children_leaves && pointing_to_me_ == 0 && first_block_->start_index_ < start_index_ &&
        second_block_ != this) {
        BackBlock *bb = new BackBlock(parent_, start_index_, end_index_, input_, first_block_, second_block_, offset_);
        bb->level_index_                  = level_index_;
        bb->first_occurrence_level_index_ = first_occurrence_level_index_;
        bb->left_                         = true;
        bb->right_                        = true;
        parent_->replace_child(this, bb);
        delete this;
    } else { // To avoid dangling references
        first_block_  = this;
        second_block_ = nullptr;
    }
}

bool InternalBlock::is_leaf() const { return false; }

int InternalBlock::access(const int i) const {
    // The cumulative length of the string up to the current child
    int cumulative_length = 0;
    for (Block *child : children_) {
        // Scan through the children
        cumulative_length += child->length();
        // If we found the correct index, descend into the child
        if (i < cumulative_length) {
            return child->access(i - (cumulative_length - child->length()));
        }
    }
    return -1;
}

char *InternalBlock::substr(char *buf, const int index, const int len) const {

    const size_t child_len      = children_[0]->length();
    const size_t block_index    = index / child_len;
    const size_t internal_index = index - block_index * child_len;

    // If the substring is part of this block's prefix, we can just read it from there
    if (index + len <= prefix_.size()) {
        std::copy(prefix_.begin() + index, prefix_.begin() + index + len, buf);
        return buf + len;
    }

    // If the substring is part of this block's prefix, we can just read it from there
    if (index >= length() - suffix_.size()) {
        const size_t suffix_start_index = length() - suffix_.size();
        const size_t start_in_suffix    = index - suffix_start_index;
        std::copy(suffix_.begin() + start_in_suffix, suffix_.begin() + start_in_suffix + len, buf);
        return buf + len;
    }

    // Does the substring cross a child boundary?
    if (internal_index + len > child_len) {
        // the substring's prefix is part of this child's suffix
        const size_t prefix_size = child_len - internal_index;
        const auto   prefix      = children_[block_index]->suffix(prefix_size);

        // the substring's suffix is part of the following child's prefix
        const size_t suffix_size = len - prefix_size;
        const auto   suffix      = children_[block_index + 1]->prefix(suffix_size);

        std::copy(prefix.begin(), prefix.end(), buf);
        std::copy(suffix.begin(), suffix.end(), buf + prefix.size());

        return buf + prefix.size() + suffix.size();
    }

    return children_[block_index]->substr(buf, index - block_index * child_len, len);
};
