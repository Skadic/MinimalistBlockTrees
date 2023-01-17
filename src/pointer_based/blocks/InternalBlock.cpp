#include <algorithm>
#include <iostream>
#include <numeric>
#include <pointer_based/blocks/BackBlock.h>
#include <pointer_based/blocks/InternalBlock.h>
#include <pointer_based/blocks/LeafBlock.h>
#include <ranges>
#include <string_view>

InternalBlock::InternalBlock(Block *parent, int64_t start_index, int64_t end_index, const std::string &source) :
    Block(parent, start_index, end_index, source) {}

InternalBlock::~InternalBlock() {
    for (int i = children_.size() - 1; i >= 0; i--) {
        delete children_[i];
    }
}

void InternalBlock::add_fast_substring_support(const int prefix_suffix_size) {
    const size_t len = length();

    // We aim to save the prefix and the suffix, each of size prefix_suffix_size.
    // If this block is smaller than both of these together, then we might as well just save the string represented by
    // this block, as it would be smaller
    auto prefix_end = start_index_ + prefix_suffix_size > source_.length() ? source_.end() : source_.begin() + start_index_ + prefix_suffix_size;
    auto suffix_end = end_index_ > source_.length() ? source_.end() : source_.begin() + end_index_ + 1;
    if (len <= 2 * prefix_suffix_size) {
        prefix_ = std::string_view(source_.begin() + start_index_, suffix_end);
        suffix_ = std::string_view(source_.begin() + start_index_, suffix_end);
        return;
    }

    // Insert the prefix
    prefix_ = std::string_view(source_.begin() + start_index_, prefix_end);

    // Insert the suffix
    suffix_ = std::string_view(source_.begin() + end_index_ - prefix_suffix_size + 1, suffix_end);

    for (Block *child : children_) {
        child->add_fast_substring_support(prefix_suffix_size);
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
    if (children_.size() > 0) {
        return children_;
    }

    // If the children are not calculated yet, then calculate them
    const int child_block_length = length() / arity;
    for (int i = 0; i < arity; ++i) {
        const int start = start_index_ + i * child_block_length;
        const int end   = start_index_ + (i + 1) * child_block_length - 1;
        if (start < source_.size()) {
            // If the child would be small enough to be a leaf or is smaller than its arity (so it wouldn't be able to
            // split anymore) create the children as leaf blocks
            Block *child = (child_block_length <= leaf_length || child_block_length <= arity)
                               ? static_cast<Block *>(new LeafBlock(this, start, end, source_))
                               : static_cast<Block *>(new InternalBlock(this, start, end, source_));
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
        BackBlock *bb = new BackBlock(parent_, start_index_, end_index_, source_, first_block_, second_block_, offset_);
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
