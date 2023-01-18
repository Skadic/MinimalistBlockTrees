#include <iostream>
#include <pointer_based/blocks/BackBlock.h>

BackBlock::BackBlock(Block             *parent,
                     int64_t            start_index,
                     int64_t            end_index,
                     const std::string &input,
                     Block             *first_block,
                     Block             *second_block,
                     int                offset) :
    Block(parent, start_index, end_index, input) {
    first_block_ = first_block;

    // If there is a second block from which we need to copy
    if (second_block != nullptr) {
        // If the second block's range is the same as this block's then the second block is this very block
        if (second_block->start_index_ == start_index && second_block->end_index_ == end_index) {
            second_block_ = this;
        } else {
            second_block_ = second_block;
        }
    }
    offset_ = offset;
    // "Notify" the other blocks that they are now being pointed to in case they exist
    if (first_block_ != nullptr) {
        first_block_->pointing_to_me_++;
    }
    if (second_block_ != nullptr) {
        second_block_->pointing_to_me_++;
    }
}

BackBlock::~BackBlock() {
    if (first_block_ != nullptr) {
        first_block_->pointing_to_me_ = first_block_->pointing_to_me_ - 1;
    }
    if (second_block_ != nullptr) {
        second_block_->pointing_to_me_ = second_block_->pointing_to_me_ - 1;
    }
}

int BackBlock::add_rank_select_support(const int c) {
    // We want to calculate the number of times c appears in this block and store it into ranks_
    // In second_ranks_ we store the amount of times c appears in the part of this block's source which lies in the
    // first block.

    // This is the number of times c appears before this block's source
    int first_rank = first_block_->rank(c, offset_ - 1);
    // If there is no second block to copy from then this block's source is entirely contained in first block.
    // So, we calculate the rank until the end of the source and subtract it from the rank before the start of the
    // source and so, we get the number of times this character appears in the first block of the source
    // Of there is a second block, then we can only rank until the end of the first block. We do that and subtract the
    // rank before the source
    int second_rank = (second_block_ == nullptr) ? first_block_->rank(c, offset_ + length() - 1) - first_rank
                                                 : first_block_->rank(c, first_block_->length() - 1) - first_rank;
    // Insert the number of times c appears in the part of the source that lies in first_block_
    pop_counts_in_first_block_[c] = second_rank;
    // Insert the number of times c appears in this block in total
    // If there is no second block, then second_rank already contains the answer
    // If there is, the ne need to add the number of times c appears in the second block
    pop_counts_[c] = (second_block_ == nullptr)
                         ? second_rank
                         : second_rank + second_block_->rank(c, offset_ + length() - 1 - first_block_->length());
    return pop_counts_[c];
}

int BackBlock::rank(const int c, const int i) const {
    if (i + offset_ >= first_block_->length()) {
        return pop_counts_in_first_block_.at(c) +
               second_block_->rank(c, offset_ + i - first_block_->length()); // Loop if it's itself
    }
    return first_block_->rank(c, i + offset_) - (first_block_->pop_counts_[c] - pop_counts_in_first_block_.at(c));
}

int BackBlock::select(const int c, const int j) const {
    if (j > pop_counts_in_first_block_.at(c)) {
        return second_block_->select(c, j - pop_counts_in_first_block_.at(c)) + first_block_->length() - offset_;
    }
    return first_block_->select(c, j + first_block_->pop_counts_[c] - pop_counts_in_first_block_.at(c)) - offset_;
}

int BackBlock::access(const int i) const {
    if (i + offset_ >= first_block_->length()) {
        return second_block_->access(offset_ + i - first_block_->length());
    }
    return first_block_->access(i + offset_);
}

char *BackBlock::substr(char *buf, const int index, const int len) const {
    const size_t block_index           = index / length();
    const size_t internal_index        = index - block_index * length();
    const size_t source_internal_index = internal_index + offset_;

    // Is the substring entirely contained in the first source block?
    if (source_internal_index + len <= length()) {
        const size_t first_block_internal_index = source_internal_index;
        return first_block_->substr(buf, first_block_internal_index, len);
    }

    // Is the substring entirely contained in the second source block?
    if (internal_index + offset_ >= length()) {
        const size_t second_block_internal_index = source_internal_index - length();
        return second_block_->substr(buf, second_block_internal_index, len);
    }

    // This means, the substring crosses the boundary between the first and second source blocks
    const size_t prefix_size = length() - offset_ - index;
    const size_t suffix_size = len - prefix_size;

    // the substring's prefix is part of the first source block
    const auto prefix = first_block_->suffix(prefix_size);
    // the substring's suffix is part of the second source block
    const auto suffix = second_block_->prefix(suffix_size);

    std::copy(prefix.begin(), prefix.end(), buf);
    std::copy(suffix.begin(), suffix.end(), buf + prefix.size());

    return buf + prefix.size() + suffix.size();
}
