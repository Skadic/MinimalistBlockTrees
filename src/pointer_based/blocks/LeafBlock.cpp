#include <algorithm>
#include <iostream>
#include <pointer_based/blocks/LeafBlock.h>

LeafBlock::LeafBlock(Block *parent, int64_t start_index, int64_t end_index, const std::string &source) :
    Block(parent, start_index, end_index, source) {}

LeafBlock::~LeafBlock() {}

int64_t LeafBlock::size() const {
    int64_t source_end_index = source_.size() - 1;
    return (end_index_ <= source_end_index ? end_index_ : source_end_index) - start_index_ + 1;
}

int LeafBlock::add_rank_select_support(int c) {
    pop_counts_[c] = rank(c, size() - 1);
    return pop_counts_[c];
}

int LeafBlock::rank(const int c, const int i) const {
    int r = 0;
    for (int j = 0; j <= i; ++j) {
        if (source_[start_index_ + j] == c)
            ++r;
    }
    return r;
}

int LeafBlock::select(const int c, const int j) const {
    auto remaining = j;
    for (int i = 0; i < size(); ++i) {
        if (((int) (source_[start_index_ + i])) == c)
            --remaining;
        if (!remaining)
            return i;
    }
    return -1;
}

int LeafBlock::access(const int i) const { return source_[start_index_ + i]; }

char *LeafBlock::substr(char *buf, const int index, const int len) const {
    std::copy(source_.begin() + start_index_ + index, source_.begin() + start_index_ + index + len, buf);
    return buf + len;
}
