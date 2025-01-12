#include <pointer_based/blocks/Block.h>
#include <string_view>

Block::Block(Block *parent, int64_t start_index, int64_t end_index, const std::string &input) :
    parent_(parent),
    start_index_(start_index),
    end_index_(end_index),
    input_(input),
    left_(false),
    right_(false),
    first_block_(this),
    second_block_(nullptr),
    pointing_to_me_(0),
    level_index_(0),
    offset_(0),
    first_occurrence_level_index_(0),
    prefix_(std::string_view()),
    suffix_(std::string_view()) {}

Block::~Block() = default;

void Block::add_fast_substring_support(const int prefix_suffix_size) {
    const size_t len = length();

    // We aim to save the prefix and the suffix, each of size prefix_suffix_size.
    // If this block is smaller than their size, then we might as well just save the string represented by
    // this block
    const auto prefix_end = start_index_ + prefix_suffix_size > input_.length()
                                ? input_.end()
                                : input_.begin() + start_index_ + prefix_suffix_size;
    const auto suffix_end = end_index_ > input_.length() ? input_.end() : input_.begin() + end_index_ + 1;
    if (len < prefix_suffix_size) {
        prefix_ = std::string_view(input_.begin() + start_index_, suffix_end);
        suffix_ = std::string_view(input_.begin() + start_index_, suffix_end);
    } else {
        prefix_                 = std::string_view(input_.begin() + start_index_, prefix_end);
        const auto suffix_start = input_.begin() + end_index_ - prefix_suffix_size + 1;
        // Only set the suffix if its start is before its end (which is the usual case)
        // This happens if the suffix would start beyond the input string's end
        if (std::distance(suffix_start, suffix_end) >= 0) {
            suffix_ = std::string_view(suffix_start, suffix_end);
        }
    }

    for (Block *child : children_) {
        child->add_fast_substring_support(prefix_suffix_size);
    }
}

int Block::add_rank_select_support(int c) { return 0; }

int Block::rank(const int c, const int i) const { return 0; }

int Block::select(const int c, const int j) const { return -1; }

int64_t Block::length() const { return end_index_ - start_index_ + 1; }

std::string Block::represented_string() const { return input_.substr(start_index_, length()); }

std::vector<Block *> &Block::children(const int leaf_length, const int r) { return children_; }

void Block::clean_unnecessary_expansions() {}

bool Block::is_leaf() const { return true; }

int Block::access(const int i) const { return -1; }

char *Block::substr(char *buf, const int index, const int len) const { return buf; };

void Block::replace_child(Block *old_child, Block *new_child) {
    for (int i = 0; i < children_.size(); ++i) {
        if (children_[i] == old_child) {
            children_[i] = new_child;
            return;
        }
    }
}

std::string_view Block::prefix(const size_t count) const { return {prefix_.data(), std::min(count, prefix_.size())}; }

std::string_view Block::suffix(const size_t count) const {
    return {suffix_.begin() + suffix_.length() - std::min(count, suffix_.size()), suffix_.end()};
}
