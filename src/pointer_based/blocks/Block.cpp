#include <pointer_based/blocks/Block.h>

Block::Block(Block *parent, int64_t start_index, int64_t end_index, const std::string &source) :
    parent_(parent),
    start_index_(start_index),
    end_index_(end_index),
    source_(source),
    left_(false),
    right_(false),
    first_block_(this),
    second_block_(nullptr),
    pointing_to_me_(0),
    level_index_(0),
    first_occurrence_level_index_(0) {}

Block::~Block() {}

void Block::add_fast_substring_support(int prefix_suffix_size) {}

int Block::add_rank_select_support(int c) { return 0; }

int Block::rank(const int c, const int i) const { return 0; }

int Block::select(const int c, const int j) const { return -1; }

int64_t Block::length() const { return end_index_ - start_index_ + 1; }

std::string Block::represented_string() const { return source_.substr(start_index_, length()); }

std::vector<Block *> &Block::children(const int leaf_length, const int r) { return children_; }

void Block::clean_unnecessary_expansions() {}

bool Block::is_leaf() const { return true; }

int Block::access(const int i) const { return -1; }

void Block::replace_child(Block *old_child, Block *new_child) {
    for (int i = 0; i < children_.size(); ++i) {
        if (children_[i] == old_child) {
            children_[i] = new_child;
            return;
        }
    }
}
