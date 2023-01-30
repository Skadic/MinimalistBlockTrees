#include "../../includes/blocktree.utils/RabinKarp.h"

RabinKarp::RabinKarp(const std::string &s,
                     const int          start_index,
                     const int          window_size,
                     const int          range,
                     const int          sigma) :
    sigma_(sigma),
    window_size_(window_size),
    s_(s),
    hash_(0),
    current_index_(start_index),
    rm_(1),
    range_(range) {
    // Calculate the initial hash value
    for (int i = start_index; i < start_index + window_size_; ++i) {
        hash_ = (sigma_ * hash_ + s_[i] + 128) % range_; // sigma or  little prime
    }

    for (int i = 0; i < window_size_ - 1; ++i) {
        rm_ = (rm_ * sigma_) % range_;
    }
}

uint64_t RabinKarp::hash() const { return hash_; }

void RabinKarp::next() {
    // Calculate the next hash value
    hash_ = (hash_ + range_ - rm_ * (s_[current_index_] + 128) % range_) % range_;
    current_index_++;
    hash_ = (hash_ * sigma_ + s_[current_index_ + window_size_ - 1] + 128) % range_;
}
