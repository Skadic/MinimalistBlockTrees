#include <blocktree.utils/HashString.h>

HashString::HashString(const size_t hash, const std::string &input, const int init, const int end) :
    hash_(hash),
    s_(input),
    start_(init),
    end_(end) {}

HashString::~HashString() {}

bool HashString::operator==(const HashString &other) const {
    int length = end_ - start_ + 1;
    if (length != other.end_ - other.start_ + 1) {
        return false;
    }
    if (s_.size() > 0) {
        for (int i = 0; i < length; ++i) {
            if (s_[start_ + i] != other.s_[other.start_ + i]) {
                return false;
            }
        }
        return true;
    }
    return true;
}
