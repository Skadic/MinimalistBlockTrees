#ifndef BLOCKTREE_HASHSTRING_H
#define BLOCKTREE_HASHSTRING_H

#include <string>

class HashString {
  public:
    /// The string's hash value
    size_t hash_;
    /// The source string
    std::string &s_;
    /// The string's start position in the source string
    int start_;
    /// The string's end position in the source string (inclusive)
    int end_;

    HashString(size_t hash, std::string &source, int start_position, int end_position);
    ~HashString();

    bool operator==(const HashString &other) const;
};

namespace std {
template<>
struct hash<HashString> {
    ///
    /// @brief Return this string's hash value
    ///
    /// @param hS The hash string whose hash to get.
    ///
    std::size_t operator()(const HashString &hS) const { return hS.hash_; }
};
} // namespace std
#endif // BLOCKTREE_HASHSTRING_H
