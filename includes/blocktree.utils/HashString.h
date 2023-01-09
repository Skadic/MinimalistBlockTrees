#ifndef BLOCKTREE_HASHSTRING_H
#define BLOCKTREE_HASHSTRING_H

#include <string>

class HashString {
  public:
    /// The string's hash value
    const size_t hash_;
    /// The source string
    const std::string &s_;
    /// The string's start position in the source string
    const int start_;
    /// The string's end position in the source string (inclusive)
    const int end_;

    HashString(const size_t hash, const std::string &source, const int start_position, const int end_position);
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
