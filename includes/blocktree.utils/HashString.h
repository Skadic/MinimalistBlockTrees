#ifndef BLOCKTREE_HASHSTRING_H
#define BLOCKTREE_HASHSTRING_H

#include <string>

class HashString {
  public:
    /// The string's hash value
    const size_t hash_;
    /// The input string
    const std::string &s_;
    /// The string's start position in the input string
    const int start_;
    /// The string's end position in the input string (inclusive)
    const int end_;

    /// @brief Generate a new hash string that saves the hash of a given substring
    ///
    /// @param hash The hash value of the given substring
    /// @param input The original full input string
    /// @param start_position The start position of the hashed substring inside the input string
    /// @param end_position The (inclusive) end position of the hashed substring inside the input string
    HashString(size_t hash, const std::string &input, int start_position, int end_position);
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
