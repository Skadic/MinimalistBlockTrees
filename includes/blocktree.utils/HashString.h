#ifndef BLOCKTREE_HASHSTRING_H
#define BLOCKTREE_HASHSTRING_H

#include <string>

class HashString {
public:
    size_t hash_;
    std::string& s_;
    int init_;
    int end_;

    HashString(size_t hash, std::string& source, int init, int end);
    ~HashString();

    bool operator==(const HashString& other) const;
};

namespace std {
    template <> struct hash<HashString> {
        std::size_t operator()(const HashString& hS) const {
            return hS.hash_;
        }
    };
}
#endif //BLOCKTREE_HASHSTRING_H
