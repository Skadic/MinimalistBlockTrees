#ifndef BLOCKTREE_PBLOCK_H
#define BLOCKTREE_PBLOCK_H

#include <string>
#include <vector>
#include <unordered_map>

class Block {
public:
    /**
     * @parent_ The parent block
     */
    Block* parent_;

    /**
     * @start_index_ The index in the input text at which this block starts
     */
    int64_t start_index_;


    /**
     * @end_index_ The index in the input text at which this block starts
     */
    int64_t end_index_;

    /**
     * @source_ The string this block tree is built from
     */
    std::string& source_;

    std::unordered_map<int,int> ranks_;
    std::unordered_map<int,int> second_ranks_;

    Block* first_block_;
    Block* second_block_;
    int offset_;
    bool left_;
    bool right_;
    int pointing_to_me_;
    int level_index_;
    int first_occurrence_level_index_;

    std::vector<Block*> children_;

    Block(Block* parent, int64_t start_index, int64_t end_index, std::string &source);
    virtual ~Block();

    int64_t length();
    std::string represented_string();

    virtual int add_rank_select_support(int character);

    virtual int access(int i);
    virtual int rank(int character, int i);
    virtual int select(int character, int rank);

    virtual std::vector<Block*>& children(int leaf_length, int arity);
    virtual void clean_unnecessary_expansions();
    void replace_child(Block* old_child, Block* new_child);

    virtual bool is_leaf();
};

#endif //BLOCKTREE_PBLOCK_H
