#include "pointer_based/BlockTree.h"

#include <iostream>
#include <pointer_based/blocks/InternalBlock.h>
#include <pointer_based/blocks/LeafBlock.h>

#include <stack>
#include <unordered_set>


BlockTree::BlockTree(std::string& input, int r, int root_block_arity, int leaf_length, bool clean, bool rs_support): 
  r_(r), root_arity_(root_block_arity), input_(input), leaf_length_(leaf_length), rank_select_support_(rs_support) {
    // If a leaf is supposed to be greater than the text, or if the arity is greater than the text is long, just make the entire tree one leaf node
    if (input_.size() <= leaf_length_ || input_.size()<r) {
        root_block_ = new LeafBlock(nullptr, 0, input_.size() - 1, input_);
    } else {
        const int number_of_leaves = (input_.size() % leaf_length_ == 0) ? input_.size() / leaf_length_ : input_.size() / leaf_length_ + 1;
        //The current height of the trea
        int height = 0;

        // The current max. index of nodes at the current height
        // So if there are n nodes at the current level, nl is n-1
        int nl = number_of_leaves-1;
        int64_t block_length = leaf_length_;

        // In total this calculates 
        while (nl > 0){
            height++;
            block_length*=r_;
            nl/=r_;
        }

        root_block_ = new InternalBlock(nullptr, 0, block_length - 1, input_);
    }
    
    std::unordered_set<int> characters;
    if (rank_select_support_) {
      for (char c: input) {
          characters.insert(c);
      }
    }

    if (clean) {
      this->process_back_pointers();
      this->clean_unnecessary_expansions();
    }

    for (char c: characters) {
        this->add_rank_select_support(c);
    }
}


BlockTree::~BlockTree() {
    delete root_block_;
}


void BlockTree::add_rank_select_support(int c) {
    root_block_->add_rank_select_support(c);
}


int BlockTree::rank(int c, int i) {
    return root_block_->rank(c, i);
}


int BlockTree::select(int c, int j) {
    return root_block_->select(c, j);
}


std::vector<std::vector<Block*>> BlockTree::levelwise_iterator() {
    std::vector<std::vector<Block*>> result = {{root_block_}};
    while (!dynamic_cast<LeafBlock*>(result.back()[0])) {
        std::vector<Block*> next_level = {};
        for (Block *b : result.back()) {
            const auto block_arity = b == root_block_ ? root_arity_ : r_; 
            for (Block *child : b->children(leaf_length_, block_arity)) {
                next_level.push_back(child);
            }
        }
        result.push_back(next_level);
    }

    return result;
}


void BlockTree::clean_unnecessary_expansions() {
    root_block_->clean_unnecessary_expansions();
    for (std::vector<Block*> level : levelwise_iterator()) {
        for (int i = 0; i<level.size(); ++i) {
            level[i]->level_index_ = i;
            level[i]->first_occurrence_level_index_ = level[i]->first_block_->level_index_;
        }
    }
}


int BlockTree::access(int i) {
    return root_block_->access(i);
}


/**
 * @brief This takes a vector of blocks and returns a vector of their children in order.
 *
 * As this method's name says, this is supposed to be used with input vectors 
 * which consist of blocks of a single level in the tree from left to right.
 *
 * @param level A vector of pointers to blocks whose children to return.
 */
std::vector<Block*> BlockTree::next_level(std::vector<Block*>& level) {

    // The arity of the current level. If we're at the root, then the arity is of course the root arity.
    const auto block_arity = (level.size() == 1 && level[0] == root_block_) ? root_arity_ : r_;

    std::vector<Block*> next_level;
    for (int i = 0; i < level.size(); ++i) {
        Block* b = level[i];
        for (Block *child : b->children(leaf_length_, block_arity)) { // Do it in order
            child->level_index_ = next_level.size();
            child->first_occurrence_level_index_ = next_level.size();
            next_level.push_back(child);
        }
    }
    return next_level;
}


void BlockTree::forward_pair_window_block_scan(std::vector<Block*>& level, int pair_window_size, int N, std::unordered_map<HashString, std::vector<std::pair<Block*, Block*>>>& pair_hashtable) {
    for (std::vector<Block *>::iterator it = level.begin(); it != level.end();) {
        Block *b = (*it);
        b->right_ = true;
        int offset = 0;
        RabinKarp rk(input_, (*it)->start_index_ + offset, pair_window_size, N); // offset is always 0 here
        for (; it != level.end() && ((*it) == b || (*(it-1))->end_index_ == (*it)->start_index_ - 1); it++) {
            Block* current = *(it);
            bool last_block = ((it+1) == level.end() ||  current->end_index_ != (*(it+1))->start_index_ - 1);
            for (offset = 0; offset < current->length(); ++offset) {
                if (last_block && current->length() - offset < pair_window_size)  break;
                HashString hS(rk.hash(), input_, current->start_index_ + offset, current->start_index_ + offset + pair_window_size - 1);
                std::unordered_map<HashString, std::vector<std::pair<Block*, Block*>>>::const_iterator result = pair_hashtable.find(hS);
                if (result != pair_hashtable.end()) { // Here, It could be that the scanning should have finished with the penultimate, but it never should enter this ''if''
                                                        // when We're on the penultimate block and the window exceeds the last block because if that is a first occurrence should have been occured before in a pair of blocks
                                                        // maybe use a condition more like rk's condition below could work fine too
                                                        // Same logic: for when passing a window of size 2l + 2 over 2 block of length l
                    for (std::pair<Block*,Block*> p: result->second) {
                        if (current->start_index_ + offset < p.first->start_index_) {
                            p.first->left_ = true;
                            p.second->right_ = true;
                        }
                    }
                    pair_hashtable.erase(hS);
                }
                if (current->start_index_+offset+pair_window_size < input_.size()) rk.next();
            }
        }
        (*(it-1))->left_ = true;
    }
}


void BlockTree::forward_window_block_scan(std::vector<Block*>& level, int window_size, int N, std::unordered_map<HashString, std::vector<Block*>>& hashtable) {
    int i = 0;
    for (std::vector<Block *>::iterator it = level.begin(); it != level.end();) {
        Block *b = (*it);
        int offset = 0;
        RabinKarp rk(input_, (*it)->start_index_ + offset, window_size, N);
        for (; it != level.end() && ((*it) == b || (*(it-1))->end_index_ == (*it)->start_index_ - 1); it++, i++) {
            Block* current = *(it);
            bool last_block = ((it+1) == level.end() ||  current->end_index_ != (*(it+1))->start_index_ - 1);
            for (offset = 0; offset < current->length(); ++offset) {
                if (last_block && current->length() - offset < window_size)  break;
                HashString hS(rk.hash(), input_, current->start_index_ + offset, current->start_index_ + offset + window_size - 1);
                std::unordered_map<HashString, std::vector<Block *>>::const_iterator result = hashtable.find(hS);
                if (result != hashtable.end()) {
                    std::vector<Block*> blocks = result->second;
                    for (Block* b : blocks) {
                        b->first_occurrence_level_index_ = i;
                        b->first_block_ = current;
                        b->offset_ = offset;
                        if (offset + window_size > b->first_block_->length()) b->second_block_ = (*(it+1));
                        else b->second_block_ = nullptr;
                    }
                    hashtable.erase(hS);
                }
                if (current->start_index_+offset+window_size < input_.size()) rk.next();
            }
        }
    }
}


void BlockTree::block_scan(std::vector<Block *>& level, int N, std::unordered_map<HashString, std::vector<Block*>>& hashtable) {
    for (Block* b : level) {
        RabinKarp rk(input_, b->start_index_, b->length(), N);
        HashString hS(rk.hash(), input_, b->start_index_, b->end_index_);

        std::unordered_map<HashString, std::vector<Block*>>::const_iterator result = hashtable.find(hS);

        if (result == hashtable.end())
            hashtable[hS] = {b};
        else
            hashtable[hS].push_back(b);
    }
}


void BlockTree::process_level(std::vector<Block*>& level) {
  
    // Large prime
    const int N = 6700417;
    // The length of blocks on this level
    int level_length = level.front()->length();

    // Block scan
    std::unordered_map<HashString, std::vector<Block*>> hashtable;
    block_scan(level, N, hashtable);

    // Pairs of blocks scan
    std::unordered_map<HashString, std::vector<std::pair<Block *,Block*>>> pair_hashtable;
    for (std::vector<Block *>::iterator it = level.begin(); it != level.end();) {
        for (++it; (it != level.end() && (*(it-1))->end_index_ == (*it)->start_index_ - 1); ++it) {
            Block* current = (*(it - 1));
            Block* next = (*it);
            RabinKarp rk(input_, current->start_index_, current->length() + next->length(), N);
            HashString hS(rk.hash(), input_, current->start_index_, current->start_index_ + current->length() + next->length()-1); // Second parameter is next->end_index

            std::unordered_map<HashString, std::vector<std::pair<Block *,Block*>>>::const_iterator result = pair_hashtable.find(hS);

            if (result == pair_hashtable.end())
                pair_hashtable[hS] = {{current, next}};
            else
                pair_hashtable[hS].push_back({current, next});
        }
    }


    // Window block scan
    //Establishes first occurrences of blocks
    forward_window_block_scan(level, level_length, N, hashtable);



    // Window Pair of blocks scans
    if (level.size() > 1)
        forward_pair_window_block_scan(level, level_length*2, N, pair_hashtable);




    // BackBlock creation
    for (int i = 0; i < level.size(); ++i) {
        Block* b = level[i];
        if (b->left_ && b->right_ && b->first_occurrence_level_index_ < b->level_index_) {
            // This doesn't have the bug of the dangling reference fixed with first_occurrence_level_index, because it shouldn't happen that
            // A block points back to a BackBlock
            BackBlock* bb = new BackBlock(b->parent_, b->start_index_, b->end_index_, input_,
                                          level[b->first_occurrence_level_index_], (b->second_block_ ==
                                                                                                      nullptr) ? nullptr : level[b->first_occurrence_level_index_ +1], b->offset_);
            bb->level_index_ = b->level_index_;
            bb->first_occurrence_level_index_ = b->first_occurrence_level_index_;
            bb->left_ = true;
            bb->right_ = true;
            b->parent_->replace_child(b, bb);
            delete b;
            level[i] = bb;
        }
    }

}


/**
 * @brief Actually build the block tree and insert back-links for repeated blocks.
 */
void BlockTree::process_back_pointers() {
    std::vector<Block*> current_level = {root_block_};
    std::stack<Block*> none_blocks;
    // While there is a new level in the tree to process 
    while ((current_level = next_level(current_level)).size() != 0) {
        //const auto block_arity = current_level[0]->parent_ == root_block_ ? root_block_arity_ : r_; 
        // If the nodes at this levels are at leaf-size,
        // we stop, as there are no back-pointers that can exist on this level
        if (current_level[0]->length() < r_ ||  current_level[0]->length() <= leaf_length_) break;
        // If we have at least 1 block to process on this level and the last node
        // exceeds the length of the text (thus the block is not "full") we take such blocks out before processing
        // as they cannot exist anywhere further to the left on this level (all nodes to the left would be full after all)
        while (current_level.size() != 0 && current_level.back()->end_index_ >= input_.size()) {
            none_blocks.push(current_level.back());
            current_level.pop_back();
        }
        // Process the backlinks on this level
        process_level(current_level);
        // Put the blocks back that we took out before
        while (!none_blocks.empty()) {
            current_level.push_back(none_blocks.top());
            none_blocks.pop();
        }
    }
}


void BlockTree::process_level_heuristic(std::vector<Block*>& level) {

    int N = 6700417; // Large prime
    int level_length = level.front()->length();

    // Block scan
    std::unordered_map<HashString, std::vector<Block*>> hashtable;
    block_scan(level, N, hashtable);

    // Window block scan
    // This is almost the same as forward_window_block_scan, as well as the BackBlock creation
    forward_window_block_scan(level, level_length, N, hashtable);



    // BackBlock creation
    for (int i = 0; i < level.size(); ++i) {
        Block* b = level[i];
        if (b->first_occurrence_level_index_ < b->level_index_) {

            BackBlock* bb = new BackBlock(b->parent_, b->start_index_, b->end_index_, input_,
                                          level[b->first_occurrence_level_index_], (b->second_block_ ==
                                                                                                      nullptr) ? nullptr : level[b->first_occurrence_level_index_ +1], b->offset_);
            bb->level_index_ = b->level_index_;
            bb->first_occurrence_level_index_ = b->first_occurrence_level_index_;
            b->parent_->replace_child(b, bb);
            delete b;
            level[i] = bb;
        }
    }

}


void BlockTree::process_back_pointers_heuristic() {
    std::vector<Block *> current_level = {root_block_};
    std::stack<Block*> none_blocks;
    while ((current_level = next_level(current_level)).size() != 0) {
        if (current_level[0]->length() < r_ ||
            current_level[0]->length() <= leaf_length_)
            break;
        while (current_level.size() != 0 && current_level.back()->end_index_ >= input_.size()) {
            none_blocks.push(current_level.back());
            current_level.pop_back();
        }
        process_level_heuristic(current_level);
        while (!none_blocks.empty()) {
            current_level.push_back(none_blocks.top());
            none_blocks.pop();
        }
    }
}

void BlockTree::print() {
  auto it = levelwise_iterator();
  for (auto &level : it) {
    for (Block *b : level) {
      std::cout << "[" << b->start_index_ << ", " << b->end_index_ << "], ";
    }
    std::cout << std::endl;
  }
}
