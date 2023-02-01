# MinimalistBlockTrees

This repository contains an fork of [Manuel Ariel Cáceres Reyes'](https://github.com/elarielcl)
implementation for the BlockTree data structure described [here](https://ieeexplore.ieee.org/document/7149265).
Details of the implementation follow the experimental studies
conducted in Reyes' MSc. thesis available [here](https://users.dcc.uchile.cl/~gnavarro/mem/algoritmos/tesisManuel.pdf).

This fork is a refactor of the original code, with more documentation,
altered variable/function naming and support for faster substring queries
on uncompressed and compressed block trees.
The bitmap-optimized block tree does not support substring queries at this point.

## Installation Guide

First clone the repo:

```bash
 git clone https://github.com/elarielcl/MinimalistBlockTrees.git
```

This project is a CMake project.
To build this project with some runnables you should do

```bash
cd MinimalistBlockTrees
mkdir build
cd build
cmake ..
cmake .. # Issue: second cmake necessary to compile sdsl
make
```

You can add an executable by writing your file in the `executables` directory
and add its name to the `executables/CMakeLists.txt` file,
this adds the necessary libraries for you:

```cpp
set(project_EXECUTABLES
        <new_executable>
        build_bt
        read_bt
        bt_repl)
...
```

### Commitlint

To use the commitlint config, make sure to run

```bash
npm i
```

## Usage Example

Let's suppose we want to build a BlockTree, so we do:

```cpp
...
std::string input = "AACCCTGCTGCTGATCGGATCGTAGC";
int arity = 2; // The arity of the BlockTree
int root_arity = 8; // The arity of the BlockTree's root block
int leaf_len = 16; // The length of the strings represented by the leaves of the BlockTree
bool process = false; // Whether the BlockTree should be immediately constructed (default=false)
bool rs_support = false; // Whether to build the BlockTree with rank/select support (default=false)

BlockTree* bt = new BlockTree(input, arity, root_arity, leaf_len, process, rs_support); // This creates the BlockTree object, a pointer-based implementation
bt->process_back_pointers(); // This method constructs the BackPointers in the BlockTree. Unncessary if process == true
bt->clean_unnecessary_expansions(); // This method removes the expansion of InternalBlocks that are unnecesary (this is also called pruning). Unncessary if process == true
 ...
```

In case you want to build the a more compact version
but without theoretical guarantee, you can replace the last two lines with

```cpp
bt->process_back_pointers_heuristic();
```

.. now you have a proper BlockTree answering access queries(`bt->access(i)`),
if you want to give it `bt->rank(c,i) & bt->select(c,j)` support you should do:

```cpp
for (int c: characters)
    bt->add_rank_select_support(c);
```

This is only necessary if `rs_support` was set to `false` before.

The above is a pointer-based implementation of BlockTrees,
if you want to have a more compact representation
(using bitvectors to represent the tree) you should do:

```cpp
...
CBlockTree *cbt = new CBlockTree(bt); // Builds a more compact BlockTree representation
cbt->access(i);
cbt->select(c, i); // Rank/Select equires the source block tree to support rank/select!
cbt->size(); // Obtains the approximate size (in bytes)
cbt->serialize(out); // Serializes the tree to an output stream
...
CBlockTree *loaded_cbt = new CBlockTree(in); // Loads a CBlockTree from an input stream
...
```

If you'd like to build a bitmap over your input (or your input has a binary alphabet)
you can build an (even) more compact representation specialized for bitmaps.
Just supply it with a block tree, and the input symbol that should be considered to be a one.

```cpp
...
// E.g. for a balanced parantheses sequence, the opening parenthesis is considered to be a 1 here
CBitBlockTree *cbbt = new CBitBlockTree(bt, '(');

cbbt->access(i);
cbbt->rank_1(i);
cbbt->select_0(i);
...
...
```

... and never forget to delete your trash ;)

```cpp
...
delete bt;
delete cbt;
delete loaded_cbt;
delete cbbt;
...
```
