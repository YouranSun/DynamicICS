// The @Treap class implements the treap structure.
#ifndef CORE_GADGET_TREAP2_H_
#define CORE_GADGET_TREAP2_H_

#include <vector>
#include <random>
#include <chrono>
#include "Tools.h"

namespace TREAP2 {

struct Data final {
  Node v; // the vertex id and weight
  bool isKey; // is key or not
  // int sumKey; // total number of keys in the subtree
  int mnKey; // min index of a key node in the subtree
  int mxKey; // max index of a key node in the subtree
} ;

class Treap final {
 public:
  explicit Treap();
  explicit Treap(const int n);

  int TreeNode(const int x, const Node v, const bool k);
  int NewTreeNode(const Node v, const bool k);
  void Insert(const int x, const bool f, int& r);
  void InsertAfter(const int x, const int y, int& r);
  void InsertBefore(const int x, const int y, int& r);
  void Delete(const int x, int& r);
  int Merge(const int r1, const int r2);
  int Rank(const int x) const;
  int Select(const int r, const int rank) const;
  int Root(const int x) const;
  int Minimum(const int x) const;
  int Maximum(const int x) const;
  int Size(const int r) const;
  void Check(const int r) const;
  void Traverse(const int r) const;
  void TraverseToVec(const int r, std::vector<int> &tmp) const;
  void Print() const;
  void Print(const int x) const;
  int FindNextKey(const Node x, const int r) const; // by weight
  int FindPrevKey(const int x, const int r) const; // by rank
  int FindPrev(int x) const;
  Node Get(const int x) const;
  void QueryCV(const int x, int &r, std::vector<std::pair<int, int> > &V) const;


 private:
  struct TreapNode final {
    int p;  // the parent
    int l;  // the left child
    int r;  // the right child
    int s;  // the size of the subtree rooted at this node
    int w;  // the priority
    Data d; // the data
  };

  void LeftRotate(int& x);
  void RightRotate(int& x);
  void SubCheck(const int x) const;
  void PushUp(const int x);
  int n_;
  std::vector<TreapNode> nd_;
  std::mt19937 rnd;
};
}

#endif