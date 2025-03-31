#include "treap2.h"
#include "Tools.h"

#include <random>
#include <assert.h>

#include "defs.h"

namespace TREAP2 {
Treap::Treap() {
  n_ = 0;
  nd_.clear();
  rnd = std::mt19937(2226701);
  NewTreeNode(Node(-1, 2147483647), 0);
  NewTreeNode(Node(-1, 2147483647), 0);
  nd_[0].l = nd_[0].r = nd_[0].s = 0, nd_[0].w = -2147483648;
  nd_[1].w = 2147483647;
  return;
}
Treap::Treap(const int n) {
  ASSERT(0 <= (n_ = n));
  nd_ = std::vector<TreapNode>(n_ + 2);
  // generate random priorities
  std::default_random_engine generator;
  std::uniform_int_distribution<int> distribution(-n_ + 1, n_ - 1);
  for (int i = 0; i < n_; ++i) {
    nd_[i].w = distribution(generator);
  }
  nd_[0].l = nd_[0].r = 0;
  nd_[0].s = 0;
  nd_[0].w = -n_;
  nd_[1].w = n_;
}
int Treap::TreeNode(const int x, const Node v, const bool k) {
  nd_[x].l = nd_[x].r = nd_[x].p = 0;
  nd_[x].d = Data{v, k, k ? x : 0, k ? x : 0};
  return x;
}
int Treap::NewTreeNode(const Node v, const bool k) {
  TreapNode tmp;
  tmp.l = tmp.r = 0;
  tmp.d = Data{v, k, k ? n_ : 0, k ? n_ : 0}, tmp.w = rnd();
  tmp.s = 1, tmp.p = 0;
  nd_.push_back(tmp);
  return n_++;
}
void Treap::PushUp(const int x) {
  int p = nd_[x].l, q = nd_[x].r;
  nd_[x].d.mnKey = nd_[x].d.mxKey = 0;
  if (nd_[p].d.mnKey) nd_[x].d.mnKey = nd_[p].d.mnKey;
  else if (nd_[x].d.isKey) nd_[x].d.mnKey = x;
  else if (nd_[q].d.mnKey) nd_[x].d.mnKey = nd_[q].d.mnKey;
  if (nd_[q].d.mxKey) nd_[x].d.mxKey = nd_[q].d.mxKey;
  else if (nd_[x].d.isKey) nd_[x].d.mxKey = x;
  else if (nd_[p].d.mxKey) nd_[x].d.mxKey = nd_[p].d.mxKey;
  ASSERT(nd_[x].d.mxKey == 0 || nd_[nd_[x].d.mxKey].d.isKey);
}
void Treap::Insert(const int x, const bool f, int& r) {
  if (0 == r) {
    nd_[x].l = nd_[x].r = 0;
    nd_[x].s = 1;
    r = x;
  } else {
    ++nd_[r].s;
    int& c = f ? nd_[r].l : nd_[r].r;
    Insert(x, f, c);
    PushUp(r);
    nd_[c].p = r;
    // the heap property may violate
    if (nd_[c].w > nd_[r].w) {
      if (f) {
        RightRotate(r);
      } else {
        LeftRotate(r);
      }
    }
  }
  nd_[r].p = 0;
}
void Treap::InsertAfter(const int x, const int y, int& r) {
  if (y == 0) {
    Insert(x, true, r);
    return;
  }
  Insert(x, true, nd_[y].r);
  nd_[nd_[y].r].p = y;
  int p = y;
  int c = nd_[y].r;
  while (0 != p) {
    ++nd_[p].s;
    int* gp = &r;
    if (0 != nd_[p].p) {
      if (p == nd_[nd_[p].p].l) {
        gp = &nd_[nd_[p].p].l;
      } else {
        gp = &nd_[nd_[p].p].r;
      }
    }
    if (nd_[c].w > nd_[p].w) {
      if (c == nd_[p].l) {
        RightRotate(*gp);
      } else {
        LeftRotate(*gp);
      }
    }
    PushUp(*gp);
    c = *gp;
    p = nd_[c].p;
  }
}
void Treap::InsertBefore(const int x, const int y, int& r) {
  if (y == 0) {
    Insert(x, false, r);
    return;
  }
  Insert(x, false, nd_[y].l);
  nd_[nd_[y].l].p = y;
  int p = y;
  int c = nd_[y].l;
  while (0 != p) {
    ++nd_[p].s;
    int* gp = &r;
    if (0 != nd_[p].p) {
      if (p == nd_[nd_[p].p].l) {
        gp = &nd_[nd_[p].p].l;
      } else {
        gp = &nd_[nd_[p].p].r;
      }
    }
    if (nd_[c].w > nd_[p].w) {
      if (c == nd_[p].l) {
        RightRotate(*gp);
      } else {
        LeftRotate(*gp);
      }
    }
    PushUp(*gp);
    c = *gp;
    p = nd_[c].p;
  }
}
void Treap::Delete(const int x, int& r) {
  int y = nd_[x].p;
  while (0 != y) {
    --nd_[y].s;
    y = nd_[y].p;
  }
  while (nd_[x].l != 0 || nd_[x].r != 0) {
    int* p = &r;
    if (0 != nd_[x].p) {
      if (x == nd_[nd_[x].p].l) {
        p = &nd_[nd_[x].p].l;
      } else {
        p = &nd_[nd_[x].p].r;
      }
    }
    if (nd_[nd_[x].l].w > nd_[nd_[x].r].w) {
      RightRotate(*p);
    } else {
      LeftRotate(*p);
    }
    --nd_[*p].s;
  }
  if (nd_[x].p == 0) {
    r = 0;
  } else if (nd_[nd_[x].p].l == x) {
    nd_[nd_[x].p].l = 0;
  } else {
    nd_[nd_[x].p].r = 0;
  }
  y = nd_[x].p;
  while (0 != y) {
    PushUp(y);
    y = nd_[y].p;
  }
}
int Treap::FindPrev(int x) const {
  if (nd_[x].l) {
    x = nd_[x].l;
    while (nd_[x].r) x = nd_[x].r;
    return x;
  }
  while (nd_[x].p && x == nd_[nd_[x].p].l) x = nd_[x].p;
  return nd_[x].p;
}
int Treap::Merge(const int r1, const int r2) {
  nd_[1].l = r1;
  nd_[r1].p = 1;
  nd_[1].r = r2;
  nd_[r2].p = 1;
  nd_[1].p = 0;
  nd_[1].s = nd_[r1].s + nd_[r2].s + 1;
  int r = 0;
  Delete(1, r);
  return r;
}
int Treap::Rank(const int x) const {
  int rank = nd_[nd_[x].l].s + 1;
  int y = x;
  int p = nd_[y].p;
  while (0 != p) {
    if (nd_[p].r == y) {
      rank += nd_[p].s - nd_[y].s;
    }
    y = p;
    p = nd_[p].p;
  }
  return rank;
}
int Treap::Select(const int r, const int rank) const {
  ASSERT(rank >= 1 && rank <= nd_[r].s);
  const int local_rank = nd_[nd_[r].l].s + 1;
  if (local_rank == rank) {
    return r;
  } else if (local_rank > rank) {
    return Select(nd_[r].l, rank);
  } else {
    return Select(nd_[r].r, rank - local_rank);
  }
}
int Treap::Root(const int x) const {
  int r = x;
  while (0 != nd_[r].p) {
    r = nd_[r].p;
  }
  return r;
}
int Treap::Minimum(const int x) const {
  int m = Root(x);
  while (0 != nd_[m].l) {
    m = nd_[m].l;
  }
  return m;
}
int Treap::Maximum(const int x) const {
  int m = Root(x);
  while (0 != nd_[m].r) {
    m = nd_[m].r;
  }
  return m;
}
int Treap::Size(const int r) const {
  return nd_[r].s;
}
void Treap::Check(const int r) const {
  int y = 146379;
  ASSERT(nd_[nd_[y].p].l == y || nd_[nd_[y].p].r == y);
}
void Treap::LeftRotate(int& x) {
  const int y = nd_[x].r;
  nd_[x].r = nd_[y].l;
  nd_[y].l = x;
  nd_[y].p = nd_[x].p;
  nd_[nd_[x].r].p = x;
  nd_[x].p = y;
  nd_[y].s = nd_[x].s;
  nd_[x].s = nd_[nd_[x].l].s + nd_[nd_[x].r].s + 1;
  PushUp(x);
  x = y;
  PushUp(x);
}
void Treap::RightRotate(int& x) {
  const int y = nd_[x].l;
  nd_[x].l = nd_[y].r;
  nd_[y].r = x;
  nd_[y].p = nd_[x].p;
  nd_[nd_[x].l].p = x;
  nd_[x].p = y;
  nd_[y].s = nd_[x].s;
  nd_[x].s = nd_[nd_[x].l].s + nd_[nd_[x].r].s + 1;
  PushUp(x);
  x = y;
  PushUp(x);
}
void Treap::SubCheck(const int x) const {
  if (0 != nd_[x].l) {
    ASSERT(nd_[nd_[x].l].p == x);
    SubCheck(nd_[x].l);
  }
  if (0 != nd_[x].r) {
    ASSERT(nd_[nd_[x].r].p == x);
    SubCheck(nd_[x].r);
  }
  ASSERT(nd_[x].s == nd_[nd_[x].l].s + nd_[nd_[x].r].s + 1);
  ASSERT(nd_[x].w >= nd_[nd_[x].l].w && nd_[x].w >= nd_[nd_[x].r].w);
}
void Treap::Traverse(const int r) const{
  if (0 == r) return;
  Traverse(nd_[r].l);
  eprintf("%d ", nd_[r].d.v.u);
  Traverse(nd_[r].r);
  return;
}
void Treap::TraverseToVec(const int r, std::vector<int> &tmp) const{
  if (0 == r) return;
  TraverseToVec(nd_[r].l, tmp);
  tmp.push_back(nd_[r].d.v.u);
  TraverseToVec(nd_[r].r, tmp);
  return;
}
int Treap::FindNextKey(const Node x, const int r) const {
  if (0 == r) return 0;
  int res = 0;
  if (x < nd_[nd_[nd_[r].r].d.mnKey].d.v) {
    res = FindNextKey(x, nd_[r].l);
    if (res == 0 && nd_[r].d.isKey && x < nd_[r].d.v) res = r;
    if (res == 0 && nd_[nd_[r].r].d.mnKey) res = nd_[nd_[r].r].d.mnKey;
    return res;
  }
  return FindNextKey(x, nd_[r].r);
}
int Treap::FindPrevKey(const int x, const int r) const {
  if (0 == r) return 0;
  if (x <= nd_[nd_[r].l].s) return FindPrevKey(x, nd_[r].l);
  int res = FindPrevKey(x - nd_[nd_[r].l].s - 1, nd_[r].r);
  if (res == 0) {
    if (nd_[r].d.isKey) res = r;
    else res = nd_[nd_[r].l].d.mxKey;
  }
  ASSERT(res == 0 || nd_[res].d.isKey);
  return res;
}
void Treap::Print() const {
  for (int i = 0; i < nd_.size(); ++i) eprintf("(%d %d %d %d v = %d k = %d) ", i, nd_[i].p, nd_[i].l, nd_[i].r, nd_[i].d.v.u, nd_[i].d.isKey); eputs("");
}
void Treap::Print(const int x) const {
  while (0 != x) eprintf("x = %d s = %d\n", x, nd_[x].s);
}
Node Treap::Get(const int x) const { return nd_[x].d.v; }
void Treap::QueryCV(const int x, int &r, std::vector<std::pair<int, int> > &V) const{
  if (x == 0) return;
  if (nd_[x].r) QueryCV(nd_[x].r, r, V);
  if (r) {
    V.push_back(std::make_pair(nd_[x].d.v.u, nd_[x].d.isKey));
    r -= nd_[x].d.isKey;
    if (r) QueryCV(nd_[x].l, r, V);
  }
}
}  // namespace gadget
/*
namespace {
constexpr int n = 1000000;
}  // namespace

int main() {
  gadget::Treap tree(n);
  std::vector<int> roots(100, n);
  std::vector<int> counts(100, 0);
  std::vector<int> numbers(n);

  for (int i = 0; i < n; ++i) {
    numbers[i] = rand() % 100;
    ++counts[numbers[i]];
    tree.Insert(i, rand() % 2, roots[numbers[i]]);
  }
  for (int i = 0; i < 100; ++i) {
    if (roots[i] != n) {
      tree.Check(roots[i]);
      ASSERT(tree.Size(roots[i]) == counts[i]);
    }
  }
  {
    for (int i = 1; i < 100; ++i) {
      roots[0] = tree.Merge(roots[0], roots[i]);
    }
    tree.Check(roots[0]);
    ASSERT(tree.Size(roots[0]) == n);
    printf("tree size: %d\n", tree.Size(roots[0]));
  }
  /star
  for (int i = 0; i < n; ++i) {
    const int b = numbers[i];
    const int nb = rand() % 100;
    if (nb != b) {
      --counts[b];
      ++counts[nb];
      numbers[i] = nb;
      tree.Delete(i, roots[b]);
      tree.Insert(i, rand() % 2, roots[nb]);
    }
  }
  for (int i = 0; i < 100; ++i) {
    if (roots[i] != n) {
      tree.Check(roots[i]);
      ASSERT(tree.Size(roots[i]) == counts[i]);
    }
  }
  for (int i = 0; i < n; ++i) {
    const int b = numbers[i];
    tree.Delete(i, roots[b]);
    --counts[b];
  }
  for (int i = 0; i < 100; ++i) {
    ASSERT(roots[i] == n);
  }
  for (int i = 0; i < n; ++i) {
    tree.Insert(i, false, roots[0]);
  }
  ASSERT(tree.Size(roots[0]) == n);
  for (int i = 0; i < n; ++i) {
    ASSERT(tree.Rank(i, roots[0]) == i + 1);
  }
  for (int i = 0; i < n; ++i) {
    tree.Delete(i, roots[0]);
    tree.InsertAfter(i, (n - 1 + i) % n, roots[0]);
  }
  for (int i = 0; i < n; ++i) {
    ASSERT(tree.Rank(i, roots[0]) == i + 1);
  }
  tree.Check(roots[0]);
  star/
}
*/
