// Copyright

#ifndef SRC_KDTREE_HPP_
#define SRC_KDTREE_HPP_

#include <cmath>
#include <iostream>
#include <set>
#include <stdexcept>
#include <utility>
#include <vector>
#include "Point.hpp"

template <typename T> class Node;

template <size_t N, typename ElemType>
class KDTree {
 public:
  typedef std::pair<Point<N>, ElemType> value_type;

  KDTree();

  ~KDTree();

  KDTree(const KDTree &rhs);
  KDTree &operator=(const KDTree &rhs);

  size_t dimension() const;

  size_t size() const;
  bool empty() const;

  bool contains(const Point<N> &pt) const;

  void insert(const Point<N> &pt, const ElemType &value);

  ElemType &operator[](const Point<N> &pt);

  ElemType &at(const Point<N> &pt);
  const ElemType &at(const Point<N> &pt) const;

  ElemType knn_value(const Point<N> &key, size_t k) const;

  std::vector<ElemType> knn_query(const Point<N> &key, size_t k) const;

 private:
    Node<value_type*>* nodes=nullptr;
    Node<value_type*>** ptr1;
    size_t dimension_;
    size_t size_;
};

template <size_t N, typename ElemType>
KDTree<N, ElemType>::KDTree() {
  dimension_ = N;
  size_ = 0;
}

template <size_t N, typename ElemType>
KDTree<N, ElemType>::~KDTree() {
  // TODO(me): Fill this in.
}

template <size_t N, typename ElemType>
KDTree<N, ElemType>::KDTree(const KDTree &rhs) {
  // TODO(me): Fill this in.
}

template <size_t N, typename ElemType>
KDTree<N, ElemType> &KDTree<N, ElemType>::operator=(const KDTree &rhs) {
  // TODO(me): Fill this in.
  return *this;
}

template <size_t N, typename ElemType>
size_t KDTree<N, ElemType>::dimension() const {
  return dimension_;
}

template <size_t N, typename ElemType>
size_t KDTree<N, ElemType>::size() const {
  return size_;
}

template <size_t N, typename ElemType>
bool KDTree<N, ElemType>::empty() const {
  if (size_) return false;
  return true;
}

template <size_t N, typename ElemType>
bool KDTree<N, ElemType>::contains(const Point<N> &pt) const {
  Node<value_type*>** ptr;
  value_type* node = new value_type(pt,0);
  return nodes->find(node,ptr);
}

template <size_t N, typename ElemType>
void KDTree<N, ElemType>::insert(const Point<N> &pt, const ElemType &value) {
  Node<value_type*>** ptr;
  ptr=&nodes;
  //std::cout << ptr<<"\n";
  value_type* node = new value_type(pt, value);
  if (nodes == nullptr) {
    nodes = new Node<value_type*> (node);
    size_++;
    return;
  }
  //std::cout << ptr<<"\n";
  if(nodes->insert(node,ptr))size_++;
}

template <size_t N, typename ElemType>
ElemType &KDTree<N, ElemType>::operator[](const Point<N> &pt) {
  ptr1=&nodes;
  value_type* node = new value_type(pt,0);
  if (nodes == nullptr) {
    nodes = new Node<value_type*> (node);
    size_++;
    return ((*ptr1)->get_sec());
  }
  if(nodes->insert(node,ptr1))size_++;
  std::cout<<((*ptr1)->get_sec())<<"\n";
  return ((*ptr1)->get_sec());
}

template <size_t N, typename ElemType>
ElemType& KDTree<N, ElemType>::at(const Point<N> &pt) {
  Node<value_type*>** ptr;
  ptr =&nodes;
  value_type* node = new value_type(pt,0);
  if (nodes==nullptr||!(nodes->find(node,ptr))) throw std::out_of_range("could not find node");
  return ((*ptr)->get_val())->second;
}

template <size_t N, typename ElemType>
const ElemType &KDTree<N, ElemType>::at(const Point<N> &pt) const {
  Node<value_type*>**  ptr;
  Node<value_type*>* aux;
  aux =nodes;
  ptr = &aux;
  value_type* node = new value_type(pt,0);
  if (nodes==nullptr||!(nodes->find(node,ptr))) throw std::out_of_range("could not find node");
  return ((*ptr)->get_val())->second;
}

template <size_t N, typename ElemType>
ElemType KDTree<N, ElemType>::knn_value(const Point<N> &key, size_t k) const {
  // TODO(me): Fill this in.
  ElemType new_element;
  return new_element;
}

template <size_t N, typename ElemType>
std::vector<ElemType> KDTree<N, ElemType>::knn_query(const Point<N> &key,
                                                     size_t k) const {
  // TODO(me): Fill this in.
  std::vector<ElemType> values;
  return values;
}

template<typename T>
class Node{
    private:
    T value;
    int lvl;
    int index;
    Node<T>* sons[2];
    public:
    Node(T val, int level = 0) {
        value=val;
        lvl=level;
        index = lvl%((val->first).size());
        sons[0]=sons[1]= nullptr;
        }
    Node();
    T& get_val(){
      return value;
    };
    size_t& get_sec(){
      return value->second;
    }
    void update_second(size_t temp){
      value->second=temp;
    }
    bool find(T val, Node<T>**& pointer) {
        if (value->first == val->first) return true;
        for (int i = 0; i < (value->first).size(); i++){
          //std::cout << (value->first)[i]<<" ";
        }
        //std::cout << "\n";
        //std::cout << ((val->first)[index])<<" "<<(value->first)[index]<<"\n";
        if((val->first)[index]<(value->first)[index]){
            pointer = &sons[0];
            //std::cout << 0<<"\n";
            //std::cout <<"puntero"<<pointer<<"\n";
            if (sons[0]!=nullptr)return (sons[0])->find(val,pointer);
            }
        else {
            pointer = &sons[1];
            //std::cout << 1<<"\n";
            //std::cout <<"puntero"<<pointer<<"\n";
            if (sons[1]!=nullptr)return (sons[1])->find(val,pointer);
            }
        return false;
    }
    bool insert(T val, Node<T>**& pointer) {
        if(find(val,pointer)) {
          (*pointer)->update_second(val->second);
          return false;
        }
        *pointer = new Node(val, lvl+1);
        return true;
    }
};

// TODO(me): finish the implementation of the rest of the KDTree class

#endif  // SRC_KDTREE_HPP_
