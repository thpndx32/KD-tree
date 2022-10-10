// Copyright

#ifndef SRC_KDTREE_HPP_
#define SRC_KDTREE_HPP_

#include <cmath>
#include <iostream>
#include <set>
#include <stdexcept>
#include <utility>
#include <vector>
#include <map>
#include <iterator>
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
    Node<value_type>* nodes=nullptr;
    Node<value_type>** ptr1;
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
  dimension_ = rhs.dimension_;
  size_ = rhs.size_;
  if (rhs.nodes!=nullptr) nodes = new Node<value_type>(*(rhs.nodes));
}

template <size_t N, typename ElemType>
KDTree<N, ElemType> &KDTree<N, ElemType>::operator=(const KDTree &rhs) {
  dimension_ = rhs.dimension_;
  size_ = rhs.size_;
  nodes = new Node<value_type>(*(rhs.nodes));
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
  Node<value_type>**  ptr;
  Node<value_type>* aux;
  aux =nodes;
  ptr = &aux;
  value_type node (pt,0);
  //std::cout<<"no\n";
  return nodes->find(node,ptr);
}

template <size_t N, typename ElemType>
void KDTree<N, ElemType>::insert(const Point<N> &pt, const ElemType &value) {
  Node<value_type>** ptr;
  ptr=&nodes;
  //std::cout << ptr<<"\n";
  value_type node (pt, value);
  if (nodes == nullptr) {
    nodes = new Node<value_type> (node);
    size_++;
    return;
  }
  //std::cout << ptr<<"\n";
  if(nodes->insert(node,ptr))size_++;
}

template <size_t N, typename ElemType>
ElemType &KDTree<N, ElemType>::operator[](const Point<N> &pt) {
  ptr1=&nodes;
  value_type node (pt,0);
  if (nodes == nullptr) {
    nodes = new Node<value_type> (node);
    size_++;
    return ((*ptr1)->get_sec());
  }
  if(nodes->insert_n(node,ptr1))size_++;
  //std::cout<<((*ptr1)->get_sec())<<"\n";
  return ((*ptr1)->get_sec());
}

template <size_t N, typename ElemType>
ElemType& KDTree<N, ElemType>::at(const Point<N> &pt) {
  Node<value_type>** ptr;
  ptr =&nodes;
  value_type node (pt,0);
  if (nodes==nullptr||!(nodes->find(node,ptr))) throw std::out_of_range("could not find node");
  return ((*ptr)->get_sec());
}

template <size_t N, typename ElemType>
const ElemType &KDTree<N, ElemType>::at(const Point<N> &pt) const {
  Node<value_type>**  ptr;
  Node<value_type>* aux;
  aux =nodes;
  ptr = &aux;
  value_type node (pt,0);
  if (nodes==nullptr||!(nodes->find(node,ptr))) throw std::out_of_range("could not find node");
  return ((*ptr)->get_sec());
}

template <size_t N, typename ElemType>
ElemType KDTree<N, ElemType>::knn_value(const Point<N> &key, size_t k) const {
  Node<value_type>**  ptr;
  Node<value_type>* aux;
  aux =nodes;
  ptr = &aux;
  value_type node (key,0);
  return nodes->knn(node,k,ptr);
}

template <size_t N, typename ElemType>
std::vector<ElemType> KDTree<N, ElemType>::knn_query(const Point<N> &key,
                                                     size_t k) const {
  // TODO(me): Fill this in.
  std::vector<ElemType> values;
  return values;
}

template <typename T> class sorted_list {
private:
  int size_;
  std::vector<T> vec;
  std::vector<double> dis;
  size_t dim;
  T value;

public:
  typedef decltype(value.second)& var;
  sorted_list<T>(T val,int size, size_t dime) { 
      value = val;
      size_ = size; 
      dim = dime;
    }
  void insert(T val) {
    double dist = distance(value,val);
    dis.push_back(dist);
    vec.push_back(val);
    for(int i=dis.size()-2;i>=0;i--){
      if(dis[i]>dis[i+1]){
        T temp = vec[i];
        vec[i]=vec[i+1];
        vec[i+1]=temp;
        double temp1 = dis[i];
        dis[i]=dis[i+1];
        dis[i+1]=temp1;
      } else break;
    }
    if (vec.size() > size_){
      vec.pop_back();
      dis.pop_back();
    }
  }
  std::map<decltype(value.second),int> map_knn(){
    std::map<decltype(value.second),int> map_k;
    for(int i =0;i<vec.size();i++){
    typename std::map<decltype(value.second),int>::iterator it;
      it = map_k.find(vec[i]->second);
      if (it!=map_k.end()){
        it->second++;
      } else map_k.insert(std::pair<var,int>(vec[i]->second,i));
    }
    return map_k;
  }
};

template<typename T>
class Node{
    private:
    T value;
    int lvl;
    int index;
    Node<T>* sons[2];
    public:
    typedef decltype(value.second)& var;
    Node<T>(T val, int level = 0) {
        value=val;
        lvl=level;
        index = lvl%((val.first).size());
        sons[0]=sons[1]= nullptr;
        }
    Node<T>();
    Node<T>(const Node<T>& node){
      value = node.value;
      lvl = node.lvl;
      index = lvl%((value.first).size());
      if (node.sons[0]!=nullptr)sons[0]= new Node<T>(*(node.sons[0]));
      else sons[0]=nullptr;
      if (node.sons[1]!=nullptr)sons[1]= new Node<T>(*(node.sons[1]));
      else sons[1]=nullptr;
    }
    T get_val(){
      return value;
    };
    var get_sec(){
      return value.second;
    }
    void update_second(size_t temp){
      value.second=temp;
    }
    bool find(T val, Node<T>**& pointer) {
        if (value.first == val.first) return true;
        if((val.first)[index]<(value.first)[index]){
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

    bool find_p(T val, Node<T>**& pointer) {
        std::cout <<"there\n";
        if (value.first == val.first) {
          std::cout <<"true\n";
          return true;
        }
        std::cout <<"here\n";
        for (int i = 0; i < (value.first).size(); i++){
          std::cout << (value.first)[i]<<" ";
        }
        std::cout << "\n";
        std::cout << ((val.first)[index])<<" "<<(value.first)[index]<<"\n";
        if((val.first)[index]<(value.first)[index]){
            pointer = &sons[0];
            std::cout << 0<<"\n";
            std::cout <<"puntero"<<pointer<<"\n";
            if (sons[0]!=nullptr) {
              std::cout << sons[0]<<"\n";
              return ((sons[0])->find_p(val,pointer));
            }
            }
        else {
            pointer = &sons[1];
            std::cout << 1<<"\n";
            std::cout <<"puntero"<<pointer<<"\n";
            if (sons[1]!=nullptr) {
              std::cout << sons[1]<<"\n";
              return ((sons[1])->find_p(val,pointer));
              }
            }
        std::cout << "false\n";
        return false;
    }

    var max(std::map<decltype(value.second),int> map_knn){
      int temp = 0;
      decltype(value.second) maxim =map_knn.begin()->first;
      for (auto const& x : map_knn){
          if (x.second > temp){ 
            temp = x.second;
            maxim = x.first;
          }
      }
      return maxim;
    }

    var knn(T val, size_t size,Node<T>**& pointer){
      sorted_list<T> best_points(val,size,(val.first).size());
      std::vector <Node<T>*> vec_ptr;
      size_t k =0;
      while(pointer != nullptr){
        vec_ptr.push_back(*pointer);
        best_points.insert((*pointer)->get_val());
        if((val.first)[index]<(value.first)[index]){
            pointer = &sons[0];
            }
        else {
            pointer = &sons[1];
            }
      }
      return max(best_points.map_knn());
    }


    bool insert(T val, Node<T>**& pointer) {
        if(find(val,pointer)) {
          (*pointer)->update_second(val.second);
          return false;
        }
        *pointer = new Node(val, lvl+1);
        return true;
    }

    bool insert_n(T val, Node<T>**& pointer) {
        if(find(val,pointer)) {
          return false;
        }
        *pointer = new Node(val, lvl+1);
        return true;
    }
};

// TODO(me): finish the implementation of the rest of the KDTree class

#endif  // SRC_KDTREE_HPP_
