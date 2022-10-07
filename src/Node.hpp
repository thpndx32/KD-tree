template<typename T>
class Node{
    private:
    T value;
    int lvl;
    Node<T>* sons[2] = {nullptr,nullptr};
    public:
    Node(T val, int level = 0) {
        value=val;
        lvl=level;
        }
    Node();
    bool find(T val, Node<T>*& pointer) {
        if (value->first == val->first) return true;
        int i = lvl%((val->first).size());
        if((val->first)[i]<(value->first)[i]){
            pointer = sons[0];
            if (sons[0]!=nullptr)return (sons[0])->find(val,pointer);
            }
        else {
            pointer = sons[1];
            if (sons[1]!=nullptr)return (sons[1])->find(val,pointer);
            }
        return false;
    }
    bool insert(T val, Node<T>*& pointer) {
        if(find(val,pointer)) return false;
        pointer = new Node(val, lvl+1);
        return true;
    }
};