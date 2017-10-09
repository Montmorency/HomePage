var swap = function (a, i, j){
  var tmp = a[i];
  a[i] = a[j]; 
  a[j] = tmp;
};

function compare(x, y){
  if (x < y) {
    return true;}
  else {
    return false;
  }
}

function parent_node(d){
  return Math.floor((d-1)/2);}

function left_child(d){
  return 2*d + 1;}

function right_child(d){
  return 2*d + 2;}

function down_heap(a, root, end){
  var swap_index = root;
  var left_child_ = left_child(root);
  var right_child_ = right_child(root);
  if ((left_child_ <= end) && (compare(a[swap_index], a[left_child_]))){
    swap_index = left_child_;
    }
  if ((right_child_ <= end) && (compare(a[swap_index], a[right_child_]))){
    swap_index = right_child_;
    }
  if (root != swap_index) {
    swap(a, swap_index, root);
    down_heap(a, swap_index, end)
  }
}

function heapify(a, N){
  parent_ = parent_node(N);
  start = parent_;
  left_child_ = left_child(parent_);
  right_child_ = right_child(parent_);
  while (start >= 0){
    down_heap(a, start, N-1);
    start -= 1;
  }
}

function heap_sort(a){
  heapify(a, a.length); 
  N = a.length - 1;
  while (N >= 0) {
    console.log(N);
    swap(a, 0, N);
    N -= 1;
    down_heap(a, 0, N);
  }
  return a;
}


