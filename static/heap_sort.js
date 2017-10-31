//https://stackoverflow.com/questions/951021/what-is-the-javascript-version-of-sleep/39914235#39914235
function sleep(ms) {
  return new Promise(resolve => setTimeout(resolve, ms));
  }

var swap = async function (a, i, j){
  var tmp = a[i];
  var tmp_i = i;
  var tmp_j = j;

  a[i] = a[j]; 
  a[j] = tmp;

  //var x_pos =[x(tmp_i), x(tmp_j)];
  var swap_bars = svg.selectAll(".bar")
                     .filter(function(d){return (d == a[tmp_i]) || (d == a[tmp_j]);});

  swap_bars.transition()
           .duration(50)
           .style("fill", "red")

  await sleep(500)

  move_bars = svg.selectAll(".bar")
                 .data(a, function(d){return d;})

  move_bars.transition()
           .duration(125)
           .attr("x", function(d, i) {return x(i);});

  await sleep(500)
  swap_bars.transition()
           .duration(50)
           .style("fill", "steelblue");
  await sleep(500)
  return Promise.resolve(true); 
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

async function down_heap(a, root, end){
  var swap_index = root;
  var left_child_ = left_child(root);
  var right_child_ = right_child(root);

  d3.select("#sort_phase")
    .text("Heaping")

  if ((left_child_ <= end) && (compare(a[swap_index], a[left_child_]))){
    swap_index = left_child_;
    }
  if ((right_child_ <= end) && (compare(a[swap_index], a[right_child_]))){
    swap_index = right_child_;
    }
  if (root != swap_index) {
    await swap(a, swap_index, root)
    await down_heap(a, swap_index, end)
  }
  else {return Promise.resolve(true);}
}

async function heapify(a, N){
  d3.select("#sort_phase")
    .text("Building Heap")
  parent_ = parent_node(N);
  start = parent_;
  left_child_ = left_child(parent_);
  right_child_ = right_child(parent_);
  while (start >= 0){
    await down_heap(a, start, N-1);
    start -= 1;
  }
  return Promise.resolve(true);
}

async function heap_sort(a){
  await heapify(a, a.length); 
  N = a.length - 1;
  while (N >= 0) {
    console.log(N);
    d3.select("#sort_phase")
      .text("Swapping")
    await swap(a, 0, N);
    N -= 1;
    await down_heap(a, 0, N);
  }
  d3.select("#sort_phase")
    .text("Sorted")
  return a;
}

