// TODO: unit tests for `compact` (?)
// TODO: unit tests for `rebalance` functions? (proper signaling upward, need non-root rotations)
// TODO: unit tests for mutations with duplicate/missing nodes
// TODO: unit tests for `peek_first`
// TODO: additional checks for non-AVL invariants?



use crate::primitives::int::IntNonnegFitsInUsize;



pub trait AvlTreeItem {
    type Key: Ord;

    fn key(&self) -> Self::Key;
}



/// An AVL tree.
/// 
/// Uses growing-array-based storage to minimize memory allocations and also
/// hopefully improve cache performance. With the current implementation, cache
/// peformance may not be materially improved if the tree becomes very large
/// with a large number of alternating deletions and insertions, as the access
/// within it will become highly random. The intent is to experiment with
/// reshuffling strategies in the future to address this.
/// 
/// Values are lazy-deleted and compacted if the fraction deleted exceeds half,
/// in order to maintain amortized-logarithmic deletion of individual items.
/// Insertion will fill pending "deleted" slots, to minimize reallocation.
/// 
/// The current deletion logic requires items to implement `Default` due to the
/// use of `std::mem::take`; the intent is to remove this requirement in the
/// future; doing so will require introducing custom heap array management
/// because `Vec` does not support marking an item as uninitialized.
/// 
/// In the future, the default strategy may be changed to *not* re-fill deleted
/// slots and insert strictly at the end, exploiting the tendency for less
/// recently added nodes to be near the top of the tree and thereby attempting
/// to keep the upper levels of the tree approximately in contiguous memory
/// after compacting, rather than becoming more or less randomly ordered over
/// time. This may be complemented by a periodic reshuffling to force the nodes
/// into level order fully or approximately (with frequency somehow maintaining
/// amortized constant time operations, e.g. at the same time as compacting).
/// 
/// Invariants:
/// - 01: length of `items` must fit within max for the given int type
/// - 02: `child_left`, `child_right`, and `bal_factor` same length as `items`
/// - 03: `root`, 'child_left`, and `child_right` point to valid items
/// - 04: `deleted` points to lazy-deleted items currently in `items`
/// - 05: `deleted` contains no duplicates
/// - 06: structure describes a valid AVL tree 
pub struct AvlTree<T,I=usize>
    where T: AvlTreeItem + Default, I: IntNonnegFitsInUsize {
    
    root: I,
    items: Vec<T>,
    child_left: Vec<I>,
    child_right: Vec<I>,
    bal_factor: Vec<i8>,
    deleted: Vec<I>,
}

impl<T,I> AvlTree<T,I>
    where T: AvlTreeItem + Default, I: IntNonnegFitsInUsize {

    // use `I::MAX` to indicate a child is absent; this can't be a valid index,
    // per restriction that length must fit in I
    const NIL: I = I::MAX;

    pub fn new() -> Self {
        Self {
            root: Self::NIL,
            items: vec![],
            child_left: vec![],
            child_right: vec![],
            bal_factor: vec![],
            deleted: vec![],
        }
    }

    pub fn compact(&mut self) {
        let mut is_deleted = vec![false; self.items.len()];
        for d in &self.deleted {
            is_deleted[d.to_usize_unchecked()] = true;
        }
        
        // shrink arrays by number of deleted items; move non-deleted items forward as needed
        let mut d: usize = 0;
        let mut s: usize = self.items.len() - 1;
        while d < s {
            if !is_deleted[d] {
                d += 1;
                continue;
            }
            if is_deleted[s] {
                s -= 1;
                continue;
            }
            self.items.swap(d, s);
            self.child_left.swap(d, s);
            self.child_right.swap(d, s);
            self.bal_factor.swap(d, s);
            if self.root.to_usize_unchecked() == s {
                self.root = I::from_usize_unchecked(d);
            }

            // temporarily steal `child_left` in the truncated range to map new indices
            self.child_left[s] = I::from_usize_unchecked(d);

            d += 1;
            s -= 1;
        }

        // update child pointers as needed
        let new_len = self.items.len() - self.deleted.len();
        for i in 0..new_len {
            let cl = self.child_left[i];
            let cl_usize = cl.to_usize_unchecked();
            if cl != Self::NIL && cl_usize >= new_len {
                self.child_left[i] = self.child_left[cl_usize];  // temporary mapping
            }
            let cr = self.child_right[i];
            let cr_usize = cr.to_usize_unchecked();
            if cr != Self::NIL && cr_usize >= new_len {
                self.child_right[i] = self.child_left[cr_usize];  // temporary mapping
            }
        }

        self.items.truncate(new_len);
        self.child_left.truncate(new_len);
        self.child_right.truncate(new_len);
        self.bal_factor.truncate(new_len);
        self.deleted.clear();
    }

    pub fn delete(&mut self, key: T::Key) -> Option<T> {
        // code below assumes the root exists; if not, nothing to delete
        if self.root == Self::NIL {
            return None;
        }

        // search for the item
        let d: I;
        let mut search_stack: Vec<(I,i8)> = vec![(self.root, 0)];
        loop {
            let i = search_stack.last().unwrap().0;
            let i_usize = i.to_usize_unchecked();
            let i_key = self.items[i_usize].key();
            if key < i_key {
                let cl = self.child_left[i_usize];
                if cl == Self::NIL {
                    return None;  // item not found
                }
                else {
                    search_stack.push((cl, -1));  // continue search to left
                }
            }
            else if key > i_key {
                let cr = self.child_right[i_usize];
                if cr == Self::NIL {
                    return None;  // item not found
                }
                else {
                    search_stack.push((cr, 1)); // continue search to right
                }
            }
            else {
                // found matching item
                d = i;
                break;
            }
        }

        // delete the item
        let d_usize = d.to_usize_unchecked();
        let cl = self.child_left[d_usize];
        let cr = self.child_right[d_usize];
        if cl == Self::NIL && cr == Self::NIL {
            self.delete_replace_simple(&mut search_stack, Self::NIL);
        }
        else if cr == Self::NIL {
            self.delete_replace_simple(&mut search_stack, cl);
        }
        else if cl == Self::NIL {
            self.delete_replace_simple(&mut search_stack, cr);
        }
        else {
            self.delete_replace_compound(&mut search_stack);
        }

        self.rebalance_after_delete(search_stack);

        // mark as deleted and extract item
        self.deleted.push(d);
        let item = std::mem::take(&mut self.items[d_usize]);

        // compact if it's time to do so
        if self.deleted.len() >= (self.items.len() + 1) / 2 {
            self.compact();
        }

        return Some(item);
    }

    /// Inserts an item into this AVL tree.
    /// If an existing item has the same key, replaces and returns it.
    pub fn insert(&mut self, item: T) -> Option<T> {
        // validate invariant 01: number of items must fit within max for the given int type
        if I::CHECK_CARDINALITY {
            assert!(self.items.len() < I::MAX.to_usize_unchecked() || self.deleted.len() > 0);
        }

        // code below assumes the root exists; if not, set this node as root
        if self.root == Self::NIL {
            self.root = self.insert_item(item);
            return None;
        }

        // insert item (no rebalancing yet)
        let mut search_stack: Vec<(I,i8)> = vec![(self.root, 0)];
        let key = item.key();
        loop {
            let i = search_stack.last().unwrap().0;
            let i_usize = i.to_usize_unchecked();
            let i_key = self.items[i_usize].key();
            if key < i_key {
                let cl = self.child_left[i_usize];
                if cl == Self::NIL {
                    // insert as left child
                    let new_idx = self.insert_item(item); 
                    self.child_left[i_usize] = new_idx;
                    self.bal_factor[i_usize] -= 1;
                    break;
                }
                else {
                    // continue searching to left of this node
                    search_stack.push((cl, -1));
                }
            }
            else if key > i_key {
                let cr = self.child_right[i_usize];
                if cr == Self::NIL {
                    // insert as right child
                    let new_idx = self.insert_item(item);
                    self.child_right[i_usize] = new_idx;
                    self.bal_factor[i_usize] += 1;
                    break;
                }
                else {
                    // continue searching to right of this node
                    search_stack.push((cr, 1));
                }
            }
            else {
                // found a matching key; replace it in-place and terminate
                let existing_value = std::mem::replace(&mut self.items[i_usize], item);
                return Some(existing_value);
            }
        }

        self.rebalance_after_insert(search_stack);

        return None;  // did not replace a value with same key
    }

    pub fn len(&self) -> I {
        // unchecked conversion valid per invariants 01, 04, 05
        I::from_usize_unchecked(self.items.len() - self.deleted.len())
    }

    pub fn peek_first(&self) -> Option<&T> {
        if self.root == Self::NIL {
            return None;
        }

        // get leftmost node
        let mut node = self.root.to_usize_unchecked();
        while self.child_left[node] != Self::NIL {
            node = self.child_left[node].to_usize_unchecked();
        }

        Some(&self.items[node])
    }

    pub fn pop_first(&mut self) -> Option<T> {
        // code below assumes the root exists; if not, nothing to pop
        if self.root == Self::NIL {
            return None;
        }

        // dive left
        let mut search_stack: Vec<(I,i8)> = vec![(self.root, 0)];
        loop {
            let i = search_stack.last().unwrap().0;
            let i_usize = i.to_usize_unchecked();
            let cl = self.child_left[i_usize];
            if cl == Self::NIL {
                break;
            }
            else {
                search_stack.push((cl, -1));
            }
        }
        let d = search_stack.last().unwrap().0;

        // delete the item (simple deletion: known to have no left child)
        let d_usize = d.to_usize_unchecked();
        self.delete_replace_simple(&mut search_stack, self.child_right[d_usize]);
        self.rebalance_after_delete(search_stack);

        // mark as deleted and extract item
        self.deleted.push(d);
        let item = std::mem::take(&mut self.items[d_usize]);

        // compact if it's time to do so
        if self.deleted.len() >= (self.items.len() + 1) / 2 {
            self.compact();
        }

        return Some(item);
    }

    // Private methods

    /// Compound replacement when deleting the node at the end of the search
    /// stack, for cases when it has both children.
    /// Modifies search stack so that end is where rebalancing should start.
    fn delete_replace_compound(&mut self, search_stack: &mut Vec<(I,i8)>) {
        let ss_len = search_stack.len();
        let d = search_stack[ss_len - 1].0;
        let d_usize = d.to_usize_unchecked();
        let cr = self.child_right[d_usize];

        // locate immediate successor of deleted node
        let mut search_stack_successor: Vec<(I,i8)> = vec![(cr, 1)];
        let mut s = cr;  // here known non-NIL
        loop {
            let s_usize = s.to_usize_unchecked();
            let cl2 = self.child_left[s_usize];
            if cl2 == Self::NIL {
                break;
            }
            else {
                search_stack_successor.push((cl2, -1));  // TEST: this is wrong! should be -1
                s = cl2;
            }
        }
        let s_usize = s.to_usize_unchecked();

        // replace d with s
        if ss_len > 1 {
            let p = search_stack[ss_len - 2].0;
            let p_usize = p.to_usize_unchecked();
            let dir = search_stack[ss_len - 1].1;
            match dir {
                -1 => self.child_left[p_usize] = s,
                1 => self.child_right[p_usize] = s,
                _ => unreachable!(),
            }
        }
        else {
            // d was root
            self.root = s;
        }

        // replace d with s in the stack
        search_stack[ss_len - 1].0 = s;

        // s inherits d's balance factor tentatively; might still need
        // update/rebalancing but only if rebalancing reaches this height
        self.bal_factor[s_usize] = self.bal_factor[d_usize];

        if search_stack_successor.len() > 1 {
            // s inherits d's children
            let sr_prev = self.child_right[s_usize];
            self.child_left[s_usize] = self.child_left[d_usize];
            self.child_right[s_usize] = self.child_right[d_usize];

            // update s's old parent to point to s's old right child
            let sp = search_stack_successor[search_stack_successor.len() - 2].0;
            let sp_usize = sp.to_usize_unchecked();
            self.child_left[sp_usize] = sr_prev;
            self.bal_factor[sp_usize] += 1;  // sp's left branch shrank by 1

            // rebalancing should start at s's old parent; extend search_stack
            search_stack.extend(search_stack_successor[0..(search_stack_successor.len()-1)].iter());
        }
        else {
            // s was d's right child; if s has a right child it should come along
            // (note s is known to have no left child before)
            self.child_left[s_usize] = self.child_left[d_usize];
            self.bal_factor[s_usize] -= 1;  // was d's old factor; right branch shrank by 1

            // don't extend search_stack, rebalancing should start at s
        }
    }

    /// Simple replacement when deleting the node at the end of the search
    /// stack, for cases when it has one child or zero (child = NIL).
    /// Modifies search stack so that end is where rebalancing should start.
    fn delete_replace_simple(&mut self, search_stack: &mut Vec<(I,i8)>, child: I) {
        if search_stack.len() > 1 {
            let p = search_stack[search_stack.len() - 2].0;
            let p_usize = p.to_usize_unchecked();
            let dir = search_stack[search_stack.len() - 1].1;
            match dir {
                -1 => {
                    self.child_left[p_usize] = child;
                    self.bal_factor[p_usize] += 1;
                },
                1 => {
                    self.child_right[p_usize] = child;
                    self.bal_factor[p_usize] -= 1;
                },
                _ => unreachable!(),
            }
        }
        else {            
            self.root = child;
            // becoming the root does not affect balance factor of `child`
        }

        // balance of child branch is unaffected; begin rebalancing at parent
        search_stack.pop();
    }

    fn insert_item(&mut self, item: T) -> I {
        if self.deleted.len() > 0 {
            self.insert_item_overwrite_deleted(item)
        }
        else {
            self.insert_item_extend_array(item)
        }  
    }


    fn insert_item_extend_array(&mut self, item: T) -> I {
        // maintains invariants 02, 03, 04, 05;
        // caller responsible for other invariants
        let new_idx = I::from_usize_unchecked(self.items.len());
        self.items.push(item);
        self.child_left.push(Self::NIL);
        self.child_right.push(Self::NIL);
        self.bal_factor.push(0);

        new_idx
    }

    fn insert_item_overwrite_deleted(&mut self, item: T) -> I {
        // maintains invariants 02, 03, 04, 05;
        // caller responsible for other invariants
        let new_idx = self.deleted.pop().unwrap();
        let n_usize = new_idx.to_usize_unchecked();
        self.items[n_usize] = item;
        self.child_left[n_usize] = Self::NIL;
        self.child_right[n_usize] = Self::NIL;
        self.bal_factor[n_usize] = 0;

        new_idx
    }

    fn rebalance_after_delete(&mut self, mut search_stack: Vec<(I,i8)>) {
        let mut prev_branch_dir: i8 = 0;
        let mut prev_branch_height_decr = true;
        let mut prev_branch_new_head: Option<I> = None;

        while search_stack.len() > 0 {
            let (i, dir) = search_stack.pop().unwrap();
            let i_usize = i.to_usize_unchecked();

            // update child pointer if child branch was rotated
            if let Some(c) = prev_branch_new_head {
                match prev_branch_dir {
                    -1 => self.child_left[i_usize] = c,
                    1 => self.child_right[i_usize] = c,
                    _ => unreachable!(),
                }
            }

            // update balance factor if height of child branch changed
            if prev_branch_height_decr {
                match prev_branch_dir {
                    -1 => self.bal_factor[i_usize] += 1,
                    1 => self.bal_factor[i_usize] -= 1,
                    0 => (),  // parent of deleted node, already had b.f. updated
                    _ => unreachable!(),
                }
            }
            else {
                // child branch did not change height; done rebalancing
                return;
            }

            // rebalance this node's branch
            let (branch_new_head, height_decr) = self.rebalance_after_delete_node(i);

            prev_branch_dir = dir;
            prev_branch_new_head = branch_new_head;
            prev_branch_height_decr = height_decr;
        }

        // reassign root pointer if rotation happened at root
        if let Some(n) = prev_branch_new_head {
            self.root = n;
        }
    }

    fn rebalance_after_delete_node(&mut self, node: I) -> (Option<I>, bool) {
        let n_usize = node.to_usize_unchecked();
        let n_bal = self.bal_factor[n_usize];

        if n_bal == -2 {
            let cl = self.child_left[n_usize];
            let cl_usize = cl.to_usize_unchecked();
            let cl_bal_before = self.bal_factor[cl_usize];
            if cl_bal_before == -1 || cl_bal_before == 0 {
                self.rotate_right(node);
                return (Some(cl), cl_bal_before == -1);
            }
            else if cl_bal_before == 1 {
                let new_head = self.child_right[cl_usize];
                self.rotate_leftright(node);
                return (Some(new_head), true);
            }
            else {
                unreachable!();
            }
        }
        else if n_bal == 2 {
            let cr = self.child_right[n_usize];
            let cr_usize = cr.to_usize_unchecked();
            let cr_bal_before = self.bal_factor[cr_usize];
            if cr_bal_before == 1 || cr_bal_before == 0 {
                self.rotate_left(node);
                return (Some(cr), cr_bal_before == 1);
            }
            else if self.bal_factor[cr_usize] == -1 {
                let new_head = self.child_left[cr_usize];
                self.rotate_rightleft(node);
                return (Some(new_head), true);
            }
            else {
                unreachable!();
            }
        }
        else if n_bal == -1 || n_bal == 1 {
            // no rotation needed and height of this branch has not changed
            // (balance factor must have been zero before, because it must have
            // changed if rebalancing has reached this node)
            return (None, false);
        }
        else if n_bal == 0 {
            // no rotation needed but height has decreased (balance factor must
            // have been +/- 1 before, and taller branch must have shrunk)
            return (None, true);
        }
        else {
            unreachable!();
        }
    }

    fn rebalance_after_insert(&mut self, mut search_stack: Vec<(I,i8)>) {
        let mut prev_branch_dir: i8 = 0;
        let mut prev_branch_height_incr= true;
        let mut prev_branch_new_head: Option<I> = None;

        while search_stack.len() > 0 {
            let (i, dir) = search_stack.pop().unwrap();
            let i_usize = i.to_usize_unchecked();

            // update child pointer if child branch was rotated
            if let Some(c) = prev_branch_new_head {
                match prev_branch_dir {
                    -1 => self.child_left[i_usize] = c,
                    1 => self.child_right[i_usize] = c,
                    _ => unreachable!(),
                }
            }

            // update balance factor if height of child branch changed
            if prev_branch_height_incr {
                match prev_branch_dir {
                    -1 => self.bal_factor[i_usize] -= 1,
                    1 => self.bal_factor[i_usize] += 1,
                    0 => (),  // parent of inserted node, already had b.f. updated
                    _ => unreachable!(),
                }
            }
            else {
                // child branch did not change height; done rebalancing
                return;
            }

            // rebalance this node's branch
            let (branch_new_head, height_incr) = self.rebalance_after_insert_node(i);

            prev_branch_dir = dir;
            prev_branch_new_head = branch_new_head;
            prev_branch_height_incr = height_incr;
        }

        // reassign root pointer if rotation happened at root
        if let Some(n) = prev_branch_new_head {
            self.root = n;
        }
    }

    /// Returns the index of the new head of this branch if a rotation was
    /// performed, or `None` if not. Also returns whether the height increased.
    /// Precondition: `self.bal_factor[node]` must be updated already.
    fn rebalance_after_insert_node(&mut self, node: I) -> (Option<I>, bool) {
        let n_usize = node.to_usize_unchecked();
        let n_bal = self.bal_factor[n_usize];

        // after rotation, the height of the branch becomes the same as before
        // (return `false` for second output)
        if n_bal == -2 {
            let cl = self.child_left[n_usize];
            let cl_usize = cl.to_usize_unchecked();
            if self.bal_factor[cl_usize] == -1 {
                self.rotate_right(node);
                return (Some(cl), false);
            }
            else if self.bal_factor[cl_usize] == 1 {
                let new_head = self.child_right[cl_usize];
                self.rotate_leftright(node);
                return (Some(new_head), false);
            }
            else {
                unreachable!();
            }
        }
        else if n_bal == 2 {
            let cr = self.child_right[n_usize];
            let cr_usize = cr.to_usize_unchecked();
            if self.bal_factor[cr_usize] == 1 {
                self.rotate_left(node);
                return (Some(cr), false);
            }
            else if self.bal_factor[cr_usize] == -1 {
                let new_head = self.child_left[cr_usize];
                self.rotate_rightleft(node);
                return (Some(new_head), false);
            }
            else {
                unreachable!();
            }
        }
        else if n_bal == -1 || n_bal == 1 {
            // no rotation needed here but the height of this branch increased
            // (balance factor must have been zero before, because it must have
            // changed if rebalancing has reached this node); need to signal
            // the height increase to this node's parent, update the parent's
            // balance factor, and assess whether the parent needs rotation
            return (None, true);
        }
        else if n_bal == 0 {
            // no height increase (balance factor must have been +/- 1 before,
            // height of shorter sub-branch must have been increased to match)
            return (None, false);
        }
        else {
            unreachable!();
        }
    }

    fn rotate_left(&mut self, node: I) {
        let n_usize = node.to_usize_unchecked();
        let cr = self.child_right[n_usize];
        let cr_usize = cr.to_usize_unchecked();

        self.child_right[n_usize] = self.child_left[cr_usize];
        self.child_left[cr_usize] = node;
        if self.bal_factor[cr_usize] == 1 {
            self.bal_factor[n_usize] = 0;
            self.bal_factor[cr_usize] = 0;
        }
        else if self.bal_factor[cr_usize] == 0 {
            self.bal_factor[n_usize] = 1;
            self.bal_factor[cr_usize] = -1;
        }
        else {
            unreachable!();
        }
    }

    fn rotate_leftright(&mut self, node: I) {
        let n_usize = node.to_usize_unchecked();
        let cl = self.child_left[n_usize];
        let cl_usize = cl.to_usize_unchecked();
        let clr = self.child_right[cl_usize];
        let clr_usize = clr.to_usize_unchecked();

        self.child_right[cl_usize] = self.child_left[clr_usize];
        self.child_left[clr_usize] = cl;
        self.child_left[n_usize] = self.child_right[clr_usize];
        self.child_right[clr_usize] = node;
        if self.bal_factor[clr_usize] == 0 {
            // all three nodes balanced after rotation (clr still = 0)
            self.bal_factor[n_usize] = 0;
            self.bal_factor[cl_usize] = 0;
        }
        else if self.bal_factor[clr_usize] == 1 {
            self.bal_factor[n_usize] = 0;
            self.bal_factor[cl_usize] = -1;
            self.bal_factor[clr_usize] = 0;
        }
        else if self.bal_factor[clr_usize] == -1 {
            self.bal_factor[n_usize] = 1;
            self.bal_factor[cl_usize] = 0;
            self.bal_factor[clr_usize] = 0;
        }
        else {
            unreachable!();
        }
    }

    fn rotate_right(&mut self, node: I) {
        let n_usize = node.to_usize_unchecked();
        let cl = self.child_left[n_usize];
        let cl_usize = cl.to_usize_unchecked();

        self.child_left[n_usize] = self.child_right[cl_usize];
        self.child_right[cl_usize] = node;
        if self.bal_factor[cl_usize] == -1 {
            self.bal_factor[n_usize] = 0;
            self.bal_factor[cl_usize] = 0;
        }
        else if self.bal_factor[cl_usize] == 0 {
            self.bal_factor[n_usize] = -1;
            self.bal_factor[cl_usize] = 1;
        }
        else {
            unreachable!();
        }
    }

    fn rotate_rightleft(&mut self, node: I) {
        let n_usize = node.to_usize_unchecked();
        let cr = self.child_right[n_usize];
        let cr_usize = cr.to_usize_unchecked();
        let crl = self.child_left[cr_usize];
        let crl_usize = crl.to_usize_unchecked();

        self.child_left[cr_usize] = self.child_right[crl_usize];
        self.child_right[crl_usize] = cr;
        self.child_right[n_usize] = self.child_left[crl_usize];
        self.child_left[crl_usize] = node;
        if self.bal_factor[crl_usize] == 0 {
            // all three nodes balanced after rotation (crl still = 0)
            self.bal_factor[n_usize] = 0;
            self.bal_factor[cr_usize] = 0;
        }
        else if self.bal_factor[crl_usize] == 1 {
            self.bal_factor[n_usize] = -1;
            self.bal_factor[cr_usize] = 0;
            self.bal_factor[crl_usize] = 0;
        }
        else if self.bal_factor[crl_usize] == -1 {
            self.bal_factor[n_usize] = 0;
            self.bal_factor[cr_usize] = 1;
            self.bal_factor[crl_usize] = 0;
        }
        else {
            unreachable!();
        }
    }

}



#[cfg(test)]
mod tests {

    use super::*;

    impl AvlTreeItem for i32 {
        type Key = i32;

        fn key(&self) -> i32 {
            *self
        }
    }

    fn util_delete_sequence<T>(tree: &mut AvlTree<T>, keys: Vec<T::Key>)
        where T: AvlTreeItem + Default {
        
        for x in keys {
            tree.delete(x);
        }
    }

    fn util_insert_sequence<T>(tree: &mut AvlTree<T>, items: Vec<T>)
        where T: AvlTreeItem + Default {

        for x in items {
            tree.insert(x);
        }
    }

    fn util_node_info<T>(tree: &AvlTree<T>, key: T::Key) -> (Option<T>, Option<T>, i8)
        where T: AvlTreeItem + Default + Copy {

        // find the node
        let nil = AvlTree::<T>::NIL;
        let mut n = tree.root;
        loop {
            let n_usize = n.to_usize_unchecked();
            let n_key = tree.items[n_usize].key();
            if key < n_key {
                n = tree.child_left[n_usize];
                assert_ne!(n, nil);
            }
            else if key > n_key {
                n = tree.child_right[n_usize];
                assert_ne!(n, nil);
            }
            else {
                break;
            }
        }

        // extract node info
        let n_usize = n.to_usize_unchecked();
        let cl = tree.child_left[n_usize];
        let cr = tree.child_right[n_usize];
        let cl = if cl == nil {None} else {Some(tree.items[cl.to_usize_unchecked()])};
        let cr = if cr == nil {None} else {Some(tree.items[cr.to_usize_unchecked()])};
        let bf = tree.bal_factor[n_usize];

        (cl, cr, bf)
    }

    fn util_pop_all<T>(tree: &mut AvlTree<T>) -> Vec<T>
        where T: AvlTreeItem + Default {

        let mut result = vec![];
        loop {
            let next = tree.pop_first();
            if let Some(x) = next {
                result.push(x);
            }
            else {
                break;
            }
        }
        result
    }

    fn util_verify_tree<T>(tree: &AvlTree<T>, root: T, nodes: Vec<(T::Key, Option<T>, Option<T>, i8)>)
        where T: AvlTreeItem + Default + Copy + PartialEq + std::fmt::Debug + std::fmt::Display, T::Key: Copy + std::fmt::Display {

        assert_eq!(tree.len(), nodes.len());
        assert_eq!(tree.items[tree.root], root, "Expected root to be {}, was {}", root, tree.items[tree.root]);

        for i in 0..nodes.len() {
            let info = util_node_info(&tree, nodes[i].0);
            
            // check left child
            let cl = nodes[i].1;
            if cl.is_none() {
                assert!(info.0.is_none(), "Expected {} cl to be None, was {}", nodes[i].0, info.0.unwrap());
            }
            else {
                assert_eq!(info.0.unwrap(), cl.unwrap(), "Expected {} cl to be {}, was {}", nodes[i].0, cl.unwrap(), info.0.unwrap());
            }

            // check right child
            let cr = nodes[i].2;
            if cr.is_none() {
                assert!(info.1.is_none(), "Expected {} cr to be None, was {}", nodes[i].0, info.1.unwrap());
            }
            else {
                assert_eq!(info.1.unwrap(), cr.unwrap(), "Expected {} cr to be {}, was {}", nodes[i].0, cr.unwrap(), info.1.unwrap());
            }

            // check balance factor
            assert_eq!(info.2, nodes[i].3, "Expected {} bf to be {}, was {}", nodes[i].0, nodes[i].3, info.2);
        }
    }

    #[test]
    fn test_000() {
        let add_seq = vec![4, 7, 1, 3, 9, 8, 2, 0, 5, 6];
        let del_seq = vec![4, 7, 2];

        let mut tree: AvlTree<i32> = AvlTree::new();
        util_insert_sequence(&mut tree, add_seq);
        util_delete_sequence(&mut tree, del_seq);
        let result = util_pop_all(&mut tree);

        assert_eq!(result, vec![0, 1, 3, 5, 6, 8, 9]);
    }

    #[test]
    fn test_001() {
        let add_seq1 = vec![12, 18, 2, 4, 8, 14, 10, 16, 0, 6];
        let del_seq1 = vec![10, 2, 8, 6, 12];
        let add_seq2 = vec![13, 9, 1, 7, 17, 19, 11, 5, 3, 15];
        let del_seq2 = vec![1, 11, 19, 15, 9, 18, 0, 16];

        let mut tree: AvlTree<i32> = AvlTree::new();
        util_insert_sequence(&mut tree, add_seq1);
        util_delete_sequence(&mut tree, del_seq1);
        util_insert_sequence(&mut tree, add_seq2);
        util_delete_sequence(&mut tree, del_seq2);
        let result = util_pop_all(&mut tree);

        assert_eq!(result, vec![3, 4, 5, 7, 13, 14, 17]);
    }

    #[test]
    fn test_002() {
        // alternating insertion/deletion without compacting
        let add_seq1 = vec![100, 300, 600, 500, 200, 400, 700];
        let del_seq1 = vec![300, 200];
        let add_seq2 = vec![50, 30, 70, 10];
        let del_seq2 = vec![500, 70];
        let add_seq3 = vec![450, 750, 250, 150, 350, 650, 550];
        let del_seq3 = vec![400, 30, 10, 250];

        let mut tree: AvlTree<i32> = AvlTree::new();
        util_insert_sequence(&mut tree, add_seq1);
        util_delete_sequence(&mut tree, del_seq1);
        util_insert_sequence(&mut tree, add_seq2);
        util_delete_sequence(&mut tree, del_seq2);
        util_insert_sequence(&mut tree, add_seq3);
        util_delete_sequence(&mut tree, del_seq3);
        let result = util_pop_all(&mut tree);

        assert_eq!(result, vec![50, 100, 150, 350, 450, 550, 600, 650, 700, 750]);
    }

    #[test]
    fn test_003() {
        // alternating insertion/deletion without compacting,
        // repeated deletion/insertion of same items

        let add_seq1 = vec![100, 300, 600, 500, 200, 400, 700];
        let del_seq1 = vec![300, 200];
        let add_seq2 = vec![50, 30, 70, 10];
        let del_seq2 = vec![500, 70];
        let add_seq3 = vec![450, 750, 250, 150, 350, 650, 550];
        let del_seq3 = vec![400, 30, 10, 250];

        let mut tree: AvlTree<i32> = AvlTree::new();
        util_insert_sequence(&mut tree, add_seq1.clone());
        util_insert_sequence(&mut tree, add_seq1.clone());
        util_delete_sequence(&mut tree, del_seq1.clone());
        util_delete_sequence(&mut tree, del_seq1.clone());
        util_insert_sequence(&mut tree, add_seq2.clone());
        util_insert_sequence(&mut tree, add_seq2.clone());
        util_delete_sequence(&mut tree, del_seq2.clone());
        util_delete_sequence(&mut tree, del_seq2.clone());
        util_insert_sequence(&mut tree, add_seq3.clone());
        util_insert_sequence(&mut tree, add_seq3.clone());
        util_delete_sequence(&mut tree, del_seq3.clone());
        util_insert_sequence(&mut tree, add_seq1.clone());
        util_delete_sequence(&mut tree, del_seq3.clone());
        let result = util_pop_all(&mut tree);

        assert_eq!(result, vec![50, 100, 150, 200, 300, 350, 450, 500, 550, 600, 650, 700, 750]);
    }

    #[test]
    fn test_004() {
        // testing with deleting a few items out of a big tree
        let start = 0;
        let end = 10000;

        let del_seq = vec![13, 7828, 2322, 853, 4390, 2958, 6530, 3893, 1098, 5928, 9237, 9999, 0, 1, 2, 3, 4, 5];

        let mut tree: AvlTree<i32> = AvlTree::new();
        for x in start..=end {
            tree.insert(x);
        }
        for x in &del_seq {
            tree.delete(*x);
        }

        let mut next = start;
        while next <= end {
            if !del_seq.contains(&next) {
                assert_eq!(next, tree.pop_first().unwrap());
            }
            next += 1;
        }

        assert_eq!(tree.len(), 0);
    }

    #[test]
    fn test_005() {
        // interleaved deletion sequences in a large tree
        
        let start = 1;
        let end = 100000;

        let mut tree: AvlTree<i32> = AvlTree::new();
        for x in start..(end/2) {
            let duplicate = tree.insert(x);
            assert!(duplicate.is_none());
        }
        for x in ((end/2)..=end).rev() {
            let duplicate = tree.insert(x);
            assert!(duplicate.is_none());
        }
        for x in start..=end {
            let duplicate = tree.insert(x);
            assert!(duplicate.is_some());
        }
        for x in start..=end {
            if x % 3 == 0 {
                tree.delete(x).unwrap();
            }
        }
        for x in start..=end {
            if x % 3 == 1 {
                tree.delete(x).unwrap();
            }
        }
        for x in start..=end {
            if x % 3 == 2 {
                tree.delete(x).unwrap();
            }
        }

        assert_eq!(tree.len(), 0);
    }

    #[test]
    fn test_006() {
        // delete multiples of 7 and 13 out of large tree, ensure rest intact
        
        let size = 100000;
        
        let mut tree: AvlTree<i32> = AvlTree::new();
        for x in 1..=size {
            tree.insert(x);
        }
        for x in (7..=size).step_by(7) {
            tree.delete(x);
        }
        for x in (13..=size).step_by(13) {
            tree.delete(x);
        }

        for x in 1..=size {
            if x % 7 != 0 && x % 13 != 0 {
                assert_eq!(x, tree.pop_first().unwrap());
            }
        }
    }

    #[test]
    fn test_007() {
        // interleaved insertion followed by interleaved deletion
        
        let size = 100000;

        let mut tree: AvlTree<i32> = AvlTree::new();
        for x in (1..=size).step_by(3) {
            tree.insert(x);
        }
        for x in (2..=size).step_by(3) {
            tree.insert(x);
        }
        for x in (3..=size).step_by(3) {
            tree.insert(x);
        }
        for x in (1..=size).step_by(4) {
            tree.delete(x);
        }
        for x in (2..=size).step_by(4) {
            tree.delete(x);
        }
        for x in (4..=size).step_by(4) {
            tree.delete(x);
        }

        let mut next = 3;
        while next <= size {
            assert_eq!(next, tree.pop_first().unwrap());
            next += 4;
        }
        assert_eq!(tree.len(), 0);
    }

    #[test]
    fn test_delete_replace_compound_000() {
        // case 1/6: at root, successor is right child

        let mut tree: AvlTree<i32> = AvlTree::new();
        util_insert_sequence(&mut tree, vec![8, 4, 12, 14]);
        util_verify_tree(&tree, 8, vec![
            (8, Some(4), Some(12), 1),
            (4, None, None, 0),
            (12, None, Some(14), 1),
            (14, None, None, 0),
        ]);

        tree.delete(8);

        util_verify_tree(&tree, 12, vec![
            (12, Some(4), Some(14), 0),
            (4, None, None, 0),
            (14, None, None, 0),
        ]);
    }

    #[test]
    fn test_delete_replace_compound_001() {
        // case 2/6: at root, successor is not right child

        let mut tree: AvlTree<i32> = AvlTree::new();
        util_insert_sequence(&mut tree, vec![8, 4, 12, 2, 10, 14, 11]);
        util_verify_tree(&tree, 8, vec![
            (8, Some(4), Some(12), 1),
            (4, Some(2), None, -1),
            (12, Some(10), Some(14), -1),
            (2, None, None, 0),
            (10, None, Some(11), 1),
            (14, None, None, 0),
            (11, None, None, 0),
        ]);

        tree.delete(8);

        util_verify_tree(&tree, 10, vec![
            (10, Some(4), Some(12), 0),
            (4, Some(2), None, -1),
            (12, Some(11), Some(14), 0),
            (2, None, None, 0),
            (11, None, None, 0),
            (14, None, None, 0),
        ]);
    }

    #[test]
    fn test_delete_replace_compound_002() {
        // case 3/6: at root, successor is not child or grandchild
        // (some retracing information in the successor's search stack
        // is not used unless more than 2 levels down)

        let mut tree: AvlTree<i32> = AvlTree::new();
        util_insert_sequence(&mut tree, vec![16, 8, 24, 4, 12, 20, 28, 6, 18, 22, 30, 19]);
        util_verify_tree(&tree, 16, vec![
            (16, Some(8), Some(24), 1),
            (8, Some(4), Some(12), -1),
            (24, Some(20), Some(28), -1),
            (4, None, Some(6), 1),
            (12, None, None, 0),
            (20, Some(18), Some(22), -1),
            (28, None, Some(30), 1),
            (6, None, None, 0),
            (18, None, Some(19), 1),
            (22, None, None, 0),
            (30, None, None, 0),
            (19, None, None, 0),
        ]);

        tree.delete(16);

        util_verify_tree(&tree, 18, vec![
            (18, Some(8), Some(24), 0),
            (8, Some(4), Some(12), -1),
            (24, Some(20), Some(28), 0),
            (4, None, Some(6), 1),
            (12, None, None, 0),
            (20, Some(19), Some(22), 0),
            (28, None, Some(30), 1),
            (6, None, None, 0),
            (19, None, None, 0),
            (22, None, None, 0),
            (30, None, None, 0),
        ]);
    }

    #[test]
    fn test_delete_replace_compound_003() {
        // case 4/6: not at root, successor is right child

        let mut tree: AvlTree<i32> = AvlTree::new();
        util_insert_sequence(&mut tree, vec![16, 8, 24, 4, 12, 28, 14]);
        util_verify_tree(&tree, 16, vec![
            (16, Some(8), Some(24), -1),
            (8, Some(4), Some(12), 1),
            (24, None, Some(28), 1),
            (4, None, None, 0),
            (12, None, Some(14), 1),
            (28, None, None, 0),
            (14, None, None, 0),
        ]);

        tree.delete(8);

        util_verify_tree(&tree, 16, vec![
            (16, Some(12), Some(24), 0),
            (12, Some(4), Some(14), 0),
            (24, None, Some(28), 1),
            (4, None, None, 0),
            (14, None, None, 0),
            (28, None, None, 0),
        ]);
    }

    #[test]
    fn test_delete_replace_compound_004() {
        // case 5/6: not at root, successor is not right child

        let mut tree: AvlTree<i32> = AvlTree::new();
        util_insert_sequence(&mut tree, vec![16, 8, 24, 4, 12, 20, 28, 6, 18, 26, 30, 27]);
        util_verify_tree(&tree, 16, vec![
            (16, Some(8), Some(24), 1),
            (8, Some(4), Some(12), -1),
            (24, Some(20), Some(28), 1),
            (4, None, Some(6), 1),
            (12, None, None, 0),
            (20, Some(18), None, -1),
            (28, Some(26), Some(30), -1),
            (6, None, None, 0),
            (18, None, None, 0),
            (26, None, Some(27), 1),
            (30, None, None, 0),
            (27, None, None, 0),
        ]);

        tree.delete(24);

        util_verify_tree(&tree, 16, vec![
            (16, Some(8), Some(26), 0),
            (8, Some(4), Some(12), -1),
            (26, Some(20), Some(28), 0),
            (4, None, Some(6), 1),
            (12, None, None, 0),
            (20, Some(18), None, -1),
            (28, Some(27), Some(30), 0),
            (6, None, None, 0),
            (18, None, None, 0),
            (27, None, None, 0),
            (30, None, None, 0),
        ]);
    }

    #[test]
    fn test_delete_replace_compound_005() {
        // case 6/6: not at root, successor is not child or grandchild
        // (some retracing information in the successor's search stack
        // is not used unless more than 2 levels down)

        let mut tree: AvlTree<i32> = AvlTree::new();
        util_insert_sequence(&mut tree, vec![16, 8, 24, 4, 12, 20, 28, 6, 18, 26, 30, 25]);
        util_verify_tree(&tree, 16, vec![
            (16, Some(8), Some(24), 1),
            (8, Some(4), Some(12), -1),
            (24, Some(20), Some(28), 1),
            (4, None, Some(6), 1),
            (12, None, None, 0),
            (20, Some(18), None, -1),
            (28, Some(26), Some(30), -1),
            (6, None, None, 0),
            (18, None, None, 0),
            (26, Some(25), None, -1),
            (30, None, None, 0),
            (25, None, None, 0),
        ]);

        tree.delete(24);

        util_verify_tree(&tree, 16, vec![
            (16, Some(8), Some(25), 0),
            (8, Some(4), Some(12), -1),
            (25, Some(20), Some(28), 0),
            (4, None, Some(6), 1),
            (12, None, None, 0),
            (20, Some(18), None, -1),
            (28, Some(26), Some(30), 0),
            (6, None, None, 0),
            (18, None, None, 0),
            (26, None, None, 0),
            (30, None, None, 0),
        ]);
    }

    #[test]
    fn test_rotate_left_delete_000() {
        // case 1/2 for deletion: right balance factor = 1

        let mut tree: AvlTree<i32> = AvlTree::new();
        util_insert_sequence(&mut tree, vec![8, 4, 12, 2, 10, 14, 13]);
        util_verify_tree(&tree, 8, vec![
            (8, Some(4), Some(12), 1),
            (4, Some(2), None, -1),
            (12, Some(10), Some(14), 1),
            (2, None, None, 0),
            (10, None, None, 0),
            (14, Some(13), None, -1),
            (13, None, None, 0),
        ]);
        
        tree.delete(4);

        util_verify_tree(&tree, 12, vec![
            (12, Some(8), Some(14), 0),
            (8, Some(2), Some(10), 0),
            (14, Some(13), None, -1),
            (2, None, None, 0),
            (10, None, None, 0),
            (13, None, None, 0),
        ]);
    }

    #[test]
    fn test_rotate_left_delete_001() {
        // case 2/2 for deletion: right balance factor = 0

        let mut tree: AvlTree<i32> = AvlTree::new();
        util_insert_sequence(&mut tree, vec![8, 4, 12, 2, 10, 14, 9, 13]);
        util_verify_tree(&tree, 8, vec![
            (8, Some(4), Some(12), 1),
            (4, Some(2), None, -1),
            (12, Some(10), Some(14), 0),
            (2, None, None, 0),
            (10, Some(9), None, -1),
            (14, Some(13), None, -1),
            (9, None, None, 0),
            (13, None, None, 0),
        ]);
        
        tree.delete(4);

        util_verify_tree(&tree, 12, vec![
            (12, Some(8), Some(14), -1),
            (8, Some(2), Some(10), 1),
            (14, Some(13), None, -1),
            (2, None, None, 0),
            (10, Some(9), None, -1),
            (13, None, None, 0),
            (9, None, None, 0),
        ]);
    }

    #[test]
    fn test_rotate_left_insert_000() {
        // only case for insertion: right balance factor must be 1
        
        let mut tree: AvlTree<i32> = AvlTree::new();
        util_insert_sequence(&mut tree, vec![8, 4, 12, 10, 14]);
        util_verify_tree(&tree, 8, vec![
            (4, None, None, 0),
            (8, Some(4), Some(12), 1),
            (10, None, None, 0),
            (12, Some(10), Some(14), 0),
            (14, None, None, 0),
        ]);

        tree.insert(13);

        util_verify_tree(&tree, 12, vec![
            (4, None, None, 0),
            (8, Some(4), Some(10), 0),
            (10, None, None, 0),
            (12, Some(8), Some(14), 0),
            (13, None, None, 0),
            (14, Some(13), None, -1),
        ]);
    }

    #[test]
    fn test_rotate_leftright_delete_000() {
        // case 1/3 for deletion: left-right balance factor = -1

        let mut tree: AvlTree<i32> = AvlTree::new();
        util_insert_sequence(&mut tree, vec![8, 4, 12, 2, 6, 14, 5]);
        util_verify_tree(&tree, 8, vec![
            (8, Some(4), Some(12), -1),
            (4, Some(2), Some(6), 1),
            (12, None, Some(14), 1),
            (2, None, None, 0),
            (6, Some(5), None, -1),
            (14, None, None, 0),
            (5, None, None, 0),
        ]);

        tree.delete(12);
        
        util_verify_tree(&tree, 6, vec![
            (6, Some(4), Some(8), 0),
            (4, Some(2), Some(5), 0),
            (8, None, Some(14), 1),
            (2, None, None, 0),
            (5, None, None, 0),
            (14, None, None, 0),
        ]);
    }

    #[test]
    fn test_rotate_leftright_delete_001() {
        // case 2/3 for deletion: left-right balance factor = 0

        let mut tree: AvlTree<i32> = AvlTree::new();
        util_insert_sequence(&mut tree, vec![8, 4, 12, 2, 6, 14, 5, 7]);
        util_verify_tree(&tree, 8, vec![
            (8, Some(4), Some(12), -1),
            (4, Some(2), Some(6), 1),
            (12, None, Some(14), 1),
            (2, None, None, 0),
            (6, Some(5), Some(7), 0),
            (14, None, None, 0),
            (5, None, None, 0),
            (7, None, None, 0),
        ]);

        tree.delete(14);

        util_verify_tree(&tree, 6, vec![
            (6, Some(4), Some(8), 0),
            (4, Some(2), Some(5), 0),
            (8, Some(7), Some(12), 0),
            (2, None, None, 0),
            (5, None, None, 0),
            (7, None, None, 0),
            (12, None, None, 0),
        ]);
    }

    #[test]
    fn test_rotate_leftright_delete_002() {
        // case 3/3 for deletion: left-right balance factor = 1

        let mut tree: AvlTree<i32> = AvlTree::new();
        util_insert_sequence(&mut tree, vec![8, 4, 12, 2, 6, 10, 7]);
        util_verify_tree(&tree, 8, vec![
            (8, Some(4), Some(12), -1),
            (4, Some(2), Some(6), 1),
            (12, Some(10), None, -1),
            (2, None, None, 0),
            (6, None, Some(7), 1),
            (10, None, None, 0),
            (7, None, None, 0),
        ]);

        tree.delete(12);

        util_verify_tree(&tree, 6, vec![
            (6, Some(4), Some(8), 0),
            (4, Some(2), None, -1),
            (8, Some(7), Some(10), 0),
            (2, None, None, 0),
            (7, None, None, 0),
            (10, None, None, 0),
        ]);
    }

    #[test]
    fn test_rotate_leftright_insert_000() {
        // case 1/2 for insertion: left-right balance factor = -1

        let mut tree: AvlTree<i32> = AvlTree::new();
        util_insert_sequence(&mut tree, vec![8, 4, 12, 2, 6]);
        util_verify_tree(&tree, 8, vec![
            (2, None, None, 0),
            (4, Some(2), Some(6), 0),
            (6, None, None, 0),
            (8, Some(4), Some(12), -1),
            (12, None, None, 0),
        ]);

        tree.insert(5);

        util_verify_tree(&tree, 6, vec![
            (2, None, None, 0),
            (4, Some(2), Some(5), 0),
            (5, None, None, 0),
            (6, Some(4), Some(8), 0),
            (8, None, Some(12), 1),
            (12, None, None, 0),
        ]);
    }

    #[test]
    fn test_rotate_leftright_insert_001() {
        // case 2/2 for insertion: left-right balance factor = 1

        let mut tree: AvlTree<i32> = AvlTree::new();
        util_insert_sequence(&mut tree, vec![8, 4, 12, 2, 6]);
        util_verify_tree(&tree, 8, vec![
            (2, None, None, 0),
            (4, Some(2), Some(6), 0),
            (6, None, None, 0),
            (8, Some(4), Some(12), -1),
            (12, None, None, 0),
        ]);

        tree.insert(7);

        util_verify_tree(&tree, 6, vec![
            (2, None, None, 0),
            (4, Some(2), None, -1),
            (6, Some(4), Some(8), 0),
            (7, None, None, 0),
            (8, Some(7), Some(12), 0),
            (12, None, None, 0),
        ]);
    }

    #[test]
    fn test_rotate_right_delete_000() {
        // case 1/2 for deletion: left balance factor = -1

        let mut tree: AvlTree<i32> = AvlTree::new();
        util_insert_sequence(&mut tree, vec![8, 4, 12, 2, 6, 14, 3]);
        util_verify_tree(&tree, 8, vec![
            (8, Some(4), Some(12), -1),
            (4, Some(2), Some(6), -1),
            (12, None, Some(14), 1),
            (2, None, Some(3), 1),
            (6, None, None, 0),
            (14, None, None, 0),
            (3, None, None, 0),
        ]);

        tree.delete(12);

        util_verify_tree(&tree, 4, vec![
            (4, Some(2), Some(8), 0),
            (2, None, Some(3), 1),
            (8, Some(6), Some(14), 0),
            (3, None, None, 0),
            (6, None, None, 0),
            (14, None, None, 0),
        ]);
    }

    #[test]
    fn test_rotate_right_delete_001() {
        // case 2/2 for deletion: right balance factor = 0

        let mut tree: AvlTree<i32> = AvlTree::new();
        util_insert_sequence(&mut tree, vec![8, 4, 12, 2, 6, 14, 3, 7]);
        util_verify_tree(&tree, 8, vec![
            (8, Some(4), Some(12), -1),
            (4, Some(2), Some(6), 0),
            (12, None, Some(14), 1),
            (2, None, Some(3), 1),
            (6, None, Some(7), 1),
            (14, None, None, 0),
            (3, None, None, 0),
            (7, None, None, 0),
        ]);
        
        tree.delete(12);

        util_verify_tree(&tree, 4, vec![
            (4, Some(2), Some(8), 1),
            (2, None, Some(3), 1),
            (8, Some(6), Some(14), -1),
            (3, None, None, 0),
            (6, None, Some(7), 1),
            (14, None, None, 0),
            (7, None, None, 0),
        ]);
    }

    #[test]
    fn test_rotate_right_insert_000() {
        // only case for insertion: left balance factor must be -1

        let mut tree: AvlTree<i32> = AvlTree::new();
        util_insert_sequence(&mut tree, vec![8, 4, 12, 2, 6]);
        util_verify_tree(&tree, 8, vec![
            (2, None, None, 0),
            (4, Some(2), Some(6), 0),
            (6, None, None, 0),
            (8, Some(4), Some(12), -1),
            (12, None, None, 0),
        ]);

        tree.insert(3);

        util_verify_tree(&tree, 4, vec![
            (2, None, Some(3), 1),
            (3, None, None, 0),
            (4, Some(2), Some(8), 0),
            (6, None, None, 0),
            (8, Some(6), Some(12), 0),
            (12, None, None, 0),
        ]);
    }

    #[test]
    fn test_rotate_rightleft_delete_000() {
        // case 1/3 for deletion: right-left balance factor = 1

        let mut tree: AvlTree<i32> = AvlTree::new();
        util_insert_sequence(&mut tree, vec![8, 4, 12, 2, 10, 14, 11]);
        util_verify_tree(&tree, 8, vec![
            (8, Some(4), Some(12), 1),
            (4, Some(2), None, -1),
            (12, Some(10), Some(14), -1),
            (2, None, None, 0),
            (10, None, Some(11), 1),
            (14, None, None, 0),
            (11, None, None, 0),
        ]);

        tree.delete(4);

        util_verify_tree(&tree, 10, vec![
            (10, Some(8), Some(12), 0),
            (8, Some(2), None, -1),
            (12, Some(11), Some(14), 0),
            (2, None, None, 0),
            (11, None, None, 0),
            (14, None, None, 0),
        ]);
    }

    #[test]
    fn test_rotate_rightleft_delete_001() {
        // case 2/3 for deletion: right-left balance factor = 0

        let mut tree: AvlTree<i32> = AvlTree::new();
        util_insert_sequence(&mut tree, vec![8, 4, 12, 2, 10, 14, 9, 11]);
        util_verify_tree(&tree, 8, vec![
            (8, Some(4), Some(12), 1),
            (4, Some(2), None, -1),
            (12, Some(10), Some(14), -1),
            (2, None, None, 0),
            (10, Some(9), Some(11), 0),
            (14, None, None, 0),
            (9, None, None, 0),
            (11, None, None, 0),
        ]);

        tree.delete(2);

        util_verify_tree(&tree, 10, vec![
            (10, Some(8), Some(12), 0),
            (8, Some(4), Some(9), 0),
            (12, Some(11), Some(14), 0),
            (4, None, None, 0),
            (9, None, None, 0),
            (11, None, None, 0),
            (14, None, None, 0),
        ]);
    }

    #[test]
    fn test_rotate_rightleft_delete_002() {
        // case 3/3 for deletion: right-left balance factor = -1

        let mut tree: AvlTree<i32> = AvlTree::new();
        util_insert_sequence(&mut tree, vec![8, 4, 12, 6, 10, 14, 9]);
        util_verify_tree(&tree, 8, vec![
            (8, Some(4), Some(12), 1),
            (4, None, Some(6), 1),
            (12, Some(10), Some(14), -1),
            (6, None, None, 0),
            (10, Some(9), None, -1),
            (14, None, None, 0),
            (9, None, None, 0),
        ]);

        tree.delete(4);

        util_verify_tree(&tree, 10, vec![
            (10, Some(8), Some(12), 0),
            (8, Some(6), Some(9), 0),
            (12, None, Some(14), 1),
            (6, None, None, 0),
            (9, None, None, 0),
            (14, None, None, 0),
        ]);
    }

    #[test]
    fn test_rotate_rightleft_insert_000() {
        // case 1/2 for insertion: right-left balance factor = 1
        
        let mut tree: AvlTree<i32> = AvlTree::new();
        util_insert_sequence(&mut tree, vec![8, 4, 12, 10, 14]);
        util_verify_tree(&tree, 8, vec![
            (4, None, None, 0),
            (8, Some(4), Some(12), 1),
            (10, None, None, 0),
            (12, Some(10), Some(14), 0),
            (14, None, None, 0),
        ]);

        tree.insert(11);

        util_verify_tree(&tree, 10, vec![
            (4, None, None, 0),
            (8, Some(4), None, -1),
            (10, Some(8), Some(12), 0),
            (11, None, None, 0),
            (12, Some(11), Some(14), 0),
            (14, None, None, 0),
        ]);
    }

    #[test]
    fn test_rotate_rightleft_insert_001() {
        // case 2/2 for insertion: right-left balance factor = -1
        
        let mut tree: AvlTree<i32> = AvlTree::new();
        util_insert_sequence(&mut tree, vec![8, 4, 12, 10, 14]);
        util_verify_tree(&tree, 8, vec![
            (4, None, None, 0),
            (8, Some(4), Some(12), 1),
            (10, None, None, 0),
            (12, Some(10), Some(14), 0),
            (14, None, None, 0),
        ]);

        tree.insert(9);

        util_verify_tree(&tree, 10, vec![
            (4, None, None, 0),
            (8, Some(4), Some(9), 0),
            (9, None, None, 0),
            (10, Some(8), Some(12), 0),
            (12, None, Some(14), 1),
            (14, None, None, 0),
        ]);
    }

}
