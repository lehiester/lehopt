use crate::primitives::int::IntNonnegFitsInUsize;
use crate::primitives::vector::DenseVector;
use crate::primitives::vector::SparseVector;
use crate::primitives::scalar::Scalar;

use std::collections::BinaryHeap;



/// A matrix in compressed sparse column format.
/// 
/// The invariants must be very carefully maintained because this is intended
/// for use with unsafe code which relies on them without checking.
/// 
/// Invariants:
/// - 01: positive row count
/// - 02: positive column count
/// - 03: number of columns must be at least 1 less than max for int type
/// - 04: column boundaries array is (num_cols + 1) long
/// - 05: values and row indices are same length
/// - 06: last "column boundary" is length of values
/// - 07: value indices fit within max for the given int type
/// - 08: column boundaries in range of value indices (0 to length, inclusive)
/// - 09: column boundaries are nondecreasing
/// - 10: row indices in range of row count (0 to num_rows - 1)
/// - 11: row indices are sorted ascending for each column
/// - 12: row indices have no duplicates in a given column
/// 
/// Invariant 07 is implied by invariant 06 in conjunction with the struct
/// definition but is listed separately for cases in which it's easier to check
/// it first, and as a reminder to check for integer overflow.
pub struct CSCMatrix<S=f64,I=usize> where S: Scalar, I: IntNonnegFitsInUsize {
    num_rows: I,
    num_cols: I,
    col_boundaries: Vec<I>,
    row_indices: Vec<I>,
    values: Vec<S>,
}

impl<S,I> CSCMatrix<S,I> where S: Scalar, I: IntNonnegFitsInUsize {

    // Constructors

    pub fn new(num_rows: I, num_cols: I, col_boundaries: Vec<I>, row_indices: Vec<I>, values: Vec<S>) -> Self {
        // validate invariants 01 and 02: positive row and column counts
        assert!(num_rows > I::ZERO);
        assert!(num_cols > I::ZERO);

        // validate invariant 03: number of columns must be at least 1 less than max for int type
        assert!(num_cols < I::MAX);

        // validate invariant 04: column boundaries array is (num_cols + 1) long
        assert_eq!(col_boundaries.len(), (num_cols + I::ONE).to_usize_unchecked());

        // validate invariant 05: values and row indices are same length
        let value_len = values.len();
        assert_eq!(value_len, row_indices.len());

        // validate invariant 06: last "column boundary" is length of values
        assert_eq!(col_boundaries[col_boundaries.len() - 1].to_usize_unchecked(), value_len);

        // validate invariant 07: value indices fit within max for the given int type
        // (redundant with current implementation of checking last "column boundary",
        // but included here as a guardrail for future changes, and also as a reminder
        // for other constructors which may need to explicitly check this)
        assert!(value_len <= I::MAX.to_usize_unchecked());
        let value_len: I = value_len.try_into().unwrap();

        // validate invariant 08: column boundaries in range of value indices (0 to length, inclusive)
        assert!(col_boundaries.iter().all(|x| *x >= I::ZERO));
        assert!(col_boundaries.iter().all(|x| *x <= value_len));

        // validate invariant 09: column boundaries are nondecreasing
        for c in 1..col_boundaries.len() {
            assert!(col_boundaries[c] >= col_boundaries[c - 1]);
        }

        // validate invariant 10: row indices in range of row count (0 to num_rows - 1)
        assert!(row_indices.iter().all(|x| *x >= I::ZERO));
        assert!(row_indices.iter().all(|x| *x < num_rows));

        // validate invariants 11 and 12: row indices are ascending with no duplicates in each column
        for c in 0..(col_boundaries.len() - 1) {
            let i_start = col_boundaries[c].to_usize_unchecked() + 1;
            let i_end = col_boundaries[c + 1].to_usize_unchecked();
            for i in i_start..i_end {
                assert!(row_indices[i] > row_indices[i - 1]);
            }
        }

        Self {
            num_rows,
            num_cols,
            col_boundaries,
            row_indices,
            values,
        }
    }

    pub fn new_from_columns(columns: Vec<SparseVector<S,I>>) -> Self {
        // validate invariants 01 and 02: positive row and column counts
        assert!(columns.len() > 0);
        assert!(columns[0].length() > I::ZERO);

        // validate invariant 03: number of columns must be at least 1 less than max for int type
        assert!(columns.len() < I::MAX.to_usize_unchecked());

        let num_cols = columns.len().try_into().unwrap();
        let num_rows = columns[0].length();

        // validate inputs: columns must all be of the same length
        for c in 1..columns.len() {
            assert_eq!(columns[c].length(), num_rows);
        }

        // validate invariant 07: value indices fit within max for the given int type
        let mut total_values: usize = 0;
        for col in &columns {
            total_values = total_values.checked_add(col.values().len().to_usize_unchecked()).unwrap();  // panic on overflow
        }
        assert!(total_values <= I::MAX.to_usize_unchecked());

        // pack the columns into CSC format
        let mut col_boundaries: Vec<I> = vec![I::ZERO];
        let mut row_indices: Vec<I> = vec![];
        let mut values: Vec<S> = vec![];
        let mut element_counter = I::ZERO;
        for col in columns {
            // `col.values` length fits within `I` per SparseVector invariant 03;
            // overflow already checked above
            element_counter += col.values().len().try_into().unwrap();
            col_boundaries.push(element_counter);

            // destructively move indices and values out of `col` (empty afterward)
            let (mut col_indices, mut col_values) = col.to_vecs();
            row_indices.append(&mut col_indices);
            values.append(&mut col_values);
        }

        // sanity checks for invariants 05 and 06
        assert_eq!(col_boundaries.last().unwrap().to_usize_unchecked(), total_values);
        assert_eq!(row_indices.len(), total_values);
        assert_eq!(values.len(), total_values);

        // invariants satisfied by construction:
        // - 04: column boundaries array is (num_cols + 1) long
        // - 06: last "column boundary" is length of values
        // - 08: column boundaries in range of value indices (0 to length, inclusive)
        // - 09: column boundaries are nondecreasing

        // invariants satisfied by construction and SparseVector invariants:
        // - 05: values and row indices are same length
        // - 10: row indices in range of row count (0 to num_rows - 1)
        // - 11: row indices are sorted ascending for each column
        // - 12: row indices have no duplicates in a given column

        Self {
            num_rows,
            num_cols,
            col_boundaries,
            row_indices,
            values,
        }
    }

    // Static methods

    fn check_num_cols<I2,I3>(num_cols_before: I2, num_new_cols: I3) -> I
        where I2: IntNonnegFitsInUsize, I3: IntNonnegFitsInUsize {
        
        let num_cols_after: usize = num_cols_before.to_usize_unchecked()
            .checked_add(num_new_cols.to_usize_unchecked())
            .expect("Too many matrix columns for this architecture!");
        let num_cols_after: I = num_cols_after.try_into()
            .expect(format!("Too many matrix columns; max is {}, attempted {}", I::MAX - I::ONE, num_cols_after).as_str());
        assert!(num_cols_after <= I::MAX - I::ONE, "Too many matrix columns; max is {}, attempted {}", I::MAX - I::ONE, num_cols_after);
        num_cols_after
    }

    fn check_num_rows<I2,I3>(&self, num_rows_before: I2, num_new_rows: I3) -> I
        where I2: IntNonnegFitsInUsize, I3: IntNonnegFitsInUsize {

        let num_rows_after = num_rows_before.to_usize_unchecked()
            .checked_add(num_new_rows.to_usize_unchecked())
            .expect("Too many matrix rows for this architecture!");
        let num_rows_after: I = num_rows_after.try_into()
            .expect(format!("Too many matrix rows; max is {}, attempted {}", I::MAX, num_rows_after).as_str());
        num_rows_after
    }

    fn check_num_values<C,I2,I3>(num_values_before: I2, new_stuff: C) -> I
        where C: Iterator<Item=I3>, I2: IntNonnegFitsInUsize, I3: IntNonnegFitsInUsize {
        
        let mut total: usize = num_values_before.to_usize_unchecked();
        for size in new_stuff {
            total = total.checked_add(size.to_usize_unchecked())
                .expect(format!("Too many matrix nonzeros for this architecture!").as_str());
        }
        let total: I = total.try_into()
            .expect(format!("Too many matrix nonzeros; max is {}, attempted {}", I::MAX, total).as_str());
        total
    }

    // Getters

    #[inline]
    pub fn col_boundaries(&self) -> &Vec<I> {
        &self.col_boundaries
    }

    #[inline]
    pub fn num_cols(&self) -> I {
        self.num_cols
    }

    #[inline]
    pub fn num_rows(&self) -> I {
        self.num_rows
    }

    #[inline]
    pub fn row_indices(&self) -> &Vec<I> {
        &self.row_indices
    }

    #[inline]
    pub fn values(&self) -> &Vec<S> {
        &self.values
    }

    // Instance methods

    pub fn append_cols_inplace(&mut self, cols: Vec<SparseVector<S,I>>) {
        let num_cols_after: I = Self::check_num_cols(self.num_cols, cols.len());
        let _ = Self::check_num_values(
            self.values.len(),
            cols.iter().map(|col| col.values().len())
        );

        let mut num_values: usize = self.values.len();
        for col in cols {
            assert_eq!(col.length(), self.num_rows);

            num_values += col.indices().len();
            self.col_boundaries.push(num_values.try_into().unwrap());

            // move values out of SparseVector (destructively)
            let (mut col_indices, mut col_values) = col.to_vecs();
            self.row_indices.append(&mut col_indices);
            self.values.append(&mut col_values);
        }
        
        self.num_cols = num_cols_after;
    }

    pub fn append_col_multiples_inplace(&mut self, col_indices: Vec<I>, multiplier: S) {
        // invariants satisfied by existing `self`: 01, 02, 10, 11, 12
        // invariants satisfied by construction: 04, 05, 06, 08, 09

        // validate inputs: `col_indices` must all be in range
        for c in &col_indices {
            assert!(*c >= I::ZERO);
            assert!(*c < self.num_cols);
        }

        // validate invariant 03: number of columns must be at least 1 less than max for int type
        let new_num_cols = Self::check_num_cols(self.num_cols, col_indices.len());

        // validate invariant 07: value indices fit within max for the given int type
        Self::check_num_values(self.values.len(), col_indices.iter().map(|c| self.values_in_col(*c)));

        // append a multiple of each designated column
        let mut total_values = self.values.len();
        for c in &col_indices {
            let c_usize = c.to_usize_unchecked();
            let col_start = self.col_boundaries[c_usize].to_usize_unchecked();
            let col_end = self.col_boundaries[c_usize + 1].to_usize_unchecked();
            for i in col_start..col_end {
                self.row_indices.push(self.row_indices[i]);
                self.values.push(self.values[i] * multiplier);
            }
            total_values += col_end - col_start;
            self.col_boundaries.push(total_values.try_into().unwrap());
        }

        self.num_cols = new_num_cols;
    }

    pub fn append_rows_inplace(&mut self, rows: Vec<SparseVector<S,I>>) {
        // invariants satisfied by existing `self`: 01, 02, 03, 04
        // invariants satisfied by construction: 05, 06, 08, 09, 10, 11
        // invariant 12 satisfied due to SparseVector invariant 05 (no duplicate indices)

        // IDEA: test speed/memory tradeoff--this low-memory-footprint implementation
        //       doesn't add a lot of runtime overhead, right? or does it?

        // TODO: unit tests for maintaining invariants
        // TODO: unit tests with interleaved rows to test the heap management really well

        // validate inputs: row count cannot exceed capacity of int type `I`
        let num_rows_before = self.num_rows;
        let num_rows_after = self.check_num_rows(self.num_rows, rows.len());
        self.num_rows = num_rows_after;

        // validate invariant 07: value indices fit within max for the given int type
        let _ = Self::check_num_values(
            self.values.len(),
            rows.iter().map(|row| row.values().len())
        );

        // extend `row_indices` and `values` to hold the new items
        let num_elems_before = self.values.len();
        let num_new_elems: usize = rows.iter().map(|row| row.values().len()).sum();
        for _ in 0..num_new_elems {
            self.row_indices.push(I::ZERO);
            self.values.push(S::ZERO);
        }

        // construct heap of "next" column for each row, columns in reverse order:
        // elements are (col_idx, row_idx, row_indices_idx)
        let mut next_elem_by_row: BinaryHeap<(usize,usize,usize)> = BinaryHeap::new();  // this is a max heap!
        for r in 0..rows.len() {
            let row = &rows[r];
            let len = row.indices().len();
            if len > 0 {
                next_elem_by_row.push((row.indices()[len - 1].to_usize_unchecked(), r, len - 1));
            }
        }

        // update and shift columns in reverse order
        let mut next_elem = next_elem_by_row.pop();
        let mut next_csc_idx = num_elems_before + num_new_elems - 1;
        for c in (0..self.num_cols.to_usize_unchecked()).rev() {
            let col_end_new = next_csc_idx + 1;

            // add new rows for this column
            while next_elem.is_some_and(|t| t.0 == c) {
                // add next element to the CSC structure
                let (_, r, r_ind_idx) = next_elem.unwrap();
                self.row_indices[next_csc_idx] = num_rows_before + r.try_into().unwrap();
                self.values[next_csc_idx] = rows[r].values()[r_ind_idx];
                next_csc_idx = next_csc_idx.wrapping_sub(1);
                // add the next element in this row into the heap
                if r_ind_idx > 0 {
                    let next_c = rows[r].indices()[r_ind_idx - 1].to_usize_unchecked();
                    next_elem_by_row.push((next_c, r, r_ind_idx - 1));
                }
                // pop next element/row off of the heap
                next_elem = next_elem_by_row.pop();
            }

            // once no new elements remain, all existing elements are where they should be
            if next_elem.is_none() {
                self.col_boundaries[c+1] = col_end_new.try_into().unwrap();
                break;
            }

            // shift existing rows for this column (note, column pointers are
            // not updated yet and still point to existing elements)
            let col_start = self.col_boundaries[c].to_usize_unchecked();
            let col_end = self.col_boundaries[c+1].to_usize_unchecked();
            for i in (col_start..col_end).rev() {
                self.row_indices[next_csc_idx] = self.row_indices[i];
                self.values[next_csc_idx] = self.values[i];
                next_csc_idx = next_csc_idx.wrapping_sub(1);
            }

            // update column pointer (do not update `col_boundaries[c]` yet,
            // the next loop iteration needs it and will update it afterward)
            self.col_boundaries[c+1] = col_end_new.try_into().unwrap();
        }
    }

    pub fn col_dense(&self, col_idx: I) -> DenseVector<S,I> {
        assert!(col_idx >= I::ZERO);

        let c_usize = col_idx.to_usize_unchecked();
        let col_start = self.col_boundaries[c_usize].to_usize_unchecked();
        let col_end = self.col_boundaries[c_usize+1].to_usize_unchecked();

        let indices: Vec<I> = self.row_indices[col_start..col_end].to_vec();
        let values: Vec<S> = self.values[col_start..col_end].to_vec();

        DenseVector::new_from_sparse_lists(self.num_rows, indices, values)
    }

    pub fn col_dot(&self, col_idx: I, vec: &DenseVector<S,I>) -> S {
        assert_eq!(vec.len(), self.num_rows);
        assert!(col_idx >= I::ZERO);

        let c_usize = col_idx.to_usize_unchecked();
        let col_start = self.col_boundaries[c_usize].to_usize_unchecked();
        let col_end = self.col_boundaries[c_usize + 1].to_usize_unchecked();

        let mut total = S::ZERO;
        for i in col_start..col_end {
            let r = self.row_indices[i].to_usize_unchecked();
            let value = self.values[i];
            total += value * vec.values()[r];
        }

        total
    }

    pub fn col_submatrix(&self, col_indices: &Vec<I>) -> Self {
        // IDEA: strict variant that requires `col_indices` to be a subset (all invariants satisfied)
        // IDEA (cont.): ^ checked version with HashSet

        // invariants satisfied by `self`: 01, 10, 11, 12
        // invariants satisfied by construction: 04, 05, 06, 08, 09

        // validate invariants 02, 03, 07
        assert!(col_indices.len() > 0);
        let num_cols = Self::check_num_cols(0, col_indices.len());
        let _ = Self::check_num_values(0, col_indices.iter().map(|c| self.values_in_col(*c)));

        let mut value_counter: I = I::ZERO;
        let mut col_boundaries: Vec<I> = vec![I::ZERO];
        let mut row_indices: Vec<I> = vec![];
        let mut values: Vec<S> = vec![];

        for c in col_indices {
            let c_usize = c.to_usize_unchecked();
            let c_start = self.col_boundaries[c_usize];
            let c_end = self.col_boundaries[c_usize + 1];
            let c_start_usize = c_start.to_usize_unchecked();
            let c_end_usize = c_end.to_usize_unchecked();
            row_indices.extend((c_start_usize..c_end_usize).map(|i| self.row_indices[i]));
            values.extend((c_start_usize..c_end_usize).map(|i| self.values[i]));
            value_counter += c_end - c_start;
            col_boundaries.push(value_counter);
        }

        Self {
            num_rows: self.num_rows,
            num_cols,
            col_boundaries,
            row_indices,
            values,
        }
    }

    pub fn pop_cols_inplace<I2>(&mut self, num_to_pop: I2)
        where I2: IntNonnegFitsInUsize {
        
        // invariants satisfied by `self`: 01, 03, 07, 09, 10, 11, 12
        // invariants satisfied by construction: 04, 05, 06, 08

        // assumption below: at least one col popped
        assert!(num_to_pop >= I2::ZERO, "CSCMatrix.pop_cols_inplace got negative num_to_pop");
        if num_to_pop == I2::ZERO {
            return;
        }

        // validate invariant 02: positive column count
        let num_cols_usize = self.num_cols.to_usize_unchecked();
        let num_to_pop_usize = num_to_pop.to_usize_unchecked();
        assert!(num_to_pop_usize < num_cols_usize, "CSCMatrix.pop_cols_inplace can't pop all columns--zero columns not allowed!");
        
        // truncate all arrays appropriately (does not reallocate!)
        let num_to_keep = num_cols_usize - num_to_pop_usize;
        let values_to_keep = self.col_boundaries[num_to_keep].to_usize_unchecked();
        self.col_boundaries.truncate(num_to_keep + 1);
        self.row_indices.truncate(values_to_keep);
        self.values.truncate(values_to_keep);

        self.num_cols = num_to_keep.try_into().unwrap();
    }

    pub fn pop_rows_inplace<I2>(&mut self, num_to_pop: I2)
        where I2: IntNonnegFitsInUsize {

        // invariants satisfied by `self`: 02, 03, 04, 07, 11, 12
        // invariants satisfied by construction: 05, 06, 08, 09, 10

        // assumption below: at least one row popped
        assert!(num_to_pop >= I2::ZERO, "CSCMatrix.pop_rows_inplace got negative num_to_pop");
        if num_to_pop == I2::ZERO {
            return;
        }
        
        // validate invariant 01: positive row count
        let num_rows_usize = self.num_rows.to_usize_unchecked();
        let num_to_pop_usize = num_to_pop.to_usize_unchecked();
        assert!(num_to_pop_usize < num_rows_usize, "CSCMatrix.pop_rows_inplace can't pop all rows--zero rows not allowed!");
        let num_to_keep: I = (num_rows_usize - num_to_pop_usize).try_into().unwrap();

        let mut next_idx = 0;
        for c in 0..self.num_cols.to_usize_unchecked() {
            let c_start = self.col_boundaries[c].to_usize_unchecked();
            let c_end = self.col_boundaries[c+1].to_usize_unchecked();

            // leave unaffected columns at beginning untouched (relies on invariant 11!)
            if next_idx == c_start {
                if c_start == c_end {
                    continue;
                }
                else if self.row_indices[c_end - 1] < num_to_keep {
                    next_idx = c_end;
                    continue;
                }
            }

            // shift elements for affected columns
            self.col_boundaries[c] = next_idx.try_into().unwrap();  // already have previous value: c_start
            for i in c_start..c_end {
                if self.row_indices[i] < num_to_keep {
                    self.row_indices[next_idx] = self.row_indices[i];
                    self.values[next_idx] = self.values[i];
                    next_idx += 1;
                }
                else {
                    break;
                }
            }
        }

        // update last col_boundary (invariant 06) and number of rows
        self.col_boundaries[self.num_cols.to_usize_unchecked()] = next_idx.try_into().unwrap();
        self.num_rows = num_to_keep;

        // truncate row_indices and values (this does not reallocate!)
        self.row_indices.truncate(next_idx);
        self.values.truncate(next_idx);

        // TODO: unit testing
    }

    pub fn remove_cols_inplace<I2>(&mut self, col_indices: &Vec<I2>)
        where I2: IntNonnegFitsInUsize {

        // at least one column assumed below
        if col_indices.len() == 0 {
            return;
        }

        let mut col_indices = col_indices.clone();
        col_indices.sort();

        // validate arguments: col_indices must be in bounds, no duplicates, can't be all columns
        assert!(col_indices[0] >= I2::ZERO);
        assert!(col_indices[col_indices.len() - 1].to_usize_unchecked() < self.num_cols.to_usize_unchecked());
        for i in 1..col_indices.len() {
            assert_ne!(col_indices[i], col_indices[i-1], "CSCMatrix.remove_cols_inplace: got duplicate indices: {}", col_indices[i]);
        }
        assert!(col_indices.len() < self.num_cols.to_usize_unchecked(), "CSCMatrix.remove_cols_inplace: can't remove all columns!");

        // remove and shift (first column that may need to be shifted is
        // immediately after the first removed)
        let num_cols_usize = self.num_cols.to_usize_unchecked();
        let c_first_removed = col_indices[0].to_usize_unchecked();
        let mut col_indices_idx_next = 1;
        let mut next_val_idx = self.col_boundaries[col_indices[0].to_usize_unchecked()].to_usize_unchecked();
        let mut next_col_idx = c_first_removed;
        for c in (c_first_removed + 1)..num_cols_usize {
            // skip over each column that should be removed
            if col_indices_idx_next != col_indices.len() && c == col_indices[col_indices_idx_next].to_usize_unchecked() {
                col_indices_idx_next += 1;  // now watch for the next one
                continue;
            }

            // shift this column
            let col_start = self.col_boundaries[c].to_usize_unchecked();
            let col_end = self.col_boundaries[c+1].to_usize_unchecked();
            for i in col_start..col_end {
                self.row_indices[next_val_idx] = self.row_indices[i];
                self.values[next_val_idx] = self.values[i];
                next_val_idx += 1;
            }
            next_col_idx += 1;
            self.col_boundaries[next_col_idx] = next_val_idx.try_into().unwrap();
        }

        // truncate (this does not resize/reallocate!)
        self.col_boundaries.truncate(next_col_idx + 1);
        self.row_indices.truncate(next_val_idx);
        self.values.truncate(next_val_idx);

        self.num_cols = next_col_idx.try_into().unwrap();
    }

    pub fn scale_col_inplace<I2>(&mut self, col_index: I2, multiplier: S)
        where I2: IntNonnegFitsInUsize {
        
        // all invariants unaffected
        let c = col_index.to_usize_unchecked();
        let col_start = self.col_boundaries[c].to_usize_unchecked();
        let col_end = self.col_boundaries[c+1].to_usize_unchecked();
        for i in col_start..col_end {
            self.values[i] *= multiplier;
        }
    }

    fn values_in_col(&self, c: I) -> I {
        let c_usize = c.to_usize_unchecked();
        self.col_boundaries[c_usize + 1] - self.col_boundaries[c_usize]
    }

}

impl<S,I> Clone for CSCMatrix<S,I> where S: Scalar, I: IntNonnegFitsInUsize {

    fn clone(&self) -> Self {
        Self {
            num_rows: self.num_rows,
            num_cols: self.num_cols,
            col_boundaries: self.col_boundaries.clone(),
            row_indices: self.row_indices.clone(),
            values: self.values.clone(),
        }
    }

}



#[cfg(test)]
mod tests {

    use super::*;

    #[test]
    fn test_append_rows_inplace_000() {
        let mat: CSCMatrix = CSCMatrix::new(
            1,
            5,
            vec![0, 0, 1, 1, 2, 2],
            vec![0, 0],
            vec![1.0, 2.0]
        );

        let mut m1 = mat.clone();
        m1.append_rows_inplace(vec![
            SparseVector::new(5, vec![2], vec![3.0])
        ]);
        assert_eq!(m1.num_rows, 2);
        assert_eq!(m1.num_cols, 5);
        assert_eq!(m1.col_boundaries, vec![0, 0, 1, 2, 3, 3]);
        assert_eq!(m1.row_indices, vec![0, 1, 0]);
        assert_eq!(m1.values, vec![1.0, 3.0, 2.0]);

        let mut m2 = m1.clone();
        m2.append_rows_inplace(vec![
            SparseVector::new(5, vec![4], vec![4.0]),
            SparseVector::new(5, vec![4], vec![5.0])
        ]);
        assert_eq!(m2.num_rows, 4);
        assert_eq!(m2.num_cols, 5);
        assert_eq!(m2.col_boundaries, vec![0, 0, 1, 2, 3, 5]);
        assert_eq!(m2.row_indices, vec![0, 1, 0, 2, 3]);
        assert_eq!(m2.values, vec![1.0, 3.0, 2.0, 4.0, 5.0]);

        let mut m3 = m2.clone();
        m3.append_rows_inplace(vec![
            SparseVector::new(5, vec![0], vec![6.0])
        ]);
        assert_eq!(m3.num_rows, 5);
        assert_eq!(m3.num_cols, 5);
        assert_eq!(m3.col_boundaries, vec![0, 1, 2, 3, 4, 6]);
        assert_eq!(m3.row_indices, vec![4, 0, 1, 0, 2, 3]);
        assert_eq!(m3.values, vec![6.0, 1.0, 3.0, 2.0, 4.0, 5.0]);

        let mut m4 = m3.clone();
        m4.append_rows_inplace(vec![
            SparseVector::new(5, vec![0, 1, 2, 3, 4], vec![9.0, 9.0, 9.0, 9.0, 9.0])
        ]);
        assert_eq!(m4.num_rows, 6);
        assert_eq!(m4.num_cols, 5);
        assert_eq!(m4.col_boundaries, vec![0, 2, 4, 6, 8, 11]);
        assert_eq!(m4.row_indices, vec![4, 5, 0, 5, 1, 5, 0, 5, 2, 3, 5]);
        assert_eq!(m4.values, vec![6.0, 9.0, 1.0, 9.0, 3.0, 9.0, 2.0, 9.0, 4.0, 5.0, 9.0]);
    }

    #[test]
    fn test_col_dense_000() {
        let mat: CSCMatrix = CSCMatrix::new(
            5,
            5,
            vec![0, 2, 3, 5, 5, 9],
            vec![0, 1, 0, 0, 3, 0, 1, 3, 4],
            vec![1.0, 2.0, 4.0, 8.0, 16.0, 128.0, 256.0, 512.0, 1024.0],
        );

        assert_eq!(*mat.col_dense(0).values(), vec![1.0, 2.0, 0.0, 0.0, 0.0]);
        assert_eq!(*mat.col_dense(1).values(), vec![4.0, 0.0, 0.0, 0.0, 0.0]);
        assert_eq!(*mat.col_dense(2).values(), vec![8.0, 0.0, 0.0, 16.0, 0.0]);
        assert_eq!(*mat.col_dense(3).values(), vec![0.0, 0.0, 0.0, 0.0, 0.0]);
        assert_eq![*mat.col_dense(4).values(), vec![128.0, 256.0, 0.0, 512.0, 1024.0]];
    }

    #[test]
    fn test_col_dot_000() {
        let mat: CSCMatrix = CSCMatrix::new(
            5,
            5,
            vec![0, 2, 2, 4, 6, 10],
            vec![0, 1, 0, 3, 0, 2, 0, 1, 3, 4],
            vec![1.0, 2.0, 8.0, 16.0, 32.0, 64.0, 128.0, 256.0, 512.0, 1024.0],
        );

        let v1: DenseVector = DenseVector::new(vec![1.0, -1.0, 1.0, -1.0, 1.0]);
        assert_eq!(mat.col_dot(0, &v1), -1.0);
        assert_eq!(mat.col_dot(1, &v1), 0.0);
        assert_eq!(mat.col_dot(2, &v1), -8.0);
        assert_eq!(mat.col_dot(3, &v1), 96.0);
        assert_eq!(mat.col_dot(4, &v1), 384.0);
    }

    #[test]
    fn test_col_submatrix_000() {
        let mat: CSCMatrix = CSCMatrix::new(
            5,
            6,
            vec![0, 2, 3, 5, 5, 7, 11],
            vec![0, 1, 0, 0, 3, 0, 2, 0, 1, 3, 4],
            vec![2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0, 22.0],
        );

        let sm1 = mat.col_submatrix(&vec![2, 5, 0]);
        assert_eq!(sm1.num_rows, 5);
        assert_eq!(sm1.num_cols, 3);
        assert_eq!(sm1.col_boundaries, vec![0, 2, 6, 8]);
        assert_eq!(sm1.row_indices, vec![0, 3, 0, 1, 3, 4, 0, 1]);
        assert_eq!(sm1.values, vec![8.0, 10.0, 16.0, 18.0, 20.0, 22.0, 2.0, 4.0]);

        let sm2 = mat.col_submatrix(&vec![4, 1]);
        assert_eq!(sm2.num_rows, 5);
        assert_eq!(sm2.num_cols, 2);
        assert_eq!(sm2.col_boundaries, vec![0, 2, 3]);
        assert_eq!(sm2.row_indices, vec![0, 2, 0]);
        assert_eq!(sm2.values, vec![12.0, 14.0, 6.0]);
    }

    #[test]
    fn test_new_000() {
        let m1_num_rows = 5;
        let m1_num_cols = 5;
        let m1_col_boundaries = vec![0, 2, 2, 4, 6, 10];
        let m1_row_indices = vec![0, 1, 0, 3, 0, 2, 0, 1, 3, 4];
        let m1_values = vec![1.0, 2.0, 8.0, 16.0, 32.0, 64.0, 128.0, 256.0, 512.0, 1024.0];
        let m1: CSCMatrix = CSCMatrix::new(
            m1_num_rows,
            m1_num_cols,
            m1_col_boundaries.clone(),
            m1_row_indices.clone(),
            m1_values.clone(),
        );
        assert_eq!(m1.num_rows, m1_num_rows);
        assert_eq!(m1.num_cols, m1_num_cols);
        assert_eq!(m1.col_boundaries, m1_col_boundaries);
        assert_eq!(m1.row_indices, m1_row_indices);
        assert_eq!(m1.values, m1_values);

        let m2_num_rows = 5;
        let m2_num_cols = 1;
        let m2_col_boundaries = vec![0, 2];
        let m2_row_indices = vec![1, 4];
        let m2_values = vec![4.0, 9.0];
        let m2: CSCMatrix = CSCMatrix::new(
            m2_num_rows,
            m2_num_cols,
            m2_col_boundaries.clone(),
            m2_row_indices.clone(),
            m2_values.clone(),
        );
        assert_eq!(m2.num_rows, m2_num_rows);
        assert_eq!(m2.num_cols, m2_num_cols);
        assert_eq!(m2.col_boundaries, m2_col_boundaries);
        assert_eq!(m2.row_indices, m2_row_indices);
        assert_eq!(m2.values, m2_values);

        let m3_num_rows = 1;
        let m3_num_cols = 5;
        let m3_col_boundaries = vec![0, 0, 1, 1, 2, 2];
        let m3_row_indices = vec![0, 0];
        let m3_values = vec![1.0, 2.0];
        let m3: CSCMatrix = CSCMatrix::new(
            m3_num_rows,
            m3_num_cols,
            m3_col_boundaries.clone(),
            m3_row_indices.clone(),
            m3_values.clone(),
        );
        assert_eq!(m3.num_rows, m3_num_rows);
        assert_eq!(m3.num_cols, m3_num_cols);
        assert_eq!(m3.col_boundaries, m3_col_boundaries);
        assert_eq!(m3.row_indices, m3_row_indices);
        assert_eq!(m3.values, m3_values);
    }

    #[test]
    fn new_from_columns_000() {
        let m1: CSCMatrix = CSCMatrix::new_from_columns(vec![
            SparseVector::new(5, vec![0, 1], vec![1.0, 2.0]),
            SparseVector::new(5, vec![], vec![]),
            SparseVector::new(5, vec![0, 3], vec![8.0, 16.0]),
            SparseVector::new(5, vec![0, 2], vec![32.0, 64.0]),
            SparseVector::new(5, vec![0, 1, 3, 4], vec![128.0, 256.0, 512.0, 1024.0]),
        ]);
        assert_eq!(m1.num_rows, 5);
        assert_eq!(m1.num_cols, 5);
        assert_eq!(m1.col_boundaries, vec![0, 2, 2, 4, 6, 10]);
        assert_eq!(m1.row_indices, vec![0, 1, 0, 3, 0, 2, 0, 1, 3, 4]);
        assert_eq!(m1.values, vec![1.0, 2.0, 8.0, 16.0, 32.0, 64.0, 128.0, 256.0, 512.0, 1024.0]);
        // todo!();
    }

    #[test]
    fn test_pop_rows_inplace_000() {
        let mat: CSCMatrix = CSCMatrix::new(
            6,
            5,
            vec![0, 2, 4, 6, 8, 11],
            vec![4, 5, 0, 5, 1, 5, 0, 5, 2, 3, 5],
            vec![6.0, 9.0, 1.0, 9.0, 3.0, 9.0, 2.0, 9.0, 4.0, 5.0, 9.0],
        );

        let mut m0 = mat.clone();
        m0.pop_rows_inplace(0);
        assert_eq!(m0.num_rows, mat.num_rows);
        assert_eq!(m0.num_cols, mat.num_cols);
        assert_eq!(m0.col_boundaries, mat.col_boundaries);
        assert_eq!(m0.row_indices, mat.row_indices);
        assert_eq!(m0.values, mat.values);

        let mut m1 = mat.clone();
        m1.pop_rows_inplace(1);
        assert_eq!(m1.num_rows, 5);
        assert_eq!(m1.num_cols, 5);
        assert_eq!(m1.col_boundaries, vec![0, 1, 2, 3, 4, 6]);
        assert_eq!(m1.row_indices, vec![4, 0, 1, 0, 2, 3]);
        assert_eq!(m1.values, vec![6.0, 1.0, 3.0, 2.0, 4.0, 5.0]);

        let mut m2 = mat.clone();
        m2.pop_rows_inplace(2);
        assert_eq!(m2.num_rows, 4);
        assert_eq!(m2.num_cols, 5);
        assert_eq!(m2.col_boundaries, vec![0, 0, 1, 2, 3, 5]);
        assert_eq!(m2.row_indices, vec![0, 1, 0, 2, 3]);
        assert_eq!(m2.values, vec![1.0, 3.0, 2.0, 4.0, 5.0]);

        let mut m3 = mat.clone();
        m3.pop_rows_inplace(3);
        assert_eq!(m3.num_rows, 3);
        assert_eq!(m3.num_cols, 5);
        assert_eq!(m3.col_boundaries, vec![0, 0, 1, 2, 3, 4]);
        assert_eq!(m3.row_indices, vec![0, 1, 0, 2]);
        assert_eq!(m3.values, vec![1.0, 3.0, 2.0, 4.0]);

        let mut m4 = mat.clone();
        m4.pop_rows_inplace(4);
        assert_eq!(m4.num_rows, 2);
        assert_eq!(m4.num_cols, 5);
        assert_eq!(m4.col_boundaries, vec![0, 0, 1, 2, 3, 3]);
        assert_eq!(m4.row_indices, vec![0, 1, 0]);
        assert_eq!(m4.values, vec![1.0, 3.0, 2.0]);

        let mut m5 = mat.clone();
        m5.pop_rows_inplace(5);
        assert_eq!(m5.num_rows, 1);
        assert_eq!(m5.num_cols, 5);
        assert_eq!(m5.col_boundaries, vec![0, 0, 1, 1, 2, 2]);
        assert_eq!(m5.values, vec![1.0, 2.0]);
    }

    #[test]
    fn test_remove_cols_inplace_000() {
        let mat: CSCMatrix = CSCMatrix::new(
            5,
            6,
            vec![0, 2, 3, 3, 5, 7, 11],
            vec![0, 1, 0, 0, 3, 0, 2, 0, 1, 3, 4],
            vec![2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0, 22.0],
        );

        let mut mat1 = mat.clone();
        mat1.remove_cols_inplace(&vec![0]);
        assert_eq!(mat1.num_rows, 5);
        assert_eq!(mat1.num_cols, 5);
        assert_eq!(mat1.col_boundaries, vec![0, 1, 1, 3, 5, 9]);
        assert_eq!(mat1.row_indices, vec![0, 0, 3, 0, 2, 0, 1, 3, 4]);
        assert_eq!(mat1.values, vec![6.0, 8.0, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0, 22.0]);

        let mut mat2 = mat.clone();
        mat2.remove_cols_inplace(&vec![5]);
        assert_eq!(mat2.num_rows, 5);
        assert_eq!(mat2.num_cols, 5);
        assert_eq!(mat2.col_boundaries, vec![0, 2, 3, 3, 5, 7]);
        assert_eq!(mat2.row_indices, vec![0, 1, 0, 0, 3, 0, 2]);
        assert_eq!(mat2.values, vec![2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0]);

        let mut mat3 = mat.clone();
        mat3.remove_cols_inplace(&vec![1, 4]);
        assert_eq!(mat3.num_rows, 5);
        assert_eq!(mat3.num_cols, 4);
        assert_eq!(mat3.col_boundaries, vec![0, 2, 2, 4, 8]);
        assert_eq!(mat3.row_indices, vec![0, 1, 0, 3, 0, 1, 3, 4]);
        assert_eq!(mat3.values, vec![2.0, 4.0, 8.0, 10.0, 16.0, 18.0, 20.0, 22.0]);

        let mut mat4 = mat.clone();
        mat4.remove_cols_inplace::<i32>(&vec![]);
        assert_eq!(mat4.num_rows, mat.num_rows);
        assert_eq!(mat4.num_cols, mat.num_cols);
        assert_eq!(mat4.col_boundaries, mat.col_boundaries);
        assert_eq!(mat4.row_indices, mat.row_indices);
        assert_eq!(mat4.values, mat.values);

        // TODO: unit tests for panic conditions
    }

}
