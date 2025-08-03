// TODO: store BOUNDS name, error on more than one
// TODO: store RANGES name, error on more than one
// TODO: store RHS name, error on more than one



use crate::constr::milp::MILPConstrSense;
use crate::constr::milp::MILPModel;
use crate::primitives::int::IntNonnegFitsInUsize;
use crate::primitives::matrix::CSCMatrix;
use crate::primitives::vector::DenseVector;
use crate::primitives::vector::SparseVector;

use std::collections::HashMap;
use std::fs::File;
use std::io::BufRead;  // for BufReader.lines()
use std::path::Path;



pub fn read_mps_from_path<P: AsRef<Path>>(file_path: P) -> std::io::Result<MILPModel> {
    let file = File::open(file_path)?;
    read_mps(file)
}



const RECOGNIZED_HEADERS: [&str; 5] = ["ROWS", "COLUMNS", "RHS", "RANGES", "BOUNDS"];  // besides NAME and ENDATA

struct MPSParsingState where {
    current_section: String,
    current_is_int: bool,
    endata: bool,
    data: MPSData,
}

struct MPSData where {
    model_name: String,
    mps_rows: Vec<MPSDataRow>,
    mps_row_idx_by_name: HashMap<String,usize>,
    obj_row_name: Option<String>,
    vars: Vec<MPSDataVar>,
    var_idx_by_name: HashMap<String,usize>,
}

struct MPSDataBound {
    bound_type: MPSDataBoundType,
    value: Option<f64>,
}

#[derive(Clone, Copy)]
enum MPSDataBoundType {
    FX,  // fixed value
    UP,  // upper bound
    LO,  // lower bound
    PL,  // "plus" (upper bound of +inf)
    MI,  // "minus" (lower bound of -inf)
    FR,  // "free" (-inf to +inf)
    BV,  // binary variable
    LI,  // lower bound for int variable (not in original format)
    UI,  // upper bound for int variable (not in original format)
}

impl MPSDataBoundType {
    fn from_str(string: &str) -> MPSDataBoundType {
        match string {
            "FX" => MPSDataBoundType::FX,
            "UP" => MPSDataBoundType::UP,
            "LO" => MPSDataBoundType::LO,
            "PL" => MPSDataBoundType::PL,
            "MI" => MPSDataBoundType::MI,
            "FR" => MPSDataBoundType::FR,
            "BV" => MPSDataBoundType::BV,
            "LI" => MPSDataBoundType::LI,
            "UI" => MPSDataBoundType::UI,
            other => panic!("Unrecognized bound type: {}", other)
        }
    }
}

struct MPSDataRow {
    name: String,
    row_type: MPSDataRowType,
    rhs: Option<f64>,
    range: Option<f64>,
    var_coef_indices: Vec<usize>,
    var_coef_values: Vec<f64>,
}

#[derive(Clone, Copy, PartialEq)]
enum MPSDataRowType {
    N,  // non-constraint row (first one is objective)
    E,  // equal to right-hand side
    G,  // greater than or equal to right-hand side
    L,  // less than or equal to right-hand side
}

impl MPSDataRowType {
    fn from_str(string: &str) -> MPSDataRowType {
        match string {
            "N" => MPSDataRowType::N,
            "E" => MPSDataRowType::E,
            "G" => MPSDataRowType::G,
            "L" => MPSDataRowType::L,
            other => panic!("Unrecognized row type: {}", other)
        }
    }
}

struct MPSDataVar {
    name: String,
    is_int: bool,
    bound_records: Vec<MPSDataBound>,
}



// TODO: error on duplicate/conflicting BOUNDS records
// TODO FUTURE: replace panics with more graceful error handling (new Result type that wraps either an io error or parsing error)
fn read_mps(file: File) -> std::io::Result<MILPModel> {
    // initialize parsing state
    let mut state = MPSParsingState {
        current_section: "".to_string(),
        current_is_int: false,
        endata: false,
        data: MPSData {
            model_name: "".to_string(),
            mps_rows: vec![],
            mps_row_idx_by_name: HashMap::new(),
            obj_row_name: None,
            vars: vec![],
            var_idx_by_name: HashMap::new(),
        }
    };

    // loop over each line in the file (standard ordering of sections enables a
    // single pass; no tolerance here for out-of-order sections)
    for line in std::io::BufReader::new(file).lines() {
        let line = line?;

        // skip comment lines and blank lines
        if line.starts_with("*") || line.trim().len() == 0 {
            continue;
        }

        // lines that *don't* start with a space are headers
        if !line.starts_with(" ") {
            read_mps_process_header_line(&mut state, &line);
            if state.endata {
                break;
            }
        }
        else {
            let line_trimmed = line.trim();
            match state.current_section.as_str() {
                "ROWS" => read_mps_process_section_rows(&mut state, &line_trimmed),
                "COLUMNS" => read_mps_process_section_columns(&mut state, &line_trimmed),
                "RHS" => read_mps_process_section_rhs(&mut state, &line_trimmed),
                "RANGES" => read_mps_process_section_ranges(&mut state, &line_trimmed),
                "BOUNDS" => read_mps_process_section_bounds(&mut state, &line_trimmed),
                _ => panic!("Unrecognized section: {}", state.current_section),
            }
        }
    }
    
    if !state.endata {
        panic!("End of MPS file not reached! A read error or corrupted file");
    }

    // done reading data; construct and return the model
    Ok(read_mps_construct_model(state.data))
}



/// `read_mps` child routine:
/// Constructs the `MILPModel` instance from the data read from the MPS file.
fn read_mps_construct_model(mut data: MPSData) -> MILPModel {
    // TODO FUTURE: reorder variables so integers are contiguous?

    let num_rows_with_range = data.mps_rows.iter()
        .filter(|r| r.row_type != MPSDataRowType::N)
        .filter(|r| r.range.is_some_and(|x| x != 0.0)).count();
    let num_rows_n_type = data.mps_rows.iter()
        .filter(|r| r.row_type == MPSDataRowType::N).count();
    let num_mps_rows_total = data.mps_rows.len() + num_rows_with_range;
    let num_constr_rows = num_mps_rows_total - num_rows_n_type;

    let num_cols = data.vars.len();
    
    // construct objective vector
    let obj_coefs: DenseVector<f64>;
    let obj_offset: f64;
    match data.obj_row_name {
        None => {
            // no objective row; make all coefficients zero
            obj_coefs = DenseVector::new_zeros(num_cols);
            obj_offset = 0.0;
        },
        Some(obj_row_name) => {
            // put objective row into DenseVector (destructively)
            let obj_row_idx = *data.mps_row_idx_by_name.get(&obj_row_name).unwrap();
            let obj_row = &mut data.mps_rows[obj_row_idx.to_usize_unchecked()];
            let obj_coef_indices = std::mem::take(&mut obj_row.var_coef_indices);
            let obj_coef_values = std::mem::take(&mut obj_row.var_coef_values);

            obj_coefs = DenseVector::new_from_sparse_lists(
                num_cols,
                obj_coef_indices,
                obj_coef_values,
            );
            obj_offset = -obj_row.rhs.unwrap_or(0.0);
        }
    }

    // construct coefficient matrix and rhs vector
    let mut columns: Vec<SparseVector<f64>> = (0..(data.vars.len()))
        .map(|_| SparseVector::new_zeros(num_constr_rows)).collect();
    let mut constr_rhs_indices: Vec<usize> = vec![];
    let mut constr_rhs_values: Vec<f64> = vec![];
    let mut constr_sense: Vec<MILPConstrSense> = vec![];
    let mut constr_idx_next: usize = 0;
    let mut add_constr_info = |idx, sense, rhs| {
        constr_sense.push(sense);
        if rhs != 0.0 {
            constr_rhs_indices.push(idx);
            constr_rhs_values.push(rhs);
        }
    };
    for mut mps_row in data.mps_rows {
        // skip non-constraint rows (already processed objective row above)
        if let MPSDataRowType::N = mps_row.row_type {
            continue;
        }

        // add two "constraint rows" if this row has a nonzero "range"
        let range = mps_row.range.unwrap_or(0.0);
        let constr_idx = constr_idx_next;
        constr_idx_next += 1;
        if range != 0.0 {
            constr_idx_next += 1;
        }

        // store constraint sense and right-hand side
        // (also for a second constraint row if indicated by RANGES)
        let rhs = mps_row.rhs.unwrap_or(0.0);
        match mps_row.row_type {
            MPSDataRowType::E => {
                if range == 0.0 {
                    add_constr_info(constr_idx, MILPConstrSense::Equal, rhs);
                }
                else if range > 0.0 {
                    add_constr_info(constr_idx, MILPConstrSense::Greater, rhs);
                    add_constr_info(constr_idx, MILPConstrSense::Less, rhs + range);
                }
                else {
                    add_constr_info(constr_idx, MILPConstrSense::Less, rhs);
                    add_constr_info(constr_idx, MILPConstrSense::Greater, rhs + range);
                }
            },
            MPSDataRowType::L => {
                add_constr_info(constr_idx, MILPConstrSense::Less, rhs);
                if range != 0.0 {
                    add_constr_info(constr_idx + 1, MILPConstrSense::Greater, rhs - range.abs());
                }
            },
            MPSDataRowType::G => {
                add_constr_info(constr_idx, MILPConstrSense::Greater, rhs);
                if range != 0.0 {
                    add_constr_info(constr_idx + 1, MILPConstrSense::Less, rhs + range.abs());
                }
            },
            MPSDataRowType::N => unreachable!(),
        }

        // move row elements to the appropriate columns (destructively)
        let var_coef_indices = std::mem::take(&mut mps_row.var_coef_indices);
        let var_coef_values = std::mem::take(&mut mps_row.var_coef_values);
        assert_eq!(var_coef_indices.len(), var_coef_values.len());
        for i in 0..var_coef_indices.len() {
            let c = var_coef_indices[i].to_usize_unchecked();
            let value = var_coef_values[i];
            columns[c].set_value_in_order(constr_idx, value);

            // add duplicate row if a nonzero "range" is specified
            if range != 0.0 {
                columns[c].set_value_in_order(constr_idx + 1, value);
            }
        }
    }
    assert_eq!(constr_idx_next, num_constr_rows);  // guardrail for future changes

    let constr_rhs: DenseVector<f64> = DenseVector::new_from_sparse_lists(num_constr_rows, constr_rhs_indices, constr_rhs_values);
    let coef_matrix: CSCMatrix<f64> = CSCMatrix::new_from_columns(columns);

    // process variable bounds and type
    let mut var_lb: Vec<Option<f64>> = vec![];
    let mut var_ub: Vec<Option<f64>> = vec![];
    let mut var_is_int: Vec<bool> = vec![];
    for var in data.vars {
        var_is_int.push(var.is_int);

        let mut lbs = vec![];
        let mut ubs = vec![];
        let expect_value = |t, v: Option<f64>| v.expect(format!("{} bound without value for var {}", t, var.name).as_str());
        let expect_no_value = |t, v: Option<f64>| assert!(v.is_none(), "{} bound with value for var {}", t, var.name);
        let expect_int_var = |t| assert!(var.is_int, "{} bound for non-int var {}", t, var.name);
        for bound_record in var.bound_records {
            match bound_record.bound_type {
                MPSDataBoundType::FX => {
                    let value = expect_value("FX", bound_record.value);
                    lbs.push(Some(value));
                    ubs.push(Some(value));
                },
                MPSDataBoundType::UP => {
                    let value = expect_value("UP", bound_record.value);
                    ubs.push(Some(value));
                },
                MPSDataBoundType::LO => {
                    let value = expect_value("LO", bound_record.value);
                    lbs.push(Some(value));
                },
                MPSDataBoundType::PL => {
                    expect_no_value("PL", bound_record.value);
                    ubs.push(None);
                },
                MPSDataBoundType::MI => {
                    expect_no_value("MI", bound_record.value);
                    lbs.push(None);
                },
                MPSDataBoundType::FR => {
                    expect_no_value("FR", bound_record.value);
                    lbs.push(None);
                    ubs.push(None);
                },
                MPSDataBoundType::BV => {
                    expect_int_var("BV");
                    expect_no_value("BV", bound_record.value);
                    lbs.push(Some(0.0));
                    ubs.push(Some(1.0));
                },
                MPSDataBoundType::LI => {
                    expect_int_var("LI");
                    let value = expect_value("LI", bound_record.value);
                    lbs.push(Some(value));
                },
                MPSDataBoundType::UI => {
                    expect_int_var("UI");
                    let value = expect_value("UI", bound_record.value);
                    ubs.push(Some(value));
                },
            }
        }

        // get bounds, default to [0, +inf)
        let lb = if lbs.len() > 0 {lbs[0]} else {Some(0.0)};
        let ub = if ubs.len() > 0 {ubs[0]} else {None};

        // check for conflicting bound records
        let lbdescr = |v: Option<f64>| match v { None => "-inf".to_string(), Some(x) => x.to_string() };
        if lbs.len() > 1 {
            for i in 1..lbs.len() {
                assert_eq!(lb, lbs[i], "Conflicting lower bounds for var {}: {}, {}", var.name, lbdescr(lb), lbdescr(lbs[i]));
            }
        }
        let ubdescr = |v: Option<f64>| match v { None => "+inf".to_string(), Some(x) => x.to_string() };
        if ubs.len() > 1 {
            for i in 1..ubs.len() {
                assert_eq!(ub, ubs[i], "Conflicting upper bounds for var {}: {}, {}", var.name, ubdescr(ub), ubdescr(ubs[i]));
            }
        }
        if lb.is_some() && ub.is_some() && lb.unwrap() > ub.unwrap() {
            panic!("Var {} has lower bound greater than upper bound: {}, {}", var.name, lb.unwrap(), ub.unwrap());
        }

        var_lb.push(lb);
        var_ub.push(ub);
    }

    MILPModel::new(
        data.model_name,
        obj_coefs,
        obj_offset,
        coef_matrix,
        constr_rhs,
        constr_sense,
        var_lb,
        var_ub,
        var_is_int,
    )
}



/// `read_mps` child routine:
/// Process a header line.
fn read_mps_process_header_line(state: &mut MPSParsingState, line: &str) {
    let header_line_trimmed = line.trim();

    if header_line_trimmed.starts_with("NAME ") {
        // TODO: are spaces allowed in the name?
        let parts: Vec<&str> = header_line_trimmed.split(" ").filter(|s| !s.is_empty()).collect();
        if parts.len() != 2 {
            panic!("Unexpected extra fields on NAME line after name ({} tokens)", parts.len());
        }
        state.data.model_name = parts[1].to_string();
    }
    else if header_line_trimmed == "NAME" {
        // no name given, keep as blank
    }
    else if header_line_trimmed == "ENDATA" {
        // stop processing; this is specified as the "last" line in the file but there are
        // some files in the wild that contain nonstandard content after the ENDATA line
        state.endata = true;
    }
    else if RECOGNIZED_HEADERS.contains(&header_line_trimmed) {
        state.current_section = header_line_trimmed.to_string();
    }
    else {
        panic!("Unrecognized header line: {}", header_line_trimmed);
    }
}



/// `read_mps` child routine:
/// Process a line in the "BOUNDS" section.
fn read_mps_process_section_bounds(state: &mut MPSParsingState, line_trimmed: &str) {
    let fields: Vec<&str> = line_trimmed.split(" ").filter(|s| !s.is_empty()).collect();
    let fields_len = fields.len();
    if fields_len != 3 && fields_len != 4 {
        panic!("Lines in BOUNDS section should have 3 or 4 fields, found one with {}: {}", fields_len, line_trimmed);
    }

    // TODO: multiple BOUNDS names? currently ignoring name in `fields[1]`
    let bound_type = MPSDataBoundType::from_str(fields[0]);
    let var_name = fields[2];

    // read bound value if provided
    // (require it to be a valid f64 if present)
    let bound_value: Option<f64>;
    if fields_len == 4 {
        let parsed_value = fields[3].parse()
            .expect(format!("Bound value for variable \"{}\" is not a number: {}", var_name, fields[3]).as_str());
        bound_value = Some(parsed_value);
    }
    else {
        bound_value = None;
    }

    // look up variable
    let var_idx = *state.data.var_idx_by_name.get(var_name)
        .expect(format!("Bound in BOUNDS section for variable not in COLUMNS section: {}", var_name).as_str());

    // store bound record
    let bound_data = MPSDataBound {
        bound_type: bound_type,
        value: bound_value,
    };
    state.data.vars[var_idx.to_usize_unchecked()].bound_records.push(bound_data);
}



/// `read_mps` child routine:
/// Process a line in the "COLUMNS" section.
fn read_mps_process_section_columns(state: &mut MPSParsingState, line_trimmed: &str) {
    let fields: Vec<&str> = line_trimmed.split(" ").filter(|s| !s.is_empty()).collect();
    if fields.len() == 3 {
        if fields[1] == "'MARKER'" {
            // this is a marker line
            // note: `fields[0]` unused because leftmost field is marker name which has no significance
            if fields[2] == "'INTORG'" {
                if state.current_is_int {
                    panic!("Malformed INTORG/INTEND markers: found second INTORG before an INTEND");
                }
                state.current_is_int = true;
            }
            else if fields[2] == "'INTEND'" {
                if !state.current_is_int {
                    panic!("Malformed INTORG/INTEND markers: found an INTEND without a prior INTORG");
                }
                state.current_is_int = false;
            }
            else {
                panic!("Unrecognized marker type in COLUMNS section: {}", fields[2]);
            }
        }
        else {
            // this is a coefficient line, with one coefficient
            let var_name = fields[0];
            let row_name = fields[1];
            let coef_text = fields[2];
            read_mps_process_section_columns_coef(state, var_name, row_name, coef_text);
        }
    }
    else if fields.len() == 5 {
        // this is a coefficient line, with two coefficients
        // (no special semantics except that they're for the same variable;
        // historically doubled up because this was a punch-card format and
        // this could save a lot of cards)
        let var_name = fields[0];
        let row1_name = fields[1];
        let coef1_text = fields[2];
        read_mps_process_section_columns_coef(state, var_name, row1_name, coef1_text);
        let row2_name = fields[3];
        let coef2_text = fields[4];
        read_mps_process_section_columns_coef(state, var_name, row2_name, coef2_text);
    }
    else {
        panic!("Lines in COLUMNS section should have 3 or 5 fields, found one with {}: {}", fields.len(), line_trimmed);
    }
}



/// `read_mps` child routine:
/// Process a single coefficient listed in the "COLUMNS" section.
fn read_mps_process_section_columns_coef(state: &mut MPSParsingState, var_name: &str, row_name: &str, coef_text: &str) {
    // get MPS row index
    let row_idx = *state.data.mps_row_idx_by_name.get(row_name)
        .expect(format!("Coefficient in COLUMNS section for row not in ROWS section: {}", row_name).as_str());

    // get var index, add variable if it hasn't been already
    let var_idx: usize;
    match state.data.var_idx_by_name.get(var_name) {
        Some(idx) => var_idx = *idx,
        None => {
            var_idx = state.data.vars.len();

            let var_data = MPSDataVar {
                name: var_name.to_string(),
                is_int: state.current_is_int,
                bound_records: vec![],
            };
            state.data.vars.push(var_data);
            state.data.var_idx_by_name.insert(var_name.to_string(), var_idx);
        }
    }
    
    // parse and store coefficient
    let coef: f64 = coef_text.parse()
        .expect(format!("Could not parse coefficient: {}", coef_text).as_str());
    state.data.mps_rows[row_idx.to_usize_unchecked()].var_coef_indices.push(var_idx);
    state.data.mps_rows[row_idx.to_usize_unchecked()].var_coef_values.push(coef);
}



/// `read_mps` child routine:
/// Process a line in the "RANGES" section.
fn read_mps_process_section_ranges(state: &mut MPSParsingState, line_trimmed: &str) {
    // TODO: multiple RANGE vectors? is there a consensus? (currently ignoring range name in leftmost field)
    let fields: Vec<&str> = line_trimmed.split(" ").filter(|s| !s.is_empty()).collect();
    if fields.len() == 3 {
        let row_name = fields[1];
        let range_value_text = fields[2];
        read_mps_process_section_ranges_value(state, row_name, range_value_text);
    }
    else if fields.len() == 5 {
        let row1_name = fields[1];
        let range1_value_text = fields[2];
        read_mps_process_section_ranges_value(state, row1_name, range1_value_text);
        let row2_name = fields[3];
        let range2_value_text = fields[4];
        read_mps_process_section_ranges_value(state, row2_name, range2_value_text);
    }
    else {
        panic!("Lines in RANGES section should have 3 or 5 fields, found one wih {}: {}", fields.len(), line_trimmed);
    }
}



/// `read_mps` child routine:
/// Process a single range value in the "RANGES" section.
fn read_mps_process_section_ranges_value(state: &mut MPSParsingState, row_name: &str, range_value_text: &str) {
    // get MPS row
    let row_idx = *state.data.mps_row_idx_by_name.get(row_name)
        .expect(format!("Value in RANGES section for row not in ROWS section: {}", row_name).as_str());
    let row_idx_usize = row_idx.to_usize_unchecked();
    let row = &mut state.data.mps_rows[row_idx_usize];

    // validation: do not allow multiple RANGES records for the same row
    if let Some(_) = row.range {
        panic!("Row has multiple values listed in RANGES section: {}", row_name);
    }

    // parse and store range value
    let range_value: f64 = range_value_text.parse()
        .expect(format!("Could not parse range value: {}", range_value_text).as_str());
    row.range = Some(range_value);
}



/// `read_mps` child routine:
/// Process a line in the "RHS" section.
fn read_mps_process_section_rhs(state: &mut MPSParsingState, line_trimmed: &str) {
    
    // TODO: multiple RHS vectors? is there a consensus? (currently ignoring RHS name in leftmost field)
    let fields: Vec<&str> = line_trimmed.split(" ").filter(|s| !s.is_empty()).collect();
    if fields.len() == 3 {
        let row_name = fields[1];
        let rhs_value_text = fields[2];
        read_mps_process_section_rhs_value(state, row_name, rhs_value_text);
    }
    else if fields.len() == 5 {
        let row1_name = fields[1];
        let rhs1_value_text = fields[2];
        read_mps_process_section_rhs_value(state, row1_name, rhs1_value_text);
        let row2_name = fields[3];
        let rhs2_value_text = fields[4];
        read_mps_process_section_rhs_value(state, row2_name, rhs2_value_text);
    }
    else {
        panic!("Lines in RHS section should have 3 or 5 fields, found one with {}: {}", fields.len(), line_trimmed);
    }
}



/// `read_mps` child routine:
/// Process a single right-hand side value in the "RHS" section.
fn read_mps_process_section_rhs_value(state: &mut MPSParsingState, row_name: &str, rhs_value_text: &str) {
    // get MPS row
    let row_idx = *state.data.mps_row_idx_by_name.get(row_name)
        .expect(format!("Value in RHS section for row not in ROWS section: {}", row_name).as_str());
    let row = &mut state.data.mps_rows[row_idx.to_usize_unchecked()];

    // validation: do not allow multiple RHS records for the same row
    if let Some(_) = row.rhs {
        panic!("Row has multiple values listed in RHS section: {}", row_name);
    }

    // parse and store right-hand side value
    let rhs_value: f64 = rhs_value_text.parse()
        .expect(format!("Could not parse right-hand side value: {}", rhs_value_text).as_str());
    row.rhs = Some(rhs_value);
}



/// `read_mps` child routine:
/// Process a line in the "ROWS" section.
fn read_mps_process_section_rows(state: &mut MPSParsingState, line_trimmed: &str) {
    // parse fields
    let fields: Vec<&str> = line_trimmed.split(" ").filter(|s| !s.is_empty()).collect();
    if fields.len() != 2 {
        panic!("Lines in ROWS section should have 2 fields, found one with {}: {}", fields.len(), line_trimmed);
    }
    let row_type = MPSDataRowType::from_str(fields[0]);
    let row_name = fields[1].to_string();

    // validation: no duplicate rows
    if state.data.mps_row_idx_by_name.contains_key(&row_name) {
        panic!("Duplicate lines in ROWS section for row: {}", row_name);
    }

    let row_idx = state.data.mps_rows.len();

    // first "N" type row is the objective
    // (seems to be a consensus to ignore other "N" rows as the standard behavior)
    if let MPSDataRowType::N = row_type {
        if let None = state.data.obj_row_name {
            state.data.obj_row_name = Some(row_name.clone());
        }
    }

    // store row data
    let row_data = MPSDataRow {
        name: row_name.clone(),
        row_type: row_type,
        range: None,
        rhs: None,
        var_coef_indices: vec![],
        var_coef_values: vec![],
    };
    state.data.mps_rows.push(row_data);
    state.data.mps_row_idx_by_name.insert(row_name, row_idx);
}



#[cfg(test)]
mod tests {

    use super::*;
    
    use std::path::PathBuf;

    #[test]
    fn test_read_mps_from_path_000() {
        let mps_path: PathBuf = [env!("CARGO_MANIFEST_DIR"), "models", "miplib2017", "gen-ip054.mps"].iter().collect();
        // let mps_path: PathBuf = [env!("CARGO_MANIFEST_DIR"), "models", "miplib2017", "gen-ip002.mps"].iter().collect();
        // let mps_path: PathBuf = [env!("CARGO_MANIFEST_DIR"), "models", "miplib2017", "mad.mps"].iter().collect();

        let model = read_mps_from_path(mps_path).unwrap();
        crate::constr::milp::solve_milp_bnb(&model);
    }

}
