use crate::tree;

pub enum TreePriorParams { 
    Uniform { root: f64 },
}


pub struct TreeParams {
    pub prior: TreePriorParams,
    pub tree: tree::Tree,
}

pub enum SubstitutionModelParams {
    BinaryGTR { pi_one: f64 },
}

pub struct ASRVParams {
    pub shape: f64,
    pub rates: Vec<f64>,
}

pub struct ABRVParams {
    pub shape: f64,
}

pub struct TraitsParams {
    pub num_traits: i32,
    pub subst: SubstitutionModelParams,
    pub base: f64,
    pub asrv: ASRVParams,
    pub abrv: ABRVParams,
}


pub struct Parameters {
    pub tree: TreeParams,
    pub traits: TraitsParams,
}


