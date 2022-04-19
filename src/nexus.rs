// =-=-=-=-= nexus.rs =-=-=-=-=
// A NEXUS file is a standardised phylogenetic data format
// We parse a subset of NEXUS to allow us to load phylolingistic data

use std::collections::HashMap;

#[derive(Debug)]
pub struct DataBlock<'src> {
    pub n_taxa: usize,
    pub n_chars: usize,
    pub datatype: &'src str,
    pub gap: &'src str,
    pub missing: &'src str,
    pub matrix: HashMap<&'src str, Vec<char>>,
}

impl<'src> DataBlock<'src> {

    pub fn num_traits(&self) -> Option<i32> {
        match self.matrix.values().next() {
            Some(chars) => { Some(chars.len() as i32) },
            None => None,
        }
    }

    pub fn binary_into_partials(&self) -> Option<Vec<(String, Vec<f64>)>> {
        let mut tips = vec![];
        for (tag, chars) in &self.matrix {
            let mut partials = vec![];
            for c in chars {
                match c {
                    '0' => {
                        partials.extend_from_slice(&[1.0, 0.0]);
                    },
                    '1' => {
                        partials.extend_from_slice(&[0.0, 1.0]);
                    },
                    '?' => {
                        partials.extend_from_slice(&[1.0, 1.0]);
                    },
                    _ => {
                        return None;
                    },
                }
            }
            tips.push((tag.to_string(), partials));
        }
        Some(tips)
    }
}

enum Block<'src> {
    Unrecognised,
    Data(DataBlock<'src>),
}

#[derive(Debug)]
pub struct Nexus<'src> {
    pub data: Option<DataBlock<'src>>, 
}

#[derive(Debug)]
enum Token<'src> {
    Word(&'src str),
    Semi,
    Newline,
}

fn consume_newlines<'src>(toks: &mut Vec<Token<'src>>) {
    while toks.len() > 0 {
        if let Token::Newline = toks[0] {
            toks.drain(..1);
        }
        else {
            return;
        }
    }
}

fn expect_semi<'src>(toks: &mut Vec<Token<'src>>) -> Result<(), &'static str> {
    if toks.len() == 0 { return Err("No remaining tokens (expected ';')"); }

    match toks[0] {
        Token::Newline => { Err("Expected ';', found newline") },
        Token::Word(_) => { Err("Expected ';', found word") },
        Token::Semi => { toks.drain(..1); Ok(()) }, 
    }
}

fn expect_word<'src>(toks: &mut Vec<Token<'src>>) -> Result<&'src str, &'static str> {
    if toks.len() == 0 { return Err("No remaining tokens (expected word)"); }

    match toks[0] {
        Token::Newline => { Err("Expected word, found newline") },
        Token::Word(s) => { toks.drain(..1); Ok(s) },
        Token::Semi => { Err("Expected word, found ';'") }, 
    }
}

fn parse_data_block<'src>(toks: &mut Vec<Token<'src>>) -> Result<Block<'src>, &'static str> {
    consume_newlines(toks);
    let mut block = DataBlock {n_taxa: 0, n_chars: 0, datatype: "",
                               gap: "", missing: "", matrix: HashMap::new()};
    let mut state = ('-', ' ');
    let mut taxbuf = "";
    let mut rowbuf = vec![];
    while toks.len() > 0 {
        if state.0 == '-' {
            // invariant: newlines consumed
            if let Token::Word(com) = toks[0] {
                toks.drain(..1);
                match com.to_uppercase().as_str() {
                    "DIMENSIONS" => {
                        state = ('c', 'd');
                    },
                    "FORMAT" => {
                        state = ('c', 'f');
                    },
                    "MATRIX" => {
                        state = ('m', '-');
                    },
                    "END" => {
                        expect_semi(toks)?;
                        return Ok(Block::Data(block));
                    },
                    _ => {
                        state = ('c', '?');
                    }
                }
            }
            else {
                return Err("Expected command found semi");
            }
        }
        else if state.0 == 'c' {
            match toks[0] {
                Token::Semi => {
                    toks.drain(..1);
                    consume_newlines(toks);
                    state = ('-', ' ');
                },
                Token::Word(arg) => {
                    match state.1 {
                        'd' => { 
                            toks.drain(..1);
                            let parts = arg.split("=").collect::<Vec<&'src str>>();
                            if parts.len() != 2 {
                                return Err("Invalid argument to DIMENSIONS"); 
                            }
                            match parts[0].to_uppercase().as_str() {
                                "NTAX" => {
                                    match parts[1].parse::<usize>() {
                                         Ok(n) => { block.n_taxa = n; },
                                         Err(_) => { return Err("Couldn't parse NTAX"); },
                                    }
                                },
                                "NCHAR" => {
                                    match parts[1].parse::<usize>() {
                                         Ok(n) => { block.n_chars = n; },
                                         Err(_) => { return Err("Couldn't parse NCHAR"); },
                                    }
                                },
                                _ => {},
                            }
                        },
                        'f' => {
                            toks.drain(..1);
                            let parts = arg.split("=").collect::<Vec<&'src str>>();
                            if parts.len() != 2 {
                                return Err("Invalid argument to FORMAT"); 
                            }
                            match parts[0].to_uppercase().as_str() {
                                "DATATYPE" => {
                                    block.datatype = parts[1]; 
                                },
                                "GAP" => {
                                    block.gap = parts[1];
                                },
                                "MISSING" => {
                                    block.missing = parts[1];
                                },
                                _ => {},
                            }
                        },
                        '?' => {
                            toks.drain(..1);
                        },
                        _ => { unreachable!(); },
                    }
                },
                Token::Newline => {
                    return Err("Newline within single-line command");
                }
            }
        }
        else if state.0 == 'm' {
            if state.1 == '-' {
                match toks[0] {
                    Token::Semi => {
                        toks.drain(..1);
                        consume_newlines(toks);
                        state = ('-', ' ');
                    },
                    Token::Word(t) => {
                        toks.drain(..1);
                        taxbuf = t;
                        state = ('m', 'r');
                    },
                    Token::Newline => {
                        toks.drain(..1);
                    }
                }
            }
            else if state.1 == 'r' {
                match toks[0] {
                    Token::Semi => {
                        toks.drain(..1);

                        block.matrix.insert(taxbuf, rowbuf);
                        taxbuf = "";
                        rowbuf = vec![];

                        consume_newlines(toks);
                        state = ('-', ' ');
                    },
                    Token::Word(t) => {
                        toks.drain(..1);
                        for c in t.chars() {
                            rowbuf.push(c);
                        }
                    },
                    Token::Newline => {
                        toks.drain(..1);

                        block.matrix.insert(taxbuf, rowbuf);
                        taxbuf = "";
                        rowbuf = vec![];

                        state = ('m', '-');
                    }
                }
            }
            
        }

    }
    Err("Did not reach end of block")
}

fn parse_unrecognised_block<'src>(toks: &mut Vec<Token<'src>>) -> Result<Block<'src>, &'static str> {
    let mut exp_com = false;
    loop {
        if toks.len() == 0 { return Err("No remaining tokens (in unrecognised block)"); }

        match toks[0] {
            Token::Newline => {
                toks.drain(..1);
            },
            Token::Semi => {
                toks.drain(..1); 
                exp_com = true;
            }
            Token::Word(s) => {
                if exp_com {
                    if s.to_uppercase().as_str() == "END" {
                        toks.drain(..1);
                        expect_semi(toks)?;
                        break;
                    }
                }
                else { toks.drain(..1); }
            },
        }
    }
    Ok(Block::Unrecognised)
}

fn parse_block<'src>(toks: &mut Vec<Token<'src>>) -> Result<Block<'src>, &'static str> {
    consume_newlines(toks);
    match expect_word(toks)?.to_uppercase().as_str() {
        "BEGIN" => {},
        _ => { return Err("Did not see block BEGIN"); },
    }
    match expect_word(toks)?.to_uppercase().as_str() {
        "DATA" => { expect_semi(toks)?; parse_data_block(toks) },
        _ => { expect_semi(toks)?; parse_unrecognised_block(toks) },
    }
}

fn parse_nexus<'src>(toks: &mut Vec<Token<'src>>) -> Result<Nexus<'src>, &'static str> {
    consume_newlines(toks);
    match expect_word(toks)? {
        "#NEXUS" => {},
        _ => { return Err("Did not see nexus initial token '#NEXUS'"); },
    }

    let mut nexus = Nexus { data: None };
    while toks.len() > 0 { 
        match parse_block(toks)? {
            Block::Data(b) => { nexus.data = Some(b); },
            Block::Unrecognised => { },
        }
        consume_newlines(toks);
    }
    Ok(nexus)
}

pub fn parse<'src>(file: &'src str) -> Result<Nexus<'src>, &'static str> {
    let mut toks: Vec<Token<'src>> = vec![];
    let mut lptr = 0;
    let mut state = '-'; 
    for (i, c) in file.chars().enumerate() {
        if c == ';' {
            if state == 'w' {
                toks.push(Token::Word(&file[lptr..i]));
            }
            toks.push(Token::Semi); 
            state = '-';
        }
        else if c.is_whitespace() {
            if state == 'w' {
                toks.push(Token::Word(&file[lptr..i]));
            }
            if c == '\n' {
                toks.push(Token::Newline);         
            }
            state = '-';
        }
        else {
           if state == '-' {
               lptr = i;
               state = 'w';
           }
        }
    }

    parse_nexus(&mut toks)
}
