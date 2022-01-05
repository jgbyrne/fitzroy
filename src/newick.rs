use crate::tree::*;

#[derive(Debug)]
enum Token {
    LPAREN,
    RPAREN,
    COMMA,
    COLON,
    NUMBER(f64),
    NAME(String),
    SEMI,
}

fn tokenise(s: String) -> Result<Vec<Token>, String> {
    let mut toks = vec![];
    let mut state = '-';
    let mut buf = String::new();
    for c in s.chars() {
        // lookahead states
        match state {
            'n' => {
                if c == '.' || c.is_ascii_digit() {
                    buf.push(c);
                    continue;
                }
                else {
                    match buf.parse::<f64>() {
                        Ok(num) => toks.push(Token::NUMBER(num)),
                        Err(_) => return Err(String::from("Could not parse float")),
                    }
                    buf = String::new();
                    state = '-';
                }
            },
            's' => {
                if c.is_alphabetic() {
                    buf.push(c);
                    continue;
                }
                else {
                    toks.push(Token::NAME(buf));
                    buf = String::new();
                    state = '-';
                }
            },
            _ => {},
        }

        match state {
            '-' => {
                match c {
                    '(' => toks.push(Token::LPAREN),
                    ')' => toks.push(Token::RPAREN),
                    ',' => toks.push(Token::COMMA),
                    ':' => toks.push(Token::COLON),
                    ';' => toks.push(Token::SEMI),
                    '\'' => {
                        state = 'q';
                    },
                    c if c.is_whitespace() => continue,
                    c if c.is_ascii_digit() => {
                        buf.push(c);
                        state = 'n';
                    },
                    c if c.is_alphabetic() => {
                        buf.push(c);
                        state = 's';
                    },
                    _ => {
                        return Err(format!("Invalid character: {}", c));
                    },
                }
            },
            'q' => {
                if c == '\\' {
                    state = 'e';
                    continue;
                }
                if c == '\'' {
                    toks.push(Token::NAME(buf));
                    buf = String::new();
                    state = '-';
                    continue;
                }
                buf.push(c);
            },
            'e' => {
                buf.push(c);
                state = 'q';
                continue;
            },
            _ => unreachable!(),
        }
    }

    if state != '-' {
        return Err(String::from("Unfinished description"));
    }

    Ok(toks)
}

#[derive(Debug)]
enum NodeValue {
    Tree { sub: usize },
    Leaf { name: usize },
    Internal { branches: Vec<usize>, name: usize },
    Branch {sub: usize, length: usize},
    Name(Option<String>),
    Length(Option<f64>),
}

#[derive(Debug)]
struct Node {
    id: usize,
    val: NodeValue,
}

#[derive(Debug)]
struct Nodes(Vec<Node>);

impl Nodes {
    fn add(&mut self, val: NodeValue) -> usize {
        let id = self.0.len();
        self.0.push(Node {id, val});
        id
    }
}

fn parse_branch(toks: &mut Vec<Token>, nodes: &mut Nodes) -> Result<usize, String> {
    let sub = parse_subtree(toks, nodes)?;

    if toks.len() < 2 { // invariant since ) and ;
        return Err(String::from("Cannot parse branch: ran out of tokens"));
    }

    let len = match &toks[0] {
        Token::COLON =>  {
            toks.drain(..1);
            if let Token::NUMBER(num) = toks[0] {
                toks.drain(..1);
                Some(num)
            }
            else {
                return Err(String::from("Cannot parse branch: length not number"));
            }
        },
        _ => None,
    };

    let length = nodes.add(NodeValue::Length(len));
    Ok(nodes.add(NodeValue::Branch{ sub, length }))
}

fn parse_internal(toks: &mut Vec<Token>, nodes: &mut Nodes) -> Result<usize, String> {
    // invariant toks = [(, .....]

    toks.drain(..1); // LPAREN

    // branchset
    let mut branches: Vec<usize> = vec![];
    loop {
        // branch
        let branch = parse_branch(toks, nodes)?;
        branches.push(branch);

        if toks.len() < 1 {
            return Err(String::from("Cannot parse internal: not enough tokens"));
        }

        if let Token::RPAREN = toks[0] {
            toks.drain(..1);
            break;
        }
        else if let Token::COMMA = toks[0] {
            toks.drain(..1);
            continue;
        }
        else {
            return Err(String::from("Cannot parse internal: unexpected token"));
        }
    }

    let mut named = false;
    let name = if let Token::NAME(name) = &toks[0] {
        named = true;
        nodes.add(NodeValue::Name ( Some(name.clone()) ))
    }
    else {
        nodes.add(NodeValue::Name (None))
    };
    if named { toks.drain(..1); }

    let internal = nodes.add(NodeValue::Internal { branches, name });
    Ok(internal)
}


fn parse_subtree(toks: &mut Vec<Token>, nodes: &mut Nodes) -> Result<usize, String> {
    if toks.len() < 1 {
        return Err(String::from("Cannot parse subtree: not enough tokens"));
    }
    if let Token::NAME(ref name) = toks[0] {
        let name = nodes.add(NodeValue::Name(Some(name.clone())));
        toks.drain(..1);
        Ok(nodes.add(NodeValue::Leaf {name}))
    }
    else if let Token::LPAREN = toks[0] {
        let internal = parse_internal(toks, nodes)?;
        Ok(internal)
    }
    else {
        Err(String::from("Cannot parse subtree: no good token"))
    }
}

fn parse_tree(toks: &mut Vec<Token>) -> Result<(Tree, Labels), String> {
    let mut pnodes = Nodes(vec![]); 
    let sub = parse_subtree(toks, &mut pnodes)?;
    let tree = pnodes.add(NodeValue::Tree { sub });

    if toks.len() != 1 {
        return Err(String::from("Cannot parse tree: malformed (missing ';'?)"))
    }

    if !matches!(toks[0], Token::SEMI) {
        return Err(String::from("Cannot parse tree: malformed"));
    }

    let mut nodes = vec![];
    let mut labels: Labels = vec![];

    let pnodes = pnodes.0;

    let LCHILD = 1;
    let RCHILD = 2;

    let mut nstack = vec![(tree, 0, 0)];
    while let Some((node_pid, parent_id, loc)) = nstack.pop() {
        match &pnodes[node_pid].val {
            NodeValue::Tree { sub } => {
                if let NodeValue::Internal {branches, name} = &pnodes[*sub].val {
                    nodes.push(TreeNode { id: 0, parent: 0, lchild: 0, rchild: 0, length: 0.0 });
                    match &pnodes[*name].val {
                        NodeValue::Name(Some(s)) => { labels.push(Label::Interior {name: s.clone()}); },
                        NodeValue::Name(None) => { labels.push(Label::Interior {name: String::new()}); },
                        _ => { panic!("TreeGen: `name` not NodeValue::Name"); },
                    }
                    if branches.len() != 2 { return Err(String::from("TreeGen: Not a binary tree")); }
                    nstack.push((branches[0], 0, LCHILD));
                    nstack.push((branches[1], 0, RCHILD));
                }
                else {
                    panic!("Tree does not contain NodeValue::Internal");
                }
            },
            NodeValue::Branch { sub, length } => {
                let id = nodes.len();
                let length = match &pnodes[*length].val {
                    NodeValue::Length(Some(l)) => *l,
                    NodeValue::Length(None) => 0.0,
                    _ => { panic!("TreeGen: `length` not NodeValue::Length"); },

                };

                nodes.push(TreeNode { id, parent: parent_id, lchild: 0, rchild: 0, length });
                if loc == LCHILD {
                    nodes[parent_id].lchild = id;
                }
                else if loc == RCHILD {
                     nodes[parent_id].rchild = id;
                }

                match &pnodes[*sub].val { 
                    NodeValue::Leaf { name } => {
                        match &pnodes[*name].val {
                            NodeValue::Name(Some(s)) => { labels.push(Label::Tip {name: s.clone(), data: vec![], calibration: None}); },
                            NodeValue::Name(None) => { labels.push(Label::Tip {name: String::new(), data: vec![], calibration: None}); },
                            _ => { panic!("TreeGen: `name` not NodeValue::Name"); },
                        }
                    },
                    NodeValue::Internal { branches, name } => {
                        match &pnodes[*name].val {
                            NodeValue::Name(Some(s)) => { labels.push(Label::Interior {name: s.clone()}); },
                            NodeValue::Name(None) => { labels.push(Label::Interior {name: String::new()}); },
                            _ => { panic!("TreeGen: `name` not NodeValue::Name"); },
                        }

                        if branches.len() != 2 { return Err(String::from("TreeGen: Not a binary tree")); }
                        nstack.push((branches[0], id, LCHILD));
                        nstack.push((branches[1], id, RCHILD));
                    },
                    _ => { panic!("TreeGen: `sub` of Branch neither Leaf nor Internal") }
                };
            }
            t @ _ => { panic!("TreeGen: popped neither Tree nor Branch :: {:?}", t) },
        }
    }

    let tree = Tree { nodes: nodes };
    Ok((tree, labels))
}

pub fn parse(s: String) -> Result<(Tree, Labels), String> {
    let mut toks = tokenise(s)?;
    parse_tree(&mut toks)
}

