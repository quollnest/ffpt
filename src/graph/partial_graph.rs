use crate::ff::fastframe::FastFrames;

use hashbrown::HashSet;

use bitvec::vec::BitVec;

use core::ops::BitOrAssign;

///
///
///
pub struct GraphBuffer {
    pub graph: Vec<Vec<(usize, Vec<usize>)>>
}

impl GraphBuffer {
    pub fn new() -> Self {
        GraphBuffer {
            graph: vec![Vec::new()]
        }
    }

    pub fn into_partial_order_graph(self) -> PartialOrderGraph {
        PartialOrderGraph { graph: self.graph }
    } 
}

///
///
///
pub struct GraphTraversalBuffer {
    pub remaining: Vec<(usize, Vec<usize>, Vec<usize>)>
}


impl GraphTraversalBuffer {
    pub fn new() -> Self {
        GraphTraversalBuffer {
            remaining: Vec::new()
        }
    }
}

///
///
///
#[derive(Debug)]
pub struct PartialOrderGraph { 
    pub graph: Vec<Vec<(usize, Vec<usize>)>>
}

impl PartialOrderGraph {

    pub fn len(&self) -> usize {
        self.graph.len()
    }
}




pub fn get_order(ffs: &FastFrames, 
    map: &[usize]) -> PartialOrderGraph {

    let width = ffs.get_stacks()
        .get_width();
    
    //Phase 1

    let mut gbuf = GraphBuffer::new();
    let mut rembuf = GraphTraversalBuffer::new(); 

    for qidx in 0..ffs.len() {
        let (xvec, zvec) = ffs.get_stacks().get_vecs();   
        
        let qloc = qidx * width;
        let qwidth = width;
        let qloc_end = qloc+qwidth;

        let z_stack = &zvec[qloc..qloc_end];
        let x_stack = &xvec[qloc..qloc_end];

        let vmax = z_stack.len().max(x_stack.len());
        
        let mut new_zvec = 
            BitVec::repeat(false, vmax);

        let mut new_xvec = 
            BitVec::repeat(false, vmax);
       
        new_zvec.copy_from_bitslice(z_stack); 
        new_xvec.copy_from_bitslice(x_stack);
        
        new_zvec.bitor_assign(new_xvec);


        let mut deps: HashSet<usize> = HashSet::new();
        for (dep, flag) in 
            new_zvec.iter().enumerate() {
                if *flag {
                    deps.insert(map[dep]);
                }
        }
        let deps = Vec::from_iter(deps.into_iter());

        if deps.is_empty() {
            gbuf.graph[0].push((qidx, deps));
        } else {
            rembuf.remaining.push((qidx, Vec::new(), deps));
        }

    }


    //Phase 2
    let mut layer_idx = 0;

    while !rembuf.remaining.is_empty() {
        let mut new_layer = Vec::new();
        for (known, deps) in 
            gbuf.graph
            .get(layer_idx).unwrap().iter() {

            let mut register = Vec::new();
            for (bit, 
                (_, resolved, open)) in 
                rembuf.remaining 
                .iter_mut().enumerate() {

                if let Some(resolved_idx) = 
                    open.iter().position(|&dep| dep == *known) {
                    let redundent_deps: Vec<usize> = resolved
                        .iter()
                        .enumerate()
                        .filter_map(|(i, dep)| {
                            if deps.contains(dep) {
                                Some(i)
                            } else {
                                None
                            }
                        })
                        .collect();
                    // want to remove the redundent deps; 
                        // the swap_remove works, because
                    // redundent_deps is sorted with increasing order
                    for redundent in redundent_deps.iter().rev() {
                        resolved.swap_remove(*redundent);
                    }
                    resolved.push(open.swap_remove(resolved_idx));
                    if open.is_empty() {
                        register.push(bit);
                    }
                }
            }
            for fully_resolved in register.iter().rev() {
                let (bit, deps, _) = 
                    rembuf.remaining.swap_remove(*fully_resolved);
                new_layer.push((bit, deps));
            }
        }

        assert!(
            !new_layer.is_empty(),
            "couldn't find qubit with resolved dependencies in layer {}",
            layer_idx + 1
        );

        gbuf.graph.push(new_layer);
        layer_idx += 1;
    }

    gbuf.into_partial_order_graph()
}
