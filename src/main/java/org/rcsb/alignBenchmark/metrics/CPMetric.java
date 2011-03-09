package org.rcsb.alignBenchmark.metrics;


import java.util.HashMap;
import java.util.Map;

import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.align.model.AFPChain;
import org.rcsb.alignBenchmark.metrics.Metric;
import org.rcsb.alignBenchmark.MultipleAlignment;
import org.rcsb.alignBenchmark.PDBResidue;

public abstract class CPMetric {
	/**
	 * Calculates the length of the reference alignment
	 * @author Spencer Bliven
	 *
	 */
	public static class Reference extends Metric {

		@Override
		public double calculate(MultipleAlignment reference, AFPChain align, Atom[] ca1, Atom[] ca2, Map<String, Object> metaData) {
			return isSequential(reference) ? 0. : 1.;
		}
		
		/**
		 * Checks whether all residues in this MultipleAlignment are strictly increasing.
		 * 
		 * This returns false for circular permutants and other topological rearrangements.
		 * 
		 * Sequentiality is defined with in a chain and for a given insertion code.
		 * 
		 * @param ma
		 * @return
		 */
		public boolean isSequential(MultipleAlignment ma ) {
			for(PDBResidue[] sequence : ma.getAlignmentResidues() ) {
				// Map "Chain"+"insertionCode" -> residue number of last residue.
				Map<String,Integer> prevRes = new HashMap<String,Integer>();
				for(PDBResidue res: sequence) {
					String key = res.getChainId()+res.getInsCode();
					Integer next = res.getSeqNum();
					Integer prev = prevRes.get(key);
					
					if(prev != null && prev >= next) {
						// strict ordering not preserved.
						return false;
					}
					// in order, or first residue with this key
					prevRes.put(key, next);
				}
			}
	
			// Every residue was ordered.
			return true;
		}

		@Override
		public String getName() {
			return "Ref_CP";
		}

		@Override
		public String format(double result) {
			return Integer.toString((int)Math.round(result));
		}
	}
	
	/**
	 * Calculates the length of the test alignment
	 * @author Spencer Bliven
	 *
	 */
	public static class Alignment extends Metric {

		@Override
		public double calculate(MultipleAlignment reference, AFPChain align, Atom[] ca1, Atom[] ca2, Map<String, Object> metaData) {
			return align.isSequentialAlignment() ? 0. : 1.;
		}

		@Override
		public String getName() {
			return "Aln_CP";
		}
		
		@Override
		public String format(double result) {
			return Integer.toString((int)Math.round(result));
		}
	}
}
