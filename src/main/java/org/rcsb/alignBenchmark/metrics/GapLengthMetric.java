/**
 * 
 */
package org.rcsb.alignBenchmark.metrics;

import java.util.Arrays;
import java.util.Map;

import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.StructureException;
import org.biojava.bio.structure.align.model.AFPChain;
import org.rcsb.alignBenchmark.MultipleAlignment;

/**
 * Calculates the gap penalty of the alignment.
 * 
 * Gap penalties are linear functions of the number and size of the gaps. The
 * score of a consists of a constant gap initiation penalty (i) and a linear 
 * gap extension penalty (e), which is applied for each residue in the gap.
 * The total score for a gap of lenght l is o+(l-1)*e.
 * 
 * If an alignment consists of multiple blocks, gaps between blocks are not included.
 * 
 * @author Spencer Bliven <sbliven@ucsd.edu>
 *
 */
public abstract class GapLengthMetric extends Metric{
	protected double gapExtensionPenalty;
	protected double gapInitiationPenalty;
	
	/**
	 * Create a metric with the specified gap initiation and extension values.
	 * @param gapInitiationPenalty Cost of the first residue in a gap. Generally negative.
	 * @param gapExtensionPenalty Cost of the second residue in a gap. Generally negative and smaller magnitude than the initiation.
	 */
	protected GapLengthMetric(double gapInitiationPenalty, double gapExtensionPenalty) {
		this.gapInitiationPenalty = gapInitiationPenalty;
		this.gapExtensionPenalty = gapExtensionPenalty;
	}
	
	/**
	 * Use a gapInitiation=1, gapExtension=0. This counts the number of gaps.
	 */
	protected GapLengthMetric() {
		this(1, 0);
	}
	
	/**
	 * score a single gap
	 * @param len Length of the gap
	 * @return gap penalty
	 */
	protected double scoreGap(int len) {
		if( len == 0 ) {
			return 0;
		}
		return gapInitiationPenalty + (len-1)*gapExtensionPenalty;
	}
	
	public static class Reference extends GapLengthMetric {
		public Reference(double gapInitiationPenalty, double gapExtensionPenalty) {
			super(gapInitiationPenalty,gapExtensionPenalty);
		}
		public Reference() {
			super();
		}
		/**
		 * @param reference
		 * @param align
		 * @param ca1
		 * @param ca2
		 * @param metaData
		 * @return the total gap penalty over all gaps, or NaN on error
		 * @see org.rcsb.alignBenchmark.metrics.Metric#calculate(org.rcsb.alignBenchmark.MultipleAlignment, org.biojava.bio.structure.align.model.AFPChain, org.biojava.bio.structure.Atom[], org.biojava.bio.structure.Atom[], java.util.Map)
		 */
		@Override
		public double calculate(MultipleAlignment reference, AFPChain ignored,
				Atom[] ca1, Atom[] ca2, Map<String, Object> metaData) {
			if(reference.size()<=1) {
				return 0.;
			}
			int[][] mat;
			double score = 0;
			try {
				// get alignment as 0-based residue indices
				mat = reference.getAlignmentMatrix(Arrays.asList(ca1,ca2));
			} catch (StructureException e) {
				// One of the aligned residues doesn't appear in the structure.
				return Double.NaN;
			}
			
			assert(mat.length == 2);
			assert(mat[0].length == reference.size());
			assert(mat[1].length == reference.size());


			for(int pos=1;pos<reference.size();pos++) {
				int gapLen=0;
				for(int prot=0; prot<mat.length; prot++) {
					if(mat[prot][pos] <= mat[prot][pos-1] ) {
						//circular permutation
						int[] protLens = new int[] {ca1.length, ca2.length};
						//gap is the remaining res plus the leading res
						gapLen = Math.max(gapLen, protLens[prot]-mat[prot][pos-1]-1 + mat[prot][pos]);
					}
					else if(mat[prot][pos] > mat[prot][pos-1]+1) {
						gapLen = Math.max(gapLen, mat[prot][pos]-mat[prot][pos-1]-1);
					}
				}
				//gapLen is now the maximum gap over all structures.
				score += scoreGap(gapLen);
			}
			
			return score;
		}

		/**
		 * @return
		 * @see org.rcsb.alignBenchmark.metrics.Metric#getName()
		 */
		@Override
		public String getName() {
			return String.format("Ref_gaps(%d,%d)",gapInitiationPenalty, gapExtensionPenalty);
		}
	}
	public static class Alignment extends GapLengthMetric {
		public Alignment(double gapInitiationPenalty, double gapExtensionPenalty) {
			super(gapInitiationPenalty,gapExtensionPenalty);
		}
		public Alignment() {
			super();
		}
		/**
		 * @param reference
		 * @param align
		 * @param ca1
		 * @param ca2
		 * @param metaData
		 * @return the total gap penalty over all gaps, or NaN on error
		 * @see org.rcsb.alignBenchmark.metrics.Metric#calculate(org.rcsb.alignBenchmark.MultipleAlignment, org.biojava.bio.structure.align.model.AFPChain, org.biojava.bio.structure.Atom[], org.biojava.bio.structure.Atom[], java.util.Map)
		 */
		@Override
		public double calculate(MultipleAlignment ignored, AFPChain align,
				Atom[] ca1, Atom[] ca2, Map<String, Object> metaData) {
			
			if(align.getOptLength()<=1) {
				return 0.;
			}
			
			int[][][] mat = align.getOptAln();
			int[] blockLen = align.getOptLen();
			double score = 0;
			
			for(int block=0;block<align.getBlockNum();block++) {
				for(int pos=1;pos<blockLen[block];pos++) {
					int gapLen=0;
					for(int prot=0; prot<2; prot++) {
						if(mat[block][prot][pos] <= mat[block][prot][pos-1] ) {
							//circular permutation
							int[] protLens = new int[] {ca1.length, ca2.length};
							//gap is the remaining res plus the leading res
							gapLen = Math.max(gapLen, protLens[prot]-mat[block][prot][pos-1]-1 + mat[block][prot][pos]);
						}
						else if(mat[block][prot][pos] > mat[block][prot][pos-1]+1) {
							gapLen = Math.max(gapLen, mat[block][prot][pos]-mat[block][prot][pos-1]-1);
						}
					}
					//gapLen is now the maximum gap over all structures.
					score += scoreGap(gapLen);
				}
			}
			
			return score;
		}

		/**
		 * @return
		 * @see org.rcsb.alignBenchmark.metrics.Metric#getName()
		 */
		@Override
		public String getName() {
			return String.format("Ref_gaps(%d,%d)",gapInitiationPenalty, gapExtensionPenalty);
		}
	}
	
	@Override
	public String format(double result) {
		return Integer.toString((int)result);
	}
}
