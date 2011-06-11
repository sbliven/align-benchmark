package org.rcsb.alignBenchmark.metrics;

import java.util.Map;

import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.align.model.AFPChain;
import org.rcsb.alignBenchmark.MultipleAlignment;

public abstract class JiaCEScoreMetric extends Metric{
	protected Metric lenMetric;
	protected Metric gapMetric;
	protected Metric rmsdMetric;

	protected JiaCEScoreMetric(AlignmentLengthMetric lenMetric,GapLengthMetric gapMetric,RMSDMetric rmsdMetric) {
		this.lenMetric = lenMetric;
		this.gapMetric = gapMetric;
		this.rmsdMetric = rmsdMetric;
	}
	
	//TODO this is a useful method and should live in biojava.
	/**
	 * Calculates the ce_score, as defined in
	 * 
	 * Jia et al. J Comput Biol. (2004) 11:5 p787-799
	 * 
	 * The function is:
	 * 
	 * ce_score = rmsd/aligned_length^alpha * (1+ num_gap/aligned_length^beta)
	 * 
	 * where alpha = beta = 1.
	 * 
	 * @param rmsd
	 * @param aligned_length
	 * @param num_gap
	 * @return
	 */
	public static double jiaCEScore(double rmsd, int aligned_length, int num_gap) {
		return rmsd/aligned_length * (1.+ (double)num_gap/aligned_length);
	}
	
	/**
	 * Calculates the ce_score
	 * @param reference
	 * @param align
	 * @param ca1
	 * @param ca2
	 * @param metaData
	 * @return
	 * @see org.rcsb.alignBenchmark.metrics.RMSDMetric.Alignment#calculate(org.rcsb.alignBenchmark.MultipleAlignment, org.biojava.bio.structure.align.model.AFPChain, org.biojava.bio.structure.Atom[], org.biojava.bio.structure.Atom[], java.util.Map)
	 */
	@Override
	public double calculate(MultipleAlignment reference, AFPChain align,
			Atom[] ca1, Atom[] ca2, Map<String, Object> metaData) {
		double rmsd =  rmsdMetric.calculate(reference, align, ca1, ca2, metaData);
		int alignLen = (int)lenMetric.calculate(reference, align, ca1, ca2, metaData);
		int numGap = (int)gapMetric.calculate(reference, align, ca1, ca2, metaData);
		return jiaCEScore(rmsd,alignLen,numGap);
	}
	
	@Override
	public String format(double result) {
		return Integer.toString((int)result);
	}
	
	public static class Reference extends JiaCEScoreMetric{
		
		/**
		 * @param lenMetric
		 * @param gapMetric
		 * @param rmsdMetric
		 */
		public Reference() {
			super(new AlignmentLengthMetric.Reference(),
					new GapLengthMetric.Reference(),
					new RMSDMetric.Reference() );
		}

		/**
		 * @return
		 * @see org.rcsb.alignBenchmark.metrics.RMSDMetric.Alignment#getName()
		 */
		@Override
		public String getName() {
			// TODO Auto-generated method stub
			return "Ref_CEScore";
		}
		
	}
	
	public static class Alignment extends JiaCEScoreMetric{
		Metric lenMetric;
		Metric gapMetric;
		Metric rmsdMetric;
		
		/**
		 * @param lenMetric
		 * @param gapMetric
		 * @param rmsdMetric
		 */
		public Alignment() {
			super(new AlignmentLengthMetric.Alignment(),
					new GapLengthMetric.Alignment(),
					new RMSDMetric.Alignment() );
		}

		/**
		 * @return
		 * @see org.rcsb.alignBenchmark.metrics.RMSDMetric.Alignment#getName()
		 */
		@Override
		public String getName() {
			// TODO Auto-generated method stub
			return "Aln_CEScore";
		}
		
	}

	

}
